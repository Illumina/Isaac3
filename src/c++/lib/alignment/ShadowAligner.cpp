/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file ShadowAligner.cpp
 **
 ** \brief See ShadowAligned.hh
 ** 
 ** \author Come Raczy
 **/

#include <algorithm>
#include <boost/format.hpp>

#include "alignment/Quality.hh"
#include "alignment/ShadowAligner.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "oligo/Kmer.hh"
#include "oligo/KmerGenerator.hpp"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

ShadowAligner::ShadowAligner(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const bool noSmithWaterman,
    const AlignmentCfg &alignmentCfg,
    const bool reserveBuffers)
    : gappedMismatchesMax_(gappedMismatchesMax),
      smitWatermanGapsMax_(smitWatermanGapsMax),
      noSmithWaterman_(noSmithWaterman),
      flowcellLayoutList_(flowcellLayoutList),
      ungappedAligner_(collectMismatchCycles, alignmentCfg),
      gappedAligner_(collectMismatchCycles, flowcellLayoutList, smartSmithWaterman, alignmentCfg)
{
    if (reserveBuffers)
    {
        shadowCigarBuffer_.reserve(Cigar::getMaxOperationsForReads(flowcellLayoutList) *
                       unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_);
        shadowCandidatePositions_.reserve(unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_);
    }
}

unsigned ShadowAligner::hashShadowKmers(const std::vector<char> &sequence)
{
    shadowKmerPositions_.clear();
    // initialize all k-mers to the magic value -1 (NOT_FOUND)
    shadowKmerPositions_.resize(shadowKmerCount_, -1);
    // 
    oligo::KmerGenerator<shadowKmerLength_, unsigned, std::vector<char>::const_iterator> kmerGenerator(sequence.begin(), sequence.end());
    unsigned positionsCount = 0;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    while (kmerGenerator.next(kmer, position))
    {
        if (-1 == shadowKmerPositions_[kmer])
        {
            shadowKmerPositions_[kmer] = (position - sequence.begin());
            ++positionsCount;
        }
    }
    return positionsCount;
}

void ShadowAligner::findShadowCandidatePositions(
    const reference::Contig::const_iterator referenceBegin,
    const reference::Contig::const_iterator referenceEnd,
    const std::vector<char> &shadowSequence)
{
    hashShadowKmers(shadowSequence);

    reference::Contig::const_iterator position = referenceBegin;

    while (referenceEnd != position)
    {
        // find matching positions in the reference by k-mer comparison
        oligo::KmerGenerator<shadowKmerLength_, unsigned, reference::Contig::const_iterator> kmerGenerator(position, referenceEnd);
        unsigned kmer;
        if(kmerGenerator.next(kmer, position))
        {
            if (-1 != shadowKmerPositions_[kmer])
            {
                const int64_t candidatePosition = position - referenceBegin - shadowKmerPositions_[kmer];
                // avoid spurious repetitions of start positions
                if (shadowCandidatePositions_.empty() || shadowCandidatePositions_.back() != candidatePosition)
                {
                    if (shadowCandidatePositions_.size() == shadowCandidatePositions_.capacity())
                    {
                        // too many candidate positions. Just stop here. The alignment score will be miserable anyway.
                        break;
                    }
                    shadowCandidatePositions_.push_back(candidatePosition);
//                ISAAC_THREAD_CERR << "kmer: " << oligo::Bases<2, oligo::Kmer>(kmer, shadowKmerLength_) << " " << candidatePosition << std::endl;
                }
            }
            position += std::min<std::size_t>(shadowKmerLength_, std::distance(position, referenceEnd));
        }
        else
        {
            break;
        }     
    }

    // remove duplicate positions
    if (!shadowCandidatePositions_.empty())
    {
        std::sort(shadowCandidatePositions_.begin(), shadowCandidatePositions_.end());
        shadowCandidatePositions_.erase(std::unique(shadowCandidatePositions_.begin(),
                                                    shadowCandidatePositions_.end()),
                                        shadowCandidatePositions_.end());
    }
}

/**
 * \brief if the best template is longer than the dominant template, attempt to rescue shadow
 *        within a wider range
 * \param bestTemplateLength 0 indicates there is no best template
 */
std::pair <int64_t, int64_t> calculateShadowRescueRange(
    const FragmentMetadata &orphan,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    const Cluster &cluster = orphan.getCluster();

    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const unsigned readLengths[] = {cluster[0].getLength(), cluster[1].getLength()};
    int64_t shadowMinPosition = templateLengthStatistics.mateMinPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths);
    int64_t shadowMaxPosition = templateLengthStatistics.mateMaxPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths) +
        readLengths[shadowReadIndex] - 1;

    const std::pair <int64_t, int64_t> ret(shadowMinPosition - 10, shadowMaxPosition + 10);
    return ret;
}

bool ShadowAligner::tryGapAlign(
    const flowcell::ReadMetadata& readMetadata,
    const matchSelector::FragmentSequencingAdapterClipper& adapterClipper,
    const reference::ContigList& contigList,
    const isaac::reference::ContigAnnotations& kUniqenessAnnotation,
    FragmentMetadata &fragment)
{
    // Use the gapped aligner if necessary
    if (BandedSmithWaterman::mismatchesCutoff < fragment.mismatchCount)
    {
        FragmentMetadata tmp = fragment;
        const unsigned matchCount = gappedAligner_.alignGapped(
            tmp, shadowCigarBuffer_, readMetadata, adapterClipper, contigList, kUniqenessAnnotation);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "    Rescuing:     Gap-aligned: " << tmp);
        if (matchCount && tmp.gapCount <= smitWatermanGapsMax_ && (
            tmp.smithWatermanScore < fragment.smithWatermanScore ||
            (tmp.smithWatermanScore == fragment.smithWatermanScore &&
                ISAAC_LP_LESS(fragment.logProbability, tmp.logProbability))))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "    Rescuing:     accepted: " << tmp); fragment = tmp;
            return true;
        }
    }
    return false;
}

/**
 * \return false when no reasonable placement for shadow found. If at that point the shadowList is not empty,
 * this means that the shadow falls at a repetitive region and rescuing should not be considered
 */
bool ShadowAligner::rescueShadow(
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const FragmentMetadata &orphan,
    FragmentMetadataList &shadowList,
    const flowcell::ReadMetadata &readMetadata,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    if (!templateLengthStatistics.isCoherent())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing impossible. Incoherent tls");
        return false;
    }
    shadowCigarBuffer_.reserve(Cigar::getMaxOperationsForReads(flowcellLayoutList_) * unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_);
    shadowCigarBuffer_.clear();
    ISAAC_ASSERT_MSG(2 > orphan.readIndex, "Paired reads means 2");
    const Cluster &cluster = orphan.getCluster();
    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const Read &shadowRead = cluster[shadowReadIndex];
    // identify the orientation and range of reference positions of the orphan
    const reference::Contig &contig = contigList[orphan.contigId];
    const bool shadowReverse = templateLengthStatistics.mateOrientation(orphan.readIndex, orphan.reverse);
    const std::pair<int64_t, int64_t> shadowRescueRange = calculateShadowRescueRange(orphan, templateLengthStatistics);
    if (shadowRescueRange.second < shadowRescueRange.first)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing impossible: shadowMaxPosition < shadowMinPosition "
            << shadowRescueRange.second << " < " << shadowRescueRange.first);
        return false;
    }
    if (shadowRescueRange.second + 1 + shadowRead.getLength() < 0)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing impossible: shadowMaxPosition + 1 + shadowRead.getLength() < 0 " <<
                                    (shadowRescueRange.second + 1 + shadowRead.getLength()) << "%l < 0");
        return false;
    }
    // find all the candidate positions for the shadow on the identified reference region
    shadowCandidatePositions_.reserve(unreasonablyHighDifferenceBetweenMaxAndMinInsertSizePlusFlanks_);
    shadowCandidatePositions_.clear();
    const std::vector<char> &shadowSequence = shadowReverse ? shadowRead.getReverseSequence() : shadowRead.getForwardSequence();
    const int64_t candidatePositionOffset = std::min((int64_t)contig.size(), std::max(0L, shadowRescueRange.first));
    findShadowCandidatePositions(
        contig.begin() + candidatePositionOffset,
        contig.begin() + std::min((int64_t)contig.size(), std::max(0L, shadowRescueRange.second) + 1),
        shadowSequence);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "findShadowCandidatePositions found " << shadowCandidatePositions_.size() << " positions in range [" <<
                                (candidatePositionOffset) << ";" <<
                                (std::min((int64_t)contig.size(), shadowRescueRange.second + 1)) << "]");
    // align the shadow to the candidate positions and keep the best fragment
    shadowList.clear();

    matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);

    FragmentMetadata *bestFragment = 0;
    BOOST_FOREACH(int64_t strandPosition, shadowCandidatePositions_)
    {
        if (shadowList.size() == shadowList.capacity())
        {
            return false;
        }
        strandPosition += candidatePositionOffset;
        FragmentMetadata fragment(&cluster, shadowReadIndex, 0, 0, shadowReverse, orphan.contigId, strandPosition);

        adapterClipper.checkInitStrand(fragment, contig);
        if (ungappedAligner_.alignUngapped(fragment, shadowCigarBuffer_, readMetadata, adapterClipper, contigList, kUniqenessAnnotation))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing: Aligned: " << fragment);
            shadowList.push_back(fragment);
            if (0 == bestFragment ||
                fragment.smithWatermanScore < bestFragment->smithWatermanScore ||
                (fragment.smithWatermanScore == bestFragment->smithWatermanScore &&
                    ISAAC_LP_LESS(bestFragment->logProbability, fragment.logProbability)))
            {
                bestFragment = &shadowList.back();
            }
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing: Unaligned: " << fragment);
        }
    }

    if (!bestFragment)
    {
        return false;
    }

    if(!noSmithWaterman_ && shadowList.size() >= 2 && BandedSmithWaterman::mismatchesCutoff < bestFragment->mismatchCount)
    {
        for (FragmentMetadataList::iterator firstCandidate = shadowList.begin(); shadowList.end() != firstCandidate; ++firstCandidate)
        {
            FragmentMetadataList::iterator prevCandidate = firstCandidate;
            for (FragmentMetadataList::iterator candidate = firstCandidate + 1;
                shadowList.end() != candidate  && candidate->position - prevCandidate->position < BandedSmithWaterman::WIDEST_GAP_SIZE;)
            {
                ++candidate;
                ++prevCandidate;
            }

            // if we have more than one candidate within the gap size, there is a chance we'll get a sensible gap alignment.
            if (firstCandidate != prevCandidate)
            {
                FragmentMetadata &fragment = *firstCandidate;
                if (tryGapAlign(readMetadata, adapterClipper, contigList, kUniqenessAnnotation, fragment))
                {
                    if (fragment.smithWatermanScore < bestFragment->smithWatermanScore ||
                        (fragment.smithWatermanScore == bestFragment->smithWatermanScore &&
                            ISAAC_LP_LESS(bestFragment->logProbability, fragment.logProbability)))
                    {
                        bestFragment = &fragment;
                    }
                }
            }

            // loop will advance the iterator
            firstCandidate = prevCandidate;
        }
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Not smithing: " << *bestFragment);
    }

    if (&shadowList.front() != bestFragment)
    {
        std::swap(shadowList.front(), *bestFragment);
    }
    return true;
}

} // namespace alignment
} // namespace isaac
