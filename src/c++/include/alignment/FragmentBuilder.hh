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
 ** \file FragmentBuilder.hh
 **
 ** \brief Utility classes for Fragment building and management for several reads
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH

#include "alignment/fragmentBuilder/GappedAligner.hh"
#include "alignment/fragmentBuilder/UngappedAligner.hh"
#include "alignment/fragmentBuilder/SplitReadAligner.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/HashMatchFinder.hh"
#include "alignment/Match.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/matchSelector/FragmentSequencingAdapterClipper.hh"
#include "reference/Contig.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{

class Cluster;

/**
 ** \brief Utility component creating and scoring all Fragment instances from a
 ** list Seed Matches for a single Cluster (each Read independently).
 **/
class FragmentBuilder: public boost::noncopyable
{
public:
    FragmentBuilder(
        const bool collectMismatchCycles,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const unsigned maxSeedsPerRead,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const bool avoidSmithWaterman,
        const bool noSmithWaterman,
        const bool splitAlignments,
        const AlignmentCfg &alignmentCfg,
        Cigar &cigarBuffer,
        const bool reserveBuffers);

    enum AlignmentType
    {
        // the order reflects the precedence meaning that when merging results of two read alignment attemps,
        // Normal will override any other value and so on.
        Normal,
        Rm,         // one of the seeds exactly mapped to a high repeat or too many neighbors with the same prefix
        Qc,         // all seeds contain Ns, alignment is not possible
        Nm,       // seeds have no match in the reference.
    };

    friend AlignmentType &combineAlignmentTypes(AlignmentType &left, const AlignmentType right)
    {
        left = std::min(left, right);
        return left;
    }

    friend AlignmentType &operator |= (AlignmentType &left, const AlignmentType right)
    {
        left = std::min(left, right);
        return left;
    }

    template <typename MatchFinderT>
    AlignmentType build(
        const reference::ContigList &contigList,
        const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
        const flowcell::ReadMetadata &readMetadata,
        const SeedMetadataList &seedMetadataList,
        const matchSelector::SequencingAdapterList &sequencingAdapters,
        const TemplateLengthStatistics &templateLengthStatistics,
        const MatchFinderT &matchFinder,
        const Cluster &cluster,

        bool withGaps,
        FragmentMetadataList &fragments);

    const Cigar &getCigarBuffer() const {return cigarBuffer_;}

    struct SequencingAdapterRange
    {
        SequencingAdapterRange() : defined_(false), empty_(true){}
        bool defined_;
        bool empty_;
        std::vector<char>::const_iterator adapterRangeBegin_;
        std::vector<char>::const_iterator adapterRangeEnd_;
    };

private:
    static const unsigned READS_MAX = 2;
    const unsigned repeatThreshold_;
    const unsigned gappedMismatchesMax_;
    const unsigned smitWatermanGapsMax_;
    const bool noSmithWaterman_;
    const bool splitAlignments_;

    const AlignmentCfg &alignmentCfg_;

    /**
     * \brief flag per seed indicating whether the seed matches are ignored due to
     *        a high repeat match
     */
    Cigar &cigarBuffer_;

    fragmentBuilder::UngappedAligner ungappedAligner_;
    fragmentBuilder::GappedAligner gappedAligner_;
    fragmentBuilder::SplitReadAligner splitReadAligner_;

    Matches matches_;
    typedef std::vector<unsigned char> OffsetMismatchFlags;
    OffsetMismatchFlags offsetMismatches_;

    /**
     ** \brief add a match, either by creating a new instance of
     ** FragmentMetadata or by updating an existing one
     **
     ** Initializes the FragmentMetadata in the list for the corresponding
     ** readIndex with contigId, orientation (reverse flag) and
     ** position. The fragment is initially located at the leftmost position of
     ** the read on the forward strand of the contig. This means that the
     ** position can be negative.
     **
     ** Note: spurious FragmentMetadata creation is avoided by checking if the
     ** last FragmentMetadata created for the read has same contigId, position
     ** and orientation.
     **/
    void addMatch(
        const flowcell::ReadMetadata &readMetadata,
        const SeedMetadataList &seedMetadataList, const Match &match,
        const Cluster &cluster,
        FragmentMetadataList &fragments);

    /**
     ** \brief Position of the leftmost base of a read on the forward strand,
     ** given a seed, its position and orientation.
     **
     ** For forward matches, the offset to apply is simply the seed offset
     ** indicated in the SeedMetadata. For reverse matches, the offset to apply
     ** is the remaining length of the read after the end of the seed.
     **/
    int64_t getReadPosition(
        const flowcell::ReadMetadata &readMetadata,
        const SeedMetadata &seedMetadata, const SeedId &seedId,
        const reference::ReferencePosition &seedLocation,
        const bool reverse) const;

    /// consolidate fragments with same reference position and orientation for a single read
    static void consolidateDuplicateAlignments(
        FragmentMetadataList &fragmentList,
        const bool removeUnaligned,
        std::size_t &perfectFound);
};

/**
 * \return true if at least one fragment was built.
 */
template <typename MatchFinderT>
FragmentBuilder::AlignmentType FragmentBuilder::build(
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &kUniqenessAnnotation,
    const flowcell::ReadMetadata &readMetadata,
    const SeedMetadataList &seedMetadataList,
    const matchSelector::SequencingAdapterList &sequencingAdapters,
    const TemplateLengthStatistics &templateLengthStatistics,
    const MatchFinderT &matchFinder,
    const Cluster &cluster,
    const bool withGaps,
    FragmentMetadataList &fragments)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: cluster " << cluster.getId());
    ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() > readMetadata.getIndex(), "cluster geometry must match");

    matchSelector::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);

    std::size_t repeatSeeds = 0;

    offsetMismatches_.clear();
    offsetMismatches_.resize(readMetadata.getLength(), false);
    bool firstSeed = true;

    BOOST_FOREACH(const SeedMetadata &seedMetadata, seedMetadataList)
    {
        if (seedMetadata.getReadIndex() != readMetadata.getIndex())
        {
            continue;
        }

        const OffsetMismatchFlags::const_iterator seedBegin = offsetMismatches_.begin() + seedMetadata.getOffset();
        const OffsetMismatchFlags::const_iterator seedEnd = seedBegin +
            std::min<std::size_t>(seedMetadata.getLength() * 2,
                     std::distance<OffsetMismatchFlags::const_iterator>(seedBegin, offsetMismatches_.end()));

        if(firstSeed || seedEnd != std::find(seedBegin, seedEnd, true))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    trying " << seedMetadata);

            matches_.clear();
            matchFinder.findSeedMatches(cluster, seedMetadata, readMetadata, matches_, repeatSeeds);

            BOOST_FOREACH(const Match &match, matches_)
            {
                if (!match.isTooManyMatch())
                {
                    firstSeed = false;
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    added r" << readMetadata.getIndex() << " match " << match);

                    const reference::ReferencePosition &seedLocation = match.location_;
                    const bool reverse = match.seedId_.isReverse() != match.location_.reverse();
                    const int64_t readPosition = getReadPosition(readMetadata, seedMetadata, match.getSeedId(), seedLocation, reverse);

                    {
                        FragmentMetadata fragmentMetadata(&cluster, readMetadata.getIndex(), 0, 0, reverse, seedLocation.getContigId(), readPosition);

                        // set seed anchors only for alignments where we've visited all candidates.
                        fragmentMetadata.setSeedAnchor(readMetadata, seedMetadata, match.getSeedId());

                        adapterClipper.checkInitStrand(fragmentMetadata, contigList.at(fragmentMetadata.contigId));
                        if (ungappedAligner_.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadata, adapterClipper, contigList, kUniqenessAnnotation))
                        {
                            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    Aligned    : " << fragmentMetadata);

                            fragments.push_back(fragmentMetadata);

                            if (fragmentMetadata.getMismatchCount())
                            {
                                const std::vector<char> &strandSequence = fragmentMetadata.getRead().getStrandSequence(reverse);
                                const std::vector<char>::const_iterator sequenceBegin = strandSequence.begin();

                                std::vector<char>::const_iterator sequenceIterator =
                                    sequenceBegin + fragmentMetadata.getBeginClippedLength();

                                // can't use seed offset for reverse sequenceEnd because sometimes the seed is inside the left quality trimming soft clip
                                const std::vector<char>::const_iterator sequenceEnd =
                                    strandSequence.end() - fragmentMetadata.getEndClippedLength();

                                reference::Contig::const_iterator referenceIterator =
                                    contigList.at(fragmentMetadata.contigId).begin() + fragmentMetadata.position;

                                for (std::size_t offset = std::distance(sequenceBegin, sequenceIterator);
                                    sequenceIterator != sequenceEnd; ++sequenceIterator, ++referenceIterator, ++offset)
                                {
                                    ISAAC_ASSERT_MSG(sequenceIterator >= sequenceBegin, "Sequence underrun " << fragmentMetadata << " offset:" << offset << " " << seedMetadata);
                                    ISAAC_ASSERT_MSG(sequenceIterator < strandSequence.end(), "Sequence overrun " << fragmentMetadata << " offset:" << offset << " " << seedMetadata);
                                    ISAAC_ASSERT_MSG(referenceIterator >= contigList.at(fragmentMetadata.contigId).begin(), "Reference underrun " << fragmentMetadata << " offset:" << offset << " " << seedMetadata);
                                    ISAAC_ASSERT_MSG(referenceIterator < contigList.at(fragmentMetadata.contigId).end(), "Reference overrun " << fragmentMetadata << " offset:" << offset << " " << seedMetadata);
                                    offsetMismatches_.at(reverse ? readMetadata.getLength() - offset - 1 : offset) |= !isMatch(*sequenceIterator, *referenceIterator);
    //                                ISAAC_ASSERT_MSG(isMatch(*sequenceIterator, *referenceIterator) || fragmentMetadata.getMismatchCount(), "mismatch found seq:" <<
    //                                                 common::makeFastIoString(sequenceIterator, sequenceEnd) << " ref:" <<
    //                                                 common::makeFastIoString(referenceIterator, referenceIterator + std::distance(sequenceIterator, sequenceEnd)) << " " <<
    //                                                 fragmentMetadata);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (fragments.empty())
    {
        return repeatSeeds ? FragmentBuilder::Rm : FragmentBuilder::Nm;
    }

	std::size_t perfectFound = 0;
    // make sure all positions that are in the list are unique.
    consolidateDuplicateAlignments(fragments, true, perfectFound);

    // having perfect alignments means no need to spend time on trying to improve the imperfect ones.
    if (withGaps && !perfectFound)
    {
        if (splitAlignments_)
        {
            splitReadAligner_.alignSimpleSv(cigarBuffer_, contigList, kUniqenessAnnotation, readMetadata, templateLengthStatistics, fragments);
            consolidateDuplicateAlignments(fragments, true, perfectFound);
        }

        if (withGaps && !noSmithWaterman_)
        {
            // If there are still bad alignments, try to do expensive smith-waterman on them.
            gappedAligner_.realignBadUngappedAlignments(
                gappedMismatchesMax_, smitWatermanGapsMax_, contigList, kUniqenessAnnotation, readMetadata, fragments, adapterClipper, cigarBuffer_);

            // gapped alignment and adapter trimming may have adjusted the alignment position
            consolidateDuplicateAlignments(fragments, true, perfectFound);
        }
    }

    return FragmentBuilder::Normal;
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
