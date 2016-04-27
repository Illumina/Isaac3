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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentMetadata.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/FragmentMetadata.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

static int64_t skipBackDels(const Cigar &cigarBuffer, const unsigned cigarOffset, const unsigned cigarLength, unsigned &i)
{
    int64_t ret = 0;
    for (;i < cigarLength; ++i)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[cigarOffset + i]);
        switch (cigar.second)
        {
            case Cigar::DELETE:
            {
                ret += cigar.first;
                break;
            }

            case Cigar::BACK:
            {
                ret -= cigar.first;
                break;
            }

            case Cigar::SOFT_CLIP:
            {
                return ret;
            }

            case Cigar::ALIGN:
            {
                return ret;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Reached unexpected CIGAR element before reaching ALIGN element." <<
                                 Cigar::toString(cigarBuffer, cigarOffset, cigarLength) << " got " << cigar.second << " at offset " << i);
            }
        }
    }
    ISAAC_ASSERT_MSG(false, "Reached the end of CIGAR before reaching ALIGN element." << Cigar::toString(cigarBuffer, cigarOffset, cigarLength));
    return 0;
}

unsigned FragmentMetadata::updateAnchors(
    unsigned sequenceOffset,
    const unsigned length,
    reference::Contig::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence,
    isaac::reference::ContigAnnotation::const_iterator currentAnnotation,
    Anchor& firstAnchor, Anchor& lastAnchor) const
{
    unsigned ret = 0;
    unsigned matchesInARow = 0;
    // scan backwards so that it's easier to check for k-uniqueness
    currentReference += length;
    sequenceOffset += length;
    currentSequence += length;
    currentAnnotation += length;
    for (unsigned j = 0; length > j; ++j)
    {
        --currentReference;
        --sequenceOffset;
        --currentSequence;
        --currentAnnotation;
        ISAAC_ASSERT_MSG(currentAnnotation->first, "currentAnnotation->first is 0 at sequenceOffset:" << sequenceOffset << " for " << *this)
        ISAAC_ASSERT_MSG(currentAnnotation->second, "currentAnnotation->second is 0 at sequenceOffset:" << sequenceOffset << " for " << *this)
        if (isMatch(*currentSequence, *currentReference))
        {
            ++matchesInARow;
            if (matchesInARow >= currentAnnotation->first)
            {
                if (firstAnchor.empty() || !firstAnchor.kUnique_ || firstAnchor.first > sequenceOffset)
                {
                    firstAnchor.first = sequenceOffset;
                    firstAnchor.second = sequenceOffset + currentAnnotation->first;
                    firstAnchor.kUnique_ = true;
                }
                if (lastAnchor.empty() || !lastAnchor.kUnique_ || lastAnchor.first < sequenceOffset)
                {
                    lastAnchor.first = sequenceOffset;
                    lastAnchor.second = sequenceOffset + currentAnnotation->first;
                    lastAnchor.kUnique_ = true;
                }
            }
            else if (matchesInARow >= currentAnnotation->second)
            {
                // continuously update first anchor to make it as close to start of read as possible
                if (!firstAnchor.empty() && !firstAnchor.kUnique_ &&
                    firstAnchor.first >= sequenceOffset && firstAnchor.second <= sequenceOffset + currentAnnotation->second)
                {
                    firstAnchor.first = sequenceOffset;
                    firstAnchor.second = sequenceOffset + currentAnnotation->second;
                    firstAnchor.kRepeat_ = true;
                }

                // stop updating as soon as we get one. We want the last one to be as far at the back as possible
                if (!lastAnchor.empty() && !lastAnchor.kUnique_ && !lastAnchor.kRepeat_ &&
                    lastAnchor.first >= sequenceOffset && lastAnchor.second <= sequenceOffset + currentAnnotation->second)
                {
                    lastAnchor.first = sequenceOffset;
                    lastAnchor.second = sequenceOffset + currentAnnotation->second;
                    lastAnchor.kRepeat_ = true;
                }
            }
        }
        else
        {
            ret = std::max(ret, matchesInARow);
            matchesInARow = 0;
        }
    }
    return ret;
}

double FragmentMetadata::calculateLogProbability(
    unsigned length,
    reference::Contig::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence,
    std::vector<char>::const_iterator currentQuality) const
{
    double ret = 0.0;
    while (length--)
    {
        if (isMatch(*currentSequence, *currentReference))
        {
            ret += Quality::getLogMatch(*currentQuality);
        }
        else
        {
            ret += Quality::getLogMismatch(*currentQuality);
        }
        ++currentReference;
        ++currentSequence;
        ++currentQuality;
    }
    return ret;
}

//double FragmentMetadata::calculateInsertionLogProbability(
//    unsigned length,
//    std::vector<char>::const_iterator currentQuality) const
//{
//    double ret  = std::accumulate(currentQuality, currentQuality + length, 0.0,
//                           bind(std::plus<double>(), _1, boost::bind(&Quality::getLogMatch, _2)));

// !!!!!!!!! this code is not being used but looks like it mistakenly does max instead of max_element!!!!!!!!!!!!!!!
//    const char qualityMax = *std::max(currentQuality, currentQuality + length);
//    // assume one highest-quality base mismatches
//    ret -= Quality::getLogMatch(qualityMax);
//    ret += Quality::getLogMismatch(qualityMax);
//    return ret;
//}

/**
 * \brief assume a deletion is at least worth a highest base quality mismatch, otherwise one-mismatch alignments
 *        get thrown away for zero-mismatch gapped madness
 */
double calculateDeletionLogProbability(
    std::vector<char>::const_iterator qualityBegin,
    std::vector<char>::const_iterator qualityEnd)
{
    const char qualityMax = *std::max_element(qualityBegin, qualityEnd);
    return Quality::getLogMismatch(qualityMax);
}

void FragmentMetadata::addMismatchCycles(
    std::vector<char>::const_iterator currentSequence,
    reference::Contig::const_iterator currentReference,
    unsigned sequenceOffset,
    unsigned length, bool reverse,
    const unsigned lastCycle,
    const unsigned firstCycle)
{
    while (length--)
    {
        if (!isMatch(*currentSequence, *currentReference))
        {
            addMismatchCycle(reverse ? lastCycle - sequenceOffset : firstCycle + sequenceOffset);
        }
        ++currentReference;
        ++sequenceOffset;
        ++currentSequence;
    }
}

inline void fixupAnchor(
    const unsigned beginClipOffset,
    const unsigned endClipOffset,
    Anchor& anchor)
{
    anchor.first = std::max<unsigned>(anchor.first, beginClipOffset);
    anchor.second = std::max(anchor.first, anchor.second);
    anchor.second = std::min<unsigned>(anchor.second, endClipOffset);
    anchor.first = std::min(anchor.first, anchor.second);
}

unsigned FragmentMetadata::updateAlignment(
    const bool collectMismatchCycles,
    const AlignmentCfg &cfg,
    const flowcell::ReadMetadata &readMetadata,
    const reference::ContigList &contigList,
    const isaac::reference::ContigAnnotations &contigAnnotations,
    unsigned contigId,
    const int64_t strandPosition,
    const Cigar &cigarBuffer,
    const unsigned cigarOffset)
{
    const Read &read = this->getRead();
    bool reverse = this->reverse;
    std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(reverse).begin();
    std::vector<char>::const_iterator qualityBegin = read.getStrandQuality(reverse).begin();

    ISAAC_ASSERT_MSG(!contigAnnotations.at(contigId).empty(), "Empty annotation unexpected for " << contigList.at(contigId));
    ISAAC_ASSERT_MSG(!contigList.at(contigId).empty(), "Reference contig was not loaded for " << *this);

    ISAAC_ASSERT_MSG(0 <= strandPosition, "position must be positive for CIGAR update " << *this << " strandPosition:" << strandPosition <<
                     " CIGAR: " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()));
    reference::Contig::const_iterator referenceBegin = contigList.at(contigId).begin();
    isaac::reference::ContigAnnotation::const_iterator annotationBegin = contigAnnotations.at(contigId).begin();

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " <<
//                                           Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()) << " for " << *this);


    const unsigned firstCycle = readMetadata.getFirstCycle();
    const unsigned lastCycle = readMetadata.getLastCycle();

    this->reverse = reverse;
    this->contigId = contigId;
    this->cigarBuffer = &cigarBuffer;
    this->cigarOffset = cigarOffset;
    this->cigarLength = cigarBuffer.size() - this->cigarOffset;
    // adjust cigarOffset and cigarLength
    ISAAC_ASSERT_MSG(cigarBuffer.size() > this->cigarOffset, "Expecting the new cigar is not empty");

    int64_t currentPosition = strandPosition;
    unsigned currentBase = 0;
    unsigned matchCount = 0;
    for (unsigned i = 0; this->cigarLength > i; ++i)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[this->cigarOffset + i]);
        const unsigned arg = cigar.first;
        const Cigar::OpCode opCode = cigar.second;
        if (opCode == Cigar::ALIGN)
        {
            ISAAC_ASSERT_MSG(std::size_t(currentPosition + arg) <= contigAnnotations.at(contigId).size(),
                             "overshooting the annotation contigId:" << contigId <<
                             " currentPosition:" << currentPosition <<
                             " arg:" << arg);
            // scan backwards so that it's easier to check for k-uniqueness
            this->matchesInARow = std::max(
                this->matchesInARow,
                updateAnchors(
                    currentBase,arg, referenceBegin + currentPosition, sequenceBegin + currentBase,
                    annotationBegin + currentPosition, firstAnchor_, lastAnchor_));

            this->logProbability += calculateLogProbability(
                arg, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);

            if (collectMismatchCycles)
            {
                addMismatchCycles(
                    sequenceBegin + currentBase, referenceBegin + currentPosition,
                    currentBase, arg, reverse, lastCycle, firstCycle);
            }

            const unsigned matches =
                std::inner_product(
                sequenceBegin + currentBase, sequenceBegin + currentBase + arg, referenceBegin + currentPosition,
                0, std::plus<unsigned>(), &isMatch);

            const unsigned mismatches = arg - matches;

            mismatchCount += mismatches;
            this->smithWatermanScore += cfg.normalizedMismatchScore_ * mismatches;

            matchCount += matches;

            // the edit distance includes all mismatches and ambiguous bases (Ns)
            this->editDistance += std::inner_product(
                sequenceBegin + currentBase, sequenceBegin + currentBase + arg, referenceBegin + currentPosition,
                0, std::plus<unsigned>(), std::not_equal_to<char>());

            currentPosition += arg;
            currentBase += arg;
        }
        else if (opCode == Cigar::INSERT)
        {
//            this->logProbability += calculateInsertionLogProbability(arg, qualityBegin + currentBase);
            this->logProbability += calculateLogProbability(
                arg, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);

            currentBase += arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
        }
        else if (opCode == Cigar::DELETE)
        {
            this->logProbability += calculateDeletionLogProbability(qualityBegin + currentBase, read.getStrandQuality(reverse).end());
            currentPosition += arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
            this->splitAlignment |= (arg > cfg.splitGapLength_);
        }
        else if (opCode == Cigar::BACK)
        {
            this->logProbability += calculateDeletionLogProbability(qualityBegin + currentBase, read.getStrandQuality(reverse).end());
            currentPosition -= arg;
            this->editDistance += arg;
            ++this->gapCount;
            this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (arg - 1) * cfg.normalizedGapExtendScore_);
            this->splitAlignment = true;
        }
        else if (opCode == Cigar::FLIP)
        {
            reverse = !reverse;
            sequenceBegin = read.getStrandSequence(reverse).begin();
            qualityBegin = read.getStrandQuality(reverse).begin();
            currentBase = read.getLength() - currentBase - arg;
            this->splitAlignment = true;
            // Notice, this will count flips followed by CONTIG or position adjustment as multiple gaps which is probably not
            // ideal, but so far the gapCount is not being used for anything that requires precise value.
            ++this->gapCount;
        }
        else if (opCode == Cigar::CONTIG)
        {
            ++i;
            currentPosition += skipBackDels(cigarBuffer, this->cigarOffset, this->cigarLength, i);
            ISAAC_ASSERT_MSG(0 <= currentPosition, "Unexpected negative position adjustment: " << currentPosition <<
                             " " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()))
            ISAAC_ASSERT_MSG(contigList.at(arg).size() > std::size_t(currentPosition), "Position adjustment outside the contig bounds: " << currentPosition <<
                             " " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()))
            referenceBegin = contigList.at(arg).begin();
            annotationBegin = contigAnnotations.at(arg).begin();
            contigId = arg;
            ++this->gapCount;
            --i;
            this->splitAlignment = true;
        }
        else if (opCode == Cigar::SOFT_CLIP)
        {
            // With inversions, soft clipping can occur in the middle of CIGAR
//            ISAAC_ASSERT_MSG(0 == i || i + 1 == this->cigarLength, "Soft clippings are expected to be "
//                "found only at the ends of cigar string");
            this->logProbability =
                std::accumulate(qualityBegin + currentBase, qualityBegin + currentBase + arg,
                                this->logProbability,
                                boost::bind(std::plus<double>(), _1, boost::bind(Quality::getLogMatch, _2)));

            // NOTE! Not advancing the reference for soft clips
            currentBase += arg;
        }
        else
        {
            using boost::format;
            const format message = format("Unexpected Cigar OpCode: %d") % opCode;
            BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
        }
    }
    this->observedLength = currentPosition - strandPosition;
    this->position = strandPosition;

    this->dodgy = !isWellAnchored();

    const unsigned endClipOffset = getReadLength() - getEndClippedLength();
    // Make sure empty are not crossing into the soft-clipped ends (including quality trimming).
    // This will mess up split read alignment
    fixupAnchor(getBeginClippedLength(), endClipOffset, firstAnchor_);
    fixupAnchor(getBeginClippedLength(), endClipOffset, lastAnchor_);

    // make life simpler for rest of the code.
    if (lastAnchor_.empty() && !firstAnchor_.empty())
    {
        lastAnchor_ = firstAnchor_;
    }
    if (firstAnchor_.empty() && !lastAnchor_.empty())
    {
        firstAnchor_ = lastAnchor_;
    }

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " << *this);

    return matchCount;
}

} // namespace alignment
} // namespace isaac
