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
 ** \file FragmentAccessorBamAdapter.hh
 **
 ** Generates bam SA tag out of iSAAC internal CIGAR read and reference
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_SA_TAG_MAKER_HH
#define iSAAC_BUILD_SA_TAG_MAKER_HH

#include "alignment/Cigar.hh"
#include "alignment/Mismatch.hh"
#include "common/FastIo.hh"
#include "io/Fragment.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace build
{


struct SplitInfo
{
    SplitInfo(const PackedFragmentBuffer::Index &index):
        pos_(index.pos_), cigarBegin_(index.cigarBegin_), cigarEnd_(index.cigarEnd_), editDistance_(0)
    {
        pos_.setReverse(index.reverse_);
    }
    isaac::reference::ReferencePosition pos_;
    typedef const uint32_t * CigarIterator;
    CigarIterator cigarBegin_;
    CigarIterator cigarEnd_;
    unsigned editDistance_;
};

#ifndef _GLIBCXX_DEBUG
BOOST_STATIC_ASSERT(sizeof(SplitInfo) == 32);
#endif

typedef std::vector<SplitInfo> SplitInfoList;

template<typename ContainerT>
std::size_t beginSaPart(const isaac::reference::ReferencePosition& pos,
                 const bool reverse,
                 const isaac::reference::ContigList& contigList,
                 ContainerT& result)
{
    //    beginSAPart(pos, reverse, contigList);
    const std::string& contigName = contigList.at(pos.getContigId()).name_;
    const std::size_t ret = result.size();
    std::copy(contigName.begin(), contigName.end(), std::back_inserter(result));
    result.push_back(',');
    common::appendUnsignedNumber(result, pos.getPosition() + 1);
    result.push_back(',');
    result.push_back(reverse ? '-' : '+');
    result.push_back(',');
//    ISAAC_THREAD_CERR << "beginSaPart:" << std::string(result.begin(), result.end()) << std::endl;
    return ret;
}

template<typename IteratorT, typename ContainerT>
void endSaPart(IteratorT compoundCigarBegin, IteratorT compoundCigarEnd,
                         const unsigned softClipBack, const unsigned mapQ,
                         const unsigned currentEditDistance, ContainerT& result)
{
    alignment::Cigar::toString(compoundCigarBegin, compoundCigarEnd, result);
    if (softClipBack)
    {
        common::appendUnsignedNumber(result, softClipBack);
        result.push_back('S');
    }
    result.push_back(',');
    common::appendUnsignedNumber(result, mapQ);
    result.push_back(',');
    common::appendUnsignedNumber(result, currentEditDistance);
    result.push_back(';');
}

/**
 * \param sequenceReverse direction of sequence
 * \param alignedReverse direction of this alignment
 */
inline unsigned countMismatches(
    const reference::ContigList& contigList,
    const bool sequenceReverse,
    const unsigned char* alignedBegin,
    const unsigned char* alignedRBegin,
    const bool alignedReverse,
    const unsigned alignedBases,
    const reference::ReferencePosition& pos)
{
    if (sequenceReverse == alignedReverse)
    {
        return alignment::countEditDistanceMismatches(
            contigList,
            alignedBegin,
            pos,
            alignedBases);
    }
    else
    {
        return alignment::countEditDistanceMismatches(
            contigList,
            boost::make_transform_iterator(boost::make_reverse_iterator(alignedRBegin), &oligo::getReverseBcl),
            pos,
            alignedBases);
    }
}

template<typename IteratorT>
unsigned computeEditDistance(
    const isaac::reference::ContigList& contigList,
    const alignment::CigarPosition<IteratorT>& last,
    const alignment::CigarPosition<IteratorT>& current,
    const bool sequenceReverse,
    const unsigned char* sequenceBegin,
    const unsigned char* sequenceEnd)
{
    const unsigned sequenceMoved = current.sequenceOffset_ - last.sequenceOffset_;
    const unsigned referenceMoved = current.referencePos_ - last.referencePos_;
    if (sequenceMoved == referenceMoved)
    {
        return countMismatches(
            contigList,
            sequenceReverse,
            sequenceBegin + last.sequenceOffset_, sequenceEnd - last.sequenceOffset_,
            current.reverse_,
            sequenceMoved,
            last.referencePos_);
    }
    else
    {
        ISAAC_ASSERT_MSG(
            !sequenceMoved || !referenceMoved,
            "Unexpected unequal transition in CIGAR between reference and sequence");
        if (alignment::Cigar::SOFT_CLIP
            != alignment::Cigar::decode(*last.cigarIt_).second)
        {
            return sequenceMoved + referenceMoved;
        }
    }
    return 0;
}

template<typename ExcludeCigarIteratorT, typename CurrentCigarIteratorT>
bool notTheExcludeAlignment(
    reference::ReferencePosition excludePos, const bool excludeReverse,
    const ExcludeCigarIteratorT excludeCigarBegin,
    const ExcludeCigarIteratorT excludeCigarEnd,
    const reference::ReferencePosition& currentPos,
    const CurrentCigarIteratorT currentCigarBegin,
    const CurrentCigarIteratorT currentCigarEnd)
{
    const std::size_t excludeCigarLength = std::distance(excludeCigarBegin, excludeCigarEnd);
    ISAAC_ASSERT_MSG(!excludePos.reverse(), "Unexpected use of reverse bit on index positions. excludePos:" << excludePos << " currentPos:" << currentPos);
    excludePos.setReverse(excludeReverse);
    if (excludePos != currentPos)
    {
        return true;
    }
    const unsigned currentCigarLength = std::distance(currentCigarBegin, currentCigarEnd);

    return currentCigarLength != excludeCigarLength ||
        currentCigarEnd != std::mismatch(currentCigarBegin, currentCigarEnd, excludeCigarBegin).first;
}

/**
 * \brief Serializes cigar to a container in bam SA tag format
 *
 * \return true if alignment is the primary one
 */
template <typename CigarIteratorT, typename ContainerT>
bool makeSaTagString(
    const reference::ReferencePosition &excludePos,
    const bool excludeReverse,
    const CigarIteratorT excludeCigarBegin, const CigarIteratorT excludeCigarEnd,
    const unsigned mapQ,
    const isaac::reference::ContigList &contigList,
    const SplitInfoList &splitInfoList,
    const unsigned splitInfoOffset,
    unsigned splitInfoCount,
    ContainerT &resultBuffer,
    unsigned &excludeEditDistance)
{
    int ret = -1;
//    ISAAC_THREAD_CERR << "makeSaTagString excludeCigar:" << alignment::Cigar::toString(excludeCigarBegin, excludeCigarEnd) << std::endl;
    // processing backwards just because
    while (splitInfoCount--)
    {
        const SplitInfo &splitInfo = splitInfoList.at(splitInfoOffset + splitInfoCount);
        if (notTheExcludeAlignment(
            excludePos, excludeReverse, excludeCigarBegin, excludeCigarEnd, splitInfo.pos_, splitInfo.cigarBegin_, splitInfo.cigarEnd_))
        {
            beginSaPart(splitInfo.pos_, splitInfo.pos_.reverse(), contigList, resultBuffer);
            endSaPart(splitInfo.cigarBegin_, splitInfo.cigarEnd_, 0, mapQ, splitInfo.editDistance_, resultBuffer);
        }
        else
        {
            excludeEditDistance = splitInfo.editDistance_;
            ret = splitInfoCount;
        }
    }

    ISAAC_ASSERT_MSG(-1 != ret, "No match in splitInfoList found for: " << excludePos);
    // If we are not first in the list we are the supplementary one
    return !ret;
;
}

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_SA_TAG_MAKER_HH
