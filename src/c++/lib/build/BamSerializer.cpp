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
 ** \file BamSerializer.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include "build/BamSerializer.hh"

namespace isaac
{
namespace build
{

alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> BamSerializer::makeSplit(
    const alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator>& last,
    alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> current,
    PackedFragmentBuffer& data,
    PackedFragmentBuffer::Index& index,
    const unsigned editDistance,
    BinData::IndexType& splitIndexEntries,
    alignment::Cigar& splitCigars,
    SplitInfoList &splitInfoList)
{
    // create new entry
    PackedFragmentBuffer::Index secondPart(index);
    const std::size_t before = splitCigars.size();
    if (current.reverse_ == last.reverse_)
    {
        // soft clip sequence that belongs to the previous part of the split
        splitCigars.addOperation(last.sequenceOffset_, alignment::Cigar::SOFT_CLIP);
    }
    secondPart.pos_ = current.referencePos_;
    splitCigars.addOperations(current.cigarIt_, index.cigarEnd_);
    secondPart.cigarBegin_ = &splitCigars.front() + before;
    secondPart.cigarEnd_ = &splitCigars.back() + 1;
    secondPart.reverse_ = current.reverse_;

    //patch the old one
    const PackedFragmentBuffer::Index::CigarIterator oldBegin = index.cigarBegin_;
    index.cigarBegin_ = &splitCigars.back() + 1;
    splitCigars.addOperations(oldBegin, last.cigarIt_);
    const unsigned sequenceLeftover = data.getFragment(index).readLength_ - last.sequenceOffset_;
    if (sequenceLeftover)
    {
        // soft clip sequence that belongs to the next part of the split
        splitCigars.addOperation(sequenceLeftover, alignment::Cigar::SOFT_CLIP);
    }
    index.cigarEnd_ = &splitCigars.back() + 1;
    splitIndexEntries.push_back(index);
    splitInfoList.push_back(index);
    splitInfoList.back().editDistance_ = editDistance;

    index = secondPart;
    //change iterator to travel over the CIGAR of the new entry
    current = alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator>(
        index.cigarBegin_, index.cigarEnd_, index.pos_, current.reverse_, current.readLength_);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(data.getFragment(index).clusterId_, "Split into " << index << " and " << secondPart << " editDistance:" << editDistance);
    return current;
}

/**
 * \brief split into multiple bam segments if it cannot be represented as single one.
 */
void BamSerializer::splitIfNeeded(
    const reference::ContigList &contigList,
    PackedFragmentBuffer &data,
    PackedFragmentBuffer::Index &index,
    BinData::IndexType &splitIndexEntries,
    alignment::Cigar &splitCigars,
    SplitInfoList &splitInfoList)
{
    const BinData::IndexType::iterator before = splitIndexEntries.end();
    // update all index record positions before sorting as they may have been messed up by gap realignment
    io::FragmentAccessor &fragment = data.getFragment(index);
    index.pos_ = fragment.fStrandPosition_;
    alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> last(
        index.cigarBegin_, index.cigarEnd_, index.pos_, fragment.isReverse(), fragment.readLength_);
    unsigned splits = 0;
    unsigned lastEditDistance = 0;
    alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> current = last;
    for (; !current.end(); ++current)
    {
        if (current.referencePos_.getContigId() != last.referencePos_.getContigId() ||
            current.referencePos_ < last.referencePos_ ||
            current.reverse_ != last.reverse_ ||
            ((current.sequenceOffset_ - last.sequenceOffset_) != (current.referencePos_ - last.referencePos_) &&
                (current.referencePos_ - last.referencePos_) > splitGapLength_))
        {
            ISAAC_ASSERT_MSG(splitIndexEntries.size() < splitIndexEntries.capacity(), "New entries must not cause buffer reallocation");
            ISAAC_ASSERT_MSG(splitCigars.size() + std::distance(index.cigarBegin_, index.cigarEnd_) * 2 <= splitCigars.capacity(), "New entries must not cause cigar buffer reallocation");

            current = makeSplit(last, current, data, index, lastEditDistance, splitIndexEntries, splitCigars, splitInfoList);
            lastEditDistance = 0;

            ++splits;
        }
        else
        {
            lastEditDistance += computeEditDistance(
                contigList, last, current, fragment.isReverse(), fragment.basesBegin(), fragment.basesEnd());
        }

        last = current;
    }

    if (splits)
    {
        fragment.flags_.properPair_ = false;
        fragment.flags_.splitAlignment_ = true;
        io::FragmentAccessor& mate = data.getMate(index);
        mate.flags_.properPair_ = false;

        index.splitInfoOffset_ = splitInfoList.size() - splits;
        index.splitInfoCount_ = splits + 1;
        splitInfoList.push_back(index);
        splitInfoList.back().editDistance_ = computeEditDistance(
                contigList, last, current,
                fragment.isReverse(), fragment.basesBegin(), fragment.basesBegin() + fragment.readLength_);
        BOOST_FOREACH(PackedFragmentBuffer::Index &split, std::make_pair(before, splitIndexEntries.end()))
        {
            split.splitInfoOffset_ = index.splitInfoOffset_;
            split.splitInfoCount_ = index.splitInfoCount_;
        }
    }
}

void BamSerializer::prepareForBam(
    const reference::ContigList &contigList,
    PackedFragmentBuffer &data,
    BinData::IndexType &dataIndex,
    alignment::Cigar &splitCigars,
    SplitInfoList &splitInfoList)
{
    // Caution, we will be appending to dataIndex
    BOOST_FOREACH(PackedFragmentBuffer::Index &index, std::make_pair(dataIndex.begin(), dataIndex.end()))
    {
        splitIfNeeded(contigList, data, index, dataIndex, splitCigars, splitInfoList);
    }

    std::sort(dataIndex.begin(), dataIndex.end(), boost::bind(&PackedFragmentBuffer::orderForBam, boost::ref(data), _1, _2));
}

} // namespace build
} // namespace isaac
