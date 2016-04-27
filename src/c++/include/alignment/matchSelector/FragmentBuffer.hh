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
 ** \file FragmentCollector.hh
 **
 ** \brief Buffering of fragments with indexing information for BAM generation.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH

#include <boost/filesystem.hpp>
#include <boost/integer/static_min_max.hpp>

#include "alignment/FragmentMetadata.hh"
#include "alignment/BamTemplate.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BinIndexMap.hh"
#include "io/Fragment.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief buffer of fixed-size records capable of holding all cluster fragments
 */
class FragmentBuffer
{
public:
    typedef std::vector<char>::iterator iterator;
    typedef std::vector<char>::const_iterator const_iterator;
    FragmentBuffer(const flowcell::FlowcellLayoutList &flowcellLayoutList):
        recordLength_(getRecordLength(flowcellLayoutList)),
        readOffsets_(getReadOffsets(flowcellLayoutList)), // offset of the first read is 0
        clusters_(0)
    {
        ISAAC_THREAD_CERR << "Constructed FragmentBuffer for Record length: "
            << recordLength_ << " read offsets : " << readOffsets_[0] << "," <<
            (READS_MAX == readOffsets_.size() ? readOffsets_[1] : 0)  <<  std::endl;
    }

    void resizeMin(const uint64_t clusters)
    {
        clusters_ = clusters;
        data_.resize(std::max<std::size_t>(clusters_ * recordLength_, data_.size()));
        ISAAC_THREAD_CERR << "Resized FragmentBuffer for " << clusters_ << " clusters. Record length: "
            << recordLength_ << " read offsets : " << readOffsets_[0] << "," <<
            (READS_MAX == readOffsets_.size() ? readOffsets_[1] : 0)  <<  std::endl;
    }

    void reserve(const uint64_t clusters)
    {
        // data buffer is pre-allocated as the fragments will be put in there by mutiple threads using cluster id
        // as target location
        data_.reserve(clusters * recordLength_);
        ISAAC_THREAD_CERR << "Reserved FragmentBuffer for " << clusters << " clusters. Record length: "
            << recordLength_ << " read offsets : " << readOffsets_[0] << "," <<
            (READS_MAX == readOffsets_.size() ? readOffsets_[1] : 0)  <<  std::endl;
    }

    void clear()
    {
        data_.clear();
    }

    void swap(FragmentBuffer &another)
    {
        ISAAC_ASSERT_MSG(recordLength_ == another.recordLength_, "Buffers must be formatted identically");
        ISAAC_ASSERT_MSG(readOffsets_.size() == another.readOffsets_.size(), "Buffers must be formatted identically");
        ISAAC_ASSERT_MSG(readOffsets_[0] == another.readOffsets_[0], "Read offsets must match");
        ISAAC_ASSERT_MSG(readOffsets_.size() == 1 || readOffsets_[1] == another.readOffsets_[1], "Read offsets must match");
        std::swap(clusters_, another.clusters_);
        data_.swap(another.data_);
    }

    iterator getRecordInsertIterator(const uint64_t clusterId, const unsigned readIndex)
    {
        return data_.begin() + clusterId * recordLength_ + readOffsets_[readIndex];
    }

    uint64_t getClusters() const
    {
        return clusters_;
    }

    const_iterator dataBegin() const
    {
        return data_.begin();
    }

    const_iterator dataEnd() const
    {
        return data_.begin() + clusters_ * recordLength_;
    }

    const_iterator nexCluster(std::vector<char>::const_iterator it) const
    {
        return it + recordLength_;
    }

    const_iterator advanceClusters(std::vector<char>::const_iterator it, const std::size_t clusters) const
    {
        return it + std::min<std::size_t>(std::distance(it, dataEnd()), recordLength_ * clusters);
    }

    const_iterator nextRead(std::vector<char>::const_iterator it) const
    {
        ISAAC_ASSERT_MSG(2 == readOffsets_.size(), "Invalid requries for next read in single-ended data");
        return it + readOffsets_[1];
    }

private:
    static const unsigned READS_MAX = 2;
    const unsigned recordLength_;
    const common::StaticVector<unsigned, READS_MAX> readOffsets_;
    uint64_t clusters_;
    std::vector<char> data_;

    static unsigned getRecordLength(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        const unsigned read2MaxLength = flowcell::getMaxReadLength(flowcellLayoutList, 1);
        const unsigned clusterNameMax = flowcell::getMaxClusterName(flowcellLayoutList);
        return
            io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0), clusterNameMax) +
            (read2MaxLength ? (io::FragmentHeader::getMaxTotalLength(read2MaxLength, clusterNameMax)) : 0);
    }

    /**
     * \brief First read is located at the beginning, Second is at io::FragmentHeader::getMaxTotalLength of the
     *        first read length. Only two reads are supported
     */
    static common::StaticVector<unsigned, 2> getReadOffsets(const flowcell::FlowcellLayoutList &flowcellLayoutList)
    {
        common::StaticVector<unsigned, 2> ret;
        ret.push_back(0);
        const unsigned readsMax = flowcell::getMaxReadCount(flowcellLayoutList);
        const unsigned clusterNameMax = flowcell::getMaxClusterName(flowcellLayoutList);
        if (READS_MAX == readsMax)
        {
            ret.push_back(io::FragmentHeader::getMaxTotalLength(flowcell::getMaxReadLength(flowcellLayoutList, 0), clusterNameMax));
        }
        else
        {
            ISAAC_ASSERT_MSG(1 == readsMax, "Unexpected reads count " << readsMax);
        }
        return ret;
    }

};


} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_COLLECTOR_HH
