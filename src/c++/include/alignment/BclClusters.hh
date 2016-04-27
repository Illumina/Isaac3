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
 ** \file BclClusters.hh
 **
 ** \brief In-memory representation of sequencing data
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
#define iSAAC_ALIGNMENT_BCL_CLUSTERS_HH

#include <vector>

#include "common/Debug.hh"
#include "common/Numa.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{

struct ClusterXy
{
    ClusterXy() : x_(POSITION_NOT_SET), y_(POSITION_NOT_SET){}
    ClusterXy(int x, int y) : x_(x), y_(y) {}
    ClusterXy(const std::pair<int, int> &xy) : x_(xy.first), y_(xy.second) {}
    static const int POSITION_NOT_SET = boost::integer_traits<int>::const_max;

    int x_;
    int y_;

    bool isSet() const {return POSITION_NOT_SET != x_;}
};


template <typename AllocatorT>
class BasicBclClusters : std::vector<char, AllocatorT>
{
    typedef std::vector<char, AllocatorT> BaseT;
    std::size_t clusterLength_;
    std::vector<bool> pf_;
    std::vector<ClusterXy> xy_;

public:
    typedef typename BaseT::iterator iterator;
    typedef typename BaseT::const_iterator const_iterator;
    BasicBclClusters(unsigned clusterLength):
        clusterLength_(clusterLength)
    {
    }

    BasicBclClusters(const BasicBclClusters &that):
        BaseT(that),
        clusterLength_(that.clusterLength_),
        pf_(that.pf_),
        xy_(that.xy_)
    {
        BaseT::reserve(that.capacity());
        pf_.reserve(that.pf_.capacity());
        xy_.reserve(that.xy_.capacity());
    }

    void reserveClusters(const std::size_t reserveClusters, const bool storeXy)
    {
        BaseT::reserve(clusterLength_ * reserveClusters);
        pf_.reserve(reserveClusters);
        if (storeXy)
        {
            xy_.reserve(reserveClusters);
        }
    }

    void unreserve()
    {
        BaseT().swap(*this);
        std::vector<bool>().swap(pf_);
        std::vector<ClusterXy>().swap(xy_);
    }

    std::size_t getClusterCount() const
    {
        return BaseT::size() / clusterLength_;
    }

    unsigned getClusterLength() const
    {
        return clusterLength_;
    }

    /// appends more uninitialized clusters to the end. Returns iterator to first in the new block
    typename BaseT::iterator addMoreClusters(const std::size_t moreClusters)
    {
        BaseT::resize(BaseT::size() + moreClusters * clusterLength_);
        return BaseT::begin() + (BaseT::size() - moreClusters * clusterLength_);
    }
    /**
     * \post If the size of the buffer reduces, the data already in the buffer stays there.
     **/
    void reset(const std::size_t clusterLength, const std::size_t clusters, bool clear = false)
    {
        clusterLength_ = clusterLength;
        if (!clear)
        {
            BaseT::resize(clusterLength_ * clusters);
            // those storages that support pf will set it for all clusters.
            // The ones that don't, should assume all clusters are pf
            pf_.resize(clusters, true);
            if (storeXy())
            {
                xy_.resize(clusters);
            }
        }
        else
        {
            BaseT::clear();
            pf_.clear();
            xy_.clear();
        }
    }

    using BaseT::end;
    iterator cluster(std::size_t cluster)
    {
        return BaseT::begin() + cluster * clusterLength_;
    }
    const_iterator cluster(std::size_t cluster) const
    {
        return BaseT::begin() + cluster * clusterLength_;
    }

    std::vector<bool> &pf()
    {
        return pf_;
    }

    bool pf(std::size_t cluster) const
    {
        return pf_.at(cluster);
    }

    std::vector<ClusterXy> &xy()
    {
        return xy_;
    }

    const ClusterXy &xy(const std::size_t cluster) const
    {
        static const ClusterXy unsetXy;
        return storeXy() ? xy_.at(cluster) : unsetXy;
    }

    bool storeXy() const {return xy_.capacity();}

    void swap(BasicBclClusters &that)
    {
        ISAAC_ASSERT_MSG(clusterLength_ == that.clusterLength_, "Attempting to swap with buffer of different capacity " << *this << " " << that);
        ISAAC_ASSERT_MSG(BaseT::capacity() == that.capacity(), "Attempting to swap with buffer of different capacity " << *this << " " << that);
        ISAAC_ASSERT_MSG(pf_.capacity() == that.pf_.capacity(), "Attempting to swap with buffer of different capacity " << *this << " " << that);
        ISAAC_ASSERT_MSG(xy_.capacity() == that.xy_.capacity(), "Attempting to swap with buffer of different capacity " << *this << " " << that);
        BaseT::swap(that);
        pf_.swap(that.pf_);
        xy_.swap(that.xy_);
    }

    friend std::ostream &operator <<(std::ostream &os, const BasicBclClusters &bbc)
    {
        return os << "BasicBclClusters(" <<
            bbc.clusterLength_ << "cl," <<
            bbc.capacity() << "b," <<
            bbc.pf_.capacity() << "pf," <<
            bbc.xy_.capacity() << "xy)";
    }
};

typedef BasicBclClusters<common::NumaAllocator<char, common::numa::defaultNodeInterleave> > BclClusters;

template <typename IteratorType = BclClusters::const_iterator>
class BclClusterFields
{
    typedef IteratorType iterator;
    const flowcell::ReadMetadataList &readMetadataList_;
    const unsigned barcodeLength_;

public:
    typedef std::pair<iterator, iterator> IteratorPair;
    BclClusterFields(const flowcell::ReadMetadataList &readMetadataList,
                  const unsigned barcodeLength):
        readMetadataList_(readMetadataList), barcodeLength_(barcodeLength)
    {
    }

    IteratorPair getBarcode(const iterator clusterBegin) const
    {
        return std::make_pair(clusterBegin, clusterBegin + barcodeLength_);
    }

    iterator getBclBegin(const iterator clusterBegin, const unsigned readIndex) const
    {
        const iterator begin = clusterBegin + barcodeLength_ + readMetadataList_.at(readIndex).getOffset();
        return begin;
    }

    IteratorPair getBcl(const iterator clusterBegin, const unsigned readIndex) const
    {
        const iterator begin = getBclBegin(clusterBegin, readIndex);
        return std::make_pair(begin, begin + readMetadataList_[readIndex].getLength());
    }

    IteratorPair getName(const iterator clusterBegin,
                         const unsigned nameLengthMax) const
    {
        const iterator begin = clusterBegin + barcodeLength_ + readMetadataList_.back().getOffset() + readMetadataList_.back().getLength();
        return std::make_pair(begin, std::find(begin, begin + nameLengthMax, 0));
    }
};


} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BCL_CLUSTERS_HH
