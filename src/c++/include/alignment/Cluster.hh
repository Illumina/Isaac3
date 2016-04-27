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
 ** \file Cluster.hh
 **
 ** \brief Component containing the data associated to a cluster: sequence and
 ** quality strings for all the reads in the cluster.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_CLUSTER_HH
#define iSAAC_ALIGNMENT_CLUSTER_HH

#include <vector>
#include <string>

#include "alignment/BclClusters.hh"
#include "alignment/Read.hh"
#include "common/FastIo.hh"

namespace isaac
{


namespace alignment
{

class Cluster: public std::vector<Read>
{
public:
    Cluster(const unsigned maxReadLen);

    void init(
        const flowcell::ReadMetadataList &readMetadataList,
        BclClusters::const_iterator bclData,
        const unsigned tile,
        const uint64_t id,
        const ClusterXy &xy,
        const bool pf,
        const unsigned barcodeLength,
        const unsigned nameLength);
    uint64_t getTile() const {return tile_;}
    uint64_t getId() const {return id_;}
    void setId(const uint64_t id) {id_ = id;}
    bool getPf() const {return pf_;}
    ClusterXy getXy() const {return xy_;}
    int getBarcodeLength() const { return barcodeLength_; }
    unsigned getNonEmptyReadsCount() const {return nonEmptyReads_;}
    /// Beginning of the BCL data for the indicated read
    BclClusters::const_iterator getBclData(unsigned readIndex) const;
    uint64_t getBarcodeSequence() const;

    BclClusters::const_iterator nameBegin() const
    {
        return readNameBegin_;
    }

    BclClusters::const_iterator nameEnd() const
    {
        return readNameEnd_;
    }

    template<class InpuT> friend InpuT& operator >>(InpuT &input, Cluster &cluster);
private:
    unsigned tile_;
    uint64_t id_;
    ClusterXy xy_;
    bool pf_;
    unsigned barcodeLength_;
    unsigned nonEmptyReads_;
    BclClusters::const_iterator bclData_;
    BclClusters::const_iterator readNameBegin_;
    BclClusters::const_iterator readNameEnd_;

    friend std::ostream &operator <<(std::ostream &os, const Cluster& cluster)
    {
        return os << "Cluster(" <<
            cluster.id_ << "," <<
            cluster.tile_ << "," <<
            cluster.pf_ << "," <<
            common::makeFastIoString(cluster.nameBegin(), cluster.nameEnd()) <<
            ")";
    }
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_CLUSTER_HH
