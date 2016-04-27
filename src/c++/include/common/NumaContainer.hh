/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Public License 1
 **
 ** You should have received a copy of the Illumina Public License 1
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file FileSystem.hh
 **
 ** file system helper utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_NUMA_CONTAINER_HH
#define iSAAC_COMMON_NUMA_CONTAINER_HH

#include <boost/foreach.hpp>

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA

#include "common/Debug.hh"

#include "common/Numa.hh"
#include "common/Threads.hpp"

namespace isaac
{
namespace common
{

template <typename ReplicaT>
class NumaContainerReplicas
{
    std::vector<ReplicaT> nodeContainers_;
public:
    NumaContainerReplicas(ReplicaT &&node0Container)
    {
        if (common::isNumaAvailable())
        {
            const int nodes = getNumaNodeCount();

            nodeContainers_.reserve(nodes);
            ISAAC_THREAD_CERR << "before std::move(node0Container)" << std::endl;
            nodeContainers_.push_back(std::move(node0Container));
            ISAAC_THREAD_CERR << "after std::move(node0Container)" << std::endl;
            for (int node = 1; node < nodes; ++node)
            {
                ISAAC_THREAD_CERR << "before nodeContainers_.push_back()" << std::endl;
                nodeContainers_.push_back(ReplicaT(nodeContainers_.front(), common::NumaAllocator<void, 0> (node)));
                ISAAC_THREAD_CERR << "after nodeContainers_.push_back()" << std::endl;
            }
        }
        else
        {
            nodeContainers_.push_back(std::move(node0Container));
        }
    }

    const ReplicaT &node0Container() const {return nodeContainers_.front();}
    const ReplicaT &threadNodeContainer() const {return nodeContainers_.at(common::ThreadVector::getThreadNumaNode());}
//        return hashes_[(common::ThreadVector::getThreadNumaNode()+1) % 2].findMatches(kmer);
};

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_NUMA_CONTAINER_HH
