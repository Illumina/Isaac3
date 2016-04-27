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
 ** \file FileSystem.cpp
 **
 ** file system helper utilities.
 **
 ** \author Roman Petrovski
 **/

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA

#include "common/Numa.hh"

namespace isaac
{
namespace common
{

namespace numa
{

static std::vector<unsigned> numaNodes;

/*
 * \brief both queries and
 */
static bool numaAvailable(bool set = false, bool enable = false)
{
    static bool set_ = false;
    static bool available_ = false;

    if (set)
    {
        ISAAC_VERIFY_MSG(!set_, "enabling numa is expected to be done once per lifetime of the process. available_ = " << available_);
        available_ = enable;
        set_ = true;
    }

    return available_;
}

// NB: __n is permitted to be 0.  The C++ standard says nothing
// about what the return value is when __n == 0.
void* numaAllocate(std::size_t size, const int node)
{
    if (!isNumaAvailable())
    {
        return ::operator new(size);
    }

#ifdef HAVE_NUMA
    void * ret = 0;
    if (numa::defaultNodeInterleave == node)
    {
        ret = numa_alloc_interleaved(size);
    }
    else if (numa::defaultNodeLocal == node)
    {
        ret = numa_alloc_local(size);
    }
    else
    {
//        numa_set_preferred(node);
        ret = numa_alloc_onnode(size, numa::numaNodes.at(node));
        numa_tonode_memory(ret, size, numa::numaNodes.at(node));
    }
    return ret;
#endif //HAVE_NUMA

    return ::operator new(size);
}

// __p is not permitted to be a null pointer.
void numaDeallocate(void * p, std::size_t size, const int node)
{
    if (!isNumaAvailable())
    {
        ::operator delete(p);
        return;
    }

#ifdef HAVE_NUMA
    ISAAC_THREAD_CERR << "numaDeallocate " << p << " of size " << size << "on node " <<
        ((numa::defaultNodeInterleave == node || numa::defaultNodeLocal == node) ? node : numa::numaNodes.at(node)) <<
        std::endl;
    numa_free(p, size);
    return;
#endif //HAVE_NUMA
    ::operator delete(p);
}

} // namespace numa


bool isNumaAvailable()
{
    return numa::numaAvailable();
}

bool numaInitialize(bool enable)
{
#ifdef HAVE_NUMA
    if (enable)
    {
        if (-1 == numa_available())
        {
            ISAAC_THREAD_CERR << "WARNING: numa library is unavailable while the binary is compiled to use numa. errno:" << errno << std::endl;
            return numa::numaAvailable(true, false);
        }
        else if (-1 == numa_max_node())
        {
            ISAAC_THREAD_CERR << "WARNING: numa_max_node returned -1 while the binary is compiled to use numa. errno:" << errno << std::endl;
            return numa::numaAvailable(true, false);
        }

        numa_exit_on_error = 1;
        numa_set_strict(1);

        const bitmask * nodes =  numa_get_run_node_mask();
        for (unsigned node = 0; node < nodes->size; ++node)
        {
            if (numa_bitmask_isbitset(nodes, node))
            {
                numa::numaNodes.push_back(node);
                ISAAC_THREAD_CERR << "numa allowed node " << numa::numaNodes.back() << std::endl;
            }
        }

        return numa::numaAvailable(true, true);
    }
#endif //HAVE_NUMA
    return numa::numaAvailable(true, false);
}


int getNumaNodeCount()
{
    if (!isNumaAvailable())
    {
        return 1;
    }

#ifdef HAVE_NUMA
    return numa::numaNodes.size();
#endif //HAVE_NUMA

    return 1;
}
int getThreadInterleaveNumaNode(const std::size_t threadNumber)
{
    const int runOnNode = threadNumber % getNumaNodeCount();
    ISAAC_ASSERT_MSG(8 * sizeof(uint64_t) >= unsigned(runOnNode), "numa node is too high: " << runOnNode);
    return runOnNode;
}

int bindCurrentThreadToNumaNode(const std::size_t numaNode)
{
    if (!isNumaAvailable())
    {
        return 0;
    }
#ifdef HAVE_NUMA
    const int runOnNode = numa::numaNodes.at(numaNode);
    ISAAC_VERIFY_MSG(-1 != numa_run_on_node(runOnNode), "numa_run_on_node " << runOnNode << " failed, errno: " << errno  << ":" << strerror(errno));

    ISAAC_THREAD_CERR << "bound thread to numa node " << runOnNode << std::endl;
    return runOnNode;
#endif //HAVE_NUMA

    return 0;
}


} // namespace common
} // namespace isaac
