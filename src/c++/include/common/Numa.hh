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
 ** \file FileSystem.hh
 **
 ** file system helper utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_NUMA_ALLOCATOR_HH
#define iSAAC_COMMON_NUMA_ALLOCATOR_HH

#include <iostream>

#include "common/Debug.hh"

namespace isaac
{
namespace common
{

namespace numa
{

void* numaAllocate(std::size_t size, const int node);
void numaDeallocate(void * __p, std::size_t size, const int node);

static const int defaultNodeLocal = -1;
static const int defaultNodeInterleave = -2;

} //namespace numa

/**
 * \brief attempts to initialize NUMA-aware memory management.
 *
 * @param enable    set to false if regular memory allocation should be performed instead
 * @return          isNumaAvailable()
 */
bool numaInitialize(bool enable);

/**
 * \brief   Call this once at the process startup
 * @return  Returns 'false' if compiled without NUMA support, NUMA initialization failed or numaInitialize(false) was called
 */
bool isNumaAvailable();


/**
 * @param   threadNumber application-specific thread number.
 *          Expected to be uniformly distributed and unique within the group of thread being bound
 * @return  NUMA node to which the thread would be bound by bindCurrentThreadToNumaNode or 0 when NUMA is disabled
 */
int getThreadInterleaveNumaNode(const std::size_t threadNumber);

/**
 * \brief Binds the current thread to a numa node using modulo operator to achieve even distribution.
 *
 * @param   threadNumber application-specific thread number.
 *          Expected to be uniformly distributed and unique within the group of thread being bound
 * @return  NUMA node to which the thread is bound or 0 when NUMA is disabled
 */
int bindCurrentThreadToNumaNode(const std::size_t threadNumber);

/**
 * @return number of numa nodes or 1 if numa is unavailable
 */
int getNumaNodeCount();

template<typename Tp, int defaultNode = numa::defaultNodeLocal>
class NumaAllocator
{
    int node_;

public:
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;
    typedef Tp*       pointer;
    typedef const Tp* const_pointer;
    typedef Tp&       reference;
    typedef const Tp& const_reference;
    typedef Tp        value_type;

    template<typename Tp1>
    struct rebind { typedef NumaAllocator<Tp1, defaultNode> other; };

    NumaAllocator() throw() :node_(defaultNode) { }
    explicit NumaAllocator(const int node) throw() :node_(node) { }

    NumaAllocator(const NumaAllocator& that) throw() :node_(that.node_) { }

    template<typename Tp1>
    NumaAllocator(const NumaAllocator<Tp1, defaultNode>& that) throw() :node_(that.node_) { }

    ~NumaAllocator() throw() { }

    pointer address(reference __x) const { return &__x; }

    const_pointer address(const_reference __x) const { return &__x; }

    // NB: n is permitted to be 0.  The C++ standard says nothing
    // about what the return value is when __n == 0.
    pointer allocate(size_type n, const void* = 0)
    {
        if (__builtin_expect(n > this->max_size(), false))
            std::__throw_bad_alloc();

        Tp* ret = static_cast<Tp*>(numa::numaAllocate(n * sizeof(Tp), node_));
        if (!ret)
        {
            ISAAC_THREAD_CERR << "numaAllocate failed for " << n * sizeof(Tp) << " bytes on node " << node_ << " for type " << typeid(Tp).name() << std::endl;
            std::__throw_bad_alloc();
        }
        if (isNumaAvailable())
        {
            ISAAC_THREAD_CERR << "numaAllocate allocated " << n * sizeof(Tp) << " bytes on node " << node_ << " for type " << typeid(Tp).name() << std::endl;
        }
        return ret;
    }

    // __p is not permitted to be a null pointer.
    void deallocate(pointer p, size_type n)
    {
        numa::numaDeallocate(p, n * sizeof(Tp), node_);
    }

    size_type max_size() const throw() {return size_t(-1) / sizeof(Tp);}

    void construct(pointer p, const Tp& __val) { ::new((void *)p) Tp(__val); }

    void destroy(pointer p) { p->~Tp(); }

//    NumaAllocator &operator =(const NumaAllocator &that)
//    {
//        node_ = that.node_;
//        return *this;
//    }

    bool operator == (const NumaAllocator &that) const {return that.node_ == node_;}
    bool operator != (const NumaAllocator &that) const {return that.node_ != node_;}
//    template<typename _T> friend bool operator==(const NumaAllocator<_T, defaultNode>& left, const NumaAllocator<_T, defaultNode>& right);

    template<typename Tp1, int dN> friend class NumaAllocator;
    template<typename Tp1, int dN>
        friend std::ostream operator << (std::ostream &os, const NumaAllocator<Tp1, dN> &allocator);

//    typedef std::true_type propagate_on_container_copy_assignment;
};

template<int defaultNode>
class NumaAllocator<void, defaultNode>
{
    const int node_;

public:
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;

    template<typename Tp1>
    struct rebind { typedef NumaAllocator<Tp1, defaultNode> other; };

    NumaAllocator() throw() :node_(defaultNode) { }
    explicit NumaAllocator(const int node) throw() :node_(node) { }

    NumaAllocator(const NumaAllocator&) throw() { }

    template<typename Tp1>
    NumaAllocator(const NumaAllocator<Tp1, defaultNode>& that) throw(): node_(that.node_) { }

    ~NumaAllocator() throw() { }

//    NumaAllocator &operator =(const NumaAllocator &that)
//    {
//        node_ = that.node_;
//        return *this;
//    }

    bool operator == (const NumaAllocator &that) const {return that.node_ == node_;}
    bool operator != (const NumaAllocator &that) const {return that.node_ != node_;}
//    template<typename T> friend bool operator==(const NumaAllocator<T, defaultNode>& left, const NumaAllocator<T, defaultNode>& right);
    template<typename Tp, int dN> friend class NumaAllocator;
    template<typename Tp, int dN>
        friend std::ostream operator << (std::ostream &os, const NumaAllocator<Tp, dN> &allocator);

//    typedef std::true_type propagate_on_container_copy_assignment;
};

template<typename Tp, int defaultNode>
std::ostream &operator << (std::ostream &&os, const NumaAllocator<Tp, defaultNode> &allocator)
{
    os << "NumaAllocator<" << typeid(Tp).name() <<
        ">(" << allocator.node_  << ")" << std::endl;
    return os;
}

template<int defaultNode>
std::ostream &operator << (std::ostream &os, const NumaAllocator<void, defaultNode> &allocator)
{
    os << "NumaAllocator<" << typeid(void).name() <<
        ">(" << allocator.node_  << ")" << std::endl;
    return os;
}

//template<typename Tp, int defaultNode>
//bool operator==(const NumaAllocator<Tp, defaultNode>& left, const NumaAllocator<Tp, defaultNode>& right) { return left.node_ == right.node_; }
//
//template<typename Tp, int defaultNode>
//bool operator!=(const NumaAllocator<Tp, defaultNode>& left, const NumaAllocator<Tp, defaultNode>& right) { return !(left == right); }

} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_NUMA_ALLOCATOR_HH
