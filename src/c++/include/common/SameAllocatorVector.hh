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
 ** \file FiniteCapacityVector.hh
 **
 ** Vector that passes its allocator to contained objects.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_COMMON_SAME_ALLOCATOR_VECTOR_HH
#define iSAAC_COMMON_SAME_ALLOCATOR_VECTOR_HH

#include <initializer_list>
#include <boost/foreach.hpp>

namespace isaac
{
namespace common
{

template <typename ValueT, typename AllocatorT>
class SameAllocatorVector : std::vector<ValueT, typename AllocatorT::template rebind<ValueT>::other>
{
public:
    typedef std::vector<ValueT, typename AllocatorT::template rebind<ValueT>::other> BaseT;

public:
    using typename BaseT::value_type;
    typedef typename BaseT::allocator_type allocator_type;

    explicit SameAllocatorVector(const std::size_t s)
        : BaseT(s)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const std::size_t s)" << std::endl;
    }

    SameAllocatorVector(const std::size_t s, const value_type &v)
        : BaseT(s, v)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const std::size_t s)" << std::endl;
    }

    ~SameAllocatorVector()
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">~SameAllocatorVector()" << std::endl;
    }


    SameAllocatorVector()
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">()" << std::endl;
    }

    SameAllocatorVector(const SameAllocatorVector &that)
        : BaseT(that)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const SameAllocatorVector &that)" << std::endl;
    }

    SameAllocatorVector(const SameAllocatorVector &that, const AllocatorT &allocator)
        : BaseT(allocator)
    {
        BaseT::resize(that.size());
        int i = 0;
        BOOST_FOREACH(const value_type & v, that)
        {
            (*this)[i++] = value_type(v, typename BaseT::value_type::allocator_type(allocator));
        }
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const SameAllocatorVector &that, const AllocatorT &allocator)" << std::endl;
    }

    SameAllocatorVector(SameAllocatorVector &&that)
    {
        swap(that);
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(SameAllocatorVector &&that)" << std::endl;
    }

    template<typename InputIterator>
    SameAllocatorVector(InputIterator first, InputIterator last)
    : BaseT(first, last)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(InputIterator first, InputIterator last)" << std::endl;
    }

    template<typename InputIterator>
    SameAllocatorVector(InputIterator first, InputIterator last, const AllocatorT& a)
    : BaseT(first, last, a)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">(InputIterator first, InputIterator last, const AllocatorT& a))" << std::endl;
    }

    SameAllocatorVector & operator=(SameAllocatorVector &&that)
    {
        swap(that);
//        ISAAC_THREAD_CERR << "SameAllocatorVector<" <<  demangle(typeid(ValueT).name()) <<
//            ">::operator=(SameAllocatorVector &&that)" << std::endl;
        return *this;
    }

    void swap(SameAllocatorVector &that)
    {
        BaseT::swap(that);
    }

    typedef typename BaseT::iterator iterator;
    typedef typename BaseT::const_iterator const_iterator;


    using BaseT::operator [];
    using BaseT::at;
    using BaseT::push_back;
    using BaseT::size;
    using BaseT::begin;
    using BaseT::end;
    using BaseT::front;
    using BaseT::reserve;
    using BaseT::resize;
    using BaseT::empty;
    using BaseT::clear;
};


template <typename ValueT, typename AllocatorT>
class SameAllocatorVectorEnd : std::vector<ValueT, typename AllocatorT::template rebind<ValueT>::other>
{
public:
    typedef std::vector<ValueT, typename AllocatorT::template rebind<ValueT>::other> BaseT;

public:
    using typename BaseT::value_type;
    typedef typename BaseT::allocator_type allocator_type;

    explicit SameAllocatorVectorEnd(const std::size_t s)
        : BaseT(s)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const std::size_t s)" << std::endl;
    }

    SameAllocatorVectorEnd(const std::size_t s, const value_type v)
        : BaseT(s, v)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const std::size_t s, const value_type v)" << std::endl;
    }

    ~SameAllocatorVectorEnd()
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">~SameAllocatorVectorEnd()" << std::endl;
    }


    SameAllocatorVectorEnd()
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">()" << std::endl;
    }

    SameAllocatorVectorEnd(const SameAllocatorVectorEnd &that)
        : BaseT(that)
    {
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const SameAllocatorVectorEnd &that)" << std::endl;
    }

    SameAllocatorVectorEnd(const SameAllocatorVectorEnd &that, const AllocatorT &allocator)
        : BaseT(allocator)
    {
//        static_cast<BaseT&>(*this) = that;
        BaseT::resize(that.size());
        int i = 0;
        BOOST_FOREACH(const value_type & v, that)
        {
            (*this)[i++] = v;//value_type(v, typename BaseT::value_type::allocator_type(allocator));
        }
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">(const SameAllocatorVectorEnd &that, const AllocatorT &allocator)" << std::endl;
    }

    SameAllocatorVectorEnd(SameAllocatorVectorEnd &&that)
    {
        swap(that);
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">(SameAllocatorVectorEnd &&that)" << std::endl;
    }

    SameAllocatorVectorEnd & operator=(SameAllocatorVectorEnd &&that)
    {
        swap(that);
//        ISAAC_THREAD_CERR << "SameAllocatorVectorEnd<" <<  demangle(typeid(ValueT).name()) <<
//            ">::operator=(SameAllocatorVectorEnd &&that)" << std::endl;
        return *this;
    }

    void swap(SameAllocatorVectorEnd &that)
    {
        BaseT::swap(that);
    }

    typedef typename BaseT::iterator iterator;
    typedef typename BaseT::const_iterator const_iterator;


    using BaseT::operator [];
    using BaseT::at;
    using BaseT::push_back;
    using BaseT::size;
    using BaseT::begin;
    using BaseT::end;
    using BaseT::front;
    using BaseT::reserve;
    using BaseT::resize;
    using BaseT::empty;
};


} //namespace common
} //namespace isaac

#endif // #ifndef iSAAC_COMMON_SAME_ALLOCATOR_VECTOR_HH
