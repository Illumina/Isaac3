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
 ** \file Contig.hh
 **
 ** \brief Definition of a contig
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_CONTIG_HH
#define iSAAC_REFERENCE_CONTIG_HH

#include <iostream>
#include <string>
#include <vector>

#include "common/Debug.hh"
#include "common/NumaContainer.hh"
#include "common/SameAllocatorVector.hh"

namespace isaac
{
namespace reference
{

template <typename AllocatorT>
struct BasicContig : public std::vector<char, AllocatorT>
{
    unsigned index_;
    std::string name_;
    typedef std::vector<char, AllocatorT> BaseT;

    BasicContig(const AllocatorT &allocator) : BaseT(allocator),
        index_(0), name_()
    {
//        ISAAC_THREAD_CERR << "BasicContig(const AllocatorT &allocator)" << std::endl;
    }

    ~BasicContig()
    {
//        ISAAC_THREAD_CERR << "~BasicContig(" << BaseT::size() << ")" << std::endl;
    }

    BasicContig() :
        index_(0), name_()
    {
//        ISAAC_THREAD_CERR << "BasicContig()" << std::endl;
    }

    BasicContig(const BasicContig &that) : BaseT(that),
        index_(that.index_), name_(that.name_)
    {
//        ISAAC_THREAD_CERR << "BasicContig(const BasicContig &that)" << std::endl;
    }

    BasicContig(const BasicContig &that, const AllocatorT &allocator) : BaseT(that, allocator),
        index_(that.index_), name_(that.name_)
    {
//        ISAAC_THREAD_CERR << "BasicContig(const BasicContig &that, const AllocatorT &allocator)" << std::endl;
        ;
    }

    BasicContig(BasicContig &&that) : BaseT(std::move(that)),
        index_(that.index_), name_(that.name_)
    {
//        ISAAC_THREAD_CERR << "BasicContig(BasicContig &&that)" << std::endl;
    }


    BasicContig(const unsigned index, const std::string &name) : index_(index), name_(name){;}
    template <typename ContainerT>
    BasicContig(const unsigned index, const std::string &name, const ContainerT &init) : index_(index), name_(name)
    {
        BaseT::assign(init.begin(), init.end());
    }
    size_t getLength() const {return BaseT::size();}

    BasicContig & operator=(const BasicContig &that)
    {
        index_ = that.index_;
        name_ = that.name_;
        BaseT::operator = (that);
//        ISAAC_THREAD_CERR << "operator=(const BasicContig &that)" << std::endl;
        return *this;
    }

    BasicContig & operator=(BasicContig &&that)
    {
        index_ = that.index_;
        name_ = that.name_;
        BaseT::swap(that);
//        ISAAC_THREAD_CERR << "operator=(BasicContig &&that)" << std::endl;
        return *this;
    }

    friend std::ostream &operator <<(std::ostream &os, const BasicContig<AllocatorT> &contig)
    {
        return os << "Contig(" << contig.index_ << "," << contig.name_ << "," << contig.size() << ")";
    }
};

typedef BasicContig<common::NumaAllocator<char, 0> > Contig;
typedef common::SameAllocatorVector<Contig, common::NumaAllocator<Contig, 0> > ContigList;
typedef common::SameAllocatorVector<ContigList, common::NumaAllocator<ContigList, 0> > ContigLists;

class NumaContigLists
{
    common::NumaContainerReplicas<ContigLists> replicas_;
public:

    NumaContigLists(ContigLists &&node0Lists) :replicas_(std::move(node0Lists))
    {
        ISAAC_THREAD_CERR << "NumaContigLists ContigLists constructor" << std::endl;
    }

    const ContigLists &node0Container() const {return replicas_.node0Container();}
    const ContigLists &threadNodeContainer() const {return replicas_.threadNodeContainer();}
//    operator const ContigLists &()const {return replicas_.threadNodeContainer();}
};

/// Total length of all the contigs of a genome
size_t genomeLength(const ContigList &contigList);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIG_HH
