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
 ** \file ReferenceSorter.hh
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_HASH_HH
#define iSAAC_REFERENCE_REFERENCE_HASH_HH

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>

#include "common/NumaContainer.hh"
#include "oligo/Kmer.hh"
#include "PermutatedKmerGenerator.hh"
#include "reference/ReferenceKmer.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{
template <typename KmerT> class ReferenceHasher;

template <typename KmerType, typename AllocatorT = std::allocator<void> >
class ReferenceHash
{
    typedef ReferenceHash<KmerType, AllocatorT> MyT;
    typedef typename AllocatorT::template rebind<reference::ReferencePosition> ReferencePositionAllocatorRebind;
    typedef typename ReferencePositionAllocatorRebind::other ReferencePositionAllocator;
    typedef std::vector<reference::ReferencePosition, ReferencePositionAllocator> Positions;

    typedef boost::uint32_t Offset;
    typedef typename AllocatorT::template rebind<Offset> OffsetAllocatorRebind;
    typedef typename OffsetAllocatorRebind::other OffsetAllocator;
    typedef std::vector<Offset, OffsetAllocator > Offsets;
public:
    typedef KmerType KmerT;
    typedef typename Positions::const_iterator const_iterator;
    typedef std::pair<const_iterator, const_iterator> MatchRange;
    typedef void value_type;// compatibility with std containers for numa replications


    ReferenceHash(const AllocatorT &allocator)
        : offsets_(1UL << oligo::KmerTraits<KmerT>::KMER_BITS, 0, OffsetAllocator(allocator))
        , positions_(ReferencePositionAllocator(allocator))
    {
    }

    ReferenceHash()
        : offsets_(1UL << oligo::KmerTraits<KmerT>::KMER_BITS, 0, AllocatorT())
        , positions_(AllocatorT())
    {
    }

    ReferenceHash(const ReferenceHash &that)
        : offsets_(that.offsets_)
        , positions_(that.positions_)
    {
    }

    ReferenceHash(const ReferenceHash &that, const AllocatorT &allocator)
        : offsets_(OffsetAllocator(allocator))
        , positions_(ReferencePositionAllocator(allocator))
    {
        offsets_ = that.offsets_;
        positions_ = that.positions_;
    }

    ReferenceHash(ReferenceHash &&that)
        : offsets_(std::move(that.offsets_))
        , positions_(std::move(that.positions_))
    {
    }

    ReferenceHash & operator=(const ReferenceHash &that)
    {
        offsets_ = that.offsets_;
        positions_ = that.positions_;
        return *this;
    }

    /**
     * \brief copies data over using the currently set allocator
     */
    template <typename OtherHashT>
    void assign(const OtherHashT& that)
    {
        offsets_.assign(that.offsets_.begin(), that.offsets_.end());
        positions_.assign(that.positions_.begin(), that.positions_.end());
    }

    MatchRange findMatches(const KmerT &kmer) const
    {
        boost::uint32_t positionsBegin = !kmer ? 0 : offsets_[kmer.bits_ - 1];
        boost::uint32_t positionsEnd = offsets_[kmer.bits_];
        ISAAC_ASSERT_MSG(positionsBegin <= positions_.size(), "Positions buffer overrun by positionsBegin:" << positionsBegin << " for kmer " << kmer);
        ISAAC_ASSERT_MSG(positionsBegin <= positionsEnd, "positionsEnd:" << positionsEnd << " overrun by positionsBegin:" << positionsBegin << " for kmer " << kmer);

        const MatchRange ret = std::make_pair(positions_.begin() + positionsBegin, positions_.begin() + positionsEnd);

    //    ISAAC_THREAD_CERR << "found " << std::distance(ret.first, ret.second) << " matches for " << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(kmer, oligo::KmerTraits<KmerT>::KMER_BASES) << std::endl;
    //    BOOST_FOREACH(const ReferencePosition &pos, ret)
    //    {
    //        ISAAC_THREAD_CERR << pos << std::endl;
    //    }
        return ret;
    }

private:
    Offsets offsets_;
    Positions positions_;

    friend class ReferenceHasher<MyT>;
};

template <typename ReferenceHashT>
class ReferenceHasher: PermutatedKmerGenerator<typename ReferenceHashT::KmerT, permutatedKmerGenerator::ForwardOnly>
{
    typedef typename ReferenceHashT::KmerT KmerT;
    typedef PermutatedKmerGenerator<KmerT, permutatedKmerGenerator::ForwardOnly> BaseT;
    typedef typename ReferenceHashT::Positions Positions;
    typedef typename ReferenceHashT::Offset Offset;
    typedef typename ReferenceHashT::Offsets Offsets;
    static const std::size_t THREAD_BUFFER_KMERS_MAX = 8192; // arbitrary number that reduces the cost/benefit of acquiring a mutex
public:

    ReferenceHasher(
        const SortedReferenceMetadata &sortedReferenceMetadata,
        const ContigList &contigList,
        common::ThreadVector &threads,
        const unsigned threadsMax);

    ReferenceHashT generate();
    void generate(ReferenceHashT &ret);

private:
    const reference::SortedReferenceMetadata &sortedReferenceMetadata_;
    common::ThreadVector &threads_;
    const unsigned threadsMax_;

    boost::ptr_vector<boost::mutex> mutexes_;
    typedef std::pair<typename ReferenceHashT::KmerT, reference::ReferencePosition> KmerWithPosition;
    typedef std::vector<KmerWithPosition> MutexBuffer;
    typedef std::vector<MutexBuffer> ThreadBuffer;
    typedef std::vector<ThreadBuffer> ThreadBuffers;
    ThreadBuffers threadBuffers_;

    void updateOffsets(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const KmerT &kmer);

    void countKmers(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const std::size_t threads);

    void storePosition(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const KmerT &kmer,
        const int contigId,
        const uint64_t kmerPosition,
        const bool reverse,
        const std::vector<unsigned> &karyotypes,
        const std::vector<uint64_t> &contigOffsets);

    void storePositions(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const std::size_t threads,
        const std::vector<unsigned> &karyotypes,
        const std::vector<uint64_t> &contigOffsets);

    static void sortPositions(ReferenceHashT &referenceHash, const unsigned threadNumber, const std::size_t threads);
    static void updateEmptyOffsets(Offsets& offsets);
    static Offset countsToOffsets(Offsets& offsets);
    static void dumpCounts(boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash);
    static void dumpPositions(boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash);
};

template <typename HashType>
class NumaReferenceHash
{
    common::NumaContainerReplicas<HashType> replicas_;
public:
    typedef typename HashType::KmerT KmerT;
    typedef typename HashType::MatchRange MatchRange;
    typedef typename HashType::const_iterator const_iterator;

    NumaReferenceHash(HashType &&hash) :replicas_(std::move(hash))
    {
        ISAAC_THREAD_CERR << "NumaReferenceHash copy constructor" << std::endl;
    }

    MatchRange findMatches(const KmerT &kmer) const
    {
        return replicas_.threadNodeContainer().findMatches(kmer);
    }
};
} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_HASH_HH
