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
 ** \file ReferenceHasher.cpp
 **
 ** Produces the mapping between kmer values and genomic positions
 **
 ** \author Roman Petrovski
 **/

#include <boost/assert.hpp>
#include <boost/foreach.hpp>

#include "common/Exceptions.hh"
#include "common/Numa.hh"
#include "common/ParallelSort.hpp"
#include "common/SystemCompatibility.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferencePosition.hh"
#include "reference/ReferenceHash.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

template <typename ReferenceHashT>
ReferenceHasher<ReferenceHashT>::ReferenceHasher (
    const SortedReferenceMetadata &sortedReferenceMetadata,
    const ContigList &contigList,
    common::ThreadVector &threads,
    const unsigned threadsMax)
    : BaseT(0, 0, contigList, sortedReferenceMetadata.getContigs())
    , sortedReferenceMetadata_(sortedReferenceMetadata)
    , threads_(threads)
    , threadsMax_(threadsMax)
    , mutexes_(threadsMax_ * 32) // reduce collision probability somewhat
    , threadBuffers_(threadsMax_, ThreadBuffer(mutexes_.capacity()))

{
    while (mutexes_.capacity() != mutexes_.size())
    {
        mutexes_.push_back(new boost::mutex);
    }

    for (ThreadBuffer &threadBuffer : threadBuffers_)
    {
        for (MutexBuffer &mutexBuffer : threadBuffer)
        {
            mutexBuffer.reserve(THREAD_BUFFER_KMERS_MAX);
        }
    }
}

template<typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::dumpCounts(
    boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash)
{
    boost::lock_guard<boost::mutex> lock(mutex);
    for (const KmerWithPosition& kwp : buffer)
    {
        ++referenceHash.offsets_[kwp.first.bits_];
        //    ISAAC_THREAD_CERR << pos << std::endl;
    }
    buffer.clear();
}


template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::updateOffsets(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const KmerT &kmer)
{
//    ISAAC_THREAD_CERR << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(kmer, oligo::KmerTraits<KmerT>::KMER_BASES) << " at " << kmerPosition << std::endl;
    const std::size_t mutexIndex = (kmer.bits_ / THREAD_BUFFER_KMERS_MAX) % mutexes_.size();
    MutexBuffer &buffer = threadBuffers_[threadNumber][mutexIndex];
    buffer.push_back(std::make_pair(kmer, ReferencePosition()));

    if (buffer.capacity() == buffer.size())
    {
        dumpCounts(mutexes_.at(mutexIndex), buffer, referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::countKmers(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads)
{
    BaseT::kmerThread(
        oligo::Permutate(oligo::KmerTraits<KmerT>::KMER_BASES), threadNumber, threads,
        boost::bind(&ReferenceHasher::updateOffsets, this, boost::ref(referenceHash), threadNumber, _2));

    ThreadBuffer &threadBuffer = threadBuffers_[threadNumber];
    for (std::size_t mutexIndex = 0; threadBuffer.size() > mutexIndex; ++mutexIndex)
    {
        dumpCounts(mutexes_.at(mutexIndex), threadBuffer[mutexIndex], referenceHash);
    }
}

template<typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::dumpPositions(
    boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash)
{
    boost::lock_guard<boost::mutex> lock(mutex);
    for (const KmerWithPosition& kwp : buffer)
    {
        const std::size_t offset = referenceHash.offsets_[kwp.first.bits_]++;
        referenceHash.positions_.at(offset) = kwp.second;
        //    ISAAC_THREAD_CERR << pos << std::endl;
    }
    buffer.clear();
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::storePosition(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const KmerT &kmer,
    const int contigId,
    const uint64_t kmerPosition,
    const bool reverse,
    const std::vector<unsigned> &karyotypes,
    const std::vector<uint64_t> &contigOffsets)
{
    const std::size_t mutexIndex = (kmer.bits_ / THREAD_BUFFER_KMERS_MAX) % mutexes_.size();
    MutexBuffer &buffer = threadBuffers_[threadNumber][mutexIndex];
    buffer.push_back(
        std::make_pair(kmer, ReferencePosition(contigId, kmerPosition, false, reverse).translateContig(karyotypes)));

    if (buffer.capacity() == buffer.size())
    {
        dumpPositions(mutexes_.at(mutexIndex), buffer, referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::storePositions(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads,
    const std::vector<unsigned> &karyotypes,
    const std::vector<uint64_t> &contigOffsets)
{
    BaseT::kmerThread(
        oligo::Permutate(oligo::KmerTraits<KmerT>::KMER_BASES), threadNumber, threads,
        boost::bind(&ReferenceHasher::storePosition, this,
                    boost::ref(referenceHash), _1, _2, _3, _4, _5, boost::ref(karyotypes),
                    boost::ref(contigOffsets)));
    ThreadBuffer &threadBuffer = threadBuffers_[threadNumber];
    for (std::size_t mutexIndex = 0; threadBuffer.size() > mutexIndex; ++mutexIndex)
    {
        dumpPositions(mutexes_.at(mutexIndex), threadBuffer[mutexIndex], referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::sortPositions(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads)
{
    const std::size_t blockLength = (referenceHash.offsets_.size() - 1) / threads + 1;
    const std::size_t blockBegin = blockLength * threadNumber;
    if (blockBegin < referenceHash.offsets_.size())
    {
        const std::size_t blockEnd = std::min(referenceHash.offsets_.size(), blockBegin + blockLength + 1);
        for (typename Offsets::iterator it = referenceHash.offsets_.begin() + blockBegin + 1;
            referenceHash.offsets_.begin() + blockEnd != it; ++it)
        {
            std::sort(referenceHash.positions_.begin() + *(it - 1), referenceHash.positions_.begin() + *it);
        }
    }
}

/**
 * \brief replaces 0 offsets with the offset of the previous k-mer
 */
template <typename KmerT>
void ReferenceHasher<KmerT>::updateEmptyOffsets(Offsets& offsets)
{
    Offset lastOffset = 0;
    BOOST_FOREACH(Offset &offset, offsets)
    {
        if (!offset)
        {
            offset = lastOffset;
        }
        else
        {
            lastOffset = offset;
        }
    }
}

/**
 * \brief replaces 0 offsets with the offset of the previous k-mer
 */
template <typename ReferenceHashT>
typename ReferenceHasher<ReferenceHashT>::Offset
ReferenceHasher<ReferenceHashT>::countsToOffsets(Offsets& offsets)
{
    Offset offset = 0;
    BOOST_FOREACH(Offset &count, offsets)
    {
        if (count)
        {
            using std::swap; swap(offset, count);
            offset += count;
        }
    }
    return offset;
}


template <typename ReferenceHashT>
ReferenceHashT ReferenceHasher<ReferenceHashT>::generate()
{
    ReferenceHashT ret;

    generate(ret);

    return ret;
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::generate(ReferenceHashT &ret)
{
    ISAAC_THREAD_CERR <<
        "Constructing ReferenceHasher: for " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers " << std::endl;

    const std::vector<uint64_t> contigOffsets(computeContigOffsets(sortedReferenceMetadata_.getKaryotypeOrderedContigs()));

    std::vector<unsigned> karyotypes;
    // translation from karyotyope order to the original contig index is needed as this is the order in which
    // MatchSelector expects to receive data.
    karyotypes.reserve(sortedReferenceMetadata_.getContigs().size());
    std::transform(sortedReferenceMetadata_.getContigs().begin(), sortedReferenceMetadata_.getContigs().end(), std::back_inserter(karyotypes),
                   boost::bind(&reference::SortedReferenceMetadata::Contig::karyotypeIndex_, _1));

    threads_.execute(boost::bind(&ReferenceHasher::countKmers, this, boost::ref(ret), _1, _2), threadsMax_);

    const Offset total = countsToOffsets(ret.offsets_);

    ISAAC_THREAD_CERR << " found " << total << " positions on " << threadsMax_ << " threads" << std::endl;
    ret.positions_.resize(total);
    ISAAC_THREAD_CERR << " reserving memory done for " << ret.positions_.size() << " positions" << std::endl;

    threads_.execute(
        [this, &ret, &karyotypes, &contigOffsets](const unsigned threadNumber, const std::size_t threads)
        {
            ReferenceHasher::storePositions(ret, threadNumber, threads, karyotypes, contigOffsets);
        }, threadsMax_);

    ISAAC_THREAD_CERR << " generated " << total << " positions" << std::endl;

    updateEmptyOffsets(ret.offsets_);

    threads_.execute(
        [this, &ret](const unsigned threadNumber, const std::size_t threads)
        {
            sortPositions(ret, threadNumber, threads);
        }, threadsMax_);

    ISAAC_THREAD_CERR << " sorted " << ret.offsets_.back() << " positions" << std::endl;
}
//
template class ReferenceHasher<ReferenceHash<oligo::VeryShortKmerType> >;
//template class ReferenceHasher<oligo::BasicKmerType<12> >;
////template class ReferenceHasher<oligo::BasicKmerType<14> >;
template class ReferenceHasher<ReferenceHash<oligo::ShortKmerType, common::NumaAllocator<void, 0> > >;
template class ReferenceHasher<ReferenceHash<oligo::ShortKmerType> >;
////template class ReferenceHasher<oligo::BasicKmerType<18> >;
//template class ReferenceHasher<oligo::BasicKmerType<20> >;
////template class ReferenceHasher<oligo::BasicKmerType<22> >;

} // namespace reference
} // namespace isaac
