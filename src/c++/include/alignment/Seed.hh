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
 ** \file Seed.hh
 **
 ** \brief Definition od a seed by its k-mer and identifier (SeedId).
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_HH
#define iSAAC_ALIGNMENT_SEED_HH

#include <cassert>
#include <iostream>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>

#include "oligo/Kmer.hh"
#include "alignment/SeedId.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Structured unique identifier of a seed.
 **/
template <typename KmerT>
class Seed
{
public:
    Seed() : kmer_(0), seedId_(0) {}
    Seed(KmerT kmer, SeedId seedId) : kmer_(kmer), seedId_(seedId) {}
    Seed(const Seed &seed) : kmer_(seed.kmer_), seedId_(seed.seedId_) {}
    Seed &operator=(const Seed &seed)
    {
        if (this != &seed)
        {
            kmer_ = seed.kmer_;
            seedId_ = seed.seedId_;
        }
        return *this;
    }
    KmerT &kmer() {return kmer_;}
    KmerT getKmer() const {return kmer_;}
    const SeedId &getSeedId() const {return seedId_;}
//    uint64_t getSeedIndex() const {return seedId_.getSeed();}
    void invert(){kmer_ = oligo::reverseComplement(kmer_); seedId_.invert();}
    Seed inverted() const {Seed ret = *this; ret.invert(); return ret;}
    bool isReverse() const {return getSeedId().isReverse();}
    void setForward() {seedId_.setForward();}
    void setReverse() {seedId_.setReverse();}
    void setKmer(KmerT kmer) {kmer_ = kmer;}
    void setSeedId(SeedId seedId) {seedId_ = seedId;}
private:
    KmerT kmer_;
    SeedId seedId_;
};

template <typename KmerT>
inline Seed<KmerT> makeNSeed(uint64_t tile, uint64_t barcode, uint64_t cluster, bool lowestSeedId)
{
    return Seed<KmerT>(~KmerT(0), SeedId(tile, barcode, cluster, SeedId::SEED_LENGTH_MASK, !lowestSeedId));
}

template <typename KmerT>
inline bool orderByKmerSeedIndex(const Seed<KmerT> &lhs, const Seed<KmerT> &rhs)
{
    // IMPORTANT!!!: The match finder relies on N-seeds to be at the end of the seed list after sorting by kmer.
    // N-seeds are assigned the highest possible seed index value by seed loader
    return (lhs.getKmer() < rhs.getKmer()) || (lhs.getKmer() == rhs.getKmer() && lhs.getSeedId().getSeed() < rhs.getSeedId().getSeed());
}

template <typename KmerT>
inline std::ostream &operator<<(std::ostream &os, const Seed<KmerT> &seed)
{
    return os << "Seed(" << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(seed.getKmer(), oligo::KmerTraits<KmerT>::KMER_BASES) <<
        "(" << oligo::ReverseBases<oligo::BITS_PER_BASE, KmerT>(seed.getKmer(), oligo::KmerTraits<KmerT>::KMER_BASES) << ")" <<
        "," << seed.getSeedId() << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_HH
