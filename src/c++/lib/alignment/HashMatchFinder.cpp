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
 ** \file FindHashMatchesTransition.cpp
 **
 ** \brief see FindHashMatchesTransition.hh
 **
 ** \author Roman Petrovski
 **/

#include "flowcell/Layout.hh"
#include "alignment/HashMatchFinder.hh"

namespace isaac
{
namespace alignment
{


template <typename ReferenceHash>
SeedHashMatchFinder<ReferenceHash>::SeedHashMatchFinder(
    const ReferenceHash& referenceHash,
    const unsigned seedBaseQualityMin,
    const unsigned repeatThreshold)
    : referenceHash_(referenceHash)
    , seedBaseQualityMin_(seedBaseQualityMin)
    , repeatThreshold_(repeatThreshold)
    , noExtendRepeatThreshold_(sqrt(repeatThreshold_) + 1/*FAKE_NOT_ENOUGH_REPEATS*/)
{
}

template <typename ReferenceHash>
bool SeedHashMatchFinder<ReferenceHash>::findSeedMatches(
    const Cluster& cluster,
    const alignment::SeedMetadata &seedMetadata,
    const flowcell::ReadMetadata& readMetadata,
    Matches& matches,
    std::size_t& repeatSeeds) const
{
    alignment::Seed<KmerT> seed(KmerT(0), alignment::SeedId(seedMetadata.getLength(), false));
    if (updateSeedKmer(cluster, seedMetadata, seedBaseQualityMin_, seed))
    {
        //                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(seed.getCluster(), seed);
        //                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(seed.getCluster(), "findClusterMatches " << seed << ":" << seedMetadata);
        return extendSeed(cluster, seed, seedMetadata, readMetadata, referenceHash_, matches, repeatSeeds);
    }
    return false;
}

template <typename KmerT>
bool updateSeedKmer(
    const Cluster& cluster,
    const unsigned readIndex,
    const unsigned short offset,
    const unsigned short length,
    const unsigned seedBaseQualityMin,
    alignment::Seed<KmerT> &seed)
{
    const alignment::BclClusters::const_iterator cyclesBegin = cluster.getBclData(readIndex) + offset;
    const alignment::BclClusters::const_iterator cyclesEnd = cyclesBegin + length;
    for (alignment::BclClusters::const_iterator cycle = cyclesBegin; cyclesEnd != cycle; ++cycle)
    {
        const unsigned char base = *cycle;
        if (seedBaseQualityMin > oligo::getQuality(base))
        {
            return false;
        }

        const KmerT forwardBaseValue(base & oligo::BITS_PER_BASE_MASK);

        seed.kmer() <<= oligo::BITS_PER_BASE;
        seed.kmer() |= forwardBaseValue;
    }

    if (seed.isReverse())
    {
        seed.kmer() = oligo::reverseComplement(seed.kmer());
    }

    return true;
}

template <typename KmerT>
bool updateSeedKmer(
    const Cluster& cluster,
    const alignment::SeedMetadata &seedMetadata,
    const unsigned seedBaseQualityMin,
    alignment::Seed<KmerT> &seed)
{
    return updateSeedKmer(
        cluster,
        seedMetadata.getReadIndex(), seedMetadata.getOffset(), seedMetadata.getLength(), seedBaseQualityMin,
        seed);
}

template <typename ReferenceHash>
bool SeedHashMatchFinder<ReferenceHash>::storeFwExtensionMatches(
    const typename ReferenceHash::MatchRange& oriPositions,
    const typename ReferenceHash::MatchRange& extPositions,
    const uint64_t extSeedOffset,
    const Seed<KmerT> &oriSeed,
    const Seed<KmerT> &extSeed,
    const unsigned matchesMax,
    Matches& matches) const
{
    bool neighborlessMatchFound = false;
    // store
    for (typename ReferenceHash::const_iterator ori = oriPositions.first, alt = extPositions.first;
        oriPositions.second != ori && extPositions.second != alt && matches.size() < matchesMax;)
    {
        const reference::ReferencePosition correctedPos = *ori + extSeedOffset;
        // Note: simple comparison of ReferencePosition will account of reverse flag, which not what want
        if (correctedPos.getLocation() < alt->getLocation())
        {
            ++ori;
        }
        else if (correctedPos.getLocation() > alt->getLocation())
        {
            ++alt;
        }
        else
        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oriSeed.getCluster(), "fw writeExtMatch: " << oriSeed << " " << *ori);
            reference::ReferencePosition pos(*ori);
            neighborlessMatchFound |= !pos.hasNeighbors();
            matches.push_back(Match(SeedId(oriSeed.getSeedId().getSeedLength() * 2, oriSeed.getSeedId().getReverse()),
                                    pos));

            ++ori;
            ++alt;
        }
    }

    return neighborlessMatchFound;
}

template <typename ReferenceHash>
bool SeedHashMatchFinder<ReferenceHash>::storeRvExtensionMatches(
    const typename ReferenceHash::MatchRange& oriPositions,
    const typename ReferenceHash::MatchRange& extPositions,
    const uint64_t extSeedOffset,
    const Seed<KmerT> &oriSeed,
    const Seed<KmerT> &extSeed,
    const unsigned matchesMax,
    Matches& matches) const
{
    bool neighborlessMatchFound = false;
    // store
    for (typename ReferenceHash::const_iterator ori = oriPositions.first, alt = extPositions.first;
        oriPositions.second != ori && extPositions.second != alt && matches.size() < matchesMax;)
    {
        // skip position that will get extended before the start of the contig
        if (ori->getPosition() < extSeedOffset)
        {
            ++ori;
            continue;
        }

        const reference::ReferencePosition correctedPos = *ori - extSeedOffset;
        // Note: simple comparison of ReferencePosition will account of reverse flag, which not what want
        if (correctedPos.getLocation() < alt->getLocation())
        {
            ++ori;
        }
        else if (correctedPos.getLocation() > alt->getLocation())
        {
            ++alt;
        }
        else
        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oriSeed.getCluster(), "rv writeExtMatch: " << oriSeed << " " << *alt);
            reference::ReferencePosition pos(*alt);
            neighborlessMatchFound |= !pos.hasNeighbors();
            matches.push_back(Match(SeedId(oriSeed.getSeedId().getSeedLength() * 2, oriSeed.getSeedId().getReverse()),
                                    pos));
            ++ori;
            ++alt;
        }
    }
    return neighborlessMatchFound;
}

template <typename ReferenceHash>
bool SeedHashMatchFinder<ReferenceHash>::extendSeed(
    const Cluster& cluster,
    alignment::Seed<KmerT> fwSeed,
    const alignment::SeedMetadata &seedMetadata,
    const flowcell::ReadMetadata &readMetadata,
    const ReferenceHash &referenceHash,
    Matches& matches,
    std::size_t &repeatSeeds) const
{
    // TODO: either compute these constants or have them set from command line
    const std::size_t FAKE_TOO_MANY_REPEATS = repeatThreshold_ * 5000;

    typename ReferenceHash::MatchRange fwMatchPositions = referenceHash.findMatches(fwSeed.getKmer());
    const Seed<KmerT> rvSeed = fwSeed.inverted();
    typename ReferenceHash::MatchRange rvMatchPositions = referenceHash.findMatches(rvSeed.getKmer());

    const std::size_t seedRepeatCount = std::size_t(std::distance(fwMatchPositions.first, fwMatchPositions.second) +
                                                    std::distance(rvMatchPositions.first, rvMatchPositions.second));
    if (seedRepeatCount >= FAKE_TOO_MANY_REPEATS)
    {
        ++repeatSeeds;
        return false;
    }
//    else if (seedRepeatCount < noExtendRepeatThreshold_)
//    {
//        bool neighborlessMatchFound = false;
//        BOOST_FOREACH(reference::ReferencePosition pos, fwMatchPositions)
//        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "fw writeMatch: " << fwSeed << " " << pos);
//            neighborlessMatchFound |= !pos.hasNeighbors();
//            matches.push_back(Match(fwSeed.getSeedId(), pos));
//        }
//        BOOST_FOREACH(reference::ReferencePosition pos, rvMatchPositions)
//        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "rv writeMatch: " << rvSeed << " " << pos);
//            neighborlessMatchFound |= !pos.hasNeighbors();
//            matches.push_back(Match(rvSeed.getSeedId(), pos));
//        }
//        return neighborlessMatchFound;
//    }


    const std::size_t extSeedOffset = seedMetadata.getLength();

    if (seedMetadata.getOffset() + seedMetadata.getLength() + extSeedOffset > readMetadata.getLength())
    {
        return false;
    }

    Seed<KmerT> fwExtSeed = fwSeed;
    if (!updateSeedKmer(cluster, seedMetadata.getReadIndex(), seedMetadata.getOffset() + extSeedOffset, seedMetadata.getLength(), seedBaseQualityMin_, fwExtSeed))
    {
        return false;
    }

    const Seed<KmerT> rvExtSeed = fwExtSeed.inverted();

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), fwSeed << ":" << fwExtSeed);
    typename ReferenceHash::MatchRange fwExtPositions = referenceHash.findMatches(fwExtSeed.getKmer());
    typename ReferenceHash::MatchRange rvExtPositions = referenceHash.findMatches(rvExtSeed.getKmer());

    if (std::size_t(std::distance(fwExtPositions.first, fwExtPositions.second) + std::distance(rvExtPositions.first, rvExtPositions.second)) >= FAKE_TOO_MANY_REPEATS)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "too many repeats ori:" << std::distance(fwMatchPositions.first, fwMatchPositions.second));
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "too many repeats ext:" << std::distance(fwExtPositions.first, fwExtPositions.second));
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "too many repeats ori:" << std::distance(rvMatchPositions.first, rvMatchPositions.second));
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "too many repeats ext:" << std::distance(rvExtPositions.first, rvExtPositions.second));
        // too many candidates. deal with this when annotation is available.
        ++repeatSeeds;
        return false;
    }

    const std::size_t before = matches.size();
    bool neighborlessMatchFound = storeFwExtensionMatches(fwMatchPositions, fwExtPositions, extSeedOffset, fwSeed, fwExtSeed, before + repeatThreshold_, matches);
    if (matches.size() - before < repeatThreshold_)
    {
        neighborlessMatchFound |= storeRvExtensionMatches(rvMatchPositions, rvExtPositions, extSeedOffset, rvSeed, rvExtSeed, before + repeatThreshold_, matches);
    }

    if (matches.size() - before < repeatThreshold_)
    {
        // as this is an extended seed, all alignments are based on same kmer. So, either all of them have neigbors or neither of them does.
        return neighborlessMatchFound;
    }
    else
    {
        ++repeatSeeds;
        matches.resize(before);
        return false;
    }
}


template <typename ReferenceHash>
ClusterHashMatchFinder<ReferenceHash>::ClusterHashMatchFinder(
    const ReferenceHash& referenceHash,
    const unsigned seedBaseQualityMin,
    const unsigned repeatThreshold,
    const alignment::SeedMetadataList &seedMetadataList) :
    SeedHashMatchFinder<ReferenceHash>(referenceHash, seedBaseQualityMin, repeatThreshold)
    , seedMetadataList_(seedMetadataList)
{
}

template <typename ReferenceHash>
std::size_t ClusterHashMatchFinder<ReferenceHash>::findReadMatches(
    const Cluster& cluster,
    const flowcell::ReadMetadata& readMetadata,
    Matches& matches) const
{
    std::size_t repeatSeeds = 0;

    std::size_t readSeedIndex = 0;
    bool neighborlessMatchFound = false;
    BOOST_FOREACH(const alignment::SeedMetadata &seedMetadata, seedMetadataList_)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "findClusterMatches " << readMetadata << ":" << seedMetadata);
        if (seedMetadata.getReadIndex() != readMetadata.getIndex())
        {
            continue;
        }
        ++readSeedIndex;
        neighborlessMatchFound |= BaseT::findSeedMatches(cluster, seedMetadata, readMetadata, matches, repeatSeeds);
    }
    return repeatSeeds;
}

template <typename ReferenceHash>
std::size_t ClusterHashMatchFinder<ReferenceHash>::findClusterMatches(
    const Cluster& cluster,
    const flowcell::ReadMetadataList& readMetadataList,
    Matches& matches) const
{
    std::size_t repeatSeeds = 0;
    BOOST_FOREACH(const flowcell::ReadMetadata &readMetadata, readMetadataList)
    {
        repeatSeeds += findReadMatches(cluster, readMetadata, matches);
    }

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "findClusterMatches:" << matches.size() << " matches ");
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "findClusterMatches:" << repeatSeeds << " repeatSeeds ");

    return repeatSeeds;
}

template class SeedHashMatchFinder<reference::NumaReferenceHash<reference::ReferenceHash<oligo::ShortKmerType, common::NumaAllocator<void, 0> > > >;

template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::VeryShortKmerType> >;
} // namespace alignment
} // namespace isaac
