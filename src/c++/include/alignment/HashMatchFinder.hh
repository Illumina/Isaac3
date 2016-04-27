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
 ** \file FindHashMatchesTransition.hh
 **
 ** \brief Top level component to control the analysis process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH
#define iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH

#include "alignment/Cluster.hh"
#include "alignment/Seed.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/Match.hh"
#include "common/Threads.hpp"
#include "demultiplexing/BarcodeLoader.hh"
#include "demultiplexing/BarcodeResolver.hh"
#include "demultiplexing/DemultiplexingStats.hh"
#include "flowcell/Layout.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/ReferenceHash.hh"
#include "reference/ReferenceMetadata.hh"
#include "reference/SortedReferenceMetadata.hh"

#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/DataSource.hh"
#include "workflow/alignWorkflow/FoundMatchesMetadata.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;


template <typename ReferenceHash>
class SeedHashMatchFinder
{
    typedef typename ReferenceHash::KmerT KmerT;
public:
    SeedHashMatchFinder(
        const ReferenceHash& referenceHash,
        const unsigned seedBaseQualityMin,
        const unsigned repeatThreshold);

    bool findSeedMatches(
        const Cluster& cluster,
        const alignment::SeedMetadata &seedMetadata,
        const flowcell::ReadMetadata& readMetadata,
        Matches& matches,
        std::size_t& repeatSeeds) const;

private:
    const ReferenceHash& referenceHash_;
    const unsigned seedBaseQualityMin_;
    const unsigned repeatThreshold_;
    const unsigned noExtendRepeatThreshold_;

    bool extendSeed(
        const Cluster& cluster,
        alignment::Seed<KmerT> seed,
        const alignment::SeedMetadata &seedMetadata,
        const flowcell::ReadMetadata &readMetadata,
        const ReferenceHash &referenceHash,
        Matches& matches,
        std::size_t &repeatSeeds) const;

    bool storeFwExtensionMatches(
        const typename ReferenceHash::MatchRange& oriPositions,
        const typename ReferenceHash::MatchRange& extPositions,
        const uint64_t extSeedOffset,
        const Seed<KmerT> &oriSeed,
        const Seed<KmerT> &extSeed,
        const unsigned matchesMax,
        Matches& matches) const;

    bool storeRvExtensionMatches(
        const typename ReferenceHash::MatchRange& oriPositions,
        const typename ReferenceHash::MatchRange& extPositions,
        const uint64_t extSeedOffset,
        const Seed<KmerT> &oriSeed,
        const Seed<KmerT> &extSeed,
        const unsigned matchesMax,
        Matches& matches) const;

};


template <typename ReferenceHash>
class ClusterHashMatchFinder : SeedHashMatchFinder<ReferenceHash>
{
    typedef SeedHashMatchFinder<ReferenceHash> BaseT;
public:
    ClusterHashMatchFinder(
        const ReferenceHash& referenceHash,
        const unsigned seedBaseQualityMin,
        const unsigned repeatThreshold,
        const alignment::SeedMetadataList &seedMetadataList);

    std::size_t findClusterMatches(
        const Cluster& cluster,
        const flowcell::ReadMetadataList& readMetadataList,
        Matches& matches) const;

    std::size_t findReadMatches(
        const Cluster& cluster,
        const flowcell::ReadMetadata& readMetadata,
        Matches& matches) const;

private:
    const alignment::SeedMetadataList &seedMetadataList_;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH
