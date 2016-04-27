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
 ** \file TemplateDetector.hh
 **
 ** \brief Selection the best matches among all possible candidates.
 **
 ** \author Roman Petrovski
 **/


#ifndef iSAAC_ALIGNMENT_TEMPLATE_DETECTOR_HH
#define iSAAC_ALIGNMENT_TEMPLATE_DETECTOR_HH

#include <stddef.h>

#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/matchSelector/MatchSelectorStats.hh"
#include "alignment/TemplateBuilder.hh"
#include "common/Threads.hpp"
#include "flowcell/Layout.hh"
#include "reference/Contig.hh"
#include "reference/KUniqueness.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class TemplateDetector: boost::noncopyable
{
public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;
    /// Construction of an instance for a given reference
    TemplateDetector(
        common::UnsafeThreadVector &computeThreads,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const reference::NumaContigLists &contigLists,
        const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations,
        boost::ptr_vector<TemplateBuilder> &threadTemplateBuilders,
        const int mateDriftRange,
        const TemplateLengthStatistics &defaultTemplateLengthStatistics,
        const bool perTileTls);

    template <typename MatchFinderT>
    void determineTemplateLengths(
        const flowcell::TileMetadata &tileMetadata,
        const matchFinder::ClusterInfos &clusterInfos,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        std::vector<alignment::TemplateLengthStatistics> &templateLengthStatistics,
        matchSelector::MatchSelectorStats &stats);

private:
    // The threading code in selectTileMatches can't deal with exception cleanup. Let it just crash for now.
    common::UnsafeThreadVector &computeThreads_;

    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const reference::NumaContigLists &contigLists_;
    const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations_;

    const TemplateLengthStatistics userTemplateLengthStatistics_;
    const bool perTileTls_;

    std::vector<Cluster> threadCluster_;
    boost::ptr_vector<TemplateBuilder> &threadTemplateBuilders_;
    std::vector<TemplateLengthDistribution> templateLengthDistributions_;

    mutable boost::mutex mutex_;
    unsigned unprocessedClusterId_;
    unsigned pendingClusterId_;
    boost::condition_variable stateChangedCondition_;

    template <typename MatchFinderT>
    void templateLengthThread(
        const unsigned threadNumber,
        const flowcell::TileMetadata& tileMetadata,
        const BclClusters& bclData,
        const matchFinder::ClusterInfos& clusterInfos,
        const MatchFinderT& matchFinder,
        std::size_t &statsToBuild,
        std::vector<alignment::TemplateLengthStatistics>& templateLengthStatistics);

    static const unsigned CLUSTERS_AT_A_TIME = 10000;
    typedef std::pair<unsigned, TemplateLengthDistribution::AlignmentModel> BarcodeAlignmentModel;

    template <typename MatchFinderT>
    void collectModels(
        const unsigned clusterRangeBegin,
        const unsigned clusterRangeEnd,
        const BclClusters& bclData,
        const matchFinder::ClusterInfos& clusterInfos,
        const flowcell::ReadMetadataList& tileReads,
        const flowcell::TileMetadata& tileMetadata,
        const unsigned barcodeLength,
        const unsigned readNameLength,
        const SeedMetadataList& tileSeeds,
        std::size_t& statsToBuild,
        std::vector<alignment::TemplateLengthStatistics>& templateLengthStatistics,
        Cluster& ourThreadCluster,
        TemplateBuilder& ourThreadTemplateBuilder,
        const MatchFinderT& matchFinder,
        common::StaticVector<BarcodeAlignmentModel, CLUSTERS_AT_A_TIME>& threadBarcodeModels);
};

} // namespace matchSelector
} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
