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
 ** \file MatchSelector.hh
 **
 ** \brief Selection the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_HH

#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/BclClusters.hh"
#include "alignment/Cluster.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/TemplateBuilder.hh"
#include "alignment/Match.hh"
#include "alignment/SeedMetadata.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "alignment/matchSelector/MatchSelectorStats.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"
#include "alignment/matchSelector/TemplateDetector.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/AnnotationLoader.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;

class MatchSelector: boost::noncopyable
{
public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;
    /// Construction of an instance for a given reference
    MatchSelector(
        unsigned int maxThreadCount,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const reference::NumaContigLists &contigLists,
        const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations,
        const unsigned repeatThreshold,
        const int mateDriftRange,
        const TemplateLengthStatistics &defaultTemplateLengthStatistics,
        const unsigned mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
        const bool collectCycleStats,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const bool clipOverlapping,
        const bool scatterRepeats,
        const bool rescueShadows,
        const bool anchorMate,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const bool smartSmithWaterman,
        const bool noSmithWaterman,
        const bool splitAlignments,
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned splitGapLength,
        const TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
        const bool reserveBuffers);

    /**
     * \brief frees the major memory reservations to make it safe to use dynamic memory allocations again
     */
    void unreserve()
    {
        threadTemplateBuilders_.clear();
        std::vector<Cluster>().swap(threadCluster_);
        std::vector<matchSelector::MatchSelectorStats>().swap(threadStats_);
    }

    void dumpStats(const boost::filesystem::path &statsXmlPath);
    void reserveMemory(
        const flowcell::TileMetadataList &tileMetadataList);

    template <typename MatchFinderT>
    void parallelSelect(
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const flowcell::TileMetadata &tileMetadata,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        matchSelector::FragmentStorage &fragmentStorage);

private:
    // The threading code in selectTileMatches can't deal with exception cleanup. Let it just crash for now.
    common::UnsafeThreadVector computeThreads_;

    TileMetadataList tileMetadataList_;
    /**
     * \brief threadBclFilePaths_ gets resized for every tile total readlength. If the tile read lengths
     *        changes from lower to bigger, more threadBclFilePaths_ strings get allocated which breaks the whole
     *        concept of allocating things once. For now this list contains tiles in the processing order so
     *        that the total read length goes only down. TODO: cleanup this mess for example by creating
     *        MatchSelector only for the group of tiles that have the same geometry.
     */
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const reference::NumaContigLists &contigLists_;
    const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations_;
    const unsigned repeatThreshold_;

    const unsigned mapqThreshold_;
    const bool pfOnly_;
    const bool collectCycleStats_;
    const unsigned baseQualityCutoff_;
    const bool keepUnaligned_;
    const bool clipSemialigned_;
    const bool clipOverlapping_;
    const std::vector<matchSelector::SequencingAdapterList> barcodeSequencingAdapters_;

    std::vector<matchSelector::MatchSelectorStats> allStats_;
    std::vector<matchSelector::MatchSelectorStats> threadStats_;

    std::vector<Cluster> threadCluster_;
    boost::ptr_vector<TemplateBuilder> threadTemplateBuilders_;
    std::vector<matchSelector::SemialignedEndsClipper> threadSemialignedEndsClippers_;
    std::vector<matchSelector::OverlappingEndsClipper> threadOverlappingEndsClippers_;
    // updated for barcodes relevant for the current tile
    std::vector<RestOfGenomeCorrection> restOfGenomeCorrections_;

    matchSelector::TemplateDetector templateDetector_;

    mutable boost::mutex mutex_;

    template <typename MatchFinderT>
    void alignThread(
        const unsigned threadNumber,
        const flowcell::TileMetadata & tileMetadata,
        const matchFinder::ClusterInfos &clusterInfos,
        unsigned &clusterId,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        const std::vector<TemplateLengthStatistics> & templateLengthStatistics,
        matchSelector::FragmentStorage &fragmentStorage);


    /**
     * \brief Construct the contig list from the SortedReference XML
     */
    reference::ContigList getContigList(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata) const;

    /**
     ** \brief Helper method to generate the 'rest of the genome' correction for
     ** uniquely aligned reads and fragments.
     **
     ** There is one value for each individual reads in the readMetadataList (at
     ** the corresponding location) and one additional value for cases whare all
     ** the reads match uniquely.
     **/
    std::vector<double> getRestOfGenomeCorrectionList(
        const std::vector<flowcell::ReadMetadata> &readMetadataList) const;

    template <typename MatchFinderT>
    matchSelector::TemplateAlignmentType alignCluster(
        const reference::ContigList& barcodeContigList,
        const isaac::reference::ContigAnnotations& barcodeKUniqeness,
        const flowcell::ReadMetadataList& tileReads,
        const SeedMetadataList& tileSeeds,
        const matchSelector::SequencingAdapterList& sequencingAdapters,
        const TemplateLengthStatistics& templateLengthStatistics,
        const uint64_t barcodeIndex,
        const MatchFinderT &matchFinder,
        const RestOfGenomeCorrection& restOfGenomeCorrection,
        const unsigned threadNumber, TemplateBuilder& templateBuilder,
        const Cluster& cluster, BamTemplate& bamTemplate,
        matchSelector::MatchSelectorStats& stats,
        matchSelector::FragmentStorage &fragmentStorage);
    void templateLengthThread(
        const unsigned threadNumber,
        const flowcell::TileMetadata& tileMetadata,
        const BclClusters& bclData,
        const matchFinder::ClusterInfos& clusterInfos,
        alignment::SeedHashMatchFinder<oligo::ShortKmerType>& matchFinder,
        std::size_t &statsToBuild,
        std::vector<alignment::TemplateLengthStatistics>& templateLengthStatistics);

    static const unsigned CLUSTERS_AT_A_TIME = 10000;
    typedef std::pair<unsigned, TemplateLengthDistribution::AlignmentModel> BarcodeAlignmentModel;
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
        alignment::SeedHashMatchFinder<oligo::ShortKmerType>& matchFinder,
        common::StaticVector<BarcodeAlignmentModel, CLUSTERS_AT_A_TIME>& threadBarcodeModels);
};

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
