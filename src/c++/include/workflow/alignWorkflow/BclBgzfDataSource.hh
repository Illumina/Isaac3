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
 ** \file BclBgzfDataSource.hh
 **
 ** \brief Encapsulation of BaseCalls folder with bcl.bgzf and filter files as seed and cluster data source
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH

#include "demultiplexing/BarcodeLoader.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/LocsMapper.hh"
#include "io/ClocsMapper.hh"
#include "io/FiltersMapper.hh"
#include "rta/BclBgzfTileReader.hh"
#include "rta/BclMapper.hh"
#include "rta/LaneBciMapper.hh"
#include "workflow/alignWorkflow/DataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class BclBgzfBaseCallsSource : public BarcodeSource, public TileSource
{
    const flowcell::Layout &flowcell_;
    common::ThreadVector &bclLoadThreads_;

    // temporaries to avoid memory allocations during data processing
    boost::filesystem::path bciFilePath_;
    boost::filesystem::path filterFilePath_;
    boost::filesystem::path positionsFilePath_;
    std::vector<unsigned> cycles_;
    io::FileBufWithReopen cycleBciFileBuf_;

    // References to cycleBciMappers_ and tileBciIndexMap_ get passed to threadReaders_ during construction.
    // The contents is dynamically swapped each time we switch flowcell lane
    /// Cumulative offsets of each tile in clusters. All tiles of the lane included
    std::vector<uint64_t> tileClusterOffsets_;
    std::vector<rta::CycleBciMapper> cycleBciMappers_;
    std::vector<unsigned> tileBciIndexMap_;
    const flowcell::TileMetadataList flowcellTiles_;
    flowcell::TileMetadataList::const_iterator undiscoveredTiles_;

    rta::LaneBciMapper laneBciMapper_;

    std::vector<rta::BclBgzfTileReader> threadReaders_;
    rta::ParallelBclMapper<rta::BclBgzfTileReader> bclMapper_;
    io::FiltersMapper filtersMapper_;
    io::ClocsMapper clocsMapper_;
    io::LocsMapper locsMapper_;
    demultiplexing::BarcodeLoader<rta::BclBgzfTileReader> barcodeLoader_;

    unsigned currentFlowcellIndex_;
    unsigned currentLaneNumber_;

public:
    BclBgzfBaseCallsSource(
        const flowcell::Layout &flowcell,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        common::ThreadVector &bclLoadThreads,
        const unsigned inputLoadersMax,
        const bool extractClusterXy);

    unsigned getMaxTileClusters() const { return flowcell::getMaxTileClusters(flowcellTiles_);}

    flowcell::TileMetadataList discoverTiles();

    // prepare bclData buffers to receive new tile data
    void resetBclData(
        const flowcell::TileMetadata& tileMetadata,
        alignment::BclClusters& bclData) const;

    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const flowcell::Layout &flowcell,
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        demultiplexing::Barcodes &barcodes);

private:
    /**
     * \return vector of tiles ordered by: flowcellId_, lane_, tile_
     */
    static flowcell::TileMetadataList getTiles(
        const flowcell::Layout &flowcellLayout,
        std::vector<unsigned> &tileBciIndexMap);

    void initBciMappers(
        const flowcell::Layout &flowcell,
        const flowcell::TileMetadataList &allTiles,
        const unsigned laneNumber,
        std::vector<rta::CycleBciMapper>  &cycleBciMappers,
        std::vector<uint64_t> &tileClusterOffsets,
        std::vector<unsigned> &tileBciIndexMap);

    void bclToClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData,
        const bool useLocsPositions) const;

    void initCycleBciMappers(
        const flowcell::Layout& flowcell,
        const std::vector<unsigned>& cycles,
        const unsigned laneNumber,
        std::vector<rta::CycleBciMapper>& cycleBciMappers);
};

template <>
struct DataSourceTraits<BclBgzfBaseCallsSource>
{
    static const bool SUPPORTS_XY = true;
};


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCLBGZF_DATA_SOURCE_HH
