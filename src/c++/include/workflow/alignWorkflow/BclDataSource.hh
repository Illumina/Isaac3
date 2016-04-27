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
 ** \file BclDataSource.hh
 **
 ** \brief Encapsulation of BaseCalls folder with bcl and filter files as seed and cluster data source
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH

#include "demultiplexing/BarcodeLoader.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/LocsMapper.hh"
#include "io/ClocsMapper.hh"
#include "io/FiltersMapper.hh"
#include "rta/BclReader.hh"
#include "rta/BclMapper.hh"
#include "workflow/alignWorkflow/DataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class BclTileSource: public TileSource
{
    const flowcell::TileMetadataList flowcellTiles_;
    // holds the state across multiple discoverTiles calls
    flowcell::TileMetadataList::const_iterator undiscoveredTiles_;
public:
    BclTileSource(
        const flowcell::Layout &bclFlowcellLayout);
    const flowcell::TileMetadataList &flowcellTiles() const {return flowcellTiles_;}

    flowcell::TileMetadataList discoverTiles();

    unsigned getMaxTileClusters() const
    {
        return flowcellTiles_.empty() ?
            0 : std::max_element(flowcellTiles_.begin(), flowcellTiles_.end(),
                                 boost::bind(&flowcell::TileMetadata::getClusterCount, _1)<
                                 boost::bind(&flowcell::TileMetadata::getClusterCount, _2))->getClusterCount();
    }

private:
    /**
     * \return vector of tiles ordered by: flowcellId_, lane_, tile_
     */
    flowcell::TileMetadataList getTiles(const flowcell::Layout &flowcellLayout) const;
};

class BclBaseCallsSource : public BarcodeSource
{
    const flowcell::Layout &flowcell_;
    BclTileSource tileSource_;
    common::ThreadVector &bclLoadThreads_;
    boost::filesystem::path filterFilePath_;
    boost::filesystem::path positionsFilePath_;
    std::vector<rta::BclReader> threadReaders_;
    rta::ParallelBclMapper<rta::BclReader> bclMapper_;
    io::FiltersMapper filtersMapper_;
    io::ClocsMapper clocsMapper_;
    io::LocsMapper locsMapper_;
    demultiplexing::BarcodeLoader<rta::BclReader> barcodeLoader_;


public:
    BclBaseCallsSource(
        const flowcell::Layout &flowcell,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        common::ThreadVector &bclLoadThreads,
        const unsigned inputLoadersMax,
        const bool extractClusterXy);

    unsigned getMaxTileClusters() const { return tileSource_.getMaxTileClusters();}

    flowcell::TileMetadataList discoverTiles() {return tileSource_.discoverTiles();}

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
    void bclToClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData,
        const bool useLocsPositions) const;
};

template <>
struct DataSourceTraits<BclBaseCallsSource>
{
    static const bool SUPPORTS_XY = true;
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BCL_DATA_SOURCE_HH
