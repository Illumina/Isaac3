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
 ** \file BclDataSource.cpp
 **
 ** \brief see BclDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "workflow/alignWorkflow/BclDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

BclTileSource::BclTileSource(
    const flowcell::Layout &bclFlowcellLayout) :
    flowcellTiles_(getTiles(bclFlowcellLayout)),
    undiscoveredTiles_(flowcellTiles_.begin())
{

}

flowcell::TileMetadataList BclTileSource::getTiles(const flowcell::Layout &flowcellLayout) const
{
    flowcell::TileMetadataList tileMetadataList;

    const std::string &flowcellId = flowcellLayout.getFlowcellId();
    BOOST_FOREACH(const unsigned int lane, flowcellLayout.getLaneIds())
    {
        const std::vector<unsigned int> tileList = flowcellLayout.getTileIds(lane);
        unsigned tilesFound = 0;
        BOOST_FOREACH(const unsigned int tile, tileList)
        {
            boost::filesystem::path bclFilePath;
            BOOST_FOREACH(const unsigned cycle, flowcellLayout.getDataCycles())
            {
                flowcellLayout.getLaneTileCycleAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>(
                    lane, tile, cycle, bclFilePath);
                if (boost::filesystem::exists(bclFilePath))
                {
                    const unsigned int clusterCount = rta::BclMapper::getClusterCount(bclFilePath);
                    const flowcell::TileMetadata tileMetadata(
                        flowcellId, flowcellLayout.getIndex(),
                        tile, lane,
                        clusterCount,
                        tileMetadataList.size());
                    tileMetadataList.push_back(tileMetadata);
                    ++tilesFound;
                    break;
                }
            }
        }
        if (!tilesFound)
        {
            ISAAC_THREAD_CERR << "WARNING: No tiles found for lane " << lane << std::endl;
        }
    }

    if (tileMetadataList.empty())
    {
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(std::string("No tile data found for flowcell ") + boost::lexical_cast<std::string>(flowcellLayout)));
    }

    return tileMetadataList;
}


// TileSource implementation
flowcell::TileMetadataList BclTileSource::discoverTiles()
{
    flowcell::TileMetadataList ret;
    while (flowcellTiles_.end() != undiscoveredTiles_)
    {
        ret.push_back(*undiscoveredTiles_);
        if (ret.front().getLane() != ret.back().getLane())
        {
            ret.pop_back();
            break;
        }
        ++undiscoveredTiles_;
    }
    return ret;
}


// BarcodeSource implementation
void BclBaseCallsSource::loadBarcodes(
    const flowcell::Layout &flowcell,
    const unsigned unknownBarcodeIndex,
    const flowcell::TileMetadataList &tiles,
    demultiplexing::Barcodes &barcodes)
{
    barcodeLoader_.loadBarcodes(unknownBarcodeIndex, flowcell, tiles, barcodes);
}

/////////////// BclBaseCallsSource Implementation
BclBaseCallsSource::BclBaseCallsSource(
    const flowcell::Layout &flowcell,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    common::ThreadVector &bclLoadThreads,
    const unsigned inputLoadersMax,
    const bool extractClusterXy):
    flowcell_(flowcell),
    tileSource_(flowcell_),
    bclLoadThreads_(bclLoadThreads),
    filterFilePath_(flowcell_.getLongestAttribute<flowcell::Layout::Bcl, flowcell::FiltersFilePathAttributeTag>()),
    positionsFilePath_(flowcell_.getLongestAttribute<flowcell::Layout::Bcl, flowcell::PositionsFilePathAttributeTag>()),
    threadReaders_(bclLoadThreads_.size(), rta::BclReader(
        flowcell_.getLongestAttribute<flowcell::Layout::Bcl, flowcell::BclFilePathAttributeTag>().string().size(),
        ignoreMissingBcls, tileSource_.getMaxTileClusters())),
        bclMapper_(flowcell::getTotalReadLength(flowcell_.getReadMetadataList()) + flowcell_.getBarcodeLength(),
               bclLoadThreads_, threadReaders_,
               inputLoadersMax, tileSource_.getMaxTileClusters()),
    filtersMapper_(ignoreMissingFilters),
    clocsMapper_(),
    locsMapper_(),
    barcodeLoader_(bclLoadThreads, inputLoadersMax, tileSource_.getMaxTileClusters(), threadReaders_)
{
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before filtersMapper_.reserveBuffer ")
    filtersMapper_.reserveBuffers(filterFilePath_.string().size(), tileSource_.getMaxTileClusters());
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")

    if (extractClusterXy)
    {
        clocsMapper_.reserveBuffers(positionsFilePath_.string().size(), tileSource_.getMaxTileClusters());
        locsMapper_.reserveBuffers(positionsFilePath_.string().size(), tileSource_.getMaxTileClusters());
        ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")
    }
}

void BclBaseCallsSource::loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

    bclMapper_.mapTile(flowcell_, tileMetadata);
    ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
    flowcell_.getLaneTileAttribute<flowcell::Layout::Bcl, flowcell::FiltersFilePathAttributeTag>(
        tileMetadata.getLane(), tileMetadata.getTile(), filterFilePath_);
    filtersMapper_.mapTile(filterFilePath_, tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

    bool boolUseLocsPositions = false;
    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Loading Positions data for " << tileMetadata << std::endl;
        flowcell_.getLaneTileAttribute<flowcell::Layout::Bcl, flowcell::PositionsFilePathAttributeTag>(
            tileMetadata.getLane(), tileMetadata.getTile(), positionsFilePath_);
        if (flowcell::isClocsPath(positionsFilePath_))
        {
            clocsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
        }
        else
        {
            locsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
            boolUseLocsPositions = true;
        }
        ISAAC_THREAD_CERR << "Loading Positions data done for " << tileMetadata << std::endl;
    }

    // bclToClusters mainly does transposition of bcl cycles to clusters which is a non-io operation.
    // However, the amount of CPU required is relatively low, and occurs on a single thread.
    // Avoid locking all the cores for the duration of this...
    // Also, bclMapper_ and filtersMapper_ are shared between the threads at the moment.
    bclToClusters(tileMetadata, bclData, boolUseLocsPositions);
}


void BclBaseCallsSource::resetBclData(
    const flowcell::TileMetadata& tileMetadata,
    alignment::BclClusters& bclData) const
{
    ISAAC_THREAD_CERR<< "Resetting Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.reset(bclMapper_.getCyclesCount(), tileMetadata.getClusterCount(), true);
    ISAAC_THREAD_CERR << "Resetting Bcl data done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
}

void BclBaseCallsSource::bclToClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData,
    const bool useLocsPositions) const
{

    ISAAC_THREAD_CERR << "Transposing Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    const clock_t startTranspose = clock();
    bclMapper_.transpose(bclData.addMoreClusters(tileMetadata.getClusterCount()));
    ISAAC_THREAD_CERR << "Transposing Bcl data done for " << bclData.getClusterCount() << " bcl clusters in " << (clock() - startTranspose) / 1000 << "ms" << std::endl;

    ISAAC_THREAD_CERR << "Extracting Pf values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    // gcc 4.4 has trouble figuring out which assignment implementation to use with back insert iterators
    filtersMapper_.getPf(std::back_inserter(bclData.pf()));

    ISAAC_ASSERT_MSG(bclData.pf().size() == bclData.getClusterCount(), "Mismatch between data " << bclData.getClusterCount() << " and pf " << bclData.pf().size() << "counts");
    ISAAC_THREAD_CERR << "Extracting Pf values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;

    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Extracting Positions values for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
        if (!useLocsPositions)
        {
            clocsMapper_.getPositions(std::back_inserter(bclData.xy()));
        }
        else
        {
            locsMapper_.getPositions(std::back_inserter(bclData.xy()));
        }
        ISAAC_ASSERT_MSG(bclData.xy().size() == bclData.getClusterCount(), "Mismatch between data " << bclData.getClusterCount() << " and position " << bclData.xy().size() << "counts");
        ISAAC_THREAD_CERR << "Extracting Positions values done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
    }

}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
