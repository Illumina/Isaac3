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

#include "workflow/alignWorkflow/BclBgzfDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

/////////////// BclBgzfBaseCallsSource Implementation
BclBgzfBaseCallsSource::BclBgzfBaseCallsSource(
    const flowcell::Layout &flowcell,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    common::ThreadVector &bclLoadThreads,
    const unsigned inputLoadersMax,
    const bool extractClusterXy):
    flowcell_(flowcell),
    bclLoadThreads_(bclLoadThreads),
    bciFilePath_(flowcell_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>()),
    filterFilePath_(flowcell_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::FiltersFilePathAttributeTag>()),
    positionsFilePath_(flowcell_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::PositionsFilePathAttributeTag>()),
    cycles_(flowcell::getTotalReadLength(flowcell_.getReadMetadataList()) + flowcell_.getBarcodeLength()),
    cycleBciFileBuf_(std::ios_base::binary|std::ios_base::in),
    tileClusterOffsets_(flowcell_.getAttribute<flowcell::Layout::BclBgzf, flowcell::TilesPerLaneMaxAttributeTag>()),
    cycleBciMappers_(
        flowcell::getMaxCycleNumber(flowcell_) + 1,
        rta::CycleBciMapper(tileClusterOffsets_.size())),
    tileBciIndexMap_(),
    flowcellTiles_(getTiles(flowcell, tileBciIndexMap_)),
    undiscoveredTiles_(flowcellTiles_.begin()),
    laneBciMapper_(tileClusterOffsets_.size()),
    threadReaders_(
        bclLoadThreads_.size(),
        rta::BclBgzfTileReader(
            flowcell_.getLongestAttribute<flowcell::Layout::BclBgzf, flowcell::BclFilePathAttributeTag>().string().size(),
            ignoreMissingBcls,
            flowcell::getMaxTileClusters(flowcellTiles_),
            tileBciIndexMap_,
            cycleBciMappers_)),
    bclMapper_(cycles_.size(),
               bclLoadThreads_, threadReaders_,
               inputLoadersMax, flowcell::getMaxTileClusters(flowcellTiles_)),
    filtersMapper_(ignoreMissingFilters),
    clocsMapper_(),
    locsMapper_(),
    barcodeLoader_(bclLoadThreads, inputLoadersMax, flowcell::getMaxTileClusters(flowcellTiles_), threadReaders_),
    currentFlowcellIndex_(-1U),
    currentLaneNumber_(-1U)
{
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions before filtersMapper_.reserveBuffer ")
    filtersMapper_.reserveBuffers(filterFilePath_.string().size(), flowcell::getMaxTileClusters(flowcellTiles_));
    ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")

    if (extractClusterXy)
    {
        clocsMapper_.reserveBuffers(positionsFilePath_.string().size(), flowcell::getMaxTileClusters(flowcellTiles_));
        locsMapper_.reserveBuffers(positionsFilePath_.string().size(), flowcell::getMaxTileClusters(flowcellTiles_));
        ISAAC_TRACE_STAT("SelectMatchesTransition::SelectMatchesTransitions after filtersMapper_.reserveBuffer ")
    }

}

flowcell::TileMetadataList BclBgzfBaseCallsSource::discoverTiles()
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

inline unsigned getLaneNumber(const flowcell::TileMetadataList& tiles)
{
    ISAAC_ASSERT_MSG(
        tiles.end()
            == std::find_if(
                tiles.begin(),
                tiles.end(),
                boost::bind(&flowcell::TileMetadata::getLane, _1)
                    != tiles.front().getLane()),
        "Expected all tiles to belong to the same lane");
    return tiles.front().getLane();
}

flowcell::TileMetadataList BclBgzfBaseCallsSource::getTiles(
    const flowcell::Layout &flowcellLayout,
    std::vector<unsigned> &tileBciIndexMap)
{
    flowcell::TileMetadataList tileMetadataList;

    const std::string &flowcellId = flowcellLayout.getFlowcellId();
    BOOST_FOREACH(const unsigned int lane, flowcellLayout.getLaneIds())
    {
        const std::vector<unsigned int> tileList = flowcellLayout.getTileIds(lane);

        boost::filesystem::path laneBciFilePath;
        flowcellLayout.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>(lane, laneBciFilePath);
        rta::LaneBciMapper mapper(
            flowcellLayout.getAttribute<flowcell::Layout::BclBgzf, flowcell::TilesPerLaneMaxAttributeTag>());
        mapper.mapFile(laneBciFilePath);

        BOOST_FOREACH(const unsigned int tileNumber, tileList)
        {
            const rta::LaneBciMapper::TileInfo tileInfo = mapper.getTileInfo(tileNumber);
            if (tileInfo.tileClusters_)
            {
                const flowcell::TileMetadata tileMetadata(
                    flowcellId, flowcellLayout.getIndex(),
                    tileNumber, lane,
                    tileInfo.tileClusters_,
                    tileMetadataList.size());
                tileMetadataList.push_back(tileMetadata);
                tileBciIndexMap.resize(std::max<unsigned>(tileBciIndexMap.size(), tileMetadata.getIndex() + 1));
                tileBciIndexMap[tileMetadata.getIndex()] = tileInfo.tileIndex_;
                ISAAC_THREAD_CERR << tileMetadata << std::endl;
            }
        }
    }

    return tileMetadataList;
}

void BclBgzfBaseCallsSource::initCycleBciMappers(
    const flowcell::Layout& flowcell,
    const std::vector<unsigned>& cycles,
    const unsigned laneNumber,
    std::vector<rta::CycleBciMapper>& cycleBciMappers)
{
    ISAAC_ASSERT_MSG(
        cycleBciMappers.size() > *std::max_element(cycles.begin(), cycles.end()),
        "Insufficient number of cycleBciMappers. Expected: >" << *std::max_element(cycles.begin(), cycles.end()) <<
        " have:" << cycleBciMappers.size());
    BOOST_FOREACH(const unsigned cycle, cycles)
    {
        flowcell.getLaneCycleAttribute<
        flowcell::Layout::BclBgzf,flowcell::BciFilePathAttributeTag>(laneNumber, cycle, bciFilePath_);

        std::istream is(cycleBciFileBuf_.reopen(bciFilePath_.c_str(), io::FileBufWithReopen::SequentialOnce));
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open file: " + bciFilePath_.string() + strerror(errno)));
        }

        cycleBciMappers.at(cycle).mapStream(is, bciFilePath_);
    }
}

void BclBgzfBaseCallsSource::initBciMappers(
    const flowcell::Layout &flowcell,
    const flowcell::TileMetadataList &allTiles,
    const unsigned laneNumber,
    std::vector<rta::CycleBciMapper>  &cycleBciMappers,
    std::vector<uint64_t> &tileClusterOffsets,
    std::vector<unsigned> &tileBciIndexMap)
{
    flowcell.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::BciFilePathAttributeTag>(laneNumber, bciFilePath_);
    laneBciMapper_.mapFile(bciFilePath_);
    tileBciIndexMap.resize(allTiles.size());
    BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, allTiles)
    {
        if (tileMetadata.getFlowcellIndex() == flowcell.getIndex() && laneNumber == tileMetadata.getLane())
        {
            const rta::LaneBciMapper::TileInfo tileInfo = laneBciMapper_.getTileInfo(tileMetadata.getTile());
            tileBciIndexMap.at(tileMetadata.getIndex()) = tileInfo.tileIndex_;
        }
    }

    unsigned prev = 0;
    tileClusterOffsets.resize(laneBciMapper_.getTilesCount());
    for (unsigned bciTileIndex = 0; bciTileIndex < tileClusterOffsets.size(); ++bciTileIndex)
    {
        tileClusterOffsets.at(bciTileIndex) = prev;
        prev += laneBciMapper_.getTileClusterCount(bciTileIndex);
    }

    cycles_.clear();
    cycles_.insert(cycles_.end(), flowcell.getBarcodeCycles().begin(), flowcell.getBarcodeCycles().end());
    cycles_.insert(cycles_.end(), flowcell.getDataCycles().begin(), flowcell.getDataCycles().end());

    initCycleBciMappers(flowcell, cycles_, laneNumber, cycleBciMappers);
}

// BarcodeSource implementation
void BclBgzfBaseCallsSource::loadBarcodes(
    const flowcell::Layout &flowcell,
    const unsigned unknownBarcodeIndex,
    const flowcell::TileMetadataList &tiles,
    demultiplexing::Barcodes &barcodes)
{
    initCycleBciMappers(flowcell, flowcell.getBarcodeCycles(), getLaneNumber(tiles), cycleBciMappers_);

    barcodeLoader_.loadBarcodes(unknownBarcodeIndex, flowcell, tiles, barcodes);
}

void BclBgzfBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loading Bcl data for " << tileMetadata << std::endl;

    if (currentLaneNumber_ != tileMetadata.getLane())
    {
        initBciMappers(flowcell_, flowcellTiles_, tileMetadata.getLane(), cycleBciMappers_, tileClusterOffsets_, tileBciIndexMap_);
        currentFlowcellIndex_ = flowcell_.getIndex();
        currentLaneNumber_ = tileMetadata.getLane();
    }

    bclMapper_.mapTile(flowcell_, tileMetadata);
    ISAAC_THREAD_CERR << "Loading Bcl data done for " << tileMetadata << std::endl;

    ISAAC_THREAD_CERR << "Loading Filter data for " << tileMetadata << std::endl;
    flowcell_.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::FiltersFilePathAttributeTag>(tileMetadata.getLane(), filterFilePath_);
    filtersMapper_.mapTile(filterFilePath_, tileMetadata.getClusterCount(),
                           tileClusterOffsets_.at(tileBciIndexMap_.at(tileMetadata.getIndex())));
    ISAAC_THREAD_CERR << "Loading Filter data done for " << tileMetadata << std::endl;

    bool boolUseLocsPositions = false;
    if (bclData.storeXy())
    {
        ISAAC_THREAD_CERR << "Loading Positions data for " << tileMetadata << std::endl;
        flowcell_.getLaneAttribute<flowcell::Layout::BclBgzf, flowcell::PositionsFilePathAttributeTag>(tileMetadata.getLane(), positionsFilePath_);
        if (flowcell::isClocsPath(positionsFilePath_))
        {
            ISAAC_ASSERT_MSG(false, "Clocs are not supported with bcl.bgzf data");
            clocsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount());
        }
        else
        {
            const bool patternedFlowcell = flowcell_.getAttribute<flowcell::Layout::BclBgzf, flowcell::PatternedFlowcellAttributeTag>();
            locsMapper_.mapTile(positionsFilePath_, tileMetadata.getClusterCount(),
                                patternedFlowcell ? 0 : tileClusterOffsets_.at(tileBciIndexMap_.at(tileMetadata.getIndex())));
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

void BclBgzfBaseCallsSource::resetBclData(
    const flowcell::TileMetadata& tileMetadata,
    alignment::BclClusters& bclData) const
{
    ISAAC_THREAD_CERR<< "Resetting Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
    bclData.reset(bclMapper_.getCyclesCount(), tileMetadata.getClusterCount(), true);
    ISAAC_THREAD_CERR << "Resetting Bcl data done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
}

void BclBgzfBaseCallsSource::bclToClusters(
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
            ISAAC_ASSERT_MSG(false, "Clocs are not supported with bcl.bgzf data");
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
