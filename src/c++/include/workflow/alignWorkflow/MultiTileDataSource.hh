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
 ** \file MultiTileDataSource.hh
 **
 ** \brief Wrapper for presenting more than one bcl tile as a single tile to alignment.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_MULTITILE_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_MULTITILE_DATA_SOURCE_HH

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

template <typename BaseCallsSourceT>
class MultiTileBaseCallsSource : public BarcodeSource, public TileSource
{
    const unsigned tilesPerChunk_;
    const flowcell::Layout &flowcell_;
    BaseCallsSourceT &tileDataSource_;
    const flowcell::TileMetadataList physicalTiles_;
    // [original tile index]: first is chunk index, second is the cluster id offset in the chunk
    std::vector<std::pair<unsigned, unsigned> > physicalToLogical_;
    flowcell::TileMetadataList::const_iterator undiscoveredTile_;
    unsigned undiscoveredIndex_;

public:
    MultiTileBaseCallsSource(
        const unsigned tilesPerChunk,
        const flowcell::Layout &flowcell,
        BaseCallsSourceT &baseCallsSource):
            tilesPerChunk_(tilesPerChunk),
            flowcell_(flowcell),
            tileDataSource_(baseCallsSource),
            physicalTiles_(discoverAllTiles(tileDataSource_)),
            physicalToLogical_(physicalTiles_.size(), std::make_pair(-1U, 0U)),
            undiscoveredTile_(physicalTiles_.begin()),
            undiscoveredIndex_(0)
    {
    }

    unsigned getMaxTileClusters() const { return tileDataSource_.getMaxTileClusters() * tilesPerChunk_;}

    flowcell::TileMetadataList discoverTiles()
    {
        flowcell::TileMetadataList ret;
        if (physicalTiles_.end() == undiscoveredTile_)
        {
            return ret;
        }

        do
        {
            flowcell::TileMetadata chunk(*undiscoveredTile_, undiscoveredIndex_++);
            ++undiscoveredTile_;
            for (unsigned chunkTiles = 1;
                physicalTiles_.end() != undiscoveredTile_ && chunkTiles < tilesPerChunk_ &&
                    chunk.getLane() == undiscoveredTile_->getLane() &&
                    chunk.getFlowcellIndex() == undiscoveredTile_->getFlowcellIndex();
                ++chunkTiles, ++undiscoveredTile_)
            {
                chunk.setClusterCount(chunk.getClusterCount() + undiscoveredTile_->getClusterCount());
            }
            ret.push_back(chunk);
        }
        while (physicalTiles_.end() != undiscoveredTile_ && ret.back().getLane() == undiscoveredTile_->getLane());

        return ret;
    }

    // prepare bclData buffers to receive new tile data
    void resetBclData(
        const flowcell::TileMetadata& tileMetadata,
        alignment::BclClusters& bclData) const
    {
        ISAAC_THREAD_CERR<< "Resetting Bcl data for " << tileMetadata.getClusterCount() << " bcl clusters" << std::endl;
        bclData.reset(flowcell::getTotalReadLength(flowcell_.getReadMetadataList()) + flowcell_.getBarcodeLength(), tileMetadata.getClusterCount(), true);
        ISAAC_THREAD_CERR << "Resetting Bcl data done for " << bclData.getClusterCount() << " bcl clusters" << std::endl;
    }

    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData)
    {
        const flowcell::TileMetadataList::const_iterator originalBegin =
            physicalTiles_.begin() + tileMetadata.getOriginalIndex();
        const flowcell::TileMetadataList::const_iterator originalEnd =
            originalBegin + std::min<std::size_t>(tilesPerChunk_, std::distance(originalBegin, physicalTiles_.end()));
        for (flowcell::TileMetadataList::const_iterator it = originalBegin;
            it != originalEnd && it->getLane() == originalBegin->getLane() &&
                it->getFlowcellIndex() == originalBegin->getFlowcellIndex();
            ++it)
        {
            tileDataSource_.loadClusters(*it, bclData);
        }
    }

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const flowcell::Layout &flowcell,
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        demultiplexing::Barcodes &barcodes)
    {
        flowcell::TileMetadataList physicalTiles;
        BOOST_FOREACH(const flowcell::TileMetadata &tileMetadata, tiles)
        {
            const flowcell::TileMetadataList::const_iterator originalBegin =
                physicalTiles_.begin() + tileMetadata.getOriginalIndex();
            const flowcell::TileMetadataList::const_iterator originalEnd =
                originalBegin + std::min<std::size_t>(tilesPerChunk_, std::distance(originalBegin, physicalTiles_.end()));
            unsigned clusterOffset = 0;
            for (flowcell::TileMetadataList::const_iterator it = originalBegin;
                it != originalEnd &&
                    it->getLane() == originalBegin->getLane() &&
                    it->getFlowcellIndex() == originalBegin->getFlowcellIndex();
                ++it)
            {
                // fill the mapping so that barcodes can be correctly remapped from the data source tiles into
                // the tiles that the client code sees.
                physicalToLogical_.at(it->getOriginalIndex()).first = tileMetadata.getIndex();
                physicalToLogical_.at(it->getOriginalIndex()).second = clusterOffset;
                clusterOffset += it->getClusterCount();
                physicalTiles.push_back(*it);
            }
        }

        tileDataSource_.loadBarcodes(flowcell, unknownBarcodeIndex, physicalTiles, barcodes);
        // remap the barcode tile indexes back to what the client code expects
        BOOST_FOREACH(demultiplexing::Barcode &barcode, barcodes)
        {
            ISAAC_ASSERT_MSG(-1U != physicalToLogical_.at(barcode.getTile()).first, "fail: " << barcode);
            barcode = demultiplexing::Barcode(
                barcode.getSequence(),
                demultiplexing::BarcodeId(
                    physicalToLogical_.at(barcode.getTile()).first,
                    barcode.getBarcode(),
                    physicalToLogical_.at(barcode.getTile()).second + barcode.getCluster(),
                    barcode.getMismatches()));
        }
    }

private:

    static flowcell::TileMetadataList discoverAllTiles(BaseCallsSourceT &tileSource)
    {
        flowcell::TileMetadataList ret;
        for (flowcell::TileMetadataList tiles = tileSource.discoverTiles();
            !tiles.empty(); tiles = tileSource.discoverTiles())
        {
            ret.insert(ret.end(), tiles.begin(), tiles.end());
        }
        return ret;
    }
};

template <>
template <typename BaseCallsSourceT>
struct DataSourceTraits<MultiTileBaseCallsSource<BaseCallsSourceT> > : public DataSourceTraits<BaseCallsSourceT>
{
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_MULTITILE_DATA_SOURCE_HH
