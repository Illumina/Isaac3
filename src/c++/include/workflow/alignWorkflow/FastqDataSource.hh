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
 ** \file FastqDataSource.hh
 **
 ** \brief Encapsulation of single-ended and paired data stored in fastq file(s)
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH

#include "alignment/BclClusters.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/FastqLayout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/FastqLoader.hh"
#include "workflow/alignWorkflow/DataSource.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class FastqBaseCallsSource : public TileSource, public BarcodeSource
{
    const unsigned tileClustersMax_;
    const unsigned coresMax_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::Layout &fastqFlowcellLayout_;
    const unsigned clusterLength_;
    alignment::BclClusters clusters_;
    flowcell::TileMetadataList loadedTiles_;
    const std::vector<unsigned> lanes_;
    std::vector<unsigned>::const_iterator currentLaneIterator_;
    unsigned currentTile_;
    io::FastqLoader fastqLoader_;

public:
    FastqBaseCallsSource(
        const unsigned clustersAtATimeMax,
        const bool allowVariableLength,
        const unsigned coresMax,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::Layout &fastqFlowcellLayout,
        common::ThreadVector &threads);

    // TileSource implementation
    flowcell::TileMetadataList discoverTiles();

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const flowcell::Layout &flowcell,
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        demultiplexing::Barcodes &barcodes)
    {
        ISAAC_ASSERT_MSG(false, "Barcode resolution is not implemented for Fastq data");
    }

    // prepare bclData buffers to receive new tile data
    void resetBclData(
        const flowcell::TileMetadata& tileMetadata,
        alignment::BclClusters& bclData) const
    {
        // this implementation does nothing as chunk sizes are pre-determined
    }

    void loadClusters(
        const flowcell::TileMetadata &tileMetadata,
        alignment::BclClusters &bclData);

    unsigned getMaxTileClusters() const
    {
        return tileClustersMax_;
    }

};

template <>
struct DataSourceTraits<FastqBaseCallsSource>
{
    static const bool SUPPORTS_XY = false;
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
