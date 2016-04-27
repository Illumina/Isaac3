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
 ** \file BamDataSource.hh
 **
 ** \brief Encapsulation of single-ended and paired data stored in bam file
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_BAM_DATA_SOURCE_HH

#include "alignment/BclClusters.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/BamLayout.hh"
#include "flowcell/TileMetadata.hh"
#include "io/BamLoader.hh"
#include "workflow/alignWorkflow/bamDataSource/PairedEndClusterExtractor.hh"
#include "workflow/alignWorkflow/DataSource.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

class BamClusterLoader
{
    std::string flowcellId_;

    io::BamLoader bamLoader_;
    bamDataSource::PairedEndClusterExtractor clusterExtractor_;
public:
    BamClusterLoader(
        const bool cleanupIntermediary,
        const std::size_t maxPathLength,
        common::ThreadVector &threads,
        const unsigned coresMax,
        const boost::filesystem::path &tempDirectoryPath,
        const std::size_t maxBamFileLength,
        const std::size_t maxFlowcellIdLength,
        const std::size_t minClusterLength,
        const std::size_t minReadLength) :
        bamLoader_(maxPathLength, threads, coresMax),
        clusterExtractor_(tempDirectoryPath, maxBamFileLength, maxFlowcellIdLength, minClusterLength, cleanupIntermediary,
                          // assume each uncompressed bam record is roughly sizeof(header) + (read length * 2). Double the estimate.
                          bamLoader_.BUFFER_SIZE / (sizeof(bam::BamBlockHeader) + minReadLength * 2) * 2)
    {
        flowcellId_.reserve(maxFlowcellIdLength);
    }

    void open(
        const std::string &flowcellId,
        const boost::filesystem::path &bamPath);

    template <typename ClusterInsertIt, typename PfInserIt>
    unsigned loadClusters(unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList,
        ClusterInsertIt &clusterIt, PfInserIt &pfIt);

private:
    // pessimistic assumption that each bam record is made of the header only.
    template <typename InsertIt, typename PfInserIt>
    unsigned loadPairedReads(
        unsigned clusterCount, const unsigned nameLengthMax,
        const flowcell::ReadMetadataList &readMetadataList, InsertIt &it, PfInserIt &pfIt);
    template <typename InsertIt, typename PfInserIt>
    unsigned loadSingleReads(
        unsigned clusterCount, const unsigned nameLengthMax,
        const flowcell::ReadMetadataList &readMetadataList, InsertIt &it, PfInserIt &pfIt);
};


class BamBaseCallsSource : public TileSource, public BarcodeSource
{
    const flowcell::Layout &bamFlowcellLayout_;
    const unsigned tileClustersMax_;
    const unsigned coresMax_;
    const unsigned clusterLength_;
    alignment::BclClusters clusters_;
    unsigned currentTile_;
    common::ThreadVector &threads_;
    BamClusterLoader bamClusterLoader_;

public:
    BamBaseCallsSource(
        const boost::filesystem::path &tempDirectoryPath,
        const uint64_t availableMemory,
        const unsigned clustersAtATimeMax,
        const bool cleanupIntermediary,
        const unsigned coresMax,
        const flowcell::Layout &fastqFlowcellLayout,
        common::ThreadVector &threads);

    // TileSource implementation
    flowcell::TileMetadataList discoverTiles();

    // BarcodeSource implementation
    virtual void loadBarcodes(
        const flowcell::Layout &flowcell,
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles,
        std::vector<demultiplexing::Barcode> &barcodes)
    {
        ISAAC_ASSERT_MSG(false, "Barcode resolution is not implemented for Bam data");
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

private:
    static unsigned determineMemoryCapacity(
        const uint64_t availableMemory,
        const unsigned tileClustersMax,
        const std::size_t seedsPerCluster,
        const unsigned clusterLength);
};

template <>
struct DataSourceTraits<BamBaseCallsSource>
{
    static const bool SUPPORTS_XY = false;
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FASTQ_DATA_SOURCE_HH
