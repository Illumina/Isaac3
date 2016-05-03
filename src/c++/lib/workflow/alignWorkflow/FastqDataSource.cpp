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
 ** \file FastqDataSource.cpp
 **
 ** \brief see FastqDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "workflow/alignWorkflow/FastqDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

FastqBaseCallsSource::FastqBaseCallsSource(
    const unsigned clustersAtATimeMax,
    const unsigned coresMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout &fastqFlowcellLayout,
    common::ThreadVector &threads) :
        tileClustersMax_(clustersAtATimeMax),
        coresMax_(coresMax),
        barcodeMetadataList_(barcodeMetadataList),
        fastqFlowcellLayout_(fastqFlowcellLayout),
        clusterLength_(flowcell::getTotalReadLength(fastqFlowcellLayout_.getReadMetadataList()) + fastqFlowcellLayout_.getBarcodeLength() + fastqFlowcellLayout_.getReadNameLength()),
        clusters_(clusterLength_),
        lanes_(fastqFlowcellLayout.getLaneIds()),
        currentLaneIterator_(lanes_.begin()),
        currentTile_(1),
        fastqLoader_(fastqFlowcellLayout_.getAttribute<flowcell::Layout::Fastq, flowcell::FastqVariableLengthOk>(), 0, threads, coresMax_)

{
}


flowcell::TileMetadataList FastqBaseCallsSource::discoverTiles()
{
    loadedTiles_.clear();

    // As we don't know whether the current lane has been completely loaded or we're in the
    // middle of discovering it's tiles, just attempt to load more data for it and stop only
    // when some data is loaded for this or subsequent lane or we run out of lanes.
    unsigned clustersLoaded = 0;
    while (!clustersLoaded)
    {
        if (lanes_.end() == currentLaneIterator_)
        {
            return loadedTiles_;
        }

        const unsigned clustersToLoad = tileClustersMax_;
        ISAAC_THREAD_CERR << "Resetting Fastq data for " << clustersToLoad << " clusters" << std::endl;
        clusters_.reset(clusterLength_, clustersToLoad);
        ISAAC_THREAD_CERR << "Resetting Fastq data done for " << clusters_.getClusterCount() << " clusters" << std::endl;
        // load clusters, return tile breakdown based on tileClustersMax_

        const boost::filesystem::path read1Path =
            fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
                *currentLaneIterator_, fastqFlowcellLayout_.getReadMetadataList().at(0).getNumber());
        if (1 == fastqFlowcellLayout_.getReadMetadataList().size())
        {
            // this will keep the current files open if the paths don't change
            fastqLoader_.open(read1Path, fastqFlowcellLayout_.getAttribute<flowcell::Layout::Fastq, flowcell::FastqBaseQ0>());
        }
        else // assume paired data
        {
            const boost::filesystem::path read2Path =
                fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq, flowcell::FastqFilePathAttributeTag>(
                    *currentLaneIterator_, fastqFlowcellLayout_.getReadMetadataList().at(1).getNumber());
            // this will keep the current files open if the paths don't change
            fastqLoader_.open(read1Path, read2Path, fastqFlowcellLayout_.getAttribute<flowcell::Layout::Fastq, flowcell::FastqBaseQ0>());
        }

        alignment::BclClusters::iterator clustersEnd = clusters_.cluster(0);
        clustersLoaded = fastqLoader_.loadClusters(clustersToLoad, fastqFlowcellLayout_.getReadNameLength(), fastqFlowcellLayout_.getReadMetadataList(), clustersEnd);
        ISAAC_THREAD_CERR << "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
        clusters_.reset(clusterLength_, clustersLoaded);

        if (!clustersLoaded)
        {
            // there was nothing left to load for this lane, proceed with the next one.
            ++currentLaneIterator_;
            currentTile_ = 1;
        }
    }

    const std::string &flowcellId = fastqFlowcellLayout_.getFlowcellId();
    while (clustersLoaded)
    {
        const unsigned clusterCount = std::min(clustersLoaded, tileClustersMax_);
        const flowcell::TileMetadata tileMetadata(
            flowcellId, fastqFlowcellLayout_.getIndex(),
            currentTile_++, *currentLaneIterator_,
            clusterCount,
            loadedTiles_.size());
        loadedTiles_.push_back(tileMetadata);

        if (clustersLoaded < tileClustersMax_)
        {
            break;
        }
        clustersLoaded -= tileClustersMax_;
    }

    ISAAC_ASSERT_MSG(1 == loadedTiles_.size(), "Expected to have only one tile loaded");
    return loadedTiles_;
}

void FastqBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_ASSERT_MSG(1 == loadedTiles_.size(), "Expected to have only one tile loaded");
    ISAAC_ASSERT_MSG(tileMetadata.getFlowcellIndex() == loadedTiles_.front().getFlowcellIndex(), "Unexpected tile requested");
    ISAAC_ASSERT_MSG(tileMetadata.getLane() == loadedTiles_.front().getLane(), "Unexpected tile lane requested: " << tileMetadata << " current: " << loadedTiles_.front());
    ISAAC_ASSERT_MSG(tileMetadata.getTile() == loadedTiles_.front().getTile(), "Unexpected tile tile requested: " << tileMetadata << " current: " << loadedTiles_.front());

    bclData.swap(clusters_);
}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
