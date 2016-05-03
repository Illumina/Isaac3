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
 ** \file BamDataSource.cpp
 **
 ** \brief see BamDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "oligo/Nucleotides.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/bamDataSource/SingleEndClusterExtractor.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{


/**
 * \param clusterCount  Maximum number of clusters to load
 * \param clusterIt     Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of clusterLength
 * \param pfIt          Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of pf flags
 *
 * \return Actual number of loaded clusters
 */
template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadPairedReads(
    unsigned clusterCount, const unsigned nameLengthMax,
    const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;

    while (clusterCount)
    {
        if(clusterExtractor_.extractingUnpaired())
        {
            ISAAC_THREAD_CERR << "extracting unpaired " << std::endl;
            clusterCount = clusterExtractor_.extractUnpaired(
                readMetadataList.at(0).getLength(), readMetadataList.at(1).getLength(), nameLengthMax, clusterCount,
                clusterIt, pfIt);
            // Either no room in result buffer or no more data available.
            break;
        }
        else
        {
            if (!clusterExtractor_.isEmpty())
            {
                ISAAC_THREAD_CERR << "resuming from " << clusterExtractor_.size() << " pending elements" << std::endl;
                clusterCount = clusterExtractor_.extractPairedReads(
                    nameLengthMax, clusterCount, clusterIt, pfIt, readMetadataList);

                if (!clusterCount)
                {
                    // there is no room to extract any more pending items. Resume on next call.
                    return requestedClusterCount;
                }
            }

            bamLoader_.load
            (
                boost::make_tuple(
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::append<ClusterInsertIt, PfInserIt>,
                                &clusterExtractor_, _1, _2, nameLengthMax, boost::ref(clusterCount),
                                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::removeOld,
                                &clusterExtractor_, _1, _2, boost::ref(readMetadataList)))
            );


            if (clusterCount)
            {
                // we've ran out of data in the bam file. See if unpaired items can be paired
                clusterExtractor_.startExtractingUnpaired();
            }
        }
    }

    ISAAC_THREAD_CERR << "loadPairedReads done clusterCount:" << clusterCount << std::endl;
    return requestedClusterCount - clusterCount;
}

void BamClusterLoader::open(
    const std::string &flowcellId,
    const boost::filesystem::path &bamPath)
{
    if (flowcellId_ != flowcellId)
    {
        clusterExtractor_.open(flowcellId);
        bamLoader_.open(bamPath);
        flowcellId_ = flowcellId;
    }
    else
    {
        ISAAC_THREAD_CERR << "Keeping bam stream open for flowcellId " << flowcellId_ << " " << bamPath << std::endl;
    }
}

template <typename InsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadSingleReads(
    unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList,
    InsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;
    bamDataSource::SingleEndClusterExtractor extractor;
    bamLoader_.load
    (
        boost::make_tuple(
            boost::bind(&bamDataSource::SingleEndClusterExtractor::extractSingleRead<InsertIt, PfInserIt>,
                &extractor, _1, nameLengthMax, boost::ref(clusterCount),
                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
            boost::bind(&bamDataSource::SingleEndClusterExtractor::nothing))
    );

    ISAAC_THREAD_CERR << "loadSingleReads done clusterCount:" << clusterCount << " bgzfReader_.isEof():" << std::endl;

    return requestedClusterCount - clusterCount;}

template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadClusters(
    unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    return 2 == readMetadataList.size() ?
        loadPairedReads(clusterCount, nameLengthMax, readMetadataList, clusterIt, pfIt) :
        loadSingleReads(clusterCount, nameLengthMax, readMetadataList, clusterIt, pfIt);
}

inline std::size_t getBamFileSize(const flowcell::Layout &flowcell)
{
    return common::getFileSize(flowcell.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>().c_str());
}

BamBaseCallsSource::BamBaseCallsSource(
    const boost::filesystem::path &tempDirectoryPath,
    const uint64_t availableMemory,
    const unsigned clustersAtATimeMax,
    const bool cleanupIntermediary,
    const unsigned coresMax,
    const flowcell::Layout &bamFlowcellLayout,
    common::ThreadVector &threads) :
        bamFlowcellLayout_(bamFlowcellLayout),
        tileClustersMax_(clustersAtATimeMax),
        coresMax_(coresMax),
        clusterLength_(flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList()) + bamFlowcellLayout_.getBarcodeLength() + bamFlowcellLayout_.getReadNameLength()),
        clusters_(clusterLength_),
        currentTile_(0),
        threads_(threads),
        bamClusterLoader_(
            cleanupIntermediary, 0, threads, coresMax, tempDirectoryPath,
            getBamFileSize(bamFlowcellLayout_), bamFlowcellLayout.getFlowcellId().length(),
            flowcell::getTotalReadLength(bamFlowcellLayout.getReadMetadataList()),
            flowcell::getMinReadLength(bamFlowcellLayout.getReadMetadataList()))
{
}

flowcell::TileMetadataList BamBaseCallsSource::discoverTiles()
{
    flowcell::TileMetadataList ret;

    const unsigned clustersToLoad = tileClustersMax_;
//    const unsigned clustersToLoad = clustersAtATimeMax_;
    clusters_.reset(clusterLength_, clustersToLoad);
    // load clusters, return tile breakdown based on tileClustersMax_

    const boost::filesystem::path bamPath = bamFlowcellLayout_.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>();
    // this will keep the current files open if the paths don't change
    bamClusterLoader_.open(bamFlowcellLayout_.getFlowcellId(), bamPath);

    alignment::BclClusters::iterator clustersEnd = clusters_.cluster(0);
    clusters_.pf().clear();
    std::back_insert_iterator<std::vector<bool> > pfIt(clusters_.pf());
    const unsigned clustersLoaded = bamClusterLoader_.loadClusters(
        clustersToLoad, bamFlowcellLayout_.getReadNameLength(), bamFlowcellLayout_.getReadMetadataList(), clustersEnd, pfIt);
    ISAAC_THREAD_CERR << "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
    clusters_.reset(clusterLength_, clustersLoaded);
    if (clustersLoaded)
    {
        const std::string &flowcellId = bamFlowcellLayout_.getFlowcellId();

        const flowcell::TileMetadata tileMetadata(
            flowcellId, bamFlowcellLayout_.getIndex(),
            currentTile_ + 1, 1,
            clustersLoaded,
            currentTile_);
        ++currentTile_;
        ret.push_back(tileMetadata);
        ISAAC_THREAD_CERR << "Generated bam tile: " <<
            oligo::bclToString(reinterpret_cast<const unsigned char *>(&*clusters_.cluster(0)), flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList())) << " " <<
            oligo::bclToString(reinterpret_cast<const unsigned char *>(&*clusters_.cluster(0)) +
                               flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList()) * (tileMetadata.getClusterCount() - 1),
                               flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList())) << " " << tileMetadata << std::endl;
    }

    return ret;
}

void BamBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loaded bam tile: " << tileMetadata << " with " << clusters_.getClusterCount() << " clusters" << std::endl;

    bclData.swap(clusters_);
}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
