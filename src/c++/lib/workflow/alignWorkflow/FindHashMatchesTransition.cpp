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
 ** \file FindHashMatchesTransition.cpp
 **
 ** \brief see FindHashMatchesTransition.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/ref.hpp>

#include "alignment/HashMatchFinder.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "alignment/matchSelector/DebugStorage.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/Numa.hh"
#include "common/ParallelSort.hpp"
#include "demultiplexing/DemultiplexingStatsXml.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "io/BitsetLoader.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/BclBgzfDataSource.hh"
#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/MultiTileDataSource.hh"
#include "workflow/alignWorkflow/FastqDataSource.hh"
#include "workflow/alignWorkflow/FindHashMatchesTransition.hh"



namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace findHashMatchesTransition
{

void wait(
    bool &signal,
    boost::condition_variable &condition,
    boost::unique_lock<boost::mutex> &lock,
    bool &forceTermination)
{
    if (forceTermination)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    while (signal)
    {
        if (forceTermination)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        condition.wait(lock);
    }
    signal = true;
}

void release(
    bool &signal,
    boost::condition_variable &condition,
    bool &forceTermination,
    const bool exceptionUnwinding)
{
    if (exceptionUnwinding)
    {
        forceTermination = true;
    }
    signal = false;
    condition.notify_all();
}

void Thread::binQscores(alignment::BclClusters &bclData) const
{
    ISAAC_THREAD_CERR << "Binning qscores" << std::endl;

    for(std::size_t clusterId = 0; bclData.getClusterCount() != clusterId; ++clusterId)
    {
        alignment::BclClusters::iterator clusterBegin = bclData.cluster(clusterId);
        for (const flowcell::ReadMetadata &readMetadata : flowcellLayout_.getReadMetadataList())
        {
            BclClusterFields::IteratorPair pair = bclFields_.getBcl(clusterBegin, readMetadata.getIndex());
            std::for_each(
                pair.first, pair.second, [this](char &bcl){bcl = fullBclQScoreTable_[static_cast<unsigned char>(bcl)];});
        }
    }
    ISAAC_THREAD_CERR << "Binning qscores done" << std::endl;
}

/**
 * \brief Finds matches for the lane. Updates foundMatches with match information and tile metadata identified during
 *        the processing.
 */
template <typename ReferenceHashT, typename DataSourceT>
void Thread::run(
    const flowcell::TileMetadataList &unprocessedTiles,
    unsigned &nextTile,
    unsigned &nextUnprocessedTile,
    DataSourceT &dataSource,
    alignment::SeedHashMatchFinder<ReferenceHashT> &matchFinder,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    common::ScopedMallocBlock &mallocBlock)
{
    boost::unique_lock<boost::mutex> lock(mutex_);
    while(unprocessedTiles.size() != nextTile)
    {
        const unsigned ourTile = nextTile++;
        const flowcell::TileMetadata &tileMetadata = unprocessedTiles.at(ourTile);

        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&release, boost::ref(loading_), boost::ref(stateChangedCondition_), boost::ref(forceTermination_), _1))
//        ISAAC_BLOCK_WITH_CLENAUP([this](bool){release(loading_, stateChangedCondition_);})
        {
            wait(loading_, stateChangedCondition_, lock, forceTermination_);
            {
                common::ScopedMallocBlockUnblock unblockMalloc(mallocBlock);
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);

                dataSource.resetBclData(tileMetadata, tileClusters_);
                dataSource.loadClusters(tileMetadata, tileClusters_);
                if(qScoreBin_)
                {
                    binQscores(tileClusters_);
                }
            }
        }

        ISAAC_BLOCK_WITH_CLENAUP([&](bool exceptionUnwinding)
                                 {
                                    if (exceptionUnwinding) {forceTermination_ = true;}
                                    stateChangedCondition_.notify_all();
                                 })
        {
            // make sure the order in which tiles are processed is same between different runs.
            while (nextUnprocessedTile != ourTile)
            {
                if (forceTermination_)
                {
                    BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
                }

                stateChangedCondition_.wait(lock);
            }

            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                matchSelector_.parallelSelect(tileClusterInfo, barcodeTemplateLengthStatistics, tileMetadata, matchFinder, tileClusters_, fragmentStorage_);
            }

            // swap the flush buffers while we still have compute lock
            wait(flushing_, stateChangedCondition_, lock, forceTermination_);
            {
                fragmentStorage_.prepareFlush();
            }
            ++nextUnprocessedTile;
        }

        // flush asynchronously so that other guys can load and align at the same time
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&release, boost::ref(flushing_), boost::ref(stateChangedCondition_), boost::ref(forceTermination_), _1))
//        ISAAC_BLOCK_WITH_CLENAUP([&](bool){release(flushing_, stateChangedCondition_);})
        {
            // flush slot already acquired when we had the compute slot but the state could have changed in between.
            if (forceTermination_)
            {
                BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
            }
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                fragmentStorage_.flush();
            }
        }
    }
}

} // namespace findHashMatchesTransition

FindHashMatchesTransition::FindHashMatchesTransition(
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool allowVariableFastqLength,
    const bool cleanupIntermediary,
    const unsigned bclTilesPerChunk,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    const uint64_t availableMemory,
    const unsigned clustersAtATimeMax,
    const bfs::path &tempDirectory,
    const bfs::path &demultiplexingStatsXmlPath,
    const unsigned int maxThreadCount,
    const unsigned seedBaseQualityMin,
    const unsigned repeatThreshold,
    const unsigned neighborhoodSizeThreshold,
    const bool ignoreNeighbors,
    const bool ignoreRepeats,
    const unsigned inputLoadersMax,
    const unsigned tempSaversMax,
    const common::ScopedMallocBlock::Mode memoryControl,
    const std::vector<std::size_t> &clusterIdList,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const reference::NumaContigLists &contigLists,
    const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations,
    const bool extractClusterXy,
    const int mateDriftRange,
    const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
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
    const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
    const bool qScoreBin,
    const boost::array<char, 256> &fullBclQScoreTable,
    const bool bufferBins,
    const unsigned expectedCoverage,
    const uint64_t targetBinSize,
    const double expectedBgzfCompressionRatio,
    const bool preSortBins,
    const bool preAllocateBins,
    const std::string &binRegexString
    )
    : flowcellLayoutList_(flowcellLayoutList)
    , tempDirectory_(tempDirectory)
    , demultiplexingStatsXmlPath_(demultiplexingStatsXmlPath)
    , coresMax_(maxThreadCount)
    , seedBaseQualityMin_(seedBaseQualityMin)
    , repeatThreshold_(repeatThreshold)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , barcodeMetadataList_(barcodeMetadataList)
    , allowVariableFastqLength_(allowVariableFastqLength)
    , cleanupIntermediary_(cleanupIntermediary)
    , bclTilesPerChunk_(bclTilesPerChunk)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , ignoreMissingFilters_(ignoreMissingFilters)
    , availableMemory_(availableMemory)
    , clustersAtATimeMax_(clustersAtATimeMax)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , memoryControl_(memoryControl)
    , clusterIdList_(clusterIdList)
    , sortedReferenceMetadataList_(sortedReferenceMetadataList)
    , extractClusterXy_(extractClusterXy)
    , bufferBins_(bufferBins)
    , expectedCoverage_(expectedCoverage)
    , targetBinSize_(targetBinSize)
    , expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio)
    , keepUnaligned_(keepUnaligned)
    , preSortBins_(preSortBins)
    , preAllocateBins_(preAllocateBins)
    , binRegexString_(binRegexString)

    // Have thread pool for the maximum number of threads we may potentially need.
    , threads_(std::max(inputLoadersMax_, coresMax_))

    , ioOverlapThreads_(bufferBins_ ? 3 : 2) //bin buffering needs an extra thread to flush the buffers
    , contigLists_(contigLists)
    , kUniquenessAnnotations_(kUniquenessAnnotations)

    , alignmentCfg_(gapMatchScore, gapMismatchScore, gapOpenScore, gapExtendScore, minGapExtendScore, splitGapLength)
    , matchSelector_(
        coresMax_,
        barcodeMetadataList_,
        flowcellLayoutList_,
        contigLists_,
        kUniquenessAnnotations_,
        repeatThreshold_,
        mateDriftRange,
        userTemplateLengthStatistics,
        mapqThreshold,
        perTileTls,
        pfOnly,
        collectCycleStats,
        baseQualityCutoff,
        keepUnaligned,
        clipSemialigned,
        clipOverlapping,
        scatterRepeats,
        rescueShadows,
        anchorMate,
        gappedMismatchesMax,
        smitWatermanGapsMax,
        smartSmithWaterman,
        noSmithWaterman,
        splitAlignments,
        gapMatchScore,
        gapMismatchScore,
        gapOpenScore,
        gapExtendScore,
        minGapExtendScore,
        splitGapLength,
        dodgyAlignmentScore,
        common::ScopedMallocBlock::Strict == memoryControl_),
        qScoreBin_(qScoreBin),
        fullBclQScoreTable_(fullBclQScoreTable)
{
}

static alignment::BinMetadataList buildBinPathList(
    const alignment::matchSelector::BinIndexMap &binIndexMap,
    const bfs::path &binDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const uint64_t expectedTotalReads,
    const bool preSortBins)
{
    ISAAC_THREAD_CERR << "expectedTotalReads " << expectedTotalReads << std::endl;
    ISAAC_TRACE_STAT("before buildBinPathList");
    alignment::BinMetadataList binPathList;
    ISAAC_ASSERT_MSG(!binIndexMap.empty(), "Empty binIndexMap is illegal");
    ISAAC_ASSERT_MSG(!binIndexMap.back().empty(), "Empty binIndexMap entry is illegal" << binIndexMap);
    binPathList.reserve(1 + binIndexMap.back().back());
    size_t contigIndex = 0;
    BOOST_FOREACH(const std::vector<unsigned> &contigBins, binIndexMap)
    {
        ISAAC_ASSERT_MSG(!contigBins.empty(), "Unexpected empty contigBins");
        // matchDistribution contig 0 is the first contig
        // binIndexMap contig 0  is unaligned bin
        for (unsigned i = contigBins.front(); contigBins.back() >= i; ++i)
        {
            ISAAC_ASSERT_MSG(binPathList.size() == i, "Basic sanity checking for bin numbering failed");
            const reference::ReferencePosition binStartPos = binIndexMap.getBinFirstPos(i);
            ISAAC_ASSERT_MSG(!i || binIndexMap.getBinIndex(binStartPos) == i, "BinIndexMap is broken");
            binPathList.push_back(
                alignment::BinMetadata(
                    barcodeMetadataList.size(),
                    binPathList.size(),
                    binStartPos,
                    // bin zero has length of totalReads as it contains unaligned records which are chunked by the number of reads stored
                    i ? binIndexMap.getBinFirstInvalidPos(i) - binStartPos : expectedTotalReads,
                        // Pad file names well, so that we don't have to worry about them becoming of different length.
                        // This is important for memory reservation to be stable
                    binDirectory / (boost::format("bin-%08d-%08d.dat") % contigIndex % i).str(),
                    // Normally, aim to have 1024 or less chunks.
                    // This will require about 4096*1024 (4 megabytes) of cache when pre-sorting bin during the loading in bam generator
                    preSortBins ? 1024 : 0));
        }
        ++contigIndex;
    }
    return binPathList;
}

template <typename KmerT>
void FindHashMatchesTransition::perform(
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath)
{
    align<KmerT>(foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath);
}

void FindHashMatchesTransition::resolveBarcodes(
    const flowcell::Layout &flowcell,
    const flowcell::BarcodeMetadataList &barcodeGroup,
    BarcodeSource &barcodeSource,
    // this contains tiles we are processing but they are not placed at the tile.getIndex()
    flowcell::TileMetadataList unprocessedTiles,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    demultiplexing::DemultiplexingStats &demultiplexingStats)
{
    ISAAC_ASSERT_MSG(!barcodeGroup.empty(), "At least 'none' barcode must be defined");
    if (1 == barcodeGroup.size())
    {
        const flowcell::BarcodeMetadata barcode = barcodeGroup.at(0);
        ISAAC_ASSERT_MSG(barcode.isNoIndex(), "If barcode group has only one entry it must be the 'NoIndex' barcode");
        BOOST_FOREACH(const flowcell::TileMetadata &tile, unprocessedTiles)
        {
            for (unsigned clusterId = 0; clusterId < tile.getClusterCount(); ++clusterId)
            {
                tileClusterInfo.setBarcodeIndex(tile.getIndex(), clusterId, barcode.getIndex());
                demultiplexingStats.recordBarcode(demultiplexing::BarcodeId(tile.getIndex(), barcode.getIndex(), clusterId, 0));
            }
            ISAAC_THREAD_CERR << "Forced barcode index for clusters of " << tile << " to " << barcode << std::endl;
        }
    }
    else
    {
        demultiplexing::BarcodeResolver barcodeResolver(barcodeMetadataList_, barcodeGroup);

        flowcell::TileMetadataList currentTiles; currentTiles.reserve(unprocessedTiles.size());

        while (!unprocessedTiles.empty())
        {
            currentTiles.clear();
            if (!demultiplexing::BarcodeMemoryManager::selectTiles(unprocessedTiles, currentTiles))
            {
                BOOST_THROW_EXCEPTION(common::MemoryException("Insufficient memory to load barcodes even for just one tile: " +
                    boost::lexical_cast<std::string>(unprocessedTiles.back())));
            }

            demultiplexing::Barcodes barcodes;
            // this will take at most the same amount of ram as a set of singleseeds
            ISAAC_ASSERT_MSG(barcodeGroup.size(), "Barcode list must be not empty");
            ISAAC_ASSERT_MSG(barcodeGroup.at(0).isDefault(), "The very first barcode must be the 'unknown indexes or no index' one");
            barcodeSource.loadBarcodes(flowcell, barcodeGroup.at(0).getIndex(), currentTiles, barcodes);
            barcodeResolver.resolve(barcodes, demultiplexingStats);

            BOOST_FOREACH(const demultiplexing::Barcode &barcode, barcodes)
            {
                tileClusterInfo.setBarcodeIndex(barcode.getTile(), barcode.getCluster(), barcode.getBarcode());
            }
        }
    }
}

template <typename ReferenceHashT>
ReferenceHashT buildReferenceHash(
    const reference::SortedReferenceMetadata &sortedReferenceMetadata,
    const reference::ContigList &contigList,
    common::ThreadVector &threads,
    const unsigned coresMax)
{
    reference::ReferenceHasher<ReferenceHashT> hasher(sortedReferenceMetadata, contigList, threads, coresMax);

    return hasher.generate();
}

/**
 * \brief Finds matches for the lane. Updates foundMatches with match information and tile metadata identified during
 *        the processing.
 */
template <typename ReferenceHashT, typename DataSourceT>
void FindHashMatchesTransition::findLaneMatches(
    const ReferenceHashT &referenceHash,
    const flowcell::Layout &flowcell,
    const unsigned lane,
    const flowcell::BarcodeMetadataList &laneBarcodes,
    const flowcell::TileMetadataList &unprocessedTiles,
    DataSourceT &dataSource,
    const unsigned maxTileClusters,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    alignment::matchSelector::FragmentStorage &fragmentStorage)
{
    if (!unprocessedTiles.empty())
    {
        alignment::matchFinder::TileClusterInfo tileClusterInfo(unprocessedTiles);
        ISAAC_THREAD_CERR << "Resolving barcodes for " << flowcell << " lane " << lane << std::endl;
        resolveBarcodes(
            flowcell, laneBarcodes,
            dataSource,
            unprocessedTiles, tileClusterInfo, demultiplexingStats);
        ISAAC_THREAD_CERR << "Resolving barcodes done for " << flowcell << " lane " << lane << std::endl;

        ISAAC_THREAD_CERR << "Finding hash matches for " << flowcell.getSeedMetadataList() << "with repeat threshold: " << repeatThreshold_ << std::endl;

        alignment::SeedHashMatchFinder<ReferenceHashT> matchFinder(referenceHash, seedBaseQualityMin_, repeatThreshold_);

        matchSelector_.reserveMemory(unprocessedTiles);

        boost::mutex mutex;
        boost::condition_variable stateChangedCondition;
        bool loading(false);
        bool flushing(false);
        bool forceTermination(false);

        std::vector<findHashMatchesTransition::Thread> threads(
            ioOverlapThreads_.size(),
            findHashMatchesTransition::Thread(
                mutex,
                stateChangedCondition,
                loading,
                flushing,
                forceTermination,
                flowcell,
                maxTileClusters,
                DataSourceTraits<DataSourceT>::SUPPORTS_XY && extractClusterXy_,
                qScoreBin_,
                fullBclQScoreTable_,
                matchSelector_,
                fragmentStorage));

        {
            unsigned current = 0;
            unsigned nextUnprocessed = 0;
            common::ScopedMallocBlock  mallocBlock(memoryControl_);
            ioOverlapThreads_.execute
            (
                [&](const unsigned threadNumber, const unsigned threadsTotal)
                {
                    threads.at(threadNumber).run(
                        unprocessedTiles, current, nextUnprocessed, dataSource, matchFinder, tileClusterInfo,
                        barcodeTemplateLengthStatistics, mallocBlock);
                }
            );
        }

        ISAAC_THREAD_CERR << "Finding Single-seed matches done for " << flowcell.getSeedMetadataList() << std::endl;
    }
}

/**
 * \brief Collect all barcodes belonging to the flowcell lane. Default barcode is placed at the beginning of the result list
 *
 * \return Returns subset of barcodeMetadataList or an empty list if none of the barcodes are mapped to a reference
 */
static flowcell::BarcodeMetadataList findFlowcellLaneBarcodes(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout& flowcell,
    const unsigned lane)
{
    bool allBarcodesUnmapped = true;
    // put a placeholder for the unknown barcode in the beginning of the list
    flowcell::BarcodeMetadataList ret(1);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        if (barcode.getFlowcellId() == flowcell.getFlowcellId() &&
            barcode.getLane() == lane)
        {
            ISAAC_THREAD_CERR << "adding " << barcode << std::endl;
            if (barcode.isDefault())
            {
                ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX == ret.at(0).getIndex(), "More than one explicit specification for 'default' barcode within the group.");
                ret[0] = barcode;
            }
            else
            {
                ret.push_back(barcode);
            }
            allBarcodesUnmapped &= barcode.isUnmappedReference();
        }
    }

    ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX != ret.at(0).getIndex(), "Missing default barcode specification");

    if (allBarcodesUnmapped)
    {
        ret.clear();
    }
    return ret;
}


template <typename ReferenceHashT, typename DataSourceT>
void FindHashMatchesTransition::processFlowcellTiles(
    const ReferenceHashT &referenceHash,
    const flowcell::Layout& flowcell,
    DataSourceT &dataSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    FoundMatchesMetadata &foundMatches,
    alignment::matchSelector::FragmentStorage &fragmentStorage)
{
    fragmentStorage.reserve(dataSource.getMaxTileClusters());

    for (flowcell::TileMetadataList laneTiles = dataSource.discoverTiles(); !laneTiles.empty();
        laneTiles = dataSource.discoverTiles())
    {
        const unsigned lane = laneTiles[0].getLane();
        flowcell::BarcodeMetadataList laneBarcodes = findFlowcellLaneBarcodes(barcodeMetadataList_, flowcell, lane);
        if (laneBarcodes.empty())
        {
            ISAAC_THREAD_CERR << "Skipping flowcell " << flowcell.getFlowcellId() << " lane " << lane << " as none of the barcodes map to the reference" << std::endl;
        }
        else
        {
            BOOST_FOREACH(TileMetadata &tileMetadata, laneTiles)
            {
                foundMatches.addTile(tileMetadata);
                // this fixes the tile index to be correct in the context of the global tile list.
                // it is important for findLaneMatches to use the global tile index so that it can store
                // statistics properly
                tileMetadata = foundMatches.tileMetadataList_.back();
            }
            findLaneMatches(
                referenceHash, flowcell, lane, laneBarcodes, laneTiles, dataSource, dataSource.getMaxTileClusters(),
                demultiplexingStats, barcodeTemplateLengthStatistics, fragmentStorage);
        }
    }
}

template <typename ReferenceHashT>
void FindHashMatchesTransition::alignFlowcells(
    const ReferenceHashT &referenceHash,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &foundMatches,
    alignment::matchSelector::FragmentStorage &fragmentStorage)
{

    BOOST_FOREACH(const flowcell::Layout& flowcell, flowcellLayoutList_)
    {
        switch (flowcell.getFormat())
        {
            case flowcell::Layout::Bam:
            {
                BamBaseCallsSource dataSource(
                    tempDirectory_,
                    availableMemory_,
                    clustersAtATimeMax_,
                    cleanupIntermediary_,
                    // the loading itself occurs on one thread at a time only. So, the real limit is to avoid using
                    // more cores for decompression than the system actually has.
                    // On the other hand, there might be a need to limit the io to 1 thread, while allowing
                    // for the multithreaded processing of other cpu-demanding things.
                    std::min(inputLoadersMax_, coresMax_),
                    flowcell, threads_);
                processFlowcellTiles(referenceHash, flowcell, dataSource, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::Fastq:
            {
                FastqBaseCallsSource dataSource(
                    clustersAtATimeMax_,
                    allowVariableFastqLength_,
                    coresMax_,
                    barcodeMetadataList_,
                    flowcell,
                    threads_);

                processFlowcellTiles(referenceHash, flowcell, dataSource, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::Bcl:
            {
                BclBaseCallsSource baseCalls(
                    flowcell, ignoreMissingBcls_, ignoreMissingFilters_, threads_, inputLoadersMax_, extractClusterXy_);

                MultiTileBaseCallsSource<BclBaseCallsSource> multitileBaseCalls(bclTilesPerChunk_, flowcell, baseCalls);

                processFlowcellTiles(
                    referenceHash, flowcell, multitileBaseCalls, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::BclBgzf:
            {
                BclBgzfBaseCallsSource baseCalls(
                    flowcell, ignoreMissingBcls_, ignoreMissingFilters_, threads_, inputLoadersMax_, extractClusterXy_);
                MultiTileBaseCallsSource<BclBgzfBaseCallsSource> multitileBaseCalls(
                    bclTilesPerChunk_, flowcell, baseCalls);

                processFlowcellTiles(referenceHash, flowcell, multitileBaseCalls, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Unexpected flowcell format " << flowcell.getFormat());
                break;
            }
        }
    }
}
template <typename ReferenceHashT>
void FindHashMatchesTransition::alignFlowcells(
    const ReferenceHashT &referenceHash,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &ret)
{

    // asume most fragments will have a one-component CIGAR.
    const unsigned int estimatedFragmentSize = io::FragmentHeader::getMinTotalLength(
        flowcell::getMaxReadLength(flowcellLayoutList_),
        flowcell::getMaxClusterName(flowcellLayoutList_));

    const uint64_t fragmentsPerBin = targetBinSize_
        ? targetBinSize_ / estimatedFragmentSize
        : build::Build::estimateOptimumFragmentsPerBin(
            estimatedFragmentSize, availableMemory_, expectedBgzfCompressionRatio_,
            coresMax_);

    const uint64_t expectedBinSize = targetBinSize_? targetBinSize_ : fragmentsPerBin * estimatedFragmentSize;

    ISAAC_TRACE_STAT("AlignWorkflow::selectMatches ")

    alignment::EstimatedMatchDistribution matchDistribution(
        expectedCoverage_, flowcell::getMaxReadLength(flowcellLayoutList_), sortedReferenceMetadataList_);
    alignment::matchSelector::BinIndexMap binIndexMap(
        matchDistribution, fragmentsPerBin, "skip-empty" == binRegexString_);

    binMetadataList =
        buildBinPathList(binIndexMap, tempDirectory_, barcodeMetadataList_,
                         binIndexMap.getTotalBins() * fragmentsPerBin,
                         preSortBins_);

    ISAAC_THREAD_CERR << "Selecting matches using " << fragmentsPerBin << " fragments per bin limit. expectedBinSize: " << expectedBinSize << " bytes" << std::endl;


    std::unique_ptr<alignment::matchSelector::FragmentStorage> storagePtr(!bufferBins_ ?
        static_cast<alignment::matchSelector::FragmentStorage*>(new alignment::matchSelector::BinningFragmentStorage(
                            keepUnaligned_,
                            binIndexMap, binMetadataList,
                            preAllocateBins_ ? expectedBinSize : 0)) :
        static_cast<alignment::matchSelector::FragmentStorage*>(new alignment::matchSelector::BufferingFragmentStorage(
                        keepUnaligned_, coresMax_, tempSaversMax_,
                        binIndexMap, binMetadataList,
                        flowcellLayoutList_,
                        preAllocateBins_ ? expectedBinSize : 0)));

#ifdef ISAAC_DEV_STATS_ENABLED
        alignment::matchSelector::DebugStorage debugStorage(
            contigLists_.node0Container(), kUniquenessAnnotations_.node0Container(),
            alignmentCfg_, flowcellLayoutList_, demultiplexingStatsXmlPath_.parent_path(), *storagePtr);
        alignFlowcells(referenceHash, barcodeTemplateLengthStatistics, demultiplexingStats, ret, debugStorage);
#else
        alignFlowcells(referenceHash, barcodeTemplateLengthStatistics, demultiplexingStats, ret, *storagePtr);
#endif

}

template <typename KmerT>
void FindHashMatchesTransition::align(
    FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath)
{

    typedef reference::ReferenceHash<KmerT, common::NumaAllocator<void, 0> > ReferenceHashT;
    const reference::NumaReferenceHash<ReferenceHashT> referenceHash(
        buildReferenceHash<ReferenceHashT>(
            sortedReferenceMetadataList_.front(), contigLists_.node0Container().front(), threads_, coresMax_));

    FoundMatchesMetadata ret(tempDirectory_, barcodeMetadataList_, 1, sortedReferenceMetadataList_);
    demultiplexing::DemultiplexingStats demultiplexingStats(flowcellLayoutList_, barcodeMetadataList_);

    alignFlowcells(referenceHash, binMetadataList, barcodeTemplateLengthStatistics, demultiplexingStats, ret);

    dumpStats(demultiplexingStats, ret.tileMetadataList_);
    foundMatches.swap(ret);

    matchSelector_.unreserve();

    matchSelector_.dumpStats(matchSelectorStatsXmlPath);
}

void FindHashMatchesTransition::dumpStats(
    const demultiplexing::DemultiplexingStats &demultiplexingStats,
    const flowcell::TileMetadataList &tileMetadataList) const
{
    demultiplexing::DemultiplexingStatsXml statsXml;

    const unsigned maxLaneNumber = flowcell::getMaxLaneNumber(flowcellLayoutList_);
    for (unsigned lane = 1; lane <= maxLaneNumber; ++lane)
    {
        typedef std::map<std::string, demultiplexing::LaneBarcodeStats> SampleLaneBarcodeStats;
        typedef std::map<std::string, SampleLaneBarcodeStats> ProjectSampleLaneBarcodeStats;
        typedef std::map<std::string, ProjectSampleLaneBarcodeStats> FlowcellProjectSampleLaneBarcodeStats;
        FlowcellProjectSampleLaneBarcodeStats flowcellProjectSampleStats;
        BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
        {
            if (barcode.getLane() == lane)
            {
                // put one lane stat for each unknown barcode found.
                if (barcode.isUnknown())
                {
                    const flowcell::Layout& flowcell = flowcellLayoutList_.at(barcode.getFlowcellIndex());
                    statsXml.addFlowcellLane(flowcell, lane,
                                             demultiplexingStats.getLaneUnknwonBarcodeStat(barcode.getIndex()));
                }
                const demultiplexing::LaneBarcodeStats &stat = demultiplexingStats.getLaneBarcodeStat(barcode);
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()]["all"]["all"] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats["all"]["all"]["all"] += stat;
                statsXml.addLaneBarcode(barcode.getFlowcellId(), barcode.getProject(), barcode.getSampleName(), barcode.getName(), lane, stat);
            }
        }
        BOOST_FOREACH(const FlowcellProjectSampleLaneBarcodeStats::value_type &flowcellStats, flowcellProjectSampleStats)
        {
            BOOST_FOREACH(const ProjectSampleLaneBarcodeStats::value_type &projectStats, flowcellStats.second)
            {
                BOOST_FOREACH(const SampleLaneBarcodeStats::value_type &sampleStats, projectStats.second)
                {
                    statsXml.addLaneBarcode(flowcellStats.first, projectStats.first, sampleStats.first, "all", lane, sampleStats.second);
                }
            }
        }
    }

    std::ofstream os(demultiplexingStatsXmlPath_.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + demultiplexingStatsXmlPath_.string()));
    }
    if (!(os << statsXml)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: failed to store MatchFinder statistics in : " + demultiplexingStatsXmlPath_.string()));
    }
}

template void FindHashMatchesTransition::perform<oligo::ShortKmerType>(
    FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath);


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
