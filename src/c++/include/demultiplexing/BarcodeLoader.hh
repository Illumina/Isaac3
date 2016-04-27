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
 ** \file BarcodeLoader.hh
 **
 ** Helper class for loading barcode data from Bcl files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
#define iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH

#include <boost/iterator/counting_iterator.hpp>
#include <boost/mpl/equal_to.hpp>

#include "common/Debug.hh"
#include "common/Threads.hpp"

#include "demultiplexing/Barcode.hh"

#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"

#include "io/FileBufCache.hh"
#include "rta/BclMapper.hh"

namespace isaac
{
namespace demultiplexing
{


class BarcodeMemoryManager: boost::noncopyable
{
public:
    /// Determine how many tiles can have their barcoded loaded at the same time
    static bool selectTiles(
        flowcell::TileMetadataList &unprocessedPool,
        flowcell::TileMetadataList &selectedTiles)
    {
        selectedTiles.swap(unprocessedPool);
        {
            ISAAC_THREAD_CERR << "Barcode resolution: Determining the number of tiles that can be processed simultaneously..." << std::endl;
            while(!selectedTiles.empty() && !seeIfFits(selectedTiles))
            {
                unprocessedPool.push_back(selectedTiles.back());
                selectedTiles.pop_back();
            }
            if (selectedTiles.empty())
            {
                return false;
            }

            ISAAC_THREAD_CERR << "Barcode resolution: Determining the number of tiles that can be processed simultaneously done." << std::endl;

            if (!unprocessedPool.empty())
            {
                ISAAC_THREAD_CERR << "WARNING: will resolve barcodes in parts due to the memory limit. "
                    "This pass will process only " << selectedTiles.size() << " tiles" << std::endl;
            }
        }
        return true;
    }

    static void allocate(const flowcell::TileMetadataList &tiles, Barcodes &barcodes)
    {
        const uint64_t totalClusterCount = getTotalBarcodeCount(tiles);
        ISAAC_THREAD_CERR << "Allocating barcode storage for " << totalClusterCount << " barcodes" << std::endl;

        barcodes.clear();
        barcodes.resize(totalClusterCount);

        ISAAC_THREAD_CERR << "Allocating barcode storage done for " << totalClusterCount << " seeds" << std::endl;
    }

private:
    static bool seeIfFits(const flowcell::TileMetadataList &tiles)
    {
        try
        {
            Barcodes test;
            test.reserve(getTotalBarcodeCount(tiles));
            return true;
        }
        catch (std::bad_alloc &e)
        {
            // reset errno, to prevent misleading error messages when failing code does not set errno
            errno = 0;
        }
        return false;
    }

    static unsigned getTotalBarcodeCount(const flowcell::TileMetadataList &tiles)
    {
        const uint64_t totalClusterCount = std::accumulate(tiles.begin(), tiles.end(), 0,
                                                             boost::bind(std::plus<unsigned>(), _1,
                                                                         boost::bind(&flowcell::TileMetadata::getClusterCount, _2)));
        return totalClusterCount;
    }
};


/**
 ** \brief Encapsulates the variables that are shared by all the threads while
 ** loading the barcodes.
 **/
template <typename ReaderT>
class ParallelBarcodeLoader : boost::noncopyable
{
public:
    /**
     ** \brief constructs an instance with all the required shorthands.
     **
     ** Note: all parameters are kept as references and it is the
     ** responsibility of the caller to ensure appropriate life time for the
     ** referenced variables.
     **/
    ParallelBarcodeLoader(
        const unsigned maxClustersPerTile,
        std::vector<ReaderT> &threadReaders)
            : maxClustersPerTile_(maxClustersPerTile)
            , threadBclReaders_(threadReaders)
    {

        while(threadBclMappers_.size() < threadBclReaders_.size())
        {
            threadBclMappers_.push_back(
                new rta::SingleCycleBclMapper<ReaderT>(
                    maxClustersPerTile_,
                    threadBclReaders_.at(threadBclMappers_.size())));
        }
    }


    void load(const unsigned unknownBarcodeIndex,
              Barcodes &barcodes,
              const flowcell::Layout &flowcellLayout,
              Barcodes::iterator &nextTileBarcodes,
              flowcell::TileMetadataList::const_iterator &nextTile,
              const flowcell::TileMetadataList::const_iterator tilesEnd,
              const unsigned threadNumber)
    {
        ISAAC_ASSERT_MSG(
            flowcell::Layout::Bcl == flowcellLayout.getFormat() ||
            flowcell::Layout::BclBgzf == flowcellLayout.getFormat(), "Only bcl barcode loading is supported");
        ISAAC_ASSERT_MSG(MAX_BARCODE_LENGTH >= flowcellLayout.getBarcodeLength(), "barcode cannot be longer than " << MAX_BARCODE_LENGTH << " bases");

        boost::lock_guard<boost::mutex> lock(mutex_);
        while (tilesEnd != nextTile)
        {
            flowcell::TileMetadataList::const_iterator currentTile = nextTile++;
            Barcodes::iterator destinationBegin = nextTileBarcodes;
            nextTileBarcodes += currentTile->getClusterCount();

            const Barcodes::const_iterator  destinationEnd = nextTileBarcodes;
            ISAAC_ASSERT_MSG(destinationEnd <= barcodes.end(), "Computed end is past the end of the reserved buffer");

            {
                common::unlock_guard<boost::mutex> unlock(mutex_);

                ISAAC_THREAD_CERR << "Formatting tile barcodes for " << *currentTile << std::endl;
                std::transform(boost::counting_iterator<unsigned>(0),
                               boost::counting_iterator<unsigned>(currentTile->getClusterCount()),
                               destinationBegin,
                               boost::bind(&Barcode::constructFromTileBarcodeCluster, currentTile->getIndex(), unknownBarcodeIndex, _1));
                ISAAC_THREAD_CERR << "Formatting tile barcodes done for " << *currentTile << std::endl;

                ISAAC_THREAD_CERR << "Loading tile barcodes for " << *currentTile << std::endl;
                BOOST_FOREACH(const unsigned cycle, flowcellLayout.getBarcodeCycles())
                {
                    loadTileCycle(threadBclMappers_.at(threadNumber), destinationBegin, flowcellLayout, *currentTile, cycle);
                }
                ISAAC_THREAD_CERR << "Loading tile barcodes done for " << *currentTile << std::endl;
            }
        }
    }

private:
    // The mutex used to acquire the next tile and the destination of the seeds
    boost::mutex mutex_;
    const unsigned maxClustersPerTile_;

    std::vector<ReaderT> &threadBclReaders_;
    boost::ptr_vector<rta::SingleCycleBclMapper<ReaderT> > threadBclMappers_;

    void loadTileCycle(
        rta::SingleCycleBclMapper<ReaderT> &bclMapper,
        Barcodes::iterator destination,
        const flowcell::Layout &flowcellLayout,
        const flowcell::TileMetadata &tile,
        const unsigned cycle)
    {
        bclMapper.mapTileCycle(flowcellLayout, tile, cycle);

        for (unsigned int clusterId = 0; tile.getClusterCount() > clusterId; ++clusterId)
        {
            char base = 0;
            bclMapper.get(clusterId, &base);
            Barcode &barcode = *destination++;

            const Kmer barcodeBase = (0 == base) ? 4 : (base & 3);
            barcode.setSequence((barcode.getSequence() << BITS_PER_BASE) | barcodeBase);
        }
    }
};

template <typename ReaderT>
class BarcodeLoader: boost::noncopyable
{
public:
    BarcodeLoader(
        common::ThreadVector &threads,
        const unsigned inputLoadersMax,
        const unsigned maxClustersPerTile,
        std::vector<ReaderT> &threadReaders):
            inputLoadersMax_(inputLoadersMax),
            threads_(threads),
            parallelBarcodeLoader_(maxClustersPerTile, threadReaders)
    {
    }

    /**
     * \brief resizes and fills barcodes with data.
     */
    void loadBarcodes(
        const unsigned unknownBarcodeIndex,
        const flowcell::Layout &bclFlowcellLayout,
        const flowcell::TileMetadataList &tiles,
        Barcodes &barcodes)
    {
        BarcodeMemoryManager::allocate(tiles, barcodes);
        // Start and execute the threads
        ISAAC_THREAD_CERR << "Loading data on " << inputLoadersMax_ << " threads" << std::endl;

        std::vector<flowcell::TileMetadata>::const_iterator nextTile = tiles.begin();
        Barcodes::iterator nextTileBarcodes = barcodes.begin();
        threads_.execute(boost::bind(&ParallelBarcodeLoader<ReaderT>::load, &parallelBarcodeLoader_,
                                     unknownBarcodeIndex,
                                     boost::ref(barcodes),
                                     boost::ref(bclFlowcellLayout),
                                     boost::ref(nextTileBarcodes),
                                     boost::ref(nextTile),
                                     tiles.end(), _1),
                         inputLoadersMax_);
    }

private:
    static unsigned getTotalBarcodeCount(const flowcell::TileMetadataList &tiles)
    {
        const uint64_t totalClusterCount = std::accumulate(tiles.begin(), tiles.end(), 0,
                                                             boost::bind(std::plus<unsigned>(), _1,
                                                                         boost::bind(&flowcell::TileMetadata::getClusterCount, _2)));
        return totalClusterCount;
    }

    const unsigned inputLoadersMax_;

    common::ThreadVector &threads_;
    ParallelBarcodeLoader<ReaderT> parallelBarcodeLoader_;
};

} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_BARCODE_LOADER_HH
