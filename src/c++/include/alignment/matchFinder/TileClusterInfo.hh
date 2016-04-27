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
 ** \file TileClusterBarcode.hh
 **
 ** \brief Holds the information about cluster barcode mapping.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH
#define iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/Seed.hh"
#include "alignment/SeedMetadata.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace matchFinder
{

/**
 * \brief Contains the cluster barcode index
 */
class ClusterInfo
{
public:
    static const unsigned short MAX_BARCODE_VALUE = 0xffff;
    ClusterInfo() : barcode_(MAX_BARCODE_VALUE)
    {
    }

    unsigned getBarcodeIndex() const
    {
        return barcode_;
    }

    bool isBarcodeSet() const
    {
        return MAX_BARCODE_VALUE != getBarcodeIndex();
    }

    void setBarcodeIndex(const unsigned short barcodeIndex)
    {
        ISAAC_ASSERT_MSG(barcodeIndex != MAX_BARCODE_VALUE, "Invalid barcode index: " << barcodeIndex);
        barcode_ = barcodeIndex;
    }

private:
    unsigned short barcode_;
};

typedef std::vector<ClusterInfo> ClusterInfos;

inline std::ostream& operator << (std::ostream &os, const ClusterInfo &cluster)
{
    return os << "ClusterInfo(" << cluster.getBarcodeIndex() <<
        ")";
}

/**
 * \brief Geometry: [tileIndex][clusterIndex].
 *
 * For each tile contains mapping between the tile cluster id and the index
 * of the barcode found for the cluster.
 */
struct TileClusterInfo : std::vector<ClusterInfos >
{
    TileClusterInfo(const flowcell::TileMetadataList &unprocessedTileMetadataList):
        std::vector<ClusterInfos >(
            std::max_element(
                unprocessedTileMetadataList.begin(), unprocessedTileMetadataList.end(),
                [](const flowcell::TileMetadata &left, const flowcell::TileMetadata &right)
                {return left.getIndex() < right.getIndex();})->getIndex() + 1)
    {
        BOOST_FOREACH(const flowcell::TileMetadata &unprocessedTile, unprocessedTileMetadataList)
        {
            ClusterInfos &tile = at(unprocessedTile.getIndex());
            tile.resize(unprocessedTile.getClusterCount());
        }
    }
    unsigned getBarcodeIndex(const unsigned tileIndex, const unsigned clusterIndex) const
    {
        return at(tileIndex).at(clusterIndex).getBarcodeIndex();
    }

    void setBarcodeIndex(const unsigned tileIndex, const unsigned clusterIndex, const unsigned barcodeIndex)
    {
        at(tileIndex).at(clusterIndex).setBarcodeIndex(barcodeIndex);
    }
};

} // namespace matchFinder
} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_DEMULTIPLEXING_TILE_CLUSTER_INFO_HH
