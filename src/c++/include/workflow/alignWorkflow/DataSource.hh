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
 ** \file DataSource.hh
 **
 ** \brief Abstraction of data source
 **
 ** \author Roman Petrovski
 **/


#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH

#include "alignment/BclClusters.hh"
#include "flowcell/TileMetadata.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

struct TileSource : boost::noncopyable
{
    /**
     * \brief Returns set of tiles that can be processed together.
     *
     * \return If returned set is empty, there is nothing left to process
     */
    virtual flowcell::TileMetadataList discoverTiles() = 0;
    virtual ~TileSource(){}
};

struct BarcodeSource : boost::noncopyable
{
    /// Allocate memory and load barcodes for the tiles.
    virtual void loadBarcodes(
        const flowcell::Layout &flowcell,
        const unsigned unknownBarcodeIndex,
        const flowcell::TileMetadataList &tiles, demultiplexing::Barcodes &barcodes) = 0;

    virtual ~BarcodeSource(){}
};


template <typename DataSourceT>
struct DataSourceTraits
{
//    static const bool SUPPORTS_XY = false;
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac


#endif //#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_DATA_SOURCE_HH
