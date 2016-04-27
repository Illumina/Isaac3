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
 ** \file FoundMatchesMetadata.hh
 **
 ** \brief auxiliary class
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH

#include "alignment/MatchDistribution.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

struct FoundMatchesMetadata
{
    FoundMatchesMetadata(const boost::filesystem::path &tempDirectory,
                         const flowcell::BarcodeMetadataList &barcodeMetadataList,
                         const unsigned maxIterations,
                         const reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
    {

    }

    flowcell::TileMetadataList tileMetadataList_;

    void addTile(const flowcell::TileMetadata& tile)
    {
        flowcell::TileMetadata tileWithNewIndex(tile, tileMetadataList_.size());
        tileMetadataList_.push_back(tileWithNewIndex);
    }

    void swap(FoundMatchesMetadata &another)
    {
        using std::swap;
        tileMetadataList_.swap(another.tileMetadataList_);
    }
};


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FOUND_MATCHES_METADATA_HH
