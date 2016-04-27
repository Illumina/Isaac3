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
 ** \file IndelLoader.hh
 **
 ** Loads indels from a vcf file.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_INDEL_LOADER_HH
#define iSAAC_BUILD_INDEL_LOADER_HH

#include <boost/filesystem.hpp>

#include "build/gapRealigner/RealignerGaps.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace build
{

build::gapRealigner::Gaps loadIndels(
    const boost::filesystem::path &vcfFilePath,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList);

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_INDEL_LOADER_HH
