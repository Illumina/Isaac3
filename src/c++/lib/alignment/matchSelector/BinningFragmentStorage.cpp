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
 ** \file BinningFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <cerrno>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

BinningFragmentStorage::BinningFragmentStorage(
    const bool keepUnaligned,
    const BinIndexMap &binIndexMap,
    alignment::BinMetadataList &binMetadataList,
    const uint64_t expectedBinSize):
        FragmentBinner(keepUnaligned, binIndexMap.getTotalBins(), binIndexMap, binMetadataList.size(), expectedBinSize),
        binIndexMap_(binIndexMap),
        binMetadataList_(binMetadataList)
{
    FragmentBinner::open(binMetadataList_.begin(), binMetadataList_.end());
}

BinningFragmentStorage::~BinningFragmentStorage()
{
    FragmentBinner::reclaimStorage(binMetadataList_.begin(), binMetadataList_.end());
}

void BinningFragmentStorage::store(
    const BamTemplate &bamTemplate,
    const unsigned barcodeIdx)
{
    common::StaticVector<char, READS_MAX * (sizeof(io::FragmentHeader) + FRAGMENT_BYTES_MAX)> buffer;
    if (2 == bamTemplate.getFragmentCount())
    {
        packPairedFragment(bamTemplate, 0, barcodeIdx, binIndexMap_, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment0 = reinterpret_cast<const io::FragmentAccessor&>(buffer.front());

        packPairedFragment(bamTemplate, 1, barcodeIdx, binIndexMap_, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment1 = *reinterpret_cast<const io::FragmentAccessor*>(&buffer.front() + fragment0.getTotalLength());

        storePaired(fragment0, fragment1);
    }
    else
    {
        packSingleFragment(bamTemplate, barcodeIdx, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment = reinterpret_cast<const io::FragmentAccessor&>(buffer.front());
        storeSingle(fragment);
    }
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
