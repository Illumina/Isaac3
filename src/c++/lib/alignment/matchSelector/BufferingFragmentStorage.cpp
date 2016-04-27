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
 ** \file BufferingFragmentStorage.cpp
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
#include "alignment/matchSelector/BufferingFragmentStorage.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

BufferingFragmentStorage::BufferingFragmentStorage(
    const bool keepUnaligned,
    const unsigned maxThreads,
    const unsigned maxSavers,
    const BinIndexMap &binIndexMap,
    alignment::BinMetadataList &binMetadataList,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const uint64_t expectedBinSize)
    : FragmentBinner(keepUnaligned, maxSavers, binIndexMap, binMetadataList.size(), expectedBinSize)
    , keepUnaligned_(keepUnaligned)
    , binIndexMap_(binIndexMap)
    , filesAtATime_(maxSavers)
    , threads_(maxThreads)
    , binMetadataList_(binMetadataList)
    , flushBuffer_(flowcellLayoutList)
    , storeBuffer_(flowcellLayoutList)
{
}

BufferingFragmentStorage::~BufferingFragmentStorage()
{
    FragmentBinner::reclaimStorage(binMetadataList_.begin(), binMetadataList_.end());
}

void BufferingFragmentStorage::store(const BamTemplate &bamTemplate, const unsigned barcodeIdx)
{
    const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(0);
    if (2 == bamTemplate.getFragmentCount())
    {
        FragmentPacker::packPairedFragment(bamTemplate, 0, barcodeIdx, binIndexMap_,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 0));
        FragmentPacker::packPairedFragment(bamTemplate, 1, barcodeIdx, binIndexMap_,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 1));
    }
    else
    {
        FragmentPacker::packSingleFragment(bamTemplate, barcodeIdx,
                                           storeBuffer_.getRecordInsertIterator(fragment.getCluster().getId(), 0));
    }
}

void BufferingFragmentStorage::reset(const uint64_t clusterId, const bool paired)
{

    io::FragmentAccessor &fragment1 = *reinterpret_cast<io::FragmentAccessor*>(
        &*storeBuffer_.getRecordInsertIterator(clusterId, 0));
    fragment1.flags_.initialized_ = false;

    if (paired)
    {
        io::FragmentAccessor &fragment2 = *reinterpret_cast<io::FragmentAccessor*>(
            &*storeBuffer_.getRecordInsertIterator(clusterId, 1));
        fragment2.flags_.initialized_ = false;
    }
}


void BufferingFragmentStorage::prepareFlush() noexcept
{
    storeBuffer_.swap(flushBuffer_);
}

void BufferingFragmentStorage::flush()
{
    ISAAC_THREAD_CERR << "Flushing buffer" << std::endl;

    for (alignment::BinMetadataList::iterator binIterator = binMetadataList_.begin();
        binMetadataList_.end() != binIterator; )
    {
        const unsigned binsToFlush = std::min<unsigned>(filesAtATime_, std::distance(binIterator, binMetadataList_.end()));
        open(binIterator, binIterator + binsToFlush, threads_);
        ISAAC_THREAD_CERR << "Flushing bins " << std::distance(binMetadataList_.begin(), binIterator) << "-" << std::distance(binMetadataList_.begin(), binIterator) + binsToFlush << std::endl;

        threads_.execute(
            [this](const unsigned threadNumber, const unsigned threadsTotal)
            {
                const unsigned threadsToUse = std::min<unsigned>(threadsTotal, flushBuffer_.getClusters());
                if (threadNumber >= threadsToUse)
                {
                    return ;
                }

                const std::size_t blockSize = flushBuffer_.getClusters() / threadsToUse;
                // last thread has to do the remainder
                FragmentBuffer::const_iterator it = flushBuffer_.advanceClusters(flushBuffer_.dataBegin(), blockSize * threadNumber);
                FragmentBuffer::const_iterator endIt = (threadsToUse - 1) == threadNumber ?
                    flushBuffer_.dataEnd() : flushBuffer_.advanceClusters(it, blockSize);

                for (; endIt != it; it = flushBuffer_.nexCluster(it))
                {
                    const io::FragmentAccessor &fragment0 = *reinterpret_cast<const io::FragmentAccessor*>(&*it);
                    if (fragment0.flags_.initialized_)
                    {
                        if(fragment0.flags_.paired_)
                        {
                            const FragmentBuffer::const_iterator r2It = flushBuffer_.nextRead(it);
                            ISAAC_ASSERT_MSG(flushBuffer_.dataEnd() != r2It, "Unexpected end of buffer reached");
                            const io::FragmentAccessor &fragment1 = *reinterpret_cast<const io::FragmentAccessor*>(&*r2It);
                            ISAAC_ASSERT_MSG(fragment1.flags_.initialized_, "Both reads have to be either initialized or not: " << fragment0 << " " << fragment1);
                            storePaired(fragment0, fragment1);
                        }
                        else
                        {
                            storeSingle(fragment0);
                        }
                    }
                }
            },
            std::min(threads_.size(), flushBuffer_.getClusters())
        );

        binIterator += binsToFlush;
    }

    ISAAC_THREAD_CERR << "Flushing buffer done " << std::endl;
}


} //namespace matchSelector
} // namespace alignment
} // namespace isaac
