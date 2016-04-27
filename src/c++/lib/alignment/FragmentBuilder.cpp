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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentBuilder.hh
 ** 
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/Mismatch.hh"
#include "alignment/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

FragmentBuilder::FragmentBuilder(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned maxSeedsPerRead,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const bool noSmithWaterman,
    const bool splitAlignments,
    const AlignmentCfg &alignmentCfg,
    Cigar &cigarBuffer,
    const bool reserveBuffers)
    : repeatThreshold_(repeatThreshold)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , smitWatermanGapsMax_(smitWatermanGapsMax)
    , noSmithWaterman_(noSmithWaterman)
    , splitAlignments_(splitAlignments)
    , alignmentCfg_(alignmentCfg)
    , cigarBuffer_(cigarBuffer)
    , ungappedAligner_(collectMismatchCycles, alignmentCfg_)
    , gappedAligner_(collectMismatchCycles, flowcellLayoutList, smartSmithWaterman, alignmentCfg_)
    , splitReadAligner_(collectMismatchCycles, alignmentCfg_, splitAlignments_)
    , matches_()
    , offsetMismatches_(flowcell::getMaxReadLength(flowcellLayoutList))
{
    if (reserveBuffers)
    {
        matches_.reserve(maxSeedsPerRead * repeatThreshold_);
    }
}

//inline uint64_t hashPosition(uint64_t contig, uint64_t position)
//{
//    return contig + position;
//}
//

/**
 * \brief Removes entries designating the same alignment location and strand
 *
 * This is important not only to avoid repetitive processing but also to avoid good alignments
 * scored poorly because they are present multiple times in the list.
 *
 * \param removeUnaligned if true, entries with !isAlignment() are removed
 *
 * \postcondition the fragmentList is ordered
 * \postcondition The fragmentList contains unique alignments only
 */
void FragmentBuilder::consolidateDuplicateAlignments(
    FragmentMetadataList &fragmentList,
    const bool removeUnaligned,
    std::size_t &perfectFound)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "consolidateDuplicateFragments initial size: " << fragmentList.size());
    perfectFound = 0;
    // although initially matches arrive ordered by location, gapped alignment might have moved the
    // start position of some.
    std::sort(fragmentList.begin(), fragmentList.end());
    FragmentMetadataList::iterator lastFragment = fragmentList.begin();
    while(fragmentList.end() != lastFragment && removeUnaligned && !lastFragment->isAligned())
    {
        ++lastFragment;
    }

    lastFragment = fragmentList.erase(fragmentList.begin(), lastFragment);

    if (2 <= fragmentList.size())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    keeping   " << *lastFragment);
        perfectFound += lastFragment->isPerfectMatch();

        for (FragmentMetadataList::iterator currentFragment = lastFragment + 1; fragmentList.end() != currentFragment; ++currentFragment)
        {
            if (removeUnaligned && !currentFragment->isAligned())
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    ignoring " << *currentFragment);
            }
            else if (*lastFragment == *currentFragment)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    merging " << *currentFragment);
                lastFragment->mergeAnchors(*currentFragment);
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    merged  " << *lastFragment);
            }
            else
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList[0].getCluster().getId(), "    keeping   " << *currentFragment);
                perfectFound += currentFragment->isPerfectMatch();
                ++lastFragment;
                if (lastFragment != currentFragment)
                {
                    *lastFragment = *currentFragment;
                }
            }
        }
        fragmentList.resize(1 + lastFragment - fragmentList.begin());
    }

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
        fragmentList[0].getCluster().getId(),
        "consolidateDuplicateFragments consolidated size: " << fragmentList.size() << " found " << perfectFound << " perfect alignments");
}

int64_t FragmentBuilder::getReadPosition(
    const flowcell::ReadMetadata &readMetadata,
    const SeedMetadata &seedMetadata, const SeedId &seedId,
    const reference::ReferencePosition &seedLocation, const bool reverse) const
{
    const int seedOffset = seedMetadata.getOffset();
    if (reverse)
    {
        const unsigned readLength = readMetadata.getLength();
        const unsigned seedLength = seedId.getSeedLength();
        // 'seedPosition + seedLength + seedOffset' is the first position past the end of the read
        return seedLocation.getPosition() + seedLength + seedOffset - readLength;
    }
    else
    {
        return seedLocation.getPosition() - seedOffset;
    }
}

} // namespace alignment
} // namespace isaac
