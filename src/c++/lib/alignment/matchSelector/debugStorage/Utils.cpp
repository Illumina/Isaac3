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
 ** \file DebugFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <atomic>

#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "common/Debug.hh"
#include "common/FastIo.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

std::pair<BclClusters::const_iterator, BclClusters::const_iterator> getReadName(const std::size_t readIndex, const FragmentMetadata &fragment)
{
    std::pair<BclClusters::const_iterator, BclClusters::const_iterator> ret = {fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()};
    if (ret.second == ret.first)
    {
        return ret;
    }

    if (!readIndex)
    {
        ret.first = fragment.getCluster().nameBegin();
        ret.second = std::find(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd(), '-');
    }
    else
    {
        ret.first = std::find(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd(), '-') + 1;
        ret.second = fragment.getCluster().nameEnd();
    }

    return ret;
}

reference::ReferencePosition getAlignmentPositionFromName(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    // numbers are 1-based
    const auto name = getReadName(readNumber - 1, fragment);

    if (name.second == name.first)
    {
        return reference::ReferencePosition(reference::ReferencePosition::TooManyMatch);
    }

    if ('u' == *name.first)
    {
        ISAAC_ASSERT_MSG(false, common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) << " " << fragment);
        return reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    return reference::ReferencePosition(
        std::atol(&*name.first + 2),
        std::atol(&*std::find(name.first + 2, name.second, ':') + 1),
        false,
        'r' == *name.first);
}

bool alignsCorrectly(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    const reference::ReferencePosition oriPos = getAlignmentPositionFromName(readNumber, fragment);
    if (oriPos.isTooManyMatch())
    {
        return true;
    }
//    ISAAC_THREAD_CERR << "oriPos:" << oriPos << " name " << common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) << std::endl;
    return fragment.isReverse() == oriPos.reverse() &&
        fragment.getContigId() == oriPos.getContigId() &&
        uint64_t(fragment.getPosition()) == oriPos.getPosition();
}

} // namespace debugStroage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac
