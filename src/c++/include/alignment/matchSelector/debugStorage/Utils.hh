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
 ** \file DebugFragmentStorage.hh
 **
 ** \brief Compares alignment result with truth data in the read name and produces various statistics for accuracy debugging.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_UTILS_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_UTILS_HH

#include "alignment/FragmentMetadata.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

reference::ReferencePosition getAlignmentPositionFromName(const std::size_t readNumber, const FragmentMetadata &fragment);
bool alignsCorrectly(const std::size_t readNumber, const FragmentMetadata &fragment);

} // namespace debugStorage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_UTILS_HH
