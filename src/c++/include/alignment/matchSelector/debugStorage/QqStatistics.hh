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

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_QQ_STATISTICS_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_QQ_STATISTICS_HH

#include <atomic>

#include "alignment/matchSelector/FragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

class QqStatistics
{
    const reference::ContigList &contigList_;
    const boost::filesystem::path outputFilePath_;

    std::array<std::atomic<std::size_t>, 61> basesByQScore_;
    std::array<std::atomic<std::size_t>, 61> mismatchesByQscore_;

    void updateHistogram(const FragmentMetadata& fragment);

public:
    QqStatistics(
        const reference::ContigLists &contigLists,
        const boost::filesystem::path &outputFilePath);

    ~QqStatistics();

    void updateStat(
        const unsigned readIndex,
        const BamTemplate &bamTemplate,
        bool onlyGood);
};

} // namespace debugStorage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_QQ_STATISTICS_HH
