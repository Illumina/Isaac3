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

#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "alignment/matchSelector/debugStorage/Utils.hh"
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

//std::array<std::array<std::atomic<std::size_t>, 6 >, 6> badAlignmentsByGoodAndBadSeedCount;
//std::array<std::array<std::atomic<std::size_t>, 6 >, 6> goodAlignmentsByGoodAndBadSeedCount;

//void traceMatrix(
//    const std::array<std::array<std::atomic<std::size_t>, 6>, 6>& alignmentsByGoodAndBadSeedCount,
//    std::string matrix)
//{
//    matrix += "\n";
//    for (const auto& line : alignmentsByGoodAndBadSeedCount)
//    {
//        for (const auto& cell : line)
//        {
//            matrix += "\t";
//            matrix += std::to_string(cell);
//        }
//        matrix += "\n";
//    }
//    std::cout << matrix << std::endl;
//}

template <typename AlignmentsByScore>
void traceQQ(
    std::ostream &os,
    const std::atomic<std::size_t> &unaligned,
    const std::atomic<std::size_t> &downgraded,
    const std::atomic<std::size_t> &downgradedBad,
    const AlignmentsByScore &alignmentsByScore,
    const AlignmentsByScore &badAlignmentsByScore)
{
    os << "-2\t0\t0\t" << unaligned << std::endl;
    os << "-1\t" << -10 * log10(double(downgradedBad) / double(downgraded)) << "\t" << downgradedBad << "\t" << downgraded << std::endl;
    const std::size_t toTrace = std::distance(alignmentsByScore.begin(), std::find_if(alignmentsByScore.rbegin(), alignmentsByScore.rend(), [](const std::size_t c){return c;}).base());
    for (std::size_t i= 0; i < toTrace; ++i)
    {
        os << i << "\t" << -10 * log10(double(badAlignmentsByScore[i]) / double(alignmentsByScore[i])) << "\t" << badAlignmentsByScore[i] << "\t" << alignmentsByScore[i] << std::endl;
    }
}

MapqStatistics::~MapqStatistics()
{
//    traceMatrix(badAlignmentsByGoodAndBadSeedCount, "badAlignmentsByGoodAndBadSeedCount");
//    traceMatrix(goodAlignmentsByGoodAndBadSeedCount, "goodAlignmentsByGoodAndBadSeedCount");
//
    {
        ISAAC_THREAD_CERR << "Storing SM QQ stats in " << outputSmFilePath_.c_str() << std::endl;
        std::ofstream smTsv(outputSmFilePath_.c_str());
        traceQQ(smTsv, unalignedFragments_, downgradedAlignments_, badDowngradedAlignments_, alignmentsBySm_, badAlignmentsBySm_);
    }

    {
        ISAAC_THREAD_CERR << "Storing MAPQ QQ stats in " << outputMapqFilePath_.c_str() << std::endl;
        std::ofstream mapqTsv(outputMapqFilePath_.c_str());
        traceQQ(mapqTsv, unalignedFragments_, downgradedAlignments_, badDowngradedAlignments_, alignmentsByMapq_, badAlignmentsByMapq_);
    }
}
//
//std::size_t getGoodSeedCount(const FragmentMetadata &fragment)
//{
//    return fragment.isAligned();
//}
//
//std::size_t getBadSeedCount(const FragmentMetadata &fragment)
//{
//    return !fragment.isAligned();
//}

void MapqStatistics::updateMapqHistograms(
    const BamTemplate& bamTemplate,
    const std::size_t readIndex,
    const FragmentMetadata& fragment,
    const FragmentMetadata *mate)
{
    ++alignmentsBySm_[std::min<unsigned>(fragment.getAlignmentScore(),
                                        alignmentsBySm_.size() - 1)];
    const std::size_t mapQ =
        (!fragment.isUniquelyAligned() && mate && bamTemplate.isProperPair()) ?
            // Rescue non-unique alignments only if both mate and pair are unique. This prevents from using unique pair score
            // when both fragments are non-unique as it usually results in accepting high score for unique pairing
            // wihtout having seen all of the pairings.
            std::min(bamTemplate.getAlignmentScore(), mate->getAlignmentScore()) :
                     fragment.getAlignmentScore();
    ++alignmentsByMapq_[std::min<unsigned>(mapQ, alignmentsByMapq_.size() - 1)];
    if (!alignsCorrectly(readIndex, fragment))
    {
        if (mapQ >= 40)
        {
            ISAAC_THREAD_CERR<< "Misaligned expected: " << common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) <<
            " got " << fragment << std::endl;
        }
        ++badAlignmentsBySm_[std::min<unsigned >(fragment.getAlignmentScore(), badAlignmentsBySm_.size() - 1)];
        ++badAlignmentsByMapq_[std::min<unsigned >(mapQ, alignmentsByMapq_.size() - 1)];
    }
}

bool MapqStatistics::updateStat(
    const BamTemplate& bamTemplate,
    const std::size_t readNumber,
    const FragmentMetadata &fragment,
    const FragmentMetadata *mate)
{
    if (!fragment.isAligned())
    {
        ++unalignedFragments_;
    }
    else
    {
        if (fragment.dodgy)
        {
            ++downgradedAlignments_;
            if (!alignsCorrectly(readNumber, fragment))
            {
                ++badDowngradedAlignments_;
            }
        }
        else
        {
            updateMapqHistograms(bamTemplate, readNumber, fragment, mate);
        }
    }

    return false;
}


MapqStatistics::MapqStatistics(
    const reference::ContigLists &contigLists,
    const reference::ContigAnnotationsList &contigAnnotationsList,
    const AlignmentCfg &alignmentCfg,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const boost::filesystem::path &outputMapqFilePath,
    const boost::filesystem::path &outputSmFilePath):
    contigList_(contigLists.at(0)),
    contigAnnotations_(contigAnnotationsList.at(0)),
    flowcell_(flowcellLayoutList.at(0)),
    alignmentCfg_(alignmentCfg),
    outputMapqFilePath_(outputMapqFilePath),
    outputSmFilePath_(outputSmFilePath),
    unalignedFragments_(0),
    downgradedAlignments_(0),
    badDowngradedAlignments_(0)
{
    for (auto &a : alignmentsBySm_) {a = 0;}
    for (auto &a : badAlignmentsBySm_) {a = 0;}
    for (auto &a : alignmentsByMapq_) {a = 0;}
    for (auto &a : badAlignmentsByMapq_) {a = 0;}

    ISAAC_ASSERT_MSG(1 == contigLists.size(), "multiple references are not supported");
    ISAAC_ASSERT_MSG(1 == contigAnnotationsList.size(), "multiple references are not supported");
    ISAAC_ASSERT_MSG(1 == flowcellLayoutList.size(), "multiple flowcells are not supported");
}

} // namespace debugStroage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac
