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

#include "alignment/matchSelector/debugStorage/QqStatistics.hh"
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

template <typename AlignmentsByScore>
void traceQQ(
    std::ostream &os,
    const AlignmentsByScore &basesByQScore,
    const AlignmentsByScore &mismatchesByQscore)
{
    for (std::size_t i= 0; i < basesByQScore.size(); ++i)
    {
        os << i << "\t" << -10 * log10(double(mismatchesByQscore[i]) / double(basesByQScore[i])) << "\t" << mismatchesByQscore[i] << "\t" << basesByQScore[i] << std::endl;
    }
}

QqStatistics::~QqStatistics()
{
    ISAAC_THREAD_CERR << "Storing QQ stats in " << outputFilePath_.c_str() << std::endl;
    std::ofstream qqTsv(outputFilePath_.c_str());
    traceQQ(qqTsv, basesByQScore_, mismatchesByQscore_);
}

void QqStatistics::updateHistogram(const FragmentMetadata& fragment)
{
    const std::vector<char>& sequence = fragment.getRead().getStrandSequence(fragment.isReverse());
    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    const std::vector<char>::const_iterator sequenceEnd = sequence.end();
    std::vector<char>::const_iterator qualityBegin = fragment.getRead().getStrandQuality(fragment.isReverse()).begin();
    const reference::Contig& contig = contigList_.at(fragment.getContigId());
    reference::Contig::const_iterator reference = contig.begin() + fragment.getPosition();
    //            ISAAC_THREAD_CERR << "seq: " << common::makeFastIoString(sequenceBegin, sequenceEnd) << std::endl;
    //            ISAAC_THREAD_CERR << "ref: " << common::makeFastIoString(reference, reference + std::distance(sequenceBegin, sequenceEnd)) << std::endl;
//    std::size_t cycles = 32;
    while (/*cycles && */sequenceEnd != sequenceBegin)
    {
        ++basesByQScore_.at(*qualityBegin);
        if (!alignment::isMatch(*sequenceBegin, *reference))
        {
            ++mismatchesByQscore_.at(*qualityBegin);
        }
        ++sequenceBegin;
        ++reference;
        ++qualityBegin;
//        --cycles;
    }
}

void QqStatistics::updateStat(
    const unsigned readIndex,
    const BamTemplate &bamTemplate,
    const bool onlyGood)
{
    const FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(readIndex);
    if (!onlyGood || (fragment.isUniquelyAligned() && !fragment.dodgy && !fragment.gapCount && !fragment.getBeginClippedLength() && !fragment.getEndClippedLength()))
    {
        if (2 != bamTemplate.getFragmentCount())
        {
            updateHistogram(fragment);
        }
        else if (bamTemplate.isProperPair())
        {
            const FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);
            if (!onlyGood || (mate.isUniquelyAligned() && !mate.dodgy && !mate.gapCount && !mate.getBeginClippedLength() && !mate.getEndClippedLength()))
            {
                updateHistogram(fragment);
            }
        }
    }
}


QqStatistics::QqStatistics(
    const reference::ContigLists &contigLists,
    const boost::filesystem::path &outputFilePath):
    contigList_(contigLists.at(0)),
    outputFilePath_(outputFilePath)
{
    for (auto &a : basesByQScore_) {a = 0;}
    for (auto &a : mismatchesByQscore_) {a = 0;}

    ISAAC_ASSERT_MSG(1 == contigLists.size(), "multiple references are not supported");
}

} // namespace debugStroage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac
