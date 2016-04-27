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
 ** \file MatchSelectorTileStats.hh
 **
 ** \brief Statistics aggregation helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H
#define ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H

#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

struct TileStats
{

    explicit TileStats(): fragmentCount_(0), alignedFragmentCount_(0), uniquelyAlignedFragmentCount_(0), adapterBases_(0)
    {
        std::fill(cycleBlanks_, cycleBlanks_ + MAX_CYCLES, 0);
        std::fill(cycleMismatches_, cycleMismatches_ + MAX_CYCLES, 0);
        std::fill(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES, 0);

        std::fill(cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES, 0);
        std::fill(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES, 0);
    }

    void reset()
    {
        *this = TileStats();
    }

    void reserve(const unsigned reserveClusters)
    {
    }

    static const unsigned MAX_CYCLES = 1024;
    uint64_t cycleBlanks_[MAX_CYCLES];
    uint64_t cycleUniquelyAlignedBlanks_[MAX_CYCLES];
    uint64_t cycleMismatches_[MAX_CYCLES];
    uint64_t cycleUniquelyAlignedMismatches_[MAX_CYCLES];

    int64_t cycleUniquelyAligned1MismatchFragments_[MAX_CYCLES];
    int64_t cycleUniquelyAligned2MismatchFragments_[MAX_CYCLES];
    int64_t cycleUniquelyAligned3MismatchFragments_[MAX_CYCLES];
    int64_t cycleUniquelyAligned4MismatchFragments_[MAX_CYCLES];
    int64_t cycleUniquelyAlignedMoreMismatchFragments_[MAX_CYCLES];

    int64_t cycle1MismatchFragments_[MAX_CYCLES];
    int64_t cycle2MismatchFragments_[MAX_CYCLES];
    int64_t cycle3MismatchFragments_[MAX_CYCLES];
    int64_t cycle4MismatchFragments_[MAX_CYCLES];
    int64_t cycleMoreMismatchFragments_[MAX_CYCLES];

    uint64_t fragmentCount_;
    uint64_t alignedFragmentCount_;
    uint64_t uniquelyAlignedFragmentCount_;
    uint64_t adapterBases_;

    template<typename TemplateT>
    void recordTemplate(const TemplateT &templ)
    {
    }
    template<typename FragmentT>
    void recordFragment(const bool collectCycleStats, const FragmentT &fragment, const flowcell::ReadMetadata &readMetadata)
    {
        ++fragmentCount_;
        alignedFragmentCount_ += fragment.isAligned();
        uniquelyAlignedFragmentCount_ += fragment.isUniquelyAligned();

        adapterBases_ += fragment.getAdapterClipped();

        if (collectCycleStats)
        {
            ISAAC_ASSERT_MSG(readMetadata.getFirstCycle() + readMetadata.getLength() < MAX_CYCLES,
                             "Cycle number is too great, check the MAX_CYCLES constant.");

            std::transform(fragment.getForwardSequenceBegin(), fragment.getForwardSequenceEnd(),
                           cycleBlanks_ + readMetadata.getFirstCycle(),
                           cycleBlanks_ + readMetadata.getFirstCycle(),
                           [](char first, uint64_t second){return second + (first == oligo::SEQUENCE_OLIGO_N);});
            std::for_each(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd(),
                          boost::bind(&TileStats::incrementCycleMismatches, this, _1));

            std::size_t mismatches = 0;
            BOOST_FOREACH(const unsigned short &cycle, std::make_pair(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd()))
            {
                incrementCycleXMismatchFragments(cycle, ++mismatches);
            }

            if (fragment.isUniquelyAligned())
            {
                std::transform(fragment.getForwardSequenceBegin(), fragment.getForwardSequenceEnd(),
                               cycleUniquelyAlignedBlanks_ + readMetadata.getFirstCycle(),
                               cycleUniquelyAlignedBlanks_ + readMetadata.getFirstCycle(),
                               [](char first, uint64_t second){return second + (first == oligo::SEQUENCE_OLIGO_N);});
    //                           boost::bind(std::plus<uint64_t>(), _2, boost::bind(&boost::cref<char>, _1) == oligo::SEQUENCE_OLIGO_N));
                std::for_each(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd(),
                              boost::bind(&TileStats::incrementCycleUniquelyAlignedMismatches, this, _1));

                mismatches = 0;
                BOOST_FOREACH(const unsigned short &cycle, std::make_pair(fragment.mismatchCyclesBegin(), fragment.mismatchCyclesEnd()))
                {
                    incrementCycleUniquelyAlignedXMismatchFragments(cycle, ++mismatches);
                }
            }
        }
    }

    const TileStats &operator +=(const TileStats &right)
    {
        std::transform(cycleBlanks_, cycleBlanks_ + MAX_CYCLES,
                       right.cycleBlanks_, cycleBlanks_,
                       std::plus<uint64_t>());

        std::transform(cycleMismatches_, cycleMismatches_ + MAX_CYCLES,
                       right.cycleMismatches_, cycleMismatches_,
                       std::plus<uint64_t>());

        std::transform(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES,
                       right.cycle1MismatchFragments_, cycle1MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       right.cycle2MismatchFragments_, cycle2MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       right.cycle3MismatchFragments_, cycle3MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       right.cycle4MismatchFragments_, cycle4MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       right.cycleMoreMismatchFragments_, cycleMoreMismatchFragments_,
                       std::plus<uint64_t>());


        std::transform(cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedBlanks_, cycleUniquelyAlignedBlanks_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedMismatches_, cycleUniquelyAlignedMismatches_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::plus<uint64_t>());

        std::transform(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       right.cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_,
                       std::plus<uint64_t>());

        fragmentCount_ += right.fragmentCount_;
        alignedFragmentCount_ += right.alignedFragmentCount_;
        uniquelyAlignedFragmentCount_ += right.uniquelyAlignedFragmentCount_;
        adapterBases_ += right.adapterBases_;

        return *this;
    }

    const TileStats operator +(const TileStats &right) const
    {
        TileStats ret(*this);
        ret += right;
        return ret;
    }

    void finalize()
    {
        // if at cycle 1 the fragment had 1 mismatch, this mismatch propagates to all subsequent cycles.
        std::transform(cycle1MismatchFragments_ + 1, cycle1MismatchFragments_ + MAX_CYCLES,
                       cycle1MismatchFragments_, cycle1MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycle2MismatchFragments_ + 1, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle2MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycle3MismatchFragments_ + 1, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle3MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycle4MismatchFragments_ + 1, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycle4MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycleMoreMismatchFragments_ + 1, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + 1,
                       std::plus<int64_t>());

        // once a 1-mismatch fragments gets a second mismatch, it stops beinga 1-mismatch fragment. Substract.
        std::transform(cycle1MismatchFragments_, cycle1MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle1MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle2MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycle3MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycleMoreMismatchFragments_, cycle4MismatchFragments_,
                       std::minus<int64_t>());

        // two mismatch fragments include one mismatch fragments and so on, add sequentially
        std::transform(cycle2MismatchFragments_, cycle2MismatchFragments_ + MAX_CYCLES,
                       cycle1MismatchFragments_, cycle2MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycle3MismatchFragments_, cycle3MismatchFragments_ + MAX_CYCLES,
                       cycle2MismatchFragments_, cycle3MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycle4MismatchFragments_, cycle4MismatchFragments_ + MAX_CYCLES,
                       cycle3MismatchFragments_, cycle4MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycleMoreMismatchFragments_, cycleMoreMismatchFragments_ + MAX_CYCLES,
                       cycle4MismatchFragments_, cycleMoreMismatchFragments_,
                       std::plus<int64_t>());

        // if at cycle 1 the fragment had 1 mismatch, this mismatch propagates to all subsequent cycles.
        std::transform(cycleUniquelyAligned1MismatchFragments_ + 1, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycleUniquelyAligned2MismatchFragments_ + 1, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycleUniquelyAligned3MismatchFragments_ + 1, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycleUniquelyAligned4MismatchFragments_ + 1, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + 1,
                       std::plus<int64_t>());

        std::transform(cycleUniquelyAlignedMoreMismatchFragments_ + 1, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + 1,
                       std::plus<int64_t>());

        // once a 1-mismatch fragments gets a second mismatch, it stops beinga 1-mismatch fragment. Substract.
        std::transform(cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned1MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned1MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::minus<int64_t>());
        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::minus<int64_t>());

        // two mismatch fragments include one mismatch fragments and so on, add sequentially
        std::transform(cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned2MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned1MismatchFragments_, cycleUniquelyAligned2MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned3MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned2MismatchFragments_, cycleUniquelyAligned3MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAligned4MismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned3MismatchFragments_, cycleUniquelyAligned4MismatchFragments_,
                       std::plus<int64_t>());
        std::transform(cycleUniquelyAlignedMoreMismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_ + MAX_CYCLES,
                       cycleUniquelyAligned4MismatchFragments_, cycleUniquelyAlignedMoreMismatchFragments_,
                       std::plus<int64_t>());
    }

private:
    void incrementCycleMismatches(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        ++cycleMismatches_[cycle];
    }

    void incrementCycleUniquelyAlignedMismatches(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        ++cycleUniquelyAlignedMismatches_[cycle];
    }

    void incrementCycleUniquelyAlignedXMismatchFragments(const unsigned short cycle, const unsigned mismatches)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        switch (mismatches)
        {
        case 0:
            ISAAC_ASSERT_MSG(0 != mismatches, "incrementCycleUniquelyAlignedXMismatchFragments should not be called for 0 mismatches");
            break;
        case 1:
            ++cycleUniquelyAligned1MismatchFragments_[cycle];
            break;
        case 2:
            ++cycleUniquelyAligned2MismatchFragments_[cycle];
            break;
        case 3:
            ++cycleUniquelyAligned3MismatchFragments_[cycle];
            break;
        case 4:
            ++cycleUniquelyAligned4MismatchFragments_[cycle];
            break;
        case 5:
            ++cycleUniquelyAlignedMoreMismatchFragments_[cycle];
            break;
        default:
            // don't count fragments with more than 5 mismatches;
            break;
        }
    }

    void incrementCycleXMismatchFragments(const unsigned short cycle, const unsigned mismatches)
    {
        ISAAC_ASSERT_MSG(cycle < MAX_CYCLES, "Cycle number too high.");
        switch (mismatches)
        {
        case 0:
            ISAAC_ASSERT_MSG(0 != mismatches, "incrementCycleXMismatchFragments should not be called for 0 mismatches");
            break;
        case 1:
            ++cycle1MismatchFragments_[cycle];
            break;
        case 2:
            ++cycle2MismatchFragments_[cycle];
            break;
        case 3:
            ++cycle3MismatchFragments_[cycle];
            break;
        case 4:
            ++cycle4MismatchFragments_[cycle];
            break;
        case 5:
            ++cycleMoreMismatchFragments_[cycle];
            break;
        default:
            // don't count fragments with more than 5 mismatches;
            break;
        }
    }
};

} //namespace matchSelector
} //namespace alignment
} //namespace isaac

#endif //ISAAC_ALIGNMENT_MATCH_SELECTOR_TILE_STATS_H
