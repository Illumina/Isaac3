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
 ** \file BinData.hh
 **
 ** Performs sorting and duplicate marking on a single alignment bin.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_DATA_HH
#define iSAAC_BUILD_BIN_DATA_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/FragmentAccessorBamAdapter.hh"
#include "build/FragmentIndex.hh"
#include "build/PackedFragmentBuffer.hh"
#include "build/gapRealigner/RealignerGaps.hh"
#include "io/FileBufCache.hh"

namespace isaac
{
namespace build
{

enum GapRealignerMode
{
    /// don't realign
    REALIGN_NONE,
    /// Realign against gaps found within the sample
    REALIGN_SAMPLE,
    /// Realign against gaps found in all samples of the same project
    REALIGN_PROJECT,
    /// Realign against all gaps present in the data
    REALIGN_ALL
};

struct BinData : public std::vector<PackedFragmentBuffer::Index, common::NumaAllocator<PackedFragmentBuffer::Index, common::numa::defaultNodeLocal> >
{
    typedef std::vector<PackedFragmentBuffer::Index, common::NumaAllocator<PackedFragmentBuffer::Index, common::numa::defaultNodeLocal> > BaseType;
    typedef BaseType IndexType;
    typedef std::vector<SeFragmentIndex, common::NumaAllocator<SeFragmentIndex, common::numa::defaultNodeLocal> > SeIdx;
    typedef std::vector<RStrandOrShadowFragmentIndex, common::NumaAllocator<RStrandOrShadowFragmentIndex, common::numa::defaultNodeLocal> > RIdx;
    typedef std::vector<FStrandFragmentIndex, common::NumaAllocator<FStrandFragmentIndex, common::numa::defaultNodeLocal> > FIdx;

public:
    BinData(
        const unsigned realignedGapsPerFragment,
        const BarcodeBamMapping &barcodeBamMapping,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const build::GapRealignerMode realignGaps,
        const gapRealigner::Gaps &knownIndels,
        const alignment::BinMetadata &bin,
        const unsigned binStatsIndex,
        const flowcell::TileMetadataList &tileMetadataList,
        const BuildContigMap &contigMap,
        const isaac::reference::ContigLists &contigLists,
        const unsigned maxReadLength,
        const unsigned char forcedDodgyAlignmentScore,
        const flowcell::FlowcellLayoutList &flowCellLayoutList,
        const IncludeTags includeTags,
        const bool pessimisticMapQ,
        const unsigned splitGapLength,
        const unsigned expectedCoverage) :
            bin_(bin),
            binStatsIndex_(binStatsIndex),
            barcodeBamMapping_(barcodeBamMapping),
            realignGaps_(realignGaps),
            knownIndels_(knownIndels),
            realignerGaps_(getGapGroupsCount()),
            dataDistribution_(bin_.getDataDistribution()),
            inputFileBuf_(),
            bamAdapter_(
                maxReadLength, tileMetadataList, barcodeMetadataList,
                contigMap, contigLists, forcedDodgyAlignmentScore, flowCellLayoutList, includeTags, pessimisticMapQ,
                splitGapLength, splitInfoList_)
    {
        data_.reserve(bin_);

        BaseType::reserve(bin_.getTotalElements() * 2);
        seIdx_.reserve(bin_.getSeIdxElements());
        rIdx_.reserve(bin_.getRIdxElements());
        fIdx_.reserve(bin_.getFIdxElements());
        if (REALIGN_NONE != realignGaps_)
        {
            reserveGaps(bin_, knownIndels_, barcodeMetadataList);
        }
        // assume each split read will not have any noticeable amount of extra elements above of what's listed below
        // translocation gaps with inversion require up to 7 CIGAR components: SOFT_CLIP,ALIGN,FLIP,CONTIG,DELETE,ALIGN,SOFTCLIP
        // so, when this gets broken up, each leftover will have no more than 5 components
        static const unsigned SINGLE_SPLIT_LEFTOVER_COMPONENTS = 5;
        additionalCigars_.reserve(
            // Assuming each split read will result in two separate fragments
            SINGLE_SPLIT_LEFTOVER_COMPONENTS * bin.getTotalSplitCount() * 2 +
            // assume each existing cigar gets realignedGapsPerFragment_ gaps introduced...
            (bin.getTotalCigarLength() + bin.getTotalElements() *
                // assume that to introduce k gaps one will need to have k+1 operations between the gaps
                (realignedGapsPerFragment + realignedGapsPerFragment + 1)));

        // assume each gap produces coverage of realignments and then gets split. This is extremely pessimistic and happens only for
        // --split-gap_length == 1. Now, say all variants in human are gaps. Then 3M gaps with 60 being coverage and only one bin covering
        // the entire genome, this will take about 3M*60*sizeof(SplitInfo) ~= 6G. However in reality we aim to break genome down
        // into few hundred bins. And number of gaps is at least 10 times less than total number of variants. So, for the example above
        // we will be allocating about 60M per bin and using about 1-10% of it.
        // Drawbacks: overestimation will cause us to waste RAM. Underestimation will result in buffer reallocation during splitting
        // which might result in bam generation failure (or might succeed).
        // each split produced two records
        splitInfoList_.reserve((bin_.getTotalGapCount() * expectedCoverage + bin_.getTotalSplitCount()) * 2);

        // summarize chunk sizes to get offsets
        dataDistribution_.tallyOffsets();
        if (!inputFileBuf_.open(bin_.getPathString().c_str(), std::ios_base::binary|std::ios_base::in))
        {
            BOOST_THROW_EXCEPTION(
                common::IoException(errno, (boost::format("Failed to open file %s: %s") % bin_.getPathString() % strerror(errno)).str()));
        }
    }

    void finalize();

    static uint64_t getMemoryRequirements(const alignment::BinMetadata& bin)
    {
        return PackedFragmentBuffer::getMemoryRequirements(bin) +
            bin.getSeIdxElements() * sizeof(SeFragmentIndex) +
            bin.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
            bin.getFIdxElements() * sizeof(FStrandFragmentIndex) +
            bin.getTotalElements() * sizeof(PackedFragmentBuffer::Index);
    }

    void unreserveIndexes()
    {
        SeIdx().swap(seIdx_);
        RIdx().swap(rIdx_);
        FIdx().swap(fIdx_);
    }

    unsigned getBinIndex() const
    {
        return bin_.getIndex();
    }

    bool isUnalignedBin() const {return bin_.isUnalignedBin();}
    uint64_t getUniqueRecordsCount() const {return isUnalignedBin() ? bin_.getTotalElements() : size();}

    BaseType::iterator indexBegin() {return begin();}
    BaseType::iterator indexEnd() {return end();}

    const gapRealigner::RealignerGaps &getRealignerGaps(const unsigned barcode) const;
//private:
    const alignment::BinMetadata &bin_;
    const unsigned binStatsIndex_;
    const BarcodeBamMapping &barcodeBamMapping_;
    SeIdx seIdx_;
    RIdx rIdx_;
    FIdx fIdx_;
    PackedFragmentBuffer data_;
    const GapRealignerMode realignGaps_;
    const gapRealigner::Gaps &knownIndels_;
    alignment::Cigar additionalCigars_;

    SplitInfoList splitInfoList_;
    std::vector<gapRealigner::RealignerGaps> realignerGaps_;
    alignment::BinDataDistribution dataDistribution_;
    std::filebuf inputFileBuf_;
    FragmentAccessorBamAdapter bamAdapter_;

private:
    void reserveGaps(
        const alignment::BinMetadata& bin,
        const gapRealigner::Gaps &knownIndels,
        const flowcell::BarcodeMetadataList &barcodeMetadataList);

    unsigned getGapGroupIndex(const unsigned barcode) const;
    unsigned getGapGroupsCount() const;

};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_SORTER_HH
