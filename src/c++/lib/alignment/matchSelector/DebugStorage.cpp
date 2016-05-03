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

#include "alignment/matchSelector/DebugStorage.hh"
#include "alignment/matchSelector/debugStorage/Utils.hh"
#include "common/Debug.hh"
#include "common/FastIo.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

DebugStorage::DebugStorage(
    const reference::ContigLists &contigLists,
    const reference::ContigAnnotationsList &contigAnnotationsList,
    const AlignmentCfg &alignmentCfg,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const boost::filesystem::path &outputDirectory,
    alignment::matchSelector::FragmentStorage &actualStorage):
    contigList_(contigLists.at(0)),
    contigAnnotations_(contigAnnotationsList.at(0)),
    flowcell_(flowcellLayoutList.at(0)),
    alignmentCfg_(alignmentCfg),
    outputDirectory_(outputDirectory),
    r1MapqStatistics_(
        contigLists, contigAnnotationsList,
        alignmentCfg, flowcellLayoutList,
        outputDirectory / "R1MapQ.tsv",
        outputDirectory / "R1SM.tsv"),
    r2MapqStatistics_(
        contigLists, contigAnnotationsList,
        alignmentCfg, flowcellLayoutList,
        outputDirectory / "R2MapQ.tsv",
        outputDirectory / "R2SM.tsv"),
    pairMapqStatistics_(
        contigLists, contigAnnotationsList,
        alignmentCfg, flowcellLayoutList,
        outputDirectory / "MapQ.tsv",
        outputDirectory / "SM.tsv"),
    r1QqOriginalStatistics_(contigLists, outputDirectory / "R1QQOri.tsv"),
    r2QqOriginalStatistics_(contigLists, outputDirectory / "R2QQOri.tsv"),
    r1QqStatistics_(contigLists, outputDirectory / "R1QQ.tsv"),
    r2QqStatistics_(contigLists, outputDirectory / "R2QQ.tsv"),
    actualStorage_(actualStorage)
{
    ISAAC_ASSERT_MSG(1 == contigLists.size(), "multiple references are not supported");
    ISAAC_ASSERT_MSG(1 == contigAnnotationsList.size(), "multiple references are not supported");
    ISAAC_ASSERT_MSG(1 == flowcellLayoutList.size(), "multiple flowcells are not supported");
    {
        originalCigars_[0].reserve(1);
        originalCigars_[0].addOperation(flowcell_.getReadMetadataList()[0].getLength(), Cigar::ALIGN);
    }
    if (2 == flowcell_.getReadMetadataList().size())
    {
        originalCigars_[1].reserve(1);
        originalCigars_[1].addOperation(flowcell_.getReadMetadataList()[1].getLength(), Cigar::ALIGN);
    }
}

DebugStorage::~DebugStorage()
{
}

bool DebugStorage::restoreOriginal(const std::size_t readNumber, FragmentMetadata &fragment) const
{
    const reference::ReferencePosition oriPos = debugStorage::getAlignmentPositionFromName(readNumber, fragment);
    if (oriPos.isTooManyMatch())
    {
        return false;
    }

    fragment.resetAlignment();
    if (!oriPos.isNoMatch())
    {
        fragment.updateAlignment(
            false,
            alignmentCfg_,
            flowcell_.getReadMetadataList().at(fragment.getReadIndex()),
            contigList_,
            contigAnnotations_,
            oriPos.reverse(),
            oriPos.getContigId(),
            oriPos.getPosition(),
            originalCigars_[fragment.getReadIndex()],
            0);
    }
    return true;
}

bool DebugStorage::updateMapqStats(const BamTemplate& bamTemplate)
{
    const isaac::alignment::FragmentMetadata &r0 = bamTemplate.getFragmentMetadata(0);
    const unsigned int firstReadNumber = flowcell_.getReadMetadataList().at(r0.getReadIndex()).getNumber();
    if (2 == bamTemplate.getFragmentCount())
    {
        const bool ret1 = r1MapqStatistics_.updateStat(bamTemplate, firstReadNumber, r0, 0);
        const isaac::alignment::FragmentMetadata &r1 = bamTemplate.getFragmentMetadata(1);
        unsigned int secondReadNumber = flowcell_.getReadMetadataList().at(r1.getReadIndex()).getNumber();
        const bool ret2 = r2MapqStatistics_.updateStat(bamTemplate, secondReadNumber, r1, 0);
        const bool pairRet =
            pairMapqStatistics_.updateStat(bamTemplate, firstReadNumber, r0, &r1) ||
            pairMapqStatistics_.updateStat(bamTemplate, secondReadNumber, r1, &r0);
        return ret1 || ret2 || pairRet;
    }

    return 1 == firstReadNumber ?
        r1MapqStatistics_.updateStat(bamTemplate, flowcell_.getReadMetadataList().at(r0.getReadIndex()).getNumber(), r0, 0):
        r2MapqStatistics_.updateStat(bamTemplate, flowcell_.getReadMetadataList().at(r0.getReadIndex()).getNumber(), r0, 0);
}


void DebugStorage::store(
    const BamTemplate &bamTemplate,
    const unsigned barcodeIdx)
{
    updateMapqStats(bamTemplate);
    BamTemplate originalTemplate = bamTemplate;
    isaac::alignment::FragmentMetadata &r1 = originalTemplate.getFragmentMetadata(0);
//    bool store = !r1.isAligned();
    r1QqStatistics_.updateStat(0, bamTemplate, true);
    if (restoreOriginal(flowcell_.getReadMetadataList().at(r1.getReadIndex()).getNumber(), r1))
    {
        r1QqOriginalStatistics_.updateStat(0, originalTemplate, false);
    }
    if (2 == bamTemplate.getFragmentCount())
    {
        isaac::alignment::FragmentMetadata &r2 = originalTemplate.getFragmentMetadata(1);
        r2QqStatistics_.updateStat(1, bamTemplate, true);
        unsigned int r2Number = flowcell_.getReadMetadataList().at(r2.getReadIndex()).getNumber();
//        store = !debugStorage::alignsCorrectly(r2Number, r2);
//        store = store || !r2.isAligned();
        if (restoreOriginal(r2Number, r2))
        {
            r2QqOriginalStatistics_.updateStat(1, originalTemplate, false);
        }
    }

//    if (store)
    {
        actualStorage_.store(originalTemplate, barcodeIdx);
    }
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
