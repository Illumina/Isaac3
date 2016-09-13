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
 ** \file AlignWorkflow.cpp
 **
 ** \brief see AlignWorkflow.hh
 **
 ** \author Come Raczy
 **/

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstring>
#include <cerrno>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include "alignment/MatchSelector.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FileSystem.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/ContigLoader.hh"
#include "reports/AlignmentReportGenerator.hh"
#include "vcf/VcfUtils.hh"
#include "workflow/AlignWorkflow.hh"

namespace isaac
{
namespace workflow
{


struct AllowAllContigFilter
{
    bool isMapped(unsigned, unsigned contigIndex) const {return true;}
};

AlignWorkflow::AlignWorkflow(
    const std::vector<std::string> &argv,
    const std::string &description,
    const std::vector<flowcell::Layout> &flowcellLayoutList,
    const unsigned seedLength,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool cleanupIntermediary,
    const unsigned bclTilesPerChunk,
    const bool ignoreMissingBcls,
    const bool ignoreMissingFilters,
    const unsigned expectedCoverage,
    const uint64_t matchesPerBin,
    const reference::ReferenceMetadataList &referenceMetadataList,
    const bfs::path &tempDirectory,
    const bfs::path &outputDirectory,
    const unsigned int maxThreadCount,
    const unsigned seedBaseQualityMin,
    const unsigned repeatThreshold,
    const int mateDriftRange,
    const unsigned neighborhoodSizeThreshold,
    const uint64_t availableMemory,
    const unsigned clustersAtATimeMax,
    const bool ignoreNeighbors,
    const bool ignoreRepeats,
    const unsigned mapqThreshold,
    const bool perTileTls,
    const bool pfOnly,
    const unsigned baseQualityCutoff,
    const bool keepUnaligned,
    const bool preSortBins,
    const bool preAllocateBins,
    const bool putUnalignedInTheBack,
    const bool realignGapsVigorously,
    const bool realignDodgyFragments,
    const unsigned realignedGapsPerFragment,
    const bool clipSemialigned,
    const bool clipOverlapping,
    const bool scatterRepeats,
    const bool rescueShadows,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const bool noSmithWaterman,
    const bool splitAlignments,
    const int gapMatchScore,
    const int gapMismatchScore,
    const int gapOpenScore,
    const int gapExtendScore,
    const int minGapExtendScore,
    const unsigned splitGapLength,
    const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
    const unsigned inputLoadersMax,
    const unsigned tempSaversMax,
    const unsigned tempLoadersMax,
    const unsigned outputSaversMax,
    const build::GapRealignerMode realignGaps,
    const boost::filesystem::path &knownIndelsPath,
    const int bamGzipLevel,
    const std::string &bamPuFormat,
    const bool bamProduceMd5,
    const std::vector<std::string> &bamHeaderTags,
    const double expectedBgzfCompressionRatio,
    const bool singleLibrarySamples,
    const bool keepDuplicates,
    const bool markDuplicates,
    const bool anchorMate,
    const std::string &binRegexString,
    const common::ScopedMallocBlock::Mode memoryControl,
    const std::vector<std::size_t> &clusterIdList,
    const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
    const reports::AlignmentReportGenerator::ImageFileFormat statsImageFormat,
    const bool bufferBins,
    const bool qScoreBin,
    const boost::array<char, 256> &fullBclQScoreTable,
    const OptionalFeatures optionalFeatures,
    const bool pessimisticMapQ)
    : argv_(argv)
    , description_(description)
    , flowcellLayoutList_(flowcellLayoutList)
    , seedLength_(seedLength)
    , tempDirectory_(tempDirectory)
    , statsDirectory_(outputDirectory/"Stats")
    , reportsDirectory_(outputDirectory/"Reports")
    , projectsDirectory_(outputDirectory/"Projects")
    , matchSelectorStatsXmlPath_(statsDirectory_ / "AlignmentStats.xml")
    , coresMax_(maxThreadCount)
    , seedBaseQualityMin_(seedBaseQualityMin)
    , repeatThreshold_(repeatThreshold)
    , mateDriftRange_(mateDriftRange)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , clusterIdList_(clusterIdList)
    , barcodeMetadataList_(barcodeMetadataList)
    , cleanupIntermediary_(cleanupIntermediary)
    , bclTilesPerChunk_(bclTilesPerChunk)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , ignoreMissingFilters_(ignoreMissingFilters)
    , expectedCoverage_(expectedCoverage)
    , targetBinSize_(matchesPerBin)
    , availableMemory_(availableMemory)
    , clustersAtATimeMax_(clustersAtATimeMax)
    , mapqThreshold_(mapqThreshold)
    , perTileTls_(perTileTls)
    , pfOnly_(pfOnly)
    , baseQualityCutoff_(baseQualityCutoff)
    , keepUnaligned_(keepUnaligned)
    , preSortBins_(preSortBins)
    , preAllocateBins_(preAllocateBins)
    , putUnalignedInTheBack_(putUnalignedInTheBack)
    , realignGapsVigorously_(realignGapsVigorously)
    , realignDodgyFragments_(realignDodgyFragments)
    , realignedGapsPerFragment_(realignedGapsPerFragment)
    , clipSemialigned_(clipSemialigned)
    , clipOverlapping_(clipOverlapping)
    , scatterRepeats_(scatterRepeats)
    , rescueShadows_(rescueShadows)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , smitWatermanGapsMax_(smitWatermanGapsMax)
    , smartSmithWaterman_(smartSmithWaterman)
    , noSmithWaterman_(noSmithWaterman)
    , splitAlignments_(splitAlignments)
    , gapMatchScore_(gapMatchScore)
    , gapMismatchScore_(gapMismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , minGapExtendScore_(minGapExtendScore)
    , splitGapLength_(splitGapLength)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , tempLoadersMax_(tempLoadersMax)
    , outputSaversMax_(outputSaversMax)
    , realignGaps_(realignGaps)
    , knownIndelsPath_(knownIndelsPath)
    , bamGzipLevel_(bamGzipLevel)
    , bamPuFormat_(bamPuFormat)
    , bamProduceMd5_(bamProduceMd5)
    , bamHeaderTags_(bamHeaderTags)
    , expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio)
    , singleLibrarySamples_(singleLibrarySamples)
    , keepDuplicates_(keepDuplicates)
    , markDuplicates_(markDuplicates)
    , anchorMate_(anchorMate)
    , bufferBins_(bufferBins)
    , qScoreBin_(qScoreBin)
    , fullBclQScoreTable_(fullBclQScoreTable)
    , optionalFeatures_(optionalFeatures)
    , pessimisticMapQ_(pessimisticMapQ)
    , binRegexString_(binRegexString)
    , memoryControl_(memoryControl)
    , userTemplateLengthStatistics_(userTemplateLengthStatistics)
    , demultiplexingStatsXmlPath_(statsDirectory_ / "DemultiplexingStats.xml")
    , statsImageFormat_(statsImageFormat)
    , referenceMetadataList_(referenceMetadataList)
    , sortedReferenceMetadataList_(loadSortedReferenceXml(seedLength, referenceMetadataList))
    , contigLists_(reference::loadContigs(sortedReferenceMetadataList_, AllowAllContigFilter(), common::ThreadVector(inputLoadersMax_)))
    , kUniquenessAnnotations_(reference::loadAnnotations(sortedReferenceMetadataList_))
    , state_(Start)
      // dummy initialization. Will be replaced with real object once match finding is over
    , foundMatchesMetadata_(tempDirectory_, barcodeMetadataList_, 0, sortedReferenceMetadataList_)
    , barcodeTemplateLengthStatistics_(barcodeMetadataList_.size())
{
    const std::vector<bfs::path> createList = boost::assign::list_of
        (tempDirectory_)(outputDirectory)(statsDirectory_)(reportsDirectory_)(projectsDirectory_);
    common::createDirectories(createList);

    BOOST_FOREACH(const flowcell::Layout &layout, flowcellLayoutList)
    {
        ISAAC_THREAD_CERR << "Aligner: adding base-calls path " << layout.getBaseCallsPath() << std::endl;
    }
}

reference::SortedReferenceMetadataList AlignWorkflow::loadSortedReferenceXml(
    const unsigned seedLength,
    const reference::ReferenceMetadataList &referenceMetadataList)
{
    reference::SortedReferenceMetadataList ret(referenceMetadataList.size());
    BOOST_FOREACH(const reference::ReferenceMetadata &reference, referenceMetadataList)
    {
        const unsigned referenceIndex = &reference - &referenceMetadataList.front();
        reference::SortedReferenceMetadata &ref = ret.at(referenceIndex);
        ref = reference::loadSortedReferenceXml(reference.getXmlPath());
    }
    return ret;
}

void AlignWorkflow::findMatches(
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics) const
{
    alignWorkflow::FindHashMatchesTransition findMatchesTransition(
        flowcellLayoutList_,
        barcodeMetadataList_,
        cleanupIntermediary_,
        bclTilesPerChunk_,
        ignoreMissingBcls_,
        ignoreMissingFilters_,
        availableMemory_,
        clustersAtATimeMax_,
        tempDirectory_,
        demultiplexingStatsXmlPath_,
        coresMax_,
        seedBaseQualityMin_,
        repeatThreshold_,
        neighborhoodSizeThreshold_,
        ignoreNeighbors_,
        ignoreRepeats_,
        inputLoadersMax_,
        tempSaversMax_,
        memoryControl_,
        clusterIdList_,
        sortedReferenceMetadataList_,
        contigLists_,
        kUniquenessAnnotations_,
        optionalFeatures_ & BamZX,
        mateDriftRange_,
        userTemplateLengthStatistics_, mapqThreshold_, perTileTls_, pfOnly_,
        reports::AlignmentReportGenerator::none != statsImageFormat_,
        baseQualityCutoff_,
        keepUnaligned_, clipSemialigned_, clipOverlapping_,
        scatterRepeats_, rescueShadows_, anchorMate_, gappedMismatchesMax_, smitWatermanGapsMax_, smartSmithWaterman_, noSmithWaterman_, splitAlignments_,
        gapMatchScore_, gapMismatchScore_, gapOpenScore_, gapExtendScore_, minGapExtendScore_, splitGapLength_,
        dodgyAlignmentScore_,
        qScoreBin_,
        fullBclQScoreTable_,
        bufferBins_,
        expectedCoverage_,
        targetBinSize_,
        expectedBgzfCompressionRatio_,
        preSortBins_,
        preAllocateBins_,
        binRegexString_);

    if (16 == seedLength_)
    {
        findMatchesTransition.perform<oligo::ShortKmerType>(
            foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath_);
    }
//    else if (28 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::BasicKmerType<28> >(foundMatches);
//    }
//    else if (30 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::BasicKmerType<30> >(foundMatches);
//    }
//    else if (32 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::KmerType>(foundMatches);
//    }
//    else if (34 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::BasicKmerType<34> >(foundMatches);
//    }
//    else if (36 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::BasicKmerType<36> >(foundMatches);
//    }
//    else if (64 == seedLength_)
//    {
//        findMatchesTransition.perform<oligo::LongKmerType>(foundMatches);
//    }
    else
    {
        ISAAC_ASSERT_MSG(false, "Unexpected seed length " << seedLength_);
    }
}

void AlignWorkflow::cleanupBins() const
{
    ISAAC_THREAD_CERR << "Removing intermediary bin files" << std::endl;
    unsigned removed = 0;
    BOOST_FOREACH(const alignment::BinMetadata &bin, selectedMatchesMetadata_)
    {
        removed += boost::filesystem::remove(bin.getPath());
    }
    ISAAC_THREAD_CERR << "Removing intermediary bin files done. " << removed << " files removed." << std::endl;
}

void AlignWorkflow::generateAlignmentReports() const
{
    ISAAC_THREAD_CERR << "Generating the match selector reports from " << matchSelectorStatsXmlPath_ << std::endl;
    reports::AlignmentReportGenerator reportGenerator(flowcellLayoutList_, barcodeMetadataList_,
                                                  matchSelectorStatsXmlPath_, demultiplexingStatsXmlPath_,
                                                  tempDirectory_, reportsDirectory_,
                                                  statsImageFormat_);
    reportGenerator.run();
    ISAAC_THREAD_CERR << "Generating the match selector reports done from " << matchSelectorStatsXmlPath_ << std::endl;
}

const build::BarcodeBamMapping AlignWorkflow::generateBam(
    const SelectedMatchesMetadata &binPaths,
    const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics) const
{
    ISAAC_THREAD_CERR << "Generating the BAM files" << std::endl;

    build::Build build(argv_, description_,
                       flowcellLayoutList_, foundMatchesMetadata_.tileMetadataList_, barcodeMetadataList_,
                       binPaths,
                       referenceMetadataList_,
                       barcodeTemplateLengthStatistics,
                       sortedReferenceMetadataList_,
                       contigLists_,
                       kUniquenessAnnotations_,
                       projectsDirectory_,
                       tempLoadersMax_, coresMax_, outputSaversMax_, realignGaps_, knownIndelsPath_,
                       bamGzipLevel_, bamPuFormat_, bamProduceMd5_, bamHeaderTags_, expectedCoverage_, expectedBgzfCompressionRatio_, singleLibrarySamples_,
                       keepDuplicates_, markDuplicates_, anchorMate_,
                       realignGapsVigorously_, realignDodgyFragments_, realignedGapsPerFragment_,
                       clipSemialigned_, 
                       // when splitting reads, the bin regex cannot be used to decide which 
                       // contigs to load.
                       splitAlignments_, binRegexString_,
                       alignment::TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_ ?
                           0 : boost::numeric_cast<unsigned char>(dodgyAlignmentScore_),
                       keepUnaligned_, putUnalignedInTheBack_,
                       build::IncludeTags(
                           optionalFeatures_ & BamAS,
                           optionalFeatures_ & BamBC,
                           optionalFeatures_ & BamNM,
                           optionalFeatures_ & BamOC,
                           optionalFeatures_ & BamRG,
                           optionalFeatures_ & BamSM,
                           optionalFeatures_ & BamZX,
                           optionalFeatures_ & BamZY),
                       pessimisticMapQ_, splitGapLength_);
    {
        common::ScopedMallocBlock  mallocBlock(memoryControl_);
        build.run(mallocBlock);
    }
    build.dumpStats(statsDirectory_ / "BuildStats.xml");
    ISAAC_THREAD_CERR << "Generating the BAM files done" << std::endl;
    return build.getBarcodeBamMapping();
}

void AlignWorkflow::run()
{
    ISAAC_ASSERT_MSG(Start == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(AlignDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(AlignmentReportsDone == state_, "Unexpected state");
    step();
    ISAAC_ASSERT_MSG(BamDone == state_, "Unexpected state");
}

AlignWorkflow::State AlignWorkflow::getNextState() const
{
    switch (state_)
    {
    case Start:
    {
        return AlignDone;
    }
    case AlignDone:
    {
        return AlignmentReportsDone;
    }
    case AlignmentReportsDone:
    {
        return BamDone;
    }
    case Finish:
    {
        return Finish;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state value");
        return Invalid;
    }
    }
}

AlignWorkflow::State AlignWorkflow::step()
{
    using std::swap;
    switch (state_)
    {
    case Start:
    {
        findMatches(foundMatchesMetadata_, selectedMatchesMetadata_, barcodeTemplateLengthStatistics_);
        state_ = getNextState();
        break;
    }
    case AlignDone:
    {
        generateAlignmentReports();
        state_ = getNextState();
        break;
    }
    case AlignmentReportsDone:
    {
        barcodeBamMapping_ = generateBam(selectedMatchesMetadata_, barcodeTemplateLengthStatistics_);
        state_ = getNextState();
        break;
    }
    case Finish:
    {
        ISAAC_THREAD_CERR << "Already at the Finish state" << std::endl;
        break;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state");
        break;
    }
    }
    return state_;
}

void AlignWorkflow::cleanupIntermediary()
{
    switch (state_)
    {
    case Finish:
    {
        cleanupBins();
        //fall through
    }
    case AlignmentReportsDone:
    case AlignDone:
    case Start:
    {
        break;
    }
    default:
    {
        ISAAC_ASSERT_MSG(false, "Invalid state " << state_);
        break;
    }
    }
}

/**
 * \return the initial state from which the rewind occurred
 */
AlignWorkflow::State AlignWorkflow::rewind(AlignWorkflow::State to)
{
    AlignWorkflow::State ret = state_;
    switch (to)
    {
    case Last:
    {
        // Nothing to do. We're at the Last state by definition
        break;
    }
    case Start:
    {
        // Start is always possible
        state_ = Start;
        break;
    }
    case AlignDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to MatchFinderDone is not possible"));}
        state_ = AlignDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to MatchFinderDone successful" << std::endl;
        break;
    }
    case AlignmentReportsDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to MatchSelectorReportsDone is not possible"));}
        if (AlignDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to MatchSelectorReportsDone is not possible"));}
        state_ = AlignmentReportsDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to AlignmentReportsDone successful" << std::endl;
        break;
    }
    case BamDone:
    {
        if (Start == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from Start to BamDone is not possible"));}
        if (AlignDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchFinderDone to BamDone is not possible"));}
        if (AlignmentReportsDone == state_) {BOOST_THROW_EXCEPTION(common::PreConditionException("Aligner rewind from MatchSelectorDone to BamDone is not possible"));}
        state_ = BamDone;
        ISAAC_THREAD_CERR << "Workflow state rewind to BamDone successful" << std::endl;
        break;
    }
    default:
    {
        assert(false);
        break;
    }
    }

    return ret;
}


} // namespace workflow
} // namespace isaac
