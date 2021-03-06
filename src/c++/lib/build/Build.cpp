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
 ** \file Build.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 **
 ** \author Roman Petrovski
 **/

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA
 
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/regex.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "bgzf/BgzfCompressor.hh"
#include "build/Build.hh"
#include "build/IndelLoader.hh"
#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Threads.hpp"
#include "io/Fragment.hh"
#include "reference/ContigLoader.hh"

#include "BuildStatsXml.hh"
#include "SortedReferenceXmlBamHeaderAdapter.hh"

namespace isaac
{
namespace build
{

const unsigned BuildContigMap::UNMAPPED_CONTIG;
/**
 * \return Returns the total memory in bytes required to load the bin data and indexes
 */
static uint64_t getBinTotalSize(const alignment::BinMetadata & binMetadata)
{
    return
        binMetadata.getDataSize() +
        binMetadata.getFIdxElements() * sizeof(FStrandFragmentIndex) +
        binMetadata.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
        binMetadata.getSeIdxElements() * sizeof(SeFragmentIndex);
}

uint64_t Build::estimateBinCompressedDataRequirements(
    const alignment::BinMetadata & binMetadata,
    const unsigned outputFileIndex) const
{
    // TODO: put the real number in here.
    static const uint64_t EMPTY_BGZF_BLOCK_SIZE = 1234UL;
    if (!binMetadata.getTotalElements())
    {
        return EMPTY_BGZF_BLOCK_SIZE;
    }

    uint64_t thisOutputFileBarcodeElements = 0;
    // accumulate size required to store all barcodes that map to the same output file
    BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
    {
        const unsigned barcodeOutputFileIndex = barcodeBamMapping_.getSampleIndex(barcode.getIndex());
        if (outputFileIndex == barcodeOutputFileIndex)
        {
            const unsigned barcodeIndex = barcode.getIndex();
            ISAAC_ASSERT_MSG(0 != binMetadata.getTotalElements() || 0 == binMetadata.getBarcodeElements(barcodeIndex), "Can't have empty bin with non-empty bin barcode");

            thisOutputFileBarcodeElements += binMetadata.getBarcodeElements(barcodeIndex);
        }
    }

    // assume all data will take the same fraction or less than the number derived from demultiplexed fragments.
    return EMPTY_BGZF_BLOCK_SIZE +
        ((getBinTotalSize(binMetadata) * thisOutputFileBarcodeElements +
            binMetadata.getTotalElements() - 1) / binMetadata.getTotalElements()) * expectedBgzfCompressionRatio_;
}

inline boost::filesystem::path getSampleBamPath(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadata &barcode)
{
    return outputDirectory / barcode.getProject() / barcode.getSampleName() / "sorted.bam";
}
/**
 * \brief Produces mapping so that all barcodes having the same sample name go into the same output file.
 *
 * \return Returns the pair of a vector of that maps a barcode index to a unique output file index in the
 *         second vector so that two barcodes that are supposed to go into the same file will end up
 *         having the same mapping
 */
BarcodeBamMapping mapBarcodesToFiles(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    // Map barcodes to projects
    std::vector<std::string> projects;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(projects),
                   boost::bind(&flowcell::BarcodeMetadata::getProject, _1));
    std::sort(projects.begin(), projects.end());
    projects.erase(std::unique(projects.begin(), projects.end()), projects.end());

    BarcodeBamMapping::BarcodeProjectIndexMap barcodeProject(barcodeMetadataList.size());
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        barcodeProject.at(barcode.getIndex()) =
            std::distance(projects.begin(), std::lower_bound(projects.begin(), projects.end(), barcode.getProject()));
    }

    // Map barcodes to sample paths
    std::vector<boost::filesystem::path> samples;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(samples),
                   boost::bind(&getSampleBamPath, outputDirectory, _1));
    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

    BarcodeBamMapping::BarcodeProjectIndexMap barcodeSample(barcodeMetadataList.size());
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        barcodeSample.at(barcode.getIndex()) =
            std::distance(samples.begin(), std::lower_bound(samples.begin(), samples.end(),
                                                            getSampleBamPath(outputDirectory, barcode)));
    }

    return BarcodeBamMapping(barcodeProject, barcodeSample, samples);
}

inline bool orderBySampleIndex(
    const BarcodeBamMapping &barcodeBamMapping,
    const flowcell::BarcodeMetadata &left,
    const flowcell::BarcodeMetadata &right)
{
    return barcodeBamMapping.getSampleIndex(left.getIndex()) < barcodeBamMapping.getSampleIndex(right.getIndex());
}

std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > Build::createOutputFileStreams(
    const flowcell::TileMetadataList &tileMetadataList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    boost::ptr_vector<bam::BamIndex> &bamIndexes) const
{
    unsigned sinkIndexToCreate = 0;
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > ret;
    ret.reserve(barcodeBamMapping_.getTotalSamples());

    std::vector<boost::filesystem::path> directories;
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
        directories.push_back(bamPath.parent_path());
        directories.push_back(directories.back().parent_path());
    }
    common::createDirectories(directories);

    flowcell::BarcodeMetadataList barcodesOrderedBySample(barcodeMetadataList);
    std::sort(barcodesOrderedBySample.begin(), barcodesOrderedBySample.end(),
              boost::bind(&orderBySampleIndex, boost::ref(barcodeBamMapping_), _1, _2));
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodesOrderedBySample)
    {
        if (sinkIndexToCreate == barcodeBamMapping_.getSampleIndex(barcode.getIndex()))
        {
            const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
            if (!barcode.isUnmappedReference())
            {
                ISAAC_THREAD_CERR << "Created BAM file: " << bamPath << std::endl;

                const reference::SortedReferenceMetadata &sampleReference =
                    sortedReferenceMetadataList_.at(barcode.getReferenceIndex());

                std::string compressedHeader;
                {
                    std::ostringstream oss(compressedHeader);
                    boost::iostreams::filtering_ostream bgzfStream;
                    bgzfStream.push(bgzf::BgzfCompressor(bamGzipLevel_));
                    bgzfStream.push(oss);
                    bam::serializeHeader(bgzfStream,
                                         argv_,
                                         description_,
                                         bamHeaderTags_,
                                         bamPuFormat_,
                                         makeSortedReferenceXmlBamHeaderAdapter(
                                             sampleReference,
                                             boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1),
                                             tileMetadataList, barcodeMetadataList,
                                             barcode.getSampleName()));
                    bgzfStream.strict_sync();
                    compressedHeader = oss.str();
                }

                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>(new boost::iostreams::filtering_ostream()));
                boost::iostreams::filtering_ostream &bamStream = *ret.back();
                if (bamProduceMd5_)
                {
                    bamStream.push(io::FileSinkWithMd5(bamPath.c_str(), std::ios_base::binary));
                }
                else
                {
                    bamStream.push(boost::iostreams::basic_file_sink<char>(bamPath.c_str(), std::ios_base::binary));
                }

                if (!bamStream) {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
                }

                if (!bamStream.write(compressedHeader.c_str(), compressedHeader.size()))
                {
                    BOOST_THROW_EXCEPTION(
                        common::IoException(errno, (boost::format("Failed to write %d bytes into stream %s") %
                            compressedHeader.size() % bamPath.string()).str()));
                }

                // Create BAM Indexer
                unsigned headerCompressedLength = compressedHeader.size();
                unsigned contigCount = sampleReference.getContigsCount(
                    boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1));
                bamIndexes.push_back(new bam::BamIndex(bamPath, contigCount, headerCompressedLength));
            }
            else
            {
                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>());
                bamIndexes.push_back(new bam::BamIndex());
                ISAAC_THREAD_CERR << "Skipped BAM file due to unmapped barcode reference: " << bamPath << " " << barcode << std::endl;
            }
            ++sinkIndexToCreate;
        }
    }
    ISAAC_ASSERT_MSG(barcodeBamMapping_.getTotalSamples() == sinkIndexToCreate, "must create all output file sinks");

    return ret;
}

const alignment::BinMetadataCRefList filterBins(
    const alignment::BinMetadataList& bins,
    const std::string &binRegexString)
{
    alignment::BinMetadataCRefList ret;

    if ("all" == binRegexString)
    {
        std::transform(bins.begin(), bins.end(), std::back_inserter(ret),
                       [](const alignment::BinMetadata &bm){return boost::cref(bm);});
    }
    else if ("skip-empty" == binRegexString)
    {
        BOOST_FOREACH(const alignment::BinMetadata &binMetadata, bins)
        {
            if (!binMetadata.isEmpty())
            {
                ret.push_back(boost::cref(binMetadata));
            }
        }
    }
    else // use regex to filter bins by name
    {
        std::string regexString(binRegexString);
        std::replace(regexString.begin(), regexString.end(), ',', '|');
        boost::regex re(regexString);
        BOOST_FOREACH(const alignment::BinMetadata &bin, bins)
        {
            if (!bin.isEmpty() && boost::regex_search(bin.getPath().filename().string(), re))
            {
                ret.push_back(boost::ref(bin));
            }
        }
        if (ret.empty())
        {
            ISAAC_THREAD_CERR << "WARNING: Bam files will be empty. No bins are left after applying the following regex filter: "
                << regexString << std::endl;
        }
    }
    return ret;
}

static void breakUpBin(
    const alignment::BinMetadata& bin,
    const uint64_t partsCount,
    alignment::BinMetadataList &ret)
{
    ISAAC_ASSERT_MSG(0 == bin.getIndex(), "At the moment only bin 0 is expected to be unaligned");
    const uint64_t newBinSize = bin.getDataSize() / partsCount;
    ISAAC_THREAD_CERR << "Breaking unaligned bin of " << bin.getDataSize()/1024/1024 << " megabytes into " <<
        partsCount << " bins of " << newBinSize /1024/1024 << " megabytes for better parallelization: " << bin << std::endl;

    for (uint64_t offset = 0; bin.getDataSize() > offset;)
    {
        alignment::BinMetadata part = bin.getChunks(offset, newBinSize);
        ISAAC_THREAD_CERR << " offset:" << offset << " " << part <<std::endl;
        offset += part.getDataSize();
        ret.push_back(part);
    }
}

/**
 * \brief breaks unaligned bin into about partsCount of roughly equivalent size bins
 */
static alignment::BinMetadataCRefList breakUpUnalignedBin(
    alignment::BinMetadataCRefList bins,
    const uint64_t partsCount,
    const bool keepUnaligned,
    const bool putUnalignedInTheBack,
    alignment::BinMetadataList &unalignedBinParts)
{
    if (!bins.empty())
    {
        // unaligned bins must occur at the start of the list
        const alignment::BinMetadata &bin = bins.front();
        if (bin.isUnalignedBin())
        {
            if (keepUnaligned)
            {
                breakUpBin(bin, partsCount, unalignedBinParts);
            }
            bins.erase(bins.begin());
        }

        if (putUnalignedInTheBack)
        {
            std::transform(unalignedBinParts.begin(), unalignedBinParts.end(), std::back_inserter(bins),
                           [](const alignment::BinMetadata &bm){return boost::cref(bm);});
        }
        else
        {
            std::transform(unalignedBinParts.begin(), unalignedBinParts.end(), std::inserter(bins, bins.begin()),
                           [](const alignment::BinMetadata &bm){return boost::cref(bm);});
        }
    }
    return bins;
}

Build::Build(const std::vector<std::string> &argv,
             const std::string &description,
             const flowcell::FlowcellLayoutList &flowcellLayoutList,
             const flowcell::TileMetadataList &tileMetadataList,
             const flowcell::BarcodeMetadataList &barcodeMetadataList,
             const alignment::BinMetadataList &bins,
             const reference::ReferenceMetadataList &referenceMetadataList,
             const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
             const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
             const reference::NumaContigLists &contigLists,
             const isaac::reference::NumaContigAnnotationsList &kUniquenessAnnotations,
             const boost::filesystem::path outputDirectory,
             const unsigned maxLoaders,
             const unsigned maxComputers,
             const unsigned maxSavers,
             const build::GapRealignerMode realignGaps,
             const boost::filesystem::path &knownIndelsPath,
             const int bamGzipLevel,
             const std::string &bamPuFormat,
             const bool bamProduceMd5,
             const std::vector<std::string> &bamHeaderTags,
             const unsigned expectedCoverage,
             const double expectedBgzfCompressionRatio,
             const bool singleLibrarySamples,
             const bool keepDuplicates,
             const bool markDuplicates,
             const bool anchorMate,
             const bool realignGapsVigorously,
             const bool realignDodgyFragments,
             const unsigned realignedGapsPerFragment,
             const bool clipSemialigned,
             const bool loadAllContigs,
             const std::string &binRegexString,
             const unsigned char forcedDodgyAlignmentScore,
             const bool keepUnaligned,
             const bool putUnalignedInTheBack,
             const IncludeTags includeTags,
             const bool pessimisticMapQ,
             const unsigned splitGapLength)
    :argv_(argv),
     description_(description),
     flowcellLayoutList_(flowcellLayoutList),
     tileMetadataList_(tileMetadataList),
     barcodeMetadataList_(barcodeMetadataList),
     unalignedBinParts_(),
     bins_(breakUpUnalignedBin(
         filterBins(bins, binRegexString), maxComputers, keepUnaligned, putUnalignedInTheBack, unalignedBinParts_)),
     sortedReferenceMetadataList_(sortedReferenceMetadataList),
     contigMap_(barcodeMetadataList_, bins_, sortedReferenceMetadataList_, false),//!loadAllContigs && "skip-empty" == binRegexString),
     outputDirectory_(outputDirectory),
     maxLoaders_(maxLoaders),
     maxComputers_(maxComputers),
     maxRealigners_(1),
     allocatedBins_(0),
     maxSavers_(maxSavers),
     bamGzipLevel_(bamGzipLevel),
     bamPuFormat_(bamPuFormat),
     bamProduceMd5_(bamProduceMd5),
     bamHeaderTags_(bamHeaderTags),
     forcedDodgyAlignmentScore_(forcedDodgyAlignmentScore),
     singleLibrarySamples_(singleLibrarySamples),
     keepDuplicates_(keepDuplicates),
     markDuplicates_(markDuplicates),
     anchorMate_(anchorMate),
     realignGapsVigorously_(realignGapsVigorously),
     realignDodgyFragments_(realignDodgyFragments),
     realignedGapsPerFragment_(realignedGapsPerFragment),
     clipSemialigned_(clipSemialigned),
     realignGaps_(realignGaps),
     expectedCoverage_(expectedCoverage),
     expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio),
     maxReadLength_(getMaxReadLength(flowcellLayoutList_)),
     includeTags_(includeTags),
     pessimisticMapQ_(pessimisticMapQ),
     splitGapLength_(splitGapLength),
     forceTermination_(false),
     threads_(maxComputers_ + maxLoaders_ + maxSavers_),
     contigLists_(contigLists),
     barcodeBamMapping_(mapBarcodesToFiles(outputDirectory_, barcodeMetadataList_)),
     bamIndexes_(),
     bamFileStreams_(createOutputFileStreams(tileMetadataList_, barcodeMetadataList_, bamIndexes_)),
     stats_(bins_, barcodeMetadataList_),
     threadBgzfBuffers_(threads_.size(), BgzfBuffers(bamFileStreams_.size())),
     threadBgzfStreams_(threads_.size()),
     threadBamIndexParts_(threads_.size()),
     knownIndels_((build::GapRealignerMode::REALIGN_NONE == realignGaps_ || knownIndelsPath.empty()) ?
         gapRealigner::Gaps() : loadIndels(knownIndelsPath, sortedReferenceMetadataList_)),
     gapRealigner_(threads_.size(),
         realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, clipSemialigned,
         barcodeMetadataList, barcodeTemplateLengthStatistics, contigLists_),
     binSorter_(singleLibrarySamples_, keepDuplicates_, markDuplicates_, anchorMate_,
               barcodeBamMapping_, barcodeMetadataList_, contigLists_, splitGapLength_, kUniquenessAnnotations)
{
    computeSlotWaitingBins_.reserve(threads_.size());
    while(threadBgzfStreams_.size() < threads_.size())
    {
        threadBgzfStreams_.push_back(new boost::ptr_vector<boost::iostreams::filtering_ostream>(bamFileStreams_.size()));
    }
    while(threadBamIndexParts_.size() < threads_.size())
    {
        threadBamIndexParts_.push_back(new boost::ptr_vector<bam::BamIndexPart>(bamFileStreams_.size()));
    }

    threads_.execute(boost::bind(&Build::allocateThreadData, this, _1));

    // when number of bins is smaller than number of threads, some complete tasks don't get erased before other tasks are added.
    tasks_.reserve(std::max(bins_.size(), threads_.size()));

//    testBinsFitInRam();
}

void Build::allocateThreadData(const std::size_t threadNumber)
{
//#ifdef HAVE_NUMA
//    if (common::isNumaAvailable())
//    {
//        uint64_t nodemask = 1UL << common::ThreadVector::getThreadNumaNode();
//        ISAAC_ASSERT_MSG(-1 != set_mempolicy(MPOL_BIND/*|MPOL_F_STATIC_NODES*/, &nodemask, sizeof(nodemask) * 8),
//                         "set_mempolicy for nodemask: " << nodemask <<
//                         " failed, errno: " << errno << ":" << strerror(errno));
//    }
//#endif //HAVE_NUMA
}

void Build::run(common::ScopedMallocBlock &mallocBlock)
{
    alignment::BinMetadataCRefList::const_iterator nextUnprocessedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnallocatedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnloadedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUncompressedBinIt(bins_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnsavedBinIt(bins_.begin());

    threads_.execute(boost::bind(&Build::sortBinParallel, this,
                                boost::ref(nextUnprocessedBinIt),
                                boost::ref(nextUnallocatedBinIt),
                                boost::ref(nextUnloadedBinIt),
                                boost::ref(nextUncompressedBinIt),
                                boost::ref(nextUnsavedBinIt),
                                boost::ref(mallocBlock),
                                _1));

    unsigned fileIndex = 0;
    BOOST_FOREACH(const boost::filesystem::path &bamFilePath, barcodeBamMapping_.getPaths())
    {
        // some of the streams are null_sink (that's when reference is unmapped for the sample).
        // this is the simplest way to ignore them...
        std::ostream *stm = bamFileStreams_.at(fileIndex).get();
        if (stm)
        {
            bam::serializeBgzfFooter(*stm);
            stm->flush();
            ISAAC_THREAD_CERR << "BAM file generated: " << bamFilePath.c_str() << "\n";
            bamIndexes_.at(fileIndex).flush();
            ISAAC_THREAD_CERR << "BAM index generated for " << bamFilePath.c_str() << "\n";
        }
        ++fileIndex;
    }
}

void Build::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    BuildStatsXml statsXml(sortedReferenceMetadataList_, bins_, barcodeMetadataList_, stats_);
    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }
    statsXml.serialize(os);
}

uint64_t Build::estimateOptimumFragmentsPerBin(
    const unsigned int estimatedFragmentSize,
    const uint64_t availableMemory,
    const double expectedBgzfCompressionRatio,
    const unsigned computeThreads)
{
//    const size_t maxFragmentIndexBytes = std::max(sizeof(io::RStrandOrShadowFragmentIndex),
//                                                         sizeof(io::FStrandFragmentIndex));

    const std::size_t maxFragmentDedupedIndexBytes = sizeof(PackedFragmentBuffer::Index);
    const std::size_t maxFragmentCompressedBytes = estimatedFragmentSize * expectedBgzfCompressionRatio;

    const std::size_t fragmentMemoryRequirements =
        // assume the initial indexes don't stay in memory for too long //maxFragmentIndexBytes              //index containing duplicates
        + estimatedFragmentSize                 //data
        + maxFragmentDedupedIndexBytes     //deduplicated index
        + maxFragmentCompressedBytes       //bgzf chunk
        ;

    // reasonable amount of bins-in-progress to allow for no-delay input/compute/output overlap
//    const unsigned minOverlap = 3;;
    // try to increase granularity so that the CPU gets efficiently utilized.
    const unsigned minOverlap = computeThreads;
    return availableMemory / fragmentMemoryRequirements / minOverlap;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
void Build::reserveBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber,
    boost::shared_ptr<BinData> &binDataPtr)
{
    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(threadNumber);
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts = threadBamIndexParts_.at(threadNumber);
    const alignment::BinMetadata &bin = *thisThreadBinIt;
    // bin stats have an entry per filtered bin reference.
    const unsigned binStatsIndex = std::distance(bins_.begin(), thisThreadBinIt);
    common::ScopedMallocBlockUnblock unblockMalloc(mallocBlock);
    reserveBuffers(
        bin, binStatsIndex, contigLists_.threadNodeContainer(), bgzfStreams, bamIndexParts,
        threadBgzfBuffers_.at(threadNumber), binDataPtr);
}

void Build::cleanupBinAllocationFailure(
    const alignment::BinMetadata& bin,
    boost::ptr_vector<boost::iostreams::filtering_ostream>& bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart>& bamIndexParts,
    boost::shared_ptr<BinData>& binDataPtr, BgzfBuffers& bgzfBuffers)
{
    bgzfStreams.clear();
    bamIndexParts.clear();
    // give a chance other threads to allocate what they need... TODO: this is not required anymore as allocation happens orderly
    binDataPtr.reset();
    for(bam::BgzfBuffer &bgzfBuffer : bgzfBuffers)
    {
        bam::BgzfBuffer().swap(bgzfBuffer);
    }
    // reset errno, to prevent misleading error messages when failing code does not set errno
    errno = 0;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
void Build::reserveBuffers(
    const alignment::BinMetadata &bin,
    const unsigned binStatsIndex,
    const reference::ContigLists &contigLists,
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
    BgzfBuffers &bgzfBuffers,
    boost::shared_ptr<BinData> &binDataPtr)
{
    try
    {
        binDataPtr = boost::shared_ptr<BinData>(
            new BinData(realignedGapsPerFragment_,
                        barcodeBamMapping_, barcodeMetadataList_,
                        realignGaps_, knownIndels_, bin, binStatsIndex, tileMetadataList_, contigMap_, contigLists, maxReadLength_,
                        forcedDodgyAlignmentScore_,  flowcellLayoutList_, includeTags_, pessimisticMapQ_, splitGapLength_,
                        expectedCoverage_));

        unsigned outputFileIndex = 0;
        for(bam::BgzfBuffer &bgzfBuffer : bgzfBuffers)
        {
            bgzfBuffer.reserve(estimateBinCompressedDataRequirements(bin, outputFileIndex++));
        }

        ISAAC_ASSERT_MSG(!bgzfStreams.size(), "Expecting empty pool of streams");
        while(bgzfStreams.size() < bamFileStreams_.size())
        {
            bgzfStreams.push_back(new boost::iostreams::filtering_ostream);
            bgzfStreams.back().push(bgzf::BgzfCompressor(bamGzipLevel_));
            bgzfStreams.back().push(
                boost::iostreams::back_insert_device<bam::BgzfBuffer >(
                    bgzfBuffers.at(bgzfStreams.size()-1)));
            bgzfStreams.back().exceptions(std::ios_base::badbit);
        }

        ISAAC_ASSERT_MSG(!bamIndexParts.size(), "Expecting empty pool of bam index parts");
        while(bamIndexParts.size() < bamFileStreams_.size())
        {
            bamIndexParts.push_back(new bam::BamIndexPart);
        }
    }
    catch (...)
    {
        cleanupBinAllocationFailure(bin, bgzfStreams, bamIndexParts, binDataPtr, bgzfBuffers);
        throw;
    }
}

template <typename ExceptionType>
const char *getExceptionName(ExceptionType &e)
{
    try
    {
        throw e;
    }
    catch (std::bad_alloc &)
    {
        return "std::bad_alloc";
    }
    catch (boost::iostreams::zlib_error &)
    {
        return "boost::iostreams::zlib_error";
    }
    catch (common::IoException &)
    {
        return "common::IoException";
    }
    catch (...)
    {
        return "Unknown Exception";
    }

}

template <typename ExceptionType, typename ExceptionDataT>
bool Build::handleBinAllocationFailure(
    bool warningTraced,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const ExceptionType &e,
    const ExceptionDataT &errorData)
{
    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    if (!allocatedBins_)
    {
        forceTermination_ = true;
        stateChangedCondition_.notify_all();
        // couldn't allocate our bin and not bins currently allocated. No way we will ever be able to allocate it
        BOOST_THROW_EXCEPTION(common::ThreadingException(
            (boost::format("ERROR: Failing due to: %s blocking everything with %s : %s Error data: %s ")
                        % *thisThreadBinIt % getExceptionName(e) % e.what() % errorData).str()));
    }

    if (!warningTraced)
    {
        ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
            *thisThreadBinIt << " until " << e.what() << " clears. Error data: " << errorData << std::endl;
        warningTraced = true;
    }

    return warningTraced;
}

boost::shared_ptr<BinData> Build::allocateBin(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber)
{
    bool warningTraced = false;

    boost::shared_ptr<BinData> ret;
    while(nextUnallocatedBinIt != thisThreadBinIt)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        stateChangedCondition_.wait(lock);
    }

    while(true)
    {
        try
        {
            reserveBuffers(lock, thisThreadBinIt, binsEnd, mallocBlock, threadNumber, ret);
            break;
        }
        catch (std::bad_alloc &a)
        {
            uint64_t totalBuffersNeeded = 0UL;
            for(unsigned outputFileIndex = 0; outputFileIndex < threadBgzfStreams_.at(threadNumber).size(); ++outputFileIndex)
            {
                totalBuffersNeeded += estimateBinCompressedDataRequirements(*thisThreadBinIt, outputFileIndex++);
            }
            warningTraced = handleBinAllocationFailure(
                warningTraced, thisThreadBinIt, a, BinData::getMemoryRequirements(*thisThreadBinIt) + totalBuffersNeeded);
        }
        catch (boost::iostreams::zlib_error &z)
        {
            warningTraced = handleBinAllocationFailure(warningTraced, thisThreadBinIt, z, z.error());
        }
        catch (common::IoException &io)
        {
            if (EMFILE == io.getErrorNumber())
            {
                warningTraced = handleBinAllocationFailure(warningTraced, thisThreadBinIt, io, "Increase number of open files for better efficiency");
            }
            else
            {
                throw;
            }
        }

        stateChangedCondition_.wait(lock);
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
    }

    ++allocatedBins_;
    ++nextUnallocatedBinIt;
    return ret;
}

void Build::waitForLoadSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber)
{
    bool warningTraced = false;

    while(nextUnloadedBinIt != thisThreadBinIt || !maxLoaders_)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        if (nextUnloadedBinIt == thisThreadBinIt && !warningTraced)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->get().getPath().c_str() << " until a load slot is available" << std::endl;
            warningTraced = true;
        }

        stateChangedCondition_.wait(lock);
    }

    ++nextUnloadedBinIt;
    --maxLoaders_;
}

void Build::returnLoadSlot(const bool exceptionUnwinding)
{
    ++maxLoaders_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

/**
 * @return true if this thread was the first to set task to 'complete" state
 */
bool Build::executePreemptTask(
    boost::unique_lock<boost::mutex>& lock,
    Task &task,
    const unsigned threadNumber) const
{
    //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << std::endl;
    ++task.threadsIn_;
    //    ISAAC_THREAD_CERR << "preempt " << &lock << " " << threadNumber << std::endl;
    try
    {
        task.execute(lock, threadNumber);
    }
    catch (...)
    {
        --task.threadsIn_;
        throw;
    }

    bool ret = false;

    if (!task.complete_)
    {
        //Threads don't come out of execute until there is nothing left to do.
        // stop new threads entering the task;
        task.complete_ = true;
        ret = true;
    }
    --task.threadsIn_;
    //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << " done" << std::endl;
    return ret || !task.threadsIn_;
}

bool Build::processMostUrgent(boost::unique_lock<boost::mutex> &lock, const unsigned threadNumber, Task *ownTask)
{
    // find the lowest priority incomplete and not busy task that is not higher than ownTask unless ownTask is 0
    Tasks::iterator highesttPriorityTask = std::min_element(
        tasks_.begin(), tasks_.end(), [ownTask](const Task *left, const Task *right)
        {
            if (ownTask && left->priority_ > ownTask->priority_)
            {
                return false;
            }

            if (left->complete_ || left->busy())
            {
                return false;
            }

            if (ownTask && right->priority_ > ownTask->priority_)
            {
                return true;
            }

            if (right->complete_ || right->busy())
            {
                return true;
            }
            return left->priority_ < right->priority_;
        });

    if(tasks_.end() == highesttPriorityTask)
    {
        return false;
    }

    Task *task = *highesttPriorityTask;
    if(task->complete_ || task->busy() || (ownTask && ownTask->busy() && task->priority_ > ownTask->priority_))
    {
        return false;
    }

    ISAAC_ASSERT_MSG(!ownTask || task->priority_ <= ownTask->priority_, "invalid task found");

    ISAAC_ASSERT_MSG(maxComputers_, "Unexpected maxComputers_ 0");
    --maxComputers_;

    bool ret = false;
    ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnComputeSlot, this, _1))
    {
        //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << std::endl;
        ret = executePreemptTask(lock, *task, threadNumber);
    }

    if (ret)
    {
        stateChangedCondition_.notify_all();
    }

    return ret;
}

bool Build::yieldIfPossible(
    boost::unique_lock<boost::mutex>& lock,
    const std::size_t threadNumber,
    Task *task)
{
    bool stateMightHaveChanged = false;
    if (maxComputers_)
    {
        stateMightHaveChanged = processMostUrgent(lock, threadNumber, task);
    }
    return stateMightHaveChanged;
}

template <typename OperationT>
void Build::preemptComputeSlot(
    boost::unique_lock<boost::mutex> &lock,
    const std::size_t maxThreads,
    const std::size_t priority,
    OperationT operation,
    const unsigned threadNumber)
{
    struct OperationTask : public Task
    {
        OperationT operation_;
        OperationTask(const std::size_t maxThreads, const std::size_t priority, OperationT operation):
            Task(maxThreads, priority), operation_(operation) {}
        virtual void execute(boost::unique_lock<boost::mutex> &l, const unsigned tn)
        {
//            ISAAC_THREAD_CERR << "Task::execute " << &l << " " << tn << std::endl;
            operation_(l, tn);
        }
    };
    OperationTask ourTask(maxThreads, priority, operation);
    ISAAC_ASSERT_MSG(tasks_.size() < tasks_.capacity(), "Unexpected high number of concurrent tasks. capacity: " << tasks_.capacity());
    tasks_.push_back(&ourTask);

    stateChangedCondition_.notify_all();

//    ISAAC_THREAD_CERR << "preemptComputeSlot " << &lock << " " << threadNumber << std::endl;
    // keep working until our task is complete.
    while (!ourTask.complete_)
    {
        if (forceTermination_)
        {
            // don't admit new threads
            ourTask.complete_ = true;
        }
        else if (!yieldIfPossible(lock, threadNumber, &ourTask))
        {
            stateChangedCondition_.wait(lock);
        }
    }

    // don't leave while there are some threads still in.
    while (ourTask.threadsIn_)
    {
        stateChangedCondition_.wait(lock);
    }

    ISAAC_ASSERT_MSG(tasks_.end() != std::find(tasks_.begin(), tasks_.end(), &ourTask), "Our task is gone from the list.");
    tasks_.erase(std::find(tasks_.begin(), tasks_.end(), &ourTask));

    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }
}

void Build::returnComputeSlot(const bool exceptionUnwinding)
{
    ++maxComputers_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::waitForSaveSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt)
{
    while(nextUnsavedBinIt != thisThreadBinIt)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        stateChangedCondition_.wait(lock);
    }
}

void Build::returnSaveSlot(
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
    const bool exceptionUnwinding)
{
    ++nextUnsavedBinIt;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::sortBinParallel(alignment::BinMetadataCRefList::const_iterator &nextUnprocessedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUncompressedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
                            common::ScopedMallocBlock &mallocBlock,
                            const std::size_t threadNumber)
{
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while(bins_.end() != nextUnprocessedBinIt)
    {
        alignment::BinMetadataCRefList::const_iterator thisThreadBinIt = nextUnprocessedBinIt++;

        // wait and allocate memory required for loading and compressing this bin
        boost::shared_ptr<BinData> binDataPtr =
            allocateBin(lock, thisThreadBinIt, bins_.end(), nextUnallocatedBinIt, mallocBlock, threadNumber);
        waitForLoadSlot(lock, thisThreadBinIt, bins_.end(), nextUnloadedBinIt, mallocBlock, threadNumber);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnLoadSlot, this, _1))
        {
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                BinLoader binLoader;
                binLoader.loadData(*binDataPtr);
            }
        }

        {
            preemptComputeSlot(
                lock, 1, std::distance(bins_.begin(), thisThreadBinIt),
                [this, &binDataPtr](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                {
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(l);
                    binSorter_.resolveDuplicates(*binDataPtr, stats_);
                },
                threadNumber);

            if (!binDataPtr->isUnalignedBin() && REALIGN_NONE != realignGaps_)
            {
                BinData::iterator nextUnprocessed = binDataPtr->indexBegin();
                int threadsIn = 0;
                preemptComputeSlot(
                    lock, -1, std::distance(bins_.begin(), thisThreadBinIt),
                    [this, &threadsIn, &binDataPtr, &nextUnprocessed](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                    {
                        ++threadsIn;
                        if (nextUnprocessed  == binDataPtr->indexBegin())
                        {
                            ISAAC_THREAD_CERR << "Realigning against " << getTotalGapsCount(binDataPtr->realignerGaps_) <<
                                " unique gaps. " << binDataPtr->bin_ << std::endl;
                        }
                        gapRealigner_.threadRealignGaps(l, *binDataPtr, nextUnprocessed, tn);
                        if (!--threadsIn)
                        {
                            ISAAC_THREAD_CERR << "Realigning gaps done. " << binDataPtr->bin_ << std::endl;
                        }
                    },
                    threadNumber);

            }

            preemptComputeSlot(
                lock, 1, std::distance(bins_.begin(), thisThreadBinIt),
                [this, &binDataPtr, threadNumber](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                {
                    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(l);
                    // Don't use tn!!! the streams have been allocated for the threadNumber.
                    binSorter_.serialize(
                        *binDataPtr, threadBgzfStreams_.at(threadNumber), threadBamIndexParts_.at(threadNumber));
                    threadBgzfStreams_.at(threadNumber).clear();
                },
                threadNumber);
        }
        // give back some memory to allow other threads to load
        // data while we're waiting for our turn to save
        binDataPtr.reset();

        // wait for our turn to store bam data
        waitForSaveSlot(lock, thisThreadBinIt, nextUnsavedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnSaveSlot, this, boost::ref(nextUnsavedBinIt), _1))
        {
            saveAndReleaseBuffers(lock, thisThreadBinIt->get().getPath(), threadNumber);
        }
    }

    // Don't release thread until all saving is done. Use threads that don't get anything to process for preemptive tasks such as realignment.
    while(!forceTermination_ && bins_.end() != nextUnsavedBinIt)
    {
        if (!yieldIfPossible(lock, threadNumber, 0))
        {
            stateChangedCondition_.wait(lock);
        }
    }
}

/**
 * \brief Save bgzf compressed buffers into corresponding sample files and and release associated memory
 */
void Build::saveAndReleaseBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const boost::filesystem::path &filePath,
    const std::size_t threadNumber)
{
    unsigned index = 0;
    BOOST_FOREACH(bam::BgzfBuffer &bgzfBuffer, threadBgzfBuffers_.at(threadNumber))
    {
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            std::ostream *stm = bamFileStreams_.at(index).get();
            if (!stm)
            {
                ISAAC_ASSERT_MSG(bgzfBuffer.empty(), "Unexpected data for bam file belonging to a sample with unmapped reference");
            }
            else
            {
                saveBuffer(bgzfBuffer, *stm, threadBamIndexParts_.at(threadNumber).at(index), bamIndexes_.at(index), filePath);
            }
        }
        // release rest of the memory that was reserved for this bin
        bam::BgzfBuffer().swap(bgzfBuffer);
        ++index;
    }
    --allocatedBins_;
    threadBamIndexParts_.at(threadNumber).clear();
}

void Build::saveBuffer(
    const bam::BgzfBuffer &bgzfBuffer,
    std::ostream &bamStream,
    const bam::BamIndexPart &bamIndexPart,
    bam::BamIndex &bamIndex,
    const boost::filesystem::path &filePath)
{
    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath.c_str() << std::endl;
    const clock_t start = clock();
    if(!bgzfBuffer.empty() && !bamStream.write(&bgzfBuffer.front(), bgzfBuffer.size())/* ||
        !bamStream.strict_sync()*/){
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to write bgzf block of %d bytes into bam stream") % bgzfBuffer.size()).str()));
    }
    bamIndex.processIndexPart( bamIndexPart, bgzfBuffer );

    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath.c_str() << " done in " << (clock() - start) / 1000 << "ms\n";
}

} // namespace build
} // namespace isaac
