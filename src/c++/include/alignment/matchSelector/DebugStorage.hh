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

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH

#include "alignment/matchSelector/FragmentStorage.hh"
#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "alignment/matchSelector/debugStorage/QqStatistics.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class DebugStorage: public FragmentStorage
{
    static const unsigned READS_MAX = 2;
    const reference::ContigList &contigList_;
    const reference::ContigAnnotations &contigAnnotations_;
    const flowcell::Layout &flowcell_;
    std::array<Cigar, READS_MAX> originalCigars_;
    const AlignmentCfg &alignmentCfg_;
    const boost::filesystem::path outputDirectory_;
    debugStorage::MapqStatistics r1MapqStatistics_;
    debugStorage::MapqStatistics r2MapqStatistics_;
    debugStorage::MapqStatistics pairMapqStatistics_;
    debugStorage::QqStatistics r1QqOriginalStatistics_;
    debugStorage::QqStatistics r2QqOriginalStatistics_;
    debugStorage::QqStatistics r1QqStatistics_;
    debugStorage::QqStatistics r2QqStatistics_;
    alignment::matchSelector::FragmentStorage &actualStorage_;
public:
    DebugStorage(
        const reference::ContigLists &contigLists,
        const reference::ContigAnnotationsList &contigAnnotationsList,
        const AlignmentCfg &alignmentCfg,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const boost::filesystem::path &outputDirectory,
        alignment::matchSelector::FragmentStorage &actualStorage);

    ~DebugStorage();

private:
    virtual void store(
        const BamTemplate &bamTemplate,
        const unsigned barcodeIdx);

    virtual void reset(const uint64_t clusterId, const bool paired)
    {
        actualStorage_.reset(clusterId, paired);
    }

    virtual void prepareFlush() noexcept
    {
        actualStorage_.prepareFlush();
    }
    virtual void flush()
    {
        actualStorage_.flush();
    }
    virtual void resize(const uint64_t clusters)
    {
        actualStorage_.resize(clusters);
    }
    virtual void reserve(const uint64_t clusters)
    {
        actualStorage_.reserve(clusters);
    }

private:
    bool updateMapqStats(const BamTemplate& bamTemplate);
    bool restoreOriginal(const std::size_t readNumber, FragmentMetadata &fragment) const;
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH
