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
 ** \file MatchDistribution.hh
 **
 ** \brief Tracking of the distribution of the matches across the reference
 ** genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH
#define iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH

#include <vector>

#include "common/Debug.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace alignment
{

struct MatchDistribution : public std::vector<std::vector<unsigned> >
{
    MatchDistribution(uint64_t binSize = 0) : binSize_(binSize)
    {

    }

    uint64_t getBinSize() const {return binSize_;}

protected:
    uint64_t binSize_;
};

/**
 ** \brief Pretends it knows the match distribution by using the target coverage and given read length.
 **/
class EstimatedMatchDistribution: public MatchDistribution
{
public:
    /**
     ** \brief Initialize the MatchDistribution for the specified references.
     **
     ** Creates one vector for each contig, using the max contig length of each contig
     ** from all the provided references to infer the number of bins,
     ** and initialize all bin counts to a number that would achieve desired coverage.
     **
     ** \param readLength  maximum read length out of all the input data reads.
     **/
    EstimatedMatchDistribution(
        const unsigned expectedCoverage,
        const unsigned readLength,
        const isaac::reference::SortedReferenceMetadataList &sortedReferenceMetadataList):
            MatchDistribution(std::max(1UL, isaac::reference::getLongestGenomeLength(sortedReferenceMetadataList) / BINS_PER_GENOME))
    {
        ISAAC_THREAD_CERR << "EstimatedMatchDistribution binSize_=" << binSize_ << std::endl;
        BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata &sortedReferenceMetadata, sortedReferenceMetadataList)
        {
            const std::vector<isaac::reference::SortedReferenceMetadata::Contig> xmlContigs = sortedReferenceMetadata.getContigs();
            resize(std::max(size(), xmlContigs.size()));
            BOOST_FOREACH(const isaac::reference::SortedReferenceMetadata::Contig &xmlContig, xmlContigs)
            {
                const uint64_t binCount = (xmlContig.totalBases_ + getBinSize() - 1) / getBinSize();
                at(xmlContig.karyotypeIndex_).resize(std::max(binCount, at(xmlContig.karyotypeIndex_).size()), 0);
            }
        }
        // normally this is just a rounded up genome length. However in multiplexed cases this is a number
        // representing the maximum possible length of the genome that would contain the longest of corresponding contigs
        uint64_t estimatedGenomeLength = 0;
        std::size_t bins = 0;
        BOOST_FOREACH(const std::vector<unsigned> &contigBins, *this)
        {
            estimatedGenomeLength += binSize_ * contigBins.size();
            bins += contigBins.size();
        }

        const uint64_t fragmentsPerBin = estimatedGenomeLength * expectedCoverage / readLength  / bins;

        BOOST_FOREACH(std::vector<unsigned> &contigBins, *this)
        {
            std::fill(contigBins.begin(), contigBins.end(), fragmentsPerBin);
        }
    }
    bool isEmptyContig(const size_t contigIndex) const
    {
        return false;
    }
private:
    /**
     * Bins must be granular enough to allow MatchSelector binning flexibility for highly-covered tiny genomes (e.g. PhiX).
     * This presents a bit of a memory challenge when finding matches for large genomes on big number of threads.
     * For example:
     *  64 cores, 2^32 bases long human genome, 2^8 Bin size:
     *      64 * 2^32 / 2^8 * sizeof(unsigned) = 4 gigabytes. And we need those gigabytes for seeds, reference, etc...
     * maintain large number of bins, use largest genome given to decide on bin size
     */
    static const uint64_t BINS_PER_GENOME = 1024 * 100;
    /// convert a position on a contig into a bin index
    size_t getBinIndex(uint64_t position) const {return position / getBinSize();}
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_DISTRIBUTION_HH
