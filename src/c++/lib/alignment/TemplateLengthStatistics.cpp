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
 ** \file TemplateLengthStatistics.cpp
 **
 ** \brief See TemplateLengthStatistics.hh
 ** 
 ** \author Come Raczy
 **/


#include "alignment/Quality.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/Cigar.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

const double TemplateLengthDistribution::STANDARD_DEVIATIONS_MAX = 3.0;
const double TemplateLengthDistribution::FRAGMENT_LENGTH_CONFIDENCE_INTERVAL =
    boost::math::erf(STANDARD_DEVIATIONS_MAX/boost::math::constants::root_two<double>());
const double TemplateLengthDistribution::FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z = boost::math::erf(1.0/boost::math::constants::root_two<double>());
const double TemplateLengthDistribution::LOWER_PERCENT = (1.0 - FRAGMENT_LENGTH_CONFIDENCE_INTERVAL) / 2.0;
const double TemplateLengthDistribution::UPPER_PERCENT = (1.0 + FRAGMENT_LENGTH_CONFIDENCE_INTERVAL) / 2.0;
const double TemplateLengthDistribution::LOWER_PERCENT_1Z = (1.0 - FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z) / 2.0;
const double TemplateLengthDistribution::UPPER_PERCENT_1Z = (1.0 + FRAGMENT_LENGTH_CONFIDENCE_INTERVAL_1Z) / 2.0;

TemplateLengthStatistics::TemplateLengthStatistics()
    : min_(-1U)
    , max_(-1U)
    , median_(-1U)
    , lowStdDev_(-1U)
    , highStdDev_(-1U)
    , stable_(false)
    , mateMin_(-1U)
    , mateMax_(-1U)
{
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;

}

void TemplateLengthStatistics::clear()
{
    min_ = -1U;
    max_ = -1U;
    median_ = -1U;
    lowStdDev_ = -1U;
    highStdDev_ = -1U;
    stable_ = false;
    bestModels_[0] = InvalidAlignmentModel;
    bestModels_[1] = InvalidAlignmentModel;
}

bool TemplateLengthStatistics::matchModel(const FragmentMetadata &f1, const FragmentMetadata &f2) const
{
    if (f1.getContigId() == f2.getContigId())
    {
        const uint64_t length = getLength(f1, f2);
        const AlignmentModel model = alignmentModel(f1, f2);
        const bool ret = (length <= max_ + TEMPLATE_LENGTH_THRESHOLD) && ((model == bestModels_[0]) || (model == bestModels_[1]));
        ISAAC_THREAD_CERR_DEV_TRACE((boost::format("TemplateLengthStatistics::matchModel %d: %s-%s (len: %u)") % ret %
                                     f1 % f2 % length).str());
        return ret;
    }
    return false;

}

TemplateLengthStatistics::AlignmentClass TemplateLengthStatistics::alignmentClass(const AlignmentModel model)
{
    assert(model < 8);
    return static_cast<AlignmentClass>((model < 4) ? model : ((~model) & 3));
}

const static std::vector<std::string> classNames = boost::assign::list_of("F+")("R+")("R-")("F-")("unknown");
const std::string &TemplateLengthStatistics::alignmentClassName(const AlignmentClass alignmentClass)
{
    return (classNames.size() > size_t(alignmentClass)) ? classNames.at(alignmentClass) : classNames.back();
}

const static std::vector<std::string> modelNames = boost::assign::list_of
        ("FF+")("FR+")("RF+")("RR+")("FF-")("FR-")("RF-")("RR-")("unknown");
const std::string &TemplateLengthStatistics::alignmentModelName(const AlignmentModel alignmentModel)
{
    return (modelNames.size() > size_t(alignmentModel)) ? modelNames.at(alignmentModel) : modelNames.back();
}

void TemplateLengthDistribution::updateStatistics()
{
    // save the old statistics
    const TemplateLengthStatistics oldStats = stats_;
    // find the two best alignment models
    stats_.setBestModel(histograms_[1].samples() <= histograms_[0].samples() ? TemplateLengthStatistics::FFp : TemplateLengthStatistics::FRp, 0);
    stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>((stats_.getBestModel(0) + 1) % 2), 1);
    for (size_t i = 2; histograms_.size() > i; ++i)
    {
        if (histograms_[i].samples() > histograms_[stats_.getBestModel(0)].samples())
        {
            stats_.setBestModel(stats_.getBestModel(0), 1);
            stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>(i), 0);
        }
        else if (histograms_[i].samples() > histograms_[stats_.getBestModel(1)].samples())
        {
            stats_.setBestModel(static_cast<TemplateLengthStatistics::AlignmentModel>(i), 1);
        }
    }

    // TODO: report incoherent models
    if (TemplateLengthStatistics::alignmentClass(stats_.getBestModel(0)) != TemplateLengthStatistics::alignmentClass(stats_.getBestModel(1)))
    {
        ISAAC_THREAD_CERR << "Incoherent alignment models:"
            << stats_.getBestModel(0) << " (" << TemplateLengthStatistics::alignmentModelName(stats_.getBestModel(0)) << "->"
                << TemplateLengthStatistics::alignmentClassName(TemplateLengthStatistics::alignmentClass(stats_.getBestModel(0))) << "):"
            << stats_.getBestModel(1) << " (" << TemplateLengthStatistics::alignmentModelName(stats_.getBestModel(1)) << "->"
                << TemplateLengthStatistics::alignmentClassName(TemplateLengthStatistics::alignmentClass(stats_.getBestModel(1))) << ")" << std::endl;
    }

    Histogram twoBestSum = histograms_[stats_.getBestModel(0)];
    twoBestSum += histograms_[stats_.getBestModel(1)];

    unsigned minLength = 0;
    unsigned lowStdDev = 0;
    unsigned median = TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD / 2;
    unsigned highStdDev = 0;
    unsigned maxLength = TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD;

    std::size_t pos = 0;
    for (std::size_t i = 0; i < TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD; ++i)
    {
        if (pos <= std::size_t(twoBestSum.samples() * LOWER_PERCENT))
        {
            minLength = i;
        }
        if (pos <= std::size_t(twoBestSum.samples() * LOWER_PERCENT_1Z))
        {
            lowStdDev = i;
        }
        if (pos <= std::size_t(twoBestSum.samples() * 0.5))
        {
            median = i;
        }
        if (pos <= std::size_t(twoBestSum.samples() * UPPER_PERCENT_1Z))
        {
            highStdDev = i;
        }
        if (pos <= std::size_t(twoBestSum.samples() * UPPER_PERCENT))
        {
            maxLength = i;
        }
        pos += twoBestSum[i];
    }

    lowStdDev = median - lowStdDev;
    highStdDev = highStdDev - median;

    setMin(minLength);
    setMedian(median);
    setMax(maxLength);
    stats_.setLowStdDev(lowStdDev);
    stats_.setHighStdDev(highStdDev);

//    ISAAC_THREAD_CERR << "updateStatistics oldStats:" << oldStats <<  " stats_:" << stats_ << std::endl;
//    ISAAC_THREAD_CERR << "updateStatistics templateCount_:" << templateCount_ <<  " uniqueCount_:" << uniqueCount_ << std::endl;
    // check if we reached stability
    if (oldStats.getMin() == stats_.getMin() &&
        oldStats.getMedian() == stats_.getMedian() &&
        oldStats.getMax() == stats_.getMax() &&
        oldStats.getLowStdDev() == stats_.getLowStdDev() &&
        oldStats.getHighStdDev() == stats_.getHighStdDev() &&
        oldStats.getBestModel(0) == stats_.getBestModel(0) &&
        oldStats.getBestModel(1) == stats_.getBestModel(1))
    {
        stats_.setStable(true);
    }
}

/// Check that one of the two best models has the given orientation for the readIndex
bool TemplateLengthStatistics::isValidModel(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(READS_MAX > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    return (reverse == ((bestModels_[0] >> shift) & 1)) || (reverse == ((bestModels_[1] >> shift) & 1));
}

bool TemplateLengthStatistics::firstFragment(const bool reverse, const unsigned readIndex) const
{
    ISAAC_ASSERT_MSG(READS_MAX > readIndex, "Invalid read index");
    // bit 0 is for read 2; bit 1 is for read 1
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return ((bestModels_[i] >> 2) & 1)== readIndex;
        }
    }
    assert(false); // shouldn't get there
    return false;
}

bool TemplateLengthStatistics::mateOrientation(unsigned readIndex, bool reverse) const
{
    // bit 0 is for readIndex 1 and vice-versa
    const unsigned shift = (readIndex + 1) % 2;
    for (unsigned i = 0; 2 > i; ++i)
    {
        if (reverse == ((bestModels_[i] >> shift) & 1))
        {
            return (bestModels_[i] >> readIndex) & 1;
        }
    }
    // for discorant models, return the expected orientation in the best model
    return (bestModels_[0] >> readIndex) & 1;
}

int64_t TemplateLengthStatistics::mateMinPosition(
    const unsigned readIndex,
    const bool reverse,
    const int64_t position,
    const unsigned *readLengths) const
{
    ISAAC_ASSERT_MSG(isValidModel(reverse, readIndex), "Invalid method call");
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMin_/* / 2*/ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMax_/* * 2*/ + readLengths[readIndex];
    }
}

int64_t TemplateLengthStatistics::mateMaxPosition(
    const unsigned readIndex,
    const bool reverse,
    const int64_t position,
    const unsigned *readLengths) const
{
    ISAAC_ASSERT_MSG(isValidModel(reverse, readIndex), "Invalid method call");
    if (firstFragment(reverse, readIndex))
    {
        return position + mateMax_/* * 2*/ - readLengths[(readIndex + 1) % 2];
    }
    else
    {
        return position - mateMin_/* / 2*/ + readLengths[readIndex];
    }
}



TemplateLengthDistribution::TemplateLengthDistribution(int mateDriftRange)
    : mateDriftRange_(mateDriftRange)
    , count_(0)
{
}

void TemplateLengthDistribution::clear()
{
    stats_.clear();
    count_ = 0;
    std::for_each(histograms_.begin(), histograms_.end(), std::mem_fun_ref(&Histogram::clear));
}


inline FragmentMetadataList::const_iterator findUniqueOrUniquePerfectAlignment(
    const FragmentMetadataList &list)
{
    // if the alignment is unique, assume it is good enough.
    if (1 == list.size())
    {
        return list.begin();
    }

    // if there are too many candidates, just skip it
    if (list.size() > 10)
    {
        return list.end();
    }

    // if there are multiple candidates, only look at the perfect ones and only if the perfect ones are unique.
    
    // Notice that there is still a problem

    const FragmentMetadataList::const_iterator ret = std::find_if(list.begin(), list.end(), !boost::bind(&FragmentMetadata::getEditDistance, _1));
    if (list.end() == ret)
    {
        return list.end();
    }

    const FragmentMetadataList::const_iterator it = std::find_if(ret + 1, list.end(), !boost::bind(&FragmentMetadata::getEditDistance, _1));
    if (list.end() == it)
    {
        return ret;
    }

    return list.end();
}

bool TemplateLengthDistribution::appendTemplate(const FragmentMetadataList &r0Fragments, const FragmentMetadataList &r1Fragments)
{
    const AlignmentModel model = getAlignmentModel(r0Fragments, r1Fragments);
    if (TemplateLengthStatistics::InvalidAlignmentModel != model.first)
    {
        appendModel(model.first, model.second);
    }
    return stats_.isStable();
}

void TemplateLengthDistribution::appendModel(
    const TemplateLengthStatistics::AlignmentModel am,
    const uint64_t length)
{
    if (histograms_[am].append(length))
    {
        ++count_;
    }
    // calculate the alignment statistics if appropriate
    if (0 == (count_ % UPDATE_FREQUENCY))
    {
        // save the old statistics
        const TemplateLengthStatistics oldStats = stats_;
        // get the new statistics
        updateStatistics();
        // check if we reached stability
        if (oldStats.getMin() == stats_.getMin()
            && oldStats.getMedian() == stats_.getMedian()
            && oldStats.getMax() == stats_.getMax()
            && oldStats.getLowStdDev() == stats_.getLowStdDev()
            && oldStats.getHighStdDev() == stats_.getHighStdDev())
        {
            stats_.setStable(true);
        }
    }
}

bool TemplateLengthDistribution::appendUniqueTemplate(const FragmentMetadata &perfect0, const FragmentMetadata &perfect1)
{
    const AlignmentModel model = getAlignmentModel(perfect0, perfect1);
    if (TemplateLengthStatistics::InvalidAlignmentModel != model.first)
    {
        appendModel(model.first, model.second);
    }
    return stats_.isStable();
}

TemplateLengthDistribution::AlignmentModel
TemplateLengthDistribution::getAlignmentModel(const FragmentMetadataList &r0Fragments, const FragmentMetadataList &r1Fragments)
{
    // discard templates where at least on fragment didn't align
    if (r0Fragments.empty() || r1Fragments.empty())
    {
        return AlignmentModel(TemplateLengthStatistics::InvalidAlignmentModel, 0U);
    }
    // discard templates where the alignment is not unique on both fragments
    const FragmentMetadataList::const_iterator perfect0 = findUniqueOrUniquePerfectAlignment(r0Fragments);
    if (r0Fragments.end() == perfect0)
    {
        return AlignmentModel(TemplateLengthStatistics::InvalidAlignmentModel, 0U);
    }

    const FragmentMetadataList::const_iterator perfect1 = findUniqueOrUniquePerfectAlignment(r1Fragments);
    if (r1Fragments.end() == perfect1)
    {
        return AlignmentModel(TemplateLengthStatistics::InvalidAlignmentModel, 0U);
    }
    return getAlignmentModel(*perfect0, *perfect1);
}

TemplateLengthDistribution::AlignmentModel
TemplateLengthDistribution::getAlignmentModel(const FragmentMetadata &perfect0, const FragmentMetadata &perfect1)
{
    // discard templates that span across several contigs
    if (perfect0.contigId != perfect1.contigId)
    {
        return AlignmentModel(TemplateLengthStatistics::InvalidAlignmentModel, 0U);
    }

    // calculate the lenght of the template
    const uint64_t length = TemplateLengthStatistics::getLength(perfect0, perfect1);
    // discard excessively long templates
    if (length > TemplateLengthStatistics::TEMPLATE_LENGTH_THRESHOLD)
    {
        return AlignmentModel(TemplateLengthStatistics::InvalidAlignmentModel, 0U);
    }
    // update the histogram for the appropriate alignment model
    const TemplateLengthStatistics::AlignmentModel am = TemplateLengthStatistics::alignmentModel(perfect0, perfect1);
    return std::make_pair(am, length);
}

bool TemplateLengthDistribution::finalize()
{
    // save the old statistics
    TemplateLengthStatistics oldStats = stats_;
    // get the new statistics
    updateStatistics();
    // check if we reached stability
    if (oldStats.getMin() == stats_.getMin() &&
        oldStats.getMedian() == stats_.getMedian() &&
        oldStats.getMax() == stats_.getMax() &&
        oldStats.getLowStdDev() == stats_.getLowStdDev() &&
        oldStats.getHighStdDev() == stats_.getHighStdDev())
    {
        stats_.setStable(true);
    }
    return stats_.isStable();
}


} // namespace alignment
} // namespace isaac
