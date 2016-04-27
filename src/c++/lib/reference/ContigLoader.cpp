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
 ** \file ContigLoader.cpp
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "reference/ContigLoader.hh"

namespace isaac
{
namespace reference
{

void loadContig(
    const reference::SortedReferenceMetadata::Contig &xmlContig,
    Contig &contig)
{
    contig.clear();
    contig.reserve(xmlContig.totalBases_);
    std::ifstream is(xmlContig.filePath_.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open reference file " + xmlContig.filePath_.string()));
    }
    if (!is.seekg(xmlContig.offset_))
    {
        using common::IoException;
        using boost::format;
        const format message = (boost::format("Failed to reach offset %d in reference file % s") % xmlContig.offset_ % xmlContig.filePath_);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
//        ISAAC_THREAD_CERR << (boost::format("Contig seek %s (%3d:%8d): %s") % xmlContig.name_ % xmlContig.index_ % xmlContig.totalBases_ % xmlContig.filePath_).str() << std::endl;
    static const oligo::Translator<true> translator;
    char base = 0;
    std::size_t acgtBases = 0;
    while(is && (contig.size() < xmlContig.totalBases_) && is.get(base))
    {
        if ('\r' != base && '\n' != base)
        {
            ISAAC_ASSERT_MSG(std::isalpha(base), "Invalid base read from " << xmlContig << " : " << base);
            {
                contig.push_back(oligo::getBase(translator[base], true));
                acgtBases += oligo::REFERENCE_OLIGO_N != contig.back();
            }
        }
    }
    if (xmlContig.totalBases_ != contig.size())
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d bases from reference file % s: %d") % xmlContig.totalBases_ % xmlContig.filePath_.string() % contig.size());
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    if (xmlContig.acgtBases_ != acgtBases)
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d ACGT bases from reference file % s: %d") %
            xmlContig.acgtBases_ % xmlContig.filePath_.string() % acgtBases);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    ISAAC_TRACE_STAT("Loaded contig " << xmlContig << " ");
}

struct DummyFilter {bool operator() (const unsigned contigIdx) const {return true;}} dummyFilter;
/**
 * \brief loads all the fasta file contigs into memory on multiple threads
 */
reference::ContigList loadContigs(
    const reference::SortedReferenceMetadata::Contigs &xmlContigs,
    common::ThreadVector &loadThreads)
{
    return loadContigs(xmlContigs, dummyFilter, loadThreads);
}


} // namespace reference
} // namespace isaac
