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
 ** \file AnnotationLoader.cpp
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/
#include <errno.h>

#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "reference/AnnotationLoader.hh"

namespace isaac
{
namespace reference
{

template <typename MemFn>
void loadContigAnnotations(
    const SortedReferenceMetadata::Contigs& contigs,
    std::istream &is, const boost::filesystem::path& path,
    ContigAnnotations &ret,
    MemFn memFn)
{
    BOOST_FOREACH(const SortedReferenceMetadata::Contig &contig, contigs)
    {
        ContigAnnotation &contigAnnotation = ret.at(contig.karyotypeIndex_);
        contigAnnotation.resize(contig.totalBases_);
        std::size_t currentOffset = 0;
        for (std::size_t pos = 0; pos < contigAnnotation.size(); ++pos)
        {
            if (!is.read(reinterpret_cast<char *>(&memFn(contigAnnotation[pos])), sizeof(AnnotationValue)) ||
                sizeof(AnnotationValue) != std::size_t(is.gcount()))
            {
                const boost::format message = boost::format("Failed to read %d bytes at offset %d from annotation file %s. Read %d") %
                    (sizeof(AnnotationValue)) % (currentOffset + pos) % path.string() % is.gcount();
                BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
            }
        }
        currentOffset += contigAnnotation.size();
    }
}

template <typename MemFn>
void loadAnnotations(
    const boost::filesystem::path& path,
    const isaac::reference::SortedReferenceMetadata::Contigs& contigs,
    ContigAnnotations &ret,
    MemFn memFn)
{
    boost::iostreams::filtering_istream is;
    if (common::isDotGzPath(path))
    {
        is.push(boost::iostreams::gzip_decompressor());
    }
    is.push(boost::iostreams::file_source(path.string()));
    if (!is)
    {
        BOOST_THROW_EXCEPTION(
            common::IoException(errno, "Failed to open annotation file " + path.string()));
    }
    loadContigAnnotations(contigs, is, path, ret, memFn);
}

/**
 * \brief split single reference annotation blob into per-contig annotations
 */
ContigAnnotations loadAnnotations(const reference::SortedReferenceMetadata &sortedReferenceMetadata)
{
    isaac::reference::SortedReferenceMetadata::Contigs contigs = sortedReferenceMetadata.getKaryotypeOrderedContigs();

    ContigAnnotations ret(contigs.size());
    loadAnnotations(sortedReferenceMetadata.getKUniquenessAnnotation().path_, contigs, ret, boost::mem_fn(&ContigAnnotation::value_type::first));
    loadAnnotations(sortedReferenceMetadata.getKRepeatnessAnnotation().path_, contigs, ret, boost::mem_fn(&ContigAnnotation::value_type::second));
    return ret;
}

/**
 * \brief load contig annotations for each reference in the list
 */
ContigAnnotationsList loadAnnotations(const reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
{
    ISAAC_TRACE_STAT("loadAnnotations ");

    ContigAnnotationsList ret(sortedReferenceMetadataList.size());
    std::size_t referenceIndex = 0;
    BOOST_FOREACH(const reference::SortedReferenceMetadata &ref, sortedReferenceMetadataList)
    {
        if (ref.hasKUniquenessAnnotation())
        {
            ret.at(referenceIndex) = reference::loadAnnotations(ref);
        }
        else
        {
            ISAAC_THREAD_CERR << "WARNING: No annotation is available for reference genome. Alignment scores could be misleading. " <<
                ref.getContigs().front().filePath_ << std::endl;
        }
        ++referenceIndex;
    }

    ISAAC_TRACE_STAT("loadAnnotations done ");

    return ret;
}

} // namespace reference
} // namespace isaac
