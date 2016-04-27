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
 ** \file IndelLoader.cpp
 **
 ** Reads indels from a vcf file.
 ** 
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <unordered_map>

#include "build/IndelLoader.hh"
#include "common/Debug.hh"
#include "vcf/VcfUtils.hh"


namespace isaac
{
namespace build
{

build::gapRealigner::Gaps loadIndels(
    const boost::filesystem::path& vcfFilePath,
    const reference::SortedReferenceMetadata& sortedReferenceMetadata)
{
    build::gapRealigner::Gaps ret;
    typedef std::unordered_map<std::string, unsigned> ContigLookup;
    ContigLookup contigLookup;

    const reference::SortedReferenceMetadata::Contigs &contigs = sortedReferenceMetadata.getContigs();
    std::for_each(
        contigs.begin(), contigs.end(),
        [&](const reference::SortedReferenceMetadata::Contig& contig)
        {
            contigLookup.insert(ContigLookup::value_type(contig.name_, contig.index_));
        }
    );

    std::ifstream ifs(vcfFilePath.c_str());
    if (!ifs)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno,
            (boost::format("ERROR: Unable to open known indels file: %s") % vcfFilePath.c_str()).str()));
    }

    std::string line;
    std::size_t lineNumber = 0;
    typedef vcf::ParsedVcfLine<std::string::const_iterator> ParsedLine;
    ParsedLine parsedLine;
    std::size_t insertions = 0;
    std::size_t deletions = 0;
    while (std::getline(ifs, line) || ifs.eof())
    {
        ++lineNumber;
        vcf::parseVcfLine<std::string::const_iterator>(line.begin(), line.end(),
                                                       parsedLine);
//        ISAAC_THREAD_CERR << parsedLine << std::endl;
        if (!parsedLine.isValid())
        {
            BOOST_THROW_EXCEPTION(vcf::VcfError(
                (boost::format("ERROR: %s:%d. Incorrect VCF syntax: %s") %
                    vcfFilePath.c_str() % lineNumber % line).str()));
        }
        if (parsedLine.isDataLine())
        {
            ContigLookup::const_iterator contigIndex = contigLookup.find(parsedLine.getChromosome());
            if (contigLookup.end() != contigIndex)
            {
                for (ParsedLine::AlternativesConstIterator it = parsedLine.alternativesBegin();
                    parsedLine.alternativesEnd() != it; ++it)
                {
                    const ParsedLine &alternative = *it;
                    if (1 != alternative.refLength() && 1 != alternative.altLength())
                    {
                        ISAAC_THREAD_CERR << "WARNING: " << vcfFilePath.c_str() << ":" << lineNumber <<
                            " MNV ignored: " << alternative << " " << line << std::endl;
                    }
                    else if (alternative.refLength() != alternative.altLength())
                    {
                        ret.push_back(
                            build::gapRealigner::Gap(
                                reference::ReferencePosition(contigIndex->second, parsedLine.getPosition()),
                                alternative.refLength() - alternative.altLength()));
                        insertions += ret.back().isInsertion();
                        deletions += ret.back().isDeletion();
                    }
                }
            }
            else
            {
                ISAAC_THREAD_CERR << "WARNING: " << vcfFilePath.c_str() << ":" << lineNumber <<
                    " Unknown CROM record ignored: " << line << std::endl;
            }
        }
//            ISAAC_THREAD_CERR << ret.back() << std::endl;
        if (ifs.eof())
        {
            break;
        }
    }

    ISAAC_THREAD_CERR << "Read " << insertions << " known insertions and " << deletions << " deletions from " << vcfFilePath.c_str() << std::endl;

    return ret;
}

build::gapRealigner::Gaps loadIndels(
    const boost::filesystem::path &vcfFilePath,
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList)
{
    ISAAC_ASSERT_MSG(1 == sortedReferenceMetadataList.size(), "Multiple references are not supported");
    const reference::SortedReferenceMetadata &sortedReferenceMetadata  = sortedReferenceMetadataList.front();
    return loadIndels(vcfFilePath, sortedReferenceMetadata);
}

} // namespace build
} // namespace isaac
