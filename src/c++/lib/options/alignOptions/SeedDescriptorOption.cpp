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
 ** \file SeedDescriptorOption.cpp
 **
 ** seeds option parsing
 **
 ** \author Roman Petrovski
 **/
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/SeedMetadata.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"

#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

unsigned parseManualSeedDescriptor(
    const std::string &descriptor,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    unsigned ret = 0;
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;

    std::vector<std::string> offsetList;
    split(offsetList, descriptor,  is_any_of(":"));
    if (offsetList.empty())
    {
        const boost::format message = boost::format(
            "\n   *** The list of seed offsets for %s is empty. At least one seed is needed for each read ***\n") %
            readMetadata;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    BOOST_FOREACH(const std::string &offsetString, offsetList)
    {
        try
        {
            const unsigned offset = boost::lexical_cast<unsigned>(offsetString);
            alignment::SeedMetadata seedMetadata(offset, seedLength, readMetadata.getIndex(), seedMetadataList.size());
            if (offset + seedLength > readMetadata.getLength())
            {
                ISAAC_THREAD_CERR << "WARNING: ignored " << seedMetadata <<
                    " as it stretches beyond the read " << readMetadata.getNumber() <<
                    " which is " << readMetadata.getLength() << " bases long" << std::endl;
            }
            else
            {
                seedMetadataList.push_back(seedMetadata);
                ++ret;
                ISAAC_THREAD_CERR << "Constructed manual " << seedMetadataList.back() << std::endl;
            }
        }
        catch(boost::bad_lexical_cast &)
        {
            const boost::format message = boost::format("\n   *** Invalid seed offset '%s' found in '%s' ***\n") %
                offsetString % descriptor;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }
    return ret;
}


unsigned makeExtremitySeeds(
    unsigned offset,
    unsigned endOffset,
    const unsigned seedLength,
    const unsigned readIndex,
    alignment::SeedMetadataList &seedMetadataList)
{
    // fullSeedLength is a pair of seeds next to each other.
    const unsigned fullSeedLength = seedLength * 2;
    unsigned length = endOffset - offset;
    if (length >= fullSeedLength)
    {
        alignment::SeedMetadata seedMetadata(offset, seedLength, readIndex, seedMetadataList.size());
        offset += fullSeedLength;
        length -= fullSeedLength;
        seedMetadataList.push_back(seedMetadata);
    }
    else
    {
        return 0;
    }

    if (length >= fullSeedLength)
    {
        endOffset -= fullSeedLength;
        length -= fullSeedLength;
        alignment::SeedMetadata seedMetadata(endOffset, seedLength, readIndex, seedMetadataList.size());
        seedMetadataList.push_back(seedMetadata);
    }
    else
    {
        return 1;
    }

    return 2 + makeExtremitySeeds(offset, endOffset, seedLength, readIndex, seedMetadataList);
}


/**
 * \return number of seeds generated
 */
unsigned parseAutoSeedDescriptor(
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    int index = seedMetadataList.size();
    const unsigned ret = makeExtremitySeeds(0, readMetadata.getLength(), seedLength, readMetadata.getIndex(), seedMetadataList);

    // alignment expects the offsets to grow (although users can override offsets in any order they want).
    std::sort(seedMetadataList.begin() + index, seedMetadataList.end(),
              boost::bind(&alignment::SeedMetadata::getOffset, _1) < boost::bind(&alignment::SeedMetadata::getOffset, _2));
    // renumber according to sort order
    BOOST_FOREACH(alignment::SeedMetadata &seedMetadata, std::make_pair(seedMetadataList.begin() + index, seedMetadataList.end()))
    {
        seedMetadata = alignment::SeedMetadata(seedMetadata.getOffset(), seedMetadata.getLength(), seedMetadata.getReadIndex(), index);
        ISAAC_THREAD_CERR << "Constructed auto " << seedMetadata << std::endl;
        ++index;
    }

    return ret;
}


/**
 * \return Number of generated seeds
 */
unsigned parseStepSeedDescriptor(
    const unsigned step,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    ISAAC_ASSERT_MSG(readMetadata.getLength() >= seedLength, "Read is too short for seed length " << seedLength << " " << readMetadata);

    unsigned i = 0;
    for (; i < readMetadata.getLength() - seedLength; i += step)
    {
        const alignment::SeedMetadata seedMetadata(i, seedLength, readMetadata.getIndex(), seedMetadataList.size());
        ISAAC_THREAD_CERR << "Constructed step " << seedMetadata << std::endl;
        seedMetadataList.push_back(seedMetadata);
    }
    return i;
}

/**
 * \return Maximum number of first pass seeds possible
 */
unsigned parseReadSeedDescriptor(
    const std::string &descriptor,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    if (descriptor.empty())
    {
        const boost::format message = boost::format(
            "\n   *** The seed descriptor for %s is empty. At least one seed is needed ***\n") %
            readMetadata;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    if ("step=" == descriptor.substr(0, 5))
    {
        return parseStepSeedDescriptor(boost::lexical_cast<unsigned>(descriptor.substr(5)), readMetadata, seedLength, seedMetadataList);
    }
    else if ("all" == descriptor)
    {
        return parseStepSeedDescriptor(1, readMetadata, seedLength, seedMetadataList);
    }
    else if ("auto" == descriptor)
    {
        return parseAutoSeedDescriptor(readMetadata, seedLength, seedMetadataList);
    }
    else
    {
        return parseManualSeedDescriptor(descriptor, readMetadata, seedLength, seedMetadataList);
    }
}

/**
 * Parses seed descriptor into seed objects.
 */
alignment::SeedMetadataList parseSeedDescriptor(
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::string &seedDescriptor,
    const unsigned seedLength)
{
    if (seedDescriptor.empty())
    {
        const boost::format message = boost::format("\n   *** The seed descriptor is empty. At least one seed is needed ***\n");
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    alignment::SeedMetadataList seedMetadataList;
    std::vector<std::string> seedDescriptorList; // split by read
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;
    split(seedDescriptorList, seedDescriptor,  is_any_of(","));
    if (readMetadataList.size() < seedDescriptorList.size())
    {
        const boost::format message = boost::format("\n   *** Too many lists-of-seeds in seed-descriptor '%s': found %d: %d reads only ***\n") %
            seedDescriptor % seedDescriptorList.size() % readMetadataList.size();
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    // extend the last list-of-seeds to all subsequent reads if needed
    seedDescriptorList.resize(readMetadataList.size(), seedDescriptorList.back());
    // create all the seeds
    std::vector<flowcell::ReadMetadata>::const_iterator readMetadata = readMetadataList.begin();
    BOOST_FOREACH(const std::string &descriptor, seedDescriptorList)
    {
        parseReadSeedDescriptor(descriptor, *readMetadata, seedLength, seedMetadataList);
        ++readMetadata;
    }

    return seedMetadataList;
}

} // namespace alignOptions
} // namespace option
} // namespace isaac
