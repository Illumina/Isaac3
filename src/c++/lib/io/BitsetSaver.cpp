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
 ** \file BitsetSaver.cpp
 **
 ** Helper for saving neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "io/BitsetSaver.hh"

namespace isaac
{
namespace io
{

BitsetSaver::BitsetSaver(const boost::filesystem::path filePath) :
    filePath_(filePath),
    os_(filePath.c_str()),
    pos_(0),
    byte_(0)
{
    if (!os_)
    {
        const boost::format message = boost::format("Failed to open file %s for writing: %s") % filePath_ % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

void BitsetSaver::save(
    std::vector<bool>::const_iterator begin,
    const std::vector<bool>::const_iterator end,
    const bool final)
{
    while (end != begin)
    {
        const char positionHasNeighbors = *begin;
        const unsigned shift = (pos_++) % 8;
        byte_ |= positionHasNeighbors << shift;
        if (7 == shift)
        {
            if (!os_.write(&byte_, sizeof(byte_)))
            {
                const boost::format message = boost::format("Failed to write bits into %s: %s") % filePath_.string() % strerror(errno);
                BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
            }
            byte_ = 0;
        }
        ++begin;
    }

    if (final && !os_.write(&byte_, sizeof(byte_)))
    {
        const boost::format message = boost::format("Failed to write final bits byte into %s: %s") % filePath_.string() % strerror(errno);
        BOOST_THROW_EXCEPTION(common::IoException(errno, message.str()));
    }
}

} // namespace reference
} // namespace isaac
