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
 ** \file BitsetSaver.hh
 **
 ** \brief Helper for saving neighbor flags and such from a binary file
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_BITSET_SAVER_HH
#define ISAAC_REFERENCE_BITSET_SAVER_HH

#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>

namespace isaac
{
namespace io
{

class BitsetSaver: boost::noncopyable
{
    boost::filesystem::path filePath_;
    std::ofstream os_;
    std::size_t pos_;
    char byte_;
public:
    BitsetSaver(const boost::filesystem::path filePath);
    void save(
        const std::vector<bool>::const_iterator begin,
        const std::vector<bool>::const_iterator end,
        const bool final);
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_BITSET_SAVER_HH
