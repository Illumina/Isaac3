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
 ** \file SeedId.hh
 **
 ** \brief Identification of a seed encoding offset, length and
 ** orientation.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SEED_ID_HH
#define iSAAC_ALIGNMENT_SEED_ID_HH

#include <cassert>
#include <iostream>
#include <boost/format.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>

#include "common/Exceptions.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Structured unique identifier of a seed.
 **
 ** From LSB to MSB:
 **   - reverse:  1 (            2)
 **   - seed   :  8 (          256)
 **   - length : 7  (          128)
 **   - barcode: 12 (        4,096)
 **   - tile   : 12 (        4,096)
 **
 **
 ** The order of the fields is important as it will definet the natural order of
 ** the seed ids when sorted.
 **
 ** \note The implementation choice at the moment is to use a tile number that
 ** is the sequential index (0-based) of the tile in the input data set. This
 ** solution has been selected instead of keeping the actual tile and lane, to
 ** make easier (more flexible) to retrieve the metadata and components
 ** associated to each tile throughout the application.
 **
 ** \sa isaac::demultiplexing::BarcodeId
 **/
class SeedId
{
public:
    // width in bits for each field
    static const unsigned REVERSE_WIDTH = 1;
    static const unsigned SEED_LENGTH_WIDTH = 7;

//    static const unsigned UNUSED_SEED_WIDTH = 8;
//    static const unsigned UNUSED_BARCODE_WIDTH = 12;
//    static const unsigned UNUSED_TILE_WIDTH = 36;
    // masks for the values in each field
    static const unsigned char REVERSE_MASK = ~(~0UL<<REVERSE_WIDTH);
    static const unsigned char SEED_LENGTH_MASK = ~(~0UL<<SEED_LENGTH_WIDTH);

//    static const uint64_t UNUSED_SEED_MASK = ~(~0UL<<UNUSED_SEED_WIDTH);
//    static const uint64_t UNUSED_BARCODE_MASK = ~(~0UL<<UNUSED_BARCODE_WIDTH);
//    static const uint64_t UNUSED_TILE_MASK = ~(~0UL<<UNUSED_TILE_WIDTH);
    // shifts in bits for each field
    static const unsigned REVERSE_SHIFT = 0;
    static const unsigned SEED_LENGTH_SHIFT = REVERSE_SHIFT + REVERSE_WIDTH;

//    static const unsigned UNUSED_SEED_SHIFT = SEED_LENGTH_SHIFT + SEED_LENGTH_WIDTH;
//    static const unsigned UNUSED_BARCODE_SHIFT = UNUSED_SEED_SHIFT + UNUSED_SEED_WIDTH;
//    static const unsigned UNUSED_TILE_SHIFT = UNUSED_BARCODE_SHIFT + UNUSED_BARCODE_WIDTH;
    explicit SeedId(uint64_t value = 0) : value_(value) {}
    SeedId(uint64_t seedLength, uint64_t reverse)
        : value_(
//            ((0 & UNUSED_TILE_MASK) << UNUSED_TILE_SHIFT) |
//            ((0 & UNUSED_BARCODE_MASK) << UNUSED_BARCODE_SHIFT) |
//            ((0 & UNUSED_SEED_MASK) << UNUSED_SEED_SHIFT) |
            ((seedLength & SEED_LENGTH_MASK) << SEED_LENGTH_SHIFT) |
            ((reverse & REVERSE_MASK) << REVERSE_SHIFT))
    {
        using boost::mpl::equal_to;
        using boost::mpl::int_;
//        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<REVERSE_WIDTH + SEED_LENGTH_WIDTH + UNUSED_SEED_WIDTH + UNUSED_BARCODE_WIDTH + UNUSED_TILE_WIDTH> >));
//        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<UNUSED_TILE_WIDTH + UNUSED_TILE_SHIFT> >));
        BOOST_MPL_ASSERT((equal_to<int_<8>, int_<REVERSE_WIDTH + SEED_LENGTH_WIDTH> >));
        using boost::format;
        using isaac::common::PreConditionException;
        if (
//            (UNUSED_TILE_MASK < 0) ||
//            (UNUSED_BARCODE_MASK < 0) ||
//            (UNUSED_SEED_MASK < seedLength) ||
            (SEED_LENGTH_MASK < 0) ||
            (REVERSE_MASK < reverse))
        {
            using boost::format;
            using isaac::common::PreConditionException;
            const format message = format(
                "SeqId(%ld, %ld): maximum values are (%ld, %ld)") %
                seedLength % reverse %
                (uint64_t)SEED_LENGTH_MASK %
                (uint64_t)REVERSE_MASK;
            BOOST_THROW_EXCEPTION(PreConditionException(message.str()));
        }
    }

    SeedId(uint64_t tile, uint64_t barcode, uint64_t seedLength, uint64_t seed, uint64_t reverse)
        : value_(
//            ((tile  & UNUSED_TILE_MASK) << UNUSED_TILE_SHIFT) |
//            ((barcode & UNUSED_BARCODE_MASK) << UNUSED_BARCODE_SHIFT) |
//            ((seed & UNUSED_SEED_MASK) << UNUSED_SEED_SHIFT) |
            ((seedLength & SEED_LENGTH_MASK) << SEED_LENGTH_SHIFT) |
            ((reverse & REVERSE_MASK) << REVERSE_SHIFT))
    {
        using boost::mpl::equal_to;
        using boost::mpl::int_;
//        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<REVERSE_WIDTH + SEED_LENGTH_WIDTH + UNUSED_SEED_WIDTH + UNUSED_BARCODE_WIDTH + UNUSED_TILE_WIDTH> >));
//        BOOST_MPL_ASSERT((equal_to<int_<64>, int_<UNUSED_TILE_WIDTH + UNUSED_TILE_SHIFT> >));
        BOOST_MPL_ASSERT((equal_to<int_<8>, int_<REVERSE_WIDTH + SEED_LENGTH_WIDTH> >));
        using boost::format;
        using isaac::common::PreConditionException;
        if (
//            (UNUSED_TILE_MASK < tile) ||
//            (UNUSED_BARCODE_MASK < barcode) ||
//            (UNUSED_SEED_MASK < seed) ||
            (SEED_LENGTH_MASK < seedLength) ||
            (REVERSE_MASK < reverse))
        {
            using boost::format;
            using isaac::common::PreConditionException;
            const format message = format(
                "SeqId("
//                "%ld, %ld, "
                "%ld, "
//                "%ld, "
                "%ld): maximum values are ("
//                "%ld, %ld, "
                "%ld, "
//                "%ld, "
                "%ld)") %
//                tile % barcode %
                seedLength %
//                seed %
                reverse %
                //(uint64_t)UNUSED_TILE_MASK % (uint64_t)UNUSED_BARCODE_MASK %
                (uint64_t)SEED_LENGTH_MASK %
                //(uint64_t)UNUSED_SEED_MASK %
                (uint64_t)REVERSE_MASK;
            BOOST_THROW_EXCEPTION(PreConditionException(message.str()));
        }
    }
//    uint64_t unusedgetTile() const {return (value_ >> UNUSED_TILE_SHIFT) & UNUSED_TILE_MASK;}
//    uint64_t unusedgetBarcode() const {return (value_ >> UNUSED_BARCODE_SHIFT) & UNUSED_BARCODE_MASK;}
//    uint64_t unusedgetSeed() const {return (value_ >> UNUSED_SEED_SHIFT) & UNUSED_SEED_MASK;}

    uint64_t getSeedLength() const {return (value_ >> SEED_LENGTH_SHIFT) & SEED_LENGTH_MASK;}
    uint64_t getReverse() const {return (value_ >> REVERSE_SHIFT) & REVERSE_MASK;}
    bool isReverse() const {return getReverse();}
    void setForward() {value_ &= ~(REVERSE_MASK << REVERSE_SHIFT);}
    void setReverse() {value_ |= (REVERSE_MASK << REVERSE_SHIFT);}
    void invert() {value_ ^= (REVERSE_MASK << REVERSE_SHIFT);}
    operator uint64_t() const {return value_;}
private:
    unsigned char value_;
};


inline std::ostream &operator<<(std::ostream &os, const SeedId &s)
{
    return os << "SeedId(" <<
        //s.getTile() << ":" << s.getBarcode() << ":"<<
        s.getSeedLength() << ":" <<
        //s.getSeed() << ":" <<
        (s.getReverse() ? 'r':'f') << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SEED_ID_HH
