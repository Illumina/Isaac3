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
 ** \file DirectFragmentStorage.hh
 **
 ** \brief Immediate storing of aligned data in the output bin files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH

#include "alignment/BamTemplate.hh"
#include "alignment/BinMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

/**
 * \brief interface for various fragment storage implementations
 *        TODO: remove this when the implementation settles
 */
class FragmentStorage
{
public:
    virtual ~FragmentStorage(){}

    virtual void store(
        const BamTemplate &bamTemplate,
        const unsigned barcodeIdx) = 0;
    virtual void reset(const uint64_t clusterId, const bool paired) = 0;

    virtual void prepareFlush() noexcept = 0;
    virtual void flush() = 0;
    virtual void resize(const uint64_t clusters) = 0;
    virtual void reserve(const uint64_t clusters) = 0;
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_STORAGE_HH
