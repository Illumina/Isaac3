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
 ** \file BinMetadata.hh
 **
 ** \brief Metadata associated to the unsorted alignment results
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
#define iSAAC_ALIGNMENT_BIN_METADATA_HH

#include <numeric>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "reference/ReferencePosition.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{

struct BinMetadata;
inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata);

struct BarcodeCounts
{
    BarcodeCounts() : elements_(0), gaps_(0), splits_(0), cigarLength_(0){}
    // total number of elements in the bin barcode
    uint64_t elements_;
    // total number of gaps in the bin barcode reads
    uint64_t gaps_;
    // total number of splits in the bin barcode reads. Splits are count normally once except for FLIPs which could
    // count twice unless FLIP is not followed by CONTIG or position adjustment
    uint64_t splits_;
    // sum of all fragment cigar lengths in the bin barcode.
    uint64_t cigarLength_;

    BarcodeCounts& operator += (const BarcodeCounts &that)
    {
        elements_ += that.elements_;
        gaps_ += that.gaps_;
        splits_ += that.splits_;
        cigarLength_ += that.cigarLength_;
        return *this;
    }

    friend BarcodeCounts operator + (const BarcodeCounts &left, const BarcodeCounts &right)
    {
        BarcodeCounts ret(left);
        return ret += right;
    }
};



struct BinChunk
{
    std::vector<BarcodeCounts> barcodeBreakdown_;
    uint64_t dataSize_;

    // required for serialization
    BinChunk(): dataSize_(0){}
    BinChunk(const unsigned barcodesCount) : barcodeBreakdown_(barcodesCount), dataSize_(0){}

    static void reset(BinChunk &chunk)
    {
        chunk.dataSize_ = 0;
        for (BarcodeCounts &barcode : chunk.barcodeBreakdown_) {barcode = BarcodeCounts(); }
    }

    void swap(BinChunk &that)
    {
        if (this != &that)
        {
            std::swap(dataSize_ , that.dataSize_);
            barcodeBreakdown_.swap(that.barcodeBreakdown_);
        }
    }

    BinChunk &operator += (const BinChunk &that)
    {
        dataSize_ += that.dataSize_;
        std::transform(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), that.barcodeBreakdown_.begin(), barcodeBreakdown_.begin(), std::plus<BarcodeCounts>());
        return *this;
    }

    uint64_t getTotalCigarLength() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::cigarLength_, _2)));
    }

    uint64_t getTotalElements() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::elements_, _2)));
    }

    uint64_t getBarcodeGapCount(const unsigned barcodeIdx) const
    {
        ISAAC_ASSERT_MSG(barcodeBreakdown_.size() > barcodeIdx, "Invalid barcode requested: " << barcodeIdx << " for size: " << barcodeBreakdown_.size());
        return barcodeBreakdown_[barcodeIdx].gaps_;
    }

    uint64_t getTotalGapCount() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::gaps_, _2)));
    }

    uint64_t getTotalSplitCount() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::splits_, _2)));
    }

    uint64_t getBarcodeElements(const unsigned barcodeIdx) const
        {return barcodeBreakdown_.at(barcodeIdx).elements_;}

    void incrementGapCount(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).gaps_ += by;}

    void incrementSplitCount(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).splits_ += by;}

    void incrementCigarLength(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).cigarLength_ += by;}

    void incrementSeIdxElements(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementRIdxElements(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementFIdxElements(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}

    void incrementNmElements(const uint64_t by, const unsigned barcodeIdx)
        {barcodeBreakdown_.at(barcodeIdx).elements_ += by;}
};

inline std::ostream &operator << (std::ostream &os, const BinChunk &chunk)
{
    return os << "BinChunk(" <<
        chunk.dataSize_ << "ds" <<
        ")";
}

class BinDataDistribution : std::vector<BinChunk>
{
    uint64_t chunkSize_;
    bool offsetsTallied_;

    friend std::ostream &operator << (std::ostream &os, const BinDataDistribution &distribution);

public:
    using std::vector<BinChunk>::reserve;
    using std::vector<BinChunk>::size;
    BinDataDistribution(
        const unsigned barcodesCount,
        const uint64_t length,
        const unsigned distributionChunksCount)
        // we need one for rounding and one more chunk so that tallyOffset produces the end offset for the last present chunk. Last chunk should be always empty
        :std::vector<BinChunk>(length / getChunkSize(length, distributionChunksCount) + 2UL, BinChunk(barcodesCount)),
         chunkSize_(getChunkSize(length, distributionChunksCount)), offsetsTallied_(false)
    {
    }

    void makeChunksBigger(const uint64_t newBinLength)
    {
        ISAAC_ASSERT_MSG(!empty(), "Unexpected empty bin " << *this);
        const uint64_t oldChunkSize = chunkSize_;
        chunkSize_ = getChunkSize(newBinLength, size() - 2);
        ISAAC_THREAD_CERR << "increasing chunkSize to " << chunkSize_ << " from " << oldChunkSize << " for newBinLength: " << newBinLength << " last:" << back() << std::endl;
        ISAAC_ASSERT_MSG(oldChunkSize < chunkSize_, "New chunk size " << chunkSize_ << " must be greater than " << oldChunkSize);
        iterator merged = begin();
        iterator current = begin();
        for (; end() - 1 != current;)
        {
            merged->swap(*current);
            ++current;
            for (uint64_t s = oldChunkSize; s < chunkSize_ && end() - 1 != current; s += oldChunkSize)
            {
                *merged += *current;
                ++current;
            }
            ++merged;
        }
        ISAAC_ASSERT_MSG(end() - 1 == current, "Did not reach the end." << std::distance(begin(), current) << "/" << size());
        ISAAC_ASSERT_MSG(!current->dataSize_, "Expected last chunk to be empty: " << *current << " for " << *this);
        std::for_each(merged, end(), &BinChunk::reset);
        ISAAC_THREAD_CERR << "increased chunkSize to " << chunkSize_ << " from " << oldChunkSize << " for newBinLength: " << newBinLength << " last:" << back() << std::endl;
    }

    uint64_t addBytes(uint64_t binGenomicOffset, unsigned count);
    uint64_t tallyOffsets();

    BinDataDistribution &operator =(const BinDataDistribution &that);

    BinChunk &getChunk(uint64_t binGenomicOffset)
        {return this->operator[](getIndex(binGenomicOffset));}

    std::size_t getIndex(uint64_t binGenomicOffset) const;
    uint64_t getChunkSize() const {return chunkSize_;}

    /// if distributionChunksCount is 0, assume extremely large chunk to ensure that all data fits into chunk 0
    static uint64_t getChunkSize(
        const uint64_t length,
        const unsigned distributionChunksCount)
        {return distributionChunksCount ?  std::max(1UL, length / distributionChunksCount) : -1UL;}

    uint64_t getChunkEndOffset(std::size_t chunk);

    uint64_t getTotalCigarLength() const;
    uint64_t getBarcodeGapCount(const unsigned barcodeIdx) const;
    uint64_t getTotalGapCount() const;
    uint64_t getTotalSplitCount() const;
    uint64_t getBarcodeElements(const unsigned barcodeIdx) const;
    uint64_t getTotalElements() const;

    void incrementCigarLength(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementCigarLength(by, barcodeIdx);}

    void incrementGapCount(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementGapCount(by, barcodeIdx);}

    void incrementSplitCount(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementSplitCount(by, barcodeIdx);}

    void incrementSeIdxElements(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementSeIdxElements(by, barcodeIdx);}

    void incrementRIdxElements(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementRIdxElements(by, barcodeIdx);}

    void incrementFIdxElements(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementFIdxElements(by, barcodeIdx);}

    void incrementNmElements(const uint64_t binGenomicOffset, const uint64_t by, const unsigned barcodeIdx)
        {getChunk(binGenomicOffset).incrementNmElements(by, barcodeIdx);}

    uint64_t removeChunksBefore(const uint64_t minOffset);
    uint64_t removeChunksAfter(const uint64_t minOffset);

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, BinDataDistribution &bm, const unsigned int version);
};

inline std::ostream &operator << (std::ostream &os, const BinDataDistribution &distribution)
{
    os << "BinDataDistribution(" <<
        distribution.chunkSize_ << "cs," <<
        distribution.offsetsTallied_ << "ot,[";
    static const unsigned PRINT_CHUNKS_MAX = 10;
    unsigned toPrint = PRINT_CHUNKS_MAX;
    const std::vector<BinChunk> &dist = distribution;
    BOOST_FOREACH(const BinChunk &chunk, dist)
    {
        os << chunk;
        if (!--toPrint)
        {
            os << "...";
            break;
        }
        os << ",";
    }
    return os << "])";
}

class BinMetadata
{
    unsigned binIndex_;
    /// first genomic position covered by the bin
    reference::ReferencePosition binStart_;
    /// bin length in bases
    uint64_t length_;
    boost::filesystem::path binFilePath_;
    // offset from the beginning of the data file.
    // Note that single file can later be broken down into multiple BinMetadata objects
    uint64_t dataOffset_;
    // number of bytes stored in binFilePath_ at dataOffset_
    uint64_t dataSize_;
    uint64_t seIdxElements_;
    uint64_t rIdxElements_;
    uint64_t fIdxElements_;
    uint64_t nmElements_;

    BinDataDistribution dataDistribution_;

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, BinMetadata &bm, const unsigned int version);

public:
    BinMetadata() :
        binIndex_(0),
        binStart_(0),
        length_(0),
        dataOffset_(0),
        dataSize_(0),
        seIdxElements_(0),
        rIdxElements_(0),
        fIdxElements_(0),
        nmElements_(0),
        dataDistribution_(0,0,0){}

    BinMetadata(
        const unsigned barcodesCount,
        const unsigned binIndex,
        const reference::ReferencePosition binStart,
        const uint64_t length,
        const boost::filesystem::path &binFilepath,
        const unsigned distributionChunksCount) :
            binIndex_(binIndex),
            binStart_(binStart),
            length_(length),
            binFilePath_(binFilepath),
            dataOffset_(0),
            dataSize_(0),
            seIdxElements_(0),
            rIdxElements_(0),
            fIdxElements_(0),
            nmElements_(0),
            dataDistribution_(barcodesCount, length_, distributionChunksCount){}

    /**
     * \return BinMedata which guarantees to have the chunks with
     *         offset: minOffset <= offset < (minOffset + minSize)
     */
    BinMetadata getChunks(
        const uint64_t minOffset,
        const uint64_t minSize) const
    {
        ISAAC_ASSERT_MSG(isUnalignedBin(), "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!rIdxElements_, "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!fIdxElements_, "Splitting bins is supported only for unaligned bin");
        ISAAC_ASSERT_MSG(!seIdxElements_, "Splitting bins is supported only for unaligned bin");

        BinMetadata ret(*this);
        ret.removeChunksBefore(minOffset);
        ret.removeChunksAfter(minSize);

        return ret;
    }

    void removeChunksBefore(const uint64_t minOffset)
    {
        const uint64_t removedBytes = dataDistribution_.removeChunksBefore(minOffset);
        dataOffset_ += removedBytes;
        dataSize_ -= removedBytes;
    }

    void removeChunksAfter(const uint64_t minOffset)
    {
        dataSize_ = dataDistribution_.removeChunksAfter(minOffset);
    }

    unsigned getIndex() const
    {
        return binIndex_;
    }

    reference::ReferencePosition getBinStart() const
    {
        return binStart_;
    }

    bool hasPosition(const reference::ReferencePosition pos) const
    {
        return binStart_ <= pos && pos < (binStart_ + length_);
    }

    reference::ReferencePosition getBinEnd() const
    {
        return isUnalignedBin() ? reference::ReferencePosition(reference::ReferencePosition::NoMatch) : binStart_ + length_;
    }

    bool coversPosition(const reference::ReferencePosition pos) const
    {
        ISAAC_ASSERT_MSG(!isUnalignedBin(), "Checking positions in unaligned bins is not allowed");
        return getBinStart() <= pos && getBinEnd() > pos;
    }

    bool isUnalignedBin() const {return binStart_.isTooManyMatch();}

    const boost::filesystem::path & getPath() const
    {
        return binFilePath_;
    }
    const std::string &getPathString() const
    {
        return binFilePath_.string();
    }

    uint64_t getDataOffset() const
    {
        return dataOffset_;
    }

    uint64_t getDataSize() const
    {
        return dataSize_;
    }

    bool isEmpty() const
    {
        return 0 == dataSize_;
    }

    uint64_t getLength() const {return length_;}

    /**
     * \brief increment the the corresponding chunk size and total data size.
     *
     * \return pair of total data size before increment and data offset of the corresponding chunk prior to being incremented
     */
    std::pair<uint64_t, uint64_t> incrementDataSize(const reference::ReferencePosition pos, const uint64_t by)
    {
        const std::pair<uint64_t, uint64_t> ret = std::make_pair(dataSize_, dataDistribution_.addBytes(getDataDistributionKey(pos), by));
        dataSize_ += by;
        return ret;
    }

    /**
     * \brief increment the the corresponding chunk size and total data size.
     *
     * \return pair of total data size before increment and data offset of the corresponding chunk prior to being incremented
     */
    std::pair<uint64_t, uint64_t> incrementDataSize(const uint64_t recordNumber, const uint64_t by)
    {
        const std::pair<uint64_t, uint64_t> ret = std::make_pair(dataSize_, dataDistribution_.addBytes(getDataDistributionKey(recordNumber), by));
        dataSize_ += by;
        return ret;
    }

    uint64_t getDataDistributionKey(const uint64_t clusterNumber)
    {
        ISAAC_ASSERT_MSG(isUnalignedBin(), "Aligned bins must use ReferencePosition as hash key." << *this);
        if (clusterNumber >= length_)
        {
            length_ *= 2;
            dataDistribution_.makeChunksBigger(length_);
        }
        return clusterNumber;
    }

    uint64_t getDataDistributionKey(const reference::ReferencePosition pos)
    {
        ISAAC_ASSERT_MSG(!isUnalignedBin(), "Unaligned bins must use cluster number as key." << *this);
        if (hasPosition(pos))
        {
            return pos - binStart_;
        }
        // We store some reads that don't belong to the bin (mates from other bins), as well as reads that belong to multiple bins (split reads)
        // Put them in the first chunk of the bin
        return 0;
    }

    uint64_t getSeIdxElements() const
    {
        return seIdxElements_;
    }

    void incrementSeIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementSeIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        seIdxElements_ += by;
    }

    uint64_t getRIdxElements() const
    {
        return rIdxElements_;
    }

    void incrementRIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementRIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        rIdxElements_ += by;
    }

    uint64_t getFIdxElements() const
    {
        return fIdxElements_;
    }

    void incrementFIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementFIdxElements(getDataDistributionKey(pos), by, barcodeIdx);
        fIdxElements_ += by;
    }

    uint64_t getNmElements() const
    {
        return nmElements_;
    }

    void incrementNmElements(const uint64_t sequenceHash, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementNmElements(getDataDistributionKey(sequenceHash), by, barcodeIdx);
        nmElements_ += by;
    }

    void incrementGapCount(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementGapCount(getDataDistributionKey(pos), by, barcodeIdx);
    }

    void incrementSplitCount(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementSplitCount(getDataDistributionKey(pos), by, barcodeIdx);
    }

    void incrementCigarLength(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        dataDistribution_.incrementCigarLength(getDataDistributionKey(pos), by, barcodeIdx);
    }

    uint64_t getTotalElements() const
    {
        //return  getSeIdxElements() + getRIdxElements() + getFIdxElements() + getNmElements();
        return dataDistribution_.getTotalElements();
    }

    uint64_t getBarcodeElements(const unsigned barcodeIdx) const
    {
        return dataDistribution_.getBarcodeElements(barcodeIdx);
    }

    uint64_t getBarcodeGapCount(const unsigned barcodeIdx) const
    {
        return dataDistribution_.getBarcodeGapCount(barcodeIdx);
    }

    uint64_t getTotalGapCount() const
    {
        return dataDistribution_.getTotalGapCount();
    }

    uint64_t getTotalSplitCount() const
    {
        return dataDistribution_.getTotalSplitCount();
    }

    uint64_t getTotalCigarLength() const
    {
        return dataDistribution_.getTotalCigarLength();
    }

    const BinDataDistribution &getDataDistribution() const
    {
        return dataDistribution_;
    }
};


struct BinMetadataList : public std::vector<alignment::BinMetadata>
{
    BinMetadataList(){}
    BinMetadataList(size_t size) : std::vector<alignment::BinMetadata>(size){}
};
typedef boost::reference_wrapper<const BinMetadata> BinMetadataCRef;
typedef std::vector<BinMetadataCRef >BinMetadataCRefList;


inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata)
{
    return os << "BinMetadata("
              << binMetadata.getIndex() << "id "
              << binMetadata.getBinStart() << "bs "
              << binMetadata.getLength() << "bl "
              << binMetadata.getDataSize() << "ds "
              << binMetadata.getDataOffset() << "do "
              << binMetadata.getSeIdxElements() << "se "
              << binMetadata.getRIdxElements() << "rs "
              << binMetadata.getFIdxElements() << "f "
              << binMetadata.getPathString() << ")";
}


/**
 * \param binOffset Offset relative to the beginning of the Bin
 *
 * \return data offset prior to incrementing
 */
inline uint64_t BinDataDistribution::addBytes(uint64_t binGenomicOffset, unsigned count)
{
    BinChunk &ref = getChunk(binGenomicOffset);
    const uint64_t ret = ref.dataSize_;
//        ISAAC_THREAD_CERR << "addBytes binOffset:" << binOffset << "getChunkSize()" << getChunkSize() <<
//            " count:" << count <<
//            " binIndex: " << binIndex << " ret:" << ret << std::endl;
    ref.dataSize_ += count;
    return ret;
}

/**
 * \brief replace contents of each chunk with sum of contents of all previous chunks
 *
 * \return the total number of bytes occupied by data
 */
inline uint64_t BinDataDistribution::tallyOffsets()
{
    ISAAC_ASSERT_MSG(!offsetsTallied_, "Offsets cannot be gallied twice." << *this);
    uint64_t offset = 0;
    std::vector<BinChunk> &v = *this;
    BOOST_FOREACH(BinChunk &chunk, v)
    {
        using std::swap;
        swap(offset, chunk.dataSize_);
//            ISAAC_THREAD_CERR << "tallyOffsets: " << chunkOffset << std::endl;
        offset += chunk.dataSize_;
    }
    offsetsTallied_ = true;
    return offset;
}

inline BinDataDistribution &BinDataDistribution::operator =(const BinDataDistribution &that)
{
    assign(that.begin(), that.end());
    offsetsTallied_ = that.offsetsTallied_;
    return *this;
}

inline std::size_t BinDataDistribution::getIndex(uint64_t key) const
{
    std::size_t ret = key / getChunkSize();
    ISAAC_ASSERT_MSG(size() > ret, "chunk index " << ret << " for offset " << key << "is too big " << *this);
    return ret;
}

inline uint64_t BinDataDistribution::getChunkEndOffset(std::size_t chunk)
{
    ISAAC_ASSERT_MSG(offsetsTallied_, "getChunkEndOffset for untallied distribution");
    return size() > (chunk + 1) ? at(chunk+1).dataSize_ : back().dataSize_;
}

inline uint64_t BinDataDistribution::getTotalCigarLength() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getTotalCigarLength, _2)));
}

inline uint64_t BinDataDistribution::getBarcodeGapCount(const unsigned barcodeIdx) const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getBarcodeGapCount, _2, barcodeIdx)));
}

inline uint64_t BinDataDistribution::getTotalGapCount() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getTotalGapCount, _2)));
}

inline uint64_t BinDataDistribution::getTotalSplitCount() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getTotalSplitCount, _2)));
}

inline uint64_t BinDataDistribution::getBarcodeElements(const unsigned barcodeIdx) const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getBarcodeElements, _2, barcodeIdx)));
}

inline uint64_t BinDataDistribution::getTotalElements() const
{
    return std::accumulate(begin(), end(), 0,
                           boost::bind(std::plus<uint64_t>(),
                                       _1, boost::bind(&BinChunk::getTotalElements, _2)));
}

/**
 * \return number of bytes removed
 */
inline uint64_t BinDataDistribution::removeChunksBefore(const uint64_t minOffset)
{
    uint64_t offset = 0;
    std::vector<BinChunk>::iterator e = begin();
    for (; end() != e && minOffset > offset; ++e)
    {
        ISAAC_THREAD_CERR << "removeChunksBefore minOffset: " << minOffset << " " << *e << std::endl;
        offset += e->dataSize_;
    }

    erase(begin(), e);
    return offset;
}

/**
 * \return number of bytes left
 */
inline uint64_t BinDataDistribution::removeChunksAfter(const uint64_t minOffset)
{
    uint64_t offset = 0;
    std::vector<BinChunk>::iterator b = begin();
    for (; end() != b && minOffset > offset; ++b)
    {
        offset += b->dataSize_;
    }

    erase(b, end());
    return offset;
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
