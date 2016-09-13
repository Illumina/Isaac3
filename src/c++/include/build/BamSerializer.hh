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
 ** \file BamSerializer.hh
 **
 ** Helper class for converting BinSorter data to bam.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BAM_SERIALIZER_HH
#define iSAAC_BUILD_BAM_SERIALIZER_HH

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "build/BarcodeBamMapping.hh"
#include "build/BinData.hh"
#include "build/FragmentIndex.hh"
#include "build/PackedFragmentBuffer.hh"
#include "flowcell/TileMetadata.hh"


namespace isaac
{
namespace build
{

class BamSerializer
{
public:
    BamSerializer(
        const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeOutputFileIndexMap,
        const unsigned splitGapLength):
            splitGapLength_(splitGapLength),
            barcodeOutputFileIndexMap_(barcodeOutputFileIndexMap)
    {}

    void storeAligned(
        const io::FragmentAccessor &fragment,
        boost::ptr_vector<boost::iostreams::filtering_ostream> &streams,
        boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
        FragmentAccessorBamAdapter& adapter)
    {
        unsigned serializedLength = bam::serializeAlignment(
            streams.at(barcodeOutputFileIndexMap_.at(fragment.barcode_)), adapter);
        bam::BamIndexPart& bamIndexPart = bamIndexParts.at(barcodeOutputFileIndexMap_.at(fragment.barcode_));
        bamIndexPart.processFragment( adapter, serializedLength );
//        ISAAC_THREAD_CERR << "Serialized to bam pos_: " << idx.pos_ << " dataOffset_: " << idx.dataOffset_ << std::endl;
    }

    void storeUnaligned(
        const io::FragmentAccessor &fragment,
        boost::ptr_vector<boost::iostreams::filtering_ostream> &streams,
        boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
        FragmentAccessorBamAdapter& adapter)
    {
        unsigned serializedLength = bam::serializeAlignment(
            streams.at(barcodeOutputFileIndexMap_.at(fragment.barcode_)), adapter);
        bam::BamIndexPart& bamIndexPart = bamIndexParts.at(barcodeOutputFileIndexMap_.at(fragment.barcode_));
        bamIndexPart.processFragment( adapter, serializedLength );
//        ISAAC_THREAD_CERR << "Serialized unaligned pos_: " << fragment << std::endl;
    }

    void prepareForBam(
        const reference::ContigList &contigList,
        PackedFragmentBuffer &data,
        BinData::IndexType &dataIndex,
        alignment::Cigar &splitCigars,
        SplitInfoList &splitInfoList);

private:
    void splitIfNeeded(
        const reference::ContigList &contigList,
        PackedFragmentBuffer &data,
        PackedFragmentBuffer::Index &index,
        BinData::IndexType &splitIndexEntries,
        alignment::Cigar &splitCigars,
        SplitInfoList &splitInfoList);

    alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> makeSplit(
        const alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator>& last,
        alignment::CigarPosition<PackedFragmentBuffer::Index::CigarIterator> current,
        PackedFragmentBuffer& data,
        PackedFragmentBuffer::Index& index,
        const unsigned editDistance,
        BinData::IndexType& splitIndexEntries,
        alignment::Cigar& splitCigars,
        SplitInfoList &splitInfoList);

    const unsigned splitGapLength_;
    const BarcodeBamMapping::BarcodeSampleIndexMap &barcodeOutputFileIndexMap_;
};

} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BAM_SERIALIZER_HH
