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
 **/

#ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH
#define iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/SeedMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/Match.hh"

//struct TestMatchStorage : public  std::vector<isaac::alignment::Match>
//{
//    template <typename IteratorT>
//    void storeMatches(IteratorT begin, IteratorT end)
//    {
//        insert(this->end(), begin, end);
//    }
//};

typedef isaac::alignment::Matches TestMatchStorage;


class TestHashMatchFinder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestHashMatchFinder );
    CPPUNIT_TEST( testEverything );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList_;
    const isaac::alignment::SeedMetadataList seedMetadataList_;

public:
    TestHashMatchFinder();
    void setUp();
    void tearDown();
    void testEverything();

private:
    TestMatchStorage findMatches(
        const std::string& reference,
        const std::string& sequence,
        const unsigned repeatThreshold,
        const isaac::alignment::SeedMetadataList &seedMetadataList,
        const isaac::flowcell::ReadMetadataList &readMetadataList);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SEQUENCING_ADAPTER_HH

