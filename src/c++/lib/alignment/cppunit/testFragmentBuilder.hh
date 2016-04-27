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

#ifndef iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH
#define iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>
#include <vector>

#include "alignment/FragmentBuilder.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/Layout.hh"

class TestFragmentBuilder : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestFragmentBuilder );
//    CPPUNIT_TEST( testEmptyMatchList );
//    CPPUNIT_TEST( testSingleSeed );
//    CPPUNIT_TEST( testSeedOffset );
//    CPPUNIT_TEST( testMultiSeed );
//    CPPUNIT_TEST( testRepeats );
//    CPPUNIT_TEST( testMismatches );
//    CPPUNIT_TEST( testLeadingSoftClips );
//    CPPUNIT_TEST( testTrailingSoftClips );
//    CPPUNIT_TEST( testLeadingAndTrailingSoftClips );
    CPPUNIT_TEST_SUITE_END();
private:
    const isaac::flowcell::ReadMetadataList readMetadataList;
    const isaac::alignment::SeedMetadataList seedMetadataList;
    const isaac::flowcell::FlowcellLayoutList flowcells;
    const isaac::reference::ContigList contigList;
    isaac::reference::ContigAnnotations contigAnnotations;
    const isaac::alignment::BclClusters bcl0;
    const isaac::alignment::BclClusters bcl2;
    const isaac::alignment::BclClusters bcl3;
    const isaac::alignment::BclClusters bcl4l;
    const isaac::alignment::BclClusters bcl4t;
    const isaac::alignment::BclClusters bcl4lt;
    const unsigned tile0;
    const unsigned tile2;
    const unsigned clusterId0;
    const unsigned clusterId2;
    isaac::alignment::Cluster cluster0;
    isaac::alignment::Cluster cluster2;
    isaac::alignment::Cluster cluster3;
    isaac::alignment::Cluster cluster4l;
    isaac::alignment::Cluster cluster4t;
    isaac::alignment::Cluster cluster4lt;
    std::vector<isaac::alignment::Match> matchList;
    // auxiliary method for tests dedicated to a single seed on each read
    void auxSingleSeed(const unsigned s0, const unsigned s1);
public:
    TestFragmentBuilder();
    void setUp();
    void tearDown();
    void testEmptyMatchList();
    void testSingleSeed();
    void testSeedOffset();
    void testMultiSeed();
    void testRepeats();
    void testMismatches();
    void testLeadingSoftClips();
    void testTrailingSoftClips();
    void testLeadingAndTrailingSoftClips();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_FRAGMENT_BUILDER_HH
