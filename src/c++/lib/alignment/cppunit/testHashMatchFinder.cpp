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
 ** \file TestHashMatchFinder.cpp
 **
 ** More fragment builder tests.
 **
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/construct.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testHashMatchFinder.hh"

#include "alignment/HashMatchFinder.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "BuilderInit.hh"
#include "common/Threads.hpp"
#include "reference/SortedReferenceMetadata.hh"
#include "reference/ReferenceHash.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestHashMatchFinder, registryName("HashMatchFinder"));

TestHashMatchFinder::TestHashMatchFinder() :
    readMetadataList_(boost::assign::list_of
                     (isaac::flowcell::ReadMetadata(1, 50, 0, 0)).convert_to_container<std::vector<isaac::flowcell::ReadMetadata> >()),
    seedMetadataList_(boost::assign::list_of
                     (isaac::alignment::SeedMetadata( 0, 8, 0, 0))
                     (isaac::alignment::SeedMetadata(34, 8, 0, 1)).convert_to_container<isaac::alignment::SeedMetadataList>()
                     )
{

}

void TestHashMatchFinder::setUp()
{
}

void TestHashMatchFinder::tearDown()
{
}


TestMatchStorage TestHashMatchFinder::findMatches(
    const std::string& reference, const std::string& sequence,
    const unsigned repeatThreshold,
    const isaac::alignment::SeedMetadataList &seedMetadataList,
    const isaac::flowcell::ReadMetadataList &readMetadataList)
{
    isaac::reference::ContigList contigList(1, isaac::reference::Contig(0, "contig0", vectorFromString(reference)));

    isaac::reference::SortedReferenceMetadata sortedReferenceMetadata;
    sortedReferenceMetadata.putContig(0, contigList.front().name_, "blah.txt",
                                      0, contigList.front().getLength(),
                                      contigList.front().getLength(),
                                      contigList.front().getLength(), 0, 0, "",
                                      "", "");
    const isaac::reference::SortedReferenceMetadataList sortedReferenceMetadataList(1, sortedReferenceMetadata);

    isaac::common::ThreadVector threads(1);
    isaac::reference::ReferenceHasher<isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType> > referenceHasher(
        sortedReferenceMetadata, contigList, threads, threads.size());

    const isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType> referenceHash = referenceHasher.generate();

    isaac::flowcell::FlowcellLayoutList flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, isaac::flowcell::FastqFlowcellData(false, '!', false), 8, 0, std::vector<unsigned>(),
                                         readMetadataList, seedMetadataList, "blah"));

    isaac::flowcell::TileMetadata tileMetadata(
        flowcells.front().getFlowcellId(), flowcells.front().getIndex(), 0, 1, 1, 0);
    const isaac::flowcell::TileMetadataList tileMetadataList(std::vector<isaac::flowcell::TileMetadata>(1, tileMetadata));

    const unsigned clusterLength = isaac::flowcell::getTotalReadLength(readMetadataList);
    isaac::alignment::BclClusters tileClusters(clusterLength);

    tileClusters.reset(clusterLength, 1);

//    const std::vector<char>& bcls = getBcl(readMetadataList, contigList, 0, 0, readMetadataList.front().getLength());
    const std::vector<char>& bcls = getBcl(sequence.substr(0, clusterLength));

    std::copy(bcls.begin(), bcls.end(), tileClusters.cluster(0));

    isaac::alignment::Cluster cluster(readMetadataList.front().getLength());
    cluster.init(flowcells.front().getReadMetadataList(), tileClusters.cluster(0),
                 0, 0, isaac::alignment::ClusterXy(0,0), true, 0, 0);


    TestMatchStorage matches;
//    isaac::alignment::ClusterHashMatchFinder<isaac::oligo::VeryShortKmerType> matchFinder(
//        referenceHash, flowcells, isaac::flowcell::BarcodeMetadataList(), 0, repeatThreshold, std::vector<std::size_t>(),
//        sortedReferenceMetadataList, seedMetadataList, seedMetadataList.size());
    isaac::alignment::ClusterHashMatchFinder< isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType> > matchFinder(
        referenceHash, 0, repeatThreshold, seedMetadataList);

    isaac::alignment::matchFinder::TileClusterInfo tileClusterInfo(tileMetadataList);
    tileClusterInfo.setBarcodeIndex(0, 0, 0);
    matchFinder.findClusterMatches(cluster, flowcells.front().getReadMetadataList(), matches);
    return matches;
}

bool isProperMatch(const isaac::alignment::Match &match)
{
    return !match.isTooManyMatch() && !match.isNoMatch();
}

void TestHashMatchFinder::testEverything()
{
    ISAAC_SCOPE_BLOCK_CERR
    {
    {
        std::string reference("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(reference, sequence, 10, seedMetadataList_, readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(2UL, seedMetadataList_.size());
        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
    }

    {
        std::string reference("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAA"
                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(reference, sequence, 10, seedMetadataList_, readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(2UL, seedMetadataList_.size());
        CPPUNIT_ASSERT_EQUAL(3UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
    }

    {
        std::string reference("TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(reference, sequence, 10, seedMetadataList_, readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(2UL, seedMetadataList_.size());
        CPPUNIT_ASSERT_EQUAL(3UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
    }

    {
        std::string reference("TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(
            reference, sequence, 10,
            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
    }

    {
        std::string reference(//"TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
                              "ATTAAAAAAATAAAGA"
                              "ATTAAAAAAATAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(
            reference, sequence, 2,
            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(1UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
    }

    {
        std::string reference("TTTTTTGTATTTTCTTCTTGTTA"
            // seed extension disambiguates between this one:;;;
                              "TCTTTACATTTTTAAT"
            // and this one:
                              "ATTAAAAATGTAAAGT"
                              "TAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAATGTAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(
            reference, sequence, 2,
            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(1UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 23, false, false), matches.front().location_);
        CPPUNIT_ASSERT_EQUAL(isaac::alignment::SeedId(16, true), matches.front().seedId_);
    }

    //seed at the end does not get extended
//    {
//        std::string reference("TCCCAACCAATAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 3,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( sequence.length() - 8, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(),
//            boost::assign::list_of(isaac::flowcell::ReadMetadata( 1, sequence.length() + 1, 0, 0)).convert_to_container<isaac::flowcell::ReadMetadataList>());
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
//    }

    {
        std::string reference("GGGGTTTTCCCAACCTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(
            reference, sequence, 2,
            boost::assign::list_of(isaac::alignment::SeedMetadata( sequence.length() - 8, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(),
            boost::assign::list_of(isaac::flowcell::ReadMetadata( 1, sequence.length() + 1, 0, 0)).convert_to_container<isaac::flowcell::ReadMetadataList>());

        CPPUNIT_ASSERT_EQUAL(0UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
    }

    {
        std::string reference("ATTAAAAAATTAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTATTAAAAA");
        std::string sequence ("ATTAAAAAATTAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        TestMatchStorage matches = findMatches(
            reference, sequence, 2,
            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);

        CPPUNIT_ASSERT_EQUAL(1UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
    }
    }
}
