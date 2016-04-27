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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testMatchFinderClusterInfo.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestMatchFinderClusterInfo, registryName("MatchFinderClusterInfo"));

void TestMatchFinderClusterInfo::setUp()
{
}

void TestMatchFinderClusterInfo::tearDown()
{
}


void TestMatchFinderClusterInfo::testFields()
{
    using isaac::alignment::matchFinder::ClusterInfo;
    ClusterInfo none;
    ClusterInfo all;
    ClusterInfo other;
    other.setBarcodeIndex(0U);
    ClusterInfo barcode;
    ClusterInfo r1;
    ClusterInfo r2;
    // barcode
    CPPUNIT_ASSERT(!barcode.isBarcodeSet());
    // r1
    CPPUNIT_ASSERT(!r1.isBarcodeSet());
    // r2
    CPPUNIT_ASSERT(!r2.isBarcodeSet());
    // none
    CPPUNIT_ASSERT(!none.isBarcodeSet());
    // all
    CPPUNIT_ASSERT(!all.isBarcodeSet());
    // OTHER
    CPPUNIT_ASSERT_EQUAL(0U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());

    other.setBarcodeIndex(1234U);
    CPPUNIT_ASSERT_EQUAL(1234U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());

    other.setBarcodeIndex(482U);
    CPPUNIT_ASSERT_EQUAL(482U, other.getBarcodeIndex());
    CPPUNIT_ASSERT(other.isBarcodeSet());
}

