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

#include <vector>
#include <string>
#include <boost/foreach.hpp>

#include "alignment/BclClusters.hh"
#include "alignment/SeedMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"
#include "reference/Contig.hh"

inline isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned l0 = 100, const unsigned l1 = 100)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

inline isaac::alignment::SeedMetadataList getSeedMetadataList()
{
    std::vector<isaac::alignment::SeedMetadata> ret =
        boost::assign::list_of
        (isaac::alignment::SeedMetadata( 0, 32, 0, 0))
        (isaac::alignment::SeedMetadata(32, 32, 0, 1))
        (isaac::alignment::SeedMetadata(64, 32, 0, 2))
        (isaac::alignment::SeedMetadata( 0, 32, 1, 3))
        (isaac::alignment::SeedMetadata(32, 32, 1, 4))
        (isaac::alignment::SeedMetadata(64, 32, 1, 5))
        ;
    return ret;
}

inline isaac::reference::Contig getContig(const std::string name, const unsigned length)
{
    char bases[] = {'A', 'C', 'G', 'T'};
    isaac::reference::Contig contig(0, name);
    contig.resize(length);
    BOOST_FOREACH(char &base, contig)
    {
        base = bases[rand() % 4];
    }
    return contig;
}

inline void show(const std::vector<char> &s)
{
    BOOST_FOREACH(char c, s) std::cerr << c;
}

template <typename ContainerT>
inline std::vector<char> reverseComplement(const ContainerT &forward)
{
    std::vector<char> tr(256, 'N');
    tr['A'] = 'T';
    tr['C'] = 'G';
    tr['G'] = 'C';
    tr['T'] = 'A';
    std::vector<char> reverse;
    reverse.reserve(forward.size());
    for (typename ContainerT::const_reverse_iterator b = forward.rbegin(); forward.rend() != b; ++b)
    {
        reverse.push_back(tr[*b]);
    }
    return reverse;
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

inline std::vector<char> vectorFromString(const std::string &str, bool ignoreSpaces)
{
    if (ignoreSpaces)
    {
        std::vector<char> ret;
        std::remove_copy(str.begin(), str.end(), std::back_inserter(ret), ' ');
        return ret;
    }
    return std::vector<char>(str.begin(), str.end());
}

template <typename ContainerT>
inline std::string substr(const ContainerT &from, std::string::size_type __pos = 0,
                   std::string::size_type __n = std::string::npos)
{
    return std::string(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

template <typename ContainerT>
inline std::vector<char> subv(const ContainerT &from, std::string::size_type __pos = 0,
                       std::string::size_type __n = std::string::npos)
{
    return std::vector<char>(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

//inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos,
//                       std::string::size_type __n)
//{
//    return std::vector<char>(from.begin() + __pos, from.begin() + __pos + __n);
//}

//template <typename ContainerT>
//inline std::vector<char> subv(const ContainerT &from, std::string::size_type __pos)
//{
//    return std::vector<char>(from.begin() + __pos, from.end());
//}

inline std::vector<char> operator +(const std::vector<char> &right, const isaac::reference::Contig &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}


inline std::vector<char> operator +(const std::vector<char> &right, const std::vector<char> &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline std::vector<char> operator +(const std::vector<char> &right, const std::string &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline isaac::reference::ContigList getContigList(const unsigned l0 = 210, const unsigned l1 = 220, const unsigned l4 = 60)
{
    const isaac::reference::Contig c2 = getContig("c2", 230);
    std::vector<char> v = vectorFromString(std::string("AAAAA"));
    v.insert(v.end(), c2.begin(), c2.end());
    isaac::reference::Contig c3(0, "c3", v);
    const isaac::reference::Contig c4 = getContig("c4", l4);
    const isaac::reference::Contig c1 = getContig("c1", l1);
    const isaac::reference::Contig c0 = getContig("c0", l0);

    return boost::assign::list_of(c0)(c1)(c2)(c3)(c4).convert_to_container<isaac::reference::ContigList>();
}

template <typename ContainerT>
std::vector<char> getBcl(const ContainerT &bases)
{
    std::vector<char> bcl;
    bcl.reserve(bases.size());
    BOOST_FOREACH(char b, bases)
    {
        using isaac::oligo::getValue;
        bcl.push_back((40 << 2) | getValue(b));
    }
    return bcl;
}

inline std::vector<char> getBclVector(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const isaac::reference::ContigList &contigList,
    const unsigned contigId, const int offset0, const int offset1,
    const bool reverse0 = false, const bool reverse1 = true)
{
    const isaac::reference::Contig &contig = contigList[contigId];
    const unsigned length0 = readMetadataList[0].getLength();
    const unsigned length1 = readMetadataList[1].getLength();

    isaac::reference::Contig reverse;
    reverse.reserve(contig.size());
    BOOST_REVERSE_FOREACH(const char base, contig)
    {
        using isaac::oligo::getReverseBase;
        using isaac::oligo::getValue;
        reverse.push_back(getReverseBase(getValue(base)));
    }

    const isaac::reference::Contig &s0 = reverse0 ? reverse : contig;
    const isaac::reference::Contig &s1 = reverse1 ? reverse : contig;
    std::string bases(s0.begin() + offset0, s0.begin() + offset0 + length0);
    bases += std::string(s1.begin() + offset1, s1.begin() + offset1 + length1);
    return getBcl(bases);
}

inline isaac::alignment::BclClusters getBclClusters(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const std::vector<char> &clusters)
{
    isaac::alignment::BclClusters ret(isaac::flowcell::getTotalReadLength(readMetadataList));
    ret.reserveClusters(1, false);
    std::copy(clusters.begin(), clusters.end(), ret.cluster(0));
    return ret;
}

inline isaac::alignment::BclClusters getBcl(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const isaac::reference::ContigList &contigList,
    const unsigned contigId, const int offset0, const int offset1,
    const bool reverse0 = false, const bool reverse1 = true)
{
    const std::vector<char> clusters = getBclVector(
        readMetadataList, contigList, contigId, offset0, offset1, reverse0, reverse1);
    return getBclClusters(readMetadataList, clusters);
}
