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
 ** \file BandedSmithWaterman.cpp
 **
 ** \brief See BandedSmithWaterman.hh
 **
 ** \author Come Raczy
 **/

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <boost/format.hpp>
#include <cstdint>

#include "alignment/BandedSmithWaterman.hh"

namespace isaac
{
namespace alignment
{

BandedSmithWaterman::BandedSmithWaterman(const int matchScore, const int mismatchScore,
                                         const int gapOpenScore, const int gapExtendScore,
                                         const int maxReadLength)
    : matchScore_(matchScore)
    , mismatchScore_(mismatchScore)
    , gapOpenScore_(gapOpenScore)
    , gapExtendScore_(gapExtendScore)
    , maxReadLength_(maxReadLength)
    , initialValue_(static_cast<int>(std::numeric_limits<short>::min()) + gapOpenScore_)
    , T_(new char[maxReadLength_ * 3 * WIDEST_GAP_SIZE * sizeof(int16_t)])
{
    // check that there won't be any overflows in the matrices
    const int maxScore = std::max(std::max(std::max(abs(matchScore_), abs(mismatchScore_)), abs(gapOpenScore_)), abs(gapExtendScore_));
    if ((maxReadLength_ * maxScore) >= abs(static_cast<int>(initialValue_)))
    {
        const std::string message = (boost::format("BandedSmithWaterman: unsupported read length (%i) for these scores (%i): use smaller scores or shorter reads") % maxReadLength_ % maxScore).str();
        BOOST_THROW_EXCEPTION(isaac::common::InvalidParameterException(message));
    }
}

BandedSmithWaterman::~BandedSmithWaterman()
{
    free(T_);
}

void BandedSmithWaterman::cp(int16_t source[WIDEST_GAP_SIZE], int16_t destination[WIDEST_GAP_SIZE]) const
{
    // AV
    for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
    destination[i] = source[i];
    }
}

unsigned BandedSmithWaterman::align(
    const std::vector<char> &query,
    const reference::Contig::const_iterator databaseBegin,
    const reference::Contig::const_iterator databaseEnd,
    Cigar &cigar) const
{
    return align(query.begin(), query.end(), databaseBegin, databaseEnd, cigar);
}

unsigned BandedSmithWaterman::trimTailIndels(Cigar& cigar, const size_t beginOffset) const
{
    unsigned ret = 0;
    unsigned long extend = 0;
    for (Cigar::Component component = Cigar::decode(cigar.back());
        cigar.size() != beginOffset; component = Cigar::decode(cigar.back()))
    {
        if (Cigar::DELETE == component.second)
        {
            //CASAVA does not like CIGAR beginning with a deletion in the data
            cigar.pop_back();
            ISAAC_ASSERT_MSG(
                Cigar::DELETE != Cigar::decode(cigar.back()).second,
                "two Cigar::DELETE cannot be next to each other");
            ret += component.first;
        }
        else if (Cigar::INSERT == component.second)
        {
            // tail and head insertions are biasing best alignment choice. remove them and extend the adjacent align operation
            extend += component.first;
            ret -= component.first;
            cigar.pop_back();
        }
        else
        {
            break;
        }
    }

    if (extend)
    {
        if(cigar.size() == beginOffset)
        {
            // this was a really bad s-w alignment. Something like this: 15D7I7D15I15D15I15D15I15D15I15D
            cigar.addOperation(extend, Cigar::ALIGN);
        }
        else
        {
            const Cigar::Component component = Cigar::decode(cigar.back());
            ISAAC_ASSERT_MSG(Cigar::ALIGN == component.second, "Unexpected operation at the end of cigar " <<
                             Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << " : "  << component.second);
            cigar.updateOperation(cigar.size() - 1, component.first + extend, Cigar::ALIGN);
        }
    }
    return ret;
}

void BandedSmithWaterman::removeAdjacentIndels(Cigar& cigar, const size_t beginOffset) const
{
//    ISAAC_THREAD_CERR << "BandedSmithWaterman::removeAdjacentIndels " << Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << std::endl;
    const Cigar::iterator begin = cigar.begin() + beginOffset;
    for (Cigar::iterator it = begin + 1; cigar.end() != it;)
    {
        const Cigar::iterator last = it - 1;
        const Cigar::Component lastcomponent = Cigar::decode(*last);
        const Cigar::Component component = Cigar::decode(*it);
        if (Cigar::ALIGN == lastcomponent.second)
        {
            if (Cigar::ALIGN == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::ALIGN);
                it = cigar.erase(it);
            }
            else
            {
                ++it;
            }
        }
        else if (Cigar::DELETE == lastcomponent.second)
        {
            if (Cigar::DELETE == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::DELETE);
                it = cigar.erase(it);
            }
            else if (Cigar::INSERT == component.second)
            {
                if (component.first < lastcomponent.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first - component.first, Cigar::DELETE);
                    cigar.updateOperation(std::distance(cigar.begin(), it), component.first, Cigar::ALIGN);
                }
                else if (lastcomponent.first < component.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first - lastcomponent.first, Cigar::INSERT);
                    cigar.updateOperation(std::distance(cigar.begin(), it), lastcomponent.first, Cigar::ALIGN);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
                else
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first, Cigar::ALIGN);
                    it = cigar.erase(it);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
            }
            else
            {
                ++it;
            }
        }
        else if (Cigar::INSERT == lastcomponent.second)
        {
            if (Cigar::INSERT == component.second)
            {
                cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first + component.first, Cigar::INSERT);
                it = cigar.erase(it);
            }
            else if (Cigar::DELETE == component.second)
            {
                if (component.first < lastcomponent.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), lastcomponent.first - component.first, Cigar::INSERT);
                    cigar.updateOperation(std::distance(cigar.begin(), it), component.first, Cigar::ALIGN);
                }
                else if (lastcomponent.first < component.first)
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first - lastcomponent.first, Cigar::DELETE);
                    cigar.updateOperation(std::distance(cigar.begin(), it), lastcomponent.first, Cigar::ALIGN);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
                else
                {
                    cigar.updateOperation(std::distance(cigar.begin(), last), component.first, Cigar::ALIGN);
                    it = cigar.erase(it);
                    if (begin != it - 1)
                    {
                        --it;
                    }
                }
            }
            else
            {
                ++it;
            }
        }
        else
        {
            ++it;
        }
    }
    //ISAAC_THREAD_CERR << "BandedSmithWaterman::removeAdjacentIndels done " << Cigar::toString(cigar.begin() + beginOffset, cigar.end()) << std::endl;
}


unsigned BandedSmithWaterman::align(
    const std::vector<char>::const_iterator queryBegin,
    const std::vector<char>::const_iterator queryEnd,
    const reference::Contig::const_iterator databaseBegin,
    const reference::Contig::const_iterator databaseEnd,
    Cigar &cigar) const
{
    assert(databaseEnd > databaseBegin);
    const size_t querySize = std::distance(queryBegin, queryEnd);
    ISAAC_ASSERT_MSG(querySize + WIDEST_GAP_SIZE - 1 == (unsigned long)(databaseEnd - databaseBegin), "q:" << std::string(queryBegin, queryEnd) << " db:" << std::string(databaseBegin, databaseEnd));
    assert(querySize <= size_t(maxReadLength_));
    const size_t originalCigarSize = cigar.size();
    int16_t *t = (int16_t*)T_;

    int16_t GapOpenScore[WIDEST_GAP_SIZE], GapExtendScore[WIDEST_GAP_SIZE];
    for(unsigned i = 0; i < WIDEST_GAP_SIZE; i++) {
        GapOpenScore[i] = gapOpenScore_;
        GapExtendScore[i] = gapExtendScore_;
    }
    // Initialize E, F and G
    int16_t D[WIDEST_GAP_SIZE], E[WIDEST_GAP_SIZE], F[WIDEST_GAP_SIZE], G[WIDEST_GAP_SIZE];
    for(unsigned i = 0; i < WIDEST_GAP_SIZE; i++) {
        E[i] = initialValue_;
        F[i] = 0;
        G[i] = initialValue_;
    }
    G[0] = 0;

    for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
        D[i] = *(databaseBegin + (WIDEST_GAP_SIZE - i - 2));
    }

    // iterate over all bases in the query
    int16_t F1[WIDEST_GAP_SIZE + 1];
    int16_t cmpgtEgMask1[WIDEST_GAP_SIZE + 1], maxEg1[WIDEST_GAP_SIZE + 1];
    F1[0] = initialValue_ + gapExtendScore_;
    maxEg1[0] = initialValue_ + gapOpenScore_;
    cmpgtEgMask1[0] = 0;
    std::vector<char>::const_iterator queryCurrent = queryBegin;
    for (unsigned queryOffset = 0; queryEnd != queryCurrent; ++queryOffset, ++queryCurrent)
    {
        int16_t TE[WIDEST_GAP_SIZE], TF[WIDEST_GAP_SIZE], TG[WIDEST_GAP_SIZE];
        int16_t D1[WIDEST_GAP_SIZE + 1];
        int16_t Q[WIDEST_GAP_SIZE];
        int16_t GA[WIDEST_GAP_SIZE];

        // get F[i-1, j] - extend
        int16_t cmpgtGfMask[WIDEST_GAP_SIZE];

        int16_t *cmpgtEgMaskOff = cmpgtEgMask1 + 1;
        int16_t *maxEgOff = maxEg1 + 1;
        // AV
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            cmpgtEgMaskOff[i] = E[i] > G[i] ? 1 : 0;
        }
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            maxEgOff[i] = G[i] > E[i] ? G[i] : E[i];
        }
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            cmpgtGfMask[i] = F[i] > maxEgOff[i] ? 2 : 0;
        }
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            GA[i] = maxEgOff[i] > F[i] ? maxEgOff[i] : F[i];
        }
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            TG[i] =
                   cmpgtEgMaskOff[i] >
                cmpgtGfMask[i] ? cmpgtEgMaskOff[i] : cmpgtGfMask[i];
        }

        cp(F, F1 + 1);
        int16_t GF1[WIDEST_GAP_SIZE], maxEgSubGapOpen1[WIDEST_GAP_SIZE],
            cmpgtGfMask1[WIDEST_GAP_SIZE];
        // AV
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            GF1[i] = F1[i] - GapExtendScore[i];
            maxEgSubGapOpen1[i] = maxEg1[i] - GapOpenScore[i];
            cmpgtGfMask1[i] = GF1[i] > maxEgSubGapOpen1[i] ? 2 : 0;
            TF[i] =
                cmpgtEgMask1[i] >
                cmpgtGfMask1[i] ? cmpgtEgMask1[i] : cmpgtGfMask1[i];
            F[i] =
                maxEgSubGapOpen1[i] > GF1[i] ? maxEgSubGapOpen1[i] : GF1[i];
        }

        // add the match/mismatch score
        // load the query base in all 8 values of the register
        for(unsigned i = 0; i < WIDEST_GAP_SIZE; i++) { 
            Q[i] = *queryCurrent;
        }

        // shift the database by 1 byte to the left and add the new base

        cp(D, D1+1);
        D1[0] = *(databaseBegin + queryOffset + (WIDEST_GAP_SIZE - 1));
        cp(D1, D);

        // compare query and database. 0xff if different (that also the sign bits)
        int16_t B[WIDEST_GAP_SIZE], Match[WIDEST_GAP_SIZE],
        Mismatch[WIDEST_GAP_SIZE], W[WIDEST_GAP_SIZE];

        // lea
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            B[i] = (Q[i] == D[i]) ? 0 : 0xFFFF;
            Match[i] = (~B[i]) & matchScore_;
            Mismatch[i] = B[i] & mismatchScore_;
            W[i] = Match[i] + Mismatch[i];
            G[i] = GA[i] + (W[i] | (B[i] & 0xFF00));
        }

        // E[i,j] = max(G[i, j-1] - open, E[i, j-1] - extend, F[i, j-1] - open)
         int16_t cmpgtFgMask2[WIDEST_GAP_SIZE + 1], maxFg2[WIDEST_GAP_SIZE + 1];
           int16_t *cmpgtFgMaskOff2 = cmpgtFgMask2 + 1;
           int16_t *maxFgOff2 = maxFg2 + 1;
           // AV
           for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
               cmpgtFgMask2[i] = F[i] > G[i] ? 2 : 0;
               maxFg2[i] = F[i] > G[i] ? F[i] : G[i];
               maxFg2[i] -= GapOpenScore[i];
           }

           E[WIDEST_GAP_SIZE - 1] = initialValue_;
        maxFg2[WIDEST_GAP_SIZE] = initialValue_;

        short e = initialValue_;
        short fg = initialValue_;
        for (size_t i = WIDEST_GAP_SIZE; i > 0; i--) {
            short max = fg;
            if (e > fg) {
                max = e;
            }
           E[i - 1] = max;
           fg = maxFg2[i - 1];
    
           e = max - gapExtendScore_;
        }

        // lea
        int16_t cmpgtFgSueFgMask2[WIDEST_GAP_SIZE];
        cmpgtFgMask2[WIDEST_GAP_SIZE] = initialValue_;
        for (size_t i = 0; i < WIDEST_GAP_SIZE; i++) {
            cmpgtFgSueFgMask2[i] = E[i] > maxFgOff2[i] ? 5 : 0;
            E[i] = E[i] > maxFgOff2[i] ? E[i] : maxFgOff2[i];
            TE[i] =
                (cmpgtFgSueFgMask2[i] >
                 cmpgtFgMaskOff2[i] ? cmpgtFgSueFgMask2[i] :
                 cmpgtFgMaskOff2[i]) & 3;
        }

        TF[0] = 0;

        cp(TG, t);
        cp(TE, t + WIDEST_GAP_SIZE);
        cp(TF, t + WIDEST_GAP_SIZE * 2);
        t += WIDEST_GAP_SIZE * 3;
    }
    // find the max of E, F and G at the end
    short max = G[15] - 1;

    int ii = querySize - 1;
    int jj = ii;
    unsigned maxType = 0;


    int16_t *TT[] = {G, E, F};
    for (unsigned j = WIDEST_GAP_SIZE; j > 0; j--)
    {
        for (unsigned type = 0; 3 > type; ++type)
        {
            const short value = TT[type][j - 1];
            if (value > max)
            {
                max = value;
                jj = j - 1;
                maxType = type;
            }
        }
    }

    const int jjIncrement[] = {0, 1, -1};
    const int iiIncrement[] = {-1, 0, -1};
    const Cigar::OpCode opCodes[] = {Cigar::ALIGN, Cigar::DELETE, Cigar::INSERT};
    unsigned opLength = 0;
    if (jj > 0)
    {
        cigar.addOperation(jj, Cigar::DELETE);
    }
    while(ii >= 0 && jj >= 0 && jj <= 15)
    {
        ++opLength;
        const unsigned nextMaxType = T_[(ii * 3 + maxType)
            * WIDEST_GAP_SIZE * sizeof(int16_t) + (jj * 2)];
        if (nextMaxType != maxType)
        {
            cigar.addOperation(opLength, opCodes[maxType]);
            opLength = 0;
        }
        ii += iiIncrement[maxType];
        jj += jjIncrement[maxType];
        maxType = nextMaxType;
    }
    assert(-1 == ii);
    if (1 != maxType && opLength)
    {
        cigar.addOperation(opLength, opCodes[maxType]);
        opLength = 0;
    }
    if (15 > jj)
    {
        cigar.addOperation(opLength + 15 - jj, Cigar::DELETE);
        opLength = 0;
    }
    assert(0 == opLength);
    unsigned ret = trimTailIndels(cigar, originalCigarSize);
    std::reverse(cigar.begin() + originalCigarSize, cigar.end());
    trimTailIndels(cigar, originalCigarSize);
    removeAdjacentIndels(cigar, originalCigarSize);
    return ret;
}

} // namespace alignment
} // namespace isaac
