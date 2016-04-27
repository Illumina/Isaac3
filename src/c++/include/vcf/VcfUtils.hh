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
 ** \file VcfUtils.hh
 **
 ** Various helpers for dealign with simple vcf files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_VCF_VCF_UTILS_HH
#define iSAAC_VCF_VCF_UTILS_HH

#include <utility>

#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "common/StaticVector.hh"

namespace isaac
{
namespace vcf
{
/**
 ** \brief Exception thrown when an invalid command line option was detected.
 **
 **/
class VcfError: public std::logic_error, public common::ExceptionData
{
public:
    VcfError(const std::string &message) :
        std::logic_error(message),
        ExceptionData(EINVAL, message)
    {
    }

};


// pairs of iterators in the container pointing to start and end of parsed data elements.
template <typename IteratorT>
struct ParsedVcfLine
{
    static const int COMPONENT_CHROM = 0;
    static const int COMPONENT_POS = 1;
    static const int COMPONENT_ID = 2;
    static const int COMPONENT_REF = 3;
    static const int COMPONENT_ALT = 4;
    static const int COMPONENT_QUAL = 5;

    typedef std::pair<IteratorT, IteratorT> Component;
    common::StaticVector<Component, 8> components_;

    void clear() {components_.clear();}
    bool isEmpty() const {return components_.empty() || (1 == components_.size() && components_.front().first == components_.front().second);}
    bool isComment() const {return '#' == *components_.front().first;}
    bool isValid() const {return isEmpty() || isComment() || 6 <= components_.size();}
    bool isDataLine() const {return isValid() && !isEmpty() && !isComment();}
    std::string getChromosome() const {return std::string(components_[COMPONENT_CHROM].first, components_[COMPONENT_CHROM].second);}
    uint64_t getPosition() const {return common::getUnsignedInteger(components_[COMPONENT_POS].first, components_[COMPONENT_POS].second);}
    unsigned refLength() const {return std::distance(components_[COMPONENT_REF].first, components_[COMPONENT_REF].second);}
    unsigned altLength() const {return std::distance(components_[COMPONENT_ALT].first, components_[COMPONENT_ALT].second);}
    bool isMultiAlternative() const {return components_[COMPONENT_ALT].second !=  std::find(components_[COMPONENT_ALT].first, components_[COMPONENT_ALT].second, ',');}

    class AlternativesConstIterator
    {
        const ParsedVcfLine &line_;
        IteratorT current_;
    public:
        AlternativesConstIterator(const ParsedVcfLine &line, IteratorT it):
            line_(line), current_(it){}
        bool operator == (const AlternativesConstIterator &that) const {return &line_ == &that.line_ && current_ == that.current_;}
        bool operator != (const AlternativesConstIterator &that) const {return !(*this == that);}

        AlternativesConstIterator &operator ++ ()
        {
            current_ = std::find(current_, line_.components_[COMPONENT_ALT].second, ',');
            if (line_.components_[COMPONENT_ALT].second != current_)
            {
                ++current_;
            }
            return *this;
        }

        const ParsedVcfLine operator *() const
        {
            ParsedVcfLine ret(line_);
            ret.components_[COMPONENT_ALT].first = current_;
            ret.components_[COMPONENT_ALT].second = std::find(current_, line_.components_[COMPONENT_ALT].second, ',');
            return ret;
        }
    };
    AlternativesConstIterator alternativesBegin() const {return AlternativesConstIterator(*this, components_[COMPONENT_ALT].first);}
    AlternativesConstIterator alternativesEnd() const {return AlternativesConstIterator(*this, components_[COMPONENT_ALT].second);}

    friend std::ostream &operator <<(std::ostream &os, const ParsedVcfLine &line)
    {
        os << "ParsedVcfLine(";
        BOOST_FOREACH(const Component& component, line.components_)
        {
            os << "[" << common::makeFastIoString(component.first, component.second) << "]";
        }
        return os << ")";
    }
};

template <typename IteratorT>
IteratorT parseVcfLine(IteratorT begin, IteratorT end, ParsedVcfLine<IteratorT> &parsedLine)
{
    parsedLine.clear();

    while (end != begin)
    {
        static const char delimiters[] = {'\t', '\n', '\r'};
        static const char* delimitersEnd = delimiters + sizeof(delimiters);
        static const char lineTerminators[] = {'\n', '\r'};
        static const char* lineTerminatorsEnd = lineTerminators + sizeof(lineTerminators);

        const typename ParsedVcfLine<IteratorT>::Component component =
            std::make_pair(begin, std::find_first_of(begin, end, delimiters, delimitersEnd));

        if (parsedLine.components_.capacity() != parsedLine.components_.size())
        {
            parsedLine.components_.push_back(component);
        }
        begin = component.second;
        if (end != begin)
        {
            if (lineTerminatorsEnd != std::find(lineTerminators, lineTerminatorsEnd, *begin))
            {
                begin = std::find_if(begin, end, [](char c){return lineTerminatorsEnd == std::find(lineTerminators, lineTerminatorsEnd, c);});
            }
            else
            {
                ++begin;
            }
        }
    }

    return begin;
}


} // namespace vcf
} // namespace isaac

#endif // #ifndef iSAAC_VCF_VCF_UTILS_HH
