/**
 ** Isaac Genome Alignment Software
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
 ** \file KUniqueness.hh
 **
 ** \brief Basic declarations for computing k-uniqueness
 **
 ** \author Roman Petovski
 **/

#ifndef ISAAC_REFERENCE_K_UNIQUENESS_HH
#define ISAAC_REFERENCE_K_UNIQUENESS_HH

#include "common/NumaContainer.hh"
#include "common/SameAllocatorVector.hh"

namespace isaac
{
namespace reference
{

typedef unsigned char AnnotationValue;
//typedef std::vector<std::pair<AnnotationValue, AnnotationValue>, common::NumaAllocator<std::pair<AnnotationValue, AnnotationValue>, 0> > ContigAnnotation;
//typedef std::vector<ContigAnnotation, common::NumaAllocator<ContigAnnotation, 0> > ContigAnnotations;
//typedef std::vector<ContigAnnotations, common::NumaAllocator<ContigAnnotations, 0> > ContigAnnotationsList;

typedef common::SameAllocatorVectorEnd<std::pair<AnnotationValue, AnnotationValue>, common::NumaAllocator<std::pair<AnnotationValue, AnnotationValue>, 0> > ContigAnnotation;
typedef common::SameAllocatorVector<ContigAnnotation, common::NumaAllocator<ContigAnnotation, 0> > ContigAnnotations;
typedef common::SameAllocatorVector<ContigAnnotations, common::NumaAllocator<ContigAnnotations, 0> > ContigAnnotationsList;


typedef AnnotationValue DistanceToBeNeighborless;
static const DistanceToBeNeighborless K_UNIQUE_TOO_FAR = DistanceToBeNeighborless(0) - 1;



class NumaContigAnnotationsList
{
    common::NumaContainerReplicas<ContigAnnotationsList> replicas_;
public:

    NumaContigAnnotationsList(ContigAnnotationsList &&node0Lists) :replicas_(std::move(node0Lists))
    {
        ISAAC_THREAD_CERR << "NumaContigAnnotationsList copy constructor" << std::endl;
    }

    const ContigAnnotationsList &node0Container() const {return replicas_.node0Container();}
    const ContigAnnotationsList &threadNodeContainer() const {return replicas_.threadNodeContainer();}
//    operator const ContigAnnotationsList &()const {return replicas_.threadNodeContainer();}
};

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_K_UNIQUENESS_HH
