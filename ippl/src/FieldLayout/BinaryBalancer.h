//
// A fairly simple load balancer inspired by Dan Quinlan's MLB.
//
//   It does recursive binary subdivision of a FieldLayout domain,
//   restricting the cuts to coordinate directions, so as to balance the
//   workload.  The "workload" is given by a Field of weights passed in.
//   It decides on the cut axis by cutting the longest axis of a brick,
//   and the location of that cut by balancing the weights on each side
//   of the cut.  The resulting distribution has one vnode per processor.
//
//   This is restricted to a processor number that is a power of two. 
//
//   It performs log(P) parallel reductions.
//
//   It does nothing fancy when deciding on the splits to try to make the
//   new paritioning close to the previous.  The same set of weights will
//   always give the same repartitioning, but similar sets of weights
//   could result in quite different partitionings.
//
//   There are two functions defined here:
//
//   NDIndex<Dim>
//   CalcBinaryRepartion(FieldLayout<Dim>&, BareField<double,Dim>&);
//
//   Given a FieldLayout and a Field of weights, find the domain for this
//   processor.  This does not repartition the FieldLayout, it just
//   calculates the domain.  If you want to further subdivide these
//   domains, just cut up what this function returns.
//
//   void
//   BinaryRepartition(FieldLayout<Dim>&, BareField<double,Dim>&);
//
//   Just call the above function and then repartition the FieldLayout
//   (and all the Fields defined on it).
//
// Copyright (c) 2003 - 2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef BINARY_BALANCER_H
#define BINARY_BALANCER_H

// forward declarations
template<unsigned Dim> class FieldLayout;
template<class T, unsigned Dim> class BareField;

class BinaryRepartitionFailed {  };

// Calculate the local domain for a binary repartition.
template<unsigned Dim>
NDIndex<Dim>
CalcBinaryRepartition(FieldLayout<Dim>&, BareField<double,Dim>&);

// Calculate and apply a local domain for a binary repartition.
template<unsigned Dim>
inline void
BinaryRepartition(FieldLayout<Dim>& layout, BareField<double,Dim>& weights)
{
  layout.Repartition( CalcBinaryRepartition(layout,weights) );
}

#include "FieldLayout/BinaryBalancer.hpp"

#endif
