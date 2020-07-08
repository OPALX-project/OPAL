//
// Class BeamlineVisitor
//   The abstract class BeamlineVisitor is the base class for all visitors
//   (algorithms) that can iterator over a beam line representation.
//   A BeamlineVisitor applies itself to the representation via the
//   ``Visitor'' pattern, see
//   [p]
//   E. Gamma, R. Helm, R. Johnson, and J. Vlissides,
//   [BR]
//   Design Patterns, Elements of Reusable Object-Oriented Software.
//   [p]
//   By using only pure abstract classes as an interface between the
//   BeamlineVisitor and the beam line representation,
//   we decouple the former from the implementation details of the latter.
//   [p]
//   The interface is defined in such a way that a visitor cannot modify the
//   structure of a beam line, but it can assign special data like misalignments
//   or integrators without problems.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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

#include "AbsBeamline/BeamlineVisitor.h"

BeamlineVisitor::BeamlineVisitor()
{}


BeamlineVisitor::~BeamlineVisitor()
{}
