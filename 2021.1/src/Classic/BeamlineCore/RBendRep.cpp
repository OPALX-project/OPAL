//
// Class RBendRep
//   Representation for a rectangular bend magnet.
//   A rectangular bend magnet has a rectilinear geometry about which its
//   multipole components are specified.
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
#include "BeamlineCore/RBendRep.h"
#include "Channels/IndexedChannel.h"
#include "Channels/IndirectChannel.h"
#include <cctype>


RBendRep::RBendRep():
    RBend(),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


RBendRep::RBendRep(const RBendRep &rhs):
    RBend(rhs),
    geometry(rhs.geometry),
    field(rhs.field) {
    rEntry = rhs.rEntry;
    rExit  = rhs.hExit;
    hEntry = rhs.rEntry;
    hExit  = rhs.hExit;
}


RBendRep::RBendRep(const std::string &name):
    RBend(name),
    geometry(0.0, 0.0),
    field() {
    rEntry = rExit = hEntry = hExit = 0.0;
}


RBendRep::~RBendRep()
{}


ElementBase *RBendRep::clone() const {
    return new RBendRep(*this);
}


Channel *RBendRep::getChannel(const std::string &aKey, bool create) {
    return ElementBase::getChannel(aKey, create);
}


BMultipoleField &RBendRep::getField() {
    return field;
}

const BMultipoleField &RBendRep::getField() const {
    return field;
}


RBendGeometry &RBendRep::getGeometry() {
    return geometry;
}

const RBendGeometry &RBendRep::getGeometry() const {
    return geometry;
}


double RBendRep::getB() const {
    return field.getNormalComponent(1);
}

void RBendRep::setB(double B) {
    field.setNormalComponent(1, B);
}


double RBendRep::getEntryFaceRotation() const {
    return rEntry;
}


double RBendRep::getExitFaceRotation() const {
    return rExit;
}

void RBendRep::setEntryFaceRotation(double e1) {
    rEntry = e1;
}

void RBendRep::setExitFaceRotation(double e2) {
    rExit = e2;
}

double RBendRep::getEntryFaceCurvature() const {
    return hEntry;
}

double RBendRep::getExitFaceCurvature() const {
    return hExit;
}

void RBendRep::setEntryFaceCurvature(double h1) {
    hEntry = h1;
}

void RBendRep::setExitFaceCurvature(double h2) {
    hExit = h2;
}


double RBendRep::getSlices() const {
    return slices;
}

double RBendRep::getStepsize() const {
    return stepsize;
}

void RBendRep::setSlices(double sl) {
    slices = sl;
}

void RBendRep::setStepsize(double ds) {
    stepsize = ds;
}


void RBendRep::setField(const BMultipoleField &f) {
    field = f;
}