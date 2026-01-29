//
// Class CyclotronRep
//   Representation for a cyclotron magnet system
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
#include "BeamlineCore/CyclotronRep.h"
#include "Channels/IndirectChannel.h"

// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CyclotronRep::*get)() const;
        void (CyclotronRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "RINIT",
            0, // :FIXME: Why commented out? &CyclotronRep::getRadius, 
            0  // :FIXME: &CyclotronRep::setRadius
        },
        { 0, 0, 0 }
    };
}


CyclotronRep::CyclotronRep():
    Cyclotron(),
    geometry(0.0, 0.0),
    field() {
    rInit = 0.0;
}


CyclotronRep::CyclotronRep(const CyclotronRep &rhs):
    Cyclotron(rhs),
    geometry(rhs.geometry),
    field(rhs.field) {
    rInit = rhs.rInit;
}


CyclotronRep::CyclotronRep(const std::string &name):
    Cyclotron(name),
    geometry(0.0, 0.0),
    field() {
    rInit = 0.0;
}


CyclotronRep::~CyclotronRep()
{}


ElementBase *CyclotronRep::clone() const {
    return new CyclotronRep(*this);
}


Channel *CyclotronRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
      if(aKey == entry->name) {
        return new IndirectChannel<CyclotronRep>(*this, entry->get, entry->set);
      }
    }

    return ElementBase::getChannel(aKey, create);
}


double CyclotronRep::getSlices() const {
    return slices;
}

double CyclotronRep::getStepsize() const {
    return stepsize;
}

void CyclotronRep::setSlices(double sl) {
    slices = sl;
}

void CyclotronRep::setStepsize(double ds) {
    stepsize = ds;
}

/*
void CyclotronRep::setRadius(double r)
{
  rInit = r;
}

double CyclotronRep::getRadius() const
{
  return rInit ;
}
*/

PlanarArcGeometry &CyclotronRep::getGeometry() {
    return geometry;
}

const PlanarArcGeometry &CyclotronRep::getGeometry() const {
    return geometry;
}

BMultipoleField &CyclotronRep::getField() {
    return field;
}

const BMultipoleField &CyclotronRep::getField() const {
    return field;
}

void CyclotronRep::setField(const BMultipoleField &f) {
    field = f;
}

