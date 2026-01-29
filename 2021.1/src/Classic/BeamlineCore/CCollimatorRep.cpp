//
// Class CCollimatorRep
//   Representation for a collimator.
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
#include "BeamlineCore/CCollimatorRep.h"
#include "Channels/IndirectChannel.h"


namespace {
    struct Entry {
        const char *name;
        double(CCollimatorRep::*get)() const;
        void (CCollimatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &CCollimatorRep::getElementLength,
            &CCollimatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


CCollimatorRep::CCollimatorRep():
    CCollimator(),
    geometry(0.0)
{}


CCollimatorRep::CCollimatorRep(const CCollimatorRep &right):
    CCollimator(right),
    geometry(right.geometry)
{}


CCollimatorRep::CCollimatorRep(const std::string &name):
    CCollimator(name),
    geometry()
{}


CCollimatorRep::~CCollimatorRep()
{}


ElementBase *CCollimatorRep::clone() const {
    return new CCollimatorRep(*this);
}


Channel *CCollimatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<CCollimatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &CCollimatorRep::getField() {
    return field;
}

const NullField &CCollimatorRep::getField() const {
    return field;
}


StraightGeometry &CCollimatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &CCollimatorRep::getGeometry() const {
    return geometry;
}



/*
double CCollimatorRep::getXsize() const
{
  return xSize;
}

double CCollimatorRep::getYsize() const
{
  return ySize;
}

void CCollimatorRep::setXsize(double size)
{
  INFOMSG("void CCollimatorRep::setXsize(double size) " << xSize << endl;);
  xSize = size;
}

void CCollimatorRep::setYsize(double size)
{
  ySize = size;
}
*/