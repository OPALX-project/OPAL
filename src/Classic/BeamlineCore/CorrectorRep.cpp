//
// Class CorrectorRep
//   Representation of a closed orbit corrector.
//   The base class acts on both planes.
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
#include "BeamlineCore/CorrectorRep.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(CorrectorRep::*get)() const;
        void (CorrectorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &CorrectorRep::getElementLength,
            &CorrectorRep::setElementLength
        },
        {
            "BY",
            &CorrectorRep::getBy,
            &CorrectorRep::setBy
        },
        {
            "BX",
            &CorrectorRep::getBx,
            &CorrectorRep::setBx
        },
        { 0, 0, 0 }
    };
}


CorrectorRep::CorrectorRep():
    Corrector(), geometry(), field(), active(true)
{}


CorrectorRep::CorrectorRep(const CorrectorRep &right):
    Corrector(right), geometry(right.geometry), field(right.field), active(true)
{}


CorrectorRep::CorrectorRep(const std::string &name):
    Corrector(name), geometry(), field(), active(true)
{}


CorrectorRep::~CorrectorRep()
{}


ElementBase *CorrectorRep::clone() const {
    return new CorrectorRep(*this);
}


Channel *CorrectorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *table = entries; table->name != 0; ++table) {
        if(aKey == table->name) {
            return new IndirectChannel<CorrectorRep>(*this, table->get, table->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


StraightGeometry &CorrectorRep::getGeometry() {
    return geometry;
}


const StraightGeometry &CorrectorRep::getGeometry() const {
    return geometry;
}


Corrector::Plane CorrectorRep::getPlane() const {
    return active ? XY : OFF;
}


double CorrectorRep::getBx() const {
    return field.getBx();
}


double CorrectorRep::getBy() const {
    return field.getBy();
}


BDipoleField &CorrectorRep::getField() {
    return field;
}


const BDipoleField &CorrectorRep::getField() const {
    return field;
}


void CorrectorRep::setBx(double Bx) {
    field.setBx(Bx);
}


void CorrectorRep::setBy(double By) {
    field.setBy(By);
}

void CorrectorRep::setActive(bool flag) {
    active = flag;
}
