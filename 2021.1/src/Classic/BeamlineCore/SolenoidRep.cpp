//
// Class SolenoidRep
//   Representation for a solenoid magnet.
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
#include "BeamlineCore/SolenoidRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(SolenoidRep::*get)() const;
        void (SolenoidRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &SolenoidRep::getElementLength,
            &SolenoidRep::setElementLength
        },
        {
            "BZ",
            &SolenoidRep::getBz,
            &SolenoidRep::setBz
        },
        { 0, 0, 0 }
    };
}


SolenoidRep::SolenoidRep():
    Solenoid(),
    geometry(),
    field()
{}


SolenoidRep::SolenoidRep(const SolenoidRep &right):
    Solenoid(right),
    geometry(right.geometry),
    field(right.field)
{}


SolenoidRep::SolenoidRep(const std::string &name):
    Solenoid(name), geometry(), field()
{}


SolenoidRep::~SolenoidRep()
{}


ElementBase *SolenoidRep::clone() const {
    return new SolenoidRep(*this);
}


Channel *SolenoidRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SolenoidRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


ConstBzField &SolenoidRep::getField() {
    return field;
}

const ConstBzField &SolenoidRep::getField() const {
    return field;
}


StraightGeometry &SolenoidRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SolenoidRep::getGeometry() const {
    return geometry;
}


double SolenoidRep::getBz() const {
    return field.getBz();
}


void SolenoidRep::setBz(double Bz) {
    field.setBz(Bz);
}
