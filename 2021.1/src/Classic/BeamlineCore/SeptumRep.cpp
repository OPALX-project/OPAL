//
// Class SeptumRep
//   Representation for Septum.
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
#include "BeamlineCore/SeptumRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(SeptumRep::*get)() const;
        void (SeptumRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &SeptumRep::getElementLength,
            &SeptumRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


SeptumRep::SeptumRep():
    Septum(), field(), geometry(), active(true)
{}


SeptumRep::SeptumRep(const SeptumRep &right):
    Septum(right), field(), geometry(right.geometry), active(true)
{}


SeptumRep::SeptumRep(const std::string &name):
    Septum(name), field(), geometry(), active(true)
{}


SeptumRep::~SeptumRep()
{}


ElementBase *SeptumRep::clone() const {
    return new SeptumRep(*this);
}


Channel *SeptumRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SeptumRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &SeptumRep::getField() {
    return field;
}

const NullField &SeptumRep::getField() const {
    return field;
}

StraightGeometry &SeptumRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SeptumRep::getGeometry() const {
    return geometry;
}


void SeptumRep::setActive(bool flag) {
    active = flag;
}
