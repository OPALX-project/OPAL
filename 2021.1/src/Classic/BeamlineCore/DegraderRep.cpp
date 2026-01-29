//
// Class DegraderRep
//   Representation for a degrader.
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
#include "BeamlineCore/DegraderRep.h"
#include "Channels/IndirectChannel.h"


namespace {
    struct Entry {
        const char *name;
        double(DegraderRep::*get)() const;
        void (DegraderRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &DegraderRep::getElementLength,
            &DegraderRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


DegraderRep::DegraderRep():
    Degrader(),
    geometry(0.0)
{}


DegraderRep::DegraderRep(const DegraderRep &right):
    Degrader(right),
    geometry(right.geometry)
{}


DegraderRep::DegraderRep(const std::string &name):
    Degrader(name),
    geometry()
{}


DegraderRep::~DegraderRep()
{}


ElementBase *DegraderRep::clone() const {
    return new DegraderRep(*this);
}


Channel *DegraderRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<DegraderRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &DegraderRep::getField() {
    return field;
}

const NullField &DegraderRep::getField() const {
    return field;
}


StraightGeometry &DegraderRep::getGeometry() {
    return geometry;
}

const StraightGeometry &DegraderRep::getGeometry() const {
    return geometry;
}