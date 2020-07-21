//
// Class DriftRep
//   Representation for a drift space.
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
#include "BeamlineCore/DriftRep.h"
#include "Channels/IndirectChannel.h"


namespace {
    struct Entry {
        const char *name;
        double(DriftRep::*get)() const;
        void (DriftRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &DriftRep::getElementLength,
            &DriftRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


DriftRep::DriftRep():
    Drift(),
    geometry(0.0)
{}


DriftRep::DriftRep(const DriftRep &right):
    Drift(right),
    geometry(right.geometry)
{}


DriftRep::DriftRep(const std::string &name):
    Drift(name),
    geometry()
{}


DriftRep::~DriftRep()
{}


ElementBase *DriftRep::clone() const {
    return new DriftRep(*this);
}


Channel *DriftRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<DriftRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &DriftRep::getField() {
    return field;
}

const NullField &DriftRep::getField() const {
    return field;
}


StraightGeometry &DriftRep::getGeometry() {
    return geometry;
}

const StraightGeometry &DriftRep::getGeometry() const {
    return geometry;
}