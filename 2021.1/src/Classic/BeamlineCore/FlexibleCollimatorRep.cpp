//
// Class FlexibleCollimatorRep
//   Representation for a flexible collimator.
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
#include "BeamlineCore/FlexibleCollimatorRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(FlexibleCollimatorRep::*get)() const;
        void (FlexibleCollimatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &FlexibleCollimatorRep::getElementLength,
            &FlexibleCollimatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


FlexibleCollimatorRep::FlexibleCollimatorRep():
    FlexibleCollimator(),
    geometry(0.0)
{}


FlexibleCollimatorRep::FlexibleCollimatorRep(const FlexibleCollimatorRep &right):
    FlexibleCollimator(right),
    geometry(right.geometry)
{}


FlexibleCollimatorRep::FlexibleCollimatorRep(const std::string &name):
    FlexibleCollimator(name),
    geometry()
{}


FlexibleCollimatorRep::~FlexibleCollimatorRep()
{}


ElementBase *FlexibleCollimatorRep::clone() const {
    return new FlexibleCollimatorRep(*this);
}


Channel *FlexibleCollimatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<FlexibleCollimatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &FlexibleCollimatorRep::getField() {
    return field;
}

const NullField &FlexibleCollimatorRep::getField() const {
    return field;
}


StraightGeometry &FlexibleCollimatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &FlexibleCollimatorRep::getGeometry() const {
    return geometry;
}