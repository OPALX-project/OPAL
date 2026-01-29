//
// Class SourceRep
//   Representation for a source.
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
#include "BeamlineCore/SourceRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char *name;
        double(SourceRep::*get)() const;
        void (SourceRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &SourceRep::getElementLength,
            &SourceRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


SourceRep::SourceRep():
    Source(),
    geometry()
{}


SourceRep::SourceRep(const SourceRep &right):
    Source(right),
    geometry(right.geometry)
{}


SourceRep::SourceRep(const std::string &name):
    Source(name), geometry()
{}


SourceRep::~SourceRep()
{}


ElementBase *SourceRep::clone() const {
    return new SourceRep(*this);
}


Channel *SourceRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<SourceRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &SourceRep::getField() {
    return field;
}

const NullField &SourceRep::getField() const {
    return field;
}

StraightGeometry &SourceRep::getGeometry() {
    return geometry;
}

const StraightGeometry &SourceRep::getGeometry() const {
    return geometry;
}