//
// Class ProbeRep
//   Representation for Probe.
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
#include "BeamlineCore/ProbeRep.h"
#include "Channels/IndirectChannel.h"


namespace {
    struct Entry {
        const char *name;
        double(ProbeRep::*get)() const;
        void (ProbeRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &ProbeRep::getElementLength,
            &ProbeRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


ProbeRep::ProbeRep():
    Probe(), field(), geometry(), active(true)
{}


ProbeRep::ProbeRep(const ProbeRep &right):
    Probe(right), field(), geometry(right.geometry), active(true)
{}


ProbeRep::ProbeRep(const std::string &name):
    Probe(name), field(), geometry(), active(true)
{}


ProbeRep::~ProbeRep()
{}


ElementBase *ProbeRep::clone() const {
    return new ProbeRep(*this);
}


Channel *ProbeRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<ProbeRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField &ProbeRep::getField() {
    return field;
}

const NullField &ProbeRep::getField() const {
    return field;
}

StraightGeometry &ProbeRep::getGeometry() {
    return geometry;
}

const StraightGeometry &ProbeRep::getGeometry() const {
    return geometry;
}


void ProbeRep::setActive(bool flag) {
    active = flag;
}
