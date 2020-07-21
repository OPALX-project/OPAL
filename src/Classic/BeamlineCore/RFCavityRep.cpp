//
// Class RFCavityRep
//   Representation for a RF cavity.
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
#include "BeamlineCore/RFCavityRep.h"
#include "Channels/IndirectChannel.h"


bool RFCavityRep::ignoreCavities = false;

namespace {
    struct Entry {
        const char *name;
        double(RFCavityRep::*get)() const;
        void (RFCavityRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &RFCavityRep::getElementLength,
            &RFCavityRep::setElementLength
        },
        {
            "AMPLITUDE",
            &RFCavityRep::getAmplitude,
            &RFCavityRep::setAmplitude
        },
        {
            "FREQUENCY",
            &RFCavityRep::getFrequency,
            &RFCavityRep::setFrequency
        },
        {
            "PHASE",
            &RFCavityRep::getPhase,
            &RFCavityRep::setPhase
        },
        { 0, 0, 0 }
    };
}


RFCavityRep::RFCavityRep():
    RFCavity()
{}


RFCavityRep::RFCavityRep(const RFCavityRep &right):
    RFCavity(right),
    geometry(right.geometry)
{}


RFCavityRep::RFCavityRep(const std::string &name):
    RFCavity(name)
{}


RFCavityRep::~RFCavityRep()
{}


ElementBase *RFCavityRep::clone() const {
    return new RFCavityRep(*this);
}


Channel *RFCavityRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<RFCavityRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


AcceleratingField &RFCavityRep::getField() {
    return field;
}

const AcceleratingField &RFCavityRep::getField() const {
    return field;
}


StraightGeometry &RFCavityRep::getGeometry() {
    return geometry;
}

const StraightGeometry &RFCavityRep::getGeometry() const {
    return geometry;
}


double RFCavityRep::getAmplitude() const {
    return ignoreCavities ? 0.0 : field.getEz();
}


double RFCavityRep::getFrequency() const {
    return field.getFrequency();
}


double RFCavityRep::getPhase() const {
    return field.getPhase();
}


void RFCavityRep::setAmplitude(double amplitude) {
    field.setEz(amplitude);
}


void RFCavityRep::setFrequency(double frequency) {
    field.setFrequency(frequency);
}


void RFCavityRep::setPhase(double phase) {
    field.setPhase(phase);
}

void RFCavityRep::setIgnore(bool ignore) {
    ignoreCavities = ignore;
}