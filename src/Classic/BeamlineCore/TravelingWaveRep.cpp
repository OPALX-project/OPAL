//
// Class TravelingWaveRep
//   Representation for a traveling wave.
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
#include "BeamlineCore/TravelingWaveRep.h"
#include "Channels/IndirectChannel.h"

bool TravelingWaveRep::ignoreCavities = false;

namespace {
    struct Entry {
        const char *name;
        double(TravelingWaveRep::*get)() const;
        void (TravelingWaveRep::*set)(double);
    };

    static const Entry entries[] = {
        {
            "L",
            &TravelingWaveRep::getElementLength,
            &TravelingWaveRep::setElementLength
        },
        {
            "AMPLITUDE",
            &TravelingWaveRep::getAmplitude,
            &TravelingWaveRep::setAmplitude
        },
        {
            "FREQUENCY",
            &TravelingWaveRep::getFrequency,
            &TravelingWaveRep::setFrequency
        },
        {
            "PHASE",
            &TravelingWaveRep::getPhase,
            &TravelingWaveRep::setPhase
        },
        { 0, 0, 0 }
    };
}


TravelingWaveRep::TravelingWaveRep():
    TravelingWave()
{}


TravelingWaveRep::TravelingWaveRep(const TravelingWaveRep &right):
    TravelingWave(right)
{}


TravelingWaveRep::TravelingWaveRep(const std::string &name):
    TravelingWave(name)
{}


TravelingWaveRep::~TravelingWaveRep()
{}


ElementBase *TravelingWaveRep::clone() const {
    return new TravelingWaveRep(*this);
}


Channel *TravelingWaveRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {

            return new IndirectChannel<TravelingWaveRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


AcceleratingField &TravelingWaveRep::getField() {
    return field;
}

const AcceleratingField &TravelingWaveRep::getField() const {
    return field;
}


StraightGeometry &TravelingWaveRep::getGeometry() {
    return geometry;
}

const StraightGeometry &TravelingWaveRep::getGeometry() const {
    return geometry;
}


double TravelingWaveRep::getAmplitude() const {
    return ignoreCavities ? 0.0 : field.getEz();
}


double TravelingWaveRep::getFrequency() const {
    return field.getFrequency();
}


double TravelingWaveRep::getPhase() const {
    return field.getPhase();
}


void TravelingWaveRep::setAmplitude(double amplitude) {
    field.setEz(amplitude);
}


void TravelingWaveRep::setFrequency(double frequency) {
    field.setFrequency(frequency);
}


void TravelingWaveRep::setPhase(double phase) {
    field.setPhase(phase);
}

void TravelingWaveRep::setIgnore(bool ignore) {
    ignoreCavities = ignore;
}
