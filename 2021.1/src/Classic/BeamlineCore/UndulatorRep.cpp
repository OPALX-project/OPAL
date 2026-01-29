//
// Class UndulatorRep
//   Defines a concrete undulator/wiggler representation.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#include "BeamlineCore/UndulatorRep.h"

#include "Channels/IndirectChannel.h"

// Attribute access table.
namespace {
    struct Entry {
        const char* name;
        double (UndulatorRep::*get)() const;
        void (UndulatorRep::*set)(double);
    };

    const Entry entries[] = {
        {"L", &UndulatorRep::getElementLength, &UndulatorRep::setElementLength}, {0, 0, 0}};
}  // namespace

UndulatorRep::UndulatorRep() : Undulator(), geometry(0.0) {
}

UndulatorRep::UndulatorRep(const UndulatorRep& right) : Undulator(right), geometry(right.geometry) {
}

UndulatorRep::UndulatorRep(const std::string& name) : Undulator(name), geometry(0.0) {
}

UndulatorRep::~UndulatorRep() {
}

ElementBase* UndulatorRep::clone() const {
    return new UndulatorRep(*this);
}

Channel* UndulatorRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<UndulatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField& UndulatorRep::getField() {
    return field;
}

const NullField& UndulatorRep::getField() const {
    return field;
}

StraightGeometry& UndulatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry& UndulatorRep::getGeometry() const {
    return geometry;
}