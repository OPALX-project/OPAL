//
// Class VacuumRep
//   Representation for the vacuum conditions.
//
// Copyright (c) 2018 - 2021, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#include "BeamlineCore/VacuumRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double(VacuumRep::*get)() const;
        void (VacuumRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "p",
            &Vacuum::getPressure,
            &Vacuum::setPressure
        },
        { 0, 0, 0 }
    };
}


VacuumRep::VacuumRep():
    Vacuum(),
    geometry(0.0)
{}


VacuumRep::VacuumRep(const VacuumRep& right):
    Vacuum(right),
    geometry(right.geometry)
{}


VacuumRep::VacuumRep(const std::string& name):
    Vacuum(name),
    geometry()
{}


VacuumRep::~VacuumRep()
{}


ElementBase* VacuumRep::clone() const {
    return new VacuumRep(*this);
}


Channel* VacuumRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry *entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<VacuumRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField& VacuumRep::getField() {
    return field;
}

const NullField& VacuumRep::getField() const {
    return field;
}

StraightGeometry& VacuumRep::getGeometry() {
    return geometry;
}

const StraightGeometry& VacuumRep::getGeometry() const {
    return geometry;
}
