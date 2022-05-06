//
// Class Source
//   Defines the abstract interface for a source.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/Source.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Elements/OpalBeamline.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <iostream>
#include <fstream>


Source::Source():
    Source("")
{}

Source::Source(const Source& right):
    Component(right),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    isTransparent_m(right.isTransparent_m)
{}

Source::Source(const std::string& name):
    Component(name),
    startField_m(0.0),
    endField_m(0.0),
    isTransparent_m(false)
{}

Source::~Source() {
}

void Source::accept(BeamlineVisitor& visitor) const {
    visitor.visitSource(*this);
}

bool Source::apply(const size_t& i, const double& t, Vector_t& /*E*/, Vector_t& /*B*/) {

    if (isTransparent_m) {
        return false;
    }

    const Vector_t& R = RefPartBunch_m->R[i];
    const Vector_t& P = RefPartBunch_m->P[i];

    if (online_m && R(2) <= 0.0 && P(2) < 0.0) {
        const double& dt = RefPartBunch_m->dt[i];
        const Vector_t singleStep = Physics::c * dt * Util::getBeta(P);
        double frac = -R(2) / singleStep(2);

        lossDs_m->addParticle(OpalParticle(RefPartBunch_m->ID[i],
                                           R + frac * singleStep, P,
                                           t + frac * dt,
                                           RefPartBunch_m->Q[i], RefPartBunch_m->M[i]));

        return true;
    }

    return false;
}

void Source::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    endField = startField;
    startField -= getElementLength();

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getOutputFN(),
                                                              !Options::asciidump));
}

void Source::finalise()
{}

bool Source::bends() const {
    return false;
}

void Source::goOnline(const double&) {
    online_m = true;
}

void Source::goOffline() {
    online_m = false;
    lossDs_m->save();
}

void Source::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}

ElementType Source::getType() const {
    return ElementType::SOURCE;
}

void Source::setTransparent() {
    isTransparent_m = true;
}
