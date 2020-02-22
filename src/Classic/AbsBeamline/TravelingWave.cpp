//
// Class TravelingWave
//   Defines the abstract interface for Traveling Wave.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/TravelingWave.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

TravelingWave::TravelingWave():
    TravelingWave("")
{}


TravelingWave::TravelingWave(const TravelingWave& right):
    RFCavity(right),
    scaleCore_m(right.scaleCore_m),
    scaleCoreError_m(right.scaleCoreError_m),
    phaseCore1_m(right.phaseCore1_m),
    phaseCore2_m(right.phaseCore2_m),
    phaseExit_m(right.phaseExit_m),
    startCoreField_m(right.startCoreField_m),
    startExitField_m(right.startExitField_m),
    mappedStartExitField_m(right.mappedStartExitField_m),
    periodLength_m(right.periodLength_m),
    numCells_m(right.numCells_m),
    cellLength_m(right.cellLength_m),
    mode_m(right.mode_m)
{}


TravelingWave::TravelingWave(const std::string& name):
    RFCavity(name),
    scaleCore_m(1.0),
    scaleCoreError_m(0.0),
    phaseCore1_m(0.0),
    phaseCore2_m(0.0),
    phaseExit_m(0.0),
    startCoreField_m(0.0),
    startExitField_m(0.0),
    mappedStartExitField_m(0.0),
    periodLength_m(0.0),
    numCells_m(1),
    cellLength_m(0.0),
    mode_m(1)
{}


TravelingWave::~TravelingWave() {
}


void TravelingWave::accept(BeamlineVisitor& visitor) const {
    visitor.visitTravelingWave(*this);
}

bool TravelingWave::apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool TravelingWave::apply(const Vector_t& R,
                          const Vector_t& /*P*/,
                          const double& t,
                          Vector_t& E, Vector_t& B) {

    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength()) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * periodLength_m);
    double tmpcos, tmpsin;
    ComplexVector_t tmpE, tmpB;

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR)) return getFlagDeleteOnTransverseExit();

        double scale = scaleCore_m + scaleCoreError_m;
        tmpcos = scale * std::cos(frequency_m * t + phase_m + phaseError_m);
        tmpsin = scale * std::sin(frequency_m * t + phase_m + phaseError_m);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();
    } else if (tmpR(2) < startExitField_m) {
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;
        if (!fieldmap_m->isInside(tmpR)) return getFlagDeleteOnTransverseExit();

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        double scale = scaleCore_m + scaleCoreError_m;
        double phase1 = frequency_m * t + phaseCore1_m + phaseError_m;
        tmpcos = scale * std::cos(phase1);
        tmpsin = scale * std::sin(phase1);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        const double phase2 = frequency_m * t + phaseCore2_m + phaseError_m;
        tmpcos = scale * std::cos(phase2);
        tmpsin = scale * std::sin(phase2);
        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR))
            return true;

        ComplexVector_t tmpE, tmpB;
        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        const double phase = frequency_m * t + phaseExit_m + phaseError_m;
        const double scale = scale_m + scaleError_m;
        tmpcos = scale * std::cos(phase);
        tmpsin = scale * std::sin(phase);
        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();
    }

    return false;
}

bool TravelingWave::applyToReferenceParticle(const Vector_t& R,
                                             const Vector_t& /*P*/,
                                             const double& t, Vector_t& E,
                                             Vector_t& B) {

    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength()) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * periodLength_m);
    double tmpcos, tmpsin;

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR)) return true;

        ComplexVector_t tmpE, tmpB;
        const double phase = frequency_m * t + phase_m;
        double tmpcos = scale_m * std::cos(phase);
        double tmpsin = scale_m * std::sin(phase);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();

    } else if (tmpR(2) < startExitField_m) {
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;
        if (!fieldmap_m->isInside(tmpR)) return true;

        ComplexVector_t tmpE, tmpB;
        const double phase1 = frequency_m * t * phaseCore1_m;
        double tmpcos = scaleCore_m * std::cos(phase1);
        double tmpsin = scaleCore_m * std::sin(phase1);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        const double phase2 = frequency_m * t * phaseCore2_m;
        tmpcos = scaleCore_m * std::cos(phase2);
        tmpsin = scaleCore_m * std::sin(phase2);

        tmpE = 0.0;
        tmpB = 0.0;

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR)) return true;

        ComplexVector_t tmpE, tmpB;
        const double phase = frequency_m * t + phaseExit_m;
        tmpcos = scale_m * std::cos(phase);
        tmpsin = scale_m * std::sin(phase);

        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

        E += (std::complex<double>(tmpcos, tmpsin) * tmpE).real();
        B += (std::complex<double>(tmpcos, tmpsin) * tmpB).real();
    }

    return false;
}

void TravelingWave::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {

    if (bunch == nullptr) {
        startField = - 0.5 * periodLength_m;
        endField = startExitField_m;
        return;
    }

    Inform msg("TravelingWave ", *gmsg);

    RefPartBunch_m = bunch;
    double zBegin = 0.0, zEnd = 0.0;
    RFCavity::initialise(bunch, zBegin, zEnd);
    if (std::abs(startField_m) > 0.0) {
        throw GeneralClassicException("TravelingWave::initialise",
                                      "The field map of a traveling wave structure has to begin at 0.0");
    }

    periodLength_m = (zEnd - zBegin) / 2.0;
    cellLength_m = periodLength_m * mode_m;
    startField_m = -0.5 * periodLength_m;

    startCoreField_m = periodLength_m / 2.0;
    startExitField_m = startCoreField_m + (numCells_m - 1) * cellLength_m;
    mappedStartExitField_m = startExitField_m - 3.0 * periodLength_m / 2.0;

    startField = -periodLength_m / 2.0;
    endField = startField + startExitField_m + periodLength_m / 2.0;
    setElementLength(endField - startField);

    scaleCore_m = scale_m / std::sin(Physics::two_pi * mode_m);
    scaleCoreError_m = scaleError_m / std::sin(Physics::two_pi * mode_m);
    phaseCore1_m = phase_m + Physics::pi * mode_m / 2.0;
    phaseCore2_m = phase_m + Physics::pi * mode_m * 1.5;
    phaseExit_m = phase_m - Physics::two_pi * ((numCells_m - 1) * mode_m - std::floor((numCells_m - 1) * mode_m));
}

void TravelingWave::finalise()
{}

bool TravelingWave::bends() const {
    return false;
}

void TravelingWave::goOnline(const double&) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void TravelingWave::goOffline() {
    Fieldmap::freeMap(filename_m);
}

void TravelingWave::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = -0.5 * periodLength_m;
    zEnd = zBegin + getElementLength();
}


void TravelingWave::getElementDimensions(double& begin, double& end) const {
    begin = -0.5 * periodLength_m;
    end = begin + getElementLength();
}

ElementType TravelingWave::getType() const {
    return ElementType::TRAVELINGWAVE;
}

bool TravelingWave::isInside(const Vector_t& r) const {
    return (isInsideTransverse(r)
            && r(2) >= -0.5 * periodLength_m
            && r(2) < startExitField_m);
}