/*
 *  Copyright (c) 2014, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "AbsBeamline/VariableRFCavity.h"

#include "Physics/Physics.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/GeneralClassicException.h"

const double VariableRFCavity::lengthUnit_m = 1e3; // metres -> mm

VariableRFCavity::VariableRFCavity(const std::string &name) : Component(name) {
    initNull();  // initialise everything to NULL
}

VariableRFCavity::VariableRFCavity() : Component() {
    initNull();  // initialise everything to NULL
}

VariableRFCavity::VariableRFCavity(const VariableRFCavity& var) : Component() {
    initNull();  // initialise everything to NULL
    *this = var;
}

VariableRFCavity& VariableRFCavity::operator=(const VariableRFCavity& rhs) {
    if (&rhs == this) {
        return *this;
    }
    setName(rhs.getName());
    setPhaseModel(NULL);
    setAmplitudeModel(NULL);
    setFrequencyModel(NULL);
    if (rhs.phaseTD_m != NULL)
        setPhaseModel(std::shared_ptr<AbstractTimeDependence>(rhs.phaseTD_m->clone()));
    if (rhs.amplitudeTD_m != NULL) {
        setAmplitudeModel(std::shared_ptr<AbstractTimeDependence>(rhs.amplitudeTD_m->clone()));
    }
    if (rhs.frequencyTD_m != NULL)
        setFrequencyModel(std::shared_ptr<AbstractTimeDependence>(rhs.frequencyTD_m->clone()));
    phaseName_m = rhs.phaseName_m;
    amplitudeName_m = rhs.amplitudeName_m;
    frequencyName_m = rhs.frequencyName_m;
    halfWidth_m = rhs.halfWidth_m;
    halfHeight_m = rhs.halfHeight_m;
    setLength(rhs._length/lengthUnit_m);
    return *this;
}

VariableRFCavity::~VariableRFCavity() {
  // shared_ptr should self-destruct when they are ready
}

void VariableRFCavity::initNull() {
  _length = 0.;
  phaseName_m = "";
  amplitudeName_m = "";
  frequencyName_m = "";
  halfHeight_m = 0.;
  halfWidth_m = 0;
  RefPartBunch_m = NULL;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getAmplitudeModel() const {
    return amplitudeTD_m;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getPhaseModel() const {
    return phaseTD_m;
}

std::shared_ptr<AbstractTimeDependence> VariableRFCavity::getFrequencyModel() const {
    return frequencyTD_m;
}

void VariableRFCavity::setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitude_td) {
    amplitudeTD_m = amplitude_td;
}

void VariableRFCavity::setPhaseModel(std::shared_ptr<AbstractTimeDependence> phase_td) {
    phaseTD_m = phase_td;
}

void VariableRFCavity::setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequency_td) {
    frequencyTD_m = frequency_td;
}

StraightGeometry &VariableRFCavity::getGeometry() {
    return geometry;
}

const StraightGeometry &VariableRFCavity::getGeometry() const {
    return geometry;
}

EMField &VariableRFCavity::getField() {
  throw GeneralClassicException("VariableRFCavity",
                      "No field defined for VariableRFCavity");
}

const EMField &VariableRFCavity::getField() const {
  throw GeneralClassicException("VariableRFCavity",
                      "No field defined for VariableRFCavity");
}


bool VariableRFCavity::apply(const size_t &i, const double &t,
                             Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

// If this is too slow: a quicker implementation would be to use templates not
// inheritance (vtable lookup is removed). This is in the inner
// tracking loop, so low level optimisation is possibly worthwhile.
//
// Do I need bound checking here? I have no "radius" parameter, but I do have a
// "length".
bool VariableRFCavity::apply(const Vector_t &R, const Vector_t &/*P*/,
                             const double &t, Vector_t &E, Vector_t &/*B*/) {
    if (R[2] >= 0. && R[2] < _length) {
        if (std::abs(R[0]) > halfWidth_m || std::abs(R[1]) > halfHeight_m) {
            return true;
        }

        double E0 = amplitudeTD_m->getValue(t);
        double f = frequencyTD_m->getValue(t) * 1.0E-3; // need GHz on the element we have MHz
        double phi = phaseTD_m->getValue(t);
        E = Vector_t(0., 0., E0*sin(Physics::two_pi * f * t + phi));
        return false;
    }
    return true;
}

bool VariableRFCavity::applyToReferenceParticle(const Vector_t &R, const Vector_t &P,
                                                const double &t, Vector_t &E, Vector_t &B) {
    return apply(R, P, t, E, B);
}

void VariableRFCavity::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void VariableRFCavity::finalise() {
    RefPartBunch_m = NULL;
}

ElementBase* VariableRFCavity::clone() const {
    return new VariableRFCavity(*this);
}

void VariableRFCavity::accept(BeamlineVisitor& visitor) const {
    initialise();
    visitor.visitVariableRFCavity(*this);
}

void VariableRFCavity::initialise() const {
    VariableRFCavity* cavity = const_cast<VariableRFCavity*>(this);
    std::shared_ptr<AbstractTimeDependence> phaseTD =
                    AbstractTimeDependence::getTimeDependence(phaseName_m);
    cavity->setPhaseModel(std::shared_ptr<AbstractTimeDependence>(phaseTD->clone()));
    std::shared_ptr<AbstractTimeDependence> frequencyTD =
                    AbstractTimeDependence::getTimeDependence(frequencyName_m);
    cavity->setFrequencyModel(std::shared_ptr<AbstractTimeDependence>(frequencyTD->clone()));
        std::shared_ptr<AbstractTimeDependence> amplitudeTD =
                    AbstractTimeDependence::getTimeDependence(amplitudeName_m);
    cavity->setAmplitudeModel(std::shared_ptr<AbstractTimeDependence>(amplitudeTD->clone()));

    if (halfHeight_m < 1e-9 || halfWidth_m < 1e-9)
        throw GeneralClassicException("VariableRFCavity::accept",
                            "Height or width was not set on VariableRFCavity");
}

void VariableRFCavity::setLength(double length) {
    _length = length*lengthUnit_m;
    geometry.setElementLength(_length);
}
