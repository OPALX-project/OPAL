// ------------------------------------------------------------------------
// $RCSfile: ParallelPlate.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelPlate
//   Defines the abstract interface for Parallel Plate. For benchmarking the
//   Secondary emisssion model
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2010/10/12 15:32:31 $
// $Author: Chuan Wang $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class ParallelPlate
// ------------------------------------------------------------------------

ParallelPlate::ParallelPlate():
    ParallelPlate("")
{}


ParallelPlate::ParallelPlate(const ParallelPlate &right):
    Component(right),
    filename_m(right.filename_m),
    scale_m(right.scale_m),
    phase_m(right.phase_m),
    frequency_m(right.frequency_m),
    length_m(right.length_m),
    ptime_m(0.0) {
    setElType(isRF);
}


ParallelPlate::ParallelPlate(const std::string &name):
    Component(name),
    filename_m(""),
    scale_m(1.0),
    phase_m(0.0),
    frequency_m(0.0),
    length_m(0.0),
    ptime_m(0.0) {
    setElType(isRF);
}


ParallelPlate::~ParallelPlate() {


}


void ParallelPlate::accept(BeamlineVisitor &visitor) const {
    visitor.visitParallelPlate(*this);
}


std::string ParallelPlate::getFieldMapFN() const {
    return "";
}

void ParallelPlate::setAmplitude(double vPeak) {
    scale_m = vPeak;
}

double ParallelPlate::getAmplitude() const {
    return scale_m ;
}

void ParallelPlate::setFrequency(double freq) {
    frequency_m = freq;
}

double ParallelPlate::getFrequency()const {
    return  frequency_m;
}

void ParallelPlate::setPhase(double phase) {
    phase_m = phase;
}

double ParallelPlate::getPhase() const {
    return phase_m;
}

// void ParallelPlate::setElementLength(double length) {
//     length_m = length;
// }

// double ParallelPlate::getElementLength() const {
//     return length_m;
// }


bool ParallelPlate::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return applyToReferenceParticle(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool ParallelPlate::apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {
    return applyToReferenceParticle(R, P, t, E, B);
}

bool ParallelPlate::applyToReferenceParticle(const Vector_t &R, const Vector_t &/*P*/, const double &t, Vector_t &E, Vector_t &B) {
    const double phase = frequency_m * t + phase_m;
    Vector_t tmpE(0.0, 0.0, -1.0), tmpB(0.0, 0.0, 0.0);
    ptime_m = t;

    //Here we don't check if the particle is outof bounds, just leave the matter to the particle Boundary collision model in BoundaryGeometry
    if ((R(2) >= 0.0)
        && (R(2) < length_m)) {
        E += scale_m / length_m * sin(phase) * tmpE; //Here scale_m should be Voltage between the parallel plates(V).
        B = tmpB;       //B field is always zero for our parallel plate elements used for benchmarking.
        return false;
    }
    return true;
}

void ParallelPlate::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {

    length_m = getElementLength();
    RefPartBunch_m = bunch;
    startField = 0.0;

    endField = length_m;

}

// In current version ,not implemented yet.
void ParallelPlate::initialise(PartBunchBase<double, 3> *bunch) {
    using Physics::pi;

    Inform msg("ParallelPlate initialization for cyclotron tracker ");

    RefPartBunch_m = bunch;
    // ElementEdge_m = startField;
    msg << " Currently parallelplate initialization for cyclotron tracker is  empty! " << endl;

}



void ParallelPlate::finalise()
{}

bool ParallelPlate::bends() const {
    return false;
}

void ParallelPlate::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = length_m;
}


ElementBase::ElementType ParallelPlate::getType() const {
    return PARALLELPLATE;
}
