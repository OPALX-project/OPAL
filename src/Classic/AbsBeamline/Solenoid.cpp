// ------------------------------------------------------------------------
// $RCSfile: Solenoid.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Solenoid.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class Solenoid
// ------------------------------------------------------------------------

Solenoid::Solenoid():
    Solenoid("")
{}


Solenoid::Solenoid(const Solenoid &right):
    Component(right),
    filename_m(right.filename_m),
    myFieldmap_m(right.myFieldmap_m),
    scale_m(right.scale_m),
    scaleError_m(right.scaleError_m),
    startField_m(right.startField_m),
    length_m(right.length_m),
    fast_m(right.fast_m) {
}


Solenoid::Solenoid(const std::string &name):
    Component(name),
    filename_m(""),
    myFieldmap_m(NULL),
    scale_m(1.0),
    scaleError_m(0.0),
    startField_m(0.0),
    length_m(0.0),
    fast_m(true) {
}


Solenoid::~Solenoid() {
    //    Fieldmap::deleteFieldmap(filename_m);
}


void Solenoid::accept(BeamlineVisitor &visitor) const {
    visitor.visitSolenoid(*this);
}

void Solenoid::setFieldMapFN(std::string fn) {
    filename_m = fn;
}

void Solenoid::setFast(bool fast) {
    fast_m = fast;
}


bool Solenoid::getFast() const {
    return fast_m;
}

bool Solenoid::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool Solenoid::apply(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {
    if (R(2) >= startField_m
        && R(2) < startField_m + length_m) {
        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = myFieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) return true;

        B += (scale_m + scaleError_m) * tmpB;
    }

    return false;
}

bool Solenoid::applyToReferenceParticle(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {

    if (R(2) >= startField_m
        && R(2) < startField_m + length_m) {
        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

        const bool outOfBounds = myFieldmap_m->getFieldstrength(R, tmpE, tmpB);
        if (outOfBounds) return true;

        B += scale_m * tmpB;
    }

    return false;
}

void Solenoid::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    Inform msg("Solenoid ", *gmsg);

    RefPartBunch_m = bunch;

    myFieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);

    if(myFieldmap_m != NULL) {
        msg << level2 << getName() << " using file ";
        myFieldmap_m->getInfo(&msg);

        double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
        myFieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

        startField_m = zBegin;
        length_m = zEnd - zBegin;
        endField = startField + length_m;
    } else {
        endField = startField;
    }
}

void Solenoid::finalise()
{}

bool Solenoid::bends() const {
    return false;
}


void Solenoid::goOnline(const double &) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void Solenoid::goOffline() {
    Fieldmap::freeMap(filename_m);
    online_m = false;
}

void Solenoid::setKS(double ks) {
    scale_m = ks;
}

void Solenoid::setDKS(double ks) {
    scaleError_m = ks;
}

void Solenoid::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + length_m;
}


ElementBase::ElementType Solenoid::getType() const {
    return SOLENOID;
}

bool Solenoid::isInside(const Vector_t &r) const {
    return isInsideTransverse(r)
        && myFieldmap_m->isInside(r);
}


double Solenoid::getElementLength() const {
    return length_m;
}

void Solenoid::getElementDimensions(double &begin,
                                         double &end) const {
    begin = startField_m;
    end = begin + length_m;
}