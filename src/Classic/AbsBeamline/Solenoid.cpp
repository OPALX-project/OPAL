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
    setElType(isSolenoid);
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
    setElType(isSolenoid);
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

/**
 * \brief ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the
 * solenoid element and adds it to the K vector
*/
void Solenoid::addKR(int i, double /*t*/, Vector_t &K) {
    Inform msg("Solenoid::addKR()");

    double pz = RefPartBunch_m->getZ(i) - startField_m;
    if (pz < 0.0 ||
        pz >= length_m) return;

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    const Vector_t tmpA(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);

    myFieldmap_m->getFieldstrength(tmpA, tmpE, tmpB);
    double k = Physics::q_e * scale_m * tmpB(2) / (2.0 * Physics::EMASS * RefPartBunch_m->getGamma(i));
    k *= k;
    K += Vector_t(k, k, 0.0);
}

/**
 * ENVELOPE COMPONENT for transverse kick (only important for x0, y0)
 * Calculates the transverse kick component for the solenoid element and adds it to
 * the K vector, only important for off track tracking. Otherwise KT = 0.
*/
void Solenoid::addKT(int i, double /*t*/, Vector_t &K) {
    Inform msg("Solenoid::addKT()");
    double dbdz, emg;

    double pz = RefPartBunch_m->getZ(i) - startField_m;
    if (pz < 0.0 ||
        pz >= length_m) return;

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);
    const Vector_t tmpA(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);

    // define z direction:
    DiffDirection zdir(DZ);

    myFieldmap_m->getFieldstrength(tmpA, tmpE, tmpB);

    // get derivation of B in z-direction
    myFieldmap_m->getFieldDerivative(tmpA, tmpE_diff, tmpB_diff, zdir);

    double bz = scale_m * tmpB(2);
    double g = RefPartBunch_m->getGamma(i);

    //FIXME?: BET: dz   = z - z0,
    dbdz = scale_m * tmpB_diff(2) * RefPartBunch_m->getBeta(i) * Physics::c;
    emg  = Physics::q_e / (g * Physics::EMASS);

    /** BET:
     * Vector_t temp(emg*(bz*(Cxy*tempBunch->getPy(i) + Sy) + Cxy*dbdz*(y-y0)),
     *              -emg*(bz*(Cxy*tempBunch->getPx(i) + Sx) + Cxy*dbdz*(x-x0)),
     *              0.0);
    */

    double dx = RefPartBunch_m->getX0(i);
    double dy = RefPartBunch_m->getY0(i);

    //FIXME: Py0, Px0?
    Vector_t temp(emg * (bz * (RefPartBunch_m->getPy0(i)) + dbdz * dy),
                  -emg * (bz * (RefPartBunch_m->getPx0(i)) + dbdz * dx),
                  0.0);

    K += temp;
}

bool Solenoid::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool Solenoid::apply(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {
    const Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) >= 0.0 && tmpR(2) < length_m) {
        const bool out_of_bounds = myFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        if (out_of_bounds) return true;

        B += (scale_m + scaleError_m) * tmpB;
    }

    return false;
}

bool Solenoid::applyToReferenceParticle(const Vector_t &R, const Vector_t &/*P*/, const  double &/*t*/, Vector_t &/*E*/, Vector_t &B) {
    const Vector_t tmpR(R(0), R(1), R(2) - startField_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) >= 0.0 && tmpR(2) < length_m) {
        const bool out_of_bounds = myFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        if (out_of_bounds) return true;

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
    return (r(2) >= startField_m && r(2) < startField_m + length_m && isInsideTransverse(r) && myFieldmap_m->isInside(r));
}


double Solenoid::getElementLength() const {
    return length_m;
}

void Solenoid::getElementDimensions(double &begin,
                                         double &end) const {
    begin = startField_m;
    end = begin + length_m;
}

