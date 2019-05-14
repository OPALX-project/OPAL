// ------------------------------------------------------------------------
// $RCSfile: RBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class RBend
// ------------------------------------------------------------------------

RBend::RBend():
    Bend2D()
{
    setMessageHeader("RBend ");
}

RBend::RBend(const RBend &right):
    Bend2D(right)
{
    setMessageHeader("RBend ");
}

RBend::RBend(const std::string &name):
    Bend2D(name)
{
    setMessageHeader("RBend ");
}

RBend::~RBend() {
}

void RBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitRBend(*this);
}

/*
 * OPAL-MAP methods
 * ================
 */
double RBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}

double RBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}

void RBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}

void RBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}

/*
 * BET methods.
 */
void RBend::addKR(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - getStartField();
    const Vector_t tmpA(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);

    DiffDirection zdir(DZ);
    getFieldmap()->getFieldstrength(tmpA, tmpE, tmpB);
    getFieldmap()->getFieldDerivative(tmpA, tmpE_diff, tmpB_diff, zdir);

    double g = RefPartBunch_m->getGamma(i);

    if(std::abs(getFieldAmplitude() * tmpB_diff(2)) > 0.1) {
        double cf = Physics::q_e * tmpB(2) / (g * Physics::EMASS);
        K += Vector_t(-pow(cf * getFieldAmplitude() * tmpB(0), 2) / 3.0, -pow(cf * getFieldAmplitude() * tmpB(1), 2) / 3.0, 0.0);
    }
}

void RBend::addKT(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - getStartField();
    const Vector_t tmpA(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), pz);
    getFieldmap()->getFieldstrength(tmpA, tmpE, tmpB);

    double b = RefPartBunch_m->getBeta(i);
    double g = 1 / sqrt(1 - b * b);

    double cf = -Physics::q_e * Physics::c * b * tmpB(2) * getFieldAmplitude() / (g * Physics::EMASS);
    Vector_t temp(cf * tmpB(1), cf * tmpB(0), 0.0);

    //FIXME: K += ??
}


/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */

ElementBase::ElementType RBend::getType() const {
    return RBEND;
}

void RBend::setBendAngle(double angle) {
    Bend2D::setBendAngle(angle);
    setExitAngle(angle - getEntranceAngle());
}

void RBend::setEntranceAngle(double entranceAngle) {
    Bend2D::setEntranceAngle(entranceAngle);
    setExitAngle(getBendAngle() - entranceAngle);
}

bool RBend::findChordLength(Inform &msg,
                            double &chordLength) {

    /*
     * Find bend chord length. If this was not set by the user using the
     * L (length) attribute, infer it from the field map.
     */
    const double angle = getBendAngle();
    if (std::abs(angle) > 0.0) {
        double E1 = std::copysign(1.0, angle) * getEntranceAngle();
        chordLength = 2 * getLength() * sin(0.5 * std::abs(angle)) /
            (sin(E1) + sin(std::abs(angle) - E1));
    } else {
        double refMass  = RefPartBunch_m->getM();
        double refGamma = designEnergy_m / refMass + 1.0;
        double refBetaGamma = sqrt(pow(refGamma, 2.0) - 1.0);
        double refCharge = RefPartBunch_m->getQ();
        double amplitude = (std::abs(bY_m) > 0.0? bY_m: bX_m);
        double fieldAmplitude = refCharge * std::abs(amplitude / refCharge);
        double designRadius = std::abs(refBetaGamma * refMass / (Physics::c * fieldAmplitude));
        chordLength = sin(0.5 * (asin(getLength() / designRadius - sin(getEntranceAngle())) + getEntranceAngle())) * 2 * designRadius;
    }

    return true;
}