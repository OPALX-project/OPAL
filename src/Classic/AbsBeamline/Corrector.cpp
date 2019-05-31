// ------------------------------------------------------------------------
// $RCSfile: Corrector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Corrector
//   Defines the abstract interface for a orbit corrector.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

extern Inform *gmsg;

// Class Corrector
// ------------------------------------------------------------------------

Corrector::Corrector():
    Corrector("")
{ }


Corrector::Corrector(const Corrector &right):
    Component(right),
    kickX_m(right.kickX_m),
    kickY_m(right.kickY_m),
    designEnergy_m(right.designEnergy_m),
    designEnergyChangeable_m(right.designEnergyChangeable_m),
    kickFieldSet_m(right.kickFieldSet_m),
    kickField_m(right.kickField_m)
{ }


Corrector::Corrector(const std::string &name):
    Component(name),
    kickX_m(0.0),
    kickY_m(0.0),
    designEnergy_m(0.0),
    designEnergyChangeable_m(true),
    kickFieldSet_m(false),
    kickField_m(0.0)
{ }


Corrector::~Corrector()
{ }


void Corrector::accept(BeamlineVisitor &visitor) const {
    visitor.visitCorrector(*this);
}

bool Corrector::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    Vector_t &R = RefPartBunch_m->R[i];
    Vector_t &P = RefPartBunch_m->P[i];

    return apply(R, P, t, E, B);
}

bool Corrector::apply(const Vector_t &R,
                      const Vector_t &P,
                      const double &t,
                      Vector_t &E,
                      Vector_t &B) {

    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;

        double tau = 1.0;
        const double &dt = RefPartBunch_m->getdT();
        const double stepSize = dt * Physics::c * P(2) / Util::getGamma(P);

        if (R(2) < stepSize) {
            tau = R(2) / stepSize + 0.5;
        }
        if (getElementLength() - R(2) < stepSize) {
            tau += (getElementLength() - R(2)) / stepSize - 0.5;
        }

        B += kickField_m * tau;
    }

    return false;
}

void Corrector::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
}

void Corrector::finalise()
{ }

void Corrector::goOnline(const double &) {
    const double pathLength = getGeometry().getElementLength();
    const double minLength = Physics::c * RefPartBunch_m->getdT();
    if (pathLength < minLength) {
        throw GeneralClassicException("Corrector::goOnline",
                                      "length of corrector, L= " + std::to_string(pathLength) +
                                      ", shorter than distance covered during one time step, dS= " + std::to_string(minLength));
    }

    if (!kickFieldSet_m) {
        const double momentum = sqrt(std::pow(designEnergy_m, 2.0) + 2.0 * designEnergy_m * RefPartBunch_m->getM());
        const double magnitude = momentum / (Physics::c * pathLength);
        kickField_m = magnitude * RefPartBunch_m->getQ() * Vector_t(kickY_m, -kickX_m, 0.0);
    }

    online_m = true;
}

void Corrector::setDesignEnergy(const double& ekin, bool changeable) {
    if (designEnergyChangeable_m) {
        designEnergy_m = ekin;
        designEnergyChangeable_m = changeable;
    }
    if (RefPartBunch_m) {
        if (!kickFieldSet_m) {
            const double pathLength = getGeometry().getElementLength();
            const double momentum = sqrt(std::pow(designEnergy_m, 2.0) + 2.0 * designEnergy_m * RefPartBunch_m->getM());
            const double magnitude = momentum / (Physics::c * pathLength);
            kickField_m = magnitude * RefPartBunch_m->getQ() * Vector_t(kickY_m, -kickX_m, 0.0);
        }
    }
}

bool Corrector::bends() const {
    return false;
}

void Corrector::getDimensions(double &zBegin, double &zEnd) const
{
  zBegin = 0.0;
  zEnd = getElementLength();
}

ElementBase::ElementType Corrector::getType() const {
    return CORRECTOR;
}