//
// Class Degrader
//   Defines the abstract interface for a beam degrader.
//
// Copyright (c) 2000 - 2023, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/Degrader.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Solvers/ParticleMatterInteractionHandler.h"
#include "Utilities/Options.h"

#include <cmath>
#include <memory>

extern Inform* gmsg;


Degrader::Degrader():
    Degrader("")
{}

Degrader::Degrader(const Degrader& right):
    Component(right),
    PosX_m(right.PosX_m),
    PosY_m(right.PosY_m),
    PosZ_m(right.PosZ_m),
    MomentumX_m(right.MomentumX_m),
    MomentumY_m(right.MomentumY_m),
    MomentumZ_m(right.MomentumZ_m),
    time_m(right.time_m),
    id_m(right.id_m),
    width_m(right.width_m),
    height_m(right.height_m)
{}

Degrader::Degrader(const std::string& name):
    Component(name),
    width_m(0.0),
    height_m(0.0)
{ }

Degrader::~Degrader() {
    if(online_m)
        goOffline();
}

void Degrader::accept(BeamlineVisitor& visitor) const {
    visitor.visitDegrader(*this);
}

void Degrader::setDimensions(double xsize, double ysize) {
    width_m  = xsize;
    height_m = ysize;
}


bool Degrader::isInside(const Vector_t& R) const {
    /**
     check if the particle is in the degrader material
    */
    bool hit = false;
    if ( R(2) > 0.0 && R(2) <= getElementLength() &&
         (4 * (std::pow(R(0) / width_m, 2) + std::pow(R(1) / height_m, 2)) <= 1) ) {
        hit = true;
    }
    return hit;
}

bool Degrader::apply(const size_t& i, const double& t, Vector_t& /*E*/, Vector_t& /*B*/) {

    const Vector_t& R = RefPartBunch_m->R[i];
    const Vector_t& P = RefPartBunch_m->P[i];

    if ( isInside(R) ) {
        //if particle was already marked as -1 (it means it should have gone into degrader but didn't)
        //set the label to -2 (will not go into degrader and will be deleted when particles per core > 2)
        if (RefPartBunch_m->Bin[i] < 0) {
            RefPartBunch_m->Bin[i] = -2;
        } else {
            RefPartBunch_m->Bin[i] = -1;
        }
        const double& dt = RefPartBunch_m->dt[i];
        const double recpgamma = Physics::c * dt / std::sqrt(1.0  + dot(P, P));
        double frac = -R(2) / (P(2) * recpgamma);
        PosX_m.push_back(R(0));
        PosY_m.push_back(R(1));
        PosZ_m.push_back(R(2));
        MomentumX_m.push_back(P(0));
        MomentumY_m.push_back(P(1));
        MomentumZ_m.push_back(P(2));
        time_m.push_back(t + frac * dt);
        id_m.push_back(RefPartBunch_m->ID[i]);
    }

    return false;
}

bool Degrader::applyToReferenceParticle(const Vector_t& R,
                                        const Vector_t& P,
                                        const double& /*t*/,
                                        Vector_t& E,
                                        Vector_t& /*B*/) {
    if ( !isInside(R) ) return false;

    Vector_t updatedP = P;
    bool isDead = getParticleMatterInteraction()->computeEnergyLoss(RefPartBunch_m, updatedP, RefPartBunch_m->getdT(), false);
    double deltaP = euclidean_norm(updatedP) - euclidean_norm(P);
    E(2) += deltaP * RefPartBunch_m->getM() / (RefPartBunch_m->getdT() * RefPartBunch_m->getQ() * Physics::c);

    return isDead;
}

void Degrader::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {
    initialise(bunch);
    endField = startField + getElementLength();
}

void Degrader::initialise(PartBunchBase<double, 3>* bunch) {
    RefPartBunch_m = bunch;
}

void Degrader::finalise() {
    *gmsg << "* Finalize degrader " << getName() << endl;
}

void Degrader::goOnline(const double&) {
    Inform msg("Degrader::goOnline ");

    int maximumSize = (int)(1.1 * RefPartBunch_m->getLocalNum());

    PosX_m.reserve(maximumSize);
    PosY_m.reserve(maximumSize);
    PosZ_m.reserve(maximumSize);
    MomentumX_m.reserve(maximumSize);
    MomentumY_m.reserve(maximumSize);
    MomentumZ_m.reserve(maximumSize);
    time_m.reserve(maximumSize);
    id_m.reserve(maximumSize);
    online_m = true;
}

void Degrader::goOffline() {
    Inform msg("Degrader::goOffline ");
    online_m = false;
    msg << " done..." << endl;
}

bool Degrader::bends() const {
    return false;
}

void Degrader::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementType Degrader::getType() const {
    return ElementType::DEGRADER;
}
