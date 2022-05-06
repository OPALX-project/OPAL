//
// Class ParticleMatterInteractionHandler
//   Defines the handler for particle-matter interactions in the elements
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
#ifndef PARTICLEMATTERINTERACTIONHANDLER_HH
#define PARTICLEMATTERINTERACTIONHANDLER_HH

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"

#include <string>

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;
class Inform;

struct InsideTester {
    virtual ~InsideTester() {}

    virtual bool checkHit(const Vector_t& R) = 0;
};

class ParticleMatterInteractionHandler {

public:
    ParticleMatterInteractionHandler(std::string name, ElementBase* elref);
    virtual ~ParticleMatterInteractionHandler() { };
    virtual void apply(PartBunchBase<double, 3>* bunch,
                       const std::pair<Vector_t, double>& boundingSphere) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    virtual bool stillActive() = 0;
    virtual double getTime() = 0;
    virtual std::string getName() = 0;
    virtual size_t getParticlesInMat() = 0;
    virtual unsigned getRediffused() = 0;
    virtual unsigned int getNumEntered() = 0;
    void setFlagAllParticlesIn(bool p);
    bool getFlagAllParticlesIn() const;
    void updateElement(ElementBase* newref);
    ElementBase* getElement();

    virtual bool computeEnergyLoss(PartBunchBase<double, 3>* bunch,
                                   Vector_t& P,
                                   const double deltat,
                                   bool includeFluctuations = true) const = 0;

protected:
    ElementBase* element_ref_m;

    bool allParticleInMat_m; ///< if all particles are in matter stay inside the particle matter interaction

    const std::string name_m; // label of PARTICLEMATTERINTERACTION

    std::unique_ptr<InsideTester> hitTester_m; // tests whether particles are inside material
};

inline
ParticleMatterInteractionHandler::ParticleMatterInteractionHandler(std::string name, ElementBase* elref):
    element_ref_m(elref),
    allParticleInMat_m(false),
    name_m(name)
{}

inline
void ParticleMatterInteractionHandler::updateElement(ElementBase* newref) {
    element_ref_m = newref;
}

inline
ElementBase* ParticleMatterInteractionHandler::getElement() {
    return element_ref_m;
}

inline
void ParticleMatterInteractionHandler::setFlagAllParticlesIn(bool p) {
  allParticleInMat_m = p;
}

inline
bool ParticleMatterInteractionHandler::getFlagAllParticlesIn() const {
    return allParticleInMat_m;
}
#endif // PARTICLEMATTERINTERACTION_HH
