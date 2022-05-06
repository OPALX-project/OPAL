//
// Class Beam
//   The class for the OPAL BEAM command.
//   A BEAM definition is used by most physics commands to define the
//   particle charge and the reference momentum, together with some other data.
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
#include "Structure/Beam.h"

#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <iterator>

using namespace Expressions;

// The attributes of class Beam.
namespace {
    enum {
        PARTICLE,   // The particle name
        MASS,       // The particle rest mass in GeV
        CHARGE,     // The particle charge in proton charges
        ENERGY,     // The particle energy in GeV
        PC,         // The particle momentum in GeV/c
        GAMMA,      // ENERGY / MASS
        BCURRENT,   // Beam current in A
        BFREQ,      // Beam frequency in MHz
        NPART,      // Number of particles per bunch
        SIZE
    };
}


Beam::Beam():
    Definition(SIZE, "BEAM",
               "The \"BEAM\" statement defines data for the particles "
               "in a beam."),
    reference(1.0, Physics::m_p * Units::GeV2eV, 1.0 * Units::GeV2eV) {

    itsAttr[PARTICLE] = Attributes::makePredefinedString
                        ("PARTICLE", "Name of particle to be used",
                         {"ELECTRON",
                          "POSITRON",
                          "MUON",
                          "PROTON",
                          "ANTIPROTON",
                          "DEUTERON",
                          "HMINUS",
                          "H2P",
                          "ALPHA",
                          "CARBON",
                          "XENON",
                          "URANIUM"});

    itsAttr[MASS]     = Attributes::makeReal
                        ("MASS", "Particle rest mass [GeV]");

    itsAttr[CHARGE]   = Attributes::makeReal
                        ("CHARGE", "Particle charge in proton charges");

    itsAttr[ENERGY]   = Attributes::makeReal
                        ("ENERGY", "Particle energy [GeV]");

    itsAttr[PC]       = Attributes::makeReal
                        ("PC", "Particle momentum [GeV/c]");

    itsAttr[GAMMA]    = Attributes::makeReal
                        ("GAMMA", "ENERGY / MASS");

    itsAttr[BCURRENT] = Attributes::makeReal
                        ("BCURRENT", "Beam current [A] (all bunches)");

    itsAttr[BFREQ]    = Attributes::makeReal
                        ("BFREQ", "Beam frequency [MHz] (all bunches)");

    itsAttr[NPART]    = Attributes::makeReal
                        ("NPART", "Number of particles in bunch");

    // Set up default beam.
    Beam* defBeam = clone("UNNAMED_BEAM");
    defBeam->builtin = true;

    try {
        defBeam->update();
        OpalData::getInstance()->define(defBeam);
    } catch(...) {
        delete defBeam;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


Beam::Beam(const std::string& name, Beam* parent):
    Definition(name, parent),
    reference(parent->reference)
{}


Beam::~Beam()
{}


bool Beam::canReplaceBy(Object* object) {
    // Can replace only by another BEAM.
    return dynamic_cast<Beam*>(object) != 0;
}


Beam* Beam::clone(const std::string& name) {
    return new Beam(name, this);
}


void Beam::execute() {
    update();
    // Check if energy explicitly has been set with the BEAM command
    if (!itsAttr[GAMMA] && !(itsAttr[ENERGY]) && !(itsAttr[PC])) {
        throw OpalException("Beam::execute()",
                            "The energy hasn't been set. "
                            "Set either \"GAMMA\", \"ENERGY\" or \"PC\".");
    }

    if ( !(itsAttr[PARTICLE]) && (!itsAttr[MASS] || !(itsAttr[CHARGE])) ) {
        throw OpalException("Beam::execute()",
                            "The beam particle hasn't been set. "
                            "Set either \"PARTICLE\" or \"MASS\" and \"CHARGE\".");
    }

    if (!(itsAttr[NPART])) {
        throw OpalException("Beam::execute()", "\"NPART\" must be set.");
    }
}


Beam* Beam::find(const std::string& name) {
    Beam* beam = dynamic_cast<Beam*>(OpalData::getInstance()->find(name));

    if (beam == 0) {
        throw OpalException("Beam::find()", "Beam \"" + name + "\" not found.");
    }

    return beam;
}

size_t Beam::getNumberOfParticles() const {
    if (Attributes::getReal(itsAttr[NPART]) > 0) {
        return (size_t)Attributes::getReal(itsAttr[NPART]);
    } else {
        throw OpalException("Beam::getNumberOfParticles()",
                            "Wrong number of particles in beam!. \"NPART\" must be positive");
    }
}

const PartData& Beam::getReference() const {
    // Cast away const, to allow logically constant Beam to update.
    const_cast<Beam*>(this)->update();
    return reference;
}

double Beam::getCurrent() const {
    return Attributes::getReal(itsAttr[BCURRENT]);
}

double Beam::getCharge() const {
    return Attributes::getReal(itsAttr[CHARGE]);
}

double Beam::getMass() const {
    return Attributes::getReal(itsAttr[MASS]);
}

std::string Beam::getParticleName() const {
    return Attributes::getString(itsAttr[PARTICLE]);
}

double Beam::getFrequency() const {
    return Attributes::getReal(itsAttr[BFREQ]);
}

double Beam::getChargePerParticle() const {
    return std::copysign(1.0, getCharge()) * getCurrent()
        / (getFrequency() * Units::MHz2Hz)
        / getNumberOfParticles();
}

double Beam::getMassPerParticle() const {
    return getMass() * getChargePerParticle() / (getCharge() * Physics::q_e);
}

void Beam::update() {

    if (itsAttr[PARTICLE]) {
        std::string pName  = getParticleName();
        ParticleType pType = ParticleProperties::getParticleType(pName);
        Attributes::setReal(itsAttr[MASS], ParticleProperties::getParticleMass(pType));
        Attributes::setReal(itsAttr[CHARGE], ParticleProperties::getParticleCharge(pType));
    }

    // Set up particle reference; convert all to eV for CLASSIC.
    double mass = (itsAttr[MASS] ? getMass() : Physics::m_p) * Units::GeV2eV;
    double charge = itsAttr[CHARGE] ? getCharge() : 1.0;

    reference = PartData(charge, mass, 1.0);

    if (itsAttr[GAMMA]) {
        double gamma = Attributes::getReal(itsAttr[GAMMA]);
        if (gamma > 1.0) {
            reference.setGamma(gamma);
        } else {
            throw OpalException("Beam::update()",
                                "\"GAMMA\" should be greater than 1.");
        }
    } else if (itsAttr[ENERGY]) {
        double energy = Attributes::getReal(itsAttr[ENERGY]) * Units::GeV2eV;
        if (energy > reference.getM()) {
            reference.setE(energy);
        } else {
            throw OpalException("Beam::update()",
                                "\"ENERGY\" should be greater than \"MASS\".");
        }
    } else if (itsAttr[PC]) {
        double pc = Attributes::getReal(itsAttr[PC]) * Units::GeV2eV;
        if (pc > 0.0) {
            reference.setP(pc);
        } else {
            throw OpalException("Beam::update()",
                                "\"PC\" should be greater than 0.");
        }
    }

    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_BEAM");
}


void Beam::print(std::ostream& os) const {
    double charge = Attributes::getReal(itsAttr[CHARGE]);
    os << "* ************* B E A M ************************************************************ " << std::endl;
    os << "* BEAM        " << getOpalName() << '\n'
       << "* PARTICLE    " << Attributes::getString(itsAttr[PARTICLE]) << '\n'
       << "* REST MASS   " << Attributes::getReal(itsAttr[MASS]) << " [GeV]\n"
       << "* CHARGE      " << (charge > 0 ? '+' : '-') << "e * " << std::abs(charge) << " \n"
       << "* MOMENTUM    " << reference.getP() << " [eV/c]\n"
       << "* CURRENT     " << Attributes::getReal(itsAttr[BCURRENT]) << " [A]\n"
       << "* FREQUENCY   " << Attributes::getReal(itsAttr[BFREQ]) << " [MHz]\n"
       << "* NPART       " << Attributes::getReal(itsAttr[NPART]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}
