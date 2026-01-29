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
#include "Physics/Physics.h"
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

const double Beam::energy_scale = 1.0e9;

Beam::Beam():
    Definition(SIZE, "BEAM",
               "The \"BEAM\" statement defines data for the particles "
               "in a beam."),
    reference(1.0, Physics::m_p *energy_scale, 1.0 * energy_scale) {

    itsAttr[PARTICLE] = Attributes::makePredefinedString
                        ("PARTICLE", "Name of particle to be used",
                         {"ELECTRON",
                          "PROTON",
                          "POSITRON",
                          "ANTIPROTON",
                          "CARBON",
                          "HMINUS",
                          "URANIUM",
                          "MUON",
                          "DEUTERON",
                          "XENON",
                          "H2P",
                          "ALPHA"});

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
                            "The energy hasn't been set. Set either GAMMA, ENERGY or PC.");
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
    return (size_t)Attributes::getReal(itsAttr[NPART]);
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
        / (getFrequency() * 1.0e6)
        / getNumberOfParticles();
}

double Beam::getMassPerParticle() const {
    return getMass() * getChargePerParticle() / (getCharge() * Physics::q_e);
}

void Beam::update() {
    // Find the particle name.
    if (itsAttr[PARTICLE]) {
        static const char *names[] = {
            "ELECTRON",
            "PROTON",
            "POSITRON",
            "ANTIPROTON",
            "CARBON",
            "HMINUS",
            "URANIUM",
            "MUON",
            "DEUTERON",
            "XENON",
            "H2P",
            "ALPHA"
        };

        static const double masses[] = {
            Physics::m_e,
            Physics::m_p,
            Physics::m_e,
            Physics::m_p,
            Physics::m_c,
            Physics::m_hm,
            Physics::m_u,
            Physics::m_mu,
            Physics::m_d,
            Physics::m_xe,
            Physics::m_h2p,
            Physics::m_alpha
        };

        static const double charges[] = {
            -1.0, 1.0, 1.0, -1.0, 12.0, -1.0, 35.0, -1.0, 1.0, 20.0, 1.0, 2.0
        };
        const unsigned int numParticleNames = std::end(names) - std::begin(names);

        std::string pName  = Attributes::getString(itsAttr[PARTICLE]);
        for (unsigned int i = 0; i < numParticleNames; ++ i) {
            if (pName == names[i]) {
                Attributes::setReal(itsAttr[MASS], masses[i]);
                Attributes::setReal(itsAttr[CHARGE], charges[i]);
                break;
            }
        }
    }

    // Set up particle reference; convert all to eV for CLASSIC.
    double mass =
        (itsAttr[MASS] ? Attributes::getReal(itsAttr[MASS]) : Physics::m_p) * energy_scale;
    double charge = itsAttr[CHARGE] ? Attributes::getReal(itsAttr[CHARGE]) : 1.0;
    reference = PartData(charge, mass, 1.0);

    // Checks
    if (itsAttr[GAMMA]) {
        double gamma = Attributes::getReal(itsAttr[GAMMA]);
        if (gamma > 1.0) {
            reference.setGamma(gamma);
        } else {
            throw OpalException("Beam::update()",
                                "\"GAMMA\" should be greater than 1.");
        }
    } else if (itsAttr[ENERGY]) {
        double energy = Attributes::getReal(itsAttr[ENERGY]) * energy_scale;
        if (energy > reference.getM()) {
            reference.setE(energy);
        } else {
            throw OpalException("Beam::update()",
                                "\"ENERGY\" should be greater than \"MASS\".");
        }
    } else if (itsAttr[PC]) {
        double pc = Attributes::getReal(itsAttr[PC]) * energy_scale;
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


//ff
double Beam::getGamma() const { //obtain value for gamma
    return Attributes::getReal(itsAttr[GAMMA]);
}

//ff
double Beam::getPC() const { //obtain value for PC
    return Attributes::getReal(itsAttr[PC]);
}


void Beam::print(std::ostream& os) const {
    double charge = Attributes::getReal(itsAttr[CHARGE]);
    os << "* ************* B E A M ************************************************************ " << std::endl;
    os << "* BEAM        " << getOpalName() << '\n'
       << "* PARTICLE    " << Attributes::getString(itsAttr[PARTICLE]) << '\n'
       << "* CURRENT     " << Attributes::getReal(itsAttr[BCURRENT]) << " A\n"
       << "* FREQUENCY   " << Attributes::getReal(itsAttr[BFREQ]) << " MHz\n"
       << "* CHARGE      " << (charge > 0 ? '+' : '-') << "e * " << std::abs(charge) << " \n"
       << "* REST MASS   " << Attributes::getReal(itsAttr[MASS]) << " GeV\n"
       << "* MOMENTUM    " << Attributes::getReal(itsAttr[PC])   << '\n'
       << "* NPART       " << Attributes::getReal(itsAttr[NPART])   << '\n';
    os << "* ********************************************************************************** " << std::endl;
}