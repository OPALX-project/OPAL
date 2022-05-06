//
// Class OpalRBend
//   The parent class of all OPAL bending magnets.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalRBend3D.h"
#include "Attributes/Attributes.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"
#include "AbsBeamline/RBend3D.h"
#include "Utilities/OpalException.h"

#include <iostream>


OpalRBend3D::OpalRBend3D():
    OpalElement(SIZE, "RBEND3D", "The \"RBEND3D\" element defines an RBEND with 3D field maps"),
    owk_m(0),
    parmatint_m(nullptr) {
    itsAttr[ANGLE]  = Attributes::makeReal
                      ("ANGLE", "Upright dipole coefficient in m^(-1)");
    itsAttr[K0]     = Attributes::makeReal
                      ("K0", "Normal dipole coefficient in m^(-1)");
    itsAttr[K0S]    = Attributes::makeReal
                      ("K0S", "Skew dipole coefficient in m^(-1)");
    itsAttr[E1]     = Attributes::makeReal
                      ("E1", "Entry pole face angle in rad", 0.0);
    itsAttr[FMAPFN] = Attributes::makeString
                      ("FMAPFN", "Filename for the fieldmap");
    itsAttr[GAP]    = Attributes::makeReal
                      ("GAP", "Full gap height of the magnet (m)", 0.0);
    itsAttr[HAPERT] = Attributes::makeReal
                      ("HAPERT", "Bend plane magnet aperture (m)", 0.0);
    itsAttr[DESIGNENERGY] = Attributes::makeReal
                            ("DESIGNENERGY", "the mean energy of the particles in MeV");

    registerOwnership();

    setElement(new RBend3D("RBEND3D"));
}

OpalRBend3D::OpalRBend3D(const std::string& name, OpalRBend3D* parent):
    OpalElement(name, parent),
    owk_m(0),
    parmatint_m(nullptr)
{
    setElement(new RBend3D(name));
}

OpalRBend3D::~OpalRBend3D() {
    delete owk_m;
    delete parmatint_m;
}

OpalRBend3D* OpalRBend3D::clone(const std::string& name) {
    return new OpalRBend3D(name, this);
}


void OpalRBend3D::update() {
    OpalElement::update();

    // Define geometry.
    RBend3D* bend =
        dynamic_cast<RBend3D*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle  = Attributes::getReal(itsAttr[ANGLE]);
    double e1     = Attributes::getReal(itsAttr[E1]);
    double k0 =
        itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
        length ? 2 * sin(angle / 2) / length : angle;
    double k0s = itsAttr[K0S] ? Attributes::getReal(itsAttr[K0S]) : 0.0;

    // Set field amplitude or bend angle.
    if (itsAttr[ANGLE]) {
        if (bend->isPositioned() && angle < 0.0) {
            e1 = -e1;
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
            bend->releasePosition();
            bend->setCSTrafoGlobal2Local(CoordinateSystemTrafo(g2l.getOrigin(),
                                                               rotAboutZ * g2l.getRotation()));
            bend->fixPosition();
        }
        bend->setBendAngle(angle);
    } else {
        bend->setFieldAmplitude(k0, k0s);
    }

    if (itsAttr[FMAPFN]) {
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    } else if (bend->getName() != "RBEND3D") {
        ERRORMSG(bend->getName() << ": No filename for a field map given." << endl);
        throw OpalException("OpalRBend3D::update", bend->getName() + ": No filename for field map given");
    }

    bend->setEntranceAngle(e1);

    // Energy in eV.
    if (itsAttr[DESIGNENERGY] && Attributes::getReal(itsAttr[DESIGNENERGY]) != 0.0) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    } else if (bend->getName() != "RBEND3D") {
        throw OpalException("OpalRBend3D::update",
                            "RBend3D requires non-zero DESIGNENERGY");
    }

    bend->setFullGap(Attributes::getReal(itsAttr[GAP]));

    if (itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        bend->setAperture(ApertureType::RECTANGULAR, std::vector<double>({hapert, hapert, 1.0}));
    }

    if (itsAttr[LENGTH]) {
        bend->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    } else
        bend->setElementLength(0.0);

    if (itsAttr[WAKEF] && itsAttr[DESIGNENERGY] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*bend);
        bend->setWake(owk_m->wf_m);
    }

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*bend);
        bend->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}

void OpalRBend3D::print(std::ostream& os) const {
    OpalElement::print(os);
}
