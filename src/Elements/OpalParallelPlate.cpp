//
// Class OpalParallelPlate
//   The ParallelPlate element.
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
#include "Elements/OpalParallelPlate.h"
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/ParallelPlateRep.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

OpalParallelPlate::OpalParallelPlate():
    OpalElement(SIZE, "PARALLELPLATE",
                "The \"PARALLELPLATE\" element defines an  cavity."),
    obgeo_m(NULL)
{
    itsAttr[VOLT] = Attributes::makeReal
                    ("VOLT", " voltage in MV");
    itsAttr[FREQ] = Attributes::makeReal
                    ("FREQ", " frequency in MHz");
    itsAttr[LAG] = Attributes::makeReal
                   ("LAG", "Phase lag (rad), !!!! was before in multiples of (2*pi) !!!!");

    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for ParallelPlate");

    itsAttr[PLENGTH] = Attributes::makeReal
                       ("PLENGTH", " Gap length in Meter");


    registerRealAttribute("VOLT");
    registerRealAttribute("FREQ");
    registerRealAttribute("LAG");
    registerStringAttribute("GEOMETRY");
    registerRealAttribute("PLENGTH");

    registerOwnership();

    setElement(new ParallelPlateRep("ParallelPlate"));
}


OpalParallelPlate::OpalParallelPlate(const std::string &name, OpalParallelPlate *parent):
    OpalElement(name, parent),
    obgeo_m(NULL)
{
    setElement(new ParallelPlateRep(name));
}


OpalParallelPlate::~OpalParallelPlate() {

}


OpalParallelPlate *OpalParallelPlate::clone(const std::string &name) {
    return new OpalParallelPlate(name, this);
}


void OpalParallelPlate::fillRegisteredAttributes(const ElementBase &base) {
    OpalElement::fillRegisteredAttributes(base);

    const ParallelPlateRep *pplate =
        dynamic_cast<const ParallelPlateRep *>(&base);
    attributeRegistry["VOLT"]->setReal(pplate->getAmplitude());
    attributeRegistry["FREQ"]->setReal(pplate->getFrequency());
    attributeRegistry["LAG"]->setReal(pplate->getPhase());
    attributeRegistry["PLENGTH"]->setReal(pplate->getElementLength());
}


void OpalParallelPlate::update() {
    OpalElement::update();

    ParallelPlateRep *pplate =
        dynamic_cast<ParallelPlateRep *>(getElement());

    double vPeak  = Attributes::getReal(itsAttr[VOLT]);
    //    double phase  = two_pi * Attributes::getReal(itsAttr[LAG]);
    double phase  = Attributes::getReal(itsAttr[LAG]);
    double freq   = (1.0e6 * Physics::two_pi) * Attributes::getReal(itsAttr[FREQ]);
    double length = Attributes::getReal(itsAttr[PLENGTH]);

    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if(obgeo_m) {
            //obgeo_m->initialize();

            pplate->setBoundaryGeometry(obgeo_m);
        }
    }

    pplate->setAmplitude(1.0e6 * vPeak);
    pplate->setFrequency(freq);
    pplate->setPhase(phase);

    pplate->setElementLength(length);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(pplate);
}