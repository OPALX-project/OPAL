/*
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "Elements/OpalScalingFFAMagnet.h"

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Attributes/Attributes.h"
#include "Physics/Units.h"

extern Inform *gmsg;

OpalScalingFFAMagnet::OpalScalingFFAMagnet() :
    OpalElement(SIZE, "SCALINGFFAMAGNET",
             "The \"ScalingFFAMagnet\" element defines a FFA scaling magnet with zero or non-zero spiral angle.") {
    itsAttr[B0] = Attributes::makeReal
                              ("B0", "The nominal dipole field of the magnet [T].");
    itsAttr[R0] = Attributes::makeReal("R0", "Radial scale [m].");
    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX",
      "The scaling magnet field index.");
    itsAttr[TAN_DELTA] = Attributes::makeReal("TAN_DELTA",
      "Tangent of the spiral angle; set to 0 to make a radial sector magnet.");
    itsAttr[MAX_Y_POWER] = Attributes::makeReal("MAX_Y_POWER",
      "The maximum power in y that will be considered in the field expansion.");
    itsAttr[END_FIELD_MODEL] = Attributes::makeString("END_FIELD_MODEL",
      "Names the end field model of the ring, giving dipole field along a line of "
      "constant radius. If blank, uses the default 'END_LENGTH' and 'CENTRE_LENGTH' "
      "parameters and a tanh model. If 'END_FIELD_MODEL' is not blank, Opal will seek "
      "an END_FIELD_MODEL corresponding to the name defined in this string.");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
                                          "The end length of the spiral FFA [m].");
    itsAttr[HEIGHT] = Attributes::makeReal("HEIGHT",
                                       "Full height of the magnet. Particles moving more than height/2. off the midplane (either above or below) are out of the aperture [m].");
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
                                       "The centre length of the spiral FFA [m].");
    itsAttr[RADIAL_NEG_EXTENT] = Attributes::makeReal("RADIAL_NEG_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0-RADIAL_NEG_EXTENT [m].");
    itsAttr[RADIAL_POS_EXTENT] = Attributes::makeReal("RADIAL_POS_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0+RADIAL_POS_EXTENT [m].");
    itsAttr[MAGNET_START] = Attributes::makeReal("MAGNET_START",
                                          "Determines the position of the central portion of the magnet field relative to the element start (default is 2*end_length). [m]");
    itsAttr[MAGNET_END] = Attributes::makeReal("MAGNET_END",
                                       "Offset to the end of the magnet, i.e. placement of the next element. Default is centre_length + 4*end_length.");
    itsAttr[AZIMUTHAL_EXTENT] = Attributes::makeReal("AZIMUTHAL_EXTENT",
                                       "The field will be assumed zero if particles are more than AZIMUTHAL_EXTENT from the magnet centre (psi=0). Default is CENTRE_LENGTH/2.+5.*END_LENGTH [m].");
    registerOwnership();

    ScalingFFAMagnet* magnet = new ScalingFFAMagnet("ScalingFFAMagnet");
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet);
}


OpalScalingFFAMagnet::OpalScalingFFAMagnet(const std::string &name,
                                             OpalScalingFFAMagnet *parent) :
    OpalElement(name, parent) {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(name);
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet);
}


OpalScalingFFAMagnet::~OpalScalingFFAMagnet() {
}


OpalScalingFFAMagnet *OpalScalingFFAMagnet::clone(const std::string &name) {
    return new OpalScalingFFAMagnet(name, this);
}

void OpalScalingFFAMagnet::setupDefaultEndField() {
    ScalingFFAMagnet *magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());
    // get centre length and end length in metres
    endfieldmodel::Tanh* endField = new endfieldmodel::Tanh();
    double end_length = Attributes::getReal(itsAttr[END_LENGTH]);
    double centre_length = Attributes::getReal(itsAttr[CENTRE_LENGTH])/2.;
    endField->setLambda(end_length);
    // x0 is the distance between B=0.5*B0 and B=B0 i.e. half the centre length
    endField->setX0(centre_length);
    std::shared_ptr<endfieldmodel::EndFieldModel> efm(endField);
    std::string endName = "__opal_internal__"+getOpalName();
    endfieldmodel::EndFieldModel::setEndFieldModel(endName, efm);
    magnet->setEndFieldName(endName);
}

void OpalScalingFFAMagnet::setupNamedEndField() {
    if (!itsAttr[END_FIELD_MODEL]) {
        return;
    }
    std::string name = Attributes::getString(itsAttr[END_FIELD_MODEL]);
    ScalingFFAMagnet *magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());
    magnet->setEndFieldName(name);
}

void OpalScalingFFAMagnet::update() {
    ScalingFFAMagnet *magnet = dynamic_cast<ScalingFFAMagnet*>(getElement());

    // use L = r0*theta; we define the magnet ito length for UI but ito angles
    // internally; and use m as external default unit and mm internally
    // get r0 in m
    double r0Abs = std::abs(Attributes::getReal(itsAttr[R0]));
    double r0Signed = Attributes::getReal(itsAttr[R0]);
    magnet->setR0(r0Signed*Units::m2mm);
    // get B0 in T
    magnet->setDipoleConstant(Attributes::getReal(itsAttr[B0])*Units::T2kG);

    // dimensionless quantities
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setTanDelta(Attributes::getReal(itsAttr[TAN_DELTA]));
    int maxOrder = floor(Attributes::getReal(itsAttr[MAX_Y_POWER]));
    magnet->setMaxOrder(maxOrder);

    if (itsAttr[END_FIELD_MODEL]) {
        setupNamedEndField();
    } else {
        setupDefaultEndField();    
    }
    // internally OpalScalingFFAMagnet uses radians, so we scale all lengths to
    // radians.
    magnet->getEndField()->rescale(1/r0Abs);

    // get rmin and rmax bounding box edge in mm
    if (!itsAttr[RADIAL_NEG_EXTENT]) {
            throw OpalException("OpalScalingFFAMagnet::update()",
                                "RADIAL_NEG_EXTENT needs to be defined");    
    }
    double rmin = r0Abs-Attributes::getReal(itsAttr[RADIAL_NEG_EXTENT]);
    if (!itsAttr[RADIAL_POS_EXTENT]) {
            throw OpalException("OpalScalingFFAMagnet::update()",
                                "RADIAL_POS_EXTENT needs to be defined");    
    }
    double rmax = r0Abs+Attributes::getReal(itsAttr[RADIAL_POS_EXTENT]);
    magnet->setRMin(rmin*Units::m2mm);
    magnet->setRMax(rmax*Units::m2mm);
    Vector_t centre(r0Signed*Units::m2mm, 0, 0);
    magnet->setCentre(centre);

    // we store maximum vertical displacement (which is half the height)
    double height = Attributes::getReal(itsAttr[HEIGHT])*Units::m2mm;
    magnet->setVerticalExtent(height/2.);
    // end of the magnet marks the point at which the next element starts
    if (itsAttr[MAGNET_END]) {
        if (Attributes::getReal(itsAttr[MAGNET_END]) < 0.0) {
            throw OpalException("OpalScalingFFAMagnet::update()",
                                "MAGNET_END must be > 0.0");
        }
        double phi_end = Attributes::getReal(itsAttr[MAGNET_END])/r0Abs;
        magnet->setPhiEnd(phi_end);
    } else {
        magnet->setPhiEnd(-1); // flag for setupEndField
    }

    // get start of the magnet element in radians
    // setPhiStart sets the position of the 0 point of the endFieldModel, which
    // is typically the magnet centre
    if (itsAttr[MAGNET_START]) {
        if (Attributes::getReal(itsAttr[MAGNET_START]) < 0.0) {
            throw OpalException("OpalScalingFFAMagnet::update()",
                                "MAGNET_START must be > 0.0");
        }
        double phi_start = Attributes::getReal(itsAttr[MAGNET_START])/r0Abs;
        magnet->setPhiStart(phi_start);
    } else {
        magnet->setPhiStart(-1); // flag for setupEndField
    }
    // get azimuthal extent in radians; this is just the bounding box
    if (itsAttr[AZIMUTHAL_EXTENT]) {
        if (Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT]) < 0.0) {
            throw OpalException("OpalScalingFFAMagnet::update()",
                                "AZIMUTHAL_EXTENT must be > 0.0");
        }
        magnet->setAzimuthalExtent(Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT])/r0Abs);
    } else {
        magnet->setAzimuthalExtent(-1); // flag for setupEndField
    }
    magnet->initialise();
    setElement(magnet);

}