// Copyright (c) 2023, Chris Rogers
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

#include "Elements/OpalOutputPlane.h"
#include "Classic/AbsBeamline/OutputPlane.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

const std::string docstring =
    std::string("The \"OUTPUTPLANE\" element writes out position at which")+
    std::string("trajectories cross a given plane.");

OpalOutputPlane::OpalOutputPlane():
    OpalElement(SIZE, "OUTPUTPLANE", docstring.c_str()) {

    itsAttr[CENTRE] = Attributes::makeRealArray
                      ("CENTRE", "3-vector position of the plane centre [m]");
    itsAttr[NORMAL] = Attributes::makeRealArray
                      ("NORMAL", "3-vector normal to the plane");
    itsAttr[REFERENCE_ALIGNMENT_PARTICLE] = Attributes::makeReal
                      ("REFERENCE_ALIGNMENT_PARTICLE",
                       "Set to a particle number (usually 0, the reference particle). "
                       "The first time that the particle crosses the reference plane, then "
                       "the plane will be moved to centre on that particle and point in S of the particle.", -1);
    itsAttr[TOLERANCE] = Attributes::makeReal
                      ("TOLERANCE", "Tolerance on position of track intercept [m]", 1e-6);
    itsAttr[WIDTH] = Attributes::makeReal
                      ("WIDTH", "Full width of the output plane [m], defined in the lab coordinate system (*not* the output plane coordinate system).");
    itsAttr[HEIGHT] = Attributes::makeReal
                      ("HEIGHT", "Full height of the output plane [m], defined in the lab coordinate system (*not* the output plane coordinate system)");
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Maximum distance from centre of plane for crossings [m].");
    itsAttr[ALGORITHM] = Attributes::makePredefinedString
                      ("ALGORITHM",
                       "The algorithm used to step from the track point to the plane", {"INTERPOLATION", "RK4"}, "RK4");
    itsAttr[XSTART] = Attributes::makeReal
                      ("XSTART", "Define a plane with horizontal extent [m] from XSTART to XEND and vertical extent from YSTART to YEND", 0.0);
    itsAttr[YSTART] = Attributes::makeReal
                      ("YSTART", "Define a plane with horizontal extent [m] from XSTART to XEND and vertical extent from YSTART to YEND", 1.0);
    itsAttr[XEND] = Attributes::makeReal
                      ("XEND", "Define a plane with horizontal extent [m] from XSTART to XEND and vertical extent from YSTART to YEND", 0.0);
    itsAttr[YEND] = Attributes::makeReal
                      ("YEND", "Define a plane with horizontal extent [m] from XSTART to XEND and vertical extent from YSTART to YEND", 0.0);

    itsAttr[PLACEMENT_STYLE] = Attributes::makePredefinedString("PLACEMENT_STYLE",
                     "Set to PROBE to define the plane using XSTART, XEND, YSTART, YEND or CENTRE_NORMAL to define the plane using centre and normal", {"CENTRE_NORMAL", "PROBE"}, "PROBE");
    itsAttr[VERBOSE] = Attributes::makeReal("VERBOSE",
                     "Set to 0 for minimal output up to 4 to output diagnostics on every track step. Output is sent to OPAL", 0);
    registerOwnership();
    setElement(new OutputPlane("OUTPUTPLANE"));
}


OpalOutputPlane::OpalOutputPlane(const std::string &name, OpalOutputPlane *parent):
    OpalElement(name, parent) {
    OutputPlane* plane = dynamic_cast<OutputPlane*>(parent->getElement());
    setElement(new OutputPlane(*plane));
}


OpalOutputPlane::~OpalOutputPlane() {
}


OpalOutputPlane *OpalOutputPlane::clone(const std::string &name) {
    return new OpalOutputPlane(name, this);
}

void OpalOutputPlane::update() {
    OpalElement::update();
    std::string placementStyle = Attributes::getString(itsAttr[PLACEMENT_STYLE]);
    OutputPlane *output = dynamic_cast<OutputPlane *>(getElement());
    Vector_t centre;
    Vector_t normal;
    double tolerance = Attributes::getReal(itsAttr[TOLERANCE]);
    output->setTolerance(tolerance);

    if (placementStyle == "CENTRE_NORMAL") {
        std::vector<double> centreVec = Attributes::getRealArray(itsAttr[CENTRE]);
        if (centreVec.size() != 3) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane centre should be a 3-vector");
        }
        centre = Vector_t(centreVec[0], centreVec[1], centreVec[2]);
        std::vector<double> normalVec = Attributes::getRealArray(itsAttr[NORMAL]);
        if (normalVec.size() != 3) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane normal should be a 3-vector");
        }
        normal = Vector_t(normalVec[0], normalVec[1], normalVec[2]);
        if (euclidean_norm(normal) < tolerance) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane normal was not defined or almost zero");
        }
    } else if (placementStyle == "PROBE") {
        centre = Vector_t(
            Attributes::getReal(itsAttr[XSTART])+Attributes::getReal(itsAttr[XEND]),
            Attributes::getReal(itsAttr[YSTART])+Attributes::getReal(itsAttr[YEND]),
            0.0
        );
        centre *= 0.5; // average
        normal = Vector_t(
            Attributes::getReal(itsAttr[YSTART])-Attributes::getReal(itsAttr[YEND]),
            Attributes::getReal(itsAttr[XEND])-Attributes::getReal(itsAttr[XSTART]),
            0.0
        );
        double width = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
        if (width < tolerance) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane had very small width or size was not defined");
        }
        output->setHorizontalExtent(width/2.0);
    }
    output->setCentre(centre);
    output->setNormal(normal);


    if (itsAttr[WIDTH]) {
        double width = Attributes::getReal(itsAttr[WIDTH]);
        if (width < 0) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane had negative width");
        }
        output->setHorizontalExtent(width/2.0);
    }
    if (itsAttr[HEIGHT]) {
        double height = Attributes::getReal(itsAttr[HEIGHT]);
        if (height < 0) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane had negative height");
        }
        output->setVerticalExtent(height/2.0);
    }
    if (itsAttr[RADIUS]) {
        double radius = Attributes::getReal(itsAttr[RADIUS]);
        if (radius < 0) {
            throw OpalException("OpalOutputPlane::update()",
                "OutputPlane had negative radius");
        }
        output->setRadialExtent(radius);
    }

    if (itsAttr[REFERENCE_ALIGNMENT_PARTICLE]) {
        double ref = Attributes::getReal(itsAttr[REFERENCE_ALIGNMENT_PARTICLE]);
        int refAlign = std::floor(ref + 0.5);
        if (refAlign >= 0) { // disabled if negative
            output->setRecentre(refAlign);
        }
    } else {
        output->setRecentre(-1);
    }
    std::string algorithm = Attributes::getString(itsAttr[ALGORITHM]);
    if (algorithm == "RK4") {
        output->setAlgorithm(OutputPlane::algorithm::RK4STEP);
    } else if (algorithm == "INTERPOLATION") {
        output->setAlgorithm(OutputPlane::algorithm::INTERPOLATION);
    }
    output->setVerboseLevel(Attributes::getReal(itsAttr[VERBOSE]));
    output->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(output);
}