//
// Class RBend
//   Interface for a rectangular bend magnet.
//
//   A rectangular bend magnet physically has a rectangular shape.
//
//   The standard rectangular magnet, for purposes of definitions, has a field in
//   the y direction. This produces a bend in the horizontal (x) plane. Bends in
//   other planes can be accomplished by rotating the magnet about the axes.
//
//   A positive bend angle is defined as one that bends a beam to the right when
//   looking down (in the negative y direction) so that the beam is bent in the
//   negative x direction. (This definition of a positive bend is the same whether
//   the charge is positive or negative.)
//
//   A zero degree entrance edge angle is parallel to the x direction in an x/y/s
//   coordinate system. A positive entrance edge angle is defined as one that
//   rotates the positive edge (in x) of the angle toward the positive s axis.
//
//   Since the magnet geometry is a fixed rectangle, the exit edge angle is
//   defined by the bend angle of the magnet and the entrance edge angle. In
//   general, the exit edge angle is equal to the bend angle minus the entrance
//   edge angle.
//
//   ------------------------------------------------------------------------
//
//   This class defines two interfaces:
//
//   1) Interface for rectangular magnets for OPAL-MAP.
//
//    Here we specify multipole components about the curved magnet trajectory.
//
//   2) Interface for rectangular magnets for OPAL-T.
//
//   Here we defined the magnet as a field map.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "AbsBeamline/RBend.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"

#include <iostream>
#include <fstream>
#include <cmath>

RBend::RBend():
    RBend("")
{}

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

bool RBend::findChordLength(double &chordLength) {

    /*
     * Find bend chord length. If this was not set by the user using the
     * L (length) attribute, infer it from the field map.
     */
    const double angle = getBendAngle();
    if (std::abs(angle) > 0.0) {
        double E1 = std::copysign(1.0, angle) * getEntranceAngle();
        chordLength = 2 * getLength() * std::sin(0.5 * std::abs(angle)) /
            (std::sin(E1) + std::sin(std::abs(angle) - E1));
    } else {
        double refCharge = RefPartBunch_m->getQ();
        double amplitude = (fieldAmplitudeY_m != 0.0) ? fieldAmplitudeY_m : fieldAmplitudeX_m;
        double fieldAmplitude = std::copysign(1.0, refCharge) * std::abs(amplitude);
        double designRadius = calcDesignRadius(fieldAmplitude);
        chordLength = std::sin(0.5 * (std::asin(getLength() / designRadius - std::sin(getEntranceAngle())) + getEntranceAngle())) * 2 * designRadius;
    }

    return true;
}