//
// Class RingSection
//   Defines the component placement handler in ring geometry.
//
// Copyright (c) 2012 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "Utilities/RingSection.h"

#include "AbsBeamline/Offset.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"


RingSection::RingSection():
    component_m(nullptr),
    componentPosition_m(0.), componentOrientation_m(0.),
    startPosition_m(0.), startOrientation_m(0.),
    endPosition_m(0.), endOrientation_m(0.) {
}

RingSection::RingSection(const RingSection& rhs):
    component_m(nullptr),
    componentPosition_m(0.), componentOrientation_m(0.),
    startPosition_m(0.), startOrientation_m(0.),
    endPosition_m(0.), endOrientation_m(0.) {
    *this = rhs;
}

RingSection::~RingSection() {
    //if (component_m != nullptr)
    //    delete component_m;
}

// Assignment operator
RingSection& RingSection::operator=(const RingSection& rhs) {
    if (&rhs != this) {
        component_m = dynamic_cast<Component*>(rhs.component_m->clone());
        if (component_m == nullptr) {
            throw GeneralClassicException("RingSection::operator=",
                                          "Failed to copy RingSection");
        }
        componentPosition_m = rhs.componentPosition_m;
        componentOrientation_m = rhs.componentOrientation_m;
        startPosition_m = rhs.startPosition_m;
        startOrientation_m = rhs.startOrientation_m;
        endPosition_m = rhs.endPosition_m;
        endOrientation_m = rhs.endOrientation_m;
    }
    return *this;
}

bool RingSection::isOnOrPastStartPlane(const Vector_t& pos) const {
    Vector_t posTransformed = pos - startPosition_m;
    // check that pos-startPosition_m is in front of startOrientation_m
    double normProd = posTransformed(0) * startOrientation_m(0) +
                      posTransformed(1) * startOrientation_m(1) +
                      posTransformed(2) * startOrientation_m(2);
    // check that pos and startPosition_m are on the same side of the ring
    double posProd = pos(0) * startPosition_m(0) +
                     pos(1) * startPosition_m(1) +
                     pos(2) * startPosition_m(2);
    return normProd >= 0. && posProd >= 0.;
}

bool RingSection::isPastEndPlane(const Vector_t& pos) const {
    Vector_t posTransformed = pos - endPosition_m;
    double normProd = posTransformed(0) * endOrientation_m(0) +
                      posTransformed(1) * endOrientation_m(1)+
                      posTransformed(2) * endOrientation_m(2);
    // check that pos and endPosition_m are on the same side of the ring
    double posProd = pos(0) * endPosition_m(0) +
                     pos(1) * endPosition_m(1) +
                     pos(2) * endPosition_m(2);
    return normProd > 0. && posProd > 0.;
}

bool RingSection::getFieldValue(const Vector_t& pos,
                                const Vector_t& /*centroid*/, const double& t,
                                Vector_t& E, Vector_t& B) const {
    // transform position into local coordinate system
    Vector_t pos_local = pos-componentPosition_m;
    rotate(pos_local);
    rotateToTCoordinates(pos_local);
    bool outOfBounds = component_m->apply(pos_local, Vector_t(0.0), t, E, B);
    if (outOfBounds) {
        return true;
    }
    // rotate fields back to global coordinate system
    rotateToCyclCoordinates(E);
    rotateToCyclCoordinates(B);
    rotate_back(E);
    rotate_back(B);
    return false;
}

void RingSection::updateComponentOrientation() {
    sin2_m = std::sin(componentOrientation_m(2));
    cos2_m = std::cos(componentOrientation_m(2));
}

std::vector<Vector_t> RingSection::getVirtualBoundingBox() const {
    Vector_t startParallel({getStartNormal()(1), -getStartNormal()(0), 0});
    Vector_t endParallel({getEndNormal()(1), -getEndNormal()(0), 0});
    normalise(startParallel);
    normalise(endParallel);
    double startRadius = 0.99 * std::sqrt(getStartPosition()(0) * getStartPosition()(0) +
                                          getStartPosition()(1) * getStartPosition()(1));
    double endRadius = 0.99 * std::sqrt(getEndPosition()(0) * getEndPosition()(0) +
                                        getEndPosition()(1) * getEndPosition()(1));
    std::vector<Vector_t> bb;
    bb.push_back(getStartPosition() - startParallel * startRadius);
    bb.push_back(getStartPosition() + startParallel * startRadius);
    bb.push_back(getEndPosition() - endParallel * endRadius);
    bb.push_back(getEndPosition() + endParallel * endRadius);
    return bb;
}

//    double phi = atan2(r(1), r(0))+Physics::pi;
bool RingSection::doesOverlap(double phiStart, double phiEnd) const {
    RingSection phiVirtualORS;
    phiVirtualORS.setStartPosition(Vector_t({std::sin(phiStart),
                                            std::cos(phiStart),
                                            0.}));
    phiVirtualORS.setStartNormal(Vector_t({std::cos(phiStart),
                                          -std::sin(phiStart),
                                          0.}));
    phiVirtualORS.setEndPosition(Vector_t({std::sin(phiEnd),
                                          std::cos(phiEnd),
                                          0.}));
    phiVirtualORS.setEndNormal(Vector_t({std::cos(phiEnd),
                                        -std::sin(phiEnd),
                                        0.}));
    std::vector<Vector_t> virtualBB = getVirtualBoundingBox();
    // at least one of the bounding box coordinates is in the defined sector
    // std::cerr << "RingSection::doesOverlap " << phiStart << " " << phiEnd << " "
    //          << phiVirtualORS.getStartPosition() << " " << phiVirtualORS.getStartNormal() << " "
    //          << phiVirtualORS.getEndPosition() << " " << phiVirtualORS.getEndNormal() << " " << std::endl
    //          << "    Component " << this << " " << getStartPosition() << " " << getStartNormal() << " "
    //          << getEndPosition() << " " << getEndNormal() << " " << std::endl;
    for (size_t i = 0; i < virtualBB.size(); ++i) {
        // std::cerr << "    VBB " << i << " " << virtualBB[i] << std::endl;
        if (phiVirtualORS.isOnOrPastStartPlane(virtualBB[i]) &&
            !phiVirtualORS.isPastEndPlane(virtualBB[i]))
            return true;
    }
    // the bounding box coordinates span the defined sector and the sector
    // sits inside the bb
    bool hasBefore = false; // some elements in bb are before phiVirtualORS
    bool hasAfter = false; // some elements in bb are after phiVirtualORS
    for (size_t i = 0; i < virtualBB.size(); ++i) {
        hasBefore = hasBefore ||
                    !phiVirtualORS.isOnOrPastStartPlane(virtualBB[i]);
        hasAfter = hasAfter ||
                   phiVirtualORS.isPastEndPlane(virtualBB[i]);
        // std::cerr << "    " << hasBefore << " " << hasAfter << std::endl;
    }
    if (hasBefore && hasAfter)
        return true;
    // std::cerr << "DOES NOT OVERLAP" << std::endl;
    return false;
}


void RingSection::rotate(Vector_t& vector) const {
    const Vector_t temp(vector);
    vector(0) = +cos2_m * temp(0) + sin2_m * temp(1);
    vector(1) = +sin2_m * temp(0) - cos2_m * temp(1);
}

void RingSection::rotate_back(Vector_t& vector) const {
    const Vector_t temp(vector);
    vector(0) = +cos2_m * temp(0) + sin2_m * temp(1);
    vector(1) = +sin2_m * temp(0) - cos2_m * temp(1);
}

void RingSection::handleOffset() {
    Offset* offsetCast = dynamic_cast<Offset*>(component_m);
    if (offsetCast == nullptr) {
        return; // nothing to do, it wasn't an offset at all
    }
    offsetCast->updateGeometry(startPosition_m, startOrientation_m);
}

