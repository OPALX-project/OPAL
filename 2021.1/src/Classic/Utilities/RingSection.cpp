/*
 *  Copyright (c) 2012-2014, Chris Rogers
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

#include "Utilities/RingSection.h"

#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

RingSection::RingSection()
  : component_m(NULL),
    componentPosition_m(0.), componentOrientation_m(0.),
    startPosition_m(0.), startOrientation_m(0.),
    endPosition_m(0.), endOrientation_m(0.) {
}

RingSection::RingSection(const RingSection& rhs)
  : component_m(NULL),
    componentPosition_m(0.), componentOrientation_m(0.),
    startPosition_m(0.), startOrientation_m(0.),
    endPosition_m(0.), endOrientation_m(0.) {
    *this = rhs;
}

RingSection::~RingSection() {
    //if (component_m != NULL)
    //    delete component_m;
}

// Assignment operator
RingSection& RingSection::operator=(const RingSection& rhs) {
    if (&rhs != this) {
        component_m = dynamic_cast<Component*>(rhs.component_m->clone());
        if (component_m == NULL)
            throw GeneralClassicException("RingSection::operator=",
                                "Failed to copy RingSection");
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
    Vector_t posTransformed = pos-startPosition_m;
    // check that pos-startPosition_m is in front of startOrientation_m
    double normProd = posTransformed(0)*startOrientation_m(0)+
                      posTransformed(1)*startOrientation_m(1)+
                      posTransformed(2)*startOrientation_m(2);
    // check that pos and startPosition_m are on the same side of the ring
    double posProd = pos(0)*startPosition_m(0)+
                     pos(1)*startPosition_m(1)+
                     pos(2)*startPosition_m(2);
    return normProd >= 0. && posProd >= 0.;
}

bool RingSection::isPastEndPlane(const Vector_t& pos) const {
    Vector_t posTransformed = pos-endPosition_m;
    double normProd = posTransformed(0)*endOrientation_m(0)+
                      posTransformed(1)*endOrientation_m(1)+
                      posTransformed(2)*endOrientation_m(2);
    // check that pos and startPosition_m are on the same side of the ring
    double posProd = pos(0)*endPosition_m(0)+
                     pos(1)*endPosition_m(1)+
                     pos(2)*endPosition_m(2);
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
    sin2_m = sin(componentOrientation_m(2));
    cos2_m = cos(componentOrientation_m(2));
}

std::vector<Vector_t> RingSection::getVirtualBoundingBox() const {
    Vector_t startParallel(getStartNormal()(1), -getStartNormal()(0), 0);
    Vector_t endParallel(getEndNormal()(1), -getEndNormal()(0), 0);
    normalise(startParallel);
    normalise(endParallel);
    double startRadius = 0.99*sqrt(getStartPosition()(0)*getStartPosition()(0)+
                                   getStartPosition()(1)*getStartPosition()(1));
    double endRadius = 0.99*sqrt(getEndPosition()(0)*getEndPosition()(0)+
                                 getEndPosition()(1)*getEndPosition()(1));
    std::vector<Vector_t> bb;
    bb.push_back(getStartPosition()-startParallel*startRadius);
    bb.push_back(getStartPosition()+startParallel*startRadius);
    bb.push_back(getEndPosition()-endParallel*endRadius);
    bb.push_back(getEndPosition()+endParallel*endRadius);
    return bb;
}

//    double phi = atan2(r(1), r(0))+Physics::pi;
bool RingSection::doesOverlap(double phiStart, double phiEnd) const {
    RingSection phiVirtualORS;
    // phiStart -= Physics::pi;
    // phiEnd -= Physics::pi;
    phiVirtualORS.setStartPosition(Vector_t(sin(phiStart),
                                            cos(phiStart),
                                            0.));
    phiVirtualORS.setStartNormal(Vector_t(cos(phiStart),
                                          -sin(phiStart),
                                          0.));
    phiVirtualORS.setEndPosition(Vector_t(sin(phiEnd),
                                          cos(phiEnd),
                                          0.));
    phiVirtualORS.setEndNormal(Vector_t(cos(phiEnd),
                                        -sin(phiEnd),
                                        0.));
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
    vector(1) = -sin2_m * temp(0) + cos2_m * temp(1);
}

void RingSection::rotate_back(Vector_t& vector) const {
    const Vector_t temp(vector);
    vector(0) = +cos2_m * temp(0) - sin2_m * temp(1);
    vector(1) = +sin2_m * temp(0) + cos2_m * temp(1);
}
