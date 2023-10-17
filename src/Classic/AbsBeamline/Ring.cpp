//
// Class Ring
//   Defines the abstract interface for a ring type geometry.
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
#include "AbsBeamline/Ring.h"

#include "Utility/Inform.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Fields/EMField.h"
#include "Physics/Units.h"
#include "Structure/LossDataSink.h"


#include <limits>
#include <sstream>

// fairly generous tolerance here... maybe too generous? Maybe should be
// user parameter?
const double Ring::lengthTolerance_m = 1e-2;
const double Ring::angleTolerance_m = 1e-2;

extern Inform* gmsg;
extern Inform* gmsgALL;

Ring::Ring(const std::string& ring)
    : Component(ring), planarArcGeometry_m(1, 1),
      refPartBunch_m(nullptr),
      lossDS_m(nullptr),
      beamRInit_m(0.), beamPRInit_m(0.), beamPhiInit_m(0.),
      latticeRInit_m(0.), latticePhiInit_m(0.), latticeThetaInit_m(0.),
      isLocked_m(false), isClosed_m(true),
      symmetry_m(0), cyclHarm_m(0), rfFreq_m(0),
      phiStep_m(Physics::pi/100.),
      ringSections_m(),
      section_list_m() {
    setRefPartBunch(nullptr);
}

void Ring::accept(BeamlineVisitor& visitor) const {
    visitor.visitRing(*this);
}

Ring::Ring(const Ring& ring)
    : Component(ring.getName()),
      planarArcGeometry_m(ring.planarArcGeometry_m),
      lossDS_m(nullptr),
      beamRInit_m(ring.beamRInit_m),
      beamPRInit_m(ring.beamPRInit_m),
      beamPhiInit_m(ring.beamPhiInit_m),
      beamThetaInit_m(ring.beamThetaInit_m),
      latticeRInit_m(ring.latticeRInit_m),
      latticePhiInit_m(ring.latticePhiInit_m),
      latticeThetaInit_m(ring.latticeThetaInit_m),
      willDoAperture_m(ring.willDoAperture_m),
      minR2_m(ring.minR2_m),
      maxR2_m(ring.maxR2_m),
      isLocked_m(ring.isLocked_m),
      isClosed_m(ring.isClosed_m),
      symmetry_m(ring.symmetry_m),
      cyclHarm_m(ring.cyclHarm_m),
      phiStep_m(ring.phiStep_m),
      ringSections_m(),
      section_list_m(ring.section_list_m.size()) {
    setRefPartBunch(ring.refPartBunch_m);
    if (ring.lossDS_m != nullptr)
        throw GeneralClassicException("Ring::Ring(const Ring&)",
                                      "Can't copy construct LossDataSink so copy constructor fails");
    for (size_t i = 0; i < section_list_m.size(); ++i) {
        section_list_m[i] = new RingSection(*ring.section_list_m[i]);
    }
    buildRingSections();
}

Ring::~Ring() {
    for (size_t i = 0; i < section_list_m.size(); ++i)
        delete section_list_m[i];
}

bool Ring::apply(const size_t& id, const double& t,
                 Vector_t& E, Vector_t& B) {

    bool flagNeedUpdate = apply(refPartBunch_m->R[id], refPartBunch_m->P[id], t, E, B);

    if (flagNeedUpdate) {
        *gmsgALL << level4 << getName() << ": particle " << id
                << " at " << refPartBunch_m->R[id]
                << " m out of the field map boundary" << endl;

        lossDS_m->addParticle(OpalParticle(id,
                                           refPartBunch_m->R[id] , refPartBunch_m->P[id],
                                           t,
                                           refPartBunch_m->Q[id], refPartBunch_m->M[id]));

        refPartBunch_m->Bin[id] = -1;
    }

    return flagNeedUpdate;
}

bool Ring::apply(const Vector_t& R, const Vector_t& /*P*/,
                 const double& t, Vector_t& E, Vector_t& B) {
    B = Vector_t(0.0, 0.0, 0.0);
    E = Vector_t(0.0, 0.0, 0.0);

    std::vector<RingSection*> sections = getSectionsAt(R); // I think this doesn't actually use R -DW
    bool outOfBounds = true;
    // assume field maps don't set B, E to 0...
    if (willDoAperture_m) {
        double rSquared = R[0] * R[0] + R[1] * R[1];
        if (rSquared < minR2_m || rSquared > maxR2_m) {
            return true;
        }
    }
    for (size_t i = 0; i < sections.size(); ++i) {
        Vector_t B_temp(0.0, 0.0, 0.0);
        Vector_t E_temp(0.0, 0.0, 0.0);
        outOfBounds &= sections[i]->getFieldValue(R, refPartBunch_m->get_centroid(), t, E_temp, B_temp);
        B += (scale_m * B_temp);
        E += (scale_m * E_temp);
    }
    return outOfBounds;
}

void Ring::setLossDataSink(LossDataSink* sink) {
    // lossDS_m is a borrowed pointer - we never delete it
    lossDS_m = sink;
}

void Ring::getDimensions(double& /*zBegin*/, double& /*zEnd*/) const {
    throw GeneralClassicException("Ring::getDimensions",
                                  "Cannot get s-dimension of a ring");
}

void Ring::initialise(PartBunchBase<double, 3>* bunch) {
    online_m = true;
    setRefPartBunch(bunch);
    setLossDataSink(new LossDataSink(getName(), false));
}

void Ring::initialise(PartBunchBase<double, 3>* bunch,
                      double& /*startField*/, double& /*endField*/) {
    initialise(bunch);
}

void Ring::finalise() {
    lossDS_m->save();
    online_m = false;
    setLossDataSink(nullptr);
}

void Ring::setRefPartBunch(PartBunchBase<double, 3>* bunch) {
    RefPartBunch_m = bunch; // inherited from Component
    refPartBunch_m = bunch; // private data (obeys style guide)
}

std::vector<RingSection*> Ring::getSectionsAt(const Vector_t& /*r*/) {
    return section_list_m;
}

Rotation3D Ring::getRotationStartToEnd(Euclid3D delta) const {
    // rotation from start to end in local coordinates
    // Euclid3D/Rotation3D doesnt have a "getAngle" method so we use fairly
    // obscure technique to extract it.
    Vector3D rotationTest(1., 0., 0.);
    rotationTest = delta.getRotation() * rotationTest;
    double deltaAngle = std::atan2(rotationTest(2), rotationTest(0));
    Rotation3D elementRotation = Rotation3D::ZRotation(deltaAngle);
    return elementRotation;
}

void Ring::checkMidplane(Euclid3D delta) const {
    if (std::abs(delta.getVector()(2)) > lengthTolerance_m ||
        std::abs(delta.getRotation().getAxis()(0)) > angleTolerance_m ||
	    std::abs(delta.getRotation().getAxis()(1)) > angleTolerance_m) {
        throw GeneralClassicException("Ring::checkMidplane",
                                      std::string("Placement of elements out of the ")+
                                      "midplane is not supported by Ring");
    }
}

void Ring::rotateToCyclCoordinates(Euclid3D& delta) const {
    Vector3D v = delta.getVector();
    Vector3D r = delta.getRotation().getAxis();
    //    Euclid3D euclid(v(0), -v(2), v(1), r(0), -r(2), r(1));
    Euclid3D euclid(v(2), v(0), -v(1), r(2), r(0), -r(1));
    delta = euclid;
}


Vector_t Ring::getNextPosition() const {
    if (!section_list_m.empty()) {
        return section_list_m.back()->getEndPosition();
    } else {
        return getStartPosition();
    }
}

Vector_t Ring::getStartPosition() const {
    return Vector_t(latticeRInit_m * std::cos(latticePhiInit_m),
                    latticeRInit_m * std::sin(latticePhiInit_m),
                    0.);
}

Vector_t Ring::getNextNormal() const {
    if (!section_list_m.empty()) {
        return section_list_m.back()->getEndNormal();
    } else {
        return getStartNormal();
    }
}

Vector_t Ring::getStartNormal() const {
    return Vector_t(-std::sin(latticePhiInit_m+latticeThetaInit_m),
                    std::cos(latticePhiInit_m+latticeThetaInit_m),
                    0.);
}

void Ring::appendElement(const Component& element) {
    if (isLocked_m) {
        throw GeneralClassicException("Ring::appendElement",
                                      "Attempt to append element " + element.getName() +
                                      " when ring is locked");
    }

    RingSection* section = new RingSection();
    Vector_t startPos = getNextPosition();
    Vector_t startNorm = getNextNormal();

    section->setComponent(dynamic_cast<Component*>(element.clone()));
    section->setStartPosition(startPos);
    section->setStartNormal(startNorm);
    section->handleOffset();
    // delta is transform from start of bend to end with x, z as horizontal
    // I failed to get Rotation3D to work so use rotations written by hand.
    // Probably an error in my call to Rotation3D.
    Euclid3D delta = section->getComponent()->getGeometry().getTotalTransform();
    rotateToCyclCoordinates(delta);
    checkMidplane(delta);

    double placeF = std::atan2(startNorm(0), startNorm(1)); // angle between y axis and norm
    Vector_t endPos = Vector_t(+delta.getVector()(0) * std::sin(placeF) + delta.getVector()(1) * std::cos(placeF),
                               +delta.getVector()(0) * std::cos(placeF) - delta.getVector()(1) * std::sin(placeF),
                               0) + startPos;
    section->setEndPosition(endPos);

    double endF = delta.getRotation().getAxis()(2);//+
    //atan2(delta.getVector()(1), delta.getVector()(0));
    Vector_t endNorm = Vector_t(+startNorm(0) * std::cos(endF) - startNorm(1) * std::sin(endF),
                                +startNorm(0) * std::sin(endF) + startNorm(1) * std::cos(endF),
                                0);
    section->setEndNormal(endNorm);

    section->setComponentPosition(startPos);
    double orientation = std::atan2(startNorm(1), startNorm(0));
    section->setComponentOrientation(Vector_t(0, 0, orientation));

    section_list_m.push_back(section);
    *gmsg << "* Added " << element.getName() << " to Ring" << endl;
    *gmsg << "* Start position ("
        << section->getStartPosition()(0) << ", "
        << section->getStartPosition()(1) << ", "
        << section->getStartPosition()(2) << ") normal ("
        << section->getStartNormal()(0) << ", "
        << section->getStartNormal()(1) << ", "
        << section->getStartNormal()(2) << ")" << endl;
    *gmsg << "* End position   ("
        << section->getEndPosition()(0) << ", "
        << section->getEndPosition()(1) << ", "
        << section->getEndPosition()(2) << ") normal ("
        << section->getEndNormal()(0) << ", "
        << section->getEndNormal()(1) << ", "
        << section->getEndNormal()(2) << ")" << endl;
}

void Ring::lockRing() {
    if (isLocked_m) {
        throw GeneralClassicException("Ring::lockRing",
                                      "Attempt to lock ring when it is already locked");
    }
    // check for any elements at all
    size_t sectionListSize = section_list_m.size();
    if (sectionListSize == 0) {
        throw GeneralClassicException("Ring::lockRing",
                                      "Failed to find any elements in Ring");
    }
    // Apply symmetry properties; I think it is fastest to just duplicate
    // elements rather than try to do some angle manipulations when we do field
    // lookups because we do field lookups in Cartesian coordinates in general.
    *gmsg << "Applying symmetry to Ring" << endl;
    for (int i = 1; i < symmetry_m; ++i) {
        for (size_t j = 0; j < sectionListSize; ++j) {
            appendElement(*section_list_m[j]->getComponent());
        }
    }
    //resetAzimuths();
    // Check that the ring is closed within tolerance and close exactly
    if (isClosed_m)
        checkAndClose();
    buildRingSections();
    isLocked_m = true;
    for (size_t i = 0; i < section_list_m.size(); i++) {
    }
}

void Ring::resetAzimuths() {
    for (size_t i = 0; i < section_list_m.size(); ++i) {
        Vector_t startPos = section_list_m[i]->getEndPosition();
        Vector_t startDir = section_list_m[i]->getEndNormal();
        Vector_t endPos = section_list_m[i]->getEndPosition();
        Vector_t endDir = section_list_m[i]->getEndNormal();
        if (!section_list_m[i]->isOnOrPastStartPlane(endPos)) {
            section_list_m[i]->setEndPosition(startPos);
            section_list_m[i]->setEndNormal(startDir);
            section_list_m[i]->setStartPosition(endPos);
            section_list_m[i]->setStartNormal(endDir);
        }
    }
}

void Ring::checkAndClose() {
    Vector_t first_pos  = section_list_m[0]->getStartPosition();
    Vector_t first_norm = section_list_m[0]->getStartNormal();
    Vector_t last_pos   = section_list_m.back()->getEndPosition();
    Vector_t last_norm  = section_list_m.back()->getEndNormal();
    for (int i = 0; i < 3; ++i) {
        if (std::abs(first_pos(i) - last_pos(i)) > lengthTolerance_m ||
            std::abs(first_norm(i) - last_norm(i)) > angleTolerance_m)
            throw GeneralClassicException("Ring::lockRing",
                                          "Ring is not closed");
    }
    section_list_m.back()->setEndPosition(first_pos);
    section_list_m.back()->setEndNormal(first_norm);
}

void Ring::buildRingSections() {
    size_t nSections = Physics::two_pi/phiStep_m+1;
    ringSections_m = std::vector< std::vector<RingSection*> >(nSections);
    for (size_t i = 0; i < ringSections_m.size(); ++i) {
        double phi0 = i*phiStep_m;
        double phi1 = (i+1)*phiStep_m;
        for (size_t j = 0; j < section_list_m.size(); ++j) {
            if (section_list_m[j]->doesOverlap(phi0, phi1))
                ringSections_m[i].push_back(section_list_m[j]);
        }
    }
}

RingSection* Ring::getLastSectionPlaced() const {
    if (section_list_m.empty()) {
        return nullptr;
    }
    return section_list_m.back();
}

bool Ring::sectionCompare(RingSection const* const sec1,
                          RingSection const* const sec2) {
    return sec1 > sec2;
}

void Ring::setRingAperture(double minR, double maxR) {
    if (minR < 0 || maxR < 0) {
        throw GeneralClassicException("Ring::setRingAperture",
                                      "Could not parse negative or undefined aperture limit");
    }

    willDoAperture_m = true;
    minR2_m = minR * minR;
    maxR2_m = maxR * maxR;
}

RingSection* Ring::getSection(int i) const {
    if (i < 0 || i >= int(getNumberOfRingSections())) {
        std::stringstream err;
        err << "Attempt to get RingSection for element " << i
            << ". Should be in range 0 <= i < " << getNumberOfRingSections();
        throw GeneralClassicException("Ring::getSection(int)", err.str());
    }
    RingSection* sec = section_list_m[i];
    if (sec == nullptr) {
        throw GeneralClassicException("Ring::getSection",
                                      "Opal internal error - RingSection was null");
    }
    return sec;
}

size_t Ring::getNumberOfRingSections() const {
    return section_list_m.size();
}
