// ------------------------------------------------------------------------
// $RCSfile: Component.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Component
//   An abstract base class which defines the common interface for all
//   CLASSIC components, i.e. beamline members which are not themselves
//   beamlines.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//c
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Utilities/LogicalError.h"
#include "Algorithms/PartBunchBase.h"

// Class Component
// ------------------------------------------------------------------------
//   Represents an arbitrary component in an accelerator.  A component is
//   the basic element in the accelerator model, and can be thought of as
//   acting as a leaf in the Composite pattern.  A Component is associated
//   with an electromagnetic field.
// 2017-03-20 (Rogers) set default aperture to something huge; else we get a
//         segmentation fault by default from NULL dereference during tracking

const std::vector<double> Component::defaultAperture_m =
                                std::vector<double>({1e6, 1e6, 1.0});


Component::Component():
    Component("")
{}


Component::Component(const Component &right):
    ElementBase(right),
    exit_face_slope_m(right.exit_face_slope_m),
    RefPartBunch_m(right.RefPartBunch_m),
    online_m(right.online_m) {
}


Component::Component(const std::string &name):
    ElementBase(name),
    exit_face_slope_m(0.0),
    RefPartBunch_m(NULL),
    online_m(false) {
    setAperture(ElementBase::ELLIPTICAL, defaultAperture_m);

}


Component::~Component()
{}


const ElementBase &Component::getDesign() const {
    return *this;
}

void Component::trackBunch(PartBunchBase<double, 3> *, const PartData &, bool, bool) const {
    throw LogicalError("Component::trackBunch()",
                       "Called for component \"" + getName() + "\".");
}


void Component::
trackMap(FVps<double, 6> &, const PartData &, bool, bool) const {
    throw LogicalError("Component::trackMap()",
                       "Called for component \"" + getName() + "\".");
}

void Component::goOnline(const double &) {
    online_m = true;
}

void Component::goOffline() {
    online_m = false;
}

bool Component::Online() {
    return online_m;
}

ElementBase::ElementType Component::getType() const {
    return ElementBase::ANY;
}

bool Component::apply(const size_t &i,
                      const double &/*t*/,
                      Vector_t &/*E*/,
                      Vector_t &/*B*/) {
    const Vector_t &R = RefPartBunch_m->R[i];
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}

bool Component::apply(const Vector_t &R,
                      const Vector_t &/*P*/,
                      const double &/*t*/,
                      Vector_t &/*E*/,
                      Vector_t &/*B*/) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}

bool Component::applyToReferenceParticle(const Vector_t &R,
                                         const Vector_t &/*P*/,
                                         const double &/*t*/,
                                         Vector_t &/*E*/,
                                         Vector_t &/*B*/) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R)) return true;
    }
    return false;
}
