//
// Class ElementBase
//   The very base class for beam line representation objects.  A beam line
//   is modelled as a composite structure having a single root object
//   (the top level beam line), which contains both ``single'' leaf-type
//   elements (Components), as well as sub-lines (composites).
//
//   Interface for basic beam line object.
//   This class defines the abstract interface for all objects which can be
//   contained in a beam line. ElementBase forms the base class for two distinct
//   but related hierarchies of objects:
//   [OL]
//   [LI]
//   A set of concrete accelerator element classes, which compose the standard
//   accelerator component library (SACL).
//   [LI]
//   [/OL]
//   Instances of the concrete classes for single elements are by default
//   sharable. Instances of beam lines and integrators are by
//   default non-sharable, but they may be made sharable by a call to
//   [b]makeSharable()[/b].
//   [p]
//   An ElementBase object can return two lengths, which may be different:
//   [OL]
//   [LI]
//   The arc length along the geometry.
//   [LI]
//   The design length, often measured along a straight line.
//   [/OL]
//   Class ElementBase contains a map of name versus value for user-defined
//   attributes (see file AbsBeamline/AttributeSet.hh). The map is primarily
//   intended for processes that require algorithm-specific data in the
//   accelerator model.
//   [P]
//   The class ElementBase has as base class the abstract class RCObject.
//   Virtual derivation is used to make multiple inheritance possible for
//   derived concrete classes. ElementBase implements three copy modes:
//   [OL]
//   [LI]
//   Copy by reference: Call RCObject::addReference() and use [b]this[/b].
//   [LI]
//   Copy structure: use ElementBase::copyStructure().
//   During copying of the structure, all sharable items are re-used, while
//   all non-sharable ones are cloned.
//   [LI]
//   Copy by cloning: use ElementBase::clone().
//   This returns a full deep copy.
//   [/OL]
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
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/Channel.h"
#include <string>

#include "Structure/BoundaryGeometry.h"    // OPAL file
#include "Solvers/WakeFunction.hh"
#include "Solvers/ParticleMatterInteractionHandler.hh"


ElementBase::ElementBase():
    ElementBase("")
{}


ElementBase::ElementBase(const ElementBase &right):
    RCObject(),
    shareFlag(true),
    csTrafoGlobal2Local_m(right.csTrafoGlobal2Local_m),
    misalignment_m(right.misalignment_m),
    aperture_m(right.aperture_m),
    elementEdge_m(right.elementEdge_m),
    rotationZAxis_m(right.rotationZAxis_m),
    elementID(right.elementID),
    userAttribs(right.userAttribs),
    wake_m(right.wake_m),
    bgeometry_m(right.bgeometry_m),
    parmatint_m(right.parmatint_m),
    positionIsFixed(right.positionIsFixed),
    elementPosition_m(right.elementPosition_m),
    elemedgeSet_m(right.elemedgeSet_m)
{

    if(parmatint_m) {
        parmatint_m->updateElement(this);
    }
    if(bgeometry_m)
        bgeometry_m->updateElement(this);
}


ElementBase::ElementBase(const std::string &name):
    RCObject(),
    shareFlag(true),
    csTrafoGlobal2Local_m(),
    misalignment_m(),
    elementEdge_m(0),
    rotationZAxis_m(0.0),
    elementID(name),
    userAttribs(),
    wake_m(NULL),
    bgeometry_m(NULL),
    parmatint_m(NULL),
    positionIsFixed(false),
    elementPosition_m(0.0),
    elemedgeSet_m(false)
{}


ElementBase::~ElementBase()

{}


const std::string &ElementBase::getName() const

{
    return elementID;
}


void ElementBase::setName(const std::string &name) {
    elementID = name;
}


double ElementBase::getAttribute(const std::string &aKey) const {
    const ConstChannel *aChannel = getConstChannel(aKey);

    if(aChannel != NULL) {
        double val = *aChannel;
        delete aChannel;
        return val;
    } else {
        return 0.0;
    }
}


bool ElementBase::hasAttribute(const std::string &aKey) const {
    const ConstChannel *aChannel = getConstChannel(aKey);

    if(aChannel != NULL) {
        delete aChannel;
        return true;
    } else {
        return false;
    }
}


void ElementBase::removeAttribute(const std::string &aKey) {
    userAttribs.removeAttribute(aKey);
}


void ElementBase::setAttribute(const std::string &aKey, double val) {
    Channel *aChannel = getChannel(aKey, true);

    if(aChannel != NULL  &&  aChannel->isSettable()) {
        *aChannel = val;
        delete aChannel;
    } else
        std::cout << "Channel NULL or not Settable" << std::endl;
}


Channel *ElementBase::getChannel(const std::string &aKey, bool create) {
    return userAttribs.getChannel(aKey, create);
}


const ConstChannel *ElementBase::getConstChannel(const std::string &aKey) const {
    // Use const_cast to allow calling the non-const method GetChannel().
    // The const return value of this method will nevertheless inhibit set().
    return const_cast<ElementBase *>(this)->getChannel(aKey);
}


std::string ElementBase::getTypeString(ElementBase::ElementType type) {
    switch (type) {
    case BEAMBEAM:
        return "BeamBeam";
    case BEAMLINE:
        return "Beamline";
    case BEAMSTRIPPING:
        return "BeamStripping";
    case CCOLLIMATOR:
        return "CCollimator";
    case CORRECTOR:
        return "Corrector";
    case CYCLOTRON:
        return "Cyclotron";
    case DEGRADER:
        return "Degrader";
    case DIAGNOSTIC:
        return "Diagnostic";
    case DRIFT:
        return "Drift";
    case INTEGRATOR:
        return "Integrator";
    case LAMBERTSON:
        return "Lambertson";
    case MARKER:
        return "Marker";
    case MONITOR:
        return "Monitor";
    case MULTIPOLE:
        return "Multipole";
    case OFFSET:
        return "Offset";
    case PARALLELPLATE:
        return "ParallelPlate";
    case PATCH:
        return "Patch";
    case PROBE:
        return "Probe";
    case RBEND:
        return "RBend";
    case RFCAVITY:
        return "RFCavity";
    case RFQUADRUPOLE:
        return "RFQuadrupole";
    case RING:
        return "Ring";
    case SBEND3D:
        return "SBend3D";
    case SBEND:
        return "SBend";
    case SEPARATOR:
        return "Separator";
    case SEPTUM:
        return "Septum";
    case SOLENOID:
        return "Solenoid";
    case STRIPPER:
        return "Stripper";
    case TRAVELINGWAVE:
        return "TravelingWave";
    case VARIABLERFCAVITY:
        return "VariableRFCavity";
    case ANY:
    default:
        return "'unknown' type";
    }
}

ElementImage *ElementBase::getImage() const {
    std::string type = getTypeString();
    return new ElementImage(getName(), type, userAttribs);
}


ElementBase *ElementBase::copyStructure() {
    if(isSharable()) {
        return this;
    } else {
        return clone();
    }
}


void ElementBase::makeSharable() {
    shareFlag = true;
}


bool ElementBase::update(const AttributeSet &set) {
    for(AttributeSet::const_iterator i = set.begin(); i != set.end(); ++i) {
        setAttribute(i->first, i->second);
    }

    return true;
}

void ElementBase::setWake(WakeFunction *wk) {
    wake_m = wk;//->clone(getName() + std::string("_wake")); }
}

void ElementBase::setBoundaryGeometry(BoundaryGeometry *geo) {
    bgeometry_m = geo;//->clone(getName() + std::string("_wake")); }
}

void ElementBase::setParticleMatterInteraction(ParticleMatterInteractionHandler *parmatint) {
    parmatint_m = parmatint;
}

void ElementBase::setCurrentSCoordinate(double s) {
    if (actionRange_m.size() > 0 && actionRange_m.front().second < s) {
        actionRange_m.pop();
        if (actionRange_m.size() > 0) {
            elementEdge_m = actionRange_m.front().first;
        }
    }
}