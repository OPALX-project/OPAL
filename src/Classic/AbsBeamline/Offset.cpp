//
// Class Offset
//   Defines the abstract interface for offset of elements.
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
#include "AbsBeamline/Offset.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

#include <cmath>

double Offset::float_tolerance = 1e-12;

Offset::Offset(const std::string& name)
    : Component(name), _is_local(false), geometry_m(nullptr) {
    geometry_m = new Euclid3DGeometry(Euclid3D());
}

Offset::Offset()
    : Offset("")
{}

Offset::Offset(const std::string& name, const Offset& rhs)
    : Component(name), _is_local(false), geometry_m(nullptr) {
    *this = rhs;
}

Offset::Offset(const Offset& rhs)
    : Component(rhs.getName()), _is_local(false), geometry_m(nullptr) {
    *this = rhs;
}

Offset::~Offset() {
    delete geometry_m;
}

Offset& Offset::operator=(const Offset& rhs) {
    if (&rhs == this) {
        return *this;
    }
    setName(rhs.getName());
    _end_position = rhs._end_position;
    _end_direction = rhs._end_direction;
    _is_local = rhs._is_local;

    if (geometry_m != nullptr)
        delete geometry_m;
    if (rhs.geometry_m == nullptr) {
        geometry_m = nullptr;
    } else {
        geometry_m = new Euclid3DGeometry(rhs.geometry_m->getTotalTransform());
    }
    setAperture(rhs.getAperture().first, rhs.getAperture().second);
    return *this;
}

void Offset::accept(BeamlineVisitor& visitor) const {
    visitor.visitOffset(*this);
}

EMField& Offset::getField() {
    throw GeneralClassicException("Offset::getField()",
                        "No field defined for Offset");
}

const EMField& Offset::getField() const {
    throw GeneralClassicException("Offset::getField() const",
                        "No field defined for Offset");
}

void Offset::initialise(PartBunchBase<double, 3>* bunch, double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
}

void Offset::finalise() {
    RefPartBunch_m = nullptr;
}

ElementBase* Offset::clone() const {
    return new Offset(*this);
}

void Offset::setEndPosition(Vector_t position) {
    _end_position = position;
}

Vector_t Offset::getEndPosition() const {
    return _end_position;
}

void Offset::setEndDirection(Vector_t direction) {
    _end_direction = direction;
}

Vector_t Offset::getEndDirection() const {
    return _end_direction;
}

void Offset::setIsLocal(bool isLocal) {
    _is_local = isLocal;
}

bool Offset::getIsLocal() const {
    return _is_local;
}

Euclid3DGeometry& Offset::getGeometry() {
    return *geometry_m;
}

const Euclid3DGeometry& Offset::getGeometry() const {
    return *geometry_m;
}

// std::ostream& operator<<(std::ostream& out, const Vector_t& vec) {
//     out << "(" << vec(0) << ", " << vec(1) << ", " << vec(2) << ")";
//     return out;
// }

double Offset::getTheta(Vector_t vec1, Vector_t vec2) {
    if (std::abs(vec1(2)) > 1e-9 || std::abs(vec2(2)) > 1e-9) {
        throw GeneralClassicException("Offset::getTheta",
                                      "Rotations out of midplane are not implemented");
    }

    // probably not the most efficient, but only called at set up
    double theta = std::atan2(vec2(1), vec2(0)) - std::atan2(vec1(1), vec1(0));
    if (theta < -Physics::pi) {
        theta += Physics::two_pi; // force into domain -pi < theta < pi
    }
    return theta;
}

Vector_t Offset::rotate(Vector_t vec, double theta) {
    double s = std::sin(theta);
    double c = std::cos(theta);
    return Vector_t({+vec(0)*c-vec(1)*s,
                    +vec(0)*s+vec(1)*c,
                    0.});
}

void Offset::updateGeometry() {
    if (!_is_local) {
        throw GeneralClassicException("Offset::updateGeometry",
                                      "Global offset needs a local coordinate system");
    }

    Vector_t translation = getEndPosition();
    double length = std::sqrt(translation(0) * translation(0) +
                              translation(1) * translation(1) +
                              translation(2) * translation(2));
    double theta_in  = getTheta(Vector_t({0., 1., 0.}), translation);
    double theta_out = getTheta(Vector_t({0., 1., 0.}), getEndDirection());
    Euclid3D euclid3D(-std::sin(theta_in) * length, 0., std::cos(theta_in) * length,
                      0., -theta_out, 0.);

    if (geometry_m != nullptr) {
        delete geometry_m;
    }

    geometry_m = new Euclid3DGeometry(euclid3D);
}

void Offset::updateGeometry(Vector_t startPosition, Vector_t startDirection) {
    if (!_is_local) {
        // thetaIn is the angle between the y axis and startDirection
        double thetaIn = std::atan2(-startDirection[0], startDirection[1]); // global OPAL-CYCL coords
        // thetaOut is the angle between the y axis and endDirection
        double thetaOut = std::atan2(-_end_direction[0], _end_direction[1]); // global OPAL-CYCL coords
        // thetaRel is the angle between thetaOut and thetaIn
        double thetaRel = thetaOut-thetaIn;
        // deltaPosition is the position change in the global coordinate system
        Vector_t deltaPosition = _end_position-startPosition;
        // endPosition is the difference between end and startPosition in 
        // startDirection coordinate system
        _end_position = rotate(deltaPosition, -thetaIn);
        // end direction is the normal in the coordinate system of startDirection
        _end_direction = Vector_t({std::sin(-thetaRel), std::cos(-thetaRel), 0});
        _is_local = true;
    }
    updateGeometry();
}

bool Offset::isGeometryAllocated() const {
    return geometry_m != nullptr;
}

bool operator==(const Offset& off1, const Offset& off2) {
    const double tol = Offset::float_tolerance;
    if (off1.getName() != off2.getName() ||
        off1.getIsLocal() != off2.getIsLocal()) {
        return false;
    }
    for (int i = 0; i < 3; ++i) {
      if ( (std::abs(off1.getEndPosition()(i)-off2.getEndPosition()(i)) > tol) ||
           (std::abs(off1.getEndDirection()(i)-off2.getEndDirection()(i)) > tol))
            return false;
    }
    if ( (!off1.isGeometryAllocated() && off2.isGeometryAllocated()) ||
         (!off2.isGeometryAllocated() && off1.isGeometryAllocated()))
        return false;
    Euclid3D transform1 = off1.getGeometry().getTotalTransform();
    Euclid3D transform2 = off2.getGeometry().getTotalTransform();
    Vector3D dTranslation = transform1.getVector() - transform2.getVector();
    Vector3D dRotation = transform1.getRotation().getAxis() -
                         transform2.getRotation().getAxis();
    for (size_t i = 0; i < 3; ++i)
        if (std::abs(dTranslation(i)) > tol || std::abs(dRotation(i)) > tol)
            return false;
    return true;
}

bool operator!=(const Offset& off1, const Offset& off2) {
    return !(off1 == off2);
}

std::ostream& operator<<(std::ostream& out, const Offset& off) {
    out << "Offset " << off.getName() << " local " << off.getIsLocal()
        << " end pos: " << off.getEndPosition()
        << " end dir: " << off.getEndDirection() << std::endl;
    return out;
}

bool Offset::bends() const {
    if (geometry_m == nullptr) {
        throw GeneralClassicException("Offset::bends",
              "Try to determine if Offset bends when geometry_m not allocated");
    }
    Rotation3D rotation = geometry_m->getTotalTransform().getRotation();
    for (size_t i = 0; i < 3; ++i)
        if (std::abs(rotation.getAxis()(i)) > float_tolerance) {
            return true;
    }
    Vector3D vector = geometry_m->getTotalTransform().getVector();
    if (std::abs(vector(0)) > float_tolerance || std::abs(vector(1)) > float_tolerance) {
        return true;
    }
    return false;
}


Offset Offset::localCylindricalOffset(const std::string& name,
                                      double phi_in,
                                      double phi_out,
                                      double displacement) {
    Offset off(name);
    off.setEndPosition(Vector_t({-std::sin(phi_in)*displacement,
                                std::cos(phi_in)*displacement,
                                0.}));
    off.setEndDirection(Vector_t({-std::sin(phi_in+phi_out), std::cos(phi_in+phi_out), 0.}));
    off.setIsLocal(true);
    off.updateGeometry();
    return off;
}

Offset Offset::globalCylindricalOffset(const std::string& name,
                                       double radius_out,
                                       double phi_out,
                                       double theta_out) {
    Offset off(name);
    off.setEndPosition(Vector_t({std::cos(phi_out)*radius_out,
                                std::sin(phi_out)*radius_out,
                                0.}));
    off.setEndDirection(Vector_t({std::sin(phi_out+theta_out),
                                 std::cos(phi_out+theta_out),
                                 0.}));
    off.setIsLocal(false);
    return off;
}

Offset Offset::localCartesianOffset(const std::string& name,
                                    Vector_t end_position,
                                    Vector_t end_direction) {
    Offset off(name);
    off.setEndPosition(end_position);
    off.setEndDirection(end_direction);
    off.setIsLocal(true);
    off.updateGeometry();
    return off;
}

Offset Offset::globalCartesianOffset(const std::string& name,
                                     Vector_t end_position,
                                     Vector_t end_direction) {
    Offset off(name);
    off.setEndPosition(end_position);
    off.setEndDirection(end_direction);
    off.setIsLocal(false);
    return off;
}
