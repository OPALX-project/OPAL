/*
 *  Copyright (c) 2014, Chris Rogers
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

#include "AbsBeamline/Offset.h"

#include <cmath>

#include "Utilities/GeneralClassicException.h"
#include "BeamlineGeometry/Euclid3DGeometry.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"

const double Offset::lengthUnits_m = 1e3;
double Offset::float_tolerance = 1e-12;

Offset::Offset(const std::string& name)
    : Component(name), _is_local(false), geometry_m(NULL) {
    geometry_m = new Euclid3DGeometry(Euclid3D());
    setElType(isDrift);
}

Offset::Offset()
    : Component(), _is_local(false), geometry_m(NULL) {
    geometry_m = new Euclid3DGeometry(Euclid3D());
    setElType(isDrift);
}


Offset::Offset(std::string name, const Offset& rhs)
    : Component(name), _is_local(false), geometry_m(NULL) {
    *this = rhs;
    setName(name);
}

Offset::Offset(const Offset& rhs)
    : Component(rhs.getName()), _is_local(false), geometry_m(NULL) {
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

    if (geometry_m != NULL)
        delete geometry_m;
    if (rhs.geometry_m == NULL) {
        geometry_m = NULL;
    } else {
        geometry_m = new Euclid3DGeometry(rhs.geometry_m->getTotalTransform());
    }
    setAperture(rhs.getAperture().first, rhs.getAperture().second);
    setElType(rhs.getElType());
    return *this;
}

void Offset::accept(BeamlineVisitor & visitor) const {
    visitor.visitOffset(*this);
}

EMField &Offset::getField() {
    throw GeneralClassicException("Offset::getField()",
                        "No field defined for Offset");
}

const EMField &Offset::getField() const {
    throw GeneralClassicException("Offset::getField() const",
                        "No field defined for Offset");
}

void Offset::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
}

void Offset::finalise() {
    RefPartBunch_m = NULL;
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

Euclid3DGeometry &Offset::getGeometry() {
    return *geometry_m;
}

const Euclid3DGeometry &Offset::getGeometry() const {
    return *geometry_m;
}

// std::ostream& operator<<(std::ostream& out, const Vector_t& vec) {
//     out << "(" << vec(0) << ", " << vec(1) << ", " << vec(2) << ")";
//     return out;
// }

double Offset::getTheta(Vector_t vec1, Vector_t vec2) {
    if (fabs(vec1(2)) > 1e-9 || fabs(vec2(2)) > 1e-9)
        throw GeneralClassicException("Offset::getTheta(...)",
                            "Rotations out of midplane are not implemented");
    // probably not the most efficient, but only called at set up
    double theta = atan2(vec2(1), vec2(0))-atan2(vec1(1), vec1(0));
    if (theta < -Physics::pi)
        theta += 2.*Physics::pi; // force into domain -pi < theta < pi
    return theta;
}

Vector_t Offset::rotate(Vector_t vec, double theta) {
    double s = sin(theta);
    double c = cos(theta);
    return Vector_t(+vec(0)*c-vec(1)*s,
                    +vec(0)*s+vec(1)*c,
                    0.);
}

void Offset::updateGeometry() {
    if (!_is_local)
        throw GeneralClassicException("Offset::updateGeometry(...)",
                            "Global offset needs a local coordinate system");
    Vector_t translation = getEndPosition();
    double length = sqrt(translation(0)*translation(0)+
                         translation(1)*translation(1)+
                         translation(2)*translation(2));
    double theta_in = getTheta(Vector_t(1., 0., 0.), translation);
    double theta_out = getTheta(Vector_t(1., 0., 0.), getEndDirection());
    Euclid3D euclid3D(-sin(theta_in)*length, 0., cos(theta_in)*length,
                      0., -theta_out, 0.);
    if (geometry_m != NULL)
        delete geometry_m;
    geometry_m = new Euclid3DGeometry(euclid3D);
}

void Offset::updateGeometry(Vector_t startPosition, Vector_t startDirection) {
    if (!_is_local) {
        Vector_t translationGlobal = _end_position;
        double theta_g2l = getTheta(startDirection, Vector_t(1, 0, 0));
        _end_position = rotate(translationGlobal, theta_g2l);
        _end_direction = rotate(_end_direction, theta_g2l);
        _is_local = true;
    }
    updateGeometry();
}

bool Offset::isGeometryAllocated() const {
    return geometry_m != NULL;
}

bool operator==(const Offset& off1, const Offset& off2) {
    const double tol = Offset::float_tolerance;
    if (off1.getName() != off2.getName() ||
        off1.getIsLocal() != off2.getIsLocal()) {
        return false;
    }
    for (int i = 0; i < 3; ++i) {
      if ( (fabs(off1.getEndPosition()(i)-off2.getEndPosition()(i)) > tol) ||
           (fabs(off1.getEndDirection()(i)-off2.getEndDirection()(i)) > tol))
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
        if (fabs(dTranslation(i)) > tol || fabs(dRotation(i)) > tol)
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
    if (geometry_m == NULL) {
        throw GeneralClassicException("Offset::bends",
              "Try to determine if Offset bends when geometry_m not allocated");
    }
    Rotation3D rotation = geometry_m->getTotalTransform().getRotation();
    for (size_t i = 0; i < 3; ++i)
        if (fabs(rotation.getAxis()(i)) > float_tolerance) {
            return true;
    }
    Vector3D vector = geometry_m->getTotalTransform().getVector();
    if (fabs(vector(0)) > float_tolerance || fabs(vector(1)) > float_tolerance) {
        return true;
    }
    return false;
}


Offset Offset::localCylindricalOffset(std::string name,
                                      double phi_in,
                                      double phi_out,
                                      double displacement) {
    Offset off(name);
    displacement *= lengthUnits_m;
    off.setEndPosition(Vector_t(cos(phi_in)*displacement,
                                sin(phi_in)*displacement,
                                0.));
    off.setEndDirection(Vector_t(cos(phi_in+phi_out), sin(phi_in+phi_out), 0.));
    off.setIsLocal(true);
    off.updateGeometry();
    return off;
}

Offset Offset::globalCylindricalOffset(std::string name,
                                       double radius_out,
                                       double phi_out,
                                       double theta_out) {
    Offset off(name);
    radius_out *= lengthUnits_m;
    off.setEndPosition(Vector_t(cos(phi_out)*radius_out,
                                sin(phi_out)*radius_out,
                                0.));
    off.setEndDirection(Vector_t(sin(phi_out+theta_out),
                                 cos(phi_out+theta_out),
                                 0.));
    off.setIsLocal(false);
    return off;
}

Offset Offset::localCartesianOffset(std::string name,
                                    Vector_t end_position,
                                    Vector_t end_direction) {
    Offset off(name);
    end_position *= lengthUnits_m;
    off.setEndPosition(end_position);
    off.setEndDirection(end_direction);
    off.setIsLocal(true);
    off.updateGeometry();
    return off;
}

Offset Offset::globalCartesianOffset(std::string name,
                                     Vector_t end_position,
                                     Vector_t end_direction) {
    Offset off(name);
    end_position *= lengthUnits_m;
    off.setEndPosition(end_position);
    off.setEndDirection(end_direction);
    off.setIsLocal(false);
    return off;
}