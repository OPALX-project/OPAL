// ------------------------------------------------------------------------
// $RCSfile: Rotation3D.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Rotation3D
//   Represents a rotation in 3-dimensional space.
//
//   There are two possible external representations:
//
//   - Matrix3D R.
//     The matrix R is an element of SO(3). Its column vectors define the
//     unit vectors of the rotated coordinate system, expressed in the
//     original system. This representation is used internally and can only
//     be read from a Rotation3D object. To build such an object, use the
//     other representation.
//
//   - Vector3D V.
//     The direction of V defines the axis of rotation, and its length the
//     rotation angle. A zero vector represents the identity.
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:32:00 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Rotation3D.h"
#include <algorithm>
#include <cmath>


// Class Rotation3D
// ------------------------------------------------------------------------

Rotation3D::Rotation3D(double vx, double vy, double vz):
    R() {
    double phi = std::sqrt(vx * vx + vy * vy + vz * vz);

    if(phi != 0.0) {
        double factor = std::sin(phi / 2.0) / phi;
        double S0 = std::cos(phi / 2.0);
        double S1 = factor * vx;
        double S2 = factor * vy;
        double S3 = factor * vz;

        R(0, 0) = 2.0 * (S1 * S1 + S0 * S0) - 1.0;
        R(0, 1) = 2.0 * (S1 * S2 - S0 * S3);
        R(0, 2) = 2.0 * (S1 * S3 + S0 * S2);

        R(1, 0) = 2.0 * (S2 * S1 + S0 * S3);
        R(1, 1) = 2.0 * (S2 * S2 + S0 * S0) - 1.0;
        R(1, 2) = 2.0 * (S2 * S3 - S0 * S1);

        R(2, 0) = 2.0 * (S3 * S1 - S0 * S2);
        R(2, 1) = 2.0 * (S3 * S2 + S0 * S1);
        R(2, 2) = 2.0 * (S3 * S3 + S0 * S0) - 1.0;
    }
}


void Rotation3D::getAxis(double &vx, double &vy, double &vz) const {
    vx = (R(2, 1) - R(1, 2)) / 2.0;
    vy = (R(0, 2) - R(2, 0)) / 2.0;
    vz = (R(1, 0) - R(0, 1)) / 2.0;
    double c = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) / 2.0;
    double s = std::sqrt(vx * vx + vy * vy + vz * vz);
    double phi = std::atan2(s, c);

    if(c < 0.0) {
        // NOTE: We must avoid negative arguments to sqrt(),
        // which may occur due to rounding errors.
        vx = (vx > 0.0 ? phi : (-phi)) * std::sqrt(std::max(R(0, 0) - c, 0.0) / (1.0 - c));
        vy = (vy > 0.0 ? phi : (-phi)) * std::sqrt(std::max(R(1, 1) - c, 0.0) / (1.0 - c));
        vz = (vz > 0.0 ? phi : (-phi)) * std::sqrt(std::max(R(2, 2) - c, 0.0) / (1.0 - c));
    } else if(std::abs(s) > 1.0e-10) {
        double factor = phi / s;
        vx *= factor;
        vy *= factor;
        vz *= factor;
    }
}


Vector3D Rotation3D::getAxis() const {
    double vx, vy, vz;
    getAxis(vx, vy, vz);
    return Vector3D(vx, vy, vz);
}


Rotation3D Rotation3D::inverse() const {
    Rotation3D inv;

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            inv.R(i, j) = R(j, i);
        }
    }

    return inv;
}


bool Rotation3D::isPureXRotation() const {
    return (R(0, 1) == 0.0) && (R(1, 0) == 0.0) && (R(0, 2) == 0.0) && (R(2, 0) == 0.0);
}


bool Rotation3D::isPureYRotation() const {
    return (R(1, 2) == 0.0) && (R(2, 1) == 0.0) && (R(1, 0) == 0.0) && (R(0, 1) == 0.0);
}


bool Rotation3D::isPureZRotation() const {
    return (R(2, 0) == 0.0) && (R(0, 2) == 0.0) && (R(2, 1) == 0.0) && (R(1, 2) == 0.0);
}


Rotation3D Rotation3D::Identity() {
    return Rotation3D();
}


Rotation3D Rotation3D::XRotation(double angle) {
    Rotation3D result;
    result.R(1, 1) =    result.R(2, 2) = std::cos(angle);
    result.R(1, 2) = - (result.R(2, 1) = std::sin(angle));
    return result;
}


Rotation3D Rotation3D::YRotation(double angle) {
    Rotation3D result;
    result.R(2, 2) =    result.R(0, 0) = std::cos(angle);
    result.R(2, 0) = - (result.R(0, 2) = std::sin(angle));
    return result;
}


Rotation3D Rotation3D::ZRotation(double angle) {
    Rotation3D result;
    result.R(0, 0) =    result.R(1, 1) = std::cos(angle);
    result.R(0, 1) = - (result.R(1, 0) = std::sin(angle));
    return result;
}