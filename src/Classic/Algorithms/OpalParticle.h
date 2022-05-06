//
// Class OpalParticle
//   This class represents the canonical coordinates of a particle.
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
#ifndef CLASSIC_OpalParticle_HH
#define CLASSIC_OpalParticle_HH

#include "Vektor.h"

class OpalParticle {

public:

    // Particle coordinate numbers.
    enum { X, Y, L, INVALID };

    /// Constructor.
    //  Construct particle with the given coordinates.
    OpalParticle(int64_t id,
                 double x, double px,
                 double y, double py,
                 double z, double pz,
                 double time,
                 double q, double m);

    OpalParticle(int64_t id,
                 Vector_t const& R, Vector_t const& P,
                 double time, double q, double m);

    OpalParticle();

    /// Set the horizontal position in m.
    void setX(double) ;

    /// Set the horizontal momentum.
    void setPx(double);

    /// Set the vertical displacement in m.
    void setY(double) ;

    /// Set the vertical momentum.
    void setPy(double);

    /// Set longitudinal position in m.
    void setZ(double) ;

    /// Set the longitudinal momentum
    void setPz(double);

    /// Set position in m
    void setR(Vector_t const&);

    /// Set momentum
    void setP(Vector_t const&);

    /// Set the time
    void setTime(double t);

    /// Get the id of the particle
    int64_t getId() const;

    /// Get coordinate.
    //  Access coordinate by index for constant particle.
    double operator[](unsigned int) const;

    /// Get horizontal position in m.
    double getX() const;

    /// Get horizontal momentum (no dimension).
    double getPx() const;

    /// Get vertical displacement in m.
    double getY() const;

    /// Get vertical momentum (no dimension).
    double getPy() const;

    /// Get longitudinal displacement c*t in m.
    double getZ() const;

    /// Get relative momentum error (no dimension).
    double getPz() const;

    /// Get position in m
    const Vector_t& getR() const;

    /// Get momentum
    const Vector_t& getP() const;

    /// Get time
    double getTime() const;

    /// Get charge in Coulomb
    double getCharge() const;

    /// Get mass in GeV/c^2
    double getMass() const;

private:
    int64_t id_m;
    Vector_t R_m;
    Vector_t P_m;
    double time_m;
    double charge_m;
    double mass_m;
};

inline
void OpalParticle::setX(double val)
{
    R_m[X] = val;
}

inline
void OpalParticle::setY(double val)
{
    R_m[Y] = val;
}

inline
void OpalParticle::setZ(double val)
{
    R_m[L] = val;
}

inline
void OpalParticle::setPx(double val)
{
    P_m[X] = val;
}

inline
void OpalParticle::setPy(double val)
{
    P_m[Y] = val;
}

inline
void OpalParticle::setPz(double val)
{
    P_m[L]  = val;
}

inline
void OpalParticle::setR(Vector_t const& R)
{
    R_m = R;
}

inline
void OpalParticle::setP(Vector_t const& P)
{
    P_m = P;
}

inline
void OpalParticle::setTime(double t)
{
    time_m = t;
}

inline
int64_t OpalParticle::getId() const
{
    return id_m;
}

inline
double OpalParticle::operator[](unsigned int i) const
{
    PAssert_LT(i, 6u);
    return i % 2 == 0? R_m[i / 2]: P_m[i / 2];
}

inline
double OpalParticle::getX() const
{
    return R_m[X];
}

inline
double OpalParticle::getY() const
{
    return R_m[Y];
}

inline
double OpalParticle::getZ() const
{
    return R_m[L];
}

inline
double OpalParticle::getPx() const
{
    return P_m[X];
}

inline
double OpalParticle::getPy() const
{
    return P_m[Y];
}

inline
double OpalParticle::getPz() const
{
    return P_m[L];
}

inline
const Vector_t& OpalParticle::getR() const
{
    return R_m;
}

inline
const Vector_t& OpalParticle::getP() const
{
    return P_m;
}

inline
double OpalParticle::getTime() const
{
    return time_m;
}

inline
double OpalParticle::getCharge() const
{
    return charge_m;
}

inline
double OpalParticle::getMass() const
{
    return mass_m;
}

#endif // CLASSIC_OpalParticle_HH