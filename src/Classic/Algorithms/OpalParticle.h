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
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  Construct particle with the given coordinates.
    OpalParticle(double x, double px,
                 double y, double py,
                 double t, double pt,
                 double q, double m);

    OpalParticle();

    /// Get coordinate.
    //  Access coordinate by index.
    //  Note above order of phase space coordinates!
    double &operator[](int);

    /// Get reference to horizontal position in m.
    double &x() ;

    /// Get reference to horizontal momentum (no dimension).
    double &px();

    /// Get reference to vertical displacement in m.
    double &y() ;

    /// Get reference to vertical momentum (no dimension).
    double &py();

    /// Get reference to longitudinal displacement c*t in m.
    double &t() ;

    /// Get reference to relative momentum error (no dimension).
    double &pt();

    /// Get reference to position in m
    Vector_t &R();

    /// Get reference to momentum
    Vector_t &P();

    /// Get coordinate.
    //  Access coordinate by index for constant particle.
    double operator[](int) const;

    /// Get horizontal position in m.
    double x() const;

    /// Get horizontal momentum (no dimension).
    double px() const;

    /// Get vertical displacement in m.
    double y() const;

    /// Get vertical momentum (no dimension).
    double py() const;

    /// Get longitudinal displacement c*t in m.
    double t() const;

    /// Get relative momentum error (no dimension).
    double pt() const;

    /// Get position in m
    Vector_t R() const;

    /// Get momentum
    Vector_t P() const;

    /// Get charge in Coulomb
    double charge() const;

    /// Get mass in GeV/c^2
    double mass() const;

private:

    // The particle's phase space coordinates:
    Vector_t R_m;
    Vector_t P_m;
    double charge_m;
    double mass_m;
};


// Inline member functions.
// ------------------------------------------------------------------------

inline double &OpalParticle::operator[](int i)
{
    return i < PX? R_m[i] : P_m[i - PX];
}

inline double &OpalParticle::x()
{
    return R_m[X];
}

inline double &OpalParticle::y()
{
    return R_m[Y];
}

inline double &OpalParticle::t()
{
    return R_m[T];
}

inline double &OpalParticle::px()
{
    return P_m[X];
}

inline double &OpalParticle::py()
{
    return P_m[Y];
}

inline double &OpalParticle::pt()
{
    return P_m[T];
}

inline Vector_t &OpalParticle::R()
{
    return R_m;
}

inline Vector_t &OpalParticle::P()
{
    return P_m;
}

inline double OpalParticle::operator[](int i) const
{
    return i < PX? R_m[i]: P_m[i - PX];
}

inline double OpalParticle::x() const
{
    return R_m[X];
}

inline double OpalParticle::y() const
{
    return R_m[Y];
}

inline double OpalParticle::t() const
{
    return R_m[T];
}

inline double OpalParticle::px() const
{
    return P_m[X];
}

inline double OpalParticle::py() const
{
    return P_m[Y];
}

inline double OpalParticle::pt() const
{
    return P_m[T];
}

inline Vector_t OpalParticle::R() const
{
    return R_m;
}

inline Vector_t OpalParticle::P() const
{
    return P_m;
}

inline double OpalParticle::charge() const
{
    return charge_m;
}

inline double OpalParticle::mass() const
{
    return mass_m;
}

#endif // CLASSIC_OpalParticle_HH