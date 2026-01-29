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
#include "Algorithms/OpalParticle.h"


OpalParticle::OpalParticle()
{}


OpalParticle::OpalParticle(int64_t id,
                           double x, double px,
                           double y, double py,
                           double z, double pz,
                           double t,
                           double q, double m):
    id_m(id),
    R_m(x, y, z),
    P_m(px, py, pz),
    time_m(t),
    charge_m(q),
    mass_m(m)
{}

OpalParticle::OpalParticle(int64_t id,
                           Vector_t const& R, Vector_t const& P,
                           double t, double q, double m):
    id_m(id),
    R_m(R),
    P_m(P),
    time_m(t),
    charge_m(q),
    mass_m(m)
{}