//
// Class LF2
//   Second order Leap-Frog time integrator
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
#include "BorisPusher.h"
#include "Physics/Units.h"

template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::doAdvance_m(PartBunchBase<double, 3>* bunch,
                                                    const size_t& i,
                                                    const double& t,
                                                    const double dt,
                                                    Arguments& ... args) const
{
    bool flagNoDeletion = true;

    // push for first LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * Units::ns2s);

    flagNoDeletion = kick_m(bunch, i, t, dt * Units::ns2s, args ...);

    // push for second LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * Units::ns2s);

    return flagNoDeletion;
}


template <typename FieldFunction, typename ... Arguments>
void LF2<FieldFunction, Arguments ...>::push_m(Vector_t& R, const Vector_t& P,
                                               const double& h) const
{
    double const gamma = sqrt(1.0 + dot(P, P));
    double const c_gamma = Physics::c / gamma;
    Vector_t const v = P * c_gamma;
    R += h * v;
}


template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::kick_m(PartBunchBase<double, 3>* bunch, const size_t& i,
                                               const double& t, const double& h,
                                               Arguments& ... args) const
{
    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);

    bool outOfBound = this->fieldfunc_m(t, i, externalE, externalB, args ...);

    if ( outOfBound )
        return false;


    double const q = bunch->Q[0] / Physics::q_e; // For now all particles have the same charge
    double const M = bunch->M[0] * Units::GeV2eV; // For now all particles have the same rest energy

    BorisPusher pusher;

    pusher.kick(bunch->R[i], bunch->P[i],
                externalE, externalB,
                h, M, q);

    return true;
}