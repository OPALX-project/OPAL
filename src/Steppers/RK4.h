//
// Class RK4
//   Fourth order Runge-Kutta time integrator
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
#ifndef RK4_H
#define RK4_H

#include "Stepper.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

/// 4-th order Runnge-Kutta stepper
template <typename FieldFunction, typename ... Arguments>
class RK4 : public Stepper<FieldFunction, Arguments...> {
    
public:
    
    RK4(const FieldFunction& fieldfunc) : Stepper<FieldFunction, Arguments ...>(fieldfunc) { }

private:
    bool doAdvance_m(PartBunchBase<double, 3>* bunch,
                     const size_t& i,
                     const double& t,
                     const double dt,
                     Arguments& ... args) const;
    
    /**
     * 
     *
     * @param y
     * @param t
     * @param yp
     * @param Pindex
     *
     * @return
     */
    bool derivate_m(PartBunchBase<double, 3>* bunch,
                    double *y,
                    const double& t,
                    double* yp,
                    const size_t& i,
                    Arguments& ... args) const;
    
    
    void copyTo(const Vector_t& R, const Vector_t& P, double* x) const;
    
    void copyFrom(Vector_t& R, Vector_t& P, double* x) const;
    
    const double mass_coeff = 1.0e9 * Units::GeV2kg; // from GeV/c^2 to basic unit: GV*C*s^2/m^2, (1.0e9 converts V*C*s^2/m^2 to GV*C*s^2/m^2)
    const double c_mmtns = Physics::c * Units::m2mm / Units::s2ns;
    const double c_mtns = Physics::c / Units::s2ns;
};

#include "RK4.hpp"

#endif
