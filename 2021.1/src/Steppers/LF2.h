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
#ifndef LF2_H
#define LF2_H

#include "Stepper.h"
#include "Physics/Physics.h"

/// Leap-Frog 2nd order
template <typename FieldFunction, typename ... Arguments>
class LF2 : public Stepper<FieldFunction, Arguments...> {
    
public:
    
    LF2(const FieldFunction& fieldfunc) : Stepper<FieldFunction, Arguments ...>(fieldfunc) { }
    
private:
    bool doAdvance_m(PartBunchBase<double, 3>* bunch,
                     const size_t& i,
                     const double& t,
                     const double dt,
                     Arguments& ... args) const;
    
    
    void push_m(Vector_t& R, const Vector_t& P, const double& h) const;
    
    bool kick_m(PartBunchBase<double, 3>* bunch, const size_t& i,
                const double& t, const double& h,
                Arguments& ... args) const;
};

#include "LF2.hpp"

#endif
