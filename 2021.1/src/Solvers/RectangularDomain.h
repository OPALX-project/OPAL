//
// Class RectangularDomain
//   This class provides a rectangular beam pipe. The mesh adapts to the bunch
//   in longitudinal direction.
//
// Copyright (c) 2008,        Yves Ineichen, ETH Zürich,
//               2013 - 2015, Tülin Kaman, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2017 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the master thesis
// "A Parallel Multigrid Solver for Beam Dynamics"
// and the paper
// "A fast parallel Poisson solver on irregular domains applied to beam dynamics simulations"
// (https://doi.org/10.1016/j.jcp.2010.02.022)
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
#ifndef RECTANGULAR_DOMAIN_H
#define RECTANGULAR_DOMAIN_H

#include <vector>
#include <string>
#include "Solvers/RegularDomain.h"

class RectangularDomain : public RegularDomain {

public:
    /**
     * \param a is the longer side a of the rectangle
     * \param b is the shorter side b of the rectangle
     *
     */
    RectangularDomain(double a, double b, IntVector_t nr, Vector_t hr);

    /// calculates intersection with the beam pipe
    void compute(Vector_t hr, NDIndex<3> /*localId*/);

    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int /*z*/) const {
        double xx = (x - (nr_m[0] - 1) / 2.0) * hr_m[0];
        double yy = (y - (nr_m[1] - 1) / 2.0) * hr_m[1];
        return (xx <= getXRangeMax() && yy < getYRangeMax());
    }

private:
    /// conversion from (x,y,z) to index on the 3D grid
    int indexAccess(int x, int y, int z) const {
        return y * nr_m[0] + x + z * getNumXY();
    }

    int coordAccess(int idx) const {
        return idx % getNumXY();
    }
};

#endif

