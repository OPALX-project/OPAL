//
// Class IrregularDomain
//   Defines a common abstract interface for different types of boundaries.
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

#include "Solvers/IrregularDomain.h"

#include "Utility/PAssert.h"
#include "Utilities/OpalException.h"

IrregularDomain::IrregularDomain(const IntVector_t& nr, const Vector_t& hr,
                                 const std::string& interpl)
    : nr_m(nr)
    , hr_m(hr)
{
    if (interpl == "CONSTANT")
        interpolationMethod_m = CONSTANT;
    else if (interpl == "LINEAR")
        interpolationMethod_m = LINEAR;
    else if (interpl == "QUADRATIC")
        interpolationMethod_m = QUADRATIC;
}


void IrregularDomain::getNeighbours(int x, int y, int z, StencilIndex_t& index) const
{
    index.west  = -1;
    index.east  = -1;
    index.south = -1;
    index.north = -1;
    index.front = -1;
    index.back  = -1;

    if (x > 0 && isInside(x - 1, y, z))
        index.west = getIdx(x - 1, y, z);

    if (x < nr_m[0] - 1 && isInside(x + 1, y, z))
        index.east = getIdx(x + 1, y, z);

    if (y > 0 && isInside(x, y - 1, z))
        index.south = getIdx(x, y - 1, z);

    if (y < nr_m[1] - 1 && isInside(x, y + 1, z))
        index.north = getIdx(x, y + 1, z);

    if (z > 0 && isInside(x, y, z - 1))
        index.front = getIdx(x, y, z - 1);

    if (z < nr_m[2] - 1 && isInside(x, y, z + 1))
        index.back = getIdx(x, y, z + 1);
}


void IrregularDomain::getNeighbours(int id, StencilIndex_t& index) const {
    int x = 0, y = 0, z = 0;
    getCoord(id, x, y, z);
    getNeighbours(x, y, z, index);
}


void IrregularDomain::getCoord(int idx, int& x, int& y, int& z) const {
    int xy = coordAccess(idx);
    x = xy % nr_m[0];
    y = (xy - x) / nr_m[0];
    z = idx / getNumXY();
}

int IrregularDomain::getIdx(int x, int y, int z) const {
    if (x < 0 || y < 0 || z < 0 || !isInside(x, y, z))
        return -1;
    return indexAccess(x, y, z);
}

void IrregularDomain::getBoundaryStencil(int x, int y, int z, StencilValue_t& value,
                                         double &scaleFactor) const
{
    scaleFactor = 1.0;

    // determine which interpolation method we use for points near the boundary
    switch (interpolationMethod_m) {
        case CONSTANT:
            constantInterpolation(x, y, z, value, scaleFactor);
            break;
        case LINEAR:
            linearInterpolation(x, y, z, value, scaleFactor);
            break;
        case QUADRATIC:
            quadraticInterpolation(x, y, z, value, scaleFactor);
            break;
    }

    // stencil center value has to be positive!
    PAssert(value.center > 0);
}


void IrregularDomain::getBoundaryStencil(int id, StencilValue_t& value,
                                         double &scaleFactor) const
{
    int idx = 0, idy = 0, idz = 0;
    getCoord(id, idx, idy, idz);
    getBoundaryStencil(idx, idy, idz, value, scaleFactor);
}


void IrregularDomain::resizeMesh(Vector_t& origin, Vector_t& hr,
                                 const Vector_t& /*rmin*/, const Vector_t& /*rmax*/,
                                 double /*dh*/)
{
    origin = min_m;
    for (int i = 0; i < 3; i++)
        hr[i] = (max_m[i] - min_m[i]) / nr_m[i];
};


void IrregularDomain::constantInterpolation(int /*idx*/, int /*idy*/, int /*idz*/,
                                            StencilValue_t& /*value*/, double &/*scaleFactor*/) const
{
    throw OpalException("IrregularDomain::constantInterpolation()",
                        "No constant interpolation method Implemented!");
}


void IrregularDomain::linearInterpolation(int /*idx*/, int /*idy*/, int /*idz*/,
                                          StencilValue_t& /*value*/, double &/*scaleFactor*/) const
{
    throw OpalException("IrregularDomain::linearInterpolation()",
                        "No linear interpolation method Implemented!");
}


void IrregularDomain::quadraticInterpolation(int /*idx*/, int /*idy*/, int /*idz*/,
                                             StencilValue_t& /*value*/, double &/*scaleFactor*/) const
{
    throw OpalException("IrregularDomain::quadraticInterpolation()",
                        "No quadratic interpolation method Implemented!");
}