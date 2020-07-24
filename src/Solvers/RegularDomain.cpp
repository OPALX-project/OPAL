//
// Class RegularDomain
//   Base class for simple domains that do not change the x-y shape in
//   longitudinal direction.
//
// Copyright (c) 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Solvers/RegularDomain.h"

RegularDomain::RegularDomain(const IntVector_t& nr,
                             const Vector_t& hr,
                             const std::string& interpl)
    : IrregularDomain(nr, hr, interpl)
    , nxy_m(nr[0] * nr[1])
{ }


void RegularDomain::resizeMesh(Vector_t& origin, Vector_t& hr, const Vector_t& rmin,
                                const Vector_t& rmax, double dh)
{
    Vector_t mymax = Vector_t(0.0, 0.0, 0.0);

    // apply bounding box increment dh, i.e., "BBOXINCR" input argument
    double zsize = rmax[2] - rmin[2];

    setMinMaxZ(rmin[2] - zsize * (1.0 + dh),
               rmax[2] + zsize * (1.0 + dh));

    origin = Vector_t(getXRangeMin(), getYRangeMin(), getMinZ());
    mymax  = Vector_t(getXRangeMax(), getYRangeMax(), getMaxZ());

    for (int i = 0; i < 3; ++i)
        hr[i] = (mymax[i] - origin[i]) / nr_m[i];
}


void RegularDomain::constantInterpolation(int x, int y, int z, StencilValue_t& value,
                                          double &scaleFactor) const
{
    scaleFactor = 1.0;
    value.west  = -1.0 / (hr_m[0] * hr_m[0]);
    value.east  = -1.0 / (hr_m[0] * hr_m[0]);
    value.north = -1.0 / (hr_m[1] * hr_m[1]);
    value.south = -1.0 / (hr_m[1] * hr_m[1]);
    value.front = -1.0 / (hr_m[2] * hr_m[2]);
    value.back  = -1.0 / (hr_m[2] * hr_m[2]);

    value.center = 2.0 / (hr_m[0] * hr_m[0])
                 + 2.0 / (hr_m[1] * hr_m[1])
                 + 2.0 / (hr_m[2] * hr_m[2]);

    // we are a right boundary point
    if (!isInside(x + 1, y, z))
        value.east = 0.0;

    // we are a left boundary point
    if (!isInside(x - 1, y, z))
        value.west = 0.0;

    // we are a upper boundary point
    if (!isInside(x, y + 1, z))
        value.north = 0.0;

    // we are a lower boundary point
    if (!isInside(x, y - 1, z))
        value.south = 0.0;

    // handle boundary condition in z direction
    robinBoundaryStencil(z, value.front, value.back, value.center);
}


void RegularDomain::robinBoundaryStencil(int z, double &F, double &B, double &C) const {
    if (z == 0 || z == nr_m[2] - 1) {

        // case where we are on the Robin BC in Z-direction
        // where we distinguish two cases
        // IFF: this values should not matter because they
        // never make it into the discretization matrix
        if (z == 0)
            F = 0.0;
        else
            B = 0.0;

        // add contribution of Robin discretization to center point
        // d the distance between the center of the bunch and the boundary
        double d = 0.5 * hr_m[2] * (nr_m[2] - 1);
        C += 2.0 / (d * hr_m[2]);
    }
}
