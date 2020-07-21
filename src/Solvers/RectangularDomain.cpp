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
#ifdef HAVE_SAAMG_SOLVER

#include "Solvers/RectangularDomain.h"
#include "Utilities/OpalException.h"

RectangularDomain::RectangularDomain(Vector_t nr, Vector_t hr)
    : RectangularDomain(0.1, 0.1, nr, hr)
{}

RectangularDomain::RectangularDomain(double a, double b, Vector_t nr, Vector_t hr)
    : a_m(a)
    , b_m(b)
    , nxy_m(nr[0] * nr[1])
{
    setNr(nr);
    setHr(hr);
}

void RectangularDomain::compute(Vector_t hr){
    setHr(hr);
    nxy_m = nr[0] * nr[1];
}

void RectangularDomain::compute(Vector_t hr, NDIndex<3> /*localId*/) {
    compute(hr);
}

int RectangularDomain::getNumXY(int /*z*/) {
    return nxy_m;
}

void RectangularDomain::getBoundaryStencil(int x, int y, int z, double &W,
                                           double &E, double &S, double &N,
                                           double &F, double &B, double &C,
                                           double &scaleFactor)
{
    scaleFactor = 1.0;
    W = -1.0 / (hr[0] * hr[0]);
    E = -1.0 / (hr[0] * hr[0]);
    N = -1.0 / (hr[1] * hr[1]);
    S = -1.0 / (hr[1] * hr[1]);
    F = -1.0 / (hr[2] * hr[2]);
    B = -1.0 / (hr[2] * hr[2]);

    C = 2.0 / (hr[0] * hr[0])
      + 2.0 / (hr[1] * hr[1])
      + 2.0 / (hr[2] * hr[2]);

    if (!isInside(x + 1, y, z))
        E = 0.0; //we are a right boundary point

    if (!isInside(x - 1, y, z))
        W = 0.0; //we are a left boundary point

    if (!isInside(x, y + 1, z))
        N = 0.0; //we are a upper boundary point

    if (!isInside(x, y - 1, z))
        S = 0.0; //we are a lower boundary point

    //dirichlet
    if (z == 0)
        F = 0.0;

    if (z == nr[2] - 1)
        B = 0.0;
}

void RectangularDomain::getBoundaryStencil(int idx, double &W, double &E,
                                           double &S, double &N, double &F,
                                           double &B, double &C, double &scaleFactor)
{
    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);

}

void RectangularDomain::getNeighbours(int idx, int &W, int &E, int &S,
                                      int &N, int &F, int &B)
{
    int x = 0, y = 0, z = 0;

    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

void RectangularDomain::getNeighbours(int x, int y, int z, int &W, int &E,
                                      int &S, int &N, int &F, int &B)
{
    if (x > 0)
        W = getIdx(x - 1, y, z);
    else
        W = -1;
    if (x < nr[0] - 1)
        E = getIdx(x + 1, y, z);
    else
        E = -1;

    if (y < nr[1] - 1)
        N = getIdx(x, y + 1, z);
    else
        N = -1;
    if (y > 0)
        S = getIdx(x, y - 1, z);
    else
        S = -1;

    if (z > 0)
        F = getIdx(x, y, z - 1);
    else
        F = -1;
    if (z < nr[2] - 1)
        B = getIdx(x, y, z + 1);
    else
        B = -1;
}

void RectangularDomain::resizeMesh(Vector_t& origin, Vector_t& hr, const Vector_t& rmin,
                                const Vector_t& rmax, double dh)
{
    // apply bounding box increment, i.e., "BBOXINCR" input argument
    double zsize = rmax[2] - rmin[2];

    setMinMaxZ(rmin[2] - zsize * (1.0 + dh),
               rmax[2] + zsize * (1.0 + dh));

    origin = Vector_t(-a_m, -b_m, getMinZ());

    Vector_t mymax = Vector_t(a_m, b_m, getMaxZ());

    for (int i = 0; i < 3; ++i)
        hr[i]   = (mymax[i] - origin[i]) / nr[i];
}

#endif //#ifdef HAVE_SAAMG_SOLVER