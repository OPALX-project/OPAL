//
// Class EllipticDomain
//   This class provides an elliptic beam pipe. The mesh adapts to the bunch size
//   in the longitudinal direction. At the intersection of the mesh with the
//   beam pipe, three stencil interpolation methods are available.
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
#include "Solvers/EllipticDomain.h"

#include <map>
#include <cmath>
#include <iostream>
#include <cassert>

//FIXME: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)


EllipticDomain::EllipticDomain(Vector_t nr, Vector_t hr) {
    setNr(nr);
    setHr(hr);
}

EllipticDomain::EllipticDomain(double semimajor, double semiminor, Vector_t nr,
                               Vector_t hr, std::string interpl)
    : semiMajor_m(semimajor)
    , semiMinor_m(semiminor)
{
    setNr(nr);
    setHr(hr);

    if (interpl == "CONSTANT")
        interpolationMethod_m = CONSTANT;
    else if (interpl == "LINEAR")
        interpolationMethod_m = LINEAR;
    else if (interpl == "QUADRATIC")
        interpolationMethod_m = QUADRATIC;
}

EllipticDomain::EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr,
                               std::string interpl)
    : EllipticDomain(bgeom->getA(), bgeom->getB(), nr, hr, interpl)
{
    setMinMaxZ(bgeom->getS(), bgeom->getLength());
}

EllipticDomain::EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr, std::string interpl)
    : EllipticDomain(bgeom, nr,
                     Vector_t(2.0 * bgeom->getA() / nr[0],
                              2.0 * bgeom->getB() / nr[1],
                              (bgeom->getLength() - bgeom->getS()) / nr[2]),
                     interpl)
{ }

EllipticDomain::~EllipticDomain() {
    //nothing so far
}


// for this geometry we only have to calculate the intersection on
// one x-y-plane
// for the moment we center the ellipse around the center of the grid
void EllipticDomain::compute(Vector_t hr, NDIndex<3> localId){
    // there is nothing to be done if the mesh spacings in x and y have not changed
    if (hr[0] == getHr()[0] &&
        hr[1] == getHr()[1])
    {
        hasGeometryChanged_m = false;
        return;
    }

    setHr(hr);
    hasGeometryChanged_m = true;
    //reset number of points inside domain
    nxy_m = 0;

    // clear previous coordinate maps
    idxMap_m.clear();
    coordMap_m.clear();
    //clear previous intersection points
    intersectYDir_m.clear();
    intersectXDir_m.clear();

    // build a index and coordinate map
    int idx = 0;
    int x, y;

    /* we need to scan the full x-y-plane on all cores
     * in order to figure out the number of valid
     * grid points per plane --> otherwise we might
     * get not unique global indices in the Tpetra::CrsMatrix
     */
    for (x = 0; x < nr[0]; ++x) {
        for (y = 0; y < nr[1]; ++y) {
            if (isInside(x, y, 1)) {
                idxMap_m[toCoordIdx(x, y)] = idx;
                coordMap_m[idx++] = toCoordIdx(x, y);
                nxy_m++;
            }
        }
    }

    switch (interpolationMethod_m) {
        case CONSTANT:
            break;
        case LINEAR:
        case QUADRATIC:

            double smajsq = semiMajor_m * semiMajor_m;
            double sminsq = semiMinor_m * semiMinor_m;
            double yd = 0.0;
            double xd = 0.0;
            double pos = 0.0;

            // calculate intersection with the ellipse
            for (x = localId[0].first(); x <= localId[0].last(); x++) {
                pos = - semiMajor_m + hr[0] * (x + 0.5);
                if (std::abs(pos) >= semiMajor_m)
                {
                    intersectYDir_m.insert(std::pair<int, double>(x, 0));
                    intersectYDir_m.insert(std::pair<int, double>(x, 0));
                } else {
                    yd = std::abs(std::sqrt(sminsq - sminsq * pos * pos / smajsq));
                    intersectYDir_m.insert(std::pair<int, double>(x, yd));
                    intersectYDir_m.insert(std::pair<int, double>(x, -yd));
                }

            }

            for (y = localId[0].first(); y < localId[1].last(); y++) {
                pos = - semiMinor_m + hr[1] * (y + 0.5);
                if (std::abs(pos) >= semiMinor_m)
                {
                    intersectXDir_m.insert(std::pair<int, double>(y, 0));
                    intersectXDir_m.insert(std::pair<int, double>(y, 0));
                } else {
                    xd = std::abs(std::sqrt(smajsq - smajsq * pos * pos / sminsq));
                    intersectXDir_m.insert(std::pair<int, double>(y, xd));
                    intersectXDir_m.insert(std::pair<int, double>(y, -xd));
                }
            }
    }
}

void EllipticDomain::resizeMesh(Vector_t& origin, Vector_t& hr, const Vector_t& rmin,
                                const Vector_t& rmax, double dh)
{
    Vector_t mymax = Vector_t(0.0, 0.0, 0.0);

    // apply bounding box increment dh, i.e., "BBOXINCR" input argument
    double zsize = rmax[2] - rmin[2];

    setMinMaxZ(rmin[2] - zsize * (1.0 + dh),
               rmax[2] + zsize * (1.0 + dh));

    origin = Vector_t(-semiMajor_m, -semiMinor_m, getMinZ());
    mymax  = Vector_t( semiMajor_m,  semiMinor_m, getMaxZ());

    for (int i = 0; i < 3; ++i)
        hr[i] = (mymax[i] - origin[i]) / nr[i];
}

void EllipticDomain::getBoundaryStencil(int x, int y, int z, double &W,
                                        double &E, double &S, double &N,
                                        double &F, double &B, double &C,
                                        double &scaleFactor)
{
    scaleFactor = 1.0;

    // determine which interpolation method we use for points near the boundary
    switch (interpolationMethod_m) {
        case CONSTANT:
            constantInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
        case LINEAR:
            linearInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
        case QUADRATIC:
            quadraticInterpolation(x, y, z, W, E, S, N, F, B, C, scaleFactor);
            break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void EllipticDomain::getBoundaryStencil(int idx, double &W, double &E,
                                        double &S, double &N, double &F,
                                        double &B, double &C, double &scaleFactor)
{
    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);
}


void EllipticDomain::getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B) {

    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

void EllipticDomain::getNeighbours(int x, int y, int z, int &W, int &E,
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

void EllipticDomain::constantInterpolation(int x, int y, int z, double &WV,
                                           double &EV, double &SV, double &NV,
                                           double &FV, double &BV, double &CV,
                                           double &scaleFactor)
{
    scaleFactor = 1.0;

    WV = -1.0 / (hr[0] * hr[0]);
    EV = -1.0 / (hr[0] * hr[0]);
    NV = -1.0 / (hr[1] * hr[1]);
    SV = -1.0 / (hr[1] * hr[1]);
    FV = -1.0 / (hr[2] * hr[2]);
    BV = -1.0 / (hr[2] * hr[2]);

    CV =  2.0 / (hr[0] * hr[0])
       +  2.0 / (hr[1] * hr[1])
       +  2.0 / (hr[2] * hr[2]);

    // we are a right boundary point
    if (!isInside(x + 1, y, z))
        EV = 0.0;

    // we are a left boundary point
    if (!isInside(x - 1, y, z))
        WV = 0.0;

    // we are a upper boundary point
    if (!isInside(x, y + 1, z))
        NV = 0.0;

    // we are a lower boundary point
    if (!isInside(x, y - 1, z))
        SV = 0.0;

    // handle boundary condition in z direction
    robinBoundaryStencil(z, FV, BV, CV);
}

void EllipticDomain::linearInterpolation(int x, int y, int z, double &W,
                                         double &E, double &S, double &N,
                                         double &F, double &B, double &C,
                                         double &scaleFactor)
{
    scaleFactor = 1.0;

    double cx = - semiMajor_m + hr[0] * (x + 0.5);
    double cy = - semiMinor_m + hr[1] * (y + 0.5);

    double dx = 0.0;
    std::multimap<int, double>::iterator it = intersectXDir_m.find(y);

    if (cx < 0)
        ++it;
    dx = it->second;

    double dy = 0.0;
    it = intersectYDir_m.find(x);
    if (cy < 0)
        ++it;
    dy = it->second;


    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    C = 0.0;

    // we are a right boundary point
    if (!isInside(x + 1, y, z)) {
        C += 1.0 / ((dx - cx) * de);
        E = 0.0;
    } else {
        C += 1.0 / (de * de);
        E = -1.0 / (de * de);
    }

    // we are a left boundary point
    if (!isInside(x - 1, y, z)) {
        C += 1.0 / ((std::abs(dx) - std::abs(cx)) * dw);
        W = 0.0;
    } else {
        C += 1.0 / (dw * dw);
        W = -1.0 / (dw * dw);
    }

    // we are a upper boundary point
    if (!isInside(x, y + 1, z)) {
        C += 1.0 / ((dy - cy) * dn);
        N = 0.0;
    } else {
        C += 1.0 / (dn * dn);
        N = -1.0 / (dn * dn);
    }

    // we are a lower boundary point
    if (!isInside(x, y - 1, z)) {
        C += 1.0 / ((std::abs(dy) - std::abs(cy)) * ds);
        S = 0.0;
    } else {
        C += 1.0 / (ds * ds);
        S = -1.0 / (ds * ds);
    }

    // handle boundary condition in z direction
    F = -1.0 / (hr[2] * hr[2]);
    B = -1.0 / (hr[2] * hr[2]);
    C += 2.0 / (hr[2] * hr[2]);
    robinBoundaryStencil(z, F, B, C);
}

void EllipticDomain::quadraticInterpolation(int x, int y, int z,
                                            double &W, double &E, double &S,
                                            double &N, double &F, double &B, double &C,
                                            double &scaleFactor)
{
    scaleFactor = 1.0;

    double cx = (x - (nr[0] - 1) / 2.0) * hr[0];
    double cy = (y - (nr[1] - 1) / 2.0) * hr[1];

    // since every vector for elliptic domains has ALWAYS size 2 we
    // can catch all cases manually
    double dx = 0.0;
    std::multimap<int, double>::iterator it = intersectXDir_m.find(y);
    if (cx < 0)
        ++it;
    dx = it->second;

    double dy = 0.0;
    it = intersectYDir_m.find(x);
    if (cy < 0)
        ++it;
    dy = it->second;

    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];

    W = 0.0;
    E = 0.0;
    S = 0.0;
    N = 0.0;
    C = 0.0;

    // we are a right boundary point
    if (!isInside(x + 1, y, z)) {
        double s = dx - cx;
        C -= 2.0 * (s - 1.0) / (s * de);
        E = 0.0;
        W += (s - 1.0) / ((s + 1.0) * de);
    } else {
        C += 1.0 / (de * de);
        E = -1.0 / (de * de);
    }

    // we are a left boundary point
    if (!isInside(x - 1, y, z)) {
        double s = std::abs(dx) - std::abs(cx);
        C -= 2.0 * (s - 1.0) / (s * de);
        W = 0.0;
        E += (s - 1.0) / ((s + 1.0) * de);
    } else {
        C += 1.0 / (dw * dw);
        W = -1.0 / (dw * dw);
    }

    // we are a upper boundary point
    if (!isInside(x, y + 1, z)) {
        double s = dy - cy;
        C -= 2.0 * (s - 1.0) / (s * dn);
        N = 0.0;
        S += (s - 1.0) / ((s + 1.0) * dn);
    } else {
        C += 1.0 / (dn * dn);
        N = -1.0 / (dn * dn);
    }

    // we are a lower boundary point
    if (!isInside(x, y - 1, z)) {
        double s = std::abs(dy) - std::abs(cy);
        C -= 2.0 * (s - 1.0) / (s * dn);
        S = 0.0;
        N += (s - 1.0) / ((s + 1.0) * dn);
    } else {
        C += 1.0 / (ds * ds);
        S = -1.0 / (ds * ds);
    }

    // handle boundary condition in z direction
    F = -1.0 / (hr[2] * hr[2]);
    B = -1.0 / (hr[2] * hr[2]);
    C += 2.0 / (hr[2] * hr[2]);
    robinBoundaryStencil(z, F, B, C);
}

void EllipticDomain::robinBoundaryStencil(int z, double &F, double &B, double &C) {
    if (z == 0 || z == nr[2] - 1) {

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
        double d = 0.5 * hr[2] * (nr[2] - 1);
        C += 2.0 / (d * hr[2]);
    }
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
