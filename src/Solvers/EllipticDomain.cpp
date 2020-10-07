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

//FIXME: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)


EllipticDomain::EllipticDomain(BoundaryGeometry *bgeom, IntVector_t nr, Vector_t hr,
                               std::string interpl)
    : RegularDomain(nr, hr, interpl)
{
    Vector_t min(-bgeom->getA(), -bgeom->getB(), bgeom->getS());
    Vector_t max( bgeom->getA(),  bgeom->getB(), bgeom->getS() + bgeom->getLength());
    setRangeMin(min);
    setRangeMax(max);
    setMinMaxZ(min[2], max[2]);
}

EllipticDomain::~EllipticDomain() {
    //nothing so far
}

// for this geometry we only have to calculate the intersection on
// one x-y-plane
// for the moment we center the ellipse around the center of the grid
void EllipticDomain::compute(Vector_t hr, NDIndex<3> localId) {
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

    // clear previous coordinate maps
    idxMap_m.clear();
    coordMap_m.clear();
    //clear previous intersection points
    intersectYDir_m.clear();
    intersectXDir_m.clear();

    // build a index and coordinate map
    int idx = 0;
    int x, y;

    int nxy = 0;

    /* we need to scan the full x-y-plane on all cores
     * in order to figure out the number of valid
     * grid points per plane --> otherwise we might
     * get not unique global indices in the Tpetra::CrsMatrix
     */
    for (x = 0; x < nr_m[0]; ++x) {
        for (y = 0; y < nr_m[1]; ++y) {
            if (isInside(x, y, 1)) {
                idxMap_m[toCoordIdx(x, y)] = idx;
                coordMap_m[idx++] = toCoordIdx(x, y);
                nxy++;
            }
        }
    }

    setNumXY(nxy);

    switch (interpolationMethod_m) {
        case CONSTANT:
            break;
        case LINEAR:
        case QUADRATIC:

            double smajsq = getXRangeMax() * getXRangeMax();
            double sminsq = getYRangeMax() * getYRangeMax();
            double yd = 0.0;
            double xd = 0.0;
            double pos = 0.0;

            // calculate intersection with the ellipse
            for (x = localId[0].first(); x <= localId[0].last(); x++) {
                pos = getXRangeMin() + hr_m[0] * (x + 0.5);
                if (std::abs(pos) >= getXRangeMax())
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
                pos = getYRangeMin() + hr_m[1] * (y + 0.5);
                if (std::abs(pos) >= getYRangeMax())
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

void EllipticDomain::linearInterpolation(int x, int y, int z, StencilValue_t& value,
                                         double &scaleFactor) const
{
    scaleFactor = 1.0;

    double cx = getXRangeMin() + hr_m[0] * (x + 0.5);
    double cy = getYRangeMin() + hr_m[1] * (y + 0.5);

    double dx = 0.0;
    std::multimap<int, double>::const_iterator it = intersectXDir_m.find(y);

    if (cx < 0)
        ++it;
    dx = it->second;

    double dy = 0.0;
    it = intersectYDir_m.find(x);
    if (cy < 0)
        ++it;
    dy = it->second;


    double dw = hr_m[0];
    double de = hr_m[0];
    double dn = hr_m[1];
    double ds = hr_m[1];
    value.center = 0.0;

    // we are a right boundary point
    if (!isInside(x + 1, y, z)) {
        value.center += 1.0 / ((dx - cx) * de);
        value.east = 0.0;
    } else {
        value.center += 1.0 / (de * de);
        value.east = -1.0 / (de * de);
    }

    // we are a left boundary point
    if (!isInside(x - 1, y, z)) {
        value.center += 1.0 / ((std::abs(dx) - std::abs(cx)) * dw);
        value.west = 0.0;
    } else {
        value.center += 1.0 / (dw * dw);
        value.west = -1.0 / (dw * dw);
    }

    // we are a upper boundary point
    if (!isInside(x, y + 1, z)) {
        value.center += 1.0 / ((dy - cy) * dn);
        value.north = 0.0;
    } else {
        value.center += 1.0 / (dn * dn);
        value.north = -1.0 / (dn * dn);
    }

    // we are a lower boundary point
    if (!isInside(x, y - 1, z)) {
        value.center += 1.0 / ((std::abs(dy) - std::abs(cy)) * ds);
        value.south = 0.0;
    } else {
        value.center += 1.0 / (ds * ds);
        value.south = -1.0 / (ds * ds);
    }

    // handle boundary condition in z direction
    value.front = -1.0 / (hr_m[2] * hr_m[2]);
    value.back  = -1.0 / (hr_m[2] * hr_m[2]);
    value.center += 2.0 / (hr_m[2] * hr_m[2]);
    robinBoundaryStencil(z, value.front, value.back, value.center);
}

void EllipticDomain::quadraticInterpolation(int x, int y, int z,
                                            StencilValue_t& value,
                                            double &scaleFactor) const
{
    scaleFactor = 1.0;

    double cx = (x - (nr_m[0] - 1) / 2.0) * hr_m[0];
    double cy = (y - (nr_m[1] - 1) / 2.0) * hr_m[1];

    // since every vector for elliptic domains has ALWAYS size 2 we
    // can catch all cases manually
    double dx = 0.0;
    std::multimap<int, double>::const_iterator it = intersectXDir_m.find(y);
    if (cx < 0)
        ++it;
    dx = it->second;

    double dy = 0.0;
    it = intersectYDir_m.find(x);
    if (cy < 0)
        ++it;
    dy = it->second;

    double dw = hr_m[0];
    double de = hr_m[0];
    double dn = hr_m[1];
    double ds = hr_m[1];

    value.west = 0.0;
    value.east = 0.0;
    value.south = 0.0;
    value.north = 0.0;
    value.center = 0.0;

    // we are a right boundary point
    if (!isInside(x + 1, y, z)) {
        double s = dx - cx;
        value.center -= 2.0 * (s - 1.0) / (s * de);
        value.east = 0.0;
        value.west += (s - 1.0) / ((s + 1.0) * de);
    } else {
        value.center += 1.0 / (de * de);
        value.east = -1.0 / (de * de);
    }

    // we are a left boundary point
    if (!isInside(x - 1, y, z)) {
        double s = std::abs(dx) - std::abs(cx);
        value.center -= 2.0 * (s - 1.0) / (s * de);
        value.west = 0.0;
        value.east += (s - 1.0) / ((s + 1.0) * de);
    } else {
        value.center += 1.0 / (dw * dw);
        value.west = -1.0 / (dw * dw);
    }

    // we are a upper boundary point
    if (!isInside(x, y + 1, z)) {
        double s = dy - cy;
        value.center -= 2.0 * (s - 1.0) / (s * dn);
        value.north = 0.0;
        value.south += (s - 1.0) / ((s + 1.0) * dn);
    } else {
        value.center += 1.0 / (dn * dn);
        value.north = -1.0 / (dn * dn);
    }

    // we are a lower boundary point
    if (!isInside(x, y - 1, z)) {
        double s = std::abs(dy) - std::abs(cy);
        value.center -= 2.0 * (s - 1.0) / (s * dn);
        value.south = 0.0;
        value.north += (s - 1.0) / ((s + 1.0) * dn);
    } else {
        value.center += 1.0 / (ds * ds);
        value.south = -1.0 / (ds * ds);
    }

    // handle boundary condition in z direction
    value.front = -1.0 / (hr_m[2] * hr_m[2]);
    value.back  = -1.0 / (hr_m[2] * hr_m[2]);
    value.center += 2.0 / (hr_m[2] * hr_m[2]);
    robinBoundaryStencil(z, value.front, value.back, value.center);
}
