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
#ifndef ELLIPTICAL_DOMAIN_H
#define ELLIPTICAL_DOMAIN_H

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "Solvers/RegularDomain.h"
#include "Structure/BoundaryGeometry.h"
#include "Utilities/OpalException.h"

class EllipticDomain : public RegularDomain {

public:
    EllipticDomain(BoundaryGeometry *bgeom, IntVector_t nr,
                   Vector_t hr, std::string interpl);

    ~EllipticDomain();

    /// queries if a given (x,y,z) coordinate lies inside the domain
    bool isInside(int x, int y, int z) const override {
        double xx = getXRangeMin() + hr_m[0] * (x + 0.5);
        double yy = getYRangeMin() + hr_m[1] * (y + 0.5);

        bool isInsideEllipse = (xx * xx / (getXRangeMax() * getXRangeMax()) +
                                yy * yy / (getYRangeMax() * getYRangeMax()) < 1);

        return (isInsideEllipse && z >= 0 && z < nr_m[2]);
    }

    /// calculates intersection
    void compute(Vector_t hr, NDIndex<3> localId) override;

private:

    /// Map from a single coordinate (x or y) to a list of intersection values with
    /// boundary.
    typedef std::multimap<int, double> EllipticPointList_t;

    /// all intersection points with grid lines in X direction
    EllipticPointList_t intersectXDir_m;

    /// all intersection points with grid lines in Y direction
    EllipticPointList_t intersectYDir_m;

    /// conversion from (x,y) to index in xy plane
    int toCoordIdx(int x, int y) const { return y * nr_m[0] + x; }

    /// conversion from (x,y,z) to index on the 3D grid
    int indexAccess(int x, int y, int z) const override {
        return idxMap_m.at(toCoordIdx(x, y)) + z * getNumXY();
    }

    int coordAccess(int idx) const override {
        int ixy = idx % getNumXY();
        return coordMap_m.at(ixy);
    }

    /// different interpolation methods for boundary points
    void linearInterpolation(int x, int y, int z, StencilValue_t& value,
                             double &scaleFactor) const override;

    void quadraticInterpolation(int x, int y, int z, StencilValue_t& value,
                                double &scaleFactor) const override;
};

#endif //#ifdef ELLIPTICAL_DOMAIN_H
