//
// Class ArbitraryDomain
//   Interface to iterative solver and boundary geometry
//   for space charge calculation
//
// Copyright (c) 2008,        Yves Ineichen, ETH Zürich,
//               2013 - 2015, Tülin Kaman, Paul Scherrer Institut, Villigen PSI, Switzerland
//                      2016, Daniel Winklehner, Massachusetts Institute of Technology
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
#ifndef ARBITRARY_DOMAIN_H
#define ARBITRARY_DOMAIN_H

#include <mpi.h>
#include <hdf5.h>
#include "H5hut.h"

#include <map>
#include <string>
#include <tuple>
#include <vector>
#include "IrregularDomain.h"

class BoundaryGeometry;

class ArbitraryDomain : public IrregularDomain {

public:
    ArbitraryDomain(BoundaryGeometry *bgeom, IntVector_t nr, Vector_t hr,
                    std::string interpl);

    ~ArbitraryDomain();

    /// queries if a given (x,y,z) coordinate lies inside the domain
    bool isInside(int idx, int idy, int idz) const override {
        return isInsideMap_m.at(toCoordIdx(idx, idy, idz));
    }

    // calculates intersection with rotated and shifted geometry
    void compute(Vector_t hr, NDIndex<3> localId) override;

private:
    BoundaryGeometry *bgeom_m;

    /** PointList_t maps from an (x,z) resp. (y,z) pair to double values
     * (=intersections with boundary)
     */
    typedef std::multimap< std::tuple<int, int, int>, double > PointList_t;

    /// all intersection points with gridlines in X direction
    PointList_t intersectHiX_m, intersectLoX_m;

    /// all intersection points with gridlines in Y direction
    PointList_t intersectHiY_m, intersectLoY_m;

    /// all intersection points with gridlines in Z direction
    PointList_t intersectHiZ_m, intersectLoZ_m;

    // Here we store the number of nodes in a xy layer for a given z coordinate
    std::map<int, int> numXY_m;

    // Mapping all cells that are inside the geometry
    std::map<int, bool> isInsideMap_m;

    Vector_t globalInsideP0_m;

    // Conversion from (x,y,z) to index in xyz plane
    int toCoordIdx(int idx, int idy, int idz) const {
        return (idz * nr_m[1] + idy) * nr_m[0]  + idx;
    }

    // Conversion from (x,y,z) to index on the 3D grid
    int indexAccess(int x, int y, int z) const override {
        return idxMap_m.at(toCoordIdx(x, y, z));
    }

    int coordAccess(int idx) const override {
        return coordMap_m.at(idx);
    }

    // Different interpolation methods for boundary points
    void constantInterpolation(int idx, int idy, int idz, StencilValue_t& value,
                               double &scaleFactor) const override;

    void linearInterpolation(int idx, int idy, int idz, StencilValue_t& value,
                             double &scaleFactor) const override;
};

#endif
