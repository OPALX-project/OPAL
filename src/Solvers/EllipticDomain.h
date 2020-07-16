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
#include "IrregularDomain.h"
#include "Structure/BoundaryGeometry.h"
#include "Utilities/OpalException.h"

class EllipticDomain : public IrregularDomain {

public:

    EllipticDomain(Vector_t nr, Vector_t hr);

    EllipticDomain(double semimajor, double semiminor, Vector_t nr,
                   Vector_t hr, std::string interpl);

    EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr,
                   Vector_t hr, std::string interpl);

    EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr, std::string interpl);

    ~EllipticDomain();

    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int x, int y, int z,
                            double &W, double &E, double &S,
                            double &N, double &F, double &B,
                            double &C, double &scaleFactor);

    /// returns discretization at 3D index
    void getBoundaryStencil(int idx, double &W, double &E,
                            double &S, double &N, double &F,
                            double &B, double &C, double &scaleFactor);

    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int x, int y, int z, int &W, int &E,
                       int &S, int &N, int &F, int &B);

    /// returns index of neighbours at 3D index
    void getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B);

    /// returns type of boundary condition
    std::string getType() {return "Elliptic";}

    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int z) {
        double xx = - semiMajor_m + hr[0] * (x + 0.5);
        double yy = - semiMinor_m + hr[1] * (y + 0.5);

        bool isInsideEllipse = (xx * xx / (semiMajor_m * semiMajor_m) +
                                yy * yy / (semiMinor_m * semiMinor_m) < 1);

        return (isInsideEllipse && z > 0 && z < nr[2] - 1);
    }

    int getNumXY(int /*z*/) { return nxy_m; }


    /// calculates intersection
    void compute(Vector_t /*hr*/) {
        throw OpalException("EllipticDomain::compute()", "This function is not available.");
    }

    void compute(Vector_t hr, NDIndex<3> localId);

    double getXRangeMin() { return -semiMajor_m; }
    double getXRangeMax() { return semiMajor_m;  }
    double getYRangeMin() { return -semiMinor_m; }
    double getYRangeMax() { return semiMinor_m;  }
    double getZRangeMin() { return zMin_m; }
    double getZRangeMax() { return zMax_m; }

    bool hasGeometryChanged() { return hasGeometryChanged_m; }


    void resizeMesh(Vector_t& origin, Vector_t& hr, const Vector_t& rmin,
                    const Vector_t& rmax, double dh);

private:

    /// Map from a single coordinate (x or y) to a list of intersection values with
    /// boundary.
    typedef std::multimap<int, double> EllipticPointList_t;

    /// all intersection points with grid lines in X direction
    EllipticPointList_t intersectXDir_m;

    /// all intersection points with grid lines in Y direction
    EllipticPointList_t intersectYDir_m;

    /// mapping (x,y,z) -> idx
    std::map<int, int> idxMap_m;

    /// mapping idx -> (x,y,z)
    std::map<int, int> coordMap_m;

    /// semi-major of the ellipse
    double semiMajor_m;

    /// semi-minor of the ellipse
    double semiMinor_m;

    /// number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;

    /// interpolation type
    int interpolationMethod_m;

    /// flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    /// conversion from (x,y) to index in xy plane
    inline int toCoordIdx(int x, int y) { return y * nr[0] + x; }

    /// conversion from (x,y,z) to index on the 3D grid
    inline int getIdx(int x, int y, int z) {
        if (isInside(x, y, z))
            return idxMap_m[toCoordIdx(x, y)] + (z - 1) * nxy_m;
        else
            return -1;
    }

    /// conversion from a 3D index to (x,y,z)
    inline void getCoord(int idx, int &x, int &y, int &z) {
        int ixy = idx % nxy_m;
        int xy = coordMap_m[ixy];
        int inr = (int)nr[0];
        x = xy % inr;
        y = (xy - x) / nr[0];
        z = (idx - ixy) / nxy_m + 1;
    }

    /// different interpolation methods for boundary points
    void constantInterpolation(int x, int y, int z,
                               double &W, double &E, double &S,
                               double &N, double &F, double &B,
                               double &C, double &scaleFactor);

    void linearInterpolation(int x, int y, int z, double &W,
                             double &E, double &S, double &N,
                             double &F, double &B, double &C,
                             double &scaleFactor);

    void quadraticInterpolation(int x, int y, int z, double &W,
                                double &E, double &S, double &N,
                                double &F, double &B, double &C,
                                double &scaleFactor);

    /// function to handle the open boundary condition in longitudinal direction
    void robinBoundaryStencil(int z, double &F, double &B, double &C);
};

#endif //#ifdef ELLIPTICAL_DOMAIN_H

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
