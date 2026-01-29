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
#ifndef IRREGULAR_DOMAIN_H
#define IRREGULAR_DOMAIN_H

#include <map>
#include <string>
#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"
#include "Index/NDIndex.h"

/// enumeration corresponding to different interpolation methods at the boundary
enum {
    CONSTANT,
    LINEAR,
    QUADRATIC
};


class IrregularDomain {

public:

    template<typename T>
    struct Stencil {
        T center;   // x,   y,   z
        T west;     // x-1, y,   z
        T east;     // x+1, y,   z
        T north;    // x,   y+1, z
        T south;    // x,   y-1, z
        T front;    // x,   y,   z-1
        T back;     // x,   y,   z+1
    };

    typedef Stencil<int>    StencilIndex_t;
    typedef Stencil<double> StencilValue_t;
    typedef Vektor<int, 3>  IntVector_t;

    IrregularDomain(const IntVector_t& nr,
                    const Vector_t& hr,
                    const std::string& interpl);


    /** method to compute the intersection points with the boundary geometry
     * (stored in some appropriate data structure)
     * \param hr updated mesh spacings
     */
    virtual void compute(Vector_t hr, NDIndex<3> localId) = 0;

    /// method to calculate the stencil at a boundary points
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \param values of stencil element
    /// \param scaleFactor of stencil values
    void getBoundaryStencil(int x, int y, int z,
                            StencilValue_t& value,
                            double &scaleFactor) const;

    /// method to calculate the stencil at a boundary points
    /// \param id index of the current element in the matrix
    // \param values of stencil element
    /// \param scaleFactor of stencil values
    void getBoundaryStencil(int id, StencilValue_t& value,
                            double &scaleFactor) const;

    /// method to calculate the neighbours in the matrix of the current index (x,y,z)
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \param index stencil indices of an element
    void getNeighbours(int x, int y, int z, StencilIndex_t& index) const;

    void getNeighbours(int idx, StencilIndex_t& index) const;

    /// Conversion from a 3D index to (x,y,z)
    virtual void getCoord(int idx, int &x, int &y, int &z) const;

    /// Conversion from (x,y,z) to index on the 3D grid
    int getIdx(int x, int y, int z) const;

    virtual int getNumXY() const { return nr_m[0] * nr_m[1]; }

    /// method that checks if a given point lies inside the boundary
    /// \param x index of the current element in the matrix
    /// \param y index of the current element in the matrix
    /// \param z index of the current element in the matrix
    /// \return boolean indicating if the point lies inside the boundary
    virtual bool isInside(int x, int y, int z)  const = 0;

    IntVector_t getNr() const { return nr_m; }
    Vector_t    getHr() const { return hr_m; }

    void setNr(IntVector_t nr) { nr_m = nr; }
    void setHr(Vector_t hr)    { hr_m = hr; }

    void setMinMaxZ(double minz, double maxz) {
        zMin_m = minz;
        zMax_m = maxz;
    }

    double getMinZ() const { return zMin_m; }
    double getMaxZ() const { return zMax_m; }

    double getXRangeMin() const { return min_m(0); }
    double getXRangeMax() const { return max_m(0); }
    double getYRangeMin() const { return min_m(1); }
    double getYRangeMax() const { return max_m(1); }
    double getZRangeMin() const { return min_m(2); }
    double getZRangeMax() const { return max_m(2); }

    void setRangeMin(const Vector_t& min) { min_m = min; }
    void setRangeMax(const Vector_t& max) { max_m = max; }

    bool hasGeometryChanged() const { return hasGeometryChanged_m; }

    virtual ~IrregularDomain() {};

    virtual void resizeMesh(Vector_t& origin, Vector_t& hr,
                            const Vector_t& /*rmin*/, const Vector_t& /*rmax*/,
                            double /*dh*/);

protected:

    virtual int indexAccess(int x, int y, int z) const = 0;

    virtual int coordAccess(int idx) const = 0;

    /// different interpolation methods for boundary points
    virtual void constantInterpolation(int x, int y, int z, StencilValue_t& value,
                               double &scaleFactor) const;

    virtual void linearInterpolation(int x, int y, int z, StencilValue_t& value,
                             double &scaleFactor) const;

    virtual void quadraticInterpolation(int x, int y, int z, StencilValue_t& value,
                                double &scaleFactor) const;


    // a irregular domain is always defined on a grid
    /// number of mesh points in each direction
    IntVector_t nr_m;
    /// mesh-spacings in each direction
    Vector_t hr_m;

    /// min/max of bunch in floor coordinates
    double zMin_m;
    double zMax_m;

    Vector_t min_m;
    Vector_t max_m;

    /// flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    /// interpolation type
    int interpolationMethod_m;

    /// mapping (x,y,z) -> idx
    std::map<int, int> idxMap_m;

    /// mapping idx -> (x,y,z)
    std::map<int, int> coordMap_m;

};

#endif
