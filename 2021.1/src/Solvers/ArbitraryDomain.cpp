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

//#define DEBUG_INTERSECT_RAY_BOUNDARY

#include "Solvers/ArbitraryDomain.h"
#include "Structure/BoundaryGeometry.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include "Utilities/OpalException.h"
#include "Index/NDIndex.h"

ArbitraryDomain::ArbitraryDomain( BoundaryGeometry * bgeom,
                                  IntVector_t nr,
                                  Vector_t hr,
                                  std::string interpl)
    : IrregularDomain(nr, hr, interpl)
{
    bgeom_m  = bgeom;

    setRangeMin(bgeom->getmincoords());
    setRangeMax(bgeom->getmaxcoords());

    bool have_inside_pt = bgeom->getInsidePoint(globalInsideP0_m);
    if (have_inside_pt == false) {
        throw OpalException(
            "ArbitraryDomain::ArbitraryDomain()",
            "No point inside geometry found/set!");
    }

    throw OpalException("ArbitraryDomain::ArbitraryDomain()",
                        "This domain is currently not available.");
}

ArbitraryDomain::~ArbitraryDomain() {
    //nothing so far
}

void ArbitraryDomain::compute(Vector_t hr, NDIndex<3> localId){

    INFOMSG(level2 << "* Starting the Boundary Intersection Tests..." << endl);

    setHr(hr);

    int zGhostOffsetLeft  = (localId[2].first()== 0) ? 0 : 1;
    int zGhostOffsetRight = (localId[2].last() == nr_m[2] - 1) ? 0 : 1;
    int yGhostOffsetLeft  = (localId[1].first()== 0) ? 0 : 1;
    int yGhostOffsetRight = (localId[1].last() == nr_m[1] - 1) ? 0 : 1;
    int xGhostOffsetLeft  = (localId[0].first()== 0) ? 0 : 1;
    int xGhostOffsetRight = (localId[0].last() == nr_m[0] - 1) ? 0 : 1;

    hasGeometryChanged_m = true;

    intersectLoX_m.clear();
    intersectHiX_m.clear();
    intersectLoY_m.clear();
    intersectHiY_m.clear();
    intersectLoZ_m.clear();
    intersectHiZ_m.clear();

    // Calculate intersection
    Vector_t P, dir, I;
    Vector_t P0 = globalInsideP0_m;

    // We cannot assume that the geometry is symmetric about the xy, xz, and yz planes!
    // In my spiral inflector simulation, this is not the case for z direction for
    // example (-0.13 to +0.025). -DW
    for (int idz = localId[2].first()-zGhostOffsetLeft; idz <= localId[2].last()+zGhostOffsetRight; idz++) {

        P[2] = getZRangeMin() + (idz + 0.5) * hr[2];

        for (int idy = localId[1].first()-yGhostOffsetLeft; idy <= localId[1].last()+yGhostOffsetRight; idy++) {

            P[1] = getYRangeMin() + (idy + 0.5) * hr[1];

            for (int idx = localId[0].first()-xGhostOffsetLeft; idx <= localId[0].last()+xGhostOffsetRight; idx++) {

                P[0] = getXRangeMin() + (idx + 0.5) * hr[0];

                if (bgeom_m->fastIsInside(P0, P) % 2 == 0) {

                    // Fill the map with true or false values for very fast isInside tests
                    // during the rest of the fieldsolve.
                    isInsideMap_m[toCoordIdx(idx, idy, idz)] = true;

                    // Replace the old reference point with the new point (which we know is
                    // inside because we just tested for it. This makes the algorithm faster
                    // because fastIsInside() creates a number of segments that depends on the
                    // distance between P and P0. Using the previous P as the new P0
                    // assures the smallest possible distance in most cases. -DW
                    P0 = P;

                    std::tuple<int, int, int> pos(idx, idy, idz);

                    dir = Vector_t(0, 0, 1);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        intersectHiZ_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=+1 " << dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        intersectLoZ_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=-1 " << -dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    dir = Vector_t(0, 1, 0);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        intersectHiY_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=+1 " << dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        intersectLoY_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=-1" << -dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    dir = Vector_t(1, 0, 0);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        intersectHiX_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=+1 " << dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)){
                        intersectLoX_m.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=-1 " << -dir << " x,y,z= " << idx << "," << idy
                              << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }
                } else {
                    isInsideMap_m[toCoordIdx(idx, idy, idz)] = false;
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                    *gmsg << "OUTSIDE" << " x,y,z= " << idx << "," << idy
                          << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                }
            }
        }
    }

    INFOMSG(level2 << "* Finding number of ghost nodes to the left..." << endl);

    //number of ghost nodes to the left
    int numGhostNodesLeft = 0;
    if(localId[2].first() != 0) {
        for(int idx = 0; idx < nr_m[0]; idx++) {
            for(int idy = 0; idy < nr_m[1]; idy++) {
                if(isInside(idx, idy, localId[2].first() - zGhostOffsetLeft))
                    numGhostNodesLeft++;
            }
        }
    }

    INFOMSG(level2 << "* Finding number of xy points in each plane along z..." << endl);

    //xy points in z plane
    int numtotal = 0;
    numXY_m.clear();
    for(int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        int numxy = 0;
        for(int idx = 0; idx < nr_m[0]; idx++) {
            for(int idy = 0; idy < nr_m[1]; idy++) {
                if(isInside(idx, idy, idz))
                    numxy++;
            }
        }
        numXY_m[idz-localId[2].first()] = numxy;
        numtotal += numxy;
    }

    int startIdx = 0;
    MPI_Scan(&numtotal, &startIdx, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    startIdx -= numtotal;

    // Build up index and coord map
    idxMap_m.clear();
    coordMap_m.clear();
    int index = startIdx - numGhostNodesLeft;

    INFOMSG(level2 << "* Building up index and coordinate map..." << endl);

    for(int idz = localId[2].first() - zGhostOffsetLeft; idz <= localId[2].last() + zGhostOffsetRight; idz++) {
        for(int idy = 0; idy < nr_m[1]; idy++) {
            for(int idx = 0; idx < nr_m[0]; idx++) {
                if(isInside(idx, idy, idz)) {
                    idxMap_m[toCoordIdx(idx, idy, idz)] = index;
                    coordMap_m[index] = toCoordIdx(idx, idy, idz);
                    index++;
                }
            }
        }
    }

    INFOMSG(level2 << "* Done." << endl);
}

void ArbitraryDomain::constantInterpolation(int idx, int idy, int idz,
                                            StencilValue_t& value, double& /*scaleFactor*/) const
{
    value.west = -1/(hr_m[0]*hr_m[0]);
    value.east = -1/(hr_m[0]*hr_m[0]);
    value.north = -1/(hr_m[1]*hr_m[1]);
    value.south = -1/(hr_m[1]*hr_m[1]);
    value.front = -1/(hr_m[2]*hr_m[2]);
    value.back = -1/(hr_m[2]*hr_m[2]);
    value.center = 2/(hr_m[0]*hr_m[0]) + 2/(hr_m[1]*hr_m[1]) + 2/(hr_m[2]*hr_m[2]);

    if(!isInside(idx-1,idy,idz))
        value.west = 0.0;
    if(!isInside(idx+1,idy,idz))
        value.east = 0.0;

    if(!isInside(idx,idy+1,idz))
        value.north = 0.0;
    if(!isInside(idx,idy-1,idz))
        value.south = 0.0;

    if(!isInside(idx,idy,idz-1))
        value.front = 0.0;
    if(!isInside(idx,idy,idz+1))
        value.back = 0.0;
}

void ArbitraryDomain::linearInterpolation(int idx, int idy, int idz,
                                          StencilValue_t& value, double &scaleFactor) const
{
    scaleFactor = 1;

    double cx = (idx - (nr_m[0]-1)/2.0)*hr_m[0];
    double cy = (idy - (nr_m[1]-1)/2.0)*hr_m[1];
    double cz = (idz - (nr_m[2]-1)/2.0)*hr_m[2];

    double dx_w=hr_m[0];
    double dx_e=hr_m[0];
    double dy_n=hr_m[1];
    double dy_s=hr_m[1];
    double dz_f=hr_m[2];
    double dz_b=hr_m[2];
    value.center = 0.0;

    std::tuple<int, int, int> coordxyz(idx, idy, idz);

    if (idx == nr_m[0]-1)
        dx_e = std::abs(intersectHiX_m.find(coordxyz)->second - cx);
    if (idx == 0)
        dx_w = std::abs(intersectLoX_m.find(coordxyz)->second - cx);
    if (idy == nr_m[1]-1)
        dy_n = std::abs(intersectHiY_m.find(coordxyz)->second - cy);
    if (idy == 0)
        dy_s = std::abs(intersectLoY_m.find(coordxyz)->second - cy);
    if (idz == nr_m[2]-1)
        dz_b = std::abs(intersectHiZ_m.find(coordxyz)->second - cz);
    if (idz == 0)
        dz_f = std::abs(intersectLoZ_m.find(coordxyz)->second - cz);

    if(dx_w != 0)
        value.west = -(dz_f + dz_b) * (dy_n + dy_s) / dx_w;
    else
        value.west = 0;
    if(dx_e != 0)
        value.east = -(dz_f + dz_b) * (dy_n + dy_s) / dx_e;
    else
        value.east = 0;
    if(dy_n != 0)
        value.north = -(dz_f + dz_b) * (dx_w + dx_e) / dy_n;
    else
        value.north = 0;
    if(dy_s != 0)
        value.south = -(dz_f + dz_b) * (dx_w + dx_e) / dy_s;
    else
        value.south = 0;
    if(dz_f != 0)
        value.front = -(dx_w + dx_e) * (dy_n + dy_s) / dz_f;
    else
        value.front = 0;
    if(dz_b != 0)
        value.back = -(dx_w + dx_e) * (dy_n + dy_s) / dz_b;
    else
        value.back = 0;

    //RHS scaleFactor for current 3D index
    //0.5* comes from discretiztaion
    //scaleFactor = 0.5*(dw+de)*(dn+ds)*(df+db);
    scaleFactor = 0.5;
    if(dx_w + dx_e != 0)
        scaleFactor *= (dx_w + dx_e);
    if(dy_n + dy_s != 0)
        scaleFactor *= (dy_n + dy_s);
    if(dz_f + dz_b != 0)
        scaleFactor *= (dz_f + dz_b);

    //catch the case where a point lies on the boundary
    double m1 = dx_w * dx_e;
    double m2 = dy_n * dy_s;
    if(dx_e == 0)
        m1 = dx_w;
    if(dx_w == 0)
        m1 = dx_e;
    if(dy_n == 0)
        m2 = dy_s;
    if(dy_s == 0)
        m2 = dy_n;

    value.center = 2 / hr_m[2];
    if(dx_w != 0 || dx_e != 0)
        value.center *= (dx_w + dx_e);
    if(dy_n != 0 || dy_s != 0)
        value.center *= (dy_n + dy_s);
    if(dx_w != 0 || dx_e != 0)
        value.center += (dz_f + dz_b) * (dy_n + dy_s) * (dx_w + dx_e) / m1;
    if(dy_n != 0 || dy_s != 0)
        value.center += (dz_f + dz_b) * (dx_w + dx_e) * (dy_n + dy_s) / m2;
}
