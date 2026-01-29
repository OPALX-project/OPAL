//
// Class BoxCornerDomain
//   :FIXME: add brief description
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
#include "Solvers/BoxCornerDomain.h"
#include "Utilities/OpalException.h"

#include <map>
#include <string>
#include <cmath>
#include <iostream>

//FIXME: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)

BoxCornerDomain::BoxCornerDomain(double A, double B, double C,
                                 double L1, double L2, IntVector_t nr, Vector_t hr,
                                 std::string interpl)
    : RegularDomain(nr, hr, interpl)
{
    setRangeMin(Vector_t(-A, -B, L1));
    setRangeMax(Vector_t( A,  B, L1 + L2));
    C_m = C;

    throw OpalException("BoxCornerDomain::BoxCornerDomain()",
                        "This domain is currently not supported!");
}

BoxCornerDomain::~BoxCornerDomain() {
    //nothing so far
}


// for this geometry we only have to calculate the intersection on
// all x-y-planes
// for the moment we center the box corner geometry around the center of the grid
// hr holds the grid-spacings (boundary ellipse embedded in hr-grid)

void BoxCornerDomain::compute(Vector_t hr, NDIndex<3> /*localId*/){

    //there is nothing to be done if the mesh spacings have not changed
    //    if(hr_m[0] == getHr()[0] && hr_m[1] == getHr()[1] && hr_m[2] == getHr()[2]) {
    //      hasGeometryChanged_m = false;
    //      return;
    //  }

    setHr(hr);
    hasGeometryChanged_m = true;

    double bL= getB(getMinZ());
    double bH= getB(getMaxZ());

    actBMin_m = getYRangeMin();
    actBMax_m = std::max(bL,bH);

    //reset number of points inside domain

    // clear previous coordinate maps
    idxMap_m.clear();
    coordMap_m.clear();
    //clear previous intersection points
    IntersectYDir.clear();
    IntersectXDir.clear();

    // build a index and coordinate map
    int idx = 0;
    int x, y, z;
    for(x = 0; x < nr_m[0]; x++) {
        for(y = 0; y < nr_m[1]; y++) {
            for(z = 0; z < nr_m[2]; z++) {
                if(isInside(x, y, z)) {
                    idxMap_m[toCoordIdx(x, y, z)] = idx;
                    coordMap_m[idx++] = toCoordIdx(x, y, z);
                }
            }
        }
    }

    //XXX: calculate intersection on the fly
    /*
    switch(interpolationMethod_m) {

    case CONSTANT:
        break;
    case LINEAR:
    case QUADRATIC:

        // calculate intersection

        for(int z = 0; z < nr_m[2]; z++) {

            for(int x = 0; x < nr_m[0]; x++) {
                // the x coordinate does not change in the CornerBox geometry
                std::pair<int, int> pos(x, z);
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, 0.5*getXRangeMax()));
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, 0.5*getXRangeMin()));
            }

            for(int y = 0; y < nr_m[1]; y++) {
                std::pair<int, int> pos(y, z);
                double yt = getB(z*hr_m[2]);
                double yb = 0.5*getYRangeMin();
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, yt));
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, yb));
            }
        }
    }
    */
}

void BoxCornerDomain::linearInterpolation(int x, int y, int z, StencilValue_t& value,
                                          double &scaleFactor) const
{
    scaleFactor = 1.0;

    double cx = x * hr_m[0] - (nr_m[0] - 1) * hr_m[0] / 2.0;
    double cy = y * hr_m[1] - (nr_m[1] - 1) * hr_m[1] / 2.0;

    //XXX: calculate intersection on the fly
    /*
    multimap< pair<int, int>, double >::iterator it;
    pair< multimap< pair<int, int>, double>::iterator, multimap< pair<int, int>, double>::iterator > ret;

    double dx = 0.0;
    std::pair<int, int> coordxz(x, z);
    ret = IntersectXDir.equal_range(coordxz);
    if(cx < 0)
        it++;
    dx = it->second;

    double dy = 0.0;
    std::pair<int, int> coordyz(y, z);
    ret = IntersectYDir.equal_range(coordyz);
    if(cy < 0)
        it++;
    dy = it->second;
    */

    double dw = hr_m[0];
    double de = hr_m[0];
    double dn = hr_m[1];
    double ds = hr_m[1];
    value.center = 0.0;

    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        double dx = getXIntersection(cx, z);
        value.center += 1 / ((dx - cx) * de);
        value.east = 0.0;
    } else {
        value.center += 1 / (de * de);
        value.east = -1 / (de * de);
    }

    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        double dx = getXIntersection(cx, z);
        value.center += 1 / ((std::abs(dx) - std::abs(cx)) * dw);
        value.west = 0.0;
    } else {
        value.center += 1 / (dw * dw);
        value.west = -1 / (dw * dw);
    }

    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        double dy = getYIntersection(cy, z);
        value.center += 1 / ((dy - cy) * dn);
        value.north = 0.0;
    } else {
        value.center += 1 / (dn * dn);
        value.north = -1 / (dn * dn);
    }

    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        double dy = getYIntersection(cy, z);
        value.center += 1 / ((std::abs(dy) - std::abs(cy)) * ds);
        value.south = 0.0;
    } else {
        value.center += 1 / (ds * ds);
        value.south = -1 / (ds * ds);
    }

    value.front = -1 / (hr_m[2] * hr_m[2]);
    value.back = -1 / (hr_m[2] * hr_m[2]);
    value.center += 2 / (hr_m[2] * hr_m[2]);

    // handle boundary condition in z direction
    if(z == 0 || z == nr_m[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            value.front = 0.0;
        else
            value.back = 0.0;

        //hr_m[2]*(nr_m2[2]-1)/2 = radius
        double d = hr_m[2] * (nr_m[2] - 1) / 2;
        value.center += 2 / (d * hr_m[2]);

        value.west   /= 2.0;
        value.east   /= 2.0;
        value.north  /= 2.0;
        value.south  /= 2.0;
        value.center /= 2.0;
        scaleFactor  *= 0.5;

    }

}

//FIXME: this probably needs some cleanup/rewriting
void BoxCornerDomain::quadraticInterpolation(int x, int y, int z, StencilValue_t& value,
                                             double &scaleFactor) const
{
    double cx = (x - (nr_m[0] - 1) / 2.0) * hr_m[0];
    double cy = (y - (nr_m[1] - 1) / 2.0) * hr_m[1];

    double dx = getXIntersection(cx, z);
    double dy = getYIntersection(cy, z);

    //XXX: calculate intersection on the fly
    /*
    multimap< pair<int, int>, double >::iterator it;
    pair< multimap< pair<int, int>, double>::iterator, multimap< pair<int, int>, double>::iterator > ret;

    double dx = 0.0;
    std::pair<int, int> coordxz(x, z);
    ret = IntersectXDir.equal_range(coordxz);
    if(cx < 0)
        it++;
    dx = it->second;

    double dy = 0.0;
    std::pair<int, int> coordyz(y, z);
    ret = IntersectYDir.equal_range(coordyz);
    if(cy < 0)
        it++;
    dy = it->second;
    */

    double dw = hr_m[0];
    double de = hr_m[0];
    double dn = hr_m[1];
    double ds = hr_m[1];
    value.west   = 1.0;
    value.east   = 1.0;
    value.north  = 1.0;
    value.south  = 1.0;
    value.front  = 1.0;
    value.back   = 1.0;
    value.center = 0.0;

    //TODO: = cx+hr_m[0] > dx && cx > 0
    //if((x-nr_m[0]/2.0+1)*hr_m[0] > dx && cx > 0) {
    ////we are a right boundary point
    ////if(!isInside(x+1,y,z)) {
    //de = dx-cx;
    //}

    //if((x-nr_m[0]/2.0-1)*hr_m[0] < dx && cx < 0) {
    ////we are a left boundary point
    ////if(!isInside(x-1,y,z)) {
    //dw = std::abs(dx)-std::abs(cx);
    //}

    //if((y-nr_m[1]/2.0+1)*hr_m[1] > dy && cy > 0) {
    ////we are a upper boundary point
    ////if(!isInside(x,y+1,z)) {
    //dn = dy-cy;
    //}

    //if((y-nr_m[1]/2.0-1)*hr_m[1] < dy && cy < 0) {
    ////we are a lower boundary point
    ////if(!isInside(x,y-1,z)) {
    //ds = std::abs(dy)-std::abs(cy);
    //}

    //TODO: = cx+hr_m[0] > dx && cx > 0
    //if((x-nr_m[0]/2.0+1)*hr_m[0] > dx && cx > 0) {
    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        de = dx - cx;
        value.east = 0.0;
    }

    //if((x-nr_m[0]/2.0-1)*hr_m[0] < dx && cx < 0) {
    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        dw = std::abs(dx) - std::abs(cx);
        value.west = 0.0;
    }

    //if((y-nr_m[1]/2.0+1)*hr_m[1] > dy && cy > 0) {
    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        dn = dy - cy;
        value.north = 0.0;
    }

    //if((y-nr_m[1]/2.0-1)*hr_m[1] < dy && cy < 0) {
    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        ds = std::abs(dy) - std::abs(cy);
        value.south = 0.0;
    }

    //2/dw*(dw_de)
    value.west  *= -1.0 / (dw * (dw + de));
    value.east  *= -1.0 / (de * (dw + de));
    value.north *= -1.0 / (dn * (dn + ds));
    value.south *= -1.0 / (ds * (dn + ds));
    value.front = -1 / (hr_m[2] * (hr_m[2] + hr_m[2]));
    value.back  = -1 / (hr_m[2] * (hr_m[2] + hr_m[2]));

    //TODO: problem when de,dw,dn,ds == 0
    //is NOT a regular BOUND PT
    value.center += 1 / de * 1 / dw;
    value.center += 1 / dn * 1 / ds;
    value.center += 1 / hr_m[2] * 1 / hr_m[2];
    scaleFactor = 0.5;


    //for regular gridpoints no problem with symmetry, just boundary
    //z direction is right
    //implement isLastInside(dir)
    //we have LOCAL x,y coordinates!

    /*
       if(dw != 0 && !wIsB)
       W = -1/dw * (dn+ds) * 2*hr_m[2];
       else
       W = 0;
       if(de != 0 && !eIsB)
       E = -1/de * (dn+ds) * 2*hr_m[2];
       else
       E = 0;
       if(dn != 0 && !nIsB)
       N = -1/dn * (dw+de) * 2*hr_m[2];
       else
       N = 0;
       if(ds != 0 && !sIsB)
       S = -1/ds * (dw+de) * 2*hr_m[2];
       else
       S = 0;
       F = -(dw+de)*(dn+ds)/hr_m[2];
       B = -(dw+de)*(dn+ds)/hr_m[2];
       */

    //if(dw != 0)
    //W = -2*hr_m[2]*(dn+ds)/dw;
    //else
    //W = 0;
    //if(de != 0)
    //E = -2*hr_m[2]*(dn+ds)/de;
    //else
    //E = 0;
    //if(dn != 0)
    //N = -2*hr_m[2]*(dw+de)/dn;
    //else
    //N = 0;
    //if(ds != 0)
    //S = -2*hr_m[2]*(dw+de)/ds;
    //else
    //S = 0;
    //F = -(dw+de)*(dn+ds)/hr_m[2];
    //B = -(dw+de)*(dn+ds)/hr_m[2];

    //// RHS scaleFactor for current 3D index
    //// Factor 0.5 results from the SW/quadratic extrapolation
    //scaleFactor = 0.5*(dw+de)*(dn+ds)*(2*hr_m[2]);

    // catch the case where a point lies on the boundary
    //FIXME: do this more elegant!
    //double m1 = dw*de;
    //double m2 = dn*ds;
    //if(de == 0)
    //m1 = dw;
    //if(dw == 0)
    //m1 = de;
    //if(dn == 0)
    //m2 = ds;
    //if(ds == 0)
    //m2 = dn;
    ////XXX: dn+ds || dw+de can be 0
    ////C = 2*(dn+ds)*(dw+de)/hr_m[2];
    //C = 2/hr_m[2];
    //if(dw != 0 || de != 0)
    //C *= (dw+de);
    //if(dn != 0 || ds != 0)
    //C *= (dn+ds);
    //if(dw != 0 || de != 0)
    //C += (2*hr_m[2])*(dn+ds)*(dw+de)/m1;
    //if(dn != 0 || ds != 0)
    //C += (2*hr_m[2])*(dw+de)*(dn+ds)/m2;

    //handle Neumann case
    //if(z == 0 || z == nr_m[2]-1) {

    //if(z == 0)
    //F = 0.0;
    //else
    //B = 0.0;

    ////neumann stuff
    //W = W/2.0;
    //E = E/2.0;
    //N = N/2.0;
    //S = S/2.0;
    //C /= 2.0;

    //scaleFactor /= 2.0;
    //}

    // handle boundary condition in z direction
    if(z == 0 || z == nr_m[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            value.front = 0.0;
        else
            value.back = 0.0;

        //value.center += 2/((hr_m[2]*(nr_m[2]-1)/2.0) * hr_m[2]);
        //hr_m[2]*(nr_m2[2]-1)/2 = radius
        double d = hr_m[2] * (nr_m[2] - 1) / 2;
        value.center += 2 / (d * hr_m[2]);

        value.west   /= 2.0;
        value.east   /= 2.0;
        value.north  /= 2.0;
        value.south  /= 2.0;
        value.center /= 2.0;
        scaleFactor  /= 2.0;

    }
}
