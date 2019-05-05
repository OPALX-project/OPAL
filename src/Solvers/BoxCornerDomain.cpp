#ifdef HAVE_SAAMG_SOLVER
#include "Solvers/BoxCornerDomain.h"

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <assert.h>

//FIXME: ORDER HOW TO TRAVERSE NODES IS FIXED, THIS SHOULD BE MORE GENERIC! (PLACES MARKED)

extern Inform *gmsg;

BoxCornerDomain::BoxCornerDomain(Vector_t nr, Vector_t hr) {
    setNr(nr);
    setHr(hr);
}

BoxCornerDomain::BoxCornerDomain(double A, double B, double C, double Length, double L1, double L2, Vector_t nr, Vector_t hr, std::string interpl) {
    A_m = A;
    B_m = B;
    C_m = C;
    Length_m = Length;
    L1_m = L1;
    L2_m = L2;

    setNr(nr);
    setHr(hr);

    if(interpl == "CONSTANT")
        interpolationMethod = CONSTANT;
    else if(interpl == "LINEAR")
        interpolationMethod = LINEAR;
    else if(interpl == "QUADRATIC")
        interpolationMethod = QUADRATIC;

    if(Ippl::getNodes() == 1) {
      *gmsg << " Write BoxCorner data to file boxcorner.dat" << endl;
      std::string file("boxcorner.dat");
      os_m.open(file.c_str());
      if(os_m.bad()) {
          *gmsg << "Unable to open output file " <<  file << endl;
      }
      //os_m << "# ...." << endl;
    }
}

BoxCornerDomain::~BoxCornerDomain() {
    //nothing so far
}


// for this geometry we only have to calculate the intersection on
// all x-y-planes
// for the moment we center the box corner geometry around the center of the grid
// hr holds the grid-spacings (boundary ellipse embedded in hr-grid)

void BoxCornerDomain::compute(Vector_t hr){

    //there is nothing to be done if the mesh spacings have not changed
    //    if(hr[0] == getHr()[0] && hr[1] == getHr()[1] && hr[2] == getHr()[2]) {
    //      hasGeometryChanged_m = false;
    //      return;
    //  }

    setHr(hr);
    hasGeometryChanged_m = true;

    double bL= getB(getMinZ());
    double bH= getB(getMaxZ());

    actBMin_m = -B_m;
    actBMax_m = std::max(bL,bH);

    INFOMSG(" BoxCorner L= " << Length_m << " L1= " << L1_m << " L2= " << L2_m << " A= " << A_m << " B= " << B_m << " C= " << C_m
            << " bL= " << bL << " bH= " << bH <<  " actBMin= " << actBMin_m << " actBMax=max(bL,bH)= " << actBMax_m << endl);

    //reset number of points inside domain

    // clear previous coordinate maps
    IdxMap.clear();
    CoordMap.clear();
    //clear previous intersection points
    IntersectYDir.clear();
    IntersectXDir.clear();

    // build a index and coordinate map
    int idx = 0;
    int x, y, z;
    for(x = 0; x < nr[0]; x++) {
        for(y = 0; y < nr[1]; y++) {
            for(z = 0; z < nr[2]; z++) {
                if(isInside(x, y, z)) {
                    IdxMap[toCoordIdx(x, y, z)] = idx;
                    CoordMap[idx++] = toCoordIdx(x, y, z);
                }
            }
        }
    }

    //XXX: calculate intersection on the fly
    /*
    switch(interpolationMethod) {

    case CONSTANT:
        break;
    case LINEAR:
    case QUADRATIC:

        // calculate intersection

        for(int z = 0; z < nr[2]; z++) {

            for(int x = 0; x < nr[0]; x++) {
                // the x coordinate does not change in the CornerBox geometry
                std::pair<int, int> pos(x, z);
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, 0.5*A_m));
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, -0.5*A_m));
            }

            for(int y = 0; y < nr[1]; y++) {
                std::pair<int, int> pos(y, z);
                double yt = getB(z*hr[2]);
                double yb = -0.5*B_m;
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, yt));
                IntersectXDir.insert(std::pair< std::pair<int, int>, double >(pos, yb));
            }
        }
    }
    */
}

void BoxCornerDomain::compute(Vector_t hr, NDIndex<3> localId){
}

void BoxCornerDomain::getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod) {
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

void BoxCornerDomain::getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {


    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getBoundaryStencil(x, y, z, W, E, S, N, F, B, C, scaleFactor);
}


void BoxCornerDomain::getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B) {

    int x = 0, y = 0, z = 0;
    getCoord(idx, x, y, z);
    getNeighbours(x, y, z, W, E, S, N, F, B);

}

void BoxCornerDomain::getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B) {

    if(x > 0)
        W = getIdx(x - 1, y, z);
    else
        W = -1;
    if(x < nr[0] - 1)
        E = getIdx(x + 1, y, z);
    else
        E = -1;

    if(y < nr[1] - 1)
        N = getIdx(x, y + 1, z);
    else
        N = -1;
    if(y > 0)
        S = getIdx(x, y - 1, z);
    else
        S = -1;

    if(z > 0)
        F = getIdx(x, y, z - 1);
    else
        F = -1;
    if(z < nr[2] - 1)
        B = getIdx(x, y, z + 1);
    else
        B = -1;

}

void BoxCornerDomain::constantInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;

    W = -1 / hr[0] * 1 / hr[0];
    E = -1 / hr[0] * 1 / hr[0];
    N = -1 / hr[1] * 1 / hr[1];
    S = -1 / hr[1] * 1 / hr[1];
    F = -1 / hr[2] * 1 / hr[2];
    B = -1 / hr[2] * 1 / hr[2];
    C = 2 / hr[0] * 1 / hr[0] + 2 / hr[1] * 1 / hr[1] + 2 / hr[2] * 1 / hr[2];

    // we are a right boundary point
    if(!isInside(x + 1, y, z))
        E = 0.0;

    // we are a left boundary point
    if(!isInside(x - 1, y, z))
        W = 0.0;

    // we are a upper boundary point
    if(!isInside(x, y + 1, z))
        N = 0.0;

    // we are a lower boundary point
    if(!isInside(x, y - 1, z))
        S = 0.0;

    if(z == 1 || z == nr[2] - 2) {

        // case where we are on the Robin BC in Z-direction
        // where we distinguish two cases
        // IFF: this values should not matter because they
        // never make it into the discretization matrix
        if(z == 1)
            F = 0.0;
        else
            B = 0.0;

        // add contribution of Robin discretization to center point
        // d the distance between the center of the bunch and the boundary
        //double cx = (x-(nr[0]-1)/2)*hr[0];
        //double cy = (y-(nr[1]-1)/2)*hr[1];
        //double cz = hr[2]*(nr[2]-1);
        //double d = sqrt(cx*cx+cy*cy+cz*cz);
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);
        //C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);

        // scale all stencil-points in z-plane with 0.5 (Robin discretization)
        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;
    }

}

void BoxCornerDomain::linearInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;

    double cx = x * hr[0] - (nr[0] - 1) * hr[0] / 2.0;
    double cy = y * hr[1] - (nr[1] - 1) * hr[1] / 2.0;

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

    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    C = 0.0;

    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        double dx = getXIntersection(cx, z);
        C += 1 / ((dx - cx) * de);
        E = 0.0;
    } else {
        C += 1 / (de * de);
        E = -1 / (de * de);
    }

    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        double dx = getXIntersection(cx, z);
        C += 1 / ((std::abs(dx) - std::abs(cx)) * dw);
        W = 0.0;
    } else {
        C += 1 / (dw * dw);
        W = -1 / (dw * dw);
    }

    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        double dy = getYIntersection(cy, z);
        C += 1 / ((dy - cy) * dn);
        N = 0.0;
    } else {
        C += 1 / (dn * dn);
        N = -1 / (dn * dn);
    }

    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        double dy = getYIntersection(cy, z);
        C += 1 / ((std::abs(dy) - std::abs(cy)) * ds);
        S = 0.0;
    } else {
        C += 1 / (ds * ds);
        S = -1 / (ds * ds);
    }

    F = -1 / (hr[2] * hr[2]);
    B = -1 / (hr[2] * hr[2]);
    C += 2 / (hr[2] * hr[2]);

    // handle boundary condition in z direction
    if(z == 0 || z == nr[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        //hr[2]*(nr2[2]-1)/2 = radius
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);

        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor *= 0.5;

    }

}

//FIXME: this probably needs some cleanup/rewriting
void BoxCornerDomain::quadraticInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    double cx = (x - (nr[0] - 1) / 2.0) * hr[0];
    double cy = (y - (nr[1] - 1) / 2.0) * hr[1];

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

    double dw = hr[0];
    double de = hr[0];
    double dn = hr[1];
    double ds = hr[1];
    W = 1.0;
    E = 1.0;
    N = 1.0;
    S = 1.0;
    F = 1.0;
    B = 1.0;
    C = 0.0;

    //TODO: = cx+hr[0] > dx && cx > 0
    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    ////we are a right boundary point
    ////if(!isInside(x+1,y,z)) {
    //de = dx-cx;
    //}

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    ////we are a left boundary point
    ////if(!isInside(x-1,y,z)) {
    //dw = std::abs(dx)-std::abs(cx);
    //}

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    ////we are a upper boundary point
    ////if(!isInside(x,y+1,z)) {
    //dn = dy-cy;
    //}

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    ////we are a lower boundary point
    ////if(!isInside(x,y-1,z)) {
    //ds = std::abs(dy)-std::abs(cy);
    //}

    //TODO: = cx+hr[0] > dx && cx > 0
    //if((x-nr[0]/2.0+1)*hr[0] > dx && cx > 0) {
    //we are a right boundary point
    if(!isInside(x + 1, y, z)) {
        de = dx - cx;
        E = 0.0;
    }

    //if((x-nr[0]/2.0-1)*hr[0] < dx && cx < 0) {
    //we are a left boundary point
    if(!isInside(x - 1, y, z)) {
        dw = std::abs(dx) - std::abs(cx);
        W = 0.0;
    }

    //if((y-nr[1]/2.0+1)*hr[1] > dy && cy > 0) {
    //we are a upper boundary point
    if(!isInside(x, y + 1, z)) {
        dn = dy - cy;
        N = 0.0;
    }

    //if((y-nr[1]/2.0-1)*hr[1] < dy && cy < 0) {
    //we are a lower boundary point
    if(!isInside(x, y - 1, z)) {
        ds = std::abs(dy) - std::abs(cy);
        S = 0.0;
    }

    //2/dw*(dw_de)
    W *= -1.0 / (dw * (dw + de));
    E *= -1.0 / (de * (dw + de));
    N *= -1.0 / (dn * (dn + ds));
    S *= -1.0 / (ds * (dn + ds));
    F = -1 / (hr[2] * (hr[2] + hr[2]));
    B = -1 / (hr[2] * (hr[2] + hr[2]));

    //TODO: problem when de,dw,dn,ds == 0
    //is NOT a regular BOUND PT
    C += 1 / de * 1 / dw;
    C += 1 / dn * 1 / ds;
    C += 1 / hr[2] * 1 / hr[2];
    scaleFactor = 0.5;


    //for regular gridpoints no problem with symmetry, just boundary
    //z direction is right
    //implement isLastInside(dir)
    //we have LOCAL x,y coordinates!

    /*
       if(dw != 0 && !wIsB)
       W = -1/dw * (dn+ds) * 2*hr[2];
       else
       W = 0;
       if(de != 0 && !eIsB)
       E = -1/de * (dn+ds) * 2*hr[2];
       else
       E = 0;
       if(dn != 0 && !nIsB)
       N = -1/dn * (dw+de) * 2*hr[2];
       else
       N = 0;
       if(ds != 0 && !sIsB)
       S = -1/ds * (dw+de) * 2*hr[2];
       else
       S = 0;
       F = -(dw+de)*(dn+ds)/hr[2];
       B = -(dw+de)*(dn+ds)/hr[2];
       */

    //if(dw != 0)
    //W = -2*hr[2]*(dn+ds)/dw;
    //else
    //W = 0;
    //if(de != 0)
    //E = -2*hr[2]*(dn+ds)/de;
    //else
    //E = 0;
    //if(dn != 0)
    //N = -2*hr[2]*(dw+de)/dn;
    //else
    //N = 0;
    //if(ds != 0)
    //S = -2*hr[2]*(dw+de)/ds;
    //else
    //S = 0;
    //F = -(dw+de)*(dn+ds)/hr[2];
    //B = -(dw+de)*(dn+ds)/hr[2];

    //// RHS scaleFactor for current 3D index
    //// Factor 0.5 results from the SW/quadratic extrapolation
    //scaleFactor = 0.5*(dw+de)*(dn+ds)*(2*hr[2]);

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
    ////C = 2*(dn+ds)*(dw+de)/hr[2];
    //C = 2/hr[2];
    //if(dw != 0 || de != 0)
    //C *= (dw+de);
    //if(dn != 0 || ds != 0)
    //C *= (dn+ds);
    //if(dw != 0 || de != 0)
    //C += (2*hr[2])*(dn+ds)*(dw+de)/m1;
    //if(dn != 0 || ds != 0)
    //C += (2*hr[2])*(dw+de)*(dn+ds)/m2;

    //handle Neumann case
    //if(z == 0 || z == nr[2]-1) {

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
    if(z == 0 || z == nr[2] - 1) {

        // case where we are on the NEUMAN BC in Z-direction
        // where we distinguish two cases
        if(z == 0)
            F = 0.0;
        else
            B = 0.0;

        //C += 2/((hr[2]*(nr[2]-1)/2.0) * hr[2]);
        //hr[2]*(nr2[2]-1)/2 = radius
        double d = hr[2] * (nr[2] - 1) / 2;
        C += 2 / (d * hr[2]);

        W /= 2.0;
        E /= 2.0;
        N /= 2.0;
        S /= 2.0;
        C /= 2.0;
        scaleFactor /= 2.0;

    }
}


#endif //#ifdef HAVE_SAAMG_SOLVER