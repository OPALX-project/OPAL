// ------------------------------------------------------------------------
// $Version: 1.2.1 $
// ------------------------------------------------------------------------
// Copyright & License: See Copyright.readme in src directory
// ------------------------------------------------------------------------
// Class ArbitraryDomain
//   Interface to iterative solver and boundary geometry
//   for space charge calculation
//
// ------------------------------------------------------------------------
// $Author: kaman $
// $Date: 2014 $
// ------------------------------------------------------------------------
//#define DEBUG_INTERSECT_RAY_BOUNDARY

#ifdef HAVE_SAAMG_SOLVER
#include "Solvers/ArbitraryDomain.h"
#include "Structure/BoundaryGeometry.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <assert.h>

#include <math.h>

ArbitraryDomain::ArbitraryDomain( BoundaryGeometry * bgeom,
                                  Vector_t nr,
                                  Vector_t hr,
                                  std::string interpl){
    bgeom_m  = bgeom;
    minCoords_m = bgeom->getmincoords();
    maxCoords_m = bgeom->getmaxcoords();
    geomCentroid_m = (minCoords_m + maxCoords_m)/2.0;

    // TODO: THis needs to be made into OPTION of the geometry.
    // A user defined point that is INSIDE with 100% certainty. -DW
    globalInsideP0_m = Vector_t(0.0, 0.0, -0.13);

    setNr(nr);
    for(int i=0; i<3; i++)
        Geo_hr_m[i] = (maxCoords_m[i] - minCoords_m[i])/nr[i];
    setHr(Geo_hr_m);

    startId = 0;

    if (interpl == "CONSTANT")
        interpolationMethod = CONSTANT;
    else if(interpl == "LINEAR")
        interpolationMethod = LINEAR;
    else if(interpl == "QUADRATIC")
        interpolationMethod = QUADRATIC;
}

ArbitraryDomain::~ArbitraryDomain() {
    //nothing so far
}

void ArbitraryDomain::compute(Vector_t hr){

    setHr(hr);

    globalMeanR_m = getGlobalMeanR();

    globalToLocalQuaternion_m = getGlobalToLocalQuaternion();
    localToGlobalQuaternion_m[0] = globalToLocalQuaternion_m[0];
    for (int i=1; i<4; i++)
        localToGlobalQuaternion_m[i] = -globalToLocalQuaternion_m[i];

    hasGeometryChanged_m = true;

    IntersectLoX.clear();
    IntersectHiX.clear();
    IntersectLoY.clear();
    IntersectHiY.clear();
    IntersectLoZ.clear();
    IntersectHiZ.clear();

    //calculate intersection
    Vector_t P, saveP, dir, I;
    //Reference Point
    Vector_t P0 = geomCentroid_m;

    for (int idz = 0; idz < nr[2] ; idz++) {
        saveP[2] = (idz - (nr[2]-1)/2.0)*hr[2];
        for (int idy = 0; idy < nr[1] ; idy++) {
            saveP[1] = (idy - (nr[1]-1)/2.0)*hr[1];
            for (int idx = 0; idx <nr[0]; idx++) {
                saveP[0] = (idx - (nr[0]-1)/2.0)*hr[0];
                P = saveP;

                rotateWithQuaternion(P, localToGlobalQuaternion_m);
                P += geomCentroid_m;

                if (bgeom_m->fastIsInside(P0, P) % 2 == 0) {
                    P0 = P;

                    std::tuple<int, int, int> pos(idx, idy, idz);

                    rotateZAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    rotateYAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=-1" << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    rotateXAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)){
                        I -= geomCentroid_m;
                        rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }
                } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                    *gmsg << "OUTSIDE" << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                }
            }
        }
    }
    IdxMap.clear();
    CoordMap.clear();

    int id=0;
    int idx, idy, idz;
    for (idz = 0; idz < nr[2]; idz++) {
        for (idy = 0; idy < nr[1]; idy++) {
            for (idx = 0; idx < nr[0]; idx++) {
                if (isInside(idx, idy, idz)) {
                    IdxMap[toCoordIdx(idx, idy, idz)] = id;
                    CoordMap[id] = toCoordIdx(idx, idy, idz);
                    id++;
                }
            }
        }
    }
}

void ArbitraryDomain::compute(Vector_t hr, NDIndex<3> localId){

    INFOMSG(level2 << "* Starting the Boundary Intersection Tests..." << endl);

    setHr(hr);

    globalMeanR_m = getGlobalMeanR();

    globalToLocalQuaternion_m = getGlobalToLocalQuaternion();
    localToGlobalQuaternion_m[0] = globalToLocalQuaternion_m[0];
    for (int i=1; i<4; i++)
        localToGlobalQuaternion_m[i] = -globalToLocalQuaternion_m[i];

    int zGhostOffsetLeft  = (localId[2].first()== 0) ? 0 : 1;
    int zGhostOffsetRight = (localId[2].last() == nr[2] - 1) ? 0 : 1;
    int yGhostOffsetLeft  = (localId[1].first()== 0) ? 0 : 1;
    int yGhostOffsetRight = (localId[1].last() == nr[1] - 1) ? 0 : 1;
    int xGhostOffsetLeft  = (localId[0].first()== 0) ? 0 : 1;
    int xGhostOffsetRight = (localId[0].last() == nr[0] - 1) ? 0 : 1;

    hasGeometryChanged_m = true;

    IntersectLoX.clear();
    IntersectHiX.clear();
    IntersectLoY.clear();
    IntersectHiY.clear();
    IntersectLoZ.clear();
    IntersectHiZ.clear();

    // Calculate intersection
    Vector_t P, dir, I;
    // Vector_t saveP, saveP_old;
    Vector_t P0 = globalInsideP0_m;

    // We cannot assume that the geometry is symmetric about the xy, xz, and yz planes!
    // In my spiral inflector simulation, this is not the case for z direction for
    // example (-0.13 to +0.025). -DW
    for (int idz = localId[2].first()-zGhostOffsetLeft; idz <= localId[2].last()+zGhostOffsetRight; idz++) {

        //saveP_old[2] = (idz - (nr[2]-1)/2.0)*hr[2];
        P[2] = minCoords_m[2] + (idz + 0.5) * hr[2];

        for (int idy = localId[1].first()-yGhostOffsetLeft; idy <= localId[1].last()+yGhostOffsetRight; idy++) {

            //saveP_old[1] = (idy - (nr[1]-1)/2.0)*hr[1];
            P[1] = minCoords_m[1] + (idy + 0.5) * hr[1];

            for (int idx = localId[0].first()-xGhostOffsetLeft; idx <= localId[0].last()+xGhostOffsetRight; idx++) {

                //saveP_old[0] = (idx - (nr[0]-1)/2.0)*hr[0];
                P[0] = minCoords_m[0] + (idx + 0.5) * hr[0];

                // *gmsg << "Now working on point " << saveP << " (original was " << saveP_old << ")" << endl;

                //P = saveP;

                //rotateWithQuaternion(P, localToGlobalQuaternion_m);
                //P += geomCentroid_m; //sorry, this doesn't make sense. -DW
                //P += globalMeanR_m;

                if (bgeom_m->fastIsInside(P0, P) % 2 == 0) {

                    // Fill the map with true or false values for very fast isInside tests
                    // during the rest of the fieldsolve.
                    IsInsideMap[toCoordIdx(idx, idy, idz)] = true;

                    // Replace the old reference point with the new point (which we know is
                    // inside because we just tested for it. This makes the algorithm faster
                    // because fastIsInside() creates a number of segments that depends on the
                    // distance between P and P0. Using the previous P as the new P0
                    // assures the smallest possible distance in most cases. -DW
                    P0 = P;

                    std::tuple<int, int, int> pos(idx, idy, idz);

                    //rotateZAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    dir = Vector_t(0, 0, 1);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoZ.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[2]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "zdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    //rotateYAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    dir = Vector_t(0, 1, 0);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)) {
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoY.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[1]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "ydir=-1" << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    //rotateXAxisWithQuaternion(dir, localToGlobalQuaternion_m);
                    dir = Vector_t(1, 0, 0);

                    if (bgeom_m->intersectRayBoundary(P, dir, I)) {
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectHiX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=+1 " << dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }

                    if (bgeom_m->intersectRayBoundary(P, -dir, I)){
                        //I -= geomCentroid_m;
                        //I -= globalMeanR_m;
                        //rotateWithQuaternion(I, globalToLocalQuaternion_m);
                        IntersectLoX.insert(std::pair< std::tuple<int, int, int>, double >(pos, I[0]));
                    } else {
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                        *gmsg << "xdir=-1 " << -dir << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                    }
                } else {
                    IsInsideMap[toCoordIdx(idx, idy, idz)] = false;
#ifdef DEBUG_INTERSECT_RAY_BOUNDARY
                    *gmsg << "OUTSIDE" << " x,y,z= " << idx << "," << idy << "," << idz << " P=" << P <<" I=" << I << endl;
#endif
                }
            }
        }
    }

    INFOMSG(level2 << "* Finding number of ghost nodes to the left..." << endl);

    //number of ghost nodes to the left
    int numGhostNodesLeft = 0;
    if(localId[2].first() != 0) {
        for(int idx = 0; idx < nr[0]; idx++) {
            for(int idy = 0; idy < nr[1]; idy++) {
                if(isInside(idx, idy, localId[2].first() - zGhostOffsetLeft))
                    numGhostNodesLeft++;
            }
        }
    }

    INFOMSG(level2 << "* Finding number of xy points in each plane along z..." << endl);

    //xy points in z plane
    int numtotal = 0;
    numXY.clear();
    for(int idz = localId[2].first(); idz <= localId[2].last(); idz++) {
        int numxy = 0;
        for(int idx = 0; idx < nr[0]; idx++) {
            for(int idy = 0; idy < nr[1]; idy++) {
                if(isInside(idx, idy, idz))
                    numxy++;
            }
        }
        numXY[idz-localId[2].first()] = numxy;
        numtotal += numxy;
    }

    int startIdx = 0;
    MPI_Scan(&numtotal, &startIdx, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    startIdx -= numtotal;

    // Build up index and coord map
    IdxMap.clear();
    CoordMap.clear();
    int index = startIdx - numGhostNodesLeft;

    INFOMSG(level2 << "* Building up index and coordinate map..." << endl);

    for(int idz = localId[2].first() - zGhostOffsetLeft; idz <= localId[2].last() + zGhostOffsetRight; idz++) {
        for(int idy = 0; idy < nr[1]; idy++) {
            for(int idx = 0; idx < nr[0]; idx++) {
                if(isInside(idx, idy, idz)) {
                    IdxMap[toCoordIdx(idx, idy, idz)] = index;
                    CoordMap[index] = toCoordIdx(idx, idy, idz);
                    index++;
                }
            }
        }
    }

    INFOMSG(level2 << "* Done." << endl);
}

// Conversion from (x,y,z) to index in xyz plane
inline int ArbitraryDomain::toCoordIdx(int idx, int idy, int idz) {
    return (idz * nr[1] + idy) * nr[0]  + idx;
}

// Conversion from (x,y,z) to index on the 3D grid
int ArbitraryDomain::getIdx(int idx, int idy, int idz) {

    if(isInside(idx, idy, idz) && idx >= 0 && idy >= 0 && idz >= 0)
        return IdxMap[toCoordIdx(idx, idy, idz)];
    else
        return -1;
}

// Conversion from a 3D index to (x,y,z)
inline void ArbitraryDomain::getCoord(int idxyz, int &idx, int &idy, int &idz) {
    int id = CoordMap[idxyz];
    idx = id % (int)nr[0];
    id /= nr[0];
    idy = id % (int)nr[1];
    id /= nr[1];
    idz = id;
}

inline bool ArbitraryDomain::isInside(int idx, int idy, int idz) {

    return IsInsideMap[toCoordIdx(idx, idy, idz)];
}

/*
  inline bool ArbitraryDomain::isInside(int idx, int idy, int idz) {
  Vector_t P;

  P[0] = minCoords_m[0] + (idx + 0.5) * hr[0];
  P[1] = minCoords_m[1] + (idy + 0.5) * hr[1];
  P[2] = minCoords_m[2] + (idz + 0.5) * hr[2];

  return (bgeom_m->fastIsInside(globalInsideP0_m, P) % 2 == 0);
  }
*/

/*
  inline bool ArbitraryDomain::isInside(int idx, int idy, int idz) {
  Vector_t P;

  P[0] = (idx - (nr[0]-1)/2.0) * hr[0];
  P[1] = (idy - (nr[1]-1)/2.0) * hr[1];
  P[2] = (idz - (nr[2]-1)/2.0) * hr[2];

  bool ret = false;
  int  countH, countL;
  std::multimap < std::tuple<int, int, int>, double >::iterator itrH, itrL;
  std::tuple<int, int, int> coordxyz(idx, idy, idz);

  //check if z is inside with x,y coords
  itrH = IntersectHiZ.find(coordxyz);
  itrL = IntersectLoZ.find(coordxyz);

  countH = IntersectHiZ.count(coordxyz);
  countL = IntersectLoZ.count(coordxyz);
  if(countH == 1 && countL == 1)
  ret = (P[2] <= itrH->second) && (P[2] >= itrL->second);

  //check if y is inside with x,z coords
  itrH = IntersectHiY.find(coordxyz);
  itrL = IntersectLoY.find(coordxyz);

  countH = IntersectHiY.count(coordxyz);
  countL = IntersectLoY.count(coordxyz);
  if(countH == 1 && countL == 1)
  ret = ret && (P[1] <= itrH->second) && (P[1] >= itrL->second);

  //check if x is inside with y,z coords
  itrH = IntersectHiX.find(coordxyz);
  itrL = IntersectLoX.find(coordxyz);

  countH = IntersectHiX.count(coordxyz);
  countL = IntersectLoX.count(coordxyz);
  if(countH == 1 && countL == 1)
  ret = ret && (P[0] <= itrH->second) && (P[0] >= itrL->second);

  return ret;
  }
*/

int ArbitraryDomain::getNumXY(int z) {

    return numXY[z];
}

void ArbitraryDomain::getBoundaryStencil(int idxyz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {
    int idx = 0, idy = 0, idz = 0;

    getCoord(idxyz, idx, idy, idz);
    getBoundaryStencil(idx, idy, idz, W, E, S, N, F, B, C, scaleFactor);
}

void ArbitraryDomain::getBoundaryStencil(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor) {

    scaleFactor = 1.0;
    // determine which interpolation method we use for points near the boundary
    switch(interpolationMethod){
    case CONSTANT:
        constantInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        break;
    case LINEAR:
        linearInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        break;
    case QUADRATIC:
        //  QuadraticInterpolation(idx,idy,idz,W,E,S,N,F,B,C,scaleFactor);
        break;
    }

    // stencil center value has to be positive!
    assert(C > 0);
}

void ArbitraryDomain::constantInterpolation(int idx, int idy, int idz, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor) {

    W = -1/(hr[0]*hr[0]);
    E = -1/(hr[0]*hr[0]);
    N = -1/(hr[1]*hr[1]);
    S = -1/(hr[1]*hr[1]);
    F = -1/(hr[2]*hr[2]);
    B = -1/(hr[2]*hr[2]);
    C = 2/(hr[0]*hr[0]) + 2/(hr[1]*hr[1]) + 2/(hr[2]*hr[2]);

    if(!isInside(idx-1,idy,idz))
        W = 0.0;
    if(!isInside(idx+1,idy,idz))
        E = 0.0;

    if(!isInside(idx,idy+1,idz))
        N = 0.0;
    if(!isInside(idx,idy-1,idz))
        S = 0.0;

    if(!isInside(idx,idy,idz-1))
        F = 0.0;
    if(!isInside(idx,idy,idz+1))
        B = 0.0;
}

void ArbitraryDomain::linearInterpolation(int idx, int idy, int idz, double& W, double& E, double& S, double& N, double& F, double& B, double& C, double &scaleFactor)
{
    scaleFactor = 1;

    double cx = (idx - (nr[0]-1)/2.0)*hr[0];
    double cy = (idy - (nr[1]-1)/2.0)*hr[1];
    double cz = (idz - (nr[2]-1)/2.0)*hr[2];

    double dx_w=hr[0];
    double dx_e=hr[0];
    double dy_n=hr[1];
    double dy_s=hr[1];
    double dz_f=hr[2];
    double dz_b=hr[2];
    C = 0.0;

    std::tuple<int, int, int> coordxyz(idx, idy, idz);

    if (idx == nr[0]-1)
        dx_e = fabs(IntersectHiX.find(coordxyz)->second - cx);
    if (idx == 0)
        dx_w = fabs(IntersectLoX.find(coordxyz)->second - cx);
    if (idy == nr[1]-1)
        dy_n = fabs(IntersectHiY.find(coordxyz)->second - cy);
    if (idy == 0)
        dy_s = fabs(IntersectLoY.find(coordxyz)->second - cy);
    if (idz == nr[2]-1)
        dz_b = fabs(IntersectHiZ.find(coordxyz)->second - cz);
    if (idz == 0)
        dz_f = fabs(IntersectLoZ.find(coordxyz)->second - cz);

    if(dx_w != 0)
        W = -(dz_f + dz_b) * (dy_n + dy_s) / dx_w;
    else
        W = 0;
    if(dx_e != 0)
        E = -(dz_f + dz_b) * (dy_n + dy_s) / dx_e;
    else
        E = 0;
    if(dy_n != 0)
        N = -(dz_f + dz_b) * (dx_w + dx_e) / dy_n;
    else
        N = 0;
    if(dy_s != 0)
        S = -(dz_f + dz_b) * (dx_w + dx_e) / dy_s;
    else
        S = 0;
    if(dz_f != 0)
        F = -(dx_w + dx_e) * (dy_n + dy_s) / dz_f;
    else
        F = 0;
    if(dz_b != 0)
        B = -(dx_w + dx_e) * (dy_n + dy_s) / dz_b;
    else
        B = 0;

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

    C = 2 / hr[2];
    if(dx_w != 0 || dx_e != 0)
        C *= (dx_w + dx_e);
    if(dy_n != 0 || dy_s != 0)
        C *= (dy_n + dy_s);
    if(dx_w != 0 || dx_e != 0)
        C += (dz_f + dz_b) * (dy_n + dy_s) * (dx_w + dx_e) / m1;
    if(dy_n != 0 || dy_s != 0)
        C += (dz_f + dz_b) * (dx_w + dx_e) * (dy_n + dy_s) / m2;
}

void ArbitraryDomain::getNeighbours(int id, int &W, int &E, int &S, int &N, int &F, int &B) {

    int idx = 0, idy = 0, idz = 0;

    getCoord(id, idx, idy, idz);
    getNeighbours(idx, idy, idz, W, E, S, N, F, B);
}

void ArbitraryDomain::getNeighbours(int idx, int idy, int idz, int &W, int &E, int &S, int &N, int &F, int &B) {

    W = getIdx(idx - 1, idy, idz);
    E = getIdx(idx + 1, idy, idz);
    N = getIdx(idx, idy + 1, idz);
    S = getIdx(idx, idy - 1, idz);
    F = getIdx(idx, idy, idz - 1);
    B = getIdx(idx, idy, idz + 1);

    if(!isInside(idx+1,idy,idz))
        E = -1;

    if(!isInside(idx-1,idy,idz))
        W = -1;

    if(!isInside(idx,idy+1,idz))
        N = -1;

    if(!isInside(idx,idy-1,idz))
        S = -1;

    if(!isInside(idx,idy,idz-1))
        F = -1;

    if(!isInside(idx,idy,idz+1))
        B = -1;

}


inline void ArbitraryDomain::crossProduct(double A[], double B[], double C[]) {
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
}

inline void ArbitraryDomain::rotateWithQuaternion(Vector_t & v, Quaternion_t const quaternion) {
    // rotates a Vector_t (3 elements) using a quaternion.
    // Flip direction of rotation by quaternionVectorcomponent *= -1

    Vector_t const quaternionVectorComponent = Vector_t(quaternion(1), quaternion(2), quaternion(3));
    double const quaternionScalarComponent = quaternion(0);

    v = 2.0 * dot(quaternionVectorComponent, v) * quaternionVectorComponent
        + (quaternionScalarComponent * quaternionScalarComponent
           -  dot(quaternionVectorComponent, quaternionVectorComponent)) * v
        + 2.0 * quaternionScalarComponent * cross(quaternionVectorComponent, v);
}

inline void ArbitraryDomain::rotateXAxisWithQuaternion(Vector_t & v, Quaternion_t const quaternion) {
    // rotates the positive xaxis using a quaternion.

    v(0) = (quaternion(0) * quaternion(0)
            + quaternion(1) * quaternion(1)
            - quaternion(2) * quaternion(2)
            - quaternion(3) * quaternion(3));

    v(1) = 2.0 * (quaternion(1) * quaternion(2) + quaternion(0) * quaternion(3));
    v(2) = 2.0 * (quaternion(1) * quaternion(3) - quaternion(0) * quaternion(2));
}

inline void ArbitraryDomain::rotateYAxisWithQuaternion(Vector_t & v, Quaternion_t const quaternion) {
    // rotates the positive yaxis using a quaternion.

    v(0) = 2.0 * (quaternion(1) * quaternion(2) - quaternion(0) * quaternion(3));

    v(1) = (quaternion(0) * quaternion(0)
            - quaternion(1) * quaternion(1)
            + quaternion(2) * quaternion(2)
            - quaternion(3) * quaternion(3));

    v(2) = 2.0 * (quaternion(2) * quaternion(3) + quaternion(0) * quaternion(1));
}

inline void ArbitraryDomain::rotateZAxisWithQuaternion(Vector_t & v, Quaternion_t const quaternion) {
    // rotates the positive zaxis using a quaternion.
    v(0) = 2.0 * (quaternion(1) * quaternion(3) + quaternion(0) * quaternion(2));
    v(1) = 2.0 * (quaternion(2) * quaternion(3) - quaternion(0) * quaternion(1));

    v(2) = (quaternion(0) * quaternion(0)
            - quaternion(1) * quaternion(1)
            - quaternion(2) * quaternion(2)
            + quaternion(3) * quaternion(3));

}
#endif //#ifdef HAVE_SAAMG_SOLVER