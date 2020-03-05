//
// Class ArbitraryDomain
//   Interface to iterative solver and boundary geometry
//   for space charge calculation
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
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

    ArbitraryDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr, std::string interpl);
    ArbitraryDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr, Vector_t globalMeanR, Quaternion_t globalToLocalQuaternion, std::string interpl);

    ~ArbitraryDomain();

    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns discretization at 3D index
    void getBoundaryStencil(int idxyz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int idx, int idy, int idz, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns index of neighbours at 3D index
    void getNeighbours(int idxyz, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns type of boundary condition
    std::string getType() {return "Geometric";}
    /// queries if a given (x,y,z) coordinate lies inside the domain
    bool isInside(int idx, int idy, int idz);
    /// returns number of nodes in xy plane
    int getNumXY(int idz);
    // calculates intersection
    void compute(Vector_t hr);
    // calculates intersection with rotated and shifted geometry
    void compute(Vector_t hr, NDIndex<3> localId);

    int getStartId() {return startId;}

    double getXRangeMin(){ return minCoords_m(0); }
    double getYRangeMin(){ return minCoords_m(1); }
    double getZRangeMin(){ return minCoords_m(2); }

    double getXRangeMax(){ return maxCoords_m(0); }
    double getYRangeMax(){ return maxCoords_m(1); }
    double getZRangeMax(){ return maxCoords_m(2); }

    void setXRangeMin(double xmin){ minCoords_m(0) = xmin; }
    void setYRangeMin(double ymin){ minCoords_m(1) = ymin; }
    void setZRangeMin(double zmin){ minCoords_m(2) = zmin; }

    void setXRangeMax(double xmax){ maxCoords_m(0) = xmax; }
    void setYRangeMax(double ymax){ maxCoords_m(1) = ymax; }
    void setZRangeMax(double zmax){ maxCoords_m(2) = zmax; }


    bool hasGeometryChanged() { return hasGeometryChanged_m; }

private:
    BoundaryGeometry *bgeom_m;

    /// PointList maps from an (x,z) resp. (y,z) pair to double values (=intersections with boundary)
    typedef std::multimap< std::tuple<int, int, int>, double > PointList;

    /// all intersection points with gridlines in X direction
    PointList IntersectHiX, IntersectLoX;

    /// all intersection points with gridlines in Y direction
    PointList IntersectHiY, IntersectLoY;

    /// all intersection points with gridlines in Z direction
    PointList IntersectHiZ, IntersectLoZ;

    // meanR to shift from global to local frame
    Vector_t globalMeanR_m;
    //    Quaternion_t globalToLocalQuaternion_m;  because defined in parent class
    Quaternion_t localToGlobalQuaternion_m;

    int startId;

    // Here we store the number of nodes in a xy layer for a given z coordinate
    std::map<int, int> numXY;
    std::map<int, int> numYZ;
    std::map<int, int> numXZ;

    // Number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m[1000];
    // mapping (x,y,z) -> idxyz
    std::map<int, int> IdxMap;
    // Mapping idxyz -> (x,y,z)
    std::map<int, int> CoordMap;

    // Mapping all cells that are inside the geometry
    std::map<int, bool> IsInsideMap;

    // Interpolation type
    int interpolationMethod;

    // Flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    Vector_t Geo_nr_m;
    Vector_t Geo_hr_m;
    Vector_t geomCentroid_m;
    Vector_t minCoords_m;
    Vector_t maxCoords_m;
    Vector_t globalInsideP0_m;

    // Conversion from (x,y,z) to index in xyz plane
    inline int toCoordIdx(int idx, int idy, int idz);
    // Conversion from (x,y,z) to index on the 3D grid
    int getIdx(int idx, int idy, int idz);
    // Conversion from a 3D index to (x,y,z)
    inline void getCoord(int idxyz, int &x, int &y, int &z);

    inline void crossProduct(double A[], double B[], double C[]);
    inline double dotProduct(double v1[], double v2[]) { return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]); }

    // Different interpolation methods for boundary points
    void constantInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void linearInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void quadraticInterpolation(int idx, int idy, int idz, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);

    // Rotate positive axes with quaternion -DW
    inline void rotateWithQuaternion(Vector_t &v, Quaternion_t const quaternion);

    inline void rotateXAxisWithQuaternion(Vector_t &v, Quaternion_t const quaternion);
    inline void rotateYAxisWithQuaternion(Vector_t &v, Quaternion_t const quaternion);
    inline void rotateZAxisWithQuaternion(Vector_t &v, Quaternion_t const quaternion);
};

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
