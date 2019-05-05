#ifndef BOXCORNER_DOMAIN_H
#define BOXCORNER_DOMAIN_H
#ifdef HAVE_SAAMG_SOLVER

//#include <vector>
#include <map>
//#include <multimap>
#include <string>
#include <cmath>
#include <iostream>  // Neeeded for stream I/O
#include <fstream>   // Needed for file I/O
//#include <iomanip>   // Needed for I/O manipulators

//#include "Ippl.h"
#include "IrregularDomain.h"


/*

    A_m and B_m are the half apperture of the box


                                     / (A_m,B_m)
                                    /
                                   /
                                  /
    L1_m                         /
------------      --------------+ (-A_m,B_m)
           | L2_m |             |
        C_m|      |             |
           |------|             |      /
         .....                  |     /
(0,0)---.......-----------------+    /
         .....                  |   /
   z                            |  /
   |                            | /
--------------------------------+/ (-A_m,-B_m)

            Length_m

Test in which of the 3 parts of the geometry we are in.

    if((z < L1_m) || (z > (L1_m + L2_m)))
        b = B_m;
    else
        b = B_m-C_m;

*/

class BoxCornerDomain : public IrregularDomain {

public:

    BoxCornerDomain(Vector_t nr, Vector_t hr);
    BoxCornerDomain(double A, double B, double C, double Length, double L1, double L2, Vector_t nr, Vector_t hr,
                    std::string interpl);
    ~BoxCornerDomain();

    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B,
                            double &C, double &scaleFactor);

    /// returns discretization at 3D index
    void getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C,
                            double &scaleFactor);

    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B);

    /// returns index of neighbours at 3D index
    void getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B);

    /// returns type of boundary condition
    std::string getType() {return "BoxCorner";}


    /// we do not need this
    int getNumXY(int z) { return -1;}


    /// as a function of z, determine the hight (B) of the geometry
    inline double getB(double z) {
      if((z < L1_m) || (z > (L1_m + L2_m)))
            return B_m;
        else
            return B_m - C_m;
    }

    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int z) {
        const double xx = (x - (nr[0] - 1) / 2.0) * hr[0];
        const double yy = (y - (nr[1] - 1) / 2.0) * hr[1];
        const double b = getB(z * hr[2]);
        return (xx < A_m && yy < b && z != 0 && z != nr[2] - 1);
    }

    /// set semi-minor
    //void setSemiMinor(double sm) {SemiMinor = sm;}
    /// set semi-major
    //void setSemiMajor(double sm) {SemiMajor = sm;}

    void compute(Vector_t hr);
    void compute(Vector_t hr, NDIndex<3> localId);

    double getXRangeMin() { return -A_m; }
    double getXRangeMax() { return  A_m; }
    double getYRangeMin() { return  -B_m;} // actBMin_m; }
    double getYRangeMax() { return  B_m; } // actBMax_m; }
    double getZRangeMin() { return  L1_m;}
    double getZRangeMax() { return  L1_m+L2_m; }


    //TODO: ?
    int getStartIdx() {return 0;}

    bool hasGeometryChanged() { return hasGeometryChanged_m; }

private:

    //XXX: since the Y coorindate is dependent on the Z value we need (int,
    //int) -> intersection. To simplify things (for now) we use the same
    //structure for X...
    /// Map from a ([(x or y], z) to a list of intersection values with
    /// boundary.
    typedef std::multimap< std::pair<int, int>, double > BoxCornerPointList;

    /// all intersection points with grid lines in X direction
    BoxCornerPointList IntersectXDir;

    /// all intersection points with grid lines in Y direction
    BoxCornerPointList IntersectYDir;

    /// mapping (x,y,z) -> idx
    std::map<int, int> IdxMap;

    /// mapping idx -> (x,y,z)
    std::map<int, int> CoordMap;

    /// depth of the box
    double A_m;

    /// the maximal hight of the box
    double B_m;


    /// because the geometry can change in the y direction
    double actBMin_m;

    double actBMax_m;

    /// hight of the corner
    double C_m;

    /// lenght of the structure
    double Length_m;

    /// lenght of the first part of the structure
    double L1_m;

    /// lenght of the corner
    double L2_m;

    /// semi-major of the ellipse
    // :FIXME: unused
    //double SemiMajor;
    /// semi-minor of the ellipse
    // :FIXME: unused
    //double SemiMinor;

    /// interpolation type
    int interpolationMethod;
    /// flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    /// for debug reasons
   std::ofstream os_m;





    inline double getXIntersection(double cx, int z) {
        if(cx < 0)
            return -A_m;
        else
            return  A_m;
    }

    inline double getYIntersection(double cy, int z) {
        if(cy < 0)
            return -B_m;
        else
            return getB(z * hr[2]);
    }

    /// conversion from (x,y,z) to index in xyz plane
    inline int toCoordIdx(int x, int y, int z) {
        return (z * nr[1] + y) * nr[0] + x;
    }

    /// conversion from (x,y,z) to index on the 3D grid
    /*inline*/
    inline int getIdx(int x, int y, int z) {
        if(isInside(x, y, z) && x >= 0 && y >= 0 && z >= 0)
            return IdxMap[toCoordIdx(x, y, z)];
        else
            return -1;
    }

    /// conversion from a 3D index to (x,y,z)
    inline void getCoord(int idx, int &x, int &y, int &z) {

        int idxx = CoordMap[idx];

        x = idxx % (int)nr[0];
        idxx /= nr[0];
        y = idxx % (int)nr[1];
        idxx /= nr[1];
        z = idxx;

    }

    /// different interpolation methods for boundary points
    void constantInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void linearInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void quadraticInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);

};

#endif //#ifdef HAVE_SAAMG_SOLVER
#endif //#ifdef BOXCORNER_DOMAIN_H