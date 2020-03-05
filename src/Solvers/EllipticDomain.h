#ifndef ELLIPTICAL_DOMAIN_H
#define ELLIPTICAL_DOMAIN_H

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "IrregularDomain.h"
#include "Structure/BoundaryGeometry.h"

class EllipticDomain : public IrregularDomain {

public:

    EllipticDomain(Vector_t nr, Vector_t hr);
    EllipticDomain(double semimajor, double semiminor, Vector_t nr, Vector_t hr, std::string interpl);
    EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr, Vector_t hr, std::string interpl);
    EllipticDomain(BoundaryGeometry *bgeom, Vector_t nr, std::string interpl);

    ~EllipticDomain();

    /// returns discretization at (x,y,z)
    void getBoundaryStencil(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns discretization at 3D index
    void getBoundaryStencil(int idx, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    /// returns index of neighbours at (x,y,z)
    void getNeighbours(int x, int y, int z, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns index of neighbours at 3D index
    void getNeighbours(int idx, int &W, int &E, int &S, int &N, int &F, int &B);
    /// returns type of boundary condition
    std::string getType() {return "Elliptic";}
    /// queries if a given (x,y,z) coordinate lies inside the domain
    inline bool isInside(int x, int y, int z) {
        double xx = (x - (nr[0] - 1) / 2.0) * hr[0];
        double yy = (y - (nr[1] - 1) / 2.0) * hr[1];
        return ((xx * xx / (SemiMajor * SemiMajor) + yy * yy / (SemiMinor * SemiMinor) < 1) && z != 0 && z != nr[2] - 1);
    }

    int getNumXY(int /*z*/) { return nxy_m; }
    /// set semi-minor
    void setSemiMinor(double sm) {SemiMinor = sm;}
    /// set semi-major
    void setSemiMajor(double sm) {SemiMajor = sm;}

    /// calculates intersection 
    void compute(Vector_t hr);
    void compute(Vector_t hr, NDIndex<3> localId);

    double getXRangeMin() { return -SemiMajor; }
    double getXRangeMax() { return SemiMajor;  }
    double getYRangeMin() { return -SemiMinor; }
    double getYRangeMax() { return SemiMinor;  }
    double getZRangeMin() { return zMin_m; }
    double getZRangeMax() { return zMax_m; }


    //TODO: ?
    int getStartIdx() {return 0;}

    bool hasGeometryChanged() { return hasGeometryChanged_m; }

private:

    /// Map from a single coordinate (x or y) to a list of intersection values with
    /// boundary.
    typedef std::multimap<int, double> EllipticPointList;

    /// all intersection points with grid lines in X direction
    EllipticPointList IntersectXDir;

    /// all intersection points with grid lines in Y direction
    EllipticPointList IntersectYDir;

    /// mapping (x,y,z) -> idx
    std::map<int, int> IdxMap;

    /// mapping idx -> (x,y,z)
    std::map<int, int> CoordMap;

    /// semi-major of the ellipse
    double SemiMajor;
    /// semi-minor of the ellipse
    double SemiMinor;
    /// number of nodes in the xy plane (for this case: independent of the z coordinate)
    int nxy_m;
    /// interpolation type
    int interpolationMethod;
    /// flag indicating if geometry has changed for the current time-step
    bool hasGeometryChanged_m;

    /// conversion from (x,y) to index in xy plane
    inline int toCoordIdx(int x, int y) { return y * nr[0] + x; }
    /// conversion from (x,y,z) to index on the 3D grid
    inline int getIdx(int x, int y, int z) {
        if(isInside(x, y, z) && x >= 0 && y >= 0 && z >= 0)
            return IdxMap[toCoordIdx(x, y)] + (z - 1) * nxy_m;
        else
            return -1;
    }
    /// conversion from a 3D index to (x,y,z)
    inline void getCoord(int idx, int &x, int &y, int &z) {
        int ixy = idx % nxy_m;
        int xy = CoordMap[ixy];
        int inr = (int)nr[0];
        x = xy % inr;
        y = (xy - x) / nr[0];
        z = (idx - ixy) / nxy_m + 1;
    }

    /// different interpolation methods for boundary points
    void constantInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void linearInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);
    void quadraticInterpolation(int x, int y, int z, double &W, double &E, double &S, double &N, double &F, double &B, double &C, double &scaleFactor);

};

#endif //#ifdef ELLIPTICAL_DOMAIN_H

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
