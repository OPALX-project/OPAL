#ifndef CLASSIC_FIELDMAP3DMAGNETOSTATIC_HH
#define CLASSIC_FIELDMAP3DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

class FM3DMagnetoStatic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    virtual bool isInside(const Vector_t &r) const;
private:
    FM3DMagnetoStatic(std::string aFilename);
    ~FM3DMagnetoStatic();

    virtual void readMap();
    virtual void freeMap();

    struct IndexTriplet {
        unsigned int i;
        unsigned int j;
        unsigned int k;
        Vector_t weight;
        IndexTriplet():
            i(0),
            j(0),
            k(0),
            weight(0.0)
        { }
    };

    IndexTriplet getIndex(const Vector_t &X) const;
    unsigned long getIndex(unsigned int i, unsigned int j, unsigned int k) const;
    Vector_t interpolateTrilinearly(const Vector_t &X) const;
    double getWeightedData(double *data, const IndexTriplet &idx, unsigned short corner) const;

    enum { LX = 0,  // low X
           LY = 0,  // low Y
           LZ = 0,  // low Z
           HX = 4,  // high X
           HY = 2,  // high Y
           HZ = 1}; // high Z

    double *FieldstrengthBz_m;    /**< 3D array with Bz */
    double *FieldstrengthBx_m;    /**< 3D array with Bx */
    double *FieldstrengthBy_m;    /**< 3D array with By */

    double xbegin_m;
    double xend_m;

    double ybegin_m;
    double yend_m;

    double zbegin_m;
    double zend_m;
    double length_m;

    double hx_m;                   /**< length between points in grid, x-direction */
    double hy_m;                   /**< length between points in grid, y-direction */
    double hz_m;                   /**< length between points in grid, z-direction */
    unsigned int num_gridpx_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    unsigned int num_gridpy_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    unsigned int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    friend class Fieldmap;
};

inline
bool FM3DMagnetoStatic::isInside(const Vector_t &r) const
{
    return ((r(0) >= xbegin_m && r(0) < xend_m) &&
            (r(1) >= ybegin_m && r(1) < yend_m) &&
            (r(2) >= 0.0 && r(2) < length_m));
            // (r(2) >= zbegin_m && r(2) < zend_m));
}

inline
double FM3DMagnetoStatic::getFrequency() const {
    return 0.0;
}

inline
void FM3DMagnetoStatic::setFrequency(double /*freq*/) { }

inline
unsigned long FM3DMagnetoStatic::getIndex(unsigned int i, unsigned int j, unsigned int k) const
{
    unsigned long result = i + j * num_gridpx_m;
    result = k + result * num_gridpz_m;
    PAssert_LT(result, num_gridpx_m * num_gridpy_m * num_gridpz_m);
    return result;
}

inline
FM3DMagnetoStatic::IndexTriplet FM3DMagnetoStatic::getIndex(const Vector_t &X) const {
    IndexTriplet idx;
    idx.i = std::floor((X(0) - xbegin_m) / hx_m);
    idx.j = std::floor((X(1) - ybegin_m) / hy_m);
    idx.k = std::floor(X(2) / hz_m);
    // idx.k = std::floor((X(2) - zbegin_m) / hz_m);

    idx.weight(0) = (X(0) - xbegin_m) / hx_m - idx.i;
    idx.weight(1) = (X(1) - ybegin_m) / hy_m - idx.j;
    idx.weight(2) = (X(2) - zbegin_m) / hz_m - idx.k;

    return idx;
}

#endif
