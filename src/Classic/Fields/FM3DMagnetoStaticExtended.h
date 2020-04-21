#ifndef CLASSIC_FM2DMAGNETOSTATICEXTENDED_HH
#define CLASSIC_FM2DMAGNETOSTATICEXTENDED_HH

#include "Fields/Fieldmap.h"

class FM3DMagnetoStaticExtended: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    virtual bool isInside(const Vector_t &r) const;
private:
    FM3DMagnetoStaticExtended(std::string aFilename);
    ~FM3DMagnetoStaticExtended();

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

    enum { LX = 0,
           LY = 0,
           LZ = 0,
           HX = 4,
           HY = 2,
           HZ = 1};

    void integrateBx(unsigned j);
    void integrateBz(unsigned j);
    void integrateBy(unsigned j);

    void smoothData(double * data, unsigned j);

    void saveField(const std::string &fname, unsigned int j) const;

    double *FieldstrengthBx_m;    /**< 3D array with Bx, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
    double *FieldstrengthBy_m;    /**< 3D array with By, read in like Bx*/
    double *FieldstrengthBz_m;    /**< 3D array with Bz, read in like Bx*/

    double xbegin_m;
    double xend_m;
    double ybegin_m;
    double yend_m;
    double zbegin_m;
    double zend_m;
    double length_m;
    double hx_m;                   /**< length between points in grid, x-direction, m*/
    double hy_m;                   /**< length between points in grid, y-direction, m*/
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    unsigned int num_gridpx_m;     /**< Read in number of points after 0(not counted here) in grid, x-direction*/
    unsigned int num_gridpy_m;     /**< Read in number of points after 0(not counted here) in grid, y-direction*/
    unsigned int num_gridpz_m;     /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    friend class Fieldmap;
};

inline
bool FM3DMagnetoStaticExtended::isInside(const Vector_t &r) const
{
    return r(2) >= 0.0 && r(2) < length_m && r(0) >= xbegin_m && r(0) < xend_m && std::abs(r(1)) < yend_m;
}

inline
unsigned long FM3DMagnetoStaticExtended::getIndex(unsigned int i, unsigned int j, unsigned int k) const
{
    unsigned long result = i + j * num_gridpx_m;
    result = k + result * num_gridpz_m;
    PAssert_LT(result, num_gridpx_m * num_gridpy_m * num_gridpz_m);
    return result;
}

inline
FM3DMagnetoStaticExtended::IndexTriplet FM3DMagnetoStaticExtended::getIndex(const Vector_t &X) const {
    IndexTriplet idx;
    idx.i = std::floor((X(0) - xbegin_m) / hx_m);
    idx.j = std::floor(std::abs(X(1)) / hy_m);
    idx.k = std::floor((X(2) - zbegin_m) / hz_m);
    PAssert_LT(idx.i, num_gridpx_m - 1);
    PAssert_LT(idx.j, num_gridpy_m - 1);
    PAssert_LT(idx.k, num_gridpz_m - 1);

    idx.weight(0) = (X(0) - xbegin_m) / hx_m - idx.i;
    idx.weight(1) = std::abs(X(1)) / hy_m - idx.j;
    idx.weight(2) = (X(2) - zbegin_m) / hz_m - idx.k;

    return idx;
}

#endif
