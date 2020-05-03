#ifndef CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH
#define CLASSIC_FIELDMAP2DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

class FM2DMagnetoStatic: public Fieldmap {

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
    FM2DMagnetoStatic(std::string aFilename);
    ~FM2DMagnetoStatic();

    virtual void readMap();
    virtual void freeMap();

    double *FieldstrengthBz_m;    /**< 2D array with Ez, read in first along z0 - r0 to rN then z1 - r0 to rN until zN - r0 to rN  */
    double *FieldstrengthBr_m;    /**< 2D array with Er, read in like Ez*/

    double rbegin_m;
    double rend_m;
    double zbegin_m;
    double zend_m;
    double hz_m;                   /**< length between points in grid, z-direction, m*/
    double hr_m;                   /**< length between points in grid, r-direction, m*/
    int num_gridpr_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool swap_m;
    friend class Fieldmap;
};

inline bool FM2DMagnetoStatic::isInside(const Vector_t &r) const
{
    return r(2) >= zbegin_m && r(2) < zend_m && std::sqrt(r(0)*r(0) + r(1)*r(1)) < rend_m;
}

#endif