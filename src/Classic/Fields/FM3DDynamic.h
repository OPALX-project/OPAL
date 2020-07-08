#ifndef CLASSIC_FIELDMAP3DDYNAMIC_HH
#define CLASSIC_FIELDMAP3DDYNAMIC_HH

#include "Fields/Fieldmap.h"

class FM3DDynamic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *msg);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

    virtual bool isInside(const Vector_t &r) const;
private:
    FM3DDynamic(std::string aFilename);
    ~FM3DDynamic();

    virtual void readMap();
    virtual void freeMap();

    double *FieldstrengthEz_m;    /**< 3D array with Ez */
    double *FieldstrengthEx_m;    /**< 3D array with Ex */
    double *FieldstrengthEy_m;    /**< 3D array with Ey */
    double *FieldstrengthBz_m;    /**< 3D array with Bz */
    double *FieldstrengthBx_m;    /**< 3D array with Bx */
    double *FieldstrengthBy_m;    /**< 3D array with By */

    double frequency_m;

    double xbegin_m;
    double xend_m;

    double ybegin_m;
    double yend_m;

    double zbegin_m;
    double zend_m;

    double hx_m;                   /**< length between points in grid, x-direction */
    double hy_m;                   /**< length between points in grid, y-direction */
    double hz_m;                   /**< length between points in grid, z-direction */
    unsigned int num_gridpx_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    unsigned int num_gridpy_m;              /**< Read in number of points after 0(not counted here) in grid, r-direction*/
    unsigned int num_gridpz_m;              /**< Read in number of points after 0(not counted here) in grid, z-direction*/

    bool normalize_m;
    friend class Fieldmap;
};

inline bool FM3DDynamic::isInside(const Vector_t &r) const
{
    return ((r(0) >= xbegin_m && r(0) < xend_m) &&
            (r(1) >= ybegin_m && r(1) < yend_m) &&
            (r(2) >= zbegin_m && r(2) < zend_m));
}

#endif