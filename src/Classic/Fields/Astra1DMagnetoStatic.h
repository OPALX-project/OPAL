#ifndef CLASSIC_AstraFIELDMAP1DMAGNETOSTATIC_HH
#define CLASSIC_AstraFIELDMAP1DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"

class _Astra1DMagnetoStatic: public _Fieldmap {

public:
    virtual ~_Astra1DMagnetoStatic();

    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    virtual bool isInside(const Vector_t &r) const;
private:
    _Astra1DMagnetoStatic(const std::string& filename);

    static Astra1DMagnetoStatic create(const std::string& filename);

    virtual void readMap();
    virtual void freeMap();

    double *FourCoefs_m;

    double zbegin_m;
    double zend_m;
    double length_m;

    int accuracy_m;
    int num_gridpz_m;

    friend class _Fieldmap;
};

inline bool _Astra1DMagnetoStatic::isInside(const Vector_t &r) const
{
    return r(2) >= zbegin_m && r(2) < zend_m;
}

using Astra1DMagnetoStatic = std::shared_ptr<_Astra1DMagnetoStatic>;

#endif