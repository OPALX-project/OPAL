#ifndef CLASSIC_AstraFIELDMAP1DDYNAMIC_HH
#define CLASSIC_AstraFIELDMAP1DDYNAMIC_HH

#include "Fields/Fieldmap.h"

class Astra1DDynamic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

    virtual bool isInside(const Vector_t &r) const;
private:
    Astra1DDynamic(std::string aFilename);
    ~Astra1DDynamic();

    virtual void readMap();
    virtual void freeMap();

    double *__restrict__ FourCoefs_m;

    double frequency_m;
    double xlrep_m;

    double zbegin_m;
    double zend_m;
    double length_m;

    int accuracy_m;
    int num_gridpz_m;

    friend class Fieldmap;
};

inline bool Astra1DDynamic::isInside(const Vector_t &r) const
{
    return r(2) >= zbegin_m && r(2) < zend_m;
}

#endif
