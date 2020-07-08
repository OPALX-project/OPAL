#ifndef CLASSIC_AstraFIELDMAP1DFAST_HH
#define CLASSIC_AstraFIELDMAP1DFAST_HH

#include "Fields/Fieldmap.h"

class Astra1D_fast: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const = 0;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const = 0;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const = 0;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const = 0;
    virtual void swap() = 0;
    virtual void getInfo(Inform *) = 0;
    virtual double getFrequency() const = 0;
    virtual void setFrequency(double freq) = 0;
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

    virtual bool isInside(const Vector_t &r) const;
protected:
    Astra1D_fast(std::string aFilename);
    virtual ~Astra1D_fast();

    virtual void readMap();
    virtual void freeMap();

    bool determineNumSamplingPoints(std::ifstream &file);
    double readFieldData(std::ifstream &file);
    void normalizeFieldData(double maxEz);
    std::vector<double> getEvenlyDistributedSamplingPoints();
    std::vector<double> interpolateFieldData(std::vector<double> &samplingPoints);
    void computeFieldDerivatives(std::vector<double> & fourierComponents, int accuracy);
    std::vector<double> computeFourierCoefficients(int accuracy, std::vector<double> &evenSampling);


    double* onAxisField_m;
    double* zvals_m;
    gsl_spline *onAxisInterpolants_m[4];
    gsl_interp_accel *onAxisAccel_m[4];

    double hz_m;

    double zbegin_m;
    double zend_m;
    double length_m;

    int num_gridpz_m;
    int numHeaderLines_m;

    friend class Fieldmap;
};

inline bool Astra1D_fast::isInside(const Vector_t &r) const
{
    return r(2) >= zbegin_m && r(2) < zend_m;
}

#endif