#ifndef CLASSIC_AstraFIELDMAP1DELECTROSTATICFAST_HH
#define CLASSIC_AstraFIELDMAP1DELECTROSTATICFAST_HH

#include "Fields/Astra1D_fast.h"

class _Astra1DElectroStatic_fast: public _Astra1D_fast {

public:
    virtual ~_Astra1DElectroStatic_fast();

    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    _Astra1DElectroStatic_fast(const std::string& filename);

    static Astra1DElectroStatic_fast create(const std::string& filename);

    virtual void readMap();

    bool readFileHeader(std::ifstream &file);
    int stripFileHeader(std::ifstream &file);

    friend class _Fieldmap;
};

using Astra1DElectroStatic_fast = std::shared_ptr<_Astra1DElectroStatic_fast>;

#endif