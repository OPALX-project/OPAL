#ifndef CLASSIC_AstraFIELDMAP1DMAGNETOSTATICFAST_HH
#define CLASSIC_AstraFIELDMAP1DMAGNETOSTATICFAST_HH

#include "Fields/Astra1D_fast.h"

class Astra1DMagnetoStatic_fast: public Astra1D_fast {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

private:
    Astra1DMagnetoStatic_fast(std::string aFilename);
    virtual ~Astra1DMagnetoStatic_fast();

    virtual void readMap();

    bool readFileHeader(std::ifstream &file);
    int stripFileHeader(std::ifstream &file);

    friend class Fieldmap;
};

#endif