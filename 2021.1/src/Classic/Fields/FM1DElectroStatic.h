#ifndef CLASSIC_FIELDMAP1DELECTROSTATIC_HH
#define CLASSIC_FIELDMAP1DELECTROSTATIC_HH

#include "Fields/Fieldmap.h"

class FM1DElectroStatic: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal,
                                    double &yIni, double &yFinal,
                                    double &zIni, double &zFinal) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E,
                                    Vector_t &B, const DiffDirection &dir) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    virtual bool isInside(const Vector_t &r) const;
private:
    FM1DElectroStatic(std::string aFilename);
    ~FM1DElectroStatic();

    virtual void readMap();
    virtual void freeMap();

    bool checkFileData(std::ifstream &fieldFile, bool parsingPassed);
    void computeFieldOffAxis(const Vector_t &R, Vector_t &E, Vector_t &B,
                             std::vector<double> fieldComponents) const;
    void computeFieldOnAxis(double z,
                            std::vector<double> &fieldComponents) const;
    void computeFourierCoefficients(double maxEz, double fieldData[]);
    void convertHeaderData();
    double readFileData(std::ifstream &fieldFile, double fieldData[]);
    bool readFileHeader(std::ifstream &fieldFile);
    void stripFileHeader(std::ifstream &fieldFile);

    double rBegin_m;                        /// Minimum radius of field.
    double rEnd_m;                          /// Maximum radius of field.
    double zBegin_m;                        /// Longitudinal start of field.
    double zEnd_m;                          /// Longitudinal end of field.
    double length_m;                        /// Field length.
    int numberOfGridPoints_m;               /// Number of grid points in field input file.

    int accuracy_m;                         /// Number of Fourier coefficients to use reconstructing field.
    std::vector<double> fourierCoefs_m;     /// Fourier coefficients derived from field map.

    friend class Fieldmap;
};

inline bool FM1DElectroStatic::isInside(const Vector_t &r) const
{
    return r(2) >= zBegin_m && r(2) < zEnd_m;
}

#endif
