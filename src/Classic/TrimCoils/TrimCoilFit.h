#ifndef TRIM_COILFIT_H
#define TRIM_COILFIT_H

#include "TrimCoils/TrimCoil.h"

#include <vector>

/// Abstract TrimCoilFit class
/// General rational function fit

class TrimCoilFit : public TrimCoil {

public:
    TrimCoilFit(double bmax,
                double rmin,
                double rmax,
                const std::vector<double>& coefnum,
                const std::vector<double>& coefdenom,
                const std::vector<double>& coefnumphi,
                const std::vector<double>& coefdenomphi);

    virtual ~TrimCoilFit() {};

protected:
    enum PolynomType  {NUM, DENOM, NUMPHI, DENOMPHI};
    enum FunctionType {RADIUS=0, PHI=2};

    /// calculate rational function and its first derivative
    void calculateRationalFunction(FunctionType, double value, double& quot, double& der_quot) const;
    /// calculate rational function and its first and second derivative
    void calculateRationalFunction(FunctionType, double value, double& quot, double& der_quot, double& der2_quot) const;

private:
    TrimCoilFit() = delete;

    /// rational function coefficients
    std::vector<std::vector<double>> coefs;
};

#endif //TRIM_COILFIT_H
