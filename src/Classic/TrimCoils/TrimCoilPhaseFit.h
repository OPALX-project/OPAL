#ifndef TRIM_COILPHASEFIT_H
#define TRIM_COILPHASEFIT_H

#include "TrimCoils/TrimCoilFit.h"

#include <vector>

/// TrimCoilPhaseFit class
/// General rational function fit of the phase shift

class TrimCoilPhaseFit : public TrimCoilFit {

public:
    TrimCoilPhaseFit(double bmax,
                     double rmin,
                     double rmax,
                     const std::vector<double>& coefnum,
                     const std::vector<double>& coefdenom,
                     const std::vector<double>& coefnumphi,
                     const std::vector<double>& coefdenomphi);

    virtual ~TrimCoilPhaseFit() { };

private:
    TrimCoilPhaseFit() = delete;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz);
};

#endif //TRIM_COILPHASEFIT_H
