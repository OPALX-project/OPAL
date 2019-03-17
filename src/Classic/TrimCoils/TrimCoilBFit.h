#ifndef TRIM_COILBFIT_H
#define TRIM_COILBFIT_H

#include "TrimCoils/TrimCoilFit.h"

#include <vector>

/// TrimCoilBFit class
/// General rational function fit
/// https://gitlab.psi.ch/OPAL/src/issues/157

class TrimCoilBFit : public TrimCoilFit {

public:
    TrimCoilBFit(double bmax,
                 double rmin,
                 double rmax,
                 const std::vector<double>& coefnum,
                 const std::vector<double>& coefdenom,
                 const std::vector<double>& coefnumphi,
                 const std::vector<double>& coefdenomphi);

    virtual ~TrimCoilBFit() { };

private:
    TrimCoilBFit() = delete;

    /// @copydoc TrimCoil::doApplyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz);
};

#endif //TRIM_COILBFIT_H
