#include "TrimCoils/TrimCoilMirrored.h"

#include "Utility/IpplInfo.h"
#include "Utilities/GeneralClassicException.h"

#include <cmath>

TrimCoilMirrored::TrimCoilMirrored(double bmax,
                                   double rmin,
                                   double rmax,
                                   double bslope):
    TrimCoil(bmax, rmin, rmax)
{
    // convert to m
    const double mm2m = 0.001;
    bslope_m = bslope / mm2m;
}

void TrimCoilMirrored::doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    /// update bz and br with trim coil contributions
    // for some discussion on the formulas see
    // http://doi.org/10.1103/PhysRevSTAB.14.054402
    // https://gitlab.psi.ch/OPAL/src/issues/157
    // https://gitlab.psi.ch/OPAL/src/issues/110

    // unitless constants
    const double Amax1 = 1;
    const double Amax2 = 3;
    const double Amin  = -2;
    const double x01   = 4;
    const double x02   = 8;
    const double h1    = 0.03;
    const double h2    = 0.2;
    const double justAnotherFudgeFactor = 1 / 2.78;
    // helper variables
    const double log10  = std::log(10);
    const double const3 = -(Amax1 - Amin) * h1 * log10;
    const double const4 =  (Amax2 - Amin) * h2 * log10;

    const double &tcr1      = rmin_m;
    const double &tcr2      = rmax_m;
    const double &slope     = bslope_m;
    const double &magnitude = bmax_m;

    double part1;
    double part2;

    if (2 * r < (tcr2 + tcr1)) {
        part1 = std::pow(10.0,  ((r - tcr1) * slope - x01) * h1);
        part2 = std::pow(10.0, -((r - tcr1) * slope - x02) * h2);
    } else {
        part1 = std::pow(10.0,  ((tcr2 - r) * slope - x01) * h1);
        part2 = std::pow(10.0, -((tcr2 - r) * slope - x02) * h2);
    }

    const double part1plusinv = 1 / (1 + part1);
    const double part2plusinv = 1 / (1 + part2);

    double part3 = const3 * slope * part1 * part1plusinv * part1plusinv;
    double part4 = const4 * slope * part2 * part2plusinv * part2plusinv;

    const double constmag = magnitude * justAnotherFudgeFactor;

    double dr  = constmag * (part3 + part4);
    double btr = constmag * (Amin - 1 +
                             (Amax1 - Amin) * part1plusinv +
                             (Amax2 - Amin) * part2plusinv);

    if (std::isnan(dr)  || std::isinf(dr) ||
        std::isnan(btr) || std::isinf(btr)) {
        ERRORMSG("r = " << r << " m, tcr1 = " << tcr1 << " m, tcr2 = " << tcr2 << " m\n");
        ERRORMSG("slope = " << slope << ", magnitude = " << magnitude << " kG\n");
        ERRORMSG("part1 = " << part1 << ", part2 = " << part2 << "\n");
        ERRORMSG("part3 = " << part3 << ", part4 = " << part4 << endl);
        throw GeneralClassicException("TrimCoilMirrored::doApplyField",
                                      "dr or btr yield results that are either nan or inf");
    }

    *bz -= btr;
    *br -= dr * z;
}