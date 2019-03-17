#include "TrimCoil.h"

#include <cmath>

#include "Physics/Physics.h"

TrimCoil::TrimCoil(double bmax,
                   double rmin,
                   double rmax)
{
    // convert to m
    const double mm2m = 0.001;
    rmin_m = rmin * mm2m;
    rmax_m = rmax * mm2m;
    // convert to kG
    bmax_m = bmax * 10.0;
}

void TrimCoil::applyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    if (std::abs(bmax_m) < 1e-20) return;
    if ((phimin_m <= phimax_m && (phi_rad < phimin_m || phi_rad > phimax_m)) ||
        (phimin_m >  phimax_m && (phi_rad < phimin_m && phi_rad > phimax_m)) ) return;

    doApplyField(r,z,phi_rad,br,bz);
}

void TrimCoil::setAzimuth(const double phimin, const double phimax)
{
    // phi convert to rad in [0,two pi]
    if (phimin_m < 0) phimin_m += 360;
    if (phimax_m < 0) phimax_m += 360;

    phimin_m = phimin * Physics::deg2rad;
    phimax_m = phimax * Physics::deg2rad;
}
