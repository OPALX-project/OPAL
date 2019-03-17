#include "TrimCoils/TrimCoilBFit.h"

TrimCoilBFit::TrimCoilBFit(double bmax,
                           double rmin,
                           double rmax,
                           const std::vector<double>& coefnum,
                           const std::vector<double>& coefdenom,
                           const std::vector<double>& coefnumphi,
                           const std::vector<double>& coefdenomphi):
    TrimCoilFit(bmax, rmin, rmax, coefnum, coefdenom, coefnumphi, coefdenomphi)
{}

void TrimCoilBFit::doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double btr = 0.0, dr = 0.0;
    calculateRationalFunction(RADIUS, r, btr, dr);
    double phi = 0.0, dphi = 0.0;
    calculateRationalFunction(PHI,    phi_rad, phi, dphi);

    //std::cout << "r " << r << " dr " <<  dr << std::endl;

    *bz += bmax_m * btr *  phi;
    *br += bmax_m * (dr *  phi + btr*dphi) * z;
}