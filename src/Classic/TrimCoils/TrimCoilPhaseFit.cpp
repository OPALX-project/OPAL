#include "TrimCoils/TrimCoilPhaseFit.h"

TrimCoilPhaseFit::TrimCoilPhaseFit(double bmax,
                                   double rmin,
                                   double rmax,
                                   const std::vector<double>& coefnum,
                                   const std::vector<double>& coefdenom,
                                   const std::vector<double>& coefnumphi,
                                   const std::vector<double>& coefdenomphi):
    TrimCoilFit(bmax, rmin, rmax, coefnum, coefdenom, coefnumphi, coefdenomphi)
{}

void TrimCoilPhaseFit::doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz)
{
    // check range
    if (r < rmin_m || r > rmax_m) return;

    double phase=0.0, d_phase=0.0, d2_phase=0.0;
    calculateRationalFunction(RADIUS, r, phase, d_phase, d2_phase);
    double phi = 0.0, d_phi = 0.0, d2_phi=0.0;
    calculateRationalFunction(PHI,    phi_rad, phi, d_phi, d2_phi);

    //std::cout << "r " << r << " dr " <<  dr << std::endl;
    double derivative =  d_phase * phi + phase *  d_phi;
    double der2       = d2_phase * phi + phase * d2_phi + 2*d_phi*d_phase;

    *bz += - bmax_m * derivative;
    *br += - bmax_m * der2 * z;
}