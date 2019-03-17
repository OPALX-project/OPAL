#ifndef TRIM_COIL_H
#define TRIM_COIL_H

#include "Physics/Physics.h"

/// Abstract TrimCoil class

class TrimCoil {

public:

    TrimCoil(double bmax, double rmin, double rmax);
    /// Apply the trim coil at position r and z to Bfields br and bz
    /// Calls virtual doApplyField
    void applyField(const double r, const double z, const double phi_rad, double *br, double *bz);
    /// Set azimuthal range
    void setAzimuth(const double phimin, const double phimax);
    virtual ~TrimCoil() { };

protected:

    /// Maximum B field (kG)
    double bmax_m;
    /// Minimum radius (m)
    double rmin_m;
    /// Maximum radius (m)
    double rmax_m;
    /// Minimal azimuth (rad)
    double phimin_m = 0.0;
    /// Maximal azimuth (rad)
    double phimax_m = Physics::two_pi;

private:

    /// virtual implementation of applyField
    virtual void doApplyField(const double r, const double z, const double phi_rad, double *br, double *bz) = 0;
};

#endif //TRIM_COIL_H
