//
// Class BeamStripping
//   Defines the abstract interface environment for
//   beam stripping physics.
//
// Copyright (c) 2018-2019, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/BeamStripping.h"
#include "AbsBeamline/Cyclotron.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Solvers/BeamStrippingPhysics.hh"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Utilities/LogicalError.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"
#include "Utilities/GeneralClassicException.h"

#include <memory>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "Ippl.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

#define CHECK_BSTP_FSCANF_EOF(arg) if(arg == EOF)\
throw GeneralClassicException("BeamStripping::getPressureFromFile",\
                              "fscanf returned EOF at " #arg);

extern Inform *gmsg;

// Class BeamStripping
// ------------------------------------------------------------------------

BeamStripping::BeamStripping():
    BeamStripping("")
{}


BeamStripping::BeamStripping(const BeamStripping &right):
    Component(right),
    gas_m(right.gas_m),
    pressure_m(right.pressure_m),
    pmapfn_m(right.pmapfn_m),
    pscale_m(right.pscale_m),
    temperature_m(right.temperature_m),
    stop_m(right.stop_m),
    minr_m(right.minr_m),
    maxr_m(right.maxr_m),
    minz_m(right.minz_m),
    maxz_m(right.maxz_m),
    parmatint_m(NULL) {
}

BeamStripping::BeamStripping(const std::string &name):
    Component(name),
    gas_m(""),
    pressure_m(0.0),
    pmapfn_m(""),
    pscale_m(0.0),
    temperature_m(0.0),
    stop_m(true),
    minr_m(0.0),
    maxr_m(0.0),
    minz_m(0.0),
    maxz_m(0.0),
    parmatint_m(NULL) {
}


BeamStripping::~BeamStripping() {
    if (online_m)
        goOffline();
}

void BeamStripping::accept(BeamlineVisitor &visitor) const {
    visitor.visitBeamStripping(*this);
}

void BeamStripping::setPressure(double pressure) {
    pressure_m = pressure;
}

double BeamStripping::getPressure() const {
    if(pressure_m > 0.0)
        return pressure_m;
    else {
        throw LogicalError("BeamStripping::getPressure",
                           "Pressure must not be zero");
    }
}

void BeamStripping::setTemperature(double temperature) {
    temperature_m = temperature;
}

double BeamStripping::getTemperature() const {
    if(temperature_m > 0.0)
        return temperature_m;
    else {
        throw LogicalError("BeamStripping::getTemperature",
                           "Temperature must not be zero");
    }
}

void BeamStripping::setPressureMapFN(std::string p) {
    pmapfn_m = p;
}

std::string BeamStripping::getPressureMapFN() const {
    return pmapfn_m;
}

void BeamStripping::setPScale(double ps) {
    pscale_m = ps;
}

double BeamStripping::getPScale() const {
    if(pscale_m > 0.0)
        return pscale_m;
    else {
        throw LogicalError("BeamStripping::getPScale",
                           "PScale must be positive");
    }
}

void BeamStripping::setResidualGas(std::string gas) {
    gas_m = gas;
}

std::string BeamStripping::getResidualGas() const {
    if(gas_m == "H2" || gas_m == "AIR")
        return gas_m;
    else {
        throw GeneralClassicException("BeamStripping::getResidualGas",
                                      "Residual gas " + gas_m + " not found");
    }
}

void BeamStripping::setStop(bool stopflag) {
    stop_m = stopflag;
}

bool BeamStripping::getStop() const {
    return stop_m;
}


bool BeamStripping::checkBeamStripping(PartBunchBase<double, 3> *bunch, Cyclotron* cycl,
                                       const int /*turnnumber*/, const double /*t*/, const double /*tstep*/) {

    bool flagNeedUpdate = false;

    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);
    std::pair<Vector_t, double> boundingSphere;
    boundingSphere.first = 0.5 * (rmax + rmin);
    boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);

    maxr_m = 1000 * cycl->getMaxR();
    minr_m = 1000 * cycl->getMinR();
    maxz_m = 1000 * cycl->getMaxZ();
    minz_m = 1000 * cycl->getMinZ();

    int pflag = 0;
    size_t tempnum = bunch->getLocalNum();
    for (unsigned int i = 0; i < tempnum; ++i) {
        pflag = checkPoint(bunch->R[i](0), bunch->R[i](1), bunch->R[i](2));
        if ( (pflag != 0) && (bunch->Bin[i] != -1) )
            flagNeedUpdate = true;
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if (flagNeedUpdate && parmatint_m) {
        dynamic_cast<BeamStrippingPhysics*>(parmatint_m)->setCyclotron(cycl);
        parmatint_m->apply(bunch, boundingSphere);
    }
    return flagNeedUpdate;
}


void BeamStripping::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    initialise(bunch, pscale_m);
}

void BeamStripping::initialise(PartBunchBase<double, 3> *bunch, const double &scaleFactor) {

    RefPartBunch_m = bunch;

    parmatint_m = getParticleMatterInteraction();

    if (pmapfn_m != "") {
        getPressureFromFile(scaleFactor);
        // calculate the radii of initial grid.
        initR(PP.rmin, PP.delr, PField.nrad);
    }

    goOnline(-1e6);

    // change the mass and charge to simulate real particles
    //*gmsg << "* Mass and charge have been reseted for beam stripping " <<  endl;
    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->M[i] = bunch->getM()*1E-9;
        bunch->Q[i] = bunch->getQ() * Physics::q_e;
        if(bunch->weHaveBins())
            bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];
    }
}

void BeamStripping::finalise() {
    if (online_m)
        goOffline();
    *gmsg << "* Finalize Beam Stripping" << endl;
}

void BeamStripping::goOnline(const double &) {
}

void BeamStripping::goOffline() {
    online_m = false;
}

bool BeamStripping::bends() const {
    return false;
}

void BeamStripping::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {}

ElementBase::ElementType BeamStripping::getType() const {
    return BEAMSTRIPPING;
}

std::string BeamStripping::getBeamStrippingShape() {
    return "BeamStripping";
}

int BeamStripping::checkPoint(const double &x, const double &y, const double &z) {
    int cn;
    double rpos = std::sqrt(x * x + y * y);
    double zpos = z;
    if (zpos >= maxz_m || zpos <= minz_m || rpos >= maxr_m || rpos <= minr_m)
        cn = 0;
    else
        cn = 1;
    return (cn);  // 0 if out, and 1 if in
}

double BeamStripping::checkPressure(const double &x, const double &y) {

    double pressure = 0.0;

    if(pmapfn_m != "") {
        const double rad = std::sqrt(x * x + y * y);
        const double xir = (rad - PP.rmin) / PP.delr;

        // ir : the number of path whose radius is less than the 4 points of cell which surround the particle.
        const int ir = (int)xir;
        // wr1 : the relative distance to the inner path radius
        const double wr1 = xir - (double)ir;
        // wr2 : the relative distance to the outer path radius
        const double wr2 = 1.0 - wr1;

        const double tempv = atan(y / x);
        double tet = tempv;
        if((x < 0) && (y >= 0)) tet = Physics::pi + tempv;
        else if((x < 0) && (y <= 0)) tet = Physics::pi + tempv;
        else if((x > 0) && (y <= 0)) tet = 2.0 * Physics::pi + tempv;
        else if((x == 0) && (y > 0)) tet = Physics::pi / 2.0;
        else if((x == 0) && (y < 0)) tet = 1.5 * Physics::pi;

        // the actual angle of particle
        tet = tet / Physics::pi * 180.0;

        // the corresponding angle on the field map
        // Note: this does not work if the start point of field map does not equal zero.
        double xit = tet / PP.dtet;
        int it = (int) xit;
        const double wt1 = xit - (double)it;
        const double wt2 = 1.0 - wt1;
        // it : the number of point on the inner path whose angle is less than the particle' corresponding angle.
        // include zero degree point
        it++;
        double epsilon = 0.06;
        if  (tet > 360 - epsilon && tet < 360 + epsilon) it = 0;

        int r1t1, r2t1, r1t2, r2t2;
        // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
        // considering  the array start with index of zero, minus 1.

        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);

        if((it >= 0) && (ir >= 0) && (it < PField.ntetS) && (ir < PField.nrad)) {
            pressure = (PField.pfld[r1t1] * wr2 * wt2 +
                        PField.pfld[r2t1] * wr1 * wt2 +
                        PField.pfld[r1t2] * wr2 * wt1 +
                        PField.pfld[r2t2] * wr1 * wt1);
            if(pressure <= 0.0) {
                *gmsg << level4 << getName() << ": Pressure data from file zero." << endl;
                *gmsg << level4 << getName() << ": Take constant value through BeamStripping::getPressure" << endl;
                pressure = getPressure();
            }
        }
        else if (ir >= PField.nrad) {
            *gmsg << level4 << getName() << ": Particle out of maximum radial position of pressure field map." << endl;
            *gmsg << level4 << getName() << ": Take constant value through BeamStripping::getPressure" << endl;
            pressure = getPressure();
        }
        else {
            throw GeneralClassicException("BeamStripping::checkPressure",
                                          "Pressure data not found");
        }
    }
    else {
        pressure = getPressure();
    }

    return pressure;
}

// Calculates radius of initial grid (dimensions in [m]!)
void BeamStripping::initR(double rmin, double dr, int nrad) {
    PP.rarr.resize(nrad);
    for(int i = 0; i < nrad; i++)
        PP.rarr[i] = rmin + i * dr;

    PP.delr = dr;
}

// Read pressure map from external file.
void BeamStripping::getPressureFromFile(const double &scaleFactor) {

    FILE *f = NULL;

    *gmsg << "* " << endl;
    *gmsg << "* Reading pressure field map " << pmapfn_m << endl;

    PP.Pfact = scaleFactor;

    if((f = fopen(pmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("BeamStripping::getPressureFromFile",
                                      "failed to open file '" + pmapfn_m + "', please check if it exists");
    }

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%lf", &PP.rmin));
    *gmsg << "* --- Minimal radius of measured pressure map: " << PP.rmin << " [mm]" << endl;

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%lf", &PP.delr));
    //if the value is negative, the actual value is its reciprocal.
    if(PP.delr < 0.0) PP.delr = 1.0 / (-PP.delr);
    *gmsg << "* --- Stepsize in radial direction: " << PP.delr << " [mm]" << endl;

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%lf", &PP.tetmin));
    *gmsg << "* --- Minimal angle of measured pressure map: " << PP.tetmin << " [deg]" << endl;

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%lf", &PP.dtet));
    //if the value is negative, the actual value is its reciprocal.
    if(PP.dtet < 0.0) PP.dtet = 1.0 / (-PP.dtet);
    *gmsg << "* --- Stepsize in azimuthal direction: " << PP.dtet << " [deg]" << endl;

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%d", &PField.ntet));
    *gmsg << "* --- Grid points along azimuth (ntet): " << PField.ntet << endl;

    CHECK_BSTP_FSCANF_EOF(fscanf(f, "%d", &PField.nrad));
    *gmsg << "* --- Grid points along radius (nrad): " << PField.nrad << endl;
    *gmsg << "* --- Maximum radial position: " << PP.rmin + (PField.nrad-1)*PP.delr << " [mm]" << endl;
    PP.rmin *= 0.001;  // mm --> m
    PP.delr *= 0.001;  // mm --> m

    PField.ntetS = PField.ntet + 1;
    *gmsg << "* --- Adding a guard cell along azimuth" << endl;

    PField.ntot = PField.nrad * PField.ntetS;
    PField.pfld.resize(PField.ntot);
    *gmsg << "* --- Total stored grid point number ((ntet+1) * nrad) : " << PField.ntot << endl;
    *gmsg << "* --- Escale factor: " << PP.Pfact << endl;

    for(int i = 0; i < PField.nrad; i++) {
        for(int k = 0; k < PField.ntet; k++) {
            CHECK_BSTP_FSCANF_EOF(fscanf(f, "%16lE", &(PField.pfld[idx(i, k)])));
            PField.pfld[idx(i, k)] *= PP.Pfact;
        }
    }
    *gmsg << "*" << endl;

    fclose(f);
}

#undef CHECK_BSTP_FSCANF_EOF