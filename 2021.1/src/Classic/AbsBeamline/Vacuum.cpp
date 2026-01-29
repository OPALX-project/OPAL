//
// Class Vacuum
//   Defines the abstract interface environment for
//   beam stripping physics.
//
// Copyright (c) 2018-2021, Pedro Calvo, CIEMAT, Spain
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
#include "AbsBeamline/Vacuum.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/Cyclotron.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Solvers/BeamStrippingPhysics.h"
#include "Solvers/ParticleMatterInteractionHandler.h"
#include "Utilities/LogicalError.h"
#include "Utilities/GeneralClassicException.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>

#include "Utility/Inform.h"

#define CHECK_VAC_FSCANF_EOF(arg) if (arg == EOF)\
throw GeneralClassicException("Vacuum::getPressureFromFile",\
                              "fscanf returned EOF at " #arg);

extern Inform* gmsg;


Vacuum::Vacuum():
    Vacuum("")
{}

Vacuum::Vacuum(const Vacuum& right):
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

Vacuum::Vacuum(const std::string& name):
    Component(name),
    gas_m(ResidualGas::NOGAS),
    pressure_m(0.0),
    pmapfn_m(""),
    pscale_m(1.0),
    temperature_m(0.0),
    stop_m(true),
    minr_m(0.0),
    maxr_m(0.0),
    minz_m(0.0),
    maxz_m(0.0),
    parmatint_m(NULL) {
}

Vacuum::~Vacuum() {
    if (online_m)
        goOffline();
}


void Vacuum::accept(BeamlineVisitor& visitor) const {
    visitor.visitVacuum(*this);
}

void Vacuum::setPressure(double pressure) {
    pressure_m = pressure;
}

double Vacuum::getPressure() const {
    if (pressure_m > 0.0) {
        return pressure_m;
    } else {
        throw LogicalError("Vacuum::getPressure",
                           "Pressure must not be zero");
    }
}

void Vacuum::setTemperature(double temperature) {
    temperature_m = temperature;
}

double Vacuum::getTemperature() const {
    if (temperature_m > 0.0) {
        return temperature_m;
    } else {
        throw LogicalError("Vacuum::getTemperature",
                           "Temperature must not be zero");
    }
}

void Vacuum::setPressureMapFN(std::string p) {
    pmapfn_m = p;
}

std::string Vacuum::getPressureMapFN() const {
    return pmapfn_m;
}

void Vacuum::setPScale(double ps) {
    pscale_m = ps;
}

double Vacuum::getPScale() const {
    if (pscale_m > 0.0) {
        return pscale_m;
    } else {
        throw LogicalError("Vacuum::getPScale",
                           "PScale must be positive");
    }
}

void Vacuum::setResidualGas(std::string gas) {
    if (gas == "AIR") {
        gas_m = ResidualGas::AIR;
    } else if (gas == "H2") {
        gas_m = ResidualGas::H2;
    }
}

ResidualGas Vacuum::getResidualGas() const {
    return gas_m;
}

std::string Vacuum::getResidualGasName() {
    switch (gas_m) {
        case ResidualGas::AIR: {
            return "AIR";
        }
        case ResidualGas::H2: {
            return "H2";
        }
        default: {
           throw GeneralClassicException("Vacuum::getResidualGasName",
                                         "Residual gas not set");
        }
    }
}

void Vacuum::setStop(bool stopflag) {
    stop_m = stopflag;
}

bool Vacuum::getStop() const {
    return stop_m;
}


bool Vacuum::checkVacuum(PartBunchBase<double, 3>* bunch, Cyclotron* cycl) {

    bool flagNeedUpdate = false;

    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);
    std::pair<Vector_t, double> boundingSphere;
    boundingSphere.first = 0.5 * (rmax + rmin);
    boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);

    maxr_m = cycl->getMaxR();
    minr_m = cycl->getMinR();
    maxz_m = cycl->getMaxZ();
    minz_m = cycl->getMinZ();

    size_t tempnum = bunch->getLocalNum();
    for (unsigned int i = 0; i < tempnum; ++i) {
        int pflag = checkPoint(bunch->R[i](0), bunch->R[i](1), bunch->R[i](2));
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


void Vacuum::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {
    endField = startField + getElementLength();
    initialise(bunch, pscale_m);
}

void Vacuum::initialise(PartBunchBase<double, 3>* bunch, const double& scaleFactor) {

    RefPartBunch_m = bunch;

    parmatint_m = getParticleMatterInteraction();

    if (pmapfn_m != "") {
        getPressureFromFile(scaleFactor);
        // calculate the radii of initial grid.
        initR(PP.rmin, PP.delr, PField.nrad);
    }

    goOnline(-1e6);

    // change the mass and charge to simulate real particles
    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->M[i] = bunch->getM()*1E-9;
        bunch->Q[i] = bunch->getQ() * Physics::q_e;
        if (bunch->weHaveBins())
            bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];
    }
}

void Vacuum::finalise() {
    if (online_m)
        goOffline();
    *gmsg << "* Finalize vacuum" << endl;
}

void Vacuum::goOnline(const double &) {
}

void Vacuum::goOffline() {
    online_m = false;
}

bool Vacuum::bends() const {
    return false;
}

void Vacuum::getDimensions(double& /*zBegin*/, double& /*zEnd*/) const {}

ElementBase::ElementType Vacuum::getType() const {
    return VACUUM;
}

std::string Vacuum::getVacuumShape() {
    return "Vacuum";
}

int Vacuum::checkPoint(const double& x, const double& y, const double& z) {
    int cn;
    double rpos = std::sqrt(x * x + y * y);
    if (z >= maxz_m || z <= minz_m || rpos >= maxr_m || rpos <= minr_m) {
        cn = 0;
    } else {
        cn = 1;
    }
    return (cn);  // 0 if out, and 1 if in
}

double Vacuum::checkPressure(const double& x, const double& y) {

    double pressure = 0.0;

    if (pmapfn_m != "") {
        const double rad = std::sqrt(x * x + y * y);
        const double xir = (rad - PP.rmin) / PP.delr;

        // ir : the number of path whose radius is less than the 4 points of cell which surround the particle.
        const int ir = (int)xir;
        // wr1 : the relative distance to the inner path radius
        const double wr1 = xir - (double)ir;
        // wr2 : the relative distance to the outer path radius
        const double wr2 = 1.0 - wr1;

        const double tempv = std::atan(y / x);
        double tet = tempv;
        if      ((x < 0) && (y >= 0)) tet = Physics::pi + tempv;
        else if ((x < 0) && (y <= 0)) tet = Physics::pi + tempv;
        else if ((x > 0) && (y <= 0)) tet = Physics::two_pi + tempv;
        else if ((x == 0) && (y > 0)) tet = Physics::pi / 2.0;
        else if ((x == 0) && (y < 0)) tet = 1.5 * Physics::pi;

        // the actual angle of particle
        tet = tet * Physics::rad2deg;

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
        if (tet > 360 - epsilon && tet < 360 + epsilon) it = 0;

        int r1t1, r2t1, r1t2, r2t2;
        // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
        // considering  the array start with index of zero, minus 1.

        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);

        if ((it >= 0) && (ir >= 0) && (it < PField.ntetS) && (ir < PField.nrad)) {
            pressure = (PField.pfld[r1t1] * wr2 * wt2 +
                        PField.pfld[r2t1] * wr1 * wt2 +
                        PField.pfld[r1t2] * wr2 * wt1 +
                        PField.pfld[r2t2] * wr1 * wt1);
            if (pressure <= 0.0) {
                *gmsg << level4 << getName() << ": Pressure data from file zero." << endl;
                *gmsg << level4 << getName() << ": Take constant value through Vacuum::getPressure" << endl;
                pressure = getPressure();
            }
        } else if (ir >= PField.nrad) {
            *gmsg << level4 << getName() << ": Particle out of maximum radial position of pressure field map." << endl;
            *gmsg << level4 << getName() << ": Take constant value through Vacuum::getPressure" << endl;
            pressure = getPressure();
        } else {
            throw GeneralClassicException("Vacuum::checkPressure",
                                          "Pressure data not found");
        }
    } else {
        pressure = getPressure();
    }

    return pressure;
}

// Calculates radius of initial grid (dimensions in [m]!)
void Vacuum::initR(double rmin, double dr, int nrad) {
    PP.rarr.resize(nrad);
    for (int i = 0; i < nrad; i++)
        PP.rarr[i] = rmin + i * dr;

    PP.delr = dr;
}

// Read pressure map from external file.
void Vacuum::getPressureFromFile(const double& scaleFactor) {

    FILE* f = NULL;

    *gmsg << "* " << endl;
    *gmsg << "* Reading pressure field map " << endl;

    PP.Pfact = scaleFactor;

    if ((f = std::fopen(pmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Vacuum::getPressureFromFile",
                                      "failed to open file '" + pmapfn_m +
                                      "', please check if it exists");
    }

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP.rmin));
    *gmsg << "* --- Minimal radius of measured pressure map: " << PP.rmin << " [mm]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP.delr));
    //if the value is negative, the actual value is its reciprocal.
    if (PP.delr < 0.0) PP.delr = 1.0 / (-PP.delr);
    *gmsg << "* --- Stepsize in radial direction: " << PP.delr << " [mm]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP.tetmin));
    *gmsg << "* --- Minimal angle of measured pressure map: " << PP.tetmin << " [deg]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP.dtet));
    //if the value is negative, the actual value is its reciprocal.
    if (PP.dtet < 0.0) PP.dtet = 1.0 / (-PP.dtet);
    *gmsg << "* --- Stepsize in azimuthal direction: " << PP.dtet << " [deg]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%d", &PField.ntet));
    *gmsg << "* --- Grid points along azimuth (ntet): " << PField.ntet << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%d", &PField.nrad));
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

    for (int i = 0; i < PField.nrad; i++) {
        for (int k = 0; k < PField.ntet; k++) {
            CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%16lE", &(PField.pfld[idx(i, k)])));
            PField.pfld[idx(i, k)] *= PP.Pfact;
        }
    }

    std::fclose(f);
}

#undef CHECK_VAC_FSCANF_EOF