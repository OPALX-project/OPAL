//
// Class Vacuum
//   Defines the abstract interface for vacuum.
//
// Copyright (c) 2018 - 2021, Pedro Calvo, CIEMAT, Spain
// All rights reserved.
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
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Solvers/BeamStrippingPhysics.h"
#include "Solvers/ParticleMatterInteractionHandler.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/LogicalError.h"

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#define CHECK_VAC_FSCANF_EOF(arg) if (arg == EOF)\
throw GeneralClassicException("Vacuum::getPressureFromFile",\
                              "fscanf returned EOF at " #arg);

extern Inform *gmsg;

const boost::bimap<ResidualGas, std::string> Vacuum::bmResidualGasString_s =
    boost::assign::list_of<const boost::bimap<ResidualGas, std::string>::relation>
        (ResidualGas::NOGAS, "NOGAS")
        (ResidualGas::AIR,   "AIR")
        (ResidualGas::H2,    "H2");


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
    parmatint_m(nullptr)
{}


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
    parmatint_m(nullptr)
{}


Vacuum::~Vacuum() {
    if (online_m)
        goOffline();
}


void Vacuum::accept(BeamlineVisitor& visitor) const {
    visitor.visitVacuum(*this);
}

void Vacuum::setResidualGas(std::string gas) {
    auto it = bmResidualGasString_s.right.find(gas);
    if (it != bmResidualGasString_s.right.end()) {
        gas_m = it->second;
    } else {
        gas_m = ResidualGas::NOGAS;
    }
}

ResidualGas Vacuum::getResidualGas() const {
    return gas_m;
}

std::string Vacuum::getResidualGasName() {
    if (gas_m != ResidualGas::NOGAS) {
        return bmResidualGasString_s.left.at(gas_m);
    } else {
        throw GeneralClassicException("Vacuum::getResidualGasName",
                                      "Residual gas not set");
    }
}

void Vacuum::setPressure(double pressure) {
    pressure_m = pressure;
}

double Vacuum::getPressure() const {
    if (pressure_m > 0.0) {
        return pressure_m;
    } else {
        throw LogicalError("Vacuum::getPressure",
                           "Pressure must be positive");
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

void Vacuum::setTemperature(double temperature) {
    temperature_m = temperature;
}

double Vacuum::getTemperature() const {
    if (temperature_m > 0.0) {
        return temperature_m;
    } else {
        throw LogicalError("Vacuum::getTemperature",
                           "Temperature must be positive");
    }
}

void Vacuum::setStop(bool stopflag) {
    stop_m = stopflag;
}

bool Vacuum::getStop() const {
    return stop_m;
}

bool Vacuum::checkPoint(const Vector_t& R) {
    bool hit = false;
    if (OpalData::getInstance()->isInOPALCyclMode()) {
        double rpos = std::sqrt(R(0) * R(0) + R(1) * R(1));
        if (R(2) <= maxz_m || R(2) >= minz_m || rpos <= maxr_m || rpos >= minr_m) {
            hit = true;
        }
    } else {
        hit = ((R(2) > 0.0) && (R(2) <= getElementLength()));
    }
    return hit;
}

bool Vacuum::apply(const size_t& /*i*/, const double& /*t*/, Vector_t& /*E*/, Vector_t& /*B*/) {
    return false;
}

bool Vacuum::applyToReferenceParticle(const Vector_t& /*R*/, const Vector_t& /*P*/,
                                      const double& /*t*/, Vector_t& /*E*/, Vector_t& /*B*/) {
    return false;
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
        int pflag = checkPoint(bunch->R[i]);
        if ( (pflag != 0) && (bunch->Bin[i] != -1) )
            flagNeedUpdate = true;
    }

    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

    if (flagNeedUpdate && parmatint_m) {
        dynamic_cast<BeamStrippingPhysics*>(parmatint_m)->setCyclotron(cycl);
        parmatint_m->apply(bunch, boundingSphere);
        parmatint_m->print(*gmsg);
    }
    return flagNeedUpdate;
}

void Vacuum::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {

    endField = startField + getElementLength();

    initialise(bunch);

    print();
}

void Vacuum::initialise(PartBunchBase<double, 3>* bunch) {

    RefPartBunch_m = bunch;

    parmatint_m = getParticleMatterInteraction();

    goOnline(-1e6);

    updateParticleAttributes();

    if (boost::filesystem::exists(pmapfn_m)) {
        getPressureFromFile();
        // calculate the radii of initial grid.
        initR(PP_m.rmin_m, PP_m.delr_m, PField_m.nrad_m);
    }
}

void Vacuum::updateParticleAttributes() {
    for (size_t i = 0; i < RefPartBunch_m->getLocalNum(); ++i) {
        RefPartBunch_m->M[i] = RefPartBunch_m->getM() * Units::eV2GeV;
        RefPartBunch_m->Q[i] = RefPartBunch_m->getQ() * Physics::q_e;
    }
}

void Vacuum::finalise() {
    *gmsg << "* Finalize vacuum " << getName() << endl;
    if (online_m)
        goOffline();
}

void Vacuum::goOnline(const double&) {
    online_m = true;
}

void Vacuum::print() {
    *gmsg << level2 << "\n" << parmatint_m->getElement()->getTypeString()
          << ": " << getName() << " -> Residual gas   = "
          << getResidualGasName() << endl;
    *gmsg << level2 << parmatint_m->getElement()->getTypeString()
          << ": " << getName() << " -> Pressure level = "
          << std::scientific << pressure_m << " [mbar]" << endl;
}

void Vacuum::goOffline() {
    online_m = false;
}

void Vacuum::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementType Vacuum::getType() const {
    return ElementType::VACUUM;
}

double Vacuum::checkPressure(const Vector_t& R) {

    const double x = R(0);
    const double y = R(1);
    double pressure = 0.0;

    if (pmapfn_m.empty()) {
        pressure = getPressure();
    } else {
        const double rad = std::sqrt(x * x + y * y);
        const double xir = (rad - PP_m.rmin_m) / PP_m.delr_m;

        // ir: the number of path whose radius is less than the 4 points of cell which surround the particle.
        const int ir = (int)xir;
        // wr1: the relative distance to the inner path radius
        const double wr1 = xir - (double)ir;
        // wr2: the relative distance to the outer path radius
        const double wr2 = 1.0 - wr1;

        const double tempv = std::atan(y / x);
        double tet = tempv;
        if      ((x < 0) && (y >= 0)) tet = Physics::pi + tempv;
        else if ((x < 0) && (y <= 0)) tet = Physics::pi + tempv;
        else if ((x > 0) && (y <= 0)) tet = Physics::two_pi + tempv;
        else if ((x == 0) && (y > 0)) tet = Physics::pi / 2.0;
        else if ((x == 0) && (y < 0)) tet = 1.5 * Physics::pi;

        // the actual angle of particle
        tet = tet * Units::rad2deg;

        // the corresponding angle on the field map
        // Note: this does not work if the start point of field map does not equal zero.
        double xit = tet / PP_m.dtet_m;
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

        if ((it >= 0) && (ir >= 0) && (it < PField_m.ntetS_m) && (ir < PField_m.nrad_m)) {
            pressure = (PField_m.pfld_m[r1t1] * wr2 * wt2 +
                        PField_m.pfld_m[r2t1] * wr1 * wt2 +
                        PField_m.pfld_m[r1t2] * wr2 * wt1 +
                        PField_m.pfld_m[r2t2] * wr1 * wt1);
            if (pressure <= 0.0) {
                *gmsg << level4 << getName() << ": Pressure data from file zero." << endl;
                *gmsg << level4 << getName() << ": Take constant value through Vacuum::getPressure" << endl;
                pressure = getPressure();
            }
        } else if (ir >= PField_m.nrad_m) {
            *gmsg << level4 << getName() << ": Particle out of maximum radial position of pressure field map." << endl;
            *gmsg << level4 << getName() << ": Take constant value through Vacuum::getPressure" << endl;
            pressure = getPressure();
        } else {
            throw GeneralClassicException("Vacuum::checkPressure",
                                          "Pressure data not found");
        }
    }
    return pressure;
}

// Calculates radius of initial grid (dimensions in [m]!)
void Vacuum::initR(double rmin_m, double dr, int nrad_m) {
    PP_m.rarr_m.resize(nrad_m);
    for (int i = 0; i < nrad_m; i++) {
        PP_m.rarr_m[i] = rmin_m + i * dr;
    }
    PP_m.delr_m = dr;
}

// Read pressure map from external file.
void Vacuum::getPressureFromFile() {

    *gmsg << "* " << endl;
    *gmsg << "* Reading pressure field map " << endl;

    PP_m.Pfact_m = pscale_m;
    FILE* f = nullptr;
    if ((f = std::fopen(pmapfn_m.c_str(), "r")) == nullptr) {
        throw GeneralClassicException("Vacuum::getPressureFromFile",
                                      "failed to open file '" + pmapfn_m +
                                      "', please check if it exists");
    }

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP_m.rmin_m));
    *gmsg << "* --- Minimal radius of measured pressure map: " << PP_m.rmin_m << " [mm]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP_m.delr_m));
    //if the value is negative, the actual value is its reciprocal.
    if (PP_m.delr_m < 0.0) PP_m.delr_m = 1.0 / (-PP_m.delr_m);
    *gmsg << "* --- Stepsize in radial direction: " << PP_m.delr_m << " [mm]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP_m.tetmin_m));
    *gmsg << "* --- Minimal angle of measured pressure map: " << PP_m.tetmin_m << " [deg]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%lf", &PP_m.dtet_m));
    //if the value is negative, the actual value is its reciprocal.
    if (PP_m.dtet_m < 0.0) PP_m.dtet_m = 1.0 / (-PP_m.dtet_m);
    *gmsg << "* --- Stepsize in azimuthal direction: " << PP_m.dtet_m << " [deg]" << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%d", &PField_m.ntet_m));
    *gmsg << "* --- Grid points along azimuth (ntet_m): " << PField_m.ntet_m << endl;

    CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%d", &PField_m.nrad_m));
    *gmsg << "* --- Grid points along radius (nrad_m): " << PField_m.nrad_m << endl;
    *gmsg << "* --- Maximum radial position: " << PP_m.rmin_m + (PField_m.nrad_m-1)*PP_m.delr_m << " [mm]" << endl;
    PP_m.rmin_m *= Units::mm2m;
    PP_m.delr_m *= Units::mm2m;

    PField_m.ntetS_m = PField_m.ntet_m + 1;
    *gmsg << "* --- Adding a guard cell along azimuth" << endl;

    PField_m.ntot_m = PField_m.nrad_m * PField_m.ntetS_m;
    PField_m.pfld_m.resize(PField_m.ntot_m);
    *gmsg << "* --- Total stored grid point number ((ntet_m+1) * nrad_m) : " << PField_m.ntot_m << endl;
    *gmsg << "* --- Escale factor: " << PP_m.Pfact_m << endl;

    for (int i = 0; i < PField_m.nrad_m; i++) {
        for (int k = 0; k < PField_m.ntet_m; k++) {
            CHECK_VAC_FSCANF_EOF(std::fscanf(f, "%16lE", &(PField_m.pfld_m[idx(i, k)])));
            PField_m.pfld_m[idx(i, k)] *= PP_m.Pfact_m;
        }
    }

    std::fclose(f);
}

#undef CHECK_VAC_FSCANF_EOF
