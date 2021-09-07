//
// Class Cyclotron
//   Defines the abstract interface for a cyclotron.
//
// Copyright (c) 2007 - 2012, Jianjun Yang and Andreas Adelmann, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2013 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#include "AbsBeamline/Cyclotron.h"

#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "TrimCoils/TrimCoil.h"
#include "Utilities/Options.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>

#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>

#define CHECK_CYC_FSCANF_EOF(arg) if (arg == EOF)\
throw GeneralClassicException("Cyclotron::getFieldFromFile",\
                              "fscanf returned EOF at " #arg);

extern Inform* gmsg;


Cyclotron::Cyclotron():
    Component() {
}

Cyclotron::Cyclotron(const Cyclotron& right):
    Component(right),
    fieldType_m(right.fieldType_m),
    fmapfn_m(right.fmapfn_m),
    rffrequ_m(right.rffrequ_m),
    rfphi_m(right.rfphi_m),
    escale_m(right.escale_m),
    superpose_m(right.superpose_m),
    symmetry_m(right.symmetry_m),
    rinit_m(right.rinit_m),
    prinit_m(right.prinit_m),
    phiinit_m(right.phiinit_m),
    zinit_m(right.zinit_m),
    pzinit_m(right.pzinit_m),
    spiralFlag_m(right.spiralFlag_m),
    trimCoilThreshold_m(right.trimCoilThreshold_m),
    typeName_m(right.typeName_m),
    harm_m(right.harm_m),
    bscale_m(right.bscale_m),
    trimcoils_m(right.trimcoils_m),
    minr_m(right.minr_m),
    maxr_m(right.maxr_m),
    minz_m(right.minz_m),
    maxz_m(right.maxz_m),
    fmLowE_m(right.fmLowE_m),
    fmHighE_m(right.fmHighE_m),
    RFfilename_m(right.RFfilename_m),
    RFFCoeff_fn_m(right.RFFCoeff_fn_m),
    RFVCoeff_fn_m(right.RFVCoeff_fn_m) {
}

Cyclotron::Cyclotron(const std::string& name):
    Component(name) {
}

Cyclotron::~Cyclotron() {
}


void Cyclotron::applyTrimCoil_m(const double r, const double z,
                                const double tet_rad,
                                double* br, double* bz) {
     for (auto trimcoil : trimcoils_m) {
         trimcoil->applyField(r,tet_rad,z,br,bz);
     }
}

void Cyclotron::applyTrimCoil(const double r, const double z,
                              const double tet_rad,
                              double& br, double& bz) {
    //Changed from > to >= to include case where bz == 0 and trimCoilThreshold_m == 0 -DW
    if (std::abs(bz) >= trimCoilThreshold_m) {
        applyTrimCoil_m(r, z, tet_rad, &br, &bz);
    } else {
        // make sure to have a smooth transition
        double tmp_bz = 0.0;
        applyTrimCoil_m(r, z, tet_rad, &br, &tmp_bz);
        bz += tmp_bz * std::abs(bz) / trimCoilThreshold_m;
    }
}

void Cyclotron::accept(BeamlineVisitor &visitor) const {
    visitor.visitCyclotron(*this);
}

void Cyclotron::setRinit(double rinit) {
    rinit_m = rinit;
}

double Cyclotron::getRinit() const {
    return rinit_m;
}

void Cyclotron::setPRinit(double prinit) {
    prinit_m = prinit;
}

double Cyclotron::getPRinit() const {
    return prinit_m;
}

void Cyclotron::setPHIinit(double phiinit) {
    phiinit_m = phiinit;
}

double Cyclotron::getPHIinit() const {
    return phiinit_m;
}

void Cyclotron::setZinit(double zinit){
    zinit_m = zinit;
}

double Cyclotron::getZinit() const {
    return zinit_m;
}

void Cyclotron::setPZinit(double pzinit){
    pzinit_m = pzinit;
}

double Cyclotron::getPZinit() const {
    return pzinit_m;
}

void Cyclotron::setTrimCoilThreshold(double trimCoilThreshold) {
    trimCoilThreshold_m = 10.0 * trimCoilThreshold; // convert from Tesla to kGauss
}

double Cyclotron::getTrimCoilThreshold() const {
    return trimCoilThreshold_m;
}

void Cyclotron::setSpiralFlag(bool spiral_flag) {
    spiralFlag_m = spiral_flag;
}

bool Cyclotron::getSpiralFlag() const {
    return spiralFlag_m;
}

void Cyclotron::setFieldMapFN(std::string f) {
    fmapfn_m = f;
}

std::string Cyclotron::getFieldMapFN() const {
    if (boost::filesystem::exists(fmapfn_m)) {
        return fmapfn_m;
    } else {
        throw GeneralClassicException("Cyclotron::getFieldMapFN",
                                      "Failed to open file '" + fmapfn_m +
                                      "', please check if it exists");
    }
}

void Cyclotron::setRfFieldMapFN(std::vector<std::string> f) {
    RFfilename_m = f;
}

void Cyclotron::setRFFCoeffFN(std::vector<std::string> f) {
    RFFCoeff_fn_m = f;
}

void Cyclotron::setRFVCoeffFN(std::vector<std::string> f) {
    RFVCoeff_fn_m = f;
}

void Cyclotron::setRfPhi(std::vector<double> f) {
    rfphi_m = f;
}

double Cyclotron::getRfPhi(unsigned int i) const {
    if (i < rfphi_m.size()) {
        return rfphi_m[i];
    } else {
        throw GeneralClassicException("Cyclotron::getRfPhi",
                                      "RFPHI not defined for CYCLOTRON!");
    }
}

void Cyclotron::setRfFrequ(std::vector<double> f) {
    rffrequ_m = f;
}

double Cyclotron::getRfFrequ(unsigned int i) const {
    if (i < rffrequ_m.size()) {
        return rffrequ_m[i];
    } else {
        throw GeneralClassicException("Cyclotron::getRfFrequ",
                                      "RFFREQ not defined for CYCLOTRON!");
    }
}

void Cyclotron::setSuperpose(std::vector<bool> flag) {
  superpose_m = flag;
}

bool Cyclotron::getSuperpose(unsigned int i) const {
    if (i < superpose_m.size()) {
        return superpose_m[i];
    } else {
        throw GeneralClassicException("Cyclotron::getSuperpose",
                                      "SUPERPOSE not defined for CYCLOTRON!");
    }
}

void Cyclotron::setSymmetry(double s) {
    symmetry_m = s;
}

double Cyclotron::getSymmetry() const {
    return symmetry_m;
}

void Cyclotron::setCyclotronType(std::string t) {
    typeName_m = t;
}

const std::string &Cyclotron::getCyclotronType() const {
    return typeName_m;
}

ElementBase::ElementType Cyclotron::getType() const {
    return CYCLOTRON;
}

void Cyclotron::setCyclHarm(double h) {
    harm_m = h;
}

void Cyclotron::setBScale(double s) {
    bscale_m = s;
}

double Cyclotron::getBScale() const {
    return bscale_m;
}

void Cyclotron::setEScale(std::vector<double> s) {
    escale_m = s;
}

double Cyclotron::getEScale(unsigned int i) const {
    if (i < escale_m.size()) {
        return escale_m[i];
    } else {
        throw GeneralClassicException("Cyclotron::getEScale",
                                      "EScale not defined for CYCLOTRON!");
    }
}

unsigned int Cyclotron::getNumberOfTrimcoils() const {
  return trimcoils_m.size();
}

double Cyclotron::getCyclHarm() const {
    return harm_m;
}

double Cyclotron::getRmin() const {
    return BP_m.rmin_m;
}

double Cyclotron::getRmax() const {
    return BP_m.rmin_m + (Bfield_m.nrad_m - 1) * BP_m.delr_m;
}

void Cyclotron::setMinR(double r) {
    // DW: This is to let the user keep using mm in the input file for now
    // while switching internally to m
    minr_m = 0.001 * r;
}

void Cyclotron::setMaxR(double r) {
    // DW: This is to let the user keep using mm in the input file for now
    // while switching internally to m
    maxr_m = 0.001 * r;
}

double Cyclotron::getMinR() const {
    return minr_m;
}

double Cyclotron::getMaxR() const {
    return maxr_m;
}

void  Cyclotron::setMinZ(double z) {
    // DW: This is to let the user keep using mm in the input file for now
    // while switching internally to m
    minz_m = 0.001 * z;
}

double Cyclotron::getMinZ() const {
    return minz_m;
}

void Cyclotron::setMaxZ(double z) {
    // DW: This is to let the user keep using mm in the input file for now
    // while switching internally to m
    maxz_m = 0.001 * z;
}

double Cyclotron::getMaxZ() const {
    return maxz_m;
}

void Cyclotron::setTrimCoils(const std::vector<TrimCoil*>& trimcoils) {
    trimcoils_m = trimcoils;
}

void Cyclotron::setFMLowE(double e) {
    fmLowE_m = e;
}

double Cyclotron::getFMLowE() const {
    return fmLowE_m;
}

void Cyclotron::setFMHighE(double e) {
    fmHighE_m = e;
}

double Cyclotron::getFMHighE() const {
    return fmHighE_m;
}

void Cyclotron::setBFieldType() {
    if (typeName_m.empty()) {
        throw GeneralClassicException("Cyclotron::setBFieldType",
                                      "TYPE is not defined for CYCLOTRON!");
    } else if (typeName_m == std::string("RING")) {
        fieldType_m = BFieldType::PSIBF;
    } else if (typeName_m == std::string("CARBONCYCL")) {
        fieldType_m = BFieldType::CARBONBF;
    } else if (typeName_m == std::string("CYCIAE")) {
        fieldType_m = BFieldType::ANSYSBF;
    } else if (typeName_m == std::string("AVFEQ")) {
        fieldType_m = BFieldType::AVFEQBF;
    } else if (typeName_m == std::string("FFA")) {
        fieldType_m = BFieldType::FFABF;
    } else if (typeName_m == std::string("BANDRF")) {
        fieldType_m = BFieldType::BANDRF;
    } else if (typeName_m == std::string("SYNCHROCYCLOTRON")) {
        fieldType_m = BFieldType::SYNCHRO;
    }
}

bool Cyclotron::apply(const size_t& id, const double& t, Vector_t& E, Vector_t& B) {

    bool flagNeedUpdate = false;

    const double rpos = std::hypot(RefPartBunch_m->R[id](0), RefPartBunch_m->R[id](1));
    const double zpos = RefPartBunch_m->R[id](2);

    Inform gmsgALL("OPAL", INFORM_ALL_NODES);
    if (zpos > maxz_m || zpos < minz_m || rpos > maxr_m || rpos < minr_m) {
        flagNeedUpdate = true;
        gmsgALL << level4 << getName() << ": Particle " << id
                << " out of the global aperture of cyclotron!" << endl;
        gmsgALL << level4 << getName()
                << ": Coords: "<< RefPartBunch_m->R[id] << " m"  << endl;

    } else {
        flagNeedUpdate = apply(RefPartBunch_m->R[id], RefPartBunch_m->P[id], t, E, B);
        if (flagNeedUpdate) {
            gmsgALL << level4 << getName() << ": Particle "<< id
                    << " out of the field map boundary!" << endl;
            gmsgALL << level4 << getName()
                    << ": Coords: "<< RefPartBunch_m->R[id] << " m" << endl;
        }
    }

    if (flagNeedUpdate) {
        lossDs_m->addParticle(OpalParticle(id, RefPartBunch_m->R[id], RefPartBunch_m->P[id],
                                           t*1.0e-9, RefPartBunch_m->Q[id], RefPartBunch_m->M[id]),
                              std::make_pair(0, RefPartBunch_m->bunchNum[id]));
        RefPartBunch_m->Bin[id] = -1;
    }

    return flagNeedUpdate;
}

bool Cyclotron::apply(const Vector_t& R, const Vector_t& /*P*/,
                      const double& t, Vector_t& E, Vector_t& B) {

    const double rad   = std::hypot(R[0],R[1]);
    const double tempv = std::atan(R[1] / R[0]);
    double tet = tempv;

    /* if ((R[0] > 0) && (R[1] >= 0)) tet = tempv;
       else*/
    if      ((R[0] < 0) && (R[1] >= 0)) tet = Physics::pi + tempv;
    else if ((R[0] < 0) && (R[1] <= 0)) tet = Physics::pi + tempv;
    else if ((R[0] > 0) && (R[1] <= 0)) tet = Physics::two_pi + tempv;
    else if ((R[0] == 0) && (R[1] > 0)) tet = Physics::pi / 2.0;
    else if ((R[0] == 0) && (R[1] < 0)) tet = 1.5 * Physics::pi;

    double tet_rad = tet;

    // the actual angle of particle in degree
    tet *= Physics::rad2deg;

    // Necessary for gap phase output -DW
    if (0 <= tet && tet <= 45) waitingGap_m = 1;

    // dB_{z}/dr, dB_{z}/dtheta, B_{z}
    double brint = 0.0, btint = 0.0, bzint = 0.0;

    if ( this->interpolate(rad, tet_rad, brint, btint, bzint) ) {

        /* Br */
        double br = - brint * R[2];

        /* Btheta */
        double bt = - btint / rad * R[2];

        /* Bz */
        double bz = - bzint;

        this->applyTrimCoil(rad, R[2], tet_rad, br, bz);

        /* Br Btheta -> Bx By */
        B[0] = br * std::cos(tet_rad) - bt * std::sin(tet_rad);
        B[1] = br * std::sin(tet_rad) + bt * std::cos(tet_rad);
        B[2] = bz;
    } else {
        return true;
    }

    if (fieldType_m != BFieldType::SYNCHRO && fieldType_m != BFieldType::BANDRF) {
        return false;
    }

    //The RF field is supposed to be sampled on a cartesian grid
    std::vector<Fieldmap *>::const_iterator fi  = RFfields_m.begin();
    std::vector<double>::const_iterator rffi    = rffrequ_m.begin();
    std::vector<double>::const_iterator rfphii  = rfphi_m.begin();
    std::vector<double>::const_iterator escali  = escale_m.begin();
    std::vector<bool>::const_iterator superposei = superpose_m.begin();
    std::vector<std::vector<double>>::const_iterator rffci = rffc_m.begin();
    std::vector<std::vector<double>>::const_iterator rfvci = rfvc_m.begin();

    double xBegin(0), xEnd(0), yBegin(0), yEnd(0), zBegin(0), zEnd(0);
    int fcount = 0;

    for (; fi != RFfields_m.end(); ++fi, ++rffi, ++rfphii, ++escali, ++superposei) {
        (*fi)->getFieldDimensions(xBegin, xEnd, yBegin, yEnd, zBegin, zEnd);
        if (fcount > 0 && *superposei == false) continue;

        // Ok, this is a total patch job, but now that the internal cyclotron units are in m, we have to
        // change stuff here to match with the input units of mm in the fieldmaps. -DW
        const Vector_t temp_R = R * Vector_t(1000.0); //Keep this until we have transitioned fully to m -DW

        if (temp_R(0) < xBegin || temp_R(0) > xEnd ||
            temp_R(1) < yBegin || temp_R(1) > yEnd ||
            temp_R(2) < zBegin || temp_R(2) > zEnd) continue;

        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
        // out of bounds?
        if ((*fi)->getFieldstrength(temp_R, tmpE, tmpB)) continue;

        ++fcount;

        double frequency = (*rffi);   // frequency in [MHz]
        double ebscale = (*escali);   // E and B field scaling

        if (fieldType_m == BFieldType::SYNCHRO) {
            double powert = 1;
            for (const double fcoef : (*rffci)) {
                powert *= (t * 1e-9);
                frequency += fcoef * powert; // Add frequency ramp [MHz/s^n]
            }
            powert = 1;
            for (const double vcoef : (*rfvci)) {
                powert *= (t * 1e-9);
                ebscale += vcoef * powert; // Add frequency ramp [MHz/s^n]
            }
        }

        double phase = Physics::two_pi * 1.0E-3 * frequency * t + (*rfphii); // f in [MHz], t in [ns]

        E += ebscale * std::cos(phase) * tmpE;
        B -= ebscale * std::sin(phase) * tmpB;

        if (fieldType_m != BFieldType::SYNCHRO)
            continue;

        // Some phase output -DW
        double phase_print = phase * Physics::rad2deg;
        if (tet >= 90.0 && waitingGap_m == 1) {
            phase_print = std::fmod(phase_print, 360) - 360.0;

            *gmsg << endl << "Gap 1 phase = " << phase_print << " Deg" << endl;
            *gmsg << "Gap 1 E-Field = (" << E[0] << "/" << E[1] << "/" << E[2] << ")" << endl;
            *gmsg << "Gap 1 B-Field = (" << B[0] << "/" << B[1] << "/" << B[2] << ")" << endl;
            *gmsg << "RF Frequency = " << frequency << " MHz" << endl;

            waitingGap_m = 2;
        } else if (tet >= 270.0 && waitingGap_m == 2) {
            phase_print = std::fmod(phase_print, 360) - 360.0;

            *gmsg << endl << "Gap 2 phase = " << phase_print << " Deg" << endl;
            *gmsg << "Gap 2 E-Field = (" << E[0] << "/" << E[1] << "/" << E[2] << ")" << endl;
            *gmsg << "Gap 2 B-Field = (" << B[0] << "/" << B[1] << "/" << B[2] << ")" << endl;
            *gmsg << "RF Frequency = " << frequency << " MHz" << endl;
            waitingGap_m = 0;
        }
        if (fieldType_m == BFieldType::SYNCHRO) {
            ++rffci, ++rfvci;
        }
    }
    return false;
}

void Cyclotron::apply(const double& rad, const double& z,
                      const double& tet_rad, double& br,
                      double& bt, double& bz) {
    this->interpolate(rad, tet_rad, br, bt, bz);
    this->applyTrimCoil(rad, z, tet_rad, br, bz);
}

void Cyclotron::finalise() {
    online_m = false;
    lossDs_m->save();
    *gmsg << "* Finalize cyclotron" << endl;
}

bool Cyclotron::bends() const {
    return true;
}

// calculate derivatives with 5-point lagrange's formula.
double Cyclotron::gutdf5d(double* f, double dx, const int kor,
                          const int krl, const int lpr) {

    double C[5][5][3], FAC[3];
    double result;
    int j;
    /* CALCULATE DERIVATIVES WITH 5-POINT LAGRANGE FORMULA
     * PARAMETERS:
     * F  STARTADDRESS FOR THE 5 SUPPORT POINTS
     * DX STEPWIDTH FOR ARGUMENT
     * KOR        ORDER OF DERIVATIVE (KOR=1,2,3).
     * KRL        NUMBER OF SUPPORT POINT, WHERE THE DERIVATIVE IS TO BE CALCULATED
     *  (USUALLY 3, USE FOR BOUNDARY 1 ,2, RESP. 4, 5)
     * LPR        DISTANCE OF THE 5 STORAGE POSITIONS (=1 IF THEY ARE NEIGHBORS OR LENGTH
     * OF COLUMNLENGTH OF A MATRIX, IF THE SUPPORT POINTS ARE ON A LINE).
     * ATTENTION! THE INDICES ARE NOW IN C-FORMAT AND NOT IN FORTRAN-FORMAT.*/

    /* COEFFICIENTS FOR THE 1ST DERIVATIVE: */
    C[0][0][0] = -50.0;
    C[1][0][0] = 96.0;
    C[2][0][0] = -72.0;
    C[3][0][0] = 32.0;
    C[4][0][0] = -6.0;
    C[0][1][0] = -6.0;
    C[1][1][0] = -20.0;
    C[2][1][0] = 36.0;
    C[3][1][0] = -12.0;
    C[4][1][0] =  2.0;
    C[0][2][0] =  2.0;
    C[1][2][0] = -16.0;
    C[2][2][0] =  0.0;
    C[3][2][0] = 16.0;
    C[4][2][0] = -2.0;
    C[0][3][0] = -2.0;
    C[1][3][0] = 12.0;
    C[2][3][0] = -36.0;
    C[3][3][0] = 20.0;
    C[4][3][0] =  6.0;
    C[0][4][0] =  6.0;
    C[1][4][0] = -32.0;
    C[2][4][0] = 72.0;
    C[3][4][0] = -96.0;
    C[4][4][0] = 50.0;

    /* COEFFICIENTS FOR THE 2ND DERIVATIVE: */
    C[0][0][1] = 35.0;
    C[1][0][1] = -104;
    C[2][0][1] = 114.0;
    C[3][0][1] = -56.0;
    C[4][0][1] = 11.0;
    C[0][1][1] = 11.0;
    C[1][1][1] = -20.0;
    C[2][1][1] =  6.0;
    C[3][1][1] =  4.0;
    C[4][1][1] = -1.0;
    C[0][2][1] = -1.0;
    C[1][2][1] = 16.0;
    C[2][2][1] = -30.0;
    C[3][2][1] = 16.0;
    C[4][2][1] = -1.0;
    C[0][3][1] = -1.0;
    C[1][3][1] =  4.0;
    C[2][3][1] =  6.0;
    C[3][3][1] = -20.0;
    C[4][3][1] = 11.0;
    C[0][4][1] = 11.0;
    C[1][4][1] = -56.0;
    C[2][4][1] = 114.0;
    C[3][4][1] = -104;
    C[4][4][1] = 35.0;

    /* COEFFICIENTS FOR THE 3RD DERIVATIVE: */
    C[0][0][2] = -10.0;
    C[1][0][2] = 36.0;
    C[2][0][2] = -48.0;
    C[3][0][2] = 28.0;
    C[4][0][2] = -6.0;
    C[0][1][2] = -6.0;
    C[1][1][2] = 20.0;
    C[2][1][2] = -24.0;
    C[3][1][2] = 12.0;
    C[4][1][2] = -2.0;
    C[0][2][2] = -2.0;
    C[1][2][2] =  4.0;
    C[2][2][2] =  0.0;
    C[3][2][2] = -4.0;
    C[4][2][2] =  2.0;
    C[0][3][2] =  2.0;
    C[1][3][2] = -12.0;
    C[2][3][2] = 24.0;
    C[3][3][2] = -20.0;
    C[4][3][2] =  6.0;
    C[0][4][2] =  6.0;
    C[1][4][2] = -28.0;
    C[2][4][2] = 48.0;
    C[3][4][2] = -36.0;
    C[4][4][2] = 10.0;

    /* FACTOR: */
    FAC[0] = 24.0;
    FAC[1] = 12.0;
    FAC[2] = 4.0;

    result = 0.0;
    for (j = 0; j < 5; j++) {
        result += C[j][krl][kor] * *(f + j * lpr);
    }

    return result / (FAC[kor] * std::pow(dx, (kor + 1)));
}


bool Cyclotron::interpolate(const double& rad,
                            const double& tet_rad,
                            double& brint,
                            double& btint,
                            double& bzint) {

    const double xir = (rad - BP_m.rmin_m) / BP_m.delr_m;

    // ir : the number of path whose radius is less than the 4 points of cell which surround the particle.
    const int ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;

    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    double tet_map = std::fmod(tet_rad * Physics::rad2deg, 360.0 / symmetry_m);

    double xit = tet_map / BP_m.dtet_m;
    int it = (int) xit;

    const double wt1 = xit - (double)it;
    const double wt2 = 1.0 - wt1;

    // it : the number of point on the inner path whose angle is less than
    // the particle' corresponding angle.
    // include zero degree point
    it++;

    int r1t1, r2t1, r1t2, r2t2;
    int ntetS = Bfield_m.ntet_m + 1;

    // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
    // considering  the array start with index of zero, minus 1.

    if (fieldType_m != BFieldType::FFABF) {
        /*
          For FFA this does not work
        */
        r1t1 = it + ntetS * ir - 1;
        r1t2 = r1t1 + 1;
        r2t1 = r1t1 + ntetS;
        r2t2 = r2t1 + 1 ;

    } else {
        /*
          With this we have B-field AND this is far more
          intuitive for me ....
        */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
    }

    if ((it >= 0) && (ir >= 0) && (it < Bfield_m.ntetS_m) && (ir < Bfield_m.nrad_m)) {
        // B_{z}
        double bzf = Bfield_m.bfld_m[r1t1] * wr2 * wt2 +
                     Bfield_m.bfld_m[r2t1] * wr1 * wt2 +
                     Bfield_m.bfld_m[r1t2] * wr2 * wt1 +
                     Bfield_m.bfld_m[r2t2] * wr1 * wt1;
        bzint = /*- */bzf ;

        // dB_{z}/dr
        brint = Bfield_m.dbr_m[r1t1] * wr2 * wt2 +
                Bfield_m.dbr_m[r2t1] * wr1 * wt2 +
                Bfield_m.dbr_m[r1t2] * wr2 * wt1 +
                Bfield_m.dbr_m[r2t2] * wr1 * wt1;

        // dB_{z}/dtheta
        btint = Bfield_m.dbt_m[r1t1] * wr2 * wt2 +
                Bfield_m.dbt_m[r2t1] * wr1 * wt2 +
                Bfield_m.dbt_m[r1t2] * wr2 * wt1 +
                Bfield_m.dbt_m[r2t2] * wr1 * wt1;

        return true;
    }
    return false;
}


void Cyclotron::read(const double& scaleFactor) {

    setBFieldType();

    switch (fieldType_m) {
        case BFieldType::PSIBF: {
            *gmsg << "* Read field data from PSI format field map file" << endl;
            getFieldFromFile_Ring(scaleFactor);
            break;
        }
        case BFieldType::CARBONBF: {
            *gmsg << "* Read data from 450MeV Carbon cyclotron field file" << endl;
            getFieldFromFile_Carbon(scaleFactor);
            break;
        }
        case BFieldType::ANSYSBF: {
            *gmsg << "* Read data from 100MeV H- cyclotron CYCIAE-100 field file" << endl;
            getFieldFromFile_CYCIAE(scaleFactor);
            break;
        }
        case BFieldType::AVFEQBF: {
            *gmsg << "* Read AVFEQ data (Riken)" << endl;
            getFieldFromFile_AVFEQ(scaleFactor);
            break;
        }
        case BFieldType::FFABF: {
            *gmsg << "* Read FFA data MSU/FNAL" << endl;
            getFieldFromFile_FFA(scaleFactor);
            break;
        }
        case BFieldType::BANDRF: {
            *gmsg << "* Read both median plane B field map and 3D E field map of RF cavity for compact cyclotron" << endl;
            getFieldFromFile_BandRF(scaleFactor);
            break;
        }
        case BFieldType::SYNCHRO: {
            *gmsg << "* Read midplane B-field, 3D RF fieldmaps, and text files with RF frequency/Voltage coefficients for Synchrocyclotron" << endl;
            getFieldFromFile_Synchrocyclotron(scaleFactor);
            break;
        }
    }

    // calculate the radii of initial grid.
    initR(BP_m.rmin_m, BP_m.delr_m, Bfield_m.nrad_m);

    // calculate the remaining derivatives
    getdiffs();
}

// evaluate other derivative of magnetic field.
void Cyclotron::getdiffs() {

    Bfield_m.dbr_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbrr_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbrrr_m.resize(Bfield_m.ntot_m);

    Bfield_m.dbrt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbrrt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbrtt_m.resize(Bfield_m.ntot_m);

    Bfield_m.f2_m.resize(Bfield_m.ntot_m);
    Bfield_m.f3_m.resize(Bfield_m.ntot_m);
    Bfield_m.g3_m.resize(Bfield_m.ntot_m);

    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        for (int k = 0; k < Bfield_m.ntet_m; k++) {

            double dtheta = Physics::deg2rad * BP_m.dtet_m;

            int kEdge;

            kEdge = std::max(k - 2, 0);
            kEdge = std::min(kEdge, Bfield_m.ntet_m - 5);

            int dkFromEdge = k - kEdge;
            int index = idx(i, k);
            int indexkEdge = idx(i, kEdge);

            Bfield_m.dbt_m[index]   = gutdf5d(&Bfield_m.bfld_m[indexkEdge], dtheta, 0, dkFromEdge, 1);
            Bfield_m.dbtt_m[index]  = gutdf5d(&Bfield_m.bfld_m[indexkEdge], dtheta, 1, dkFromEdge, 1);
            Bfield_m.dbttt_m[index] = gutdf5d(&Bfield_m.bfld_m[indexkEdge], dtheta, 2, dkFromEdge, 1);
        }
    }

    for (int k = 0; k < Bfield_m.ntet_m; k++) {
        // inner loop varies R
        for (int i = 0; i < Bfield_m.nrad_m; i++) {
            double rac = BP_m.rarr_m[i];
            // define iredg, the reference index for radial interpolation
            // standard: i-2 minimal: 0 (not negative!)  maximal: nrad-4
            int iredg = std::max(i - 2, 0);
            iredg = std::min(iredg, Bfield_m.nrad_m - 5);
            int irtak = i - iredg;
            int index = idx(i, k);
            int indexredg = idx(iredg, k);

            Bfield_m.dbr_m[index]   = gutdf5d(&Bfield_m.bfld_m[indexredg], BP_m.delr_m, 0, irtak, Bfield_m.ntetS_m);
            Bfield_m.dbrr_m[index]  = gutdf5d(&Bfield_m.bfld_m[indexredg], BP_m.delr_m, 1, irtak, Bfield_m.ntetS_m);
            Bfield_m.dbrrr_m[index] = gutdf5d(&Bfield_m.bfld_m[indexredg], BP_m.delr_m, 2, irtak, Bfield_m.ntetS_m);

            Bfield_m.dbrt_m[index]  = gutdf5d(&Bfield_m.dbt_m[indexredg],  BP_m.delr_m, 0, irtak, Bfield_m.ntetS_m);
            Bfield_m.dbrrt_m[index] = gutdf5d(&Bfield_m.dbt_m[indexredg],  BP_m.delr_m, 1, irtak, Bfield_m.ntetS_m);
            Bfield_m.dbrtt_m[index] = gutdf5d(&Bfield_m.dbtt_m[indexredg], BP_m.delr_m, 0, irtak, Bfield_m.ntetS_m);

            // fehlt noch!! f2,f3,g3,
            Bfield_m.f2_m[index] = (Bfield_m.dbrr_m[index]
                                    + Bfield_m.dbr_m[index] / rac
                                    + Bfield_m.dbtt_m[index] / rac / rac) / 2.0;

            Bfield_m.f3_m[index] = (Bfield_m.dbrrr_m[index]
                                    + Bfield_m.dbrr_m[index] / rac
                                    + (Bfield_m.dbrtt_m[index] - Bfield_m.dbr_m[index]) / rac / rac
                                    - 2.0 * Bfield_m.dbtt_m[index] / rac / rac / rac) / 6.0;

            Bfield_m.g3_m[index] = (Bfield_m.dbrrt_m[index]
                                    + Bfield_m.dbrt_m[index] / rac
                                    + Bfield_m.dbttt_m[index] / rac / rac) / 6.0;
        } // Radius Loop
    } // Azimuth loop

    // copy 1st azimuth to last + 1 to always yield an interval
    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        int iend = idx(i, Bfield_m.ntet_m);
        int istart = idx(i, 0);

        Bfield_m.bfld_m[iend]  = Bfield_m.bfld_m[istart];
        Bfield_m.dbt_m[iend]   = Bfield_m.dbt_m[istart];
        Bfield_m.dbtt_m[iend]  = Bfield_m.dbtt_m[istart];
        Bfield_m.dbttt_m[iend] = Bfield_m.dbttt_m[istart];

        Bfield_m.dbr_m[iend]   = Bfield_m.dbr_m[istart];
        Bfield_m.dbrr_m[iend]  = Bfield_m.dbrr_m[istart];
        Bfield_m.dbrrr_m[iend] = Bfield_m.dbrrr_m[istart];

        Bfield_m.dbrt_m[iend]  = Bfield_m.dbrt_m[istart];
        Bfield_m.dbrtt_m[iend] = Bfield_m.dbrtt_m[istart];
        Bfield_m.dbrrt_m[iend] = Bfield_m.dbrrt_m[istart];

        Bfield_m.f2_m[iend]    = Bfield_m.f2_m[istart];
        Bfield_m.f3_m[iend]    = Bfield_m.f3_m[istart];
        Bfield_m.g3_m[iend]    = Bfield_m.g3_m[istart];
    }

    /* debug

    for (int i = 0; i< Bfield_m.nrad_m; i++){
      for (int j = 0; j< Bfield_m.ntetS_m; j++){
    int index = idx(i,j);
    double x = (BP_m.rmin_m+i*BP_m.delr_m) * std::sin(j*BP_m.dtet_m*pi/180.0);
    double y = (BP_m.rmin_m+i*BP_m.delr_m) * std::cos(j*BP_m.dtet_m*pi/180.0);
    *gmsg<<"x= "<<x<<" y= "<<y<<" B= "<<Bfield_m.bfld_m[index]<<endl;
      }
    }
    */
}


// Calculates Radii of initial grid.
// dimensions in [m]!
void Cyclotron::initR(double rmin, double dr, int nrad) {
    BP_m.rarr_m.resize(nrad);
    for (int i = 0; i < nrad; i++) {
        BP_m.rarr_m[i] = rmin + i * dr;
    }
    BP_m.delr_m = dr;
}

void Cyclotron::initialise(PartBunchBase<double, 3>* bunch, double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
    online_m = true;
}

void Cyclotron::initialise(PartBunchBase<double, 3>* bunch, const double& scaleFactor) {
    RefPartBunch_m = bunch;
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getOutputFN(), !Options::asciidump && Options::enableHDF5));

    this->read(scaleFactor);
}


// Read field map from external file.
void Cyclotron::getFieldFromFile_Ring(const double& scaleFactor) {

    FILE *f = NULL;
    int lpar;
    char fout[100];
    double dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*             READ IN RING FIELD MAP            " << endl;
    *gmsg << "*      (The first data block is useless)        " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP_m.Bfact_m = scaleFactor;

    f = std::fopen(fmapfn_m.c_str(), "r");

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.rmin_m));
    *gmsg << "* Minimal radius of measured field map: " << BP_m.rmin_m << " [mm]" << endl;
    BP_m.rmin_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.delr_m));
    //if the value is negative, the actual value is its reciprocal.
    if (BP_m.delr_m < 0.0) BP_m.delr_m = 1.0 / (-BP_m.delr_m);
    *gmsg << "* Stepsize in radial direction: " << BP_m.delr_m << " [mm]" << endl;
    BP_m.delr_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.tetmin_m));
    *gmsg << "* Minimal angle of measured field map: " << BP_m.tetmin_m << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.dtet_m));
    //if the value is negative, the actual value is its reciprocal.
    if (BP_m.dtet_m < 0.0) BP_m.dtet_m = 1.0 / (-BP_m.dtet_m);
    *gmsg << "* Stepsize in azimuth direction: " << BP_m.dtet_m << " [deg]" << endl;

    for (int i = 0; i < 13; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
    }

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.nrad_m));
    *gmsg << "* Index in radial direction: " << Bfield_m.nrad_m << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.ntet_m));
    *gmsg << "* Index in azimuthal direction: " << Bfield_m.ntet_m << endl;

    Bfield_m.ntetS_m = Bfield_m.ntet_m + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield_m.ntetS_m << endl;

    for (int i = 0; i < 5; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
    }
    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &lpar));
    // msg<< "READ"<<lpar<<" DATA ENTRIES"<<endl;

    for (int i = 0; i < 4; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
    }

    for (int i = 0; i < lpar; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &dtmp));
    }
    for (int i = 0; i < 6; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
    }
    //*gmsg << "* READ FILE DESCRIPTION..." <<endl;
    for (int i = 0; i < 10000; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
        if (std::strcmp(fout, "LREC=") == 0)break;
    }

    for (int i = 0; i < 5; i++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
    }
    Bfield_m.ntot_m = idx(Bfield_m.nrad_m - 1, Bfield_m.ntet_m) + 1;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ): " << Bfield_m.ntot_m << endl;

    Bfield_m.bfld_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbtt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbttt_m.resize(Bfield_m.ntot_m);

    *gmsg << "* Read-in loop one block per radius" << endl;
    *gmsg << "* Rescaling of the magnetic fields with factor: " << BP_m.Bfact_m << endl;
    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        if (i > 0) {
            for (int dummy = 0; dummy < 6; dummy++) {
                CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout)); // INFO-LINE
            }
        }
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.bfld_m[idx(i, k)])));
            Bfield_m.bfld_m[idx(i, k)] *= BP_m.Bfact_m;
        }
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.dbt_m[idx(i, k)])));
            Bfield_m.dbt_m[idx(i, k)] *= BP_m.Bfact_m;
        }
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.dbtt_m[idx(i, k)])));
            Bfield_m.dbtt_m[idx(i, k)] *= BP_m.Bfact_m;
        }
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.dbttt_m[idx(i, k)])));
            Bfield_m.dbttt_m[idx(i, k)] *= BP_m.Bfact_m;
        }
    }
    std::fclose(f);

    *gmsg << "* Field Map read successfully!" << endl << endl;
}


void Cyclotron::getFieldFromFile_FFA(const double& /*scaleFactor*/) {

    /*
      Field is read in from ascii file (COSY output) in the order:
      R(m) theta(Deg) x(m) y(m) Bz(T).

      Theta is the fast varying variable

      2.0000   0.0  2.0000  0.0000      0.0000000000000000
      2.0000   1.0  1.9997  0.0349      0.0000000000000000
      2.0000   2.0  1.9988  0.0698      0.0000000000000000
      2.0000   3.0  1.9973  0.1047      0.0000000000000000

      ......
      <blank line>

      2.1000   0.0  2.1000  0.0000      0.0000000000000000
      2.1000   1.0  2.0997  0.0367      0.0000000000000000
      2.1000   2.0  2.0987  0.0733      0.0000000000000000
    */

    std::vector<double> rv;
    std::vector<double> thv;
    std::vector<double> xv;
    std::vector<double> yv;
    std::vector<double> bzv;
    std::vector<double>::iterator vit;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*             READ IN FFA FIELD MAP             " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP_m.Bfact_m = -10.0; // T->kG and H- for the current FNAL FFA

    std::ifstream file_to_read(fmapfn_m.c_str());
    const int max_num_of_char_in_a_line = 128;
    const int num_of_header_lines = 1;

    // STEP2: SKIP ALL THE HEADER LINES
    for (int i = 0; i < num_of_header_lines; ++i)
        file_to_read.ignore(max_num_of_char_in_a_line, '\n');

    // TEMP for OPAL 2.0 changing this to m -DW
    while(!file_to_read.eof()) {
        double r, th, x, y, bz;
        file_to_read >> r >> th >> x >> y >> bz;
        if ((int)th != 360) {
            rv.push_back(r);
            thv.push_back(th);
            xv.push_back(x);
            yv.push_back(y);
            bzv.push_back(bz);
        }
    }

    double maxtheta = 360.0;
    BP_m.dtet_m = thv[1] - thv[0];
    BP_m.rmin_m = *(rv.begin());
    double rmax = rv.back();

    // find out dR
    for (vit = rv.begin(); *vit <= BP_m.rmin_m; ++vit) {}
    BP_m.delr_m = *vit - BP_m.rmin_m;

    BP_m.tetmin_m = thv[0];

    Bfield_m.ntet_m = (int)((maxtheta - thv[0]) / BP_m.dtet_m);
    Bfield_m.nrad_m = (int)(rmax - BP_m.rmin_m) / BP_m.delr_m + 1;
    Bfield_m.ntetS_m = Bfield_m.ntet_m + 1;
    *gmsg << "* Minimal radius of measured field map: " << 1000.0 * BP_m.rmin_m << " [mm]" << endl;
    *gmsg << "* Maximal radius of measured field map: " << 1000.0 * rmax << " [mm]" << endl;
    *gmsg << "* Stepsize in radial direction: " << 1000.0 * BP_m.delr_m << " [mm]" << endl;
    *gmsg << "* Minimal angle of measured field map: " << BP_m.tetmin_m << " [deg]" << endl;
    *gmsg << "* Maximal angle of measured field map: " << maxtheta << " [deg]" << endl;

    //if the value is negtive, the actual value is its reciprocal.
    if (BP_m.dtet_m < 0.0) BP_m.dtet_m = 1.0 / (-BP_m.dtet_m);
    *gmsg << "* Stepsize in azimuth direction: " << BP_m.dtet_m << " [deg]" << endl;
    *gmsg << "* Total grid point along azimuth:  " << Bfield_m.ntetS_m << endl;
    *gmsg << "* Total grid point along radius: " << Bfield_m.nrad_m << endl;

    Bfield_m.ntot_m = Bfield_m.ntetS_m * Bfield_m.nrad_m;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ): " << Bfield_m.ntot_m << endl;

    Bfield_m.bfld_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbtt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbttt_m.resize(Bfield_m.ntot_m);

    *gmsg << "* Rescaling of the magnetic fields with factor: " << BP_m.Bfact_m << endl;

    int count = 0;
    for (int r = 0; r < Bfield_m.nrad_m; r++) {
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            Bfield_m.bfld_m[idx(r, k)] = bzv[count] * BP_m.Bfact_m;
            count++;
        }
    }

    if ((Ippl::getNodes()) == 1 && Options::info) {
        writeOutputFieldFiles();
    }
    *gmsg << "* Number of elements read: " << count << endl;
    *gmsg << "* Field Map read successfully!" << endl << endl;
}


void Cyclotron::getFieldFromFile_AVFEQ(const double& scaleFactor) {

    FILE *f = NULL;
    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*        READ IN AVFEQ CYCLOTRON FIELD MAP      " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    /*  From Hiroki-san
        The first line tells r minimum (500mm),
                             r maximum(4150mm),
                             r step(50mm),
                             theta minimum(0deg),
                             theta maximum(90deg)
                             theta step(0.5deg).

        From the next line data repeat the block for a given r which the
        first line of the block tells. Each block consists of the data Bz
        from theta minimum (0deg) to theta maximum(90deg) with theta
        step(0.5deg).
    */

    BP_m.Bfact_m = scaleFactor / 1000.;

    f = std::fopen(fmapfn_m.c_str(), "r");

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.rmin_m));
    *gmsg << "* Minimal radius of measured field map: " << BP_m.rmin_m << " [mm]" << endl;
    BP_m.rmin_m *= 0.001;  // mm --> m

    double rmax;
    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &rmax));
    *gmsg << "* Maximal radius of measured field map: " << rmax << " [mm]" << endl;
    rmax *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.delr_m));
    *gmsg << "* Stepsize in radial direction: " << BP_m.delr_m << " [mm]" << endl;
    BP_m.delr_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.tetmin_m));
    *gmsg << "* Minimal angle of measured field map: " << BP_m.tetmin_m << " [deg]" << endl;

    double tetmax;
    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &tetmax));
    *gmsg << "* Maximal angle of measured field map: " << tetmax << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.dtet_m));
    //if the value is nagtive, the actual value is its reciprocal.

    if (BP_m.dtet_m < 0.0) BP_m.dtet_m = 1.0 / (-BP_m.dtet_m);
    *gmsg << "* Stepsize in azimuth direction: " << BP_m.dtet_m << " [deg]" << endl;

    Bfield_m.ntetS_m = (int)((tetmax - BP_m.tetmin_m) / BP_m.dtet_m + 1);
    *gmsg << "* Total grid point along azimuth:  " << Bfield_m.ntetS_m << endl;

    Bfield_m.nrad_m = (int)(rmax - BP_m.rmin_m) / BP_m.delr_m;

    int ntotidx = idx(Bfield_m.nrad_m, Bfield_m.ntetS_m) + 1;

    Bfield_m.ntot_m = Bfield_m.ntetS_m * Bfield_m.nrad_m;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ): "
          << Bfield_m.ntot_m << " ntot-idx= " << ntotidx << endl;

    Bfield_m.bfld_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbtt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbttt_m.resize(Bfield_m.ntot_m);

    *gmsg << "* Rescaling of the magnetic fields with factor: " << BP_m.Bfact_m << endl;

    double tmp;
    int count = 0;
    for (int r = 0; r < Bfield_m.nrad_m; r++) {
        CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &tmp));   // over read
        for (int k = 0; k < Bfield_m.ntetS_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.bfld_m[idx(r, k)])));
            Bfield_m.bfld_m[idx(r, k)] *= BP_m.Bfact_m;
            count++;
        }
    }

    if ((Ippl::getNodes()) == 1 && Options::info) {
        writeOutputFieldFiles();
    }

    std::fclose(f);
    *gmsg << "* Number of elements read: " << count << endl;
    *gmsg << "* Field Map read successfully!" << endl << endl;
}


void Cyclotron::getFieldFromFile_Carbon(const double& scaleFactor) {

    FILE *f = NULL;
    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*      READ IN CARBON CYCLOTRON FIELD MAP       " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP_m.Bfact_m = scaleFactor;

    f = std::fopen(fmapfn_m.c_str(), "r");

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.rmin_m));
    *gmsg << "* Minimal radius of measured field map: " << BP_m.rmin_m << " [mm]" << endl;
    BP_m.rmin_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.delr_m));
    //if the value is negative, the actual value is its reciprocal.
    if (BP_m.delr_m < 0.0) BP_m.delr_m = 1.0 / (-BP_m.delr_m);
    *gmsg << "* Stepsize in radial direction: " << BP_m.delr_m << " [mm]" << endl;
    BP_m.delr_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.tetmin_m));
    *gmsg << "* Minimal angle of measured field map: " << BP_m.tetmin_m << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.dtet_m));
    //if the value is negative, the actual value is its reciprocal.
    if (BP_m.dtet_m < 0.0) BP_m.dtet_m = 1.0 / (-BP_m.dtet_m);
    *gmsg << "* Stepsize in azimuthal direction: " << BP_m.dtet_m << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.ntet_m));
    *gmsg << "* Grid points along azimuth (ntet): " << Bfield_m.ntet_m << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.nrad_m));
    *gmsg << "* Grid points along radius (nrad): " << Bfield_m.nrad_m << endl;

    //Bfield_m.ntetS = Bfield_m.ntet;
    Bfield_m.ntetS_m = Bfield_m.ntet_m + 1;
    //*gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield_m.ntetS << endl;

    //Bfield_m.ntot = idx(Bfield_m.nrad - 1, Bfield_m.ntet) + 1;
    Bfield_m.ntot_m = Bfield_m.nrad_m * Bfield_m.ntetS_m;

    *gmsg << "* Adding a guard cell along azimuth" << endl;
    *gmsg << "* Total stored grid point number ((ntet+1) * nrad): " << Bfield_m.ntot_m << endl;
    Bfield_m.bfld_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbtt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbttt_m.resize(Bfield_m.ntot_m);

    *gmsg << "* Rescaling of the magnetic fields with factor: " << BP_m.Bfact_m << endl;

    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%16lE", &(Bfield_m.bfld_m[idx(i, k)])));
            Bfield_m.bfld_m[idx(i, k)] *= BP_m.Bfact_m;
        }
    }

    if ((Ippl::getNodes()) == 1 && Options::info) {
        writeOutputFieldFiles();
    }

    std::fclose(f);
    *gmsg << "* Field Map read successfully!" << endl << endl;
}


void Cyclotron::getFieldFromFile_CYCIAE(const double& scaleFactor) {

    FILE *f = NULL;
    char fout[100];
    int dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*    READ IN CYCIAE-100 CYCLOTRON FIELD MAP     " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP_m.Bfact_m = scaleFactor;

    f = std::fopen(fmapfn_m.c_str(), "r");

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.rmin_m));
    *gmsg << "* Minimal radius of measured field map: " << BP_m.rmin_m << " [mm]" << endl;
    BP_m.rmin_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.delr_m));
    *gmsg << "* Stepsize in radial direction: " << BP_m.delr_m << " [mm]" << endl;
    BP_m.delr_m *= 0.001;  // mm --> m

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.tetmin_m));
    *gmsg << "* Minimal angle of measured field map: " << BP_m.tetmin_m << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &BP_m.dtet_m));
    //if the value is nagtive, the actual value is its reciprocal.
    if (BP_m.dtet_m < 0.0) BP_m.dtet_m = 1.0 / (-BP_m.dtet_m);
    *gmsg << "* Stepsize in azimuth direction: " << BP_m.dtet_m << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.ntet_m));
    *gmsg << "* Index in azimuthal direction: " << Bfield_m.ntet_m << endl;

    CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &Bfield_m.nrad_m));
    *gmsg << "* Index in radial direction: " << Bfield_m.nrad_m << endl;

    Bfield_m.ntetS_m = Bfield_m.ntet_m + 1;
    *gmsg << "* Accordingly, total grid point along azimuth: " << Bfield_m.ntetS_m << endl;

    Bfield_m.ntot_m = idx(Bfield_m.nrad_m - 1, Bfield_m.ntet_m) + 1;

    *gmsg << "* Total stored grid point number ( ntetS * nrad ): " << Bfield_m.ntot_m << endl;
    Bfield_m.bfld_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbtt_m.resize(Bfield_m.ntot_m);
    Bfield_m.dbttt_m.resize(Bfield_m.ntot_m);

    *gmsg << "* Rescaling of the magnetic fields with factor: " << BP_m.Bfact_m << endl;

    int nHalfPoints = Bfield_m.ntet_m / 2.0 + 1;

    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        for (int ii = 0; ii < 13; ii++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%s", fout));
        }
        for (int k = 0; k < nHalfPoints; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(std::fscanf(f, "%lf", &(Bfield_m.bfld_m[idx(i, k)])));
            //  T --> kGs, minus for minus hydrogen
            Bfield_m.bfld_m[idx(i, k)] = Bfield_m.bfld_m[idx(i, k)] * (-10.0);
        }
        for (int k = nHalfPoints; k < Bfield_m.ntet_m; k++) {
            Bfield_m.bfld_m[idx(i, k)] = Bfield_m.bfld_m[idx(i, Bfield_m.ntet_m-k)];
        }
    }

    std::fclose(f);
    *gmsg << "* Field Map read successfully!" << endl << endl;
}


void Cyclotron::getFieldFromFile_BandRF(const double& scaleFactor) {

    // read 3D E&B field data file
    // loop over all field maps and superpose fields
    for (auto& fm: RFfilename_m) {
        Fieldmap *f = Fieldmap::getFieldmap(fm, false);
        *gmsg << "* Reading " << fm << endl;
        f->readMap();
        RFfields_m.push_back(f);
    }
    // read CARBON type B field
    getFieldFromFile_Carbon(scaleFactor);
}


void Cyclotron::getFieldFromFile_Synchrocyclotron(const double& scaleFactor) {

    // read 3D E&B field data file
    std::vector<std::string>::const_iterator fm = RFfilename_m.begin();
    std::vector<std::string>::const_iterator rffcfni = RFFCoeff_fn_m.begin();
    std::vector<std::string>::const_iterator rfvcfni = RFVCoeff_fn_m.begin();
    // loop over all field maps and superpose fields
    int fcount = 0;
    FILE *rffcf = NULL;
    FILE *rfvcf = NULL;

    *gmsg << endl;
    *gmsg << "* ------------------------------------------------------------" << endl;
    *gmsg << "*      READ IN 3D RF Fields and Frequency Coefficients        " << endl;
    *gmsg << "* ------------------------------------------------------------" << endl;

    for (; fm != RFfilename_m.end(); ++fm, ++rffcfni, ++rfvcfni, ++fcount) {
        Fieldmap *f = Fieldmap::getFieldmap(*fm, false);
        f->readMap();
        // if (IPPL::Comm->getOutputLevel() != 0) f->getInfo(gmsg);
        RFfields_m.push_back(f);

        // Read RF Frequency Coefficients from file
        *gmsg << "RF Frequency Coefficient Filename: " << (*rffcfni) << endl;

        rffcf = std::fopen((*rffcfni).c_str(), "r");

        if (rffcf == NULL) {
            throw GeneralClassicException(
                "Cyclotron::getFieldFromFile_Synchrocyclotron",
                "failed to open file '" + *rffcfni + "', please check if it exists");
        }

        std::vector<double> fcoeff;
        int nc; //Number of coefficients
        double value;

        CHECK_CYC_FSCANF_EOF(std::fscanf(rffcf, "%d", &nc));
        *gmsg << "* Number of coefficients in file: " << nc << endl;
        for (int k = 0; k < nc; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(rffcf, "%16lE", &value));
            fcoeff.push_back(value);
            //*gmsg << "* Coefficient " << k << ": " << value << endl;
        }
        rffc_m.push_back(fcoeff);

        std::fclose(rffcf);

        // Read RF Voltage Coefficients from file
        *gmsg << "RF Voltage Coefficient Filename: " << (*rfvcfni) << endl;

        rfvcf = std::fopen((*rfvcfni).c_str(), "r");
        if (rfvcf == NULL) {
            throw GeneralClassicException(
                "Cyclotron::getFieldFromFile_Synchrocyclotron",
                "failed to open file '" + *rfvcfni + "', please check if it exists");
        }

        std::vector<double> vcoeff;

        CHECK_CYC_FSCANF_EOF(std::fscanf(rfvcf, "%d", &nc));
        *gmsg << "* Number of coefficients in file: " << nc << endl;
        for (int k = 0; k < nc; k++) {
            CHECK_CYC_FSCANF_EOF(std::fscanf(rfvcf, "%16lE", &value));
            vcoeff.push_back(value);
            //*gmsg << "* Coefficient " << k << ": " << value << endl;
        }
        rfvc_m.push_back(vcoeff);

        std::fclose(rfvcf);
    }

    // read CARBON type B field for mid-plane field
    getFieldFromFile_Carbon(scaleFactor);
}

void Cyclotron::getDimensions(double& /*zBegin*/, double& /*zEnd*/) const
{ }


void Cyclotron::writeOutputFieldFiles() {
    std::fstream fp1;
    std::string fname = Util::combineFilePath({
        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
        "gnu.out"
    });
    fp1.open(fname, std::ios::out);
    for (int i = 0; i < Bfield_m.nrad_m; i++) {
        for (int k = 0; k < Bfield_m.ntet_m; k++) {
            fp1 << BP_m.rmin_m + (i * BP_m.delr_m) << " \t "
                << k * (BP_m.tetmin_m + BP_m.dtet_m) << " \t "
                << Bfield_m.bfld_m[idx(i, k)] << std::endl;
        }
    }
    fp1.close();

    if (fieldType_m == BFieldType::BANDRF  ||
        fieldType_m == BFieldType::SYNCHRO ||
        fieldType_m == BFieldType::CARBONBF) {
        std::fstream fp2;
        fname = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "eb.out"
        });
        fp2.open(fname, std::ios::out);
        for (int i = 0; i < Bfield_m.nrad_m; i++) {
            for (int k = 0; k < Bfield_m.ntet_m; k++) {
                Vector_t tmpR = Vector_t (BP_m.rmin_m + (i * BP_m.delr_m), 0.0, k * (BP_m.tetmin_m + BP_m.dtet_m));
                Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
                for (auto& fi: RFfields_m) {
                    Vector_t E(0.0, 0.0, 0.0), B(0.0, 0.0, 0.0);
                    if (!fi->getFieldstrength(tmpR, tmpE, tmpB)) {
                        tmpE += E;
                        tmpB -= B;
                    }
                }
                fp2 << tmpR << " \t E= " << tmpE << "\t B= " << tmpB << std::endl;
            }
        }
        *gmsg << "\n* Writing 'gnu.out' and 'eb.out' files of cyclotron field maps\n" << endl;
        fp2.close();
    } else {
        *gmsg << "\n* Writing 'gnu.out' file of cyclotron field map\n" << endl;
    }
}

#undef CHECK_CYC_FSCANF_EOF