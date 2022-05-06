//
// Class TravelingWave
//   Defines the abstract interface for Traveling Wave.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#include "AbsBeamline/TravelingWave.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Units.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

TravelingWave::TravelingWave():
    TravelingWave("")
{}


TravelingWave::TravelingWave(const TravelingWave& right):
    RFCavity(right),
    scaleCore_m(right.scaleCore_m),
    scaleCoreError_m(right.scaleCoreError_m),
    phaseCore1_m(right.phaseCore1_m),
    phaseCore2_m(right.phaseCore2_m),
    phaseExit_m(right.phaseExit_m),
    startCoreField_m(right.startCoreField_m),
    startExitField_m(right.startExitField_m),
    mappedStartExitField_m(right.mappedStartExitField_m),
    periodLength_m(right.periodLength_m),
    numCells_m(right.numCells_m),
    cellLength_m(right.cellLength_m),
    mode_m(right.mode_m)
{}


TravelingWave::TravelingWave(const std::string& name):
    RFCavity(name),
    scaleCore_m(1.0),
    scaleCoreError_m(0.0),
    phaseCore1_m(0.0),
    phaseCore2_m(0.0),
    phaseExit_m(0.0),
    startCoreField_m(0.0),
    startExitField_m(0.0),
    mappedStartExitField_m(0.0),
    periodLength_m(0.0),
    numCells_m(1),
    cellLength_m(0.0),
    mode_m(1)
{}


TravelingWave::~TravelingWave() {
}


void TravelingWave::accept(BeamlineVisitor& visitor) const {
    visitor.visitTravelingWave(*this);
}

bool TravelingWave::apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool TravelingWave::apply(const Vector_t& R,
                          const Vector_t& /*P*/,
                          const double& t,
                          Vector_t& E, Vector_t& B) {

    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength()) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * periodLength_m);
    double tmpcos, tmpsin;
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR)) return getFlagDeleteOnTransverseExit();

        tmpcos =  (scale_m + scaleError_m) * std::cos(frequency_m * t + phase_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(frequency_m * t + phase_m + phaseError_m);

    } else if (tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;
        if (!fieldmap_m->isInside(tmpR)) return getFlagDeleteOnTransverseExit();

        tmpcos =  (scaleCore_m + scaleCoreError_m) * std::cos(frequency_m * t + phaseCore1_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m) * std::sin(frequency_m * t + phaseCore1_m + phaseError_m);
        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  (scaleCore_m + scaleCoreError_m) * std::cos(frequency_m * t + phaseCore2_m + phaseError_m);
        tmpsin = -(scaleCore_m + scaleCoreError_m) * std::sin(frequency_m * t + phaseCore2_m + phaseError_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR)) return getFlagDeleteOnTransverseExit();
        tmpcos =  (scale_m + scaleError_m) * std::cos(frequency_m * t + phaseExit_m + phaseError_m);
        tmpsin = -(scale_m + scaleError_m) * std::sin(frequency_m * t + phaseExit_m + phaseError_m);
    }

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);

    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

bool TravelingWave::applyToReferenceParticle(const Vector_t& R,
                                             const Vector_t& /*P*/,
                                             const double& t, Vector_t& E,
                                             Vector_t& B) {

    if (R(2) < -0.5 * periodLength_m || R(2) + 0.5 * periodLength_m >= getElementLength()) return false;

    Vector_t tmpR = Vector_t(R(0), R(1), R(2) + 0.5 * periodLength_m);
    double tmpcos, tmpsin;
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if (tmpR(2) < startCoreField_m) {
        if (!fieldmap_m->isInside(tmpR)) return true;
        tmpcos =  scale_m * std::cos(frequency_m * t + phase_m);
        tmpsin = -scale_m * std::sin(frequency_m * t + phase_m);

    } else if (tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;
        if (!fieldmap_m->isInside(tmpR)) return true;

        tmpcos =  scaleCore_m * std::cos(frequency_m * t + phaseCore1_m);
        tmpsin = -scaleCore_m * std::sin(frequency_m * t + phaseCore1_m);
        fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + cellLength_m;
        tmpR(2) = tmpR(2) - periodLength_m * std::floor(tmpR(2) / periodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * std::cos(frequency_m * t + phaseCore2_m);
        tmpsin = -scaleCore_m * std::sin(frequency_m * t + phaseCore2_m);

    } else {
        tmpR(2) -= mappedStartExitField_m;
        if (!fieldmap_m->isInside(tmpR)) return true;

        tmpcos =  scale_m * std::cos(frequency_m * t + phaseExit_m);
        tmpsin = -scale_m * std::sin(frequency_m * t + phaseExit_m);
    }

    fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return false;
}

void TravelingWave::initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) {

    if (bunch == nullptr) {
        startField = - 0.5 * periodLength_m;
        endField = startExitField_m;
        return;
    }

    Inform msg("TravelingWave ", *gmsg);

    RefPartBunch_m = bunch;
    double zBegin = 0.0, zEnd = 0.0;
    RFCavity::initialise(bunch, zBegin, zEnd);
    if (std::abs(startField_m) > 0.0) {
        throw GeneralClassicException("TravelingWave::initialise",
                                      "The field map of a traveling wave structure has to begin at 0.0");
    }

    periodLength_m = (zEnd - zBegin) / 2.0;
    cellLength_m = periodLength_m * mode_m;
    startField_m = -0.5 * periodLength_m;

    startCoreField_m = periodLength_m / 2.0;
    startExitField_m = startCoreField_m + (numCells_m - 1) * cellLength_m;
    mappedStartExitField_m = startExitField_m - 3.0 * periodLength_m / 2.0;

    startField = -periodLength_m / 2.0;
    endField = startField + startExitField_m + periodLength_m / 2.0;
    setElementLength(endField - startField);

    scaleCore_m = scale_m / std::sin(Physics::two_pi * mode_m);
    scaleCoreError_m = scaleError_m / std::sin(Physics::two_pi * mode_m);
    phaseCore1_m = phase_m + Physics::pi * mode_m / 2.0;
    phaseCore2_m = phase_m + Physics::pi * mode_m * 1.5;
    phaseExit_m = phase_m - Physics::two_pi * ((numCells_m - 1) * mode_m - std::floor((numCells_m - 1) * mode_m));
}

void TravelingWave::finalise()
{}

bool TravelingWave::bends() const {
    return false;
}

void TravelingWave::goOnline(const double&) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void TravelingWave::goOffline() {
    Fieldmap::freeMap(filename_m);
}

void TravelingWave::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = -0.5 * periodLength_m;
    zEnd = zBegin + getElementLength();
}


void TravelingWave::getElementDimensions(double& begin, double& end) const {
    begin = -0.5 * periodLength_m;
    end = begin + getElementLength();
}

ElementType TravelingWave::getType() const {
    return ElementType::TRAVELINGWAVE;
}

double TravelingWave::getAutoPhaseEstimate(const double& E0,
                                           const double& t0,
                                           const double& q,
                                           const double& mass) {
    std::vector<double> t, E, t2, E2;
    std::vector<std::pair<double, double> > F;
    double phi = 0.0, tmp_phi, dphi = 0.5 * Units::deg2rad;
    double phaseC1 = phaseCore1_m - phase_m;
    double phaseC2 = phaseCore2_m - phase_m;
    double phaseE = phaseExit_m - phase_m;

    fieldmap_m->getOnaxisEz(F);
    if (F.size() == 0) return 0.0;

    int N1 = static_cast<int>(std::floor(F.size() / 4.)) + 1;
    int N2 = F.size() - 2 * N1 + 1;
    int N3 = 2 * N1 + static_cast<int>(std::floor((numCells_m - 1) * N2 * mode_m)) - 1;
    int N4 = static_cast<int>(std::round(N2 * mode_m));
    double Dz = F[N1 + N2].first - F[N1].first;

    t.resize(N3, t0);
    t2.resize(N3, t0);
    E.resize(N3, E0);
    E2.resize(N3, E0);
    for (int i = 1; i < N1; ++ i) {
        E[i] = E0 + (F[i].first + F[i - 1].first) / 2. * scale_m / mass;
        E2[i] = E[i];
    }
    for (int i = N1; i < N3 - N1 + 1; ++ i) {
        int I = (i - N1) % N2 + N1;
        double Z = (F[I].first + F[I - 1].first) / 2. + std::floor((i - N1) / N2) * Dz;
        E[i] = E0 + Z * scaleCore_m / mass;
        E2[i] = E[i];
    }
    for (int i = N3 - N1 + 1; i < N3; ++ i) {
        int I = i - N3 - 1 + 2 * N1 + N2;
        double Z = (F[I].first + F[I - 1].first) / 2. + ((numCells_m - 1) * mode_m - 1) * Dz;
        E[i] = E0 + Z * scale_m / mass;
        E2[i] = E[i];
    }

    for (int iter = 0; iter < 10; ++ iter) {
        double A = 0.0;
        double B = 0.0;
        for (int i = 1; i < N1; ++ i) {
            t[i] = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, i, t, 0.0, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, i, t, 0.0, F);
        }
        for (int i = N1; i < N3 - N1 + 1; ++ i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * (getdA(i, I, t, phaseC1, F) + getdA(i, J, t, phaseC2, F));
            B += scaleCore_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * (getdB(i, I, t, phaseC1, F) + getdB(i, J, t, phaseC2, F));
        }
        for (int i = N3 - N1 + 1; i < N3; ++ i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            A += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdA(i, I, t, phaseE, F);
            B += scale_m * (1. + frequency_m * (t2[i] - t[i]) / dphi) * getdB(i, I, t, phaseE, F);
        }

        if (std::abs(B) > 0.0000001) {
            tmp_phi = std::atan(A / B);
        } else {
            tmp_phi = Physics::pi / 2;
        }
        if (q * (A * std::sin(tmp_phi) + B * std::cos(tmp_phi)) < 0) {
            tmp_phi += Physics::pi;
        }

        if (std::abs(phi - tmp_phi) < frequency_m * (t[N3 - 1] - t[0]) / N3) {
            for (int i = 1; i < N1; ++ i) {
                E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            }
            for (int i = N1; i < N3 - N1 + 1; ++ i) {
                int I = (i - N1) % N2 + N1;
                int J = (i - N1 + N4) % N2 + N1;
                E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            }
            for (int i = N3 - N1 + 1; i < N3; ++ i) {
                int I = i - N3 - 1 + 2 * N1 + N2;
                E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            }

            const int prevPrecision = Ippl::Info->precision(8);
            INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
                    << tmp_phi * Units::rad2deg << " deg,\n"
                    << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);
            return tmp_phi;
        }
        phi = tmp_phi - std::round(tmp_phi / Physics::two_pi) * Physics::two_pi;


        for (int i = 1; i < N1; ++ i) {
            E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, i, t, phi + dphi, F); // should I use here t or t2?
            t[i] = t[i - 1] + getdT(i, i, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, i, E2, F, mass);
            E[i] = E[i - 1] + q * scale_m * getdE(i, i, t, phi, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, i, t2, phi + dphi, F);
        }
        for (int i = N1; i < N3 - N1 + 1; ++ i) {
            int I = (i - N1) % N2 + N1;
            int J = (i - N1 + N4) % N2 + N1;
            E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1 + dphi, F) + getdE(i, J, t, phi + phaseC2 + dphi, F)); //concerning t: see above
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i] = E[i - 1] + q * scaleCore_m * (getdE(i, I, t, phi + phaseC1, F) + getdE(i, J, t, phi + phaseC2, F));
            E2[i] = E2[i - 1] + q * scaleCore_m * (getdE(i, I, t2, phi + phaseC1 + dphi, F) + getdE(i, J, t2, phi + phaseC2 + dphi, F));
        }
        for (int i = N3 - N1 + 1; i < N3; ++ i) {
            int I = i - N3 - 1 + 2 * N1 + N2;
            E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE + dphi, F); //concerning t: see above
            t[i] = t[i - 1] + getdT(i, I, E, F, mass);
            t2[i] = t2[i - 1] + getdT(i, I, E2, F, mass);
            E[i] = E[i - 1] + q * scale_m * getdE(i, I, t, phi + phaseE, F);
            E2[i] = E2[i - 1] + q * scale_m * getdE(i, I, t2, phi + phaseE + dphi, F);
        }
        //             msg << ", Ekin= " << E[N3-1] << " MeV" << endl;
    }


    const int prevPrecision = Ippl::Info->precision(8);
    INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad = "
            << tmp_phi * Units::rad2deg << " deg,\n"
            << "Ekin= " << E[N3 - 1] << " MeV" << std::setprecision(prevPrecision) << "\n" << endl);

    return phi;
}

bool TravelingWave::isInside(const Vector_t& r) const {
    return (isInsideTransverse(r)
            && r(2) >= -0.5 * periodLength_m
            && r(2) < startExitField_m);
}
