//
// Class BeamStrippingPhysics
//   Defines the physical processes of residual gas
//   interactions and Lorentz stripping
//
// Copyright (c) 2018 - 2021, Pedro Calvo, CIEMAT, Spain
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
#include "Solvers/BeamStrippingPhysics.h"

#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Vacuum.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Util.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <boost/math/special_functions/chebyshev.hpp>

#include <sys/time.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace {
    struct VacuumInsideTester: public InsideTester {
        explicit VacuumInsideTester(ElementBase* el) {
            vac_m = static_cast<Vacuum*>(el);
        }
        bool checkHit(const Vector_t& R) override {
            return vac_m->checkPoint(R);
        }
    private:
        Vacuum* vac_m;
    };
}


BeamStrippingPhysics::BeamStrippingPhysics(const std::string& name,
                                           ElementBase* element):
    ParticleMatterInteractionHandler(name, element),
    pType_m(ParticleType::UNNAMED),
    T_m(0.0),
    dT_m(0.0),
    mass_m(0.0),
    pressure_m(0.0),
    temperature_m(0.0),
    nCSA_m(0.0),
    nCSB_m(0.0),
    nCSC_m(0.0),
    nCSTotal_m(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    totalPartsInMat_m(0)
{
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));

    const gsl_rng_type* T;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default; // Generator setup
    r_m = gsl_rng_alloc (T);
    gsl_rng_set(r_m, mySeed);
}


BeamStrippingPhysics::~BeamStrippingPhysics() {
    lossDs_m->save();
    gsl_rng_free(r_m);
}


void BeamStrippingPhysics::apply(PartBunchBase<double, 3>* bunch,
                                 const std::pair<Vector_t, double>& /*boundingSphere*/) {

    ParticleType particle = bunch->getPType();
    if (particle != ParticleType::PROTON   &&
        particle != ParticleType::DEUTERON &&
        particle != ParticleType::HMINUS   &&
        particle != ParticleType::H2P      &&
        particle != ParticleType::HYDROGEN) {

        throw GeneralClassicException(
                "BeamStrippingPhysics::apply",
                "Particle " + ParticleProperties::getParticleTypeString(particle) +
                " is not supported for residual stripping interactions!");
    }

    if (element_ref_m->getType() == ElementType::VACUUM) {
         hitTester_m.reset(new VacuumInsideTester(element_ref_m));
         vac_m = dynamic_cast<Vacuum*>(element_ref_m);
    } else {
        throw GeneralClassicException("BeamStrippingPhysics::apply",
                                      "Unsupported element type");
    }

    dT_m = bunch->getdT();

    rediffusedStat_m  = 0;
    stoppedPartStat_m = 0;
    totalPartsInMat_m = 0;

    if (bunch->get_sPos() != 0) {
        doPhysics(bunch);

        gatherStatistics();

        bunch->destroyT();
    }
}


void BeamStrippingPhysics::doPhysics(PartBunchBase<double, 3>* bunch) {
    /*
      Do physics if
      -- particle in material
      -- particle not dead (bunch->Bin[i] != -1)
      Delete particle i: bunch->Bin[i] != -1;
    */
    Inform gmsgALL("OPAL", INFORM_ALL_NODES);

    temperature_m = vac_m->getTemperature();
    bool stop = vac_m->getStop();
    Vector_t extE = Vector_t(0.0, 0.0, 0.0);
    Vector_t extB = Vector_t(0.0, 0.0, 0.0); //kGauss

    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        Vector_t& R = bunch->R[i];
        Vector_t& P = bunch->P[i];
        if ( bunch->Bin[i] != -1 && hitTester_m->checkHit(R) ) {

            ++totalPartsInMat_m;

            bool pdead_GS = false;
            bool pdead_LS = false;
            pressure_m = vac_m->checkPressure(R);
            mass_m = bunch->M[i];
            pType_m = bunch->PType[i];

            double energy = Util::getKineticEnergy(P, mass_m); //GeV
            double gamma = Util::getGamma(P);
            double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
            double deltas = dT_m * beta * Physics::c;

            computeCrossSection(energy);
            pdead_GS = evalGasStripping(deltas);

            if (OpalData::getInstance()->isInOPALCyclMode() &&
                bunch->PType[i] == ParticleType::HMINUS) {
                cycl_m->apply(R, P, T_m, extE, extB);
                double bField = Units::kG2T * std::sqrt(extB[0]*extB[0] + extB[1]*extB[1] + extB[2]*extB[2]); //T
                double eField = gamma * beta * Physics::c * bField;
                pdead_LS = evalLorentzStripping(gamma, eField);
            }

            if (pdead_GS == true || pdead_LS == true) {
                lossDs_m->addParticle(OpalParticle(bunch->ID[i],
                                                   R, P,
                                                   bunch->getT(),
                                                   bunch->Q[i], bunch->M[i]));
                if (stop) {
                    bunch->Bin[i] = -1;
                    ++stoppedPartStat_m;
                    gmsgALL << level4 << getName() << ": Particle " << bunch->ID[i]
                            << " is deleted by beam stripping interactions" << endl;
                } else {
                    ++rediffusedStat_m;
                    getSecondaryParticles(bunch, i, pdead_LS);
                }
            }
        }
    }
}


void BeamStrippingPhysics::computeCrossSection(double energy) {

    const ResidualGas& gas = vac_m->getResidualGas();

    energy *=Units::GeV2keV;
    double energyThreshold = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double a4 = 0.0;
    double a5 = 0.0;
    double a6 = 0.0;
    double a7 = 0.0;
    double a8 = 0.0;
    double a9 = 0.0;
    double a10 = 0.0;
    double a11 = 0.0;
    double a12 = 0.0;
    double molecularDensity[3] = {};

    switch (gas) {
        case ResidualGas::H2: {

            molecularDensity[0] = 100 * pressure_m / (Physics::kB * Physics::q_e * temperature_m);
            double energyMin = 0.0, energyMax = 0.0;
            double csA = 0.0, csB = 0.0, csC = 0.0, csTotal = 0.0;

            if (pType_m == ParticleType::HMINUS) {
                double energyChebyshevFit = energy * Units::keV2eV / (Physics::m_hm / Physics::amu);

                // Single-electron detachment - Hydrogen Production
                energyMin = csCoefSingle_Hminus_Chebyshev[0];
                energyMax = csCoefSingle_Hminus_Chebyshev[1];
                for (int i = 0; i < 9; ++i)
                    a_m[i] = {csCoefSingle_Hminus_Chebyshev[i+2]};
                csA = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

                // Double-electron detachment - Proton Production
                energyMin = csCoefDouble_Hminus_Chebyshev[0];
                energyMax = csCoefDouble_Hminus_Chebyshev[1];
                for (int i = 0; i < 9; ++i)
                    a_m[i] = {csCoefDouble_Hminus_Chebyshev[i+2]};
                csB = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

            } else if (pType_m == ParticleType::PROTON) {
                double energyChebyshevFit = energy * Units::keV2eV / (Physics::m_p / Physics::amu);

                // Single-electron capture - Hydrogen Production
                energyMin = csCoefSingle_Hplus_Chebyshev[0];
                energyMax = csCoefSingle_Hplus_Chebyshev[1];
                for (int i = 0; i < 9; ++i)
                    a_m[i] = {csCoefSingle_Hplus_Chebyshev[i+2]};
                csA = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

                // Double-electron capture - Hminus Production
                energyMin = csCoefDouble_Hplus_Chebyshev[0];
                energyMax = csCoefDouble_Hplus_Chebyshev[1];
                for (int i = 0; i < 9; ++i)
                    a_m[i] = {csCoefDouble_Hplus_Chebyshev[i+2]};
                csB = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

            } else if (pType_m == ParticleType::DEUTERON) {
               // Single-electron capture
               energyThreshold = csCoefSingle_Hplus_Tabata[0];
               a1 = csCoefSingle_Hplus_Tabata[1];
               a2 = csCoefSingle_Hplus_Tabata[2];
               a3 = csCoefSingle_Hplus_Tabata[3];
               a4 = csCoefSingle_Hplus_Tabata[4];
               a5 = csCoefSingle_Hplus_Tabata[5];
               a6 = csCoefSingle_Hplus_Tabata[6];
               a7 = csCoefSingle_Hplus_Tabata[7];
               a8 = csCoefSingle_Hplus_Tabata[8];
               a9 = csCoefSingle_Hplus_Tabata[9];
               csA = computeCrossSectionTabata(energy, energyThreshold, a1, a2, 0.0, 0.0, a3, a4) +
                   computeCrossSectionTabata(energy, energyThreshold, a5, a2, a6, a7, a8, a9);

            } else if (pType_m == ParticleType::HYDROGEN) {
                // Single-electron detachment - Proton Production
                energyThreshold = csCoefProtonProduction_H_Tabata[0];
                a1 = csCoefProtonProduction_H_Tabata[1];
                a2 = csCoefProtonProduction_H_Tabata[2];
                a3 = csCoefProtonProduction_H_Tabata[3];
                a4 = csCoefProtonProduction_H_Tabata[4];
                a5 = csCoefProtonProduction_H_Tabata[5];
                a6 = csCoefProtonProduction_H_Tabata[6];
                csA = computeCrossSectionTabata(energy, energyThreshold, a1, a2, 0.0, 0.0, a3, a4) +
                    a5 * computeCrossSectionTabata(energy/a6, energyThreshold/a6, a1, a2, 0.0, 0.0, a3, a4);

                // Single-electron capture - Hminus Production
                energyThreshold = csCoefHminusProduction_H_Tabata[0];
                a1 = csCoefHminusProduction_H_Tabata[1];
                a2 = csCoefHminusProduction_H_Tabata[2];
                a3 = csCoefHminusProduction_H_Tabata[3];
                a4 = csCoefHminusProduction_H_Tabata[4];
                a5 = csCoefHminusProduction_H_Tabata[5];
                a6 = csCoefHminusProduction_H_Tabata[6];
                a7 = csCoefHminusProduction_H_Tabata[7];
                a8 = csCoefHminusProduction_H_Tabata[8];
                a9 = csCoefHminusProduction_H_Tabata[9];
                a10 = csCoefHminusProduction_H_Tabata[10];
                a11 = csCoefHminusProduction_H_Tabata[11];
                a12 = csCoefHminusProduction_H_Tabata[12];
                csB = computeCrossSectionTabata(energy, energyThreshold, a1, a2, a3, a4, a5, a6) +
                    computeCrossSectionTabata(energy, energyThreshold, a7, a8, a9, a10, a11, a12);

            } else if (pType_m == ParticleType::H2P) {
                double energyChebyshevFit = energy * Units::keV2eV / (Physics::m_h2p / Physics::amu);
                // Proton production
                if (energy <= energyRangeH2plusinH2[0]) {
                    energyThreshold = csCoefProtonProduction_H2plus_Tabata[0];
                    a1 = csCoefProtonProduction_H2plus_Tabata[1];
                    a2 = csCoefProtonProduction_H2plus_Tabata[2];
                    a3 = csCoefProtonProduction_H2plus_Tabata[3];
                    a4 = csCoefProtonProduction_H2plus_Tabata[4];
                    a5 = csCoefProtonProduction_H2plus_Tabata[5];
                    a6 = csCoefProtonProduction_H2plus_Tabata[6];
                    a7 = csCoefProtonProduction_H2plus_Tabata[7];
                    a8 = csCoefProtonProduction_H2plus_Tabata[8];
                    a9 = csCoefProtonProduction_H2plus_Tabata[9];
                    a10 = csCoefProtonProduction_H2plus_Tabata[10];
                    csA = computeCrossSectionTabata(energy, energyThreshold, a1, a2, 0.0, 0.0, a3, a4) +
                        computeCrossSectionTabata(energy, energyThreshold, a5, a6, 0.0, 0.0, a7, a8) +
                        a9 * computeCrossSectionTabata(energy/a10, energyThreshold/a10, a5, a6, 0.0, 0.0, a7, a8);

                } else if (energy > energyRangeH2plusinH2[0] &&
                           energy < energyRangeH2plusinH2[1]) {
                    energyMin = csCoefProtonProduction_H2plus_Chebyshev[0];
                    energyMax = csCoefProtonProduction_H2plus_Chebyshev[1];
                    for (int i = 0; i < 9; ++i)
                        a_m[i] = {csCoefProtonProduction_H2plus_Chebyshev[i+2]};
                    csA = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

                } else if (energy >= energyRangeH2plusinH2[1]) {
                    int zTarget = 1;
                    double massInAmu = Physics::m_h2p / Physics::amu;
                    csA = computeCrossSectionBohr(energy, zTarget, massInAmu);
                }

                // Hydrogen production
                energyMin = csCoefHydrogenProduction_H2plus_Chebyshev[0];
                energyMax = csCoefHydrogenProduction_H2plus_Chebyshev[1];
                for (int i = 0; i < 9; ++i)
                    a_m[i] = {csCoefHydrogenProduction_H2plus_Chebyshev[i+2]};
                csB = computeCrossSectionChebyshev(energyChebyshevFit, energyMin, energyMax);

                // H3+, H
                energyThreshold = csCoefH3plusProduction_H2plus_Tabata[0];
                a1 = csCoefH3plusProduction_H2plus_Tabata[1];
                a2 = csCoefH3plusProduction_H2plus_Tabata[2];
                a3 = csCoefH3plusProduction_H2plus_Tabata[3];
                a4 = csCoefH3plusProduction_H2plus_Tabata[4];
                a5 = csCoefH3plusProduction_H2plus_Tabata[5];
                a6 = csCoefH3plusProduction_H2plus_Tabata[6];
                csC = computeCrossSectionTabata(energy, energyThreshold, a1, a2, a3, a4, a5, a6);
            }
            csTotal = csA + csB + csC;

            nCSA_m = csA * Units::cm2m * Units::cm2m * molecularDensity[0];
            nCSB_m = csB * Units::cm2m * Units::cm2m * molecularDensity[0];
            nCSC_m = csC * Units::cm2m * Units::cm2m * molecularDensity[0];
            nCSTotal_m = csTotal * Units::cm2m * Units::cm2m * molecularDensity[0];
            break;
        }

        case ResidualGas::AIR: {

            int zTarget[3] = {7, 8, 18};
            static const double fMolarFraction[3] = {0.78084, 0.20947, 0.00934}; // N2, O2, Ar
            double csSingle[3], csDouble[3], csTotal[3];
            double nCSSingle[3], nCSDouble[3], nCS[3];
            double nCSSingleSum = 0.0;
            double nCSDoubleSum = 0.0;
            double nCSTotalSum = 0.0;

            for (int i = 0; i < 3; ++i) {
                molecularDensity[i] = 100 * pressure_m * fMolarFraction[i] / (Physics::kB * Physics::q_e * temperature_m);

                if (pType_m == ParticleType::HMINUS) {
                    // Single-electron detachment - Hydrogen Production
                    energyThreshold = csCoefSingle_Hminus[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefSingle_Hminus[i][j+1];
                    csSingle[i] = computeCrossSectionNakai(energy, energyThreshold, i);

                    // Double-electron detachment - Proton Production
                    energyThreshold = csCoefDouble_Hminus[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefDouble_Hminus[i][j+1];
                    csDouble[i] = computeCrossSectionNakai(energy, energyThreshold, i);

                } else if (pType_m == ParticleType::PROTON || pType_m == ParticleType::DEUTERON) {
                    // Single-electron capture
                    energyThreshold = csCoefSingle_Hplus[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefSingle_Hplus[i][j+1];
                    csSingle[i] = computeCrossSectionNakai(energy, energyThreshold, i) +
                        b_m[i][6] * computeCrossSectionNakai(energy/b_m[i][7], energyThreshold/b_m[i][7], i);

                    // Double-electron capture
                    energyThreshold = csCoefDouble_Hplus[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefDouble_Hplus[i][j+1];
                    if (b_m[i][7] != 0) {
                        csDouble[i] = computeCrossSectionNakai(energy, energyThreshold, i) +
                            b_m[i][6] * computeCrossSectionNakai(energy/b_m[i][7], energyThreshold/b_m[i][7], i);
                    } else {
                        csDouble[i] = computeCrossSectionNakai(energy, energyThreshold, i);
                    }

                } else if (pType_m == ParticleType::HYDROGEN) {
                    // Single-electron detachment - Proton Production
                    energyThreshold = csCoefSingleLoss_H[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefSingleLoss_H[i][j+1];
                    csSingle[i] = computeCrossSectionNakai(energy, energyThreshold, i);

                    // Single-electron capture - Hminus Production
                    energyThreshold = csCoefSingleCapt_H[i][0];
                    for (int j = 0; j < 8; ++j)
                        b_m[i][j] = csCoefSingleCapt_H[i][j+1];
                    if (b_m[i][7] != 0) {
                        csDouble[i] = computeCrossSectionNakai(energy, energyThreshold, i) +
                            b_m[i][6] * computeCrossSectionNakai(energy/b_m[i][7], energyThreshold/b_m[i][7], i);
                    } else {
                        csDouble[i] = computeCrossSectionNakai(energy, energyThreshold, i);
                    }

                } else if (pType_m == ParticleType::H2P) {
                    double massInAmu = Physics::m_h2p / Physics::amu;
                    csSingle[i] = {0.0};
                    csDouble[i] = computeCrossSectionBohr(energy, zTarget[i], massInAmu);

                } else {
                    csSingle[i] = {0.0};
                    csDouble[i] = {0.0};
                }
                csTotal[i] = csSingle[i] + csDouble[i];
                nCSSingle[i] = csSingle[i] * Units::cm2m * Units::cm2m * molecularDensity[i];
                nCSDouble[i] = csDouble[i] * Units::cm2m * Units::cm2m * molecularDensity[i];
                nCS[i] = csTotal[i] * Units::cm2m * Units::cm2m * molecularDensity[i];

                nCSSingleSum += nCSSingle[i];
                nCSDoubleSum += nCSDouble[i];
                nCSTotalSum += nCS[i];
            }
            nCSA_m = nCSSingleSum;
            nCSB_m = nCSDoubleSum;
            nCSTotal_m = nCSTotalSum;
            break;
        }
        default: {
            break;
        }
    }
}

/// Analytical expression for charge-transfer cross section between hydrogen
/// ions and atoms with gaseous atoms and molecules based on semiempirically
/// functional forms.
/// Y.Nakai et.al, "Cross sections for charge transfer of hydrogen atoms and
/// ions colliding with gaseous atoms and molecules", At. Data Nucl. Data
/// Tables 37, 69 (1987).
// -------------------------------------------------------------------------
double BeamStrippingPhysics::computeCrossSectionNakai(double energy, double energyThreshold, int &i) {

    if (energy <= energyThreshold) {
        return 0.0;
    }
    const double E_R = Physics::E_ryd * Units::GeV2keV * Physics::m_h / Physics::m_e;
    const double sigma_0 = 1E-16;
    double sigma = 0.0; //cm2
    double effectiveEnergy = energy - energyThreshold; //keV
    double f1 = sigma_0 * b_m[i][0] * std::pow((effectiveEnergy / E_R), b_m[i][1]);
    if (b_m[i][2] != 0.0 && b_m[i][3] != 0.0) {
        sigma = f1 / (1 + std::pow((effectiveEnergy / b_m[i][2]), (b_m[i][1] + b_m[i][3]))
                + std::pow((effectiveEnergy / b_m[i][4]), (b_m[i][1] + b_m[i][5])));
    } else {
        sigma = f1 / (1 + std::pow((effectiveEnergy / b_m[i][4]), (b_m[i][1] + b_m[i][5])));
    }

    return sigma;
}

/// Analytical expression for cross section of collision processes between
/// hydrogen ions, atoms, and molecules with hydrogen molecules based on
/// semiempirically functional forms.
/// T.Tabata and T.Shirai, "Analytic cross sections for collisions of H+,
/// H2+, H3+, H, H2 and H− with hydrogen molecules", At. Data Nucl. Data
/// Tables 76, 1 (2000).
// -------------------------------------------------------------------------
double BeamStrippingPhysics::computeCrossSectionTabata(double energy, double energyThreshold,
                                                       double a1, double a2, double a3,
                                                       double a4, double a5, double a6) {
    if (energy <= energyThreshold) {
        return 0.0;
    }
    const double sigma_0 = 1E-16;
    double sigma = 0.0; //cm2
    double effectiveEnergy = energy - energyThreshold; //keV
    double f1 = sigma_0 * a1 * std::pow((effectiveEnergy / (Physics::E_ryd * Units::GeV2keV)),a2);
    if (a3 != 0.0 && a4 != 0.0) {
        sigma = f1 / (1 + std::pow((effectiveEnergy / a3), (a2 + a4)) + std::pow((effectiveEnergy / a5), (a2 + a6)));
    } else {
        sigma = f1 / (1 + std::pow((effectiveEnergy / a5), (a2 + a6)));
    }

    return sigma;
}

/// Analytical expression for cross section for inelastic collisions between
/// hydrogen atoms, molecules and ions from Chebyshev fitting parameters.
/// C.F.Barnett, "Atomic data for fusion. Volume 1: Collisions of H, H2, He
/// and Li atoms and ions with atoms and molecules", ORNL-6068/V1 (1990).
// -------------------------------------------------------------------------
double BeamStrippingPhysics::computeCrossSectionChebyshev(double energy,
                                                          double energyMin,
                                                          double energyMax) {
    // energy -> eV/amu
    if (energy <= energyMin || energy >= energyMax) {
        return 0.0;
    }
    double sum_aT = 0.0;
    double aT[8] = {0.0};
    double x = ((std::log(energy)-std::log(energyMin)) - (std::log(energyMax)-std::log(energy))) / (std::log(energyMax)-std::log(energyMin));
    for (int i = 0; i < 8; ++i) {
        aT[i] = (a_m[i+1] * boost::math::chebyshev_t(i+1, x));
        sum_aT += aT[i];
    }
    double sigma = std::exp(0.5*a_m[0] + sum_aT); //cm2

    return sigma;
}

/// Analytical expressions for electron loss by light and heavy ions passing
/// through gases based on Bohr theory.
/// H.-D.Betz, "Charge states and charge-changing cross sections of fast heavy
/// ions penetrating through gaseous and solid media", Rev. Mod. Phys. 44,
/// 465 (1972).
// -------------------------------------------------------------------------
double BeamStrippingPhysics::computeCrossSectionBohr(double energy, int zTarget, double massInAmu) {

    double energyAmu = energy / massInAmu;
    if (energyAmu <= 1E3 || energyAmu >= 1E5) {
        return 0.0;
    }
    double sigma = 0.0; //cm2
    constexpr double GeV2keV = Units::GeV2eV * Units::eV2keV;
    double mass = mass_m * GeV2keV;
    const double a0v0 = Physics::h_bar * Physics::c * Physics::c / Physics::m_e;
    double v = Physics::c * std::sqrt(1 - std::pow(mass/(energy+mass), 2));
    if (zTarget < 15) {
        double z = (zTarget + zTarget*zTarget);
        sigma = Units::m2cm * Units::m2cm * 4 * Physics::pi * std::pow(a0v0, 2) * z / std::pow(v, 2);
    } else {
        sigma = Units::m2cm * Units::m2cm * 4 * Physics::pi * Physics::a0 * a0v0 * std::pow(zTarget, 2.0/3.0) / v;
    }
    return sigma;
}


bool BeamStrippingPhysics::evalGasStripping(double& deltas) {

    double xi = gsl_rng_uniform(r_m);
    double fg = 1-std::exp(-nCSTotal_m*deltas);

    return (fg >= xi);
}

/// Expression for H- lifetime based on the calculation of the electric
/// dissociation rate from the formal theory of decay.
/// L.R.Scherk, "A improved value for the electron affinity of the negative
/// hydrogen ion", Can. J. Phys. 57, 558 (1979).
// -------------------------------------------------------------------------
bool BeamStrippingPhysics::evalLorentzStripping(double& gamma, double& eField) {

    double xi = gsl_rng_uniform(r_m);

    const double eps0 = 0.75419 * Physics::q_e; // electron binding energy,
    const double hbar = Physics::h_bar * Units::GeV2eV * Physics::q_e; // Js
    const double me = Physics::m_e * Units::GeV2kg;
    const double p = 0.0126; // polarization factor
    const double s0 = 0.783; // spectroscopic coefficient
    const double a = 2.01407/Physics::a0;
    const double k0 = std::sqrt(2 * me * eps0)/hbar;
    const double n = (std::sqrt(2 * k0 * (k0+a) * (2*k0+a)))/a; // normalization factor
    double zT = eps0 / (Physics::q_e * eField);
    double tau = (4 * me * zT)/(s0 * n * n * hbar * (1+p)*(1+p) * (1-1/(2*k0*zT))) * std::exp(4*k0*zT/3);
    double fL = 1 - std::exp( -dT_m / (gamma * tau));

    return (fL >= xi);
}

void BeamStrippingPhysics::getSecondaryParticles(PartBunchBase<double, 3>* bunch,
                                                 size_t& i, bool pdead_LS) {
    double r = gsl_rng_uniform(r_m);

    const ResidualGas& gas = vac_m->getResidualGas();

    if (pType_m == ParticleType::HMINUS) {
        if (pdead_LS == true) {
            transformToSecondary(bunch, i, ParticleType::HYDROGEN);
        } else {
            if (r > nCSB_m/nCSTotal_m)
                transformToSecondary(bunch, i, ParticleType::HYDROGEN);
            else
                transformToSecondary(bunch, i, ParticleType::PROTON);
        }

    } else if (pType_m == ParticleType::PROTON) {
        if (r > nCSB_m/nCSTotal_m)
            transformToSecondary(bunch, i, ParticleType::HYDROGEN);
        else
            transformToSecondary(bunch, i, ParticleType::HMINUS);

    } else if (pType_m == ParticleType::HYDROGEN) {
        if (r > nCSB_m/nCSTotal_m)
            transformToSecondary(bunch, i, ParticleType::PROTON);
        else
            transformToSecondary(bunch, i, ParticleType::HMINUS);

    } else if (pType_m == ParticleType::H2P) {
        if (gas == ResidualGas::H2) {
            if (nCSC_m>nCSB_m && nCSB_m>nCSA_m) {
                if (r > (nCSA_m+nCSB_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::H3P);
                else if (r > nCSA_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);
                else
                    transformToSecondary(bunch, i, ParticleType::PROTON);

            } else if (nCSA_m>nCSB_m && nCSB_m>nCSC_m) {
                if (r > (nCSC_m+nCSB_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::PROTON);
                else if (r > nCSC_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);
                else
                    transformToSecondary(bunch, i, ParticleType::H3P);

            } else if (nCSA_m>nCSB_m && nCSC_m>nCSA_m) {
                if (r > (nCSA_m+nCSB_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::H3P);
                else if (r > nCSB_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::PROTON);
                else
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);

            } else if (nCSA_m>nCSC_m && nCSC_m>nCSB_m) {
                if (r > (nCSC_m+nCSB_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::PROTON);
                else if (r > nCSB_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::H3P);
                else
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);

            } else if (nCSB_m>nCSC_m && nCSB_m>nCSA_m && nCSA_m>nCSC_m) {
                if (r > (nCSC_m+nCSA_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);
                else if (r > nCSC_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::PROTON);
                else
                    transformToSecondary(bunch, i, ParticleType::H3P);

            } else {
                if (r > (nCSC_m+nCSA_m)/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::HYDROGEN);
                else if (r > nCSA_m/nCSTotal_m)
                    transformToSecondary(bunch, i, ParticleType::H3P);
                else
                    transformToSecondary(bunch, i, ParticleType::PROTON);
            }
        } else if (gas == ResidualGas::AIR) {
            if (r > nCSTotal_m)
                transformToSecondary(bunch, i, ParticleType::PROTON);
        }
    } else if (pType_m == ParticleType::DEUTERON) {
        GeneralClassicException("BeamStrippingPhysics::getSecondaryParticles",
                                "Tracking secondary particles from incident "
                                "ParticleType::DEUTERON is not implemented");
    }
}

void BeamStrippingPhysics::transformToSecondary(PartBunchBase<double, 3>* bunch,
                                                size_t& i,
                                                ParticleType type) {
    bunch->POrigin[i] = ParticleOrigin::SECONDARY;
    bunch->PType[i] = type;
    bunch->M[i] = ParticleProperties::getParticleMass(type);
    bunch->Q[i] = ParticleProperties::getParticleChargeInCoulomb(type);

    Inform gmsgALL("OPAL", INFORM_ALL_NODES);
    gmsgALL << level4 << getName() << ": Particle " << bunch->ID[i]
            << " is transformed to " << ParticleProperties::getParticleTypeString(type) << endl;
}


void BeamStrippingPhysics::print(Inform& msg) {
    Inform::FmtFlags_t ff = msg.flags();
    if (totalPartsInMat_m > 0) {
        msg << level2
            << "\n"<< "--- BeamStrippingPhysics ---\n"
            << "Name: " << name_m << " - "
            << "Element: " << element_ref_m->getName() << " - "
            << "Residual gas: " << vac_m->getResidualGasName() << "\n"
            << std::setw(31) << "Particles in the material: " << Util::toStringWithThousandSep(totalPartsInMat_m) << "\n"
            << std::setw(31) << "Rediffused as secondary: " << Util::toStringWithThousandSep(rediffusedStat_m) << "\n"
            << std::setw(31) << "Stripped by the gas: " << Util::toStringWithThousandSep(stoppedPartStat_m)
            << endl;
    }
    msg.flags(ff);
}

void BeamStrippingPhysics::gatherStatistics() {

    constexpr unsigned short numItems = 3;
    unsigned int partStatistics[numItems] = {totalPartsInMat_m,
                                             rediffusedStat_m,
                                             stoppedPartStat_m};

    allreduce(&(partStatistics[0]), numItems, std::plus<unsigned int>());

    totalPartsInMat_m = partStatistics[0];
    rediffusedStat_m = partStatistics[1];
    stoppedPartStat_m = partStatistics[2];
}


/*
    Cross sections parameters for interaction with air
    -- [1] -> Nitrogen
    -- [2] -> Oxygen
    -- [3] -> Argon
*/
const double BeamStrippingPhysics::csCoefSingle_Hminus[3][9] = {
    {7.50E-04, 4.38E+02, 7.28E-01, 8.40E-01, 2.82E-01, 4.10E+01, 1.37E+00, 0.00E+00, 0.00E+00},
    {-2.00E-04, 3.45E+02, 4.80E-01, 5.30E-02, 8.40E-02, 1.00E+01, 9.67E-01, 0.00E+00, 0.00E+00},
    {1.70E-3, 2.47E+01, 3.36E-01, 7.00E+01, 5.00E-01, 9.90E+01, 7.80E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::csCoefDouble_Hminus[3][9] = {
    {1.40E-02, 1.77E+00, 4.80E-01, 0.00E+00, 0.00E+00, 1.52E+02, 1.52E+00, 0.00E+00, 0.00E+00},
    {1.30E-02, 1.90E+00, 6.20E-01, 0.00E+00, 0.00E+00, 5.20E+01, 9.93E-01, 0.00E+00, 0.00E+00},
    {1.50E-02, 1.97E+00, 8.90E-01, 0.00E+00, 0.00E+00, 5.10E+01, 9.37E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::csCoefSingle_Hplus[3][9] = {
    {2.00E-03, 1.93E+03, 1.64E+00, 1.98E+00, 6.69E-01, 2.19E+01, 4.15E+00, 3.23E-04, 1.00E+01},
    {-1.50E-03, 3.86E+05, 1.60E+00, 6.93E-02, 3.28E-01, 7.86E+00, 3.92E+00, 3.20E-04, 1.00E+01},
    {2.18E-03, 1.61E+04, 2.12E+00, 1.16E+00, 4.44E-01, 1.39E+01, 4.07E+00, 2.99E-04, 1.45E+01}
};
const double BeamStrippingPhysics::csCoefDouble_Hplus[3][9] = {
    {2.90E-02, 2.90E-01, 1.50E+00, 1.39E+01, 1.65E+00, 3.94E+01, 5.79E+00, 0.00E+00, 0.00E+00},
    {2.20E-02, 2.64E-01, 1.50E+00, 1.40E+01, 1.70E+00, 4.00E+01, 5.80E+00, 0.00E+00, 0.00E+00},
    {2.90E-02, 1.40E-01, 1.87E+00, 2.40E+01, 1.60E+00, 3.37E+01, 4.93E+00, 4.50E-01, 2.00E-01}
};
const double BeamStrippingPhysics::csCoefSingleLoss_H[3][9] = {
    {1.36E-02, 2.59E+03, 1.78E+00, 2.69E-01, -3.52E-01, 5.60E+00, 8.99E-01, 0.00E+00, 0.00E+00},
    {1.36E-02, 3.29E+03, 1.85E+00, 2.89E-01, -2.72E-01, 8.20E+00, 1.09E+00, 0.00E+00, 0.00E+00},
    {1.36E-02, 1.16E+12, 7.40E+00, 5.36E-01, -5.42E-01, 1.23E+00, 7.26E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::csCoefSingleCapt_H[3][9] = {
    {1.48E-02, 1.15E+00, 1.18E+00, 1.15E+01, 1.19E+00, 3.88E+01, 3.38E+00, 1.00E-01, 2.82E-02},
    {1.13E-02, 3.26E+00, 1.02E+00, 2.99E+00, 4.50E-01, 1.90E+01, 3.42E+00, 0.00E+00, 0.00E+00},
    {1.50E-02, 1.24E+03, 3.38E+00, 2.46E+00, 5.20E-01, 7.00E+00, 2.56E+00, 9.10E-02, 1.95E-02}
};

// Cross sections parameters for interaction with hydrogen gas (H2)
const double BeamStrippingPhysics::csCoefSingle_Hplus_Tabata[10] = {
    2.50E-03, 2.12E+02, 1.721E+00, 6.70E-04, 3.239E-01, 4.34E-03, 1.296E+00, 1.42E-01, 9.34E+00, 2.997E+00
};
const double BeamStrippingPhysics::csCoefHminusProduction_H_Tabata[13] = {
    2.1E-02, 9.73E-03, 2.38E+00, 1.39E-02, -5.51E-01, 7.7E-02, 2.12E+00, 1.97E-06, 2.051E+00, 5.5E+00, 6.62E-01, 2.02E+01, 3.62E+00
};
const double BeamStrippingPhysics::csCoefProtonProduction_H_Tabata[9] = {
    2.0E-02, 2.53E-04, 1.728E+00, 2.164E+00, 7.74E-01, 1.639E+00, 1.43E+01, 0, 0
};
const double BeamStrippingPhysics::csCoefProtonProduction_H2plus_Tabata[11] = {
    5.0E-03, 6.34E+01, 1.78E+00, 1.38E-03, 4.06E-01, 1.63E-01, 3.27E-01, 1.554E+01, 3.903E+00, 1.735, 1.02E+01
};
const double BeamStrippingPhysics::csCoefH3plusProduction_H2plus_Tabata[7] = {
    0.00E+00, 6.05E+00, -5.247E-01, 4.088E-03, 2.872E+00, 7.3E-03, 6.99E+00
};
const double BeamStrippingPhysics::csCoefSingle_Hminus_Chebyshev[11] = {
    2.3, 1.7E+07, -73.1505813599, -1.7569509745, -2.0016760826, -0.1902804971, 0.0171353221, 0.1270833164, -0.1523126215, 0, 0
};
const double BeamStrippingPhysics::csCoefDouble_Hminus_Chebyshev[11] = {
    1.00E+03, 1.00E+07, -79.0158996582, -2.1025185585, -1.2407282591, 0.174798578, 0.1062489152, -0.0004342734, -0.0465673618, 0, 0
};
const double BeamStrippingPhysics::csCoefSingle_Hplus_Chebyshev[11] = {
    2.6E+00, 4.00E+06, -82.5164, -6.70755, -6.10977, -2.6281, 0.709759, 0.639033, 0.10298, 0.26124, -0.263817
};
const double BeamStrippingPhysics::csCoefDouble_Hplus_Chebyshev[11] = {
    200, 1.00E+06, -95.8165, -7.17049, -7.48288, -1.93034, 0.761153, 0.556689, -0.0542859, -0.270184, -0.0147
};
const double BeamStrippingPhysics::csCoefHydrogenProduction_H2plus_Chebyshev[11] = {
    2.00E+03, 1.00E+05, -70.670173645, -0.632612288, -0.6065212488, -0.0915143117, -0.0121710282, 0.0168179292, 0.0104796877, 0, 0
};
const double BeamStrippingPhysics::csCoefProtonProduction_H2plus_Chebyshev[11] = {
    1.50E+03, 1.00E+07, -74.9261474609, -2.1944284439, -0.8558337688, 0.0421306863, 0.2162267119, 0.0921146944, -0.0893079266, 0, 0
};
const double BeamStrippingPhysics::energyRangeH2plusinH2[2] = {111.4, 7.8462E+04};
double BeamStrippingPhysics::a_m[9] = {};
double BeamStrippingPhysics::b_m[3][9] = {};