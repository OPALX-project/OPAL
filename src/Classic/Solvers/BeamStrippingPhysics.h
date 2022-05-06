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
#ifndef BEAMSTRIPPINGPHYSICS_H
#define BEAMSTRIPPINGPHYSICS_H

#include "AbsBeamline/Component.h"
#include "Algorithms/PartBunchBase.h"
#include "Algorithms/Vektor.h"
#include "Solvers/ParticleMatterInteractionHandler.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include <memory>
#include <string>
#include <utility>

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;
class Inform;
class Cyclotron;
class Vacuum;

class BeamStrippingPhysics: public ParticleMatterInteractionHandler {

public:
    BeamStrippingPhysics(const std::string& name, ElementBase* element);
    ~BeamStrippingPhysics();

    void setCyclotron(Cyclotron* cycl) { cycl_m = cycl; };

    virtual void apply(PartBunchBase<double, 3>* bunch,
                       const std::pair<Vector_t, double>& boundingSphere) override;

    virtual const std::string getType() const override;
    virtual void print(Inform& msg) override;
    virtual bool stillActive() override;

    virtual double getTime() override;
    virtual std::string getName() override;
    virtual size_t getParticlesInMat() override;
    virtual unsigned getRediffused() override;
    virtual unsigned int getNumEntered() override;

    void doPhysics(PartBunchBase<double, 3>* bunch);

private:

    void computeCrossSection(double energy);

    double computeCrossSectionNakai(double energy, double energyThreshold, int& i);
    
    double computeCrossSectionTabata(double energy, double energyThreshold,
                                     double a1, double a2, double a3,
                                     double a4, double a5, double a6);
                              
    double computeCrossSectionChebyshev(double energy, double energyMin, double energyMax);

    double computeCrossSectionBohr(double energy, int zTarget, double massInAmu);

    bool evalGasStripping(double& deltas);
    bool evalLorentzStripping(double& gamma, double& eField);

    void getSecondaryParticles(PartBunchBase<double, 3>* bunch, size_t& i, bool pdead_LS);

    void transformToSecondary(PartBunchBase<double, 3>* bunch, size_t& i, ParticleType type);

    bool computeEnergyLoss(PartBunchBase<double, 3>* /*bunch*/,
                           Vector_t& /*P*/,
                           const double /*deltat*/,
                           bool /*includeFluctuations*/) const {
        return false;
    }

    void gatherStatistics();

    ParticleType pType_m;

    Vacuum* vac_m;
    Cyclotron* cycl_m;

    gsl_rng* r_m;

    double T_m;  // s
    double dT_m; // s
    double mass_m; // GeV/c2
    double pressure_m; // mbar
    double temperature_m; // K

    std::unique_ptr<LossDataSink> lossDs_m;

    /// macroscopic cross sections
    double nCSA_m;
    double nCSB_m;
    double nCSC_m;
    double nCSTotal_m;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    unsigned totalPartsInMat_m;

    static const double csCoefSingle_Hminus[3][9];
    static const double csCoefDouble_Hminus[3][9];
    static const double csCoefSingle_Hplus[3][9];
    static const double csCoefDouble_Hplus[3][9];
    static const double csCoefSingleLoss_H[3][9];
    static const double csCoefSingleCapt_H[3][9];

    static const double csCoefSingle_Hplus_Tabata[10];
    static const double csCoefHminusProduction_H_Tabata[13];
    static const double csCoefProtonProduction_H_Tabata[9];
    static const double csCoefProtonProduction_H2plus_Tabata[11];
    static const double csCoefH3plusProduction_H2plus_Tabata[7];

    static const double csCoefSingle_Hminus_Chebyshev[11];
    static const double csCoefDouble_Hminus_Chebyshev[11];
    static const double csCoefSingle_Hplus_Chebyshev[11];
    static const double csCoefDouble_Hplus_Chebyshev[11];
    static const double csCoefHydrogenProduction_H2plus_Chebyshev[11];
    static const double csCoefProtonProduction_H2plus_Chebyshev[11];
    static const double energyRangeH2plusinH2[2];
    static double a_m[9];
    static double b_m[3][9];
};


inline
double BeamStrippingPhysics::getTime() {
    return T_m;
}

inline
std::string BeamStrippingPhysics::getName() {
    return (element_ref_m->getName() + "_" + name_m);
}

inline
size_t BeamStrippingPhysics::getParticlesInMat() {
    return totalPartsInMat_m;
}

inline
unsigned int BeamStrippingPhysics::getRediffused() {
    return rediffusedStat_m;
}

inline
unsigned int BeamStrippingPhysics::getNumEntered() {
    return bunchToMatStat_m;
}

inline
bool BeamStrippingPhysics::stillActive() {
    return totalPartsInMat_m != 0;
}

inline
const std::string BeamStrippingPhysics::getType() const {
    return "BeamStrippingPhysics";
}

#endif //BEAMSTRIPPINGPHYSICS_H
