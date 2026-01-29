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
#ifndef BEAMSTRIPPINGPHYSICS_HH
#define BEAMSTRIPPINGPHYSICS_HH

#include "AbsBeamline/Component.h"
#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"
#include "Solvers/ParticleMatterInteractionHandler.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

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
                       const std::pair<Vector_t, double>& boundingSphere);

    virtual const std::string getType() const;
    virtual void print(Inform& msg);
    virtual bool stillActive();

    virtual double getTime();
    virtual std::string getName();
    virtual size_t getParticlesInMat();
    virtual unsigned getRediffused();
    virtual unsigned int getNumEntered();

    inline void doPhysics(PartBunchBase<double, 3>* bunch);

private:

    void computeCrossSection(PartBunchBase<double, 3>* bunch, size_t& i, double energy);

    double computeCrossSectionNakai(double energy, double energyThreshold, int& i);
    
    double computeCrossSectionTabata(double energy, double energyThreshold,
                                    double a1, double a2, double a3,
                                    double a4, double a5, double a6);
                              
    double computeCrossSectionChebyshev(double energy, double energyMin, double energyMax);

    double computeCrossSectionBohr(double energy, int zTarget, double massInAmu);

    bool evalGasStripping(double& deltas);
    bool evalLorentzStripping(double& gamma, double& eField);

    void getSecondaryParticles(PartBunchBase<double, 3>* bunch, size_t& i, bool pdead_LS);
    void transformToProton(PartBunchBase<double, 3>* bunch, size_t& i);
    void transformToHydrogen(PartBunchBase<double, 3>* bunch, size_t& i);
    void transformToHminus(PartBunchBase<double, 3>* bunch, size_t& i);
    void transformToH3plus(PartBunchBase<double, 3>* bunch, size_t& i);

    bool computeEnergyLoss(PartBunchBase<double, 3>* /*bunch*/,
                           Vector_t& /*P*/,
                           const double /*deltat*/,
                           bool /*includeFluctuations*/) const {
        return false;
    }

    Cyclotron* cycl_m;
    Vacuum* vac_m;

    gsl_rng* r_m;

    double T_m;  // s
    double dT_m; // s
    double mass_m;
    double charge_m;
    double pressure_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    /// macroscopic cross sections
    double nCSA;
    double nCSB;
    double nCSC;
    double nCSTotal;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t totalPartsInMat_m;

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
const std::string BeamStrippingPhysics::getType() const {
    return "BeamStrippingPhysics";
}

#endif //BEAMSTRIPPINGPHYSICS_HH
