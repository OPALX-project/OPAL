//
// Class BeamStrippingPhysics
//   This class provides beam stripping physical processes as
//   particle matter interaction type.
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
#ifndef BEAMSTRIPPINGPHYSICS_HH
#define BEAMSTRIPPINGPHYSICS_HH

#include "AbsBeamline/Component.h"
#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

template <class T, unsigned Dim>

class PartBunchBase;
class LossDataSink;
class Inform;
class Cyclotron;
class BeamStripping;

class BeamStrippingPhysics: public ParticleMatterInteractionHandler {

public:

    BeamStrippingPhysics(const std::string &name, ElementBase *element);
    ~BeamStrippingPhysics();

    void setCyclotron(Cyclotron* cycl) { cycl_m = cycl; };

    virtual void apply(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere);

    virtual const std::string getType() const;
    void print(Inform& msg);
    bool stillActive();
    bool stillAlive(PartBunchBase<double, 3> *bunch);

    inline double getTime() {return T_m;}
    std::string getName() {return element_ref_m->getName();}
    size_t getParticlesInMat() {return locPartsInMat_m;}
    unsigned getRediffused() {return rediffusedStat_m;}
    unsigned int getNumEntered() {return bunchToMatStat_m;}
    inline void doPhysics(PartBunchBase<double, 3> *bunch);

private:

    void crossSection(double Eng);

    double csAnalyticFunctionNakai(double Eng, double Eth, int &i);
    
    double csAnalyticFunctionTabata(double Eng, double Eth,
                            double a1, double a2, double a3, double a4, double a5, double a6);
                              
    double csChebyshevFitting(double Eng, double Emin, double Emax);

    bool gasStripping(double &deltas);

    bool lorentzStripping(double &gamma, double &E);

    void secondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i, bool pdead_LS);
    void transformToProton(PartBunchBase<double, 3> *bunch, size_t &i);
    void transformToHydrogen(PartBunchBase<double, 3> *bunch, size_t &i);
    void transformToHminus(PartBunchBase<double, 3> *bunch, size_t &i);
    void transformToH3plus(PartBunchBase<double, 3> *bunch, size_t &i);

    bool computeEnergyLoss(Vector_t &/*P*/, const double /*deltat*/, bool /*includeFluctuations*/) const {
        return false;
    }

    Cyclotron *cycl_m;
    BeamStripping *bstp_m;
    ElementBase::ElementType bstpshape_m;

    gsl_rng * r_m;

    double T_m;    /// s
    double dT_m;   /// s

    double mass_m;
    double charge_m;

    std::string gas_m;
    double pressure_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    double NCS_a;
    double NCS_b;
    double NCS_c;
    double NCS_total;

    unsigned bunchToMatStat_m;
    unsigned stoppedPartStat_m;
    unsigned rediffusedStat_m;
    size_t locPartsInMat_m;

    static const double csCoefSingle_Hminus[3][9];
    static const double csCoefDouble_Hminus[3][9];
    static const double csCoefSingle_Hplus[3][9];
    static const double csCoefDouble_Hplus[3][9];
    static const double csCoefSingleLoss_H[3][9];
    static const double csCoefSingleCapt_H[3][9];

    static const double csCoefHminusProduction_H_Tabata[13];
    static const double csCoefProtonProduction_H_Tabata[9];
    static const double csCoefProtonProduction_H2plus_Tabata[11];
    static const double csCoefH3plusProduction_H2plus_Tabata[7];

    static const double csCoefSingle_Hminus_Chebyshev[11];
    static const double csCoefDouble_Hminus_Chebyshev[11];
    static const double csCoefSingle_Hplus_Chebyshev[11];
    static const double csCoefDouble_Hplus_Chebyshev[11];
    static const double csCoefHydrogenProduction_H2plus_Chebyshev[11];

    static double a_m[9];
    static double b_m[3][9];
};

#endif //BEAMSTRIPPINGPHYSICS_HH