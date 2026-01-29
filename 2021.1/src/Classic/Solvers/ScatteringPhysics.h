//
// Class ScatteringPhysics
//   Defines the physical processes of beam scattering
//   and energy loss by heavy charged particles
//
// Copyright (c) 2009 - 2021, Bi, Yang, Stachel, Adelmann
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef SCATTERINGPHYSICS_HH
#define SCATTERINGPHYSICS_HH

#include "Solvers/ParticleMatterInteractionHandler.h"

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"

#include <gsl/gsl_rng.h>

#include "Utility/IpplTimings.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;
class Inform;

typedef struct {             // struct for description of particle in material
    int label;               // the status of the particle (0 = in material / -1 = move back to bunch
    unsigned localID;        // not so unique identifier of the particle
    Vector_t Rincol;         // position
    Vector_t Pincol;         // momentum
    long IDincol;            // unique identifier of the particle inherited from the bunch
    int Binincol;            // bin number
    double DTincol;          // time step size
    double Qincol;           // charge
    double Mincol;           // mass
    Vector_t Bfincol;        // magnetic field
    Vector_t Efincol;        // electric field
} PART;

struct InsideTester {
    virtual ~InsideTester()
    { }

    virtual
    bool checkHit(const Vector_t &R) = 0;
};


class ScatteringPhysics: public ParticleMatterInteractionHandler {
public:
    ScatteringPhysics(const std::string& name,
                      ElementBase* element,
                      std::string& mat,
                      bool enableRutherford,
                      double lowEnergyThr);
    ~ScatteringPhysics();

    virtual void apply(PartBunchBase<double, 3>* bunch,
                       const std::pair<Vector_t, double>& boundingSphere);

    virtual const std::string getType() const;
    virtual void print(Inform& os);
    virtual bool stillActive();

    virtual double getTime();
    virtual std::string getName();
    virtual size_t getParticlesInMat();
    virtual unsigned getRediffused();
    virtual unsigned int getNumEntered();

    void computeInteraction(PartBunchBase<double, 3>* bunch);

    virtual bool computeEnergyLoss(PartBunchBase<double, 3>* bunch,
                                   Vector_t& P,
                                   const double deltat,
                                   bool includeFluctuations = true) const;

private:

    void configureMaterialParameters();
    void computeCoulombScattering(Vector_t& R,
                                  Vector_t& P,
                                  double dt);

    void applyRotation(Vector_t& P,
                       Vector_t& R,
                       double xplane,
                       double thetacou);
    void applyRandomRotation(Vector_t& P, double theta0);

    void copyFromBunch(PartBunchBase<double, 3>* bunch,
                       const std::pair<Vector_t, double>& boundingSphere);
    void addBackToBunch(PartBunchBase<double, 3>* bunch);

    void deleteParticleFromLocalVector();

    void calcStat(double Eng);
    void gatherStatistics();

    void push();
    void resetTimeStep();
    void setTimeStepForLeavingParticles();

    double  T_m;                               // own time, maybe larger than in the bunch object
    double dT_m;                               // dt from bunch

    double mass_m;                             // mass from bunch (eV)
    double charge_m;                           // charge from bunch (elementary charges)

    gsl_rng* rGen_m;                           // random number generator
    std::string material_m;                    // type of material e.g. aluminum
    std::unique_ptr<InsideTester> hitTester_m; // tests whether particles are inside material
    ElementBase::ElementType collshape_m;      // the type of element (DEGRADER, CCOLLIMATOR or FLEXIBLECOLLIMATOR)
    std::string collshapeStr_m;                // the type of element as string

    // material parameters
    double Z_m;                                // the atomic number [1]
    double A_m;                                // the atomic mass [u]
    double rho_m;                              // the volumetric mass density in [g cm^-3]
    double X0_m;                               // the radiation length in [m]
    double I_m;                                // the mean excitation energy [eV]

    /*
       coefficients to fit model to measurement data according to Andersen-Ziegler formulae.
       see ICRU-49, "Stopping Powers and Ranges for Protons  and Alpha Particles",
       chapter 'Electronic (Collision) Stopping Powers in the Low-Energy Region'
    */
    double A1_c;
    double A2_c;
    double A3_c;
    double A4_c;
    double A5_c;
    double B1_c;
    double B2_c;
    double B3_c;
    double B4_c;
    double B5_c;

    // number of particles that enter the material in current step (count for single step)
    unsigned int bunchToMatStat_m;
    // number of particles that are stopped by the material in current step
    unsigned int stoppedPartStat_m;
    // number of particles that leave the material in current step
    unsigned int rediffusedStat_m;
    // total number of particles that are in the material
    unsigned int totalPartsInMat_m;

    // some statistics
    double Eavg_m;                            // average kinetic energy
    double Emax_m;                            // maximum kinetic energy
    double Emin_m;                            // minimum kinetic energy

    std::vector<PART> locParts_m;             // local particles that are in material

    std::unique_ptr<LossDataSink> lossDs_m;

    bool enableRutherford_m;
    double lowEnergyThr_m;

    IpplTimings::TimerRef DegraderApplyTimer_m;
    IpplTimings::TimerRef DegraderLoopTimer_m;
    IpplTimings::TimerRef DegraderDestroyTimer_m;
};

inline
void ScatteringPhysics::calcStat(double Eng) {
    Eavg_m += Eng;
    if (Emin_m > Eng)
        Emin_m = Eng;
    if (Emax_m < Eng)
        Emax_m = Eng;
}

inline
double ScatteringPhysics::getTime() {
    return T_m;
}

inline
std::string ScatteringPhysics::getName() {
    return (element_ref_m->getName() + "_" + name_m);
}

inline
size_t ScatteringPhysics::getParticlesInMat() {
    return totalPartsInMat_m;
}

inline
unsigned int ScatteringPhysics::getRediffused() {
    return rediffusedStat_m;
}

inline
unsigned int ScatteringPhysics::getNumEntered() {
    return bunchToMatStat_m;
}

inline
const std::string ScatteringPhysics::getType() const {
    return "ScatteringPhysics";
}

#endif //SCATTERINGPHYSICS_HH