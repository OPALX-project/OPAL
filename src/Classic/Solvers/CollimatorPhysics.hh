//
// Class CollimatorPhysics
//
// Defines the collimator physics models
//
// Copyright (c) 2009 - 2020, Bi, Yang, Stachel, Adelmann
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

#ifndef COLLIMATORPHYSICS_HH
#define COLLIMATORPHYSICS_HH

#include "Solvers/ParticleMatterInteractionHandler.hh"

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Vektor.h"
#include <gsl/gsl_rng.h>

#include "Utility/IpplTimings.h"

#include <memory>
#include <utility>
#include <string>
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
    Vector_t Bfincol;        // magnetic field
    Vector_t Efincol;        // electric field
} PART;

struct InsideTester {
    virtual ~InsideTester()
    { }

    virtual
    bool checkHit(const Vector_t &R) = 0;
};


class CollimatorPhysics: public ParticleMatterInteractionHandler {
public:
    CollimatorPhysics(const std::string &name,
                      ElementBase *element,
                      std::string &mat,
                      bool enableRutherford);
    ~CollimatorPhysics();

    virtual void apply(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere);

    virtual const std::string getType() const;

    virtual void print(Inform& os);
    virtual bool stillActive();

    virtual double getTime();
    virtual std::string getName();
    virtual size_t getParticlesInMat();
    virtual unsigned getRediffused();
    virtual unsigned int getNumEntered();
    void computeInteraction();

    virtual bool computeEnergyLoss(Vector_t &P,
                                   const double deltat,
                                   bool includeFluctuations = true) const;

private:

    void configureMaterialParameters();
    void computeCoulombScattering(Vector_t &R,
                                  Vector_t &P,
                                  double dt);

    void applyRotation(Vector_t &P,
                       Vector_t &R,
                       double xplane,
                       double thetacou);
    void applyRandomRotation(Vector_t &P, double theta0);

    void copyFromBunch(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere);
    void addBackToBunch(PartBunchBase<double, 3> *bunch);

    void deleteParticleFromLocalVector();

    void calcStat(double Eng);
    void gatherStatistics();

    void push();
    void resetTimeStep();
    void setTimeStepForLeavingParticles();

    double  T_m;                               // own time, maybe larger than in the bunch object
    double dT_m;                               // dt from bunch

    gsl_rng *rGen_m;                           // random number generator

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

    double A2_c;                               // coefficients to fit model to measurement data
    double A3_c;                               // see e.g. page 16-20 in H. H. Andersen, J. F. Ziegler,
    double A4_c;                               // "Hydrogen Stopping Powers and Ranges in All Elements",
    double A5_c;                               // Pergamon Press, 1977

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

    IpplTimings::TimerRef DegraderApplyTimer_m;
    IpplTimings::TimerRef DegraderLoopTimer_m;
    IpplTimings::TimerRef DegraderDestroyTimer_m;
};

inline
void CollimatorPhysics::calcStat(double Eng) {
    Eavg_m += Eng;
    if (Emin_m > Eng)
        Emin_m = Eng;
    if (Emax_m < Eng)
        Emax_m = Eng;
}

inline
double CollimatorPhysics::getTime() {
    return T_m;
}

inline
size_t CollimatorPhysics::getParticlesInMat() {
    return totalPartsInMat_m;
}

inline
unsigned int CollimatorPhysics::getRediffused() {
    return rediffusedStat_m;
}

inline
unsigned int CollimatorPhysics::getNumEntered() {
    return bunchToMatStat_m;
}

inline
const std::string CollimatorPhysics::getType() const {
    return "CollimatorPhysics";
}

#endif //COLLIMATORPHYSICS_HH
