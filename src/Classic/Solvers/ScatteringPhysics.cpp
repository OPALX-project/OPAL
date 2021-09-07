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
#include "Solvers/ScatteringPhysics.h"

#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Multipole.h"
#include "Physics/Physics.h"
#include "Physics/Material.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"

#include "Utility/Inform.h"

#include <gsl/gsl_randist.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <sys/time.h>

namespace {
    struct DegraderInsideTester: public InsideTester {
        explicit DegraderInsideTester(ElementBase* el) {
            deg_m = static_cast<Degrader*>(el);
        }
        virtual
        bool checkHit(const Vector_t& R) override {
            return deg_m->isInMaterial(R(2));
        }
    private:
        Degrader* deg_m;
    };

    struct CollimatorInsideTester: public InsideTester {
        explicit CollimatorInsideTester(ElementBase* el) {
            col_m = static_cast<CCollimator*>(el);
        }
        virtual
        bool checkHit(const Vector_t& R)  override {
            return col_m->checkPoint(R(0), R(1));
        }
    private:
        CCollimator* col_m;
    };

    struct FlexCollimatorInsideTester: public InsideTester {
        explicit FlexCollimatorInsideTester(ElementBase* el) {
            col_m = static_cast<FlexibleCollimator*>(el);
        }
        virtual
        bool checkHit(const Vector_t& R)  override {
            return col_m->isStopped(R);
        }
    private:
        FlexibleCollimator* col_m;
    };
}

ScatteringPhysics::ScatteringPhysics(const std::string& name,
                                     ElementBase* element,
                                     std::string& material,
                                     bool enableRutherford,
                                     double lowEnergyThr):
    ParticleMatterInteractionHandler(name, element),
    T_m(0.0),
    dT_m(0.0),
    mass_m(0.0),
    charge_m(0.0),
    material_m(material),
    hitTester_m(nullptr),
    Z_m(0),
    A_m(0.0),
    rho_m(0.0),
    X0_m(0.0),
    I_m(0.0),
    A1_c(0.0),
    A2_c(0.0),
    A3_c(0.0),
    A4_c(0.0),
    A5_c(0.0),
    B1_c(0.0),
    B2_c(0.0),
    B3_c(0.0),
    B4_c(0.0),
    B5_c(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    totalPartsInMat_m(0),
    Eavg_m(0.0),
    Emax_m(0.0),
    Emin_m(0.0),
    enableRutherford_m(enableRutherford),
    lowEnergyThr_m(lowEnergyThr)
{

    gsl_rng_env_setup();
    rGen_m = gsl_rng_alloc(gsl_rng_default);

    size_t mySeed = Options::seed;

    if (Options::seed == -1) {
        struct timeval tv;
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec + tv.tv_usec;
    }

    gsl_rng_set(rGen_m, mySeed + Ippl::myNode());

    configureMaterialParameters();

    collshape_m = element_ref_m->getType();
    switch (collshape_m) {
    case ElementBase::DEGRADER:
        hitTester_m.reset(new DegraderInsideTester(element_ref_m));
        break;
    case ElementBase::CCOLLIMATOR:
        hitTester_m.reset(new CollimatorInsideTester(element_ref_m));
        break;
    case ElementBase::FLEXIBLECOLLIMATOR:
        hitTester_m.reset(new FlexCollimatorInsideTester(element_ref_m));
        break;
    default:
        throw GeneralClassicException("ScatteringPhysics::ScatteringPhysics",
                                      "Unsupported element type");
    }

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump && Options::enableHDF5));

    DegraderApplyTimer_m   = IpplTimings::getTimer("DegraderApply");
    DegraderLoopTimer_m    = IpplTimings::getTimer("DegraderLoop");
    DegraderDestroyTimer_m = IpplTimings::getTimer("DegraderDestroy");
}

ScatteringPhysics::~ScatteringPhysics() {
    locParts_m.clear();
    lossDs_m->save();
    if (rGen_m)
        gsl_rng_free(rGen_m);
}


/// The material of the collimator
//  ------------------------------------------------------------------------
void  ScatteringPhysics::configureMaterialParameters() {

    auto material = Physics::Material::getMaterial(material_m);
    Z_m = material->getAtomicNumber();
    A_m = material->getAtomicMass();
    rho_m = material->getMassDensity();
    X0_m = material->getRadiationLength();
    I_m = material->getMeanExcitationEnergy();
    A1_c = material->getStoppingPowerFitCoefficients(Physics::Material::A1);
    A2_c = material->getStoppingPowerFitCoefficients(Physics::Material::A2);
    A3_c = material->getStoppingPowerFitCoefficients(Physics::Material::A3);
    A4_c = material->getStoppingPowerFitCoefficients(Physics::Material::A4);
    A5_c = material->getStoppingPowerFitCoefficients(Physics::Material::A5);
    B1_c = material->getStoppingPowerFitCoefficients(Physics::Material::B1);
    B2_c = material->getStoppingPowerFitCoefficients(Physics::Material::B2);
    B3_c = material->getStoppingPowerFitCoefficients(Physics::Material::B3);
    B4_c = material->getStoppingPowerFitCoefficients(Physics::Material::B4);
    B5_c = material->getStoppingPowerFitCoefficients(Physics::Material::B5);
}

void ScatteringPhysics::apply(PartBunchBase<double, 3>* bunch,
                              const std::pair<Vector_t, double>& boundingSphere) {
    IpplTimings::startTimer(DegraderApplyTimer_m);

    /*
      Particles that have entered material are flagged as Bin[i] == -1.

      Flagged particles are copied to a local structure within Scattering Physics locParts_m.

      Particles in that structure will be pushed in the material and either come
      back to the bunch or will be fully stopped in the material.
    */

    ParticleType pType = bunch->getPType();
    if (pType != ParticleType::PROTON   &&
        pType != ParticleType::DEUTERON &&
        pType != ParticleType::HMINUS   &&
        pType != ParticleType::MUON     &&
        pType != ParticleType::H2P      &&
        pType != ParticleType::ALPHA) {

        throw GeneralClassicException(
                "ScatteringPhysics::apply",
                "Particle " + bunch->getPTypeString() +
                " is not supported for scattering interactions!");
    }

    Eavg_m = 0.0;
    Emax_m = 0.0;
    Emin_m = 0.0;

    bunchToMatStat_m  = 0;
    rediffusedStat_m   = 0;
    stoppedPartStat_m = 0;
    totalPartsInMat_m   = 0;

    dT_m = bunch->getdT();
    T_m  = bunch->getT();
    mass_m   = bunch->getM();
    charge_m = bunch->getQ();

    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    do {
        IpplTimings::startTimer(DegraderLoopTimer_m);
        push();
        setTimeStepForLeavingParticles();

        // if we are not looping copy newly arrived particles
        if (onlyOneLoopOverParticles) {
            copyFromBunch(bunch, boundingSphere);
        }
        addBackToBunch(bunch);

        computeInteraction(bunch);

        push();
        resetTimeStep();

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        T_m += dT_m; // update local time

        gatherStatistics();

        if (allParticleInMat_m) {
            onlyOneLoopOverParticles = (rediffusedStat_m > 0 || totalPartsInMat_m <= 0);
        } else {
            onlyOneLoopOverParticles = true;
        }
    } while (onlyOneLoopOverParticles == false);

    IpplTimings::stopTimer(DegraderApplyTimer_m);
}

void ScatteringPhysics::computeInteraction(PartBunchBase<double, 3>* bunch) {
    /*
        Do physics if
        -- correct type of particle
        -- particle not stopped (locParts_m[i].label != -1.0)

        Absorbed particle i: locParts_m[i].label = -1.0;
    */
    for (size_t i = 0; i < locParts_m.size(); ++i) {
        if (locParts_m[i].label != -1) {
            Vector_t &R = locParts_m[i].Rincol;
            Vector_t &P = locParts_m[i].Pincol;
            double &dt  = locParts_m[i].DTincol;

            if (hitTester_m->checkHit(R)) {
                bool pdead = computeEnergyLoss(bunch, P, dt);
                if (!pdead) {
                    /*
                      Now scatter and transport particle in material.
                      The checkHit call just above will detect if the
                      particle is rediffused from the material into vacuum.
                    */

                    computeCoulombScattering(R, P, dt);
                } else {
                    // The particle is stopped in the material, set label to -1
                    locParts_m[i].label = -1.0;
                    ++ stoppedPartStat_m;
                    lossDs_m->addParticle(OpalParticle(locParts_m[i].IDincol,
                                                       R, P, T_m,
                                                       locParts_m[i].Qincol, locParts_m[i].Mincol));
                }
            }
        }
    }

    // delete absorbed particles
    deleteParticleFromLocalVector();
}

/// Energy Loss: using the Bethe-Bloch equation.
/// In low-energy region use Andersen-Ziegler fitting (only for protons and alpha)
/// Energy straggling: For relatively thick absorbers such that the number of collisions
/// is large, the energy loss distribution is shown to be Gaussian in form.
/// See Particle Physics Booklet, chapter 'Passage of particles through matter' or
/// Review of Particle Physics, DOI: 10.1103/PhysRevD.86.010001, page 329 ff
// -------------------------------------------------------------------------
bool ScatteringPhysics::computeEnergyLoss(PartBunchBase<double, 3>* bunch,
                                          Vector_t& P,
                                          const double deltat,
                                          bool includeFluctuations) const {

    ParticleType pType = bunch->getPType();
    constexpr double m2cm = 1e2;
    constexpr double eV2keV = 1e-3;
    constexpr double GeV2keV = 1e6;
    const double mass_keV = mass_m * eV2keV;
    constexpr double massElectron_keV = Physics::m_e * GeV2keV;

    constexpr double K = 4.0 * Physics::pi * Physics::Avo * Physics::r_e * m2cm * Physics::r_e * m2cm * massElectron_keV;

    double gamma = Util::getGamma(P);
    const double gammaSqr = std::pow(gamma, 2);
    const double betaSqr = 1.0 - 1.0 / gammaSqr;
    double beta = std::sqrt(betaSqr);
    double Ekin = (gamma - 1) * mass_keV;
    double dEdx = 0.0;
    double epsilon = 0.0;

    const double deltas = deltat * beta * Physics::c;
    const double deltasrho = deltas * m2cm * rho_m;

    const double massRatio = massElectron_keV / mass_keV;
    double Tmax = (2.0 * massElectron_keV * betaSqr * gammaSqr /
                  (std::pow(gamma + massRatio, 2) - (gammaSqr - 1.0)));

    if (pType != ParticleType::ALPHA) {

        if (Ekin >= 600 && Ekin < 1e7) {
            dEdx = (-K * std::pow(charge_m, 2) * Z_m / (A_m * betaSqr) *
                    (0.5 * std::log(2 * massElectron_keV * betaSqr * gammaSqr * Tmax / std::pow(I_m * eV2keV, 2)) - betaSqr));
        } else if (pType == ParticleType::PROTON && Ekin < 600) {
            constexpr double massProton_amu = Physics::m_p / Physics::amu;
            const double Ts = Ekin / massProton_amu;
            if (Ekin > 10) {
                const double epsilon_low = A2_c * std::pow(Ts, 0.45);
                const double epsilon_high = (A3_c / Ts) * std::log(1 + (A4_c / Ts) + (A5_c * Ts));
                epsilon = (epsilon_low * epsilon_high) / (epsilon_low + epsilon_high);
            } else if (Ekin > 1) {
                epsilon = A1_c * std::pow(Ts, 0.5);
            }
            dEdx = -epsilon / (1e18 * (A_m / Physics::Avo));
        } else {
            INFOMSG(level4 << "Particle energy out of the valid range "
                              "for energy loss calculation." << endl);
        }
    } else {
        if (Ekin > 10000 && Ekin < 1e6) {
            dEdx = (-K * std::pow(charge_m, 2) * Z_m / (A_m * betaSqr) *
                    (0.5 * std::log(2 * massElectron_keV * betaSqr * gammaSqr * Tmax / std::pow(I_m * eV2keV, 2)) - betaSqr));
        } else if (Ekin > 1 && Ekin <= 10000) {
            const double T = Ekin * 1e-3; //MeV
            const double epsilon_low = B1_c * std::pow(1000*T, B2_c);
            const double epsilon_high = (B3_c / T) * std::log(1 + (B4_c / T) + (B5_c * T));
            epsilon = (epsilon_low * epsilon_high) / (epsilon_low + epsilon_high);
            dEdx = -epsilon / (1e18 * (A_m / Physics::Avo));
        } else {
            INFOMSG(level4 << "Particle energy out of the valid range "
                              "for energy loss calculation." << endl);
        }
    }

    Ekin += deltasrho * dEdx;

    if (includeFluctuations) {
        double sigma_E = std::sqrt(K * massElectron_keV * rho_m * (Z_m / A_m) * deltas * m2cm);
        Ekin += gsl_ran_gaussian(rGen_m, sigma_E);
    }

    gamma = Ekin / mass_keV + 1.0;
    beta = std::sqrt(1.0 - 1.0 / std::pow(gamma, 2));
    P = gamma * beta * P / euclidean_norm(P);

    bool stopped = (Ekin < lowEnergyThr_m * 1e3 || dEdx > 0);
    return stopped;
}


// Implement the rotation in 2 dimensions here
// For details see: J. Beringer et al. (Particle Data Group),
// "Passage of particles through matter", Phys. Rev. D 86, 010001 (2012)
// -------------------------------------------------------------------------
void  ScatteringPhysics::applyRotation(Vector_t& P,
                                       Vector_t& R,
                                       double shift,
                                       double thetacou) {
    // Calculate the angle between the transverse and longitudinal component of the momentum
    double Psixz = std::fmod(std::atan2(P(0), P(2)) + Physics::two_pi, Physics::two_pi);

    R(0) = R(0) + shift * std::cos(Psixz);
    R(2) = R(2) - shift * std::sin(Psixz);

    // Apply the rotation about the random angle thetacou
    double Px = P(0);
    P(0) =  Px * std::cos(thetacou) + P(2) * std::sin(thetacou);
    P(2) = -Px * std::sin(thetacou) + P(2) * std::cos(thetacou);
}

void ScatteringPhysics::applyRandomRotation(Vector_t& P, double theta0) {

    double thetaru = 2.5 / std::sqrt(gsl_rng_uniform(rGen_m)) * 2.0 * theta0;
    double phiru = Physics::two_pi * gsl_rng_uniform(rGen_m);

    double normPtrans = std::sqrt(P(0) * P(0) + P(1) * P(1));
    double Theta = std::atan(normPtrans / std::abs(P(2)));
    double CosT = std::cos(Theta);
    double SinT = std::sin(Theta);

    Vector_t X(std::cos(phiru)*std::sin(thetaru),
               std::sin(phiru)*std::sin(thetaru),
               std::cos(thetaru));
    X *= euclidean_norm(P);

    Vector_t W(-P(1), P(0), 0.0);
    W = W / normPtrans;

    // Rodrigues' formula for a rotation about W by angle Theta
    P = X * CosT + cross(W, X) * SinT + W * dot(W, X) * (1.0 - CosT);
}

/// Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
/// Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
//--------------------------------------------------------------------------
void  ScatteringPhysics::computeCoulombScattering(Vector_t& R,
                                                  Vector_t& P,
                                                  double dt) {

    constexpr double sqrtThreeInv = 0.57735026918962576451; // sqrt(1.0 / 3.0)
    const double normP = euclidean_norm(P);
    const double gammaSqr = std::pow(normP, 2) + 1.0;
    const double beta = std::sqrt(1.0 - 1.0 / gammaSqr);
    const double deltas = dt * beta * Physics::c;
    const double theta0 = (13.6e6 / (beta * normP * mass_m) *
                           charge_m * std::sqrt(deltas / X0_m) *
                           (1.0 + 0.038 * std::log(deltas / X0_m)));

    double phi = Physics::two_pi * gsl_rng_uniform(rGen_m);
    for (unsigned int i = 0; i < 2; ++ i) {
        CoordinateSystemTrafo randomTrafo(R, Quaternion(cos(phi), 0, 0, sin(phi)));
        P = randomTrafo.rotateTo(P);
        R = Vector_t(0.0); // corresponds to randomTrafo.transformTo(R);

        double z1 = gsl_ran_gaussian(rGen_m, 1.0);
        double z2 = gsl_ran_gaussian(rGen_m, 1.0);

        while(std::abs(z2) > 3.5) {
            z1 = gsl_ran_gaussian(rGen_m, 1.0);
            z2 = gsl_ran_gaussian(rGen_m, 1.0);
        }

        double thetacou = z2 * theta0;
        double shift = (z1 * sqrtThreeInv + z2) * deltas * theta0 * 0.5;
        applyRotation(P, R, shift, thetacou);

        P = randomTrafo.rotateFrom(P);
        R = randomTrafo.transformFrom(R);

        phi += 0.5 * Physics::pi;
    }

    if (enableRutherford_m &&
        gsl_rng_uniform(rGen_m) < 0.0047) {

        applyRandomRotation(P, theta0);
    }
}

void ScatteringPhysics::addBackToBunch(PartBunchBase<double, 3>* bunch) {

    const size_t nL = locParts_m.size();
    if (nL == 0) return;

    unsigned int numLocalParticles = bunch->getLocalNum();
    const double elementLength = element_ref_m->getElementLength();

    for (size_t i = 0; i < nL; ++ i) {
        Vector_t& R = locParts_m[i].Rincol;

        if (R[2] >= elementLength) {

            bunch->createWithID(locParts_m[i].IDincol);

            /*
              Binincol is still <0, but now the particle is rediffused
              from the material and hence this is not a "lost" particle anymore
            */
            bunch->Bin[numLocalParticles] = 1;
            bunch->R[numLocalParticles]   = R;
            bunch->P[numLocalParticles]   = locParts_m[i].Pincol;
            bunch->Q[numLocalParticles]   = locParts_m[i].Qincol;
            bunch->M[numLocalParticles]   = locParts_m[i].Mincol;
            bunch->Bf[numLocalParticles]  = 0.0;
            bunch->Ef[numLocalParticles]  = 0.0;
            bunch->dt[numLocalParticles]  = dT_m;

            /*
              This particle is back to the bunch, by set
              ting the label to -1.0
              the particle will be deleted.
            */
            locParts_m[i].label = -1.0;

            ++ rediffusedStat_m;
            ++ numLocalParticles;
        }
    }

    // delete particles that went to the bunch
    deleteParticleFromLocalVector();
}

void ScatteringPhysics::copyFromBunch(PartBunchBase<double, 3>* bunch,
                                      const std::pair<Vector_t, double>& boundingSphere)
{
    const size_t nL = bunch->getLocalNum();
    if (nL == 0) return;

    IpplTimings::startTimer(DegraderDestroyTimer_m);
    double zmax = boundingSphere.first(2) + boundingSphere.second;
    double zmin = boundingSphere.first(2) - boundingSphere.second;
    if (zmax < 0.0 || zmin > element_ref_m->getElementLength()) {
        IpplTimings::stopTimer(DegraderDestroyTimer_m);
        return;
    }

    size_t ne = 0;
    std::set<size_t> partsToDel;
    for (size_t i = 0; i < nL; ++ i) {
        if ((bunch->Bin[i] == -1 || bunch->Bin[i] == 1) &&
            hitTester_m->checkHit(bunch->R[i]))
        {
            // adjust the time step for those particles that enter the material
            // such that it corresponds to the time needed to reach the current
            // location form the edge of the material. Only use this time step
            // for the computation of the interaction with the material, not for
            // the integration of the particles. This will ensure that the momenta
            // of all particles are reduced by approximately the same amount in
            // computeEnergyLoss.
            double betaz = bunch->P[i](2) / Util::getGamma(bunch->P[i]);
            double stepWidth = betaz * Physics::c * bunch->dt[i];
            double tau = bunch->R[i](2) / stepWidth;

            PAssert_LE(tau, 1.0);
            PAssert_GE(tau, 0.0);

            PART x;
            x.localID      = i;
            x.DTincol      = bunch->dt[i] * tau;
            x.IDincol      = bunch->ID[i];
            x.Binincol     = bunch->Bin[i];
            x.Rincol       = bunch->R[i];
            x.Pincol       = bunch->P[i];
            x.Qincol       = bunch->Q[i];
            x.Mincol       = bunch->M[i];
            x.Bfincol      = bunch->Bf[i];
            x.Efincol      = bunch->Ef[i];
            x.label        = 0;            // alive in matter

            locParts_m.push_back(x);
            ne++;
            bunchToMatStat_m++;

            partsToDel.insert(i);
        }
    }

    for (auto it = partsToDel.begin(); it != partsToDel.end(); ++ it) {
        bunch->destroy(1, *it);
    }

    if (ne > 0) {
        bunch->performDestroy(true);
    }

    IpplTimings::stopTimer(DegraderDestroyTimer_m);
}

void ScatteringPhysics::print(Inform &msg) {
    Inform::FmtFlags_t ff = msg.flags();

    if (totalPartsInMat_m > 0 ||
        bunchToMatStat_m  > 0 ||
        rediffusedStat_m  > 0 ||
        stoppedPartStat_m > 0) {

        OPALTimer::Timer time;
        msg << level2
            << "--- ScatteringPhysics\n"
            << "Name: " << name_m << " - "
            << "Material: " << material_m << " - "
            << "Element: " << element_ref_m->getName() << "\n"
            << "Particle Statistics @ " << time.time() << "\n"
            << std::setw(21) << "entered: " << Util::toStringWithThousandSep(bunchToMatStat_m) << "\n"
            << std::setw(21) << "rediffused: " << Util::toStringWithThousandSep(rediffusedStat_m) << "\n"
            << std::setw(21) << "stopped: " << Util::toStringWithThousandSep(stoppedPartStat_m) << "\n"
            << std::setw(21) << "total in material: " << Util::toStringWithThousandSep(totalPartsInMat_m)
            << endl;
    }

    msg.flags(ff);
}

bool ScatteringPhysics::stillActive() {
    return totalPartsInMat_m != 0;
}

namespace {
    bool myCompF(PART x, PART y) {
      return x.label > y.label;
    }
}

void ScatteringPhysics::deleteParticleFromLocalVector() {
    /*
      the particle to be deleted (label < 0) are all
      at the end of the vector.
    */
    sort(locParts_m.begin(), locParts_m.end(), myCompF);

    // find start of particles to delete
    std::vector<PART>::iterator inv = locParts_m.begin();

    for (; inv != locParts_m.end(); ++inv) {
        if ((*inv).label == -1)
            break;
    }
    locParts_m.erase(inv, locParts_m.end());
    locParts_m.resize(inv - locParts_m.begin());

    // update statistics
    if (locParts_m.size() > 0) {
        Eavg_m /= locParts_m.size();
        Emin_m /= locParts_m.size();
        Emax_m /= locParts_m.size();
    }
}

void ScatteringPhysics::push() {
    for (size_t i = 0; i < locParts_m.size(); ++i) {
        Vector_t& R  = locParts_m[i].Rincol;
        Vector_t& P  = locParts_m[i].Pincol;
        double gamma = Util::getGamma(P);

        R += 0.5 * dT_m * Physics::c * P / gamma;
    }
}

// adjust the time step for those particles that will rediffuse to
// vacuum such that it corresponds to the time needed to reach the
// end of the material. Only use this time step for the computation
// of the interaction with the material, not for the integration of
// the particles. This will ensure that the momenta of all particles
// are reduced by  approcimately the same amount in computeEnergyLoss.
void ScatteringPhysics::setTimeStepForLeavingParticles() {
    const double elementLength = element_ref_m->getElementLength();

    for (size_t i = 0; i < locParts_m.size(); ++i) {
        Vector_t& R  = locParts_m[i].Rincol;
        Vector_t& P  = locParts_m[i].Pincol;
        double& dt   = locParts_m[i].DTincol;
        double gamma = Util::getGamma(P);
        Vector_t stepLength = dT_m * Physics::c * P / gamma;

        if (R(2) < elementLength &&
            R(2) + stepLength(2) > elementLength) {
            // particle is likely to leave material
            double distance = elementLength - R(2);
            double tau = distance / stepLength(2);

            PAssert_LE(tau, 1.0);
            PAssert_GE(tau, 0.0);

            dt *= tau;
        }
    }
}

void ScatteringPhysics::resetTimeStep() {
    for (size_t i = 0; i < locParts_m.size(); ++i) {
        double& dt = locParts_m[i].DTincol;
        dt = dT_m;
    }
}


void ScatteringPhysics::gatherStatistics() {

    unsigned int locPartsInMat;
    locPartsInMat = locParts_m.size();

    constexpr unsigned short numItems = 4;
    unsigned int partStatistics[numItems] = {locPartsInMat,
                                             bunchToMatStat_m,
                                             rediffusedStat_m,
                                             stoppedPartStat_m};

    allreduce(&(partStatistics[0]), numItems, std::plus<unsigned int>());

    totalPartsInMat_m = partStatistics[0];
    bunchToMatStat_m = partStatistics[1];
    rediffusedStat_m = partStatistics[2];
    stoppedPartStat_m = partStatistics[3];
}