// Class:CollimatorPhysics
// Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang Stachel, Adelmann$
//-------------------------------------------------------------------------
#include "Solvers/CollimatorPhysics.hh"
#include "Physics/Physics.h"
#include "Physics/Material.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Multipole.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"

#include "Ippl.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

#ifdef OPAL_DKS
const int CollimatorPhysics::numpar_m = 13;
#endif

namespace {
    struct DegraderInsideTester: public InsideTester {
        explicit DegraderInsideTester(ElementBase * el) {
            deg_m = static_cast<Degrader*>(el);
        }
        virtual
        bool checkHit(const Vector_t &R, const Vector_t &P, double dt) override {
            return deg_m->isInMaterial(R(2));
        }

    private:
        Degrader *deg_m;
    };

    struct CollimatorInsideTester: public InsideTester {
        explicit CollimatorInsideTester(ElementBase * el) {
            col_m = static_cast<CCollimator*>(el);
        }
        virtual
        bool checkHit(const Vector_t &R, const Vector_t &P, double dt)  override {
            return col_m->checkPoint(R(0), R(1));
        }

    private:
        CCollimator *col_m;
    };

    struct FlexCollimatorInsideTester: public InsideTester {
        explicit FlexCollimatorInsideTester(ElementBase * el) {
            col_m = static_cast<FlexibleCollimator*>(el);
        }
        virtual
        bool checkHit(const Vector_t &R, const Vector_t &P, double dt)  override {
            return col_m->isStopped(R, P, Physics::c * dt / sqrt(1.0  + dot(P, P)));
        }

    private:
        FlexibleCollimator *col_m;
    };
}

CollimatorPhysics::CollimatorPhysics(const std::string &name,
                                     ElementBase *element,
                                     std::string &material,
                                     bool enableRutherford):
    ParticleMatterInteractionHandler(name, element),
    T_m(0.0),
    dT_m(0.0),
    material_m(material),
    hitTester_m(nullptr),
    Z_m(0),
    A_m(0.0),
    rho_m(0.0),
    X0_m(0.0),
    I_m(0.0),
    A2_c(0.0),
    A3_c(0.0),
    A4_c(0.0),
    A5_c(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    totalPartsInMat_m(0),
    Eavg_m(0.0),
    Emax_m(0.0),
    Emin_m(0.0),
    enableRutherford_m(enableRutherford)
#ifdef OPAL_DKS
    , curandInitSet_m(0)
    , ierr_m(0)
    , maxparticles_m(0)
    , numparticles_m(0)
    , numlocalparts_m(0)
    , par_ptr(NULL)
    , mem_ptr(NULL)
#endif
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
        throw OpalException("CollimatorPhysics::CollimatorPhysics",
                            "Unsupported element type");
    }

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        dksbase_m.setAPI("Cuda", 4);
        dksbase_m.setDevice("-gpu", 4);
        dksbase_m.initDevice();
        curandInitSet_m = -1;
    }

#endif

    DegraderApplyTimer_m = IpplTimings::getTimer("DegraderApply");
    DegraderLoopTimer_m = IpplTimings::getTimer("DegraderLoop");
    DegraderDestroyTimer_m = IpplTimings::getTimer("DegraderDestroy");
}

CollimatorPhysics::~CollimatorPhysics() {
    locParts_m.clear();
    lossDs_m->save();
    if (rGen_m)
        gsl_rng_free(rGen_m);

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        clearCollimatorDKS();
#endif

}

std::string CollimatorPhysics::getName() {
    return element_ref_m->getName();
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::configureMaterialParameters() {

    auto material = Physics::Material::getMaterial(material_m);
    Z_m = material->getAtomicNumber();
    A_m = material->getAtomicMass();
    rho_m = material->getMassDensity();
    X0_m = material->getRadiationLength();
    I_m = material->getMeanExcitationEnergy();
    A2_c = material->getStoppingPowerFitCoefficients(Physics::Material::A2);
    A3_c = material->getStoppingPowerFitCoefficients(Physics::Material::A3);
    A4_c = material->getStoppingPowerFitCoefficients(Physics::Material::A4);
    A5_c = material->getStoppingPowerFitCoefficients(Physics::Material::A5);
}

void CollimatorPhysics::apply(PartBunchBase<double, 3> *bunch,
                              const std::pair<Vector_t, double> &boundingSphere,
                              size_t numParticlesInSimulation) {
    IpplTimings::startTimer(DegraderApplyTimer_m);

    /*
      Particles that have entered material are flagged as Bin[i] == -1.

      Flagged particles are copied to a local structure within Collimator Physics locParts_m.

      Particles in that structure will be pushed in the material and either come
      back to the bunch or will be fully stopped in the material.
    */

    Eavg_m = 0.0;
    Emax_m = 0.0;
    Emin_m = 0.0;

    bunchToMatStat_m  = 0;
    rediffusedStat_m   = 0;
    stoppedPartStat_m = 0;
    totalPartsInMat_m   = 0;

    dT_m = bunch->getdT();
    T_m  = bunch->getT();

#ifdef OPAL_DKS
    if (collshape_m == ElementBase::DEGRADER && IpplInfo::DKSEnabled) {
        applyDKS(bunch, boundingSphere, numParticlesInSimulation);
    } else {
        applyNonDKS(bunch, boundingSphere, numParticlesInSimulation);
    }
#else
    applyNonDKS(bunch, boundingSphere, numParticlesInSimulation);
#endif

    IpplTimings::stopTimer(DegraderApplyTimer_m);
}

void CollimatorPhysics::applyNonDKS(PartBunchBase<double, 3> *bunch,
                                    const std::pair<Vector_t, double> &boundingSphere,
                                    size_t numParticlesInSimulation) {
    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    do {
        IpplTimings::startTimer(DegraderLoopTimer_m);
        push();
        setTimeStepForLeavingParticles();

        /*
          if we are not looping copy newly arrived particles
        */
        if (onlyOneLoopOverParticles)
            copyFromBunch(bunch, boundingSphere);

        addBackToBunch(bunch);

        computeInteraction();

        push();
        resetTimeStep();

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        T_m += dT_m;              // update local time

        gatherStatistics();

        if (allParticleInMat_m) {
            onlyOneLoopOverParticles = (rediffusedStat_m > 0 || totalPartsInMat_m <= 0);
        } else {
            onlyOneLoopOverParticles = true;
        }
    } while (onlyOneLoopOverParticles == false);
}

void CollimatorPhysics::computeInteraction() {
    /***
        Do physics if
        -- particle not stopped (locParts_m[i].label != -1.0)

        Absorbed particle i: locParts_m[i].label = -1.0;
    */

    for (size_t i = 0; i < locParts_m.size(); ++i) {
        if (locParts_m[i].label != -1) {
            Vector_t &R = locParts_m[i].Rincol;
            Vector_t &P = locParts_m[i].Pincol;
            double &dt  = locParts_m[i].DTincol;

            if (hitTester_m->checkHit(R, P, dt)) {
                bool pdead = computeEnergyLoss(P, dt);
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
                    lossDs_m->addParticle(R, P, -locParts_m[i].IDincol);
                }
            }
        }
    }

    /*
      delete absorbed particles
    */
    deleteParticleFromLocalVector();
}

/// Energy Loss:  using the Bethe-Bloch equation.
/// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
/// the energy loss distribution is shown to be Gaussian in form.
/// See Particle Physics Booklet, chapter 'Passage of particles through matter' or
/// Review of Particle Physics, DOI: 10.1103/PhysRevD.86.010001, page 329 ff
// -------------------------------------------------------------------------
bool CollimatorPhysics::computeEnergyLoss(Vector_t &P,
                                          const double deltat,
                                          bool includeFluctuations) const {

    constexpr double m2cm = 1e2;
    constexpr double GeV2keV = 1e6;
    constexpr double massElectron_keV = Physics::m_e * GeV2keV;
    constexpr double massProton_keV = Physics::m_p * GeV2keV;
    constexpr double massProton_amu = Physics::m_p / Physics::amu;
    constexpr double chargeProton = Physics::z_p;
    constexpr double K = 4.0 * Physics::pi * Physics::Avo * Physics::r_e * m2cm * Physics::r_e * m2cm * massElectron_keV;

    double gamma = Util::getGamma(P);
    const double gammaSqr = std::pow(gamma, 2);
    const double betaSqr = 1.0 - 1.0 / gammaSqr;
    double beta = sqrt(betaSqr);
    double Ekin = (gamma - 1) * massProton_keV;
    double dEdx = 0.0;

    const double deltas = deltat * beta * Physics::c;
    const double deltasrho = deltas * m2cm * rho_m;

    if (Ekin >= 600) {
        constexpr double eV2keV = 1e-3;
        constexpr double massRatio = Physics::m_e / Physics::m_p;
        double Tmax = (2.0 * massElectron_keV * betaSqr * gammaSqr /
                       (std::pow(gamma + massRatio, 2) - (gammaSqr - 1.0)));

        dEdx = (-K * std::pow(chargeProton, 2) * Z_m / (A_m * betaSqr) *
                (0.5 * std::log(2 * massElectron_keV * betaSqr * gammaSqr * Tmax / std::pow(I_m * eV2keV, 2)) - betaSqr));

        Ekin += deltasrho * dEdx;
    } else if (Ekin > 10) {
        const double Ts = Ekin / massProton_amu;
        const double epsilon_low = A2_c * pow(Ts, 0.45);
        const double epsilon_high = (A3_c / Ts) * log(1 + (A4_c / Ts) + (A5_c * Ts));
        const double epsilon = (epsilon_low * epsilon_high) / (epsilon_low + epsilon_high);
        dEdx = -epsilon / (1e18 * (A_m / Physics::Avo));

        Ekin += deltasrho * dEdx;
    }

    if (includeFluctuations) {
        double sigma_E = sqrt(K * massElectron_keV * rho_m * (Z_m / A_m) * deltas * m2cm);
        Ekin += gsl_ran_gaussian(rGen_m, sigma_E);
    }

    gamma = Ekin / massProton_keV + 1.0;
    beta = sqrt(1.0 - 1.0 / std::pow(gamma, 2));
    P = gamma * beta * P / euclidean_norm(P);

    bool stopped = (Ekin < 10 || dEdx > 0);
    return stopped;
}


// Implement the rotation in 2 dimensions here
// For details see:
//  J. Beringer et al. (Particle Data Group), "Passage of particles through matter"
//  Phys. Rev. D 86, 010001 (2012),
void  CollimatorPhysics::applyRotation(Vector_t &P,
                                       Vector_t &R,
                                       double shift,
                                       double thetacou) {
    // Calculate the angle between the transverse and longitudinal component of the momentum
    double Psixz = std::fmod(std::atan2(P(0), P(2)) + Physics::two_pi, Physics::two_pi);

    R(0) = R(0) + shift * cos(Psixz);
    R(2) = R(2) - shift * sin(Psixz);

    // Apply the rotation about the random angle thetacou
    double Px = P(0);
    P(0) =  Px * cos(thetacou) + P(2) * sin(thetacou);
    P(2) = -Px * sin(thetacou) + P(2) * cos(thetacou);
}

void CollimatorPhysics::applyRandomRotation(Vector_t &P, double theta0) {

    double thetaru = 2.5 / sqrt(gsl_rng_uniform(rGen_m)) * 2.0 * theta0;
    double phiru = Physics::two_pi * gsl_rng_uniform(rGen_m);

    double normPtrans = sqrt(P(0) * P(0) + P(1) * P(1));
    double Theta = std::atan(normPtrans / std::abs(P(2)));
    double CosT = cos(Theta);
    double SinT = sin(Theta);

    Vector_t X(cos(phiru)*sin(thetaru),
               sin(phiru)*sin(thetaru),
               cos(thetaru));
    X *= euclidean_norm(P);

    Vector_t W(-P(1), P(0), 0.0);
    W = W / normPtrans;

    // Rodrigues' formula for a rotation about W by angle Theta
    P = X * CosT + cross(W, X) * SinT + W * dot(W, X) * (1.0 - CosT);
}

/// Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
/// Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
//--------------------------------------------------------------------------
void  CollimatorPhysics::computeCoulombScattering(Vector_t &R,
                                                  Vector_t &P,
                                                  double dt) {

    constexpr double GeV2eV = 1e9;
    constexpr double massProton_eV = Physics::m_p * GeV2eV;
    constexpr double chargeProton = Physics::z_p;
    constexpr double sqrtThreeInv = 0.57735026918962576451; // sqrt(1.0 / 3.0)
    const double normP = euclidean_norm(P);
    const double gammaSqr = std::pow(normP, 2) + 1.0;
    const double beta = sqrt(1.0 - 1.0 / gammaSqr);
    const double deltas = dt * beta * Physics::c;
    const double theta0 = (13.6e6 / (beta * normP * massProton_eV) *
                           chargeProton * sqrt(deltas / X0_m) *
                           (1.0 + 0.038 * log(deltas / X0_m)));

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

void CollimatorPhysics::addBackToBunch(PartBunchBase<double, 3> *bunch) {

    const size_t nL = locParts_m.size();
    if (nL == 0) return;

    unsigned int numLocalParticles = bunch->getLocalNum();
    const double elementLength = element_ref_m->getElementLength();

    for (size_t i = 0; i < nL; ++ i) {
        Vector_t &R = locParts_m[i].Rincol;

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

    /*
      delete particles that went to the bunch
    */
    deleteParticleFromLocalVector();
}

void CollimatorPhysics::copyFromBunch(PartBunchBase<double, 3> *bunch,
                                      const std::pair<Vector_t, double> &boundingSphere)
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
    const unsigned int minNumOfParticlesPerCore = bunch->getMinimumNumberOfParticlesPerCore();
    for (size_t i = 0; i < nL; ++ i) {
        if ((bunch->Bin[i] == -1 || bunch->Bin[i] == 1) &&
            ((nL - ne) > minNumOfParticlesPerCore) &&
            hitTester_m->checkHit(bunch->R[i], bunch->P[i], dT_m))
        {
            // adjust the time step for those particles that enter the material
            // such that it corresponds to the time needed to reach the curren
            // location form the edge of the material. Only use this time step
            // for the computation of the interaction with the material, not for
            // the integration of the particles. This will ensure that the momenta
            // of all particles are reduced by approcimately the same amount in
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

void CollimatorPhysics::print(Inform &msg) {
    Inform::FmtFlags_t ff = msg.flags();

    if (totalPartsInMat_m > 0 ||
        bunchToMatStat_m > 0 ||
        rediffusedStat_m > 0 ||
        stoppedPartStat_m > 0) {

        OPALTimer::Timer time;
        msg << level2
            << "--- CollimatorPhysics - Name " << element_ref_m->getName()
            << " Material " << material_m << "\n"
            << "Particle Statistics @ " << time.time() << "\n"
            << std::setw(21) << "entered: " << Util::toStringWithThousandSep(bunchToMatStat_m) << "\n"
            << std::setw(21) << "rediffused: " << Util::toStringWithThousandSep(rediffusedStat_m) << "\n"
            << std::setw(21) << "stopped: " << Util::toStringWithThousandSep(stoppedPartStat_m) << "\n"
            << std::setw(21) << "total in material: " << Util::toStringWithThousandSep(totalPartsInMat_m)
            << endl;
    }

    msg.flags(ff);
}

bool CollimatorPhysics::stillActive() {
    return totalPartsInMat_m != 0;
}

bool CollimatorPhysics::stillAlive(PartBunchBase<double, 3> *bunch) {

    bool degraderAlive = true;

    //free GPU memory in case element is degrader, it is empty and bunch has moved past it
    if (collshape_m == ElementBase::DEGRADER && totalPartsInMat_m == 0) {
        Degrader *deg = static_cast<Degrader *>(element_ref_m);

        //get the size of the degrader
        double zBegin, zEnd;
        deg->getDimensions(zBegin, zEnd);

        //get the average Z position of the bunch
        Vector_t bunchOrigin = bunch->get_origin();

        //if bunch has moved past degrader free GPU memory
        if (bunchOrigin[2] > zBegin) {
            degraderAlive = false;
#ifdef OPAL_DKS
            if (IpplInfo::DKSEnabled)
                clearCollimatorDKS();
#endif
        }
    }

    return degraderAlive;

}

namespace {
    bool myCompF(PART x, PART y) {
      return x.label > y.label;
    }
}

void CollimatorPhysics::deleteParticleFromLocalVector() {
    /*
      the particle to be deleted (label < 0) are all at the end of
      the vector.
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

void CollimatorPhysics::push() {
    for (size_t i = 0; i < locParts_m.size(); ++i) {
        Vector_t &R  = locParts_m[i].Rincol;
        Vector_t &P  = locParts_m[i].Pincol;
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
void CollimatorPhysics::setTimeStepForLeavingParticles() {
    const double elementLength = element_ref_m->getElementLength();

    for (size_t i = 0; i < locParts_m.size(); ++i) {
        Vector_t &R  = locParts_m[i].Rincol;
        Vector_t &P  = locParts_m[i].Pincol;
        double   &dt = locParts_m[i].DTincol;
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

void CollimatorPhysics::resetTimeStep() {
    for (size_t i = 0; i < locParts_m.size(); ++i) {
        double &dt = locParts_m[i].DTincol;
        dt = dT_m;
    }

}


void CollimatorPhysics::gatherStatistics() {

    unsigned int locPartsInMat;
#ifdef OPAL_DKS
    if (collshape_m == ElementBase::DEGRADER && IpplInfo::DKSEnabled)
        locPartsInMat = numparticles_m + dksParts_m.size();
    else
        locPartsInMat = locParts_m.size();
#else
    locPartsInMat = locParts_m.size();
#endif

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


#ifdef OPAL_DKS

namespace {
    bool myCompFDKS(PART_DKS x, PART_DKS y) {
        return x.label > y.label;
    }
}

void CollimatorPhysics::applyDKS(PartBunchBase<double, 3> *bunch,
                                 const std::pair<Vector_t, double> &boundingSphere,
                                 size_t numParticlesInSimulation) {
    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    //if firs call to apply setup needed accelerator resources
    setupCollimatorDKS(bunch, numParticlesInSimulation);

    do {
        IpplTimings::startTimer(DegraderLoopTimer_m);

        //write particles to GPU if there are any to write
        if (dksParts_m.size() > 0) {
            //wrtie data from dksParts_m to the end of mem_ptr (offset = numparticles_m)
            dksbase_m.writeDataAsync<PART_DKS>(mem_ptr, &dksParts_m[0],
                                               dksParts_m.size(), -1, numparticles_m);

            //update number of particles on Device
            numparticles_m += dksParts_m.size();

            //free locParts_m vector
            dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());
        }

        //execute CollimatorPhysics kernels on GPU if any particles are there
        if (numparticles_m > 0) {
            dksbase_m.callCollimatorPhysics2(mem_ptr, par_ptr, numparticles_m);
        }

        //sort device particles and get number of particles comming back to bunch
        int numaddback = 0;
        if (numparticles_m > 0) {
            dksbase_m.callCollimatorPhysicsSort(mem_ptr, numparticles_m, numaddback);
        }

        //read particles from GPU if any are comming out of material
        if (numaddback > 0) {

            //resize dksParts_m to hold particles that need to go back to bunch
            dksParts_m.resize(numaddback);

            //read particles that need to be added back to bunch
            //particles that need to be added back are at the end of Device array
            dksbase_m.readData<PART_DKS>(mem_ptr, &dksParts_m[0], numaddback,
                                         numparticles_m - numaddback);

            //add particles back to the bunch
            for (unsigned int i = 0; i < dksParts_m.size(); ++i) {
                if (dksParts_m[i].label == -2) {
                    addBackToBunchDKS(bunch, i);
                    rediffusedStat_m++;
                } else {
                    stoppedPartStat_m++;
                    lossDs_m->addParticle(dksParts_m[i].Rincol, dksParts_m[i].Pincol,
                                          -locParts_m[dksParts_m[i].localID].IDincol);
                }
            }

            //erase particles that came from device from host array
            dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());

            //update number of particles on Device
            numparticles_m -= numaddback;
        }

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        if (onlyOneLoopOverParticles)
            copyFromBunchDKS(bunch, boundingSphere);

        T_m += dT_m;

        gatherStatistics();

        //more than one loop only if all the particles are in this degrader
        if (allParticleInMat_m) {
            onlyOneLoopOverParticles = (totalPartsInMat_m <= 0);
        } else {
            onlyOneLoopOverParticles = true;
        }

    } while (onlyOneLoopOverParticles == false);
}

void CollimatorPhysics::addBackToBunchDKS(PartBunchBase<double, 3> *bunch, unsigned i) {

    int id = dksParts_m[i].localID;

    bunch->createWithID(locParts_m[id].IDincol);

    /*
      Binincol is still <0, but now the particle is rediffused
      from the material and hence this is not a "lost" particle anymore
    */
    unsigned int last = bunch->getLocalNum() - 1;
    bunch->Bin[last] = 1;

    bunch->R[last]   = dksParts_m[i].Rincol;
    bunch->P[last]   = dksParts_m[i].Pincol;

    bunch->Q[last]   = locParts_m[id].Qincol;
    bunch->Bf[last]  = locParts_m[id].Bfincol;
    bunch->Ef[last]  = locParts_m[id].Efincol;
    bunch->dt[last]  = locParts_m[id].DTincol;

    dksParts_m[i].label = -1.0;

    ++ rediffusedStat_m;

}

void CollimatorPhysics::copyFromBunchDKS(PartBunchBase<double, 3> *bunch,
                                         const std::pair<Vector_t, double> &boundingSphere)
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
    const unsigned int minNumOfParticlesPerCore = bunch->getMinimumNumberOfParticlesPerCore();

    for (unsigned int i = 0; i < nL; ++i) {
        if ((bunch->Bin[i] == -1 || bunch->Bin[i] == 1) &&
            ((nL - ne) > minNumOfParticlesPerCore) &&
            hitTester_m->checkHit(bunch->R[i], bunch->P[i], dT_m))
        {

            PART x;
            x.localID      = numlocalparts_m; //unique id for each particle
            x.DTincol      = bunch->dt[i];
            x.IDincol      = bunch->ID[i];
            x.Binincol     = bunch->Bin[i];
            x.Rincol       = bunch->R[i];
            x.Pincol       = bunch->P[i];
            x.Qincol       = bunch->Q[i];
            x.Bfincol      = bunch->Bf[i];
            x.Efincol      = bunch->Ef[i];
            x.label        = 0;            // alive in matter

            PART_DKS x_gpu;
            x_gpu.label = x.label;
            x_gpu.localID = x.localID;
            x_gpu.Rincol = x.Rincol;
            x_gpu.Pincol = x.Pincol;

            locParts_m.push_back(x);
            dksParts_m.push_back(x_gpu);

            ne++;
            bunchToMatStat_m++;
            numlocalparts_m++;

            //mark particle to be deleted from bunch as soon as it enters the material
            bunch->destroy(1, i);
        }
    }
    if (ne > 0) {
      bunch->performDestroy(true);
    }
    IpplTimings::stopTimer(DegraderDestroyTimer_m);

}

void CollimatorPhysics::setupCollimatorDKS(PartBunchBase<double, 3> *bunch,
                                           size_t numParticlesInSimulation)
{

    if (curandInitSet_m == -1) {

        //int size = bunch->getLocalNum() + 0.5 * bunch->getLocalNum();
        //int size = bunch->getTotalNum() + 0.5 * bunch->getTotalNum();
        int size = numParticlesInSimulation;

        //allocate memory for parameters
        par_ptr = dksbase_m.allocateMemory<double>(numpar_m, ierr_m);

        //allocate memory for particles
        mem_ptr = dksbase_m.allocateMemory<PART_DKS>((int)size, ierr_m);

        maxparticles_m = (int)size;
        numparticles_m = 0;
        numlocalparts_m = 0;

        //reserve space for locParts_m vector
        locParts_m.reserve(size);

        //init curand
        dksbase_m.callInitRandoms(size, Options::seed);
        curandInitSet_m = 1;

        //create and transfer parameter array
        Degrader *deg = static_cast<Degrader *>(element_ref_m);
        double zBegin, zEnd;
        deg->getDimensions(zBegin, zEnd);

        double params[numpar_m] = {zBegin, zEnd, rho_m, Z_m,
                                   A_m, A2_c, A3_c, A4_c, A5_c, X0_m, I_m, dT_m, 1e-4
                                  };
        dksbase_m.writeDataAsync<double>(par_ptr, params, numpar_m);

    }

}

void CollimatorPhysics::clearCollimatorDKS() {
    if (curandInitSet_m == 1) {
        dksbase_m.freeMemory<double>(par_ptr, numpar_m);
        dksbase_m.freeMemory<PART_DKS>(mem_ptr, maxparticles_m);
        curandInitSet_m = -1;
    }

}

void CollimatorPhysics::deleteParticleFromLocalVectorDKS() {

    /*
      the particle to be deleted (label < 0) are all at the end of
      the vector.
    */
    sort(dksParts_m.begin(), dksParts_m.end(), myCompFDKS);

    // find start of particles to delete
    std::vector<PART_DKS>::iterator inv = dksParts_m.begin() + stoppedPartStat_m + rediffusedStat_m;

    /*
      for (; inv != dksParts_m.end(); inv++) {
      if ((*inv).label == -1)
      break;
      }
    */

    dksParts_m.erase(inv, dksParts_m.end());

}

#endif