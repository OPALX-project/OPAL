#ifndef COLLIMATORPHYSICS_HH
#define COLLIMATORPHYSICS_HH
//Class:CollimatorPhysics
//  Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang, Stachel, Adelmann$
//-------------------------------------------------------------------------
#include <vector>
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Degrader.h"
#include <gsl/gsl_rng.h>

#include "Utility/IpplTimings.h"

#ifdef OPAL_DKS
#include "DKSOPAL.h"
#endif

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;
class Inform;

#ifdef OPAL_DKS
typedef struct __align__(16) {
    int label;
    unsigned localID;
    Vector_t Rincol;
    Vector_t Pincol;
    long IDincol;
    int Binincol;
    double DTincol;
    double Qincol;
    Vector_t Bfincol;
    Vector_t Efincol;
} PART;

typedef struct {
    int label;
    unsigned localID;
    Vector_t Rincol;
    Vector_t Pincol;
} PART_DKS;

#else
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
#endif

struct InsideTester {
    virtual ~InsideTester()
    { }

    virtual
    bool checkHit(const Vector_t &R, const Vector_t &P, double dt) = 0;
};


class CollimatorPhysics: public ParticleMatterInteractionHandler {
public:
    CollimatorPhysics(const std::string &name, ElementBase *element, std::string &mat);
    ~CollimatorPhysics();

    void apply(PartBunchBase<double, 3> *bunch,
               const std::pair<Vector_t, double> &boundingSphere,
               size_t numParticlesInSimulation = 0);

    virtual const std::string getType() const;

    void print(Inform& os);
    bool stillActive();
    bool stillAlive(PartBunchBase<double, 3> *bunch);

    double getTime();
    std::string getName();
    size_t getParticlesInMat();
    unsigned getRediffused();
    unsigned int getNumEntered();
    void doPhysics(PartBunchBase<double, 3> *bunch);

    virtual bool computeEnergyLoss(double &Eng, const double deltat, bool includeFluctuations = true) const;

private:

    void configureMaterialParameters();
    void computeCoulombScattering(Vector_t &R, Vector_t &P, const double &deltat);

    void applyRotation(double &px,
                       double &pz,
                       double &x,
                       double &z,
                       double xplane,
                       double Norm_P,
                       double thetacou,
                       double deltas,
                       int coord);

    void copyFromBunch(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere);
    void addBackToBunch(PartBunchBase<double, 3> *bunch, unsigned i);

    void applyNonDKS(PartBunchBase<double, 3> *bunch,
                     const std::pair<Vector_t, double> &boundingSphere,
                     size_t numParticlesInSimulation);
#ifdef OPAL_DKS
    void applyDKS(PartBunchBase<double, 3> *bunch,
                  const std::pair<Vector_t, double> &boundingSphere,
                  size_t numParticlesInSimulation);
    void copyFromBunchDKS(PartBunchBase<double, 3> *bunch,
                        const std::pair<Vector_t, double> &boundingSphere);
    void addBackToBunchDKS(PartBunchBase<double, 3> *bunch, unsigned i);

    void setupCollimatorDKS(PartBunchBase<double, 3> *bunch, size_t numParticlesInSimulation);
    void clearCollimatorDKS();

    void applyDKS();
    void deleteParticleFromLocalVectorDKS();

#endif
    void deleteParticleFromLocalVector();

    void calcStat(double Eng);

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
    // number of local particles that are in the material
    size_t locPartsInMat_m;

    // some statistics
    double Eavg_m;                            // average kinetic energy
    double Emax_m;                            // maximum kinetic energy
    double Emin_m;                            // minimum kinetic energy

    std::vector<PART> locParts_m;             // local particles that are in material

    std::unique_ptr<LossDataSink> lossDs_m;

#ifdef OPAL_DKS
    DKSOPAL dksbase;
    int curandInitSet;

    int ierr;
    int maxparticles;
    int numparticles;
    int numlocalparts;
    void *par_ptr;
    void *mem_ptr;

    std::vector<PART_DKS> dksParts_m;

    static const int numpar;
#endif

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
    return locPartsInMat_m;
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

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
