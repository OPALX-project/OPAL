#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

// Distribution class
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.

#include <fstream>
#include <string>
#include <vector>

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Vektor.h"
#include "AppTypes/SymTenzor.h"
#include "Attributes/Attributes.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class Beam;
class Beamline;

template <class T, unsigned Dim>
class PartBunchBase;

class PartBins;
class EnvelopeBunch;
class LaserProfile;
class H5PartWrapper;

namespace DistrTypeT
{
    enum DistrTypeT {NODIST,
                     FROMFILE,
                     GAUSS,
                     BINOMIAL,
                     FLATTOP,
                     MULTIGAUSS,
                     GUNGAUSSFLATTOPTH,
                     ASTRAFLATTOPTH,
                     MATCHEDGAUSS
                    };
}

namespace EmissionModelT
{
    enum EmissionModelT {NONE,
                         ASTRA,
                         NONEQUIL
                        };
}

namespace InputMomentumUnitsT
{
    enum InputMomentumUnitsT {NONE,
                              EV
                              };
}

namespace Attrib
{
    namespace Distribution
    {
        enum AttributesT {
            TYPE,
            FNAME,
            WRITETOFILE,
            WEIGHT,
            INPUTMOUNITS,
            EMITTED,
            EMISSIONSTEPS,
            EMISSIONMODEL,
            EKIN,
            ELASER,
            W,
            FE,
            CATHTEMP,
            NBIN,
            XMULT,
            YMULT,
            ZMULT,
            TMULT,
            PXMULT,
            PYMULT,
            PZMULT,
            OFFSETX,
            OFFSETY,
            OFFSETZ,
            OFFSETT,
            OFFSETPX,
            OFFSETPY,
            OFFSETPZ,
            OFFSETP,
            SIGMAX,
            SIGMAY,
            SIGMAR,
            SIGMAZ,
            SIGMAT,
            TPULSEFWHM,
            TRISE,
            TFALL,
            SIGMAPX,
            SIGMAPY,
            SIGMAPZ,
            MX,
            MY,
            MZ,
            MT,
            CUTOFFX,
            CUTOFFY,
            CUTOFFR,
            CUTOFFLONG,
            CUTOFFPX,
            CUTOFFPY,
            CUTOFFPZ,
            FTOSCAMPLITUDE,
            FTOSCPERIODS,
            SEPPEAKS,
            NPEAKS,
            R,                          // the correlation matrix (a la transport)
            CORRX,
            CORRY,
            CORRZ,
            CORRT,
            R51,
            R52,
            R61,
            R62,
            LASERPROFFN,
            IMAGENAME,
            INTENSITYCUT,
            FLIPX,
            FLIPY,
            ROTATE90,
            ROTATE180,
            ROTATE270,
            EX,                         // below is for the matched distribution
            EY,
            ET,
            SECTOR,                     // Matched-Gauss: single sector or full machine
            NSECTORS,                   // Matched-Gauss: number of sectors to average field
            NSTEPS,                     // Matched-Gauss: number of steps for closed orbit finder
            RGUESS,                     // Matched-Gauss: guess for closed orbit finder
            DENERGY,                    // Matched-Gauss: energy step for closed orbit finder
            LINE,
            RESIDUUM,
            MAXSTEPSCO,
            MAXSTEPSSI,
            ORDERMAPS,
            //            E2,
            ID1,                       // special particle that the user can set
            ID2,                       // special particle that the user can set
            SCALABLE,
            SIZE
        };
    }

    namespace Legacy
    {
        namespace Distribution
        {
            enum LegacyAttributesT {
                // DESCRIPTION OF THE DISTRIBUTION:
                DISTRIBUTION = Attrib::Distribution::SIZE,
                // DEBIN,
                SBIN,
                SIGMAPT,
                CUTOFF,
                T,
                PT,
                // ALPHAX,
                // ALPHAY,
                // BETAX,
                // BETAY,
                // DX,
                // DDX,
                // DY,
                // DDY,
                SIZE
            };
        }
    }
}

/*
 * Class Distribution
 *
 * Defines the initial beam that is injected or emitted into the simulation.
 */

class Distribution: public Definition {

public:

    Distribution();
    virtual ~Distribution();

    virtual bool canReplaceBy(Object *object);
    virtual Distribution *clone(const std::string &name);
    virtual void execute();
    virtual void update();
    size_t getNumOfLocalParticlesToCreate(size_t n);
    void createOpalCycl(PartBunchBase<double, 3> *beam,
                        size_t numberOfParticles,
                        double current, const Beamline &bl);
    void createOpalE(Beam *beam,
                     std::vector<Distribution *> addedDistributions,
                     EnvelopeBunch *envelopeBunch,
                     double distCenter,
                     double Bz0);
    void createOpalT(PartBunchBase<double, 3> *beam,
                     std::vector<Distribution *> addedDistributions,
                     size_t &numberOfParticles);
    void createOpalT(PartBunchBase<double, 3> *beam, size_t &numberOfParticles);
    void doRestartOpalT(PartBunchBase<double, 3> *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    void doRestartOpalCycl(PartBunchBase<double, 3> *p, size_t Np, int restartStep,
                        const int specifiedNumBunch, H5PartWrapper *h5wrapper);
    void doRestartOpalE(EnvelopeBunch *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    size_t emitParticles(PartBunchBase<double, 3> *beam, double eZ);
    double getPercentageEmitted() const;
    static Distribution *find(const std::string &name);

    bool getIfDistEmitting();
    int getLastEmittedEnergyBin();
    double getMaxTOrZ();
    double getMinTOrZ();
    size_t getNumberOfEmissionSteps();
    int getNumberOfEnergyBins();
    double getEmissionDeltaT();
    double getEnergyBinDeltaT();
    double getWeight();

    double getTEmission();

    Vector_t get_pmean() const;

    std::string getTypeofDistribution();

    Inform &printInfo(Inform &os) const;

    bool Rebin();
    void setDistToEmitted(bool emitted);
    void setDistType();
    void setSigmaR_m();
    void setSigmaP_m(double massIneV);
    void shiftBeam(double &maxTOrZ, double &minTOrZ);
    double getEmissionTimeShift() const;

    void setNumberOfDistributions(unsigned int n) { numberOfDistributions_m = n; }

    DistrTypeT::DistrTypeT getType() const;
private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussTest, FullSigmaTest1);
    FRIEND_TEST(GaussTest, FullSigmaTest2);
    FRIEND_TEST(BinomialTest, FullSigmaTest1);
    FRIEND_TEST(BinomialTest, FullSigmaTest2);
#endif

    Distribution(const std::string &name, Distribution *parent);

    // Not implemented.
    Distribution(const Distribution &) = delete;
    void operator=(const Distribution &) = delete;

    std::vector<double>& getXDist();
    std::vector<double>& getBGxDist();
    std::vector<double>& getYDist();
    std::vector<double>& getBGyDist();
    std::vector<double>& getTOrZDist();
    std::vector<double>& getBGzDist();
    void eraseXDist();
    void eraseBGxDist();
    void eraseYDist();
    void eraseBGyDist();
    void eraseTOrZDist();
    void eraseBGzDist();

    //    void printSigma(SigmaGenerator<double,unsigned int>::matrix_type& M, Inform& out);
    void addDistributions();
    void applyEmissionModel(double lowEnergyLimit, double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void applyEmissModelAstra(double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void applyEmissModelNone(double &pz);
    void applyEmissModelNonEquil(double eZ, double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void create(size_t &numberOfParticles, double massIneV);
    void calcPartPerDist(size_t numberOfParticles);
    void checkEmissionParameters();
    void checkIfEmitted();
    void checkParticleNumber(size_t &numberOfParticles);
    void chooseInputMomentumUnits(InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits);
    double converteVToBetaGamma(double valueIneV, double massIneV);
    size_t getNumberOfParticlesInFile(std::ifstream &inputFile);

    class BinomialBehaviorSplitter {
    public:
        virtual ~BinomialBehaviorSplitter()
        { }

        virtual double get(double rand) = 0;
    };

    class MDependentBehavior: public BinomialBehaviorSplitter {
    public:
        MDependentBehavior(const MDependentBehavior &rhs):
            ami_m(rhs.ami_m)
        {}

        MDependentBehavior(double a)
        { ami_m = 1.0 / a; }

        virtual double get(double rand);
    private:
        double ami_m;
    };

    class GaussianLikeBehavior: public BinomialBehaviorSplitter {
    public:
        virtual double get(double rand);
    };

    void createDistributionBinomial(size_t numberOfParticles, double massIneV);
    void createDistributionFlattop(size_t numberOfParticles, double massIneV);
    void createDistributionMultiGauss(size_t numberOfParticles, double massIneV);
    void createDistributionFromFile(size_t numberOfParticles, double massIneV);
    void createDistributionGauss(size_t numberOfParticles, double massIneV);
    void createMatchedGaussDistribution(size_t numberOfParticles, double massIneV);
    void sampleUniformDisk(gsl_qrng* quasiRandGen2D, double& x1, double& x2);
    void fillEBinHistogram();
    void fillParticleBins();
    size_t findEBin(double tOrZ);
    void generateAstraFlattopT(size_t numberOfParticles);
    void generateBinomial(size_t numberOfParticles);
    void generateFlattopLaserProfile(size_t numberOfParticles);
    void generateFlattopT(size_t numberOfParticles);
    void generateFlattopZ(size_t numberOfParticles);
    void generateGaussZ(size_t numberOfParticles);
    void generateLongFlattopT(size_t numberOfParticles);
    void generateTransverseGauss(size_t numberOfParticles);
    void initializeBeam(PartBunchBase<double, 3> *beam);
    void injectBeam(PartBunchBase<double, 3> *beam);
    void printDist(Inform &os, size_t numberOfParticles) const;
    void printDistBinomial(Inform &os) const;
    void printDistFlattop(Inform &os) const;
    void printDistMultiGauss(Inform &os) const;
    void printDistFromFile(Inform &os) const;
    void printDistGauss(Inform &os) const;
    void printDistMatchedGauss(Inform &os) const;
    void printEmissionModel(Inform &os) const;
    void printEmissionModelAstra(Inform &os) const;
    void printEmissionModelNone(Inform &os) const;
    void printEmissionModelNonEquil(Inform &os) const;
    void printEnergyBins(Inform &os) const;
    void adjustPhaseSpace(double massIneV);
    void reflectDistribution(size_t &numberOfParticles);
    void scaleDistCoordinates();
    /// Select and allocate gsl random number generator
    gsl_qrng* selectRandomGenerator(std::string, unsigned int dimension);
    void setAttributes();
    void setDistParametersBinomial(double massIneV);
    void setDistParametersFlattop(double massIneV);
    void setDistParametersMultiGauss(double massIneV);
    void setDistParametersGauss(double massIneV);
    void setEmissionTime(double &maxT, double &minT);
    void setupEmissionModel(PartBunchBase<double, 3> *beam);
    void setupEmissionModelAstra(PartBunchBase<double, 3> *beam);
    void setupEmissionModelNone(PartBunchBase<double, 3> *beam);
    void setupEmissionModelNonEquil();
    void setupEnergyBins(double maxTOrZ, double minTOrZ);
    void setupParticleBins(double massIneV, PartBunchBase<double, 3> *beam);
    void shiftDistCoordinates(double massIneV);
    void writeOutFileHeader();
    void writeOutFileEmission();
    void writeOutFileInjection();

    std::string distT_m;                 /// Distribution type. Declared as string
    DistrTypeT::DistrTypeT distrTypeT_m; /// and list type for switch statements.

    unsigned int numberOfDistributions_m;

    bool emitting_m;                     /// Distribution is an emitted, and is currently
                                         /// emitting, rather than an injected, beam.

    PartData particleRefData_m;          /// Reference data for particle type (charge,
                                         /// mass etc.)

    /// Vector of distributions to be added to this distribution.
    std::vector<Distribution *> addedDistributions_m;
    std::vector<size_t> particlesPerDist_m;

    /// Emission Model.
    EmissionModelT::EmissionModelT emissionModel_m;

    /// Emission parameters.
    double tEmission_m;
    static const double percentTEmission_m; /// Increase tEmission_m by twice this percentage
                                            /// to ensure that no particles fall on the leading edge of
                                            /// the first emission time step or the trailing edge of the last
                                            /// emission time step.
    double tBin_m;
    double currentEmissionTime_m;
    int currentEnergyBin_m;
    int currentSampleBin_m;
    int numberOfEnergyBins_m;       /// Number of energy bins the distribution
                                    /// is broken into. Used for an emitted beam.
    int numberOfSampleBins_m;       /// Number of samples to use per energy bin
                                    /// when emitting beam.
    PartBins *energyBins_m;         /// Distribution energy bins.
    gsl_histogram *energyBinHist_m; /// GSL histogram used to define energy bin
                                    /// structure.

    gsl_rng *randGen_m;             /// Random number generator

    // ASTRA and NONE photo emission model.
    double pTotThermal_m;           /// Total thermal momentum.
    Vector_t pmean_m;

    // NONEQUIL photo emission model.
    double cathodeWorkFunc_m;       /// Cathode material work function (eV).
    double laserEnergy_m;           /// Laser photon energy (eV).
    double cathodeFermiEnergy_m;    /// Cathode material Fermi energy (eV).
    double cathodeTemp_m;           /// Cathode temperature (K).
    double emitEnergyUpperLimit_m;  /// Upper limit on emission energy distribution (eV).

    std::vector<std::vector<double> > additionalRNs_m;

    size_t totalNumberParticles_m;
    size_t totalNumberEmittedParticles_m;

    // Beam coordinate containers.
    std::vector<double> xDist_m;
    std::vector<double> pxDist_m;
    std::vector<double> yDist_m;
    std::vector<double> pyDist_m;
    std::vector<double> tOrZDist_m;
    std::vector<double> pzDist_m;

    // Initial coordinates for file write.
    std::vector<double> xWrite_m;
    std::vector<double> pxWrite_m;
    std::vector<double> yWrite_m;
    std::vector<double> pyWrite_m;
    std::vector<double> tOrZWrite_m;
    std::vector<double> pzWrite_m;
    std::vector<size_t> binWrite_m;

    // for compatibility reasons
    double avrgpz_m;

    //Distribution parameters.
    InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits_m;
    double sigmaTRise_m;
    double sigmaTFall_m;
    double tPulseLengthFWHM_m;
    Vector_t sigmaR_m;
    Vector_t sigmaP_m;
    Vector_t cutoffR_m;
    Vector_t cutoffP_m;
    Vector_t mBinomial_m;
    SymTenzor<double, 6> correlationMatrix_m;

    // Parameters multiGauss distribution.
    double sepPeaks_m;
    unsigned nPeaks_m;

    // Laser profile.
    std::string laserProfileFileName_m;
    std::string laserImageName_m;
    double laserIntensityCut_m;
    LaserProfile *laserProfile_m;

    // AAA This is for the matched distribution
    double I_m;
    double E_m;

    /// time binned distribution with thermal energy
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;
};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.printInfo(os);
}

inline
Vector_t Distribution::get_pmean() const {
    return pmean_m;
}

inline
DistrTypeT::DistrTypeT Distribution::getType() const {
    return distrTypeT_m;
}

inline
double Distribution::getPercentageEmitted() const {
    return (double)totalNumberEmittedParticles_m / (double)totalNumberParticles_m;
}

inline
std::string Distribution::getTypeofDistribution() {
    return (std::string) Attributes::getString(itsAttr[Attrib::Distribution::TYPE]);
}

#endif // OPAL_Distribution_HH