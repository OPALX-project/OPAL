#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

// ------------------------------------------------------------------------
// $RCSfile: Distribution.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------
#include <iosfwd>
#include <fstream>
#include <string>

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include "Algorithms/Vektor.h"
#include "Beamlines/Beamline.h"
#include "Attributes/Attributes.h"

#include "Ippl.h"

#include "H5hut.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class Beam;

template <class T, unsigned Dim>
class PartBunchBase;

class PartBins;
class EnvelopeBunch;
class BoundaryGeometry;
class LaserProfile;
class H5PartWrapper;

namespace DistrTypeT
{
    enum DistrTypeT {NODIST,
                     FROMFILE,
                     GAUSS,
                     BINOMIAL,
                     FLATTOP,
                     SURFACEEMISSION,
                     SURFACERANDCREATE,
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
            NPDARKCUR,
            INWARDMARGIN,
            EINITHR,
            FNA,
            FNB,
            FNY,
            FNVYZERO,
            FNVYSECOND,
            FNPHIW,
            FNBETA,
            FNFIELDTHR,
            FNMAXEMI,
            SECONDARYFLAG,
            NEMISSIONMODE,
            VSEYZERO,                   // sey_0 in Vaughn's model.
            VEZERO,                     // Energy related to sey_0 in Vaughan's model.
            VSEYMAX,                    // sey max in Vaughan's model.
            VEMAX,                      // Emax in Vaughan's model.
            VKENERGY,                   // Fitting parameter denotes the roughness of
            // surface for impact energy in Vaughn's model.
            VKTHETA,                    // Fitting parameter denotes the roughness of
            // surface for impact angle in Vaughn's model.
            VVTHERMAL,                  // Thermal velocity of Maxwellian distribution
            // of secondaries in Vaughan's model.
            VW,
            SURFMATERIAL,               // Add material type, currently 0 for copper
            // and 1 for stainless steel.
            EX,                         // below is for the matched distribution
            EY,
            ET,
            SECTOR,                     // Matched-Gauss: single sector or full machine
            NSTEPS,                     // Matched-Gauss: number of steps for closed orbit finder
            LINE,
            RESIDUUM,
            MAXSTEPSCO,
            MAXSTEPSSI,
            ORDERMAPS,
            //            E2,
            RGUESS,
            DENERGY,
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
    void createBoundaryGeometry(PartBunchBase<double, 3> *p, BoundaryGeometry &bg);
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
    void createPriPart(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg);
    void doRestartOpalT(PartBunchBase<double, 3> *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    void doRestartOpalCycl(PartBunchBase<double, 3> *p, size_t Np, int restartStep,
                        const int specifiedNumBunch, H5PartWrapper *h5wrapper);
    void doRestartOpalE(EnvelopeBunch *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    size_t emitParticles(PartBunchBase<double, 3> *beam, double eZ);
    double getPercentageEmitted() const;
    static Distribution *find(const std::string &name);

    void eraseXDist();
    void eraseBGxDist();
    void eraseYDist();
    void eraseBGyDist();
    void eraseTOrZDist();
    void eraseBGzDist();
    bool getIfDistEmitting();
    int getLastEmittedEnergyBin();
    double getMaxTOrZ();
    double getMinTOrZ();
    size_t getNumberOfEmissionSteps();
    int getNumberOfEnergyBins();
    double getEmissionDeltaT();
    double getEnergyBinDeltaT();
    double getWeight();
    std::vector<double>& getXDist();
    std::vector<double>& getBGxDist();
    std::vector<double>& getYDist();
    std::vector<double>& getBGyDist();
    std::vector<double>& getTOrZDist();
    std::vector<double>& getBGzDist();

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;
    double getTEmission();

    Vector_t get_pmean() const;
    double getEkin() const;
    double getLaserEnergy() const;
    double getWorkFunctionRf() const;

    size_t getNumberOfDarkCurrentParticles();
    double getDarkCurrentParticlesInwardMargin();
    double getEInitThreshold();
    double getWorkFunction();
    double getFieldEnhancement();
    size_t getMaxFNemissionPartPerTri();
    double getFieldFNThreshold();
    double getFNParameterA();
    double getFNParameterB();
    double getFNParameterY();
    double getFNParameterVYZero();
    double getFNParameterVYSecond();
    int    getSecondaryEmissionFlag();
    bool   getEmissionMode() ;

    std::string getTypeofDistribution();

    double getvSeyZero();//return sey_0 in Vaughan's model
    double getvEZero();//return the energy related to sey_0 in Vaughan's model
    double getvSeyMax();//return sey max in Vaughan's model
    double getvEmax();//return Emax in Vaughan's model
    double getvKenergy();//return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
    double getvKtheta();//return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
    double getvVThermal();//return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
    double getVw();//return velocity scalar for parallel plate benchmark;
    int getSurfMaterial();//material type for Furman-Pivi's model 0 for copper, 1 for stainless steel

    Inform &printInfo(Inform &os) const;

    bool Rebin();
    void setDistToEmitted(bool emitted);
    void setDistType();
    void shiftBeam(double &maxTOrZ, double &minTOrZ);
    double getEmissionTimeShift() const;

    bool GetPreviousH5Local() {return previousH5Local_m;}

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
    void createDistributionFromFile(size_t numberOfParticles, double massIneV);
    void createDistributionGauss(size_t numberOfParticles, double massIneV);
    void createMatchedGaussDistribution(size_t numberOfParticles, double massIneV);
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
    void printDistFromFile(Inform &os) const;
    void printDistGauss(Inform &os) const;
    void printDistMatchedGauss(Inform &os) const;
    void printDistSurfEmission(Inform &os) const;
    void printDistSurfAndCreate(Inform &os) const;
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
    void setDistParametersGauss(double massIneV);
    void setEmissionTime(double &maxT, double &minT);
    void setFieldEmissionParameters();
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

    // Laser profile.
    std::string laserProfileFileName_m;
    std::string laserImageName_m;
    double laserIntensityCut_m;
    LaserProfile *laserProfile_m;

    /*
     * Dark current calculation parameters.
     */
    size_t darkCurrentParts_m;      /// Number of dark current particles.
    double darkInwardMargin_m;      /// Dark current particle initialization position.
                                    /// Inward along the triangle normal, positive.
                                    /// Inside the geometry.
    double eInitThreshold_m;        /// Field threshold (MV/m) beyond which particles
                                    /// will be initialized.
    double workFunction_m;          /// Work function of surface material (eV).
    double fieldEnhancement_m;      /// Field enhancement factor beta for Fowler-
                                    /// Nordheim emission.
    double fieldThrFN_m;            /// Field threshold for Fowler-Nordheim
                                    /// emission (MV/m).
    size_t maxFN_m;                 /// Max. number of electrons emitted from a
                                    /// single triangle for Fowler-Nordheim emission.
    double paraFNA_m;               /// Empirical constant A in Fowler-Nordheim
                                    /// emission model.
    double paraFNB_m;               /// Empirical constant B in Fowler-Nordheim
                                    /// emission model.
    double paraFNY_m;               /// Constant for image charge effect parameter y(E)
                                    /// in Fowler-Nordheim emission model.
    double paraFNVYSe_m;            /// Second order constant for v(y) function in
                                    /// Fowler-Nordheim emission model.
    double paraFNVYZe_m;            /// Zero order constant for v(y) function in
                                    /// Fowler-Nordheim emission model.
    int    secondaryFlag_m;         /// Select the secondary model type:
                                    ///     0           ==> no secondary emission.
                                    ///     1           ==> Furman-Pivi
                                    ///     2 or larger ==> Vaughan's model
    double ppVw_m;                  /// Velocity scalar for parallel plate benchmark.
    double vVThermal_m;             /// Thermal velocity of Maxwellian distribution
                                    /// of secondaries in Vaughan's model.


    // AAA This is for the matched distribution
    double I_m;
    double E_m;

    /// time binned distribution with thermal energy
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;

    // Cyclotron for restart in local mode
    bool previousH5Local_m;
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
double Distribution::getEkin() const {
    return Attributes::getReal(itsAttr[Attrib::Distribution::EKIN]);
}

inline
double Distribution::getLaserEnergy() const {
    return Attributes::getReal(itsAttr[Attrib::Distribution::ELASER]);
}

inline
double Distribution::getWorkFunctionRf() const {
    return Attributes::getReal(itsAttr[Attrib::Distribution::W]);
}

inline
size_t Distribution::getNumberOfDarkCurrentParticles() {
    return (size_t) Attributes::getReal(itsAttr[Attrib::Distribution::NPDARKCUR]);
}

inline
double Distribution::getDarkCurrentParticlesInwardMargin() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::INWARDMARGIN]);
}

inline
double Distribution::getEInitThreshold() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::EINITHR]);
}

inline
double Distribution::getWorkFunction() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNPHIW]);
}

inline
double Distribution::getFieldEnhancement() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNBETA]);
}

inline
size_t Distribution::getMaxFNemissionPartPerTri() {
    return (size_t) Attributes::getReal(itsAttr[Attrib::Distribution::FNMAXEMI]);
}

inline
double Distribution::getFieldFNThreshold() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNFIELDTHR]);
}

inline
double Distribution::getFNParameterA() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNA]);
}

inline
double Distribution::getFNParameterB() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNB]);
}

inline
double Distribution::getFNParameterY() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNY]);
}

inline
double Distribution::getFNParameterVYZero() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNVYZERO]);
}

inline
double Distribution::getFNParameterVYSecond() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::FNVYSECOND]);
}

inline
int Distribution::getSecondaryEmissionFlag() {
    return Attributes::getReal(itsAttr[Attrib::Distribution::SECONDARYFLAG]);
}

inline
bool Distribution::getEmissionMode() {
    return Attributes::getBool(itsAttr[Attrib::Distribution::NEMISSIONMODE]);
}

inline
std::string Distribution::getTypeofDistribution() {
    return (std::string) Attributes::getString(itsAttr[Attrib::Distribution::TYPE]);
}

inline
double Distribution::getvSeyZero() {
    // return sey_0 in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VSEYZERO]);
}

inline
double Distribution::getvEZero() {
    // return the energy related to sey_0 in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VEZERO]);
}

inline
double Distribution::getvSeyMax() {
    // return sey max in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VSEYMAX]);
}

inline
double Distribution::getvEmax() {
    // return Emax in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VEMAX]);
}

inline
double Distribution::getvKenergy() {
    // return fitting parameter denotes the roughness of surface for
    // impact energy in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VKENERGY]);
}

inline
double Distribution::getvKtheta() {
    // return fitting parameter denotes the roughness of surface for
    // impact angle in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VKTHETA]);
}

inline
double Distribution::getvVThermal() {
    // thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
    return Attributes::getReal(itsAttr[Attrib::Distribution::VVTHERMAL]);
}

inline
double Distribution::getVw() {
    // velocity scalar for parallel plate benchmark;
    return Attributes::getReal(itsAttr[Attrib::Distribution::VW]);
}

inline
int Distribution::getSurfMaterial() {
    // Surface material number for Furman-Pivi's Model;
    return (int)Attributes::getReal(itsAttr[Attrib::Distribution::SURFMATERIAL]);
}

#endif // OPAL_Distribution_HH
