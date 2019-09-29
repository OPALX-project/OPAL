#ifndef PART_BUNCH_BASE_H
#define PART_BUNCH_BASE_H

#include "Ippl.h"
#include "Particle/AbstractParticle.h" //TODO should be in Ippl.h
#include "Algorithms/PBunchDefs.h"
#include "Algorithms/OpalParticle.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBinsCyc.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Quaternion.h"

#include <iosfwd>
#include <vector>

#include "Structure/LossDataSink.h"
#include "Structure/FieldSolver.h"
#include "Algorithms/ListElem.h"

class Distribution;

template <class T, int, int> class FMatrix;
template <class T, int> class FVector;

template <class T, unsigned Dim>
class PartBunchBase
{
public:
    typedef typename AbstractParticle<T, Dim>::ParticlePos_t   ParticlePos_t;
    typedef typename AbstractParticle<T, Dim>::ParticleIndex_t ParticleIndex_t;
    typedef typename AbstractParticle<T, Dim>::UpdateFlags     UpdateFlags_t;
    typedef typename AbstractParticle<T, Dim>::Position_t      Position_t;

    typedef std::pair<Vector_t, Vector_t> VectorPair_t;

    static const unsigned Dimension = Dim;

    enum UnitState_t { units = 0, unitless = 1 };

public:

    explicit PartBunchBase(AbstractParticle<T, Dim>* pb);

    virtual ~PartBunchBase() { }

    PartBunchBase(AbstractParticle<T, Dim>* pb, const PartData *ref);

    PartBunchBase(const PartBunchBase &rhs) = delete; // implement if needed

    /*
     * Bunch common member functions
     */

    // This is required since we initialize the Layout and the RegionLayout with default constructor
    virtual void initialize(FieldLayout_t *fLayout) = 0;

    bool getIfBeamEmitting();

    int getLastEmittedEnergyBin();

    size_t getNumberOfEmissionSteps();

    int getNumberOfEnergyBins();

    void Rebin();

    void setEnergyBins(int numberOfEnergyBins);

    bool weHaveEnergyBins();

    //FIXME: unify methods, use convention that all particles have own dt
    void switchToUnitlessPositions(bool use_dt_per_particle = false);

    //FIXME: unify methods, use convention that all particles have own dt
    void switchOffUnitlessPositions(bool use_dt_per_particle = false);

    void setDistribution(Distribution *d,
                         std::vector<Distribution *> addedDistributions,
                         size_t &np);

    bool isGridFixed();

    bool hasBinning();


    /*
       Energy bins related functions
     */

    void setTEmission(double t);

    double getTEmission();

    bool doEmission();

    bool weHaveBins() const;

    void setPBins(PartBins *pbin);

    void setPBins(PartBinsCyc *pbin);

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */
    size_t emitParticles(double eZ);

    void updateNumTotal();

    void rebin();

    int getNumBins();

    int getLastemittedBin();

    void setLocalBinCount(size_t num, int bin);

    /** \brief Compute the gammas of all bins */
    void calcGammas();

    void calcGammas_cycl();

    /** \brief Get gamma of one bin */
    double getBinGamma(int bin);

    /** \brief Set the charge of one bin to the value of q and all other to zero */
    virtual void setBinCharge(int bin, double q);

    /** \brief Set the charge of all other the ones in bin to zero */
    virtual void setBinCharge(int bin);

    /** \brief returns the number of particles outside of a box defined by x */
    size_t calcNumPartsOutside(Vector_t x);

    void calcLineDensity(unsigned int nBins, std::vector<double> &lineDensity,
                         std::pair<double, double> &meshInfo);

    void setBeamFrequency(double v);

    /*
       Mesh and Field Layout related functions
     */

    virtual void boundp();

    /** delete particles which are too far away from the center of beam*/
    void boundp_destroy();

    /** This is only temporary in order to get the collimator and pepperpot working */
    size_t boundp_destroyT();

    size_t destroyT();

    /*
       Read out coordinates
     */

    virtual double getPx(int i);
    virtual double getPy(int i);
    virtual double getPz(int i);

    virtual double getPx0(int i);
    virtual double getPy0(int i);

    virtual double getX(int i);
    virtual double getY(int i);
    virtual double getZ(int i);

    virtual double getX0(int i);
    virtual double getY0(int i);

    virtual void setZ(int i, double zcoo);

    void get_bounds(Vector_t &rmin, Vector_t &rmax);

    void getLocalBounds(Vector_t &rmin, Vector_t &rmax);

    std::pair<Vector_t, double> getBoundingSphere();

    std::pair<Vector_t, double> getLocalBoundingSphere();


    /*
       Compatibility function push_back
     */

    void push_back(OpalParticle p);

    void set_part(FVector<double, 6> z, int ii);

    void set_part(OpalParticle p, int ii);

    OpalParticle get_part(int ii);

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].
    void maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                           double &axmax, double &aymax);

    void   setdT(double dt);
    double getdT() const;

    void   setT(double t);
    void   incrementT();
    double getT() const;

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    double get_sPos();

    void set_sPos(double s);

    double get_gamma() const;

    double get_meanKineticEnergy() const;
    Vector_t get_origin() const;
    Vector_t get_maxExtent() const;
    Vector_t get_centroid() const;
    Vector_t get_rrms() const;
    Vector_t get_rprms() const;
    Vector_t get_rmean() const;
    Vector_t get_prms() const;
    Vector_t get_pmean() const;
    Vector_t get_pmean_Distribution() const;
    Vector_t get_emit() const;
    Vector_t get_norm_emit() const;
    Vector_t get_halo() const;
    virtual Vector_t get_hr() const;

    double get_Dx() const;
    double get_Dy() const;

    double get_DDx() const;
    double get_DDy() const;

    virtual void set_meshEnlargement(double dh);

    void gatherLoadBalanceStatistics();
    size_t getLoadBalance(int p) const;

    void get_PBounds(Vector_t &min, Vector_t &max) const;

    void calcBeamParameters();

    void calcBeamParametersInitial(); // Calculate initial beam parameters before emission.

    double getCouplingConstant() const;
    void setCouplingConstant(double c);

    // set the charge per simulation particle
    void setCharge(double q);
    // set the charge per simulation particle when total particle number equals 0
    void setChargeZeroPart(double q);

    // set the mass per simulation particle
    void setMass(double mass);

    /// \brief Need Ek for the Schottky effect calculation (eV)
    double getEkin() const;

    /// Need the work function for the Schottky effect calculation (eV)
    double getWorkFunctionRf() const;

    /// Need the laser energy for the Schottky effect calculation (eV)
    double getLaserEnergy() const;

    /// get the total charge per simulation particle
    double getCharge() const;

    /// get the macro particle charge
    double getChargePerParticle() const;

    virtual void setSolver(FieldSolver *fs);

    bool hasFieldSolver();

    std::string getFieldSolverType() const;

    void setStepsPerTurn(int n);
    int getStepsPerTurn() const;

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n);
    long long getGlobalTrackStep() const;

    /// step in a TRACK command
    void setLocalTrackStep(long long n);
    void incTrackSteps();
    long long getLocalTrackStep() const;

    void setNumBunch(short n);
    short getNumBunch() const;

    // used in ParallelCyclotronTracker for multi-bunch mode
    void setTotalNumPerBunch(size_t numpart, short n);
    size_t getTotalNumPerBunch(short n) const;

    void setLocalNumPerBunch(size_t numpart, short n);
    size_t getLocalNumPerBunch(short n) const;

    /* used in initializeTracking_m of ParallelCyclotronTracker
     * for multi-bunch mode
     */
    void countTotalNumPerBunch();

    void setGlobalMeanR(Vector_t globalMeanR);
    Vector_t getGlobalMeanR();
    void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
    Quaternion_t getGlobalToLocalQuaternion();

    void setSteptoLastInj(int n);
    int getSteptoLastInj();

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi();

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
    bool resetPartBinID2(const double eta);

    bool resetPartBinBunch();

    ///@{ Access to reference data
    double getQ() const;
    double getM() const;
    double getP() const;
    double getE() const;
    ParticleType::type getPType() const;
    double getInitialBeta() const;
    double getInitialGamma() const;
    ///@}
    ///@{ Set reference data
    void resetQ(double q);
    void resetM(double m);
    void setPType(ParticleType::type);
    ///@}
    double getdE();
    virtual double getGamma(int i);
    virtual double getBeta(int i);
    virtual void actT();

    const PartData *getReference() const;

    double getEmissionDeltaT();

    Quaternion_t getQKs3D();
    void         setQKs3D(Quaternion_t q);
    Vector_t     getKs3DRefr();
    void         setKs3DRefr(Vector_t r);
    Vector_t     getKs3DRefp();
    void         setKs3DRefp(Vector_t p);

    void iterateEmittedBin(int binNumber);

    void calcEMean();

    void correctEnergy(double avrgp);

    Inform &print(Inform &os);

    /*
     * (Pure) virtual member functions
     */

    virtual void runTests();

    virtual void do_binaryRepart();

    virtual void resetInterpolationCache(bool clearCache = false);

    /** \brief calculates back the max/min of the efield on the grid */
    virtual VectorPair_t getEExtrema() = 0;

    virtual double getRho(int x, int y, int z) = 0;

    virtual void computeSelfFields() = 0;

    /** /brief used for self fields with binned distribution */
    virtual void computeSelfFields(int bin) = 0;

    virtual void computeSelfFields_cycl(double gamma) = 0;
    virtual void computeSelfFields_cycl(int bin) = 0;

    virtual void swap(unsigned int i, unsigned int j);

    /*
       Mesh and Field Layout related functions
     */

    virtual void setBCAllPeriodic();
    virtual void setBCAllOpen();

    virtual void setBCForDCBeam();


//     virtual void setMesh(Mesh_t* mesh) = 0;
//     virtual Mesh_t &getMesh() = 0;

//     virtual void setFieldLayout(FieldLayout_t* fLayout) = 0;
    virtual FieldLayout_t &getFieldLayout() = 0;

    /*
     * Wrapped member functions of IpplParticleBase
     */

    size_t getTotalNum() const;
    void setTotalNum(size_t n);
    void setLocalNum(size_t n);
    size_t getLocalNum() const;

    size_t getDestroyNum() const;
    size_t getGhostNum() const;

    unsigned int getMinimumNumberOfParticlesPerCore() const;
    void setMinimumNumberOfParticlesPerCore(unsigned int n);

    ParticleLayout<T, Dim> & getLayout();
    const ParticleLayout<T, Dim>& getLayout() const;

    bool getUpdateFlag(UpdateFlags_t f) const;
    void setUpdateFlag(UpdateFlags_t f, bool val);


    ParticleBConds<Position_t, Dimension>& getBConds() {
        return pbase->getBConds();
    }

    void setBConds(const ParticleBConds<Position_t, Dimension>& bc) {
        pbase->setBConds(bc);
    }

    bool singleInitNode() const;

    void resetID();

    void update();
    void update(const ParticleAttrib<char>& canSwap);

    void createWithID(unsigned id);
    void create(size_t M);
    void globalCreate(size_t np);

    void destroy(size_t M, size_t I, bool doNow = false);
    void performDestroy(bool updateLocalNum = false);
    void ghostDestroy(size_t M, size_t I);

protected:
    size_t calcMoments();    // Calculates bunch moments using only emitted particles.

    /* Calculates bunch moments by summing over bins
     * (not accurate when any particles have been emitted).
     */
    void calcMomentsInitial();
    /// angle range [0~2PI) degree
    double calculateAngle(double x, double y);


private:
    virtual void updateDomainLength(Vektor<int, 3>& grid) = 0;

    virtual void updateFields(const Vector_t& hr, const Vector_t& origin);

    void setup(AbstractParticle<T, Dim>* pb);

public:
    /*
     * Bunch attributes
     */


    ParticlePos_t& R;
    ParticleIndex_t& ID;


    // Particle container attributes
    ParticleAttrib< Vector_t > P;      // particle momentum //  ParticleSpatialLayout<double, 3>::ParticlePos_t P;
    ParticleAttrib< double >   Q;      // charge per simulation particle, unit: C.
    ParticleAttrib< double >   M;      // mass per simulation particle, for multi-species particle tracking, unit:GeV/c^2.
    ParticleAttrib< double >   Phi;    // the electric potential
    ParticleAttrib< Vector_t > Ef;     // e field vector
    ParticleAttrib< Vector_t > Eftmp;  // e field vector for gun simulations

    ParticleAttrib< Vector_t > Bf;    // b field vector
    ParticleAttrib< int >      Bin;   // holds the bin in which the particle is in, if zero particle is marked for deletion
    ParticleAttrib< double >   dt;   // holds the dt timestep for particle

    ParticleAttrib< short >    PType; // we can distinguish dark current particles from primary particle
    ParticleAttrib< int >      TriID; // holds the ID of triangle that the particle hit. Only for BoundaryGeometry case.
    ParticleAttrib< short >    cavityGapCrossed; ///< particle just crossed cavity gap (for ParallelCyclotronTracker)

    ParticleAttrib< short >    bunchNum; // bunch number to which particle belongs (multi-bunch mode)


    Vector_t RefPartR_m;
    Vector_t RefPartP_m;
    ParticleType::type refPType_m;
    CoordinateSystemTrafo toLabTrafo_m;


    /// avoid calls to Ippl::myNode()
    int myNode_m;

    /// avoid calls to Ippl::getNodes()
    int nodes_m;

    /// if the grid does not have to adapt
    bool fixed_grid;

    // The structure for particle binning
    PartBins *pbin_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    // save particles in case of one core
    std::unique_ptr<Inform> pmsg_m;
    std::unique_ptr<std::ofstream> f_stream;

    /// timer for IC, can not be in Distribution.h
    IpplTimings::TimerRef distrReload_m;
    IpplTimings::TimerRef distrCreate_m;

    // For AMTS integrator in OPAL-T
    double dtScInit_m, deltaTau_m;

    /// if a local node has less than 2 particles  lowParticleCount_m == true
    bool lowParticleCount_m;

    /// timer for selfField calculation
    IpplTimings::TimerRef selfFieldTimer_m;

    // get 2nd order momentum matrix
    FMatrix<double, 2 * Dim, 2 * Dim> getSigmaMatrix();

protected:
    IpplTimings::TimerRef boundpTimer_m;
    IpplTimings::TimerRef boundpBoundsTimer_m;
    IpplTimings::TimerRef boundpUpdateTimer_m;
    IpplTimings::TimerRef statParamTimer_m;

    IpplTimings::TimerRef histoTimer_m;


    const PartData *reference;


//     /*
//       Member variables starts here
//     */

    // unit state of PartBunch
    UnitState_t unit_state_;
    UnitState_t stateOfLastBoundP_;

    /// holds the centroid of the beam
    double centroid_m[2 * Dim];

    /// 6x6 matrix of the moments of the beam
    FMatrix<double, 2 * Dim, 2 * Dim> moments_m;

    /// holds the timestep in seconds
    double dt_m;
    /// holds the actual time of the integration
    double t_m;
    /// mean energy of the bunch (MeV)
    double eKin_m;
    /// energy spread of the beam in MeV
    double dE_m;
    /// the position along design trajectory
    double spos_m;

    /// Initialize the translation vector and rotation quaternion
    /// here. Cyclotron tracker will reset these values each timestep
    /// TTracker can just use 0 translation and 0 rotation (quat[1 0 0 0]).
    //Vector_t globalMeanR_m = Vector_t(0.0, 0.0, 0.0);
    //Quaternion_t globalToLocalQuaternion_m = Quaternion_t(1.0, 0.0, 0.0, 0.0);
    Vector_t globalMeanR_m;
    Quaternion_t globalToLocalQuaternion_m;

    /// maximal extend of particles
    Vector_t rmax_m;
    /// minimal extend of particles
    Vector_t rmin_m;

    /// rms beam size (m)
    Vector_t rrms_m;
    /// rms momenta
    Vector_t prms_m;
    /// mean position (m)
    Vector_t rmean_m;
    /// mean momenta
    Vector_t pmean_m;

    /// rms emittance (not normalized)
    Vector_t eps_m;

    /// rms normalized emittance
    Vector_t eps_norm_m;

    Vector_t halo_m;

    /// rms correlation
    Vector_t rprms_m;

    /// dispersion x & y
    double Dx_m;
    double Dy_m;

    /// derivative of the dispersion
    double DDx_m;
    double DDy_m;

    /// meshspacing of cartesian mesh
    Vector_t hr_m;
    /// meshsize of cartesian mesh
    Vektor<int, 3> nr_m;

    /// stores the used field solver
    FieldSolver *fs_m;

    double couplingConstant_m;

    double qi_m;

    /// counter to store the distribution dump
    int distDump_m;

    ///
    int fieldDBGStep_m;

    /// Mesh enlargement
    double dh_m; /// in % how much the mesh is enlarged

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    std::unique_ptr<double[]> bingamma_m;

    //FIXME: this should go into the Bin class!
    // holds number of emitted particles of the bin
    // jjyang: opal-cycl use *nBin_m of pbin_m
    std::unique_ptr<size_t[]> binemitted_m;

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;

    /// current bunch number
    short numBunch_m;

    /// number of particles per bunch
    std::vector<size_t> bunchTotalNum_m;
    std::vector<size_t> bunchLocalNum_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    /*
      Data structure for particle load balance information
    */

    std::unique_ptr<size_t[]> globalPartPerNode_m;


    Distribution *dist_m;

    // flag to tell if we are a DC-beam
    bool dcBeam_m;
    double periodLength_m;
    std::shared_ptr<AbstractParticle<T, Dim> > pbase;
};

#include "PartBunchBase.hpp"

#endif