#ifndef CLASSIC_OutputPlane_HH
#define CLASSIC_OutputPlane_HH

#include "AbsBeamline/PluginElement.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"

#include <memory>
#include <string>
#include <sstream>
#include <gsl/gsl_rng.h>

class Component;
class CollimatorPhysics;

/** Class OutputPlane
 *  ------------------------------------------------------------------------
 *  Class OutputPlane defines an output plane, which enables output of tracking
 *  data at an arbitrary spatial plane. The crossing of the plane in a given time
 *  step is detected by (1) checking if the step is within c*dt of the plane and
 *  then (2) checking if an RK4 step crosses the plane. 
 * 
 *  If the track crosses, then we estimate the position of the track either by
 *  doing a linear interpolation between the initial position and the position
 *  of the test step; OR we do an iterative RK4 step, reducing the RK4 step size
 *  each time until we find crossing position to some tolerance or the maximum
 *  number of iterations is reached.
 * 
 *  Additional features - we can recentre the plane on a particular track (for
 *  example, one might seek to make a plane perpendicular to a reference 
 *  trajectory); we can limit the extent of the plane so only particles crossing
 *  within a certain distance are registered; we can apply material physics model
 *  to tracks crossing the OutputPlane.
 * 
 *  BUG: if the "Global" Opal tracking is not RK4, then there is a small 
 *  possibility that the "Global" track and the OutputPlane step will be not
 *  perfectly aligned; tracks could register crossing twice or not crossing at
 *  all. An improved algorithm would be to use the "global" stepper or make some
 *  hack to consider e.g. two steps.
 *  BUG: verticalExtent_m and radialExtent_m don't appear to be used/implemented.
 */

class OutputPlane: public PluginElement {

public:
    enum algorithm {RK4, INTERPOLATION};

    /** Constructor with given name. */
    explicit OutputPlane(const std::string &name);

    /** Default constructor (initialise to empty */
    OutputPlane();

    /** Copy constructor */
    OutputPlane(const OutputPlane &);

    /** Destructor */
    virtual ~OutputPlane();

    /** Inheritable copy operation */
    ElementBase* clone() const override;
    
    /** Apply visitor to output plane. */
    virtual void accept(BeamlineVisitor &) const override;

    /** Get the field map */
    inline Component* getGlobalFieldMap() const;
    /** Set the field map */
    inline void setGlobalFieldMap(Component* field);

    /** Get the normal to the plane */
    inline Vector_t getNormal() const;
    /** Set the normal to the plane */
    inline void setNormal(Vector_t normal);

    /** Get the centre of the plane */
    inline Vector_t getCentre() const;
    /** Set the centre of the plane */
    inline void setCentre(Vector_t centre);

    /** Get the tolerance, used when finding intersection with the plane */
    inline double getTolerance() const;
    /** Set the tolerance, used when finding intersection with the plane */
    inline void setTolerance(double tolerance);

    /** Get the maximum allowed iteratiosn when finding intersection */
    inline size_t getMaxIterations() const;
    /** Set the maximum allowed iteratiosn when finding intersection */
    inline void setMaxIterations(size_t max);

    /** Return the index of the recentring particle*/
    inline int getRecentre() const;
    /** Set the index of the recentring particle*/
    inline void setRecentre(int willRecentre);

    /** Get the radial extent from the plane centre */
    inline double getRadialExtent() const;
    /** Set the radial extent from the centre */
    inline void setRadialExtent(double r);

    /** Get the vertical extent from the plane centre
     */
    inline double getVerticalExtent() const;
    /** Set the vertical extent from the plane centre */
    inline void setVerticalExtent(double z);

    /** Get the horizontal extent from the plane centre
     */
    inline double getHorizontalExtent() const;
    /** Set the horozontal extent from the plane centre */
    inline void setHorizontalExtent(double width);

    /** Get the algorithm used to interpolate from step to the plane */
    inline algorithm getAlgorithm() const;
    /** Set the algorithm used to interpolate from step to the plane */
    inline void setAlgorithm(algorithm alg);

    /** Set up the material model */
    inline void setMaterial(std::shared_ptr<CollimatorPhysics> material,
                            double thickness);


    /** Check for plane crossing. If crossed fill t, R, P with the intercept.
     *
     *  @li index: Particle number of the particle crossing
     *  @li tstep: Time step used by the global tracking
     *  @li chargeToMass: charge to mass ratio of the particle, in units of
     *                    (positron charge)/GeV/c^2
     *  @li t: time of the initial particle [ns]. Unchanged if no crossing is
     *         found, else filled with the estimated time of the particle at the
     *         plane crossing point.
     *  @li R: position three-vector of the initial particle [m]. Unchanged if 
     *         no crossing is found, else filled with the estimated position of
     *         the particle at the plane crossing point.
     *  @li P: (gamma beta) three-vector of the initial particle [].
     *         Unchanged if no crossing is found, else filled with the estimated
     *         position of the particle at the plane crossing point.
     *
     *  The routine checks for crossing by performing a full RK4 step using the
     *  user-supplied time step. The distance to the plane is calculated. If the
     *  distance changes sign, the particle is assumed to cross; if the distance
     *  does not change sign, the particle is assumed to stay on the same side.
     *  Note that particles that step across the boundary and back again will
     *  not be registered as a crossing.
     * 
     *  The intercept is found by iteratively doing RK4 steps. Successive step
     *  sizes are found by linear interpolation based on the distance from the 
     *  plane and the time step.
     *
     *  If recentre_m is positive then centre_m is set to the position of 
     *  the n^th crossing particle and normal_m is set to be parallel to the 
     *  momentum vector. recentre_m is set to -1. The idea is to generate
     *  output planes in the coordinate system of some trajectory (e.g. for
     *  studying beam transport/transfer matrix). Normally this will be the
     *  reference trajectory but it is possible to make it another trajectory if
     *  required.
     *
     *  If radialExtent_m is > 0, crossings more than radialExtent_m from the
     *  centre are ignored. If horizontalExtent_m is > 0, crossings more than
     *  horizontalExtent_m from the centre, in global x-y plane, are ignored. If
     *  verticalExtent_m is > 0, crossings more than verticalExtent_m in 
     *  direction parallel to global z axis, are ignored.
     *
     *  @returns True if the particle crosses the plane; else false.
     */
    bool checkOne(const int index, const double tstep, double chargeToMass,
                  double& t, Vector_t& R, Vector_t& P);

    /** Returns empty field */
    NullField &getField() {return nullfield_m;}
    /** Returns empty field */
    const NullField &getField() const {return nullfield_m;}

    /** Returns empty geometry */
    StraightGeometry &getGeometry() {return geom_m;}
    /** Returns empty geometry */
    const StraightGeometry &getGeometry() const {return geom_m;}
    /** Make an RK4Step */
    void RK4Step(const double& tstep,
                 const double& chargeToMass,
                 const double& t,
                 Vector_t& R,
                 Vector_t& P) const;

    void operator=(const OutputPlane &) = delete;
    
    //ElementBase::ElementType getType() const;
    ElementType getType() const;
    void recentre(Vector_t R, Vector_t P);


private:
    /// Initialise peakfinder file
    virtual void doInitialise(PartBunchBase<double, 3>* /*bunch*/);
    //virtual void doInitialise() ;

    /// Record probe hits when bunch particles pass
    virtual bool doCheck(PartBunchBase<double, 3> *bunch,
                         const int turnnumber, 
                         const double t, 
                         const double tstep) override;
    /// Hook for goOffline
    virtual void doGoOffline() override;
    void getDerivatives(const Vector_t& R,
                        const Vector_t& P,
                        const double& t,
                        const double& chargeToMass,
                        double* yp) const;
    double distanceToPlane(Vector_t point) const;
    void rk4Test(double tstep, double chargeToMass,
                 double& t, Vector_t& R, Vector_t& P);
    void interpolation(double& t, Vector_t& R, Vector_t& P);

    //void applyMaterialPhysics(Vector_t& R, Vector_t& P) const;
    
    void hardCodeMaterialPhysicsGraphite(Vector_t& P) const;
    Component* field_m = NULL; // field map - this is a borrowed pointer
    Vector_t normal_m; // normal to the output plane
    Vector_t centre_m; // centre of the output plane
    double radialExtent_m = -1.0; // maximum radial extent of the plane (circular)
    double verticalExtent_m = -1.0; // maximum vertical extent of the plane (rectangular)
    double horizontalExtent_m = -1.0; // maximum horizontal extent of the plane (rectangular)
    size_t maxIterations_m = 10; // maximum number of iterations when finding intercept
    double tolerance_m = 1e-9; // tolerance on distance from plane when finding intercept
    const double cLight_m; // m/ns;
    const double chargeToMassUnits_m; // m/ns;
    const double positionUnits_m; // mm;
    NullField nullfield_m; // dummy variable for inheritance
    StraightGeometry geom_m; // dummy variable for inheritance
    algorithm algorithm_m; // sets interpolation or RK4
    int recentre_m; // particle index for recentring
    int nHits_m = 0; // counter for number of hits on the plane
    int verbose_m = 0; // set to 0 for silent running; > 10 for verbose when step detected; > 20 for verbose when step not detected
    std::stringstream output_m;
    // mock up of a foil
    std::shared_ptr<CollimatorPhysics> material_m;
    double materialThickness_m = -1.0;
    static gsl_rng* rng_m; // for scattering
};

Component* OutputPlane::getGlobalFieldMap() const {
    return field_m;
}

void OutputPlane::setGlobalFieldMap(Component* field) {
    field_m = field;
}

Vector_t OutputPlane::getNormal() const {
    return normal_m;
}

void OutputPlane::setNormal(Vector_t normal) {
    normal_m = normal;
    normal_m /= euclidean_norm(normal);
}

Vector_t OutputPlane::getCentre() const {
    return centre_m;
}

void OutputPlane::setCentre(Vector_t centre) {
    centre_m = centre;
}

double OutputPlane::getTolerance() const {
    return tolerance_m;
}

void OutputPlane::setTolerance(double tolerance) {
    tolerance_m = tolerance;
}

size_t OutputPlane::getMaxIterations() const {
    return maxIterations_m;
}

void OutputPlane::setMaxIterations(size_t max) {
    maxIterations_m = max;
}

double OutputPlane::getHorizontalExtent() const {
    return horizontalExtent_m;
}

void OutputPlane::setHorizontalExtent(double width) {
    horizontalExtent_m = width;
}

double OutputPlane::getVerticalExtent() const {
    return verticalExtent_m;
}

void OutputPlane::setVerticalExtent(double width) {
    verticalExtent_m = width;
}

double OutputPlane::getRadialExtent() const {
    return radialExtent_m;
}

void OutputPlane::setRadialExtent(double radius) {
    radialExtent_m = radius;
}

int OutputPlane::getRecentre() const {
    return recentre_m;
}

void OutputPlane::setRecentre(int recentre) {
    recentre_m = recentre;
}

OutputPlane::algorithm OutputPlane::getAlgorithm() const {
    return algorithm_m;
}

void OutputPlane::setAlgorithm(OutputPlane::algorithm alg) {
    algorithm_m = alg;
}

void OutputPlane::setMaterial(std::shared_ptr<CollimatorPhysics> material, double thickness) {
    material_m = material;
    materialThickness_m = thickness;
}


#endif // CLASSIC_OutputPlane_HH
