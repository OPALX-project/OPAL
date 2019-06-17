//
//  Copyright & License: See Copyright.readme in src directory
//

/**
   \brief class BoundaryGeometry

   A GEOMETRY definition is used by most physics commands to define the
   particle charge and the reference momentum, together with some other
   data.

   i.e:
   G1: Geometry, FILE="input.h5"
   G2: Geometry, L=1.0, A=0.0025, B=0.0001

   :TODO: update above section
 */

#ifndef _OPAL_BOUNDARY_GEOMETRY_H
#define _OPAL_BOUNDARY_GEOMETRY_H

class OpalBeamline;
class ElementBase;

#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <array>

#include "AbstractObjects/Definition.h"
#include "Attributes/Attributes.h"
#include "Structure/SecondaryEmissionPhysics.h"

#include <gsl/gsl_rng.h>

extern Inform* gmsg;

namespace BGphysics {
    enum TPHYACTION {
        Nop = 0x01,                 // triangle is transparent to particle like beam window
        Absorption = 0x02,          // triangle has no field and secondary emission
        FNEmission = 0x04,          // triangle has field emission
        SecondaryEmission = 0x08    // triangle has secondary emission
    };
}

class BoundaryGeometry : public Definition {

public:
    BoundaryGeometry();
    virtual ~BoundaryGeometry();

    virtual bool canReplaceBy (
        Object* object);

    virtual BoundaryGeometry* clone (
        const std::string& name);

    // Check the GEOMETRY data.
    virtual void execute ();

    // Find named GEOMETRY.
    static BoundaryGeometry* find (
        const std::string& name);

    // Update the GEOMETRY data.
    virtual void update ();

    void updateElement (
        ElementBase* element);

    void initialize ();

    void createParticlesOnSurface (
        size_t n, double darkinward,
        OpalBeamline& itsOpalBeamline,
        PartBunchBase<double, 3>* itsBunch);

    void createPriPart (
        size_t n, double darkinward,
        OpalBeamline& itsOpalBeamline,
        PartBunchBase<double, 3>* itsBunch);

    int partInside (
        const Vector_t& r,
        const Vector_t& v,
        const double dt,
        int Parttype,
        const double Qloss,
        Vector_t& intecoords,
        int& triId);

    // non secondary emission version.
    int emitSecondaryNone (
        const Vector_t& intecoords,
        const int& triId);

    // call Furman-Pivi's model
    int emitSecondaryFurmanPivi (
        const Vector_t& intecoords,
        const int i,
        PartBunchBase<double, 3>* itsBunch,
        double& seyNum);

    // call Vaughan's model
    int emitSecondaryVaughan (
        const Vector_t& intecoords,
        const int i,
        PartBunchBase<double, 3>* itsBunch,
        double& seyNum);

    Inform& printInfo (
        Inform& os) const;

    void writeGeomToVtk (std::string fn);

    inline std::string getFilename () const {
        return (std::string) Attributes::getString (itsAttr[FGEOM]);
    }

    inline std::string getTopology () const {
        return (std::string) Attributes::getString (itsAttr[TOPO]);
    }

    inline std::string getDistribution () {
        return (std::string) Attributes::getString (itsAttr[DISTR]);
    }

    inline std::vector<std::string> getDistributionArray () {
        return Attributes::getStringArray (itsAttr[DISTRS]);
    }

    inline size_t getN () {
        return partsr_m.size ();
    }

    inline Vector_t getCooridinate (size_t i) {
        return partsr_m[i];
    }
    inline void clearCooridinateArray () {
        return partsr_m.clear ();
    }
    inline Vector_t getMomenta (size_t i) {
        return partsp_m[i];
    }

    inline void clearMomentaArray () {
        return partsp_m.clear ();
    }

    inline double getA() {
        return (double)Attributes::getReal(itsAttr[A]);
    }

    inline double getB() {
        return (double)Attributes::getReal(itsAttr[B]);
    }

    inline double getC() {
        return (double)Attributes::getReal(itsAttr[C]);
    }

    inline double getS() {
        return (double)Attributes::getReal(itsAttr[S]);
    }

    inline double getLength() {
        return (double)Attributes::getReal(itsAttr[LENGTH]);
    }

    inline double getL1() {
        return (double)Attributes::getReal(itsAttr[L1]);
    }

    inline double getL2() {
        return (double)Attributes::getReal(itsAttr[L2]);
    }

    inline void setNEmissionMode (bool nEmissionMode) {
        nEmissionMode_m = nEmissionMode;
    }

    inline void setWorkFunction (double workFunction) {
        workFunction_m = workFunction;
    }

    inline void setFieldEnhancement (double fieldEnhancement) {
        fieldEnhancement_m = fieldEnhancement;
    }

    inline void setMaxFN (size_t maxFNemission) {
        maxFNemission_m = maxFNemission;
    }

    inline void setFNTreshold (double fieldFNthreshold) {
        fieldFNthreshold_m = - 1.0e6 * fieldFNthreshold;
    }

    inline void setFNParameterA (double parameterFNA) {
        parameterFNA_m = parameterFNA;
    }

    inline void setFNParameterB (double parameterFNB) {
        parameterFNB_m = parameterFNB;
    }

    inline void setFNParameterY (double parameterFNY) {
        parameterFNY_m = parameterFNY;
    }

    inline void setFNParameterVYZe (double parameterFNVYZe) {
        parameterFNVYZe_m = parameterFNVYZe;
    }

    inline void setFNParameterVYSe (double parameterFNVYSe) {
        parameterFNVYSe_m = parameterFNVYSe;
    }

    inline void setBoundaryMatType (int BoundaryMatType) {
        seBoundaryMatType_m = BoundaryMatType;
    }

    inline void setEInitThreshold (double einitthreshold) {
        eInitThreshold_m = 1.0e6 * einitthreshold;
    }

    // return sey_0 in Vaughan's model
    inline void setvSeyZero (double vSeyZero) {
        vSeyZero_m = vSeyZero;
    }

    // set the energy related to sey_0 in Vaughan's model
    inline void setvEZero (double vEZero) {
        vEzero_m = vEZero;
    }

    // set sey max in Vaughan's model
    inline void setvSeyMax (double vSeyMax) {
        vSeyMax_m = vSeyMax;
    }

    // return Emax in Vaughan's model
    inline void setvEmax (double vEmax) {
        vEmax_m = vEmax;
    }

    // return fitting parameter denotes the roughness of surface for
    // impact energy in Vaughan's model
    inline void setvKenergy (double vKenergy) {
        vKenergy_m = vKenergy;
    }

    // return fitting parameter denotes the roughness of surface for impact
    // angle in Vaughan's model
    inline void setvKtheta (double vKtheta) {
        vKtheta_m = vKtheta;
    }

    // return thermal velocity Maxwellian distribution of secondaries
    // in Vaughan's model
    inline void setvVThermal (double vVThermal) {
        vVThermal_m = vVThermal;
    }

    inline void setVw (double ppVw) {
        ppVw_m = ppVw;
    }

    /**
       Return number of boundary faces.
    */
    inline int getNumBFaces () {
        return numTriangles_m;
    }

    /**
       Return the hr_m.
    */
    inline Vector_t gethr () {
        return voxelMesh_m.sizeOfVoxel;
    }
    /**
       Return the nr_m.
     */
    inline Vektor<int, 3> getnr () {
        return voxelMesh_m.nr_m;
    }

    /**
       Return the mincoords_m.
     */
    inline Vector_t getmincoords () {
        return minExtent_m;
    }
    /**
       Return the maxcoords_m.
    */
    inline Vector_t getmaxcoords () {
        return maxExtent_m;
    }

    /**
       @param  TriBarycenters_m store the coordinates of barycentric points of
       triangles, The Id number is the same with triangle Id.
    */
    std::vector<Vector_t> TriBarycenters_m;

    /**
       @param TriPrPartloss_m[i]:
       cummulative sum of primary particles charge hitting triangle 'i'
    */
    std::vector<double> TriPrPartloss_m;

    /**
       @param TriSePartloss_m store the number of secondary particles hitting the
       Id th triangle. The Id number is the same with triangle Id(not vertex ID).
    */
    std::vector<double> TriSePartloss_m;

    /**
       @param TriFEPartloss_m store the number of field emission/darkcurrent
       particles hitting the Id th triangle. The Id number is the same with
       triangle Id(not vertex ID).
    */
    std::vector<double> TriFEPartloss_m;

    /**
       @param TriBGphysicstag_m store the tags of each boundary triangle for
       proper physics action.
    */
    std::vector<short> TriBGphysicstag_m;

    inline bool isOutsideApperture(Vector_t x) {
        if (hasApperture()) {
            for (unsigned int i=0; i<apert_m.size(); i=i+3) {
                if ((apert_m[i] <= x(2)) && (x(2) < apert_m[i+1])) {
                    // yes we are inside the interval
                    const double r = apert_m[i+2] * apert_m[i+2];
                    return ((x(0)*x(0)) + (x(1)*x(1))) > r;
                }
            }
        }
        return false;
    }

    int intersectRayBoundary (
        const Vector_t& P,
        const Vector_t& v,
        Vector_t& I);

    int fastIsInside (
        const Vector_t& reference_pt,        // [in] a reference point
        const Vector_t& P                    // [in] point to test
        );

    enum DebugFlags {
        debug_isInside                         = 0x0001,
        debug_fastIsInside                     = 0x0002,
        debug_intersectRayBoundary             = 0x0004,
        debug_intersectLineSegmentBoundary     = 0x0008,
        debug_intersectTinyLineSegmentBoundary = 0x0010,
        debug_PartInside                       = 0x0020,
    };

    inline void enableDebug(enum DebugFlags flags) {
        debugFlags_m |= flags;
    }

    inline void disableDebug(enum DebugFlags flags) {
        debugFlags_m &= ~flags;
    }

private:
    int intersectTriangleVoxel (
        const int triangle_id,
        const int i,
        const int j,
        const int k);

    int intersectTinyLineSegmentBoundary (
        const Vector_t&,
        const Vector_t&,
        Vector_t&,
        int&
        );

    int intersectLineSegmentBoundary (
        const Vector_t& P0,
        const Vector_t& P1,
        Vector_t& intersection_pt,
        int& triangle_id
        );

    std::string h5FileName_m;           // H5hut filename

    std::vector<Vector_t> Points_m;     // geometry point coordinates
    std::vector<std::array<int,4>> Triangles_m;   // boundary faces defined via point IDs
                                        // please note: 4 is correct, historical reasons!
    int numTriangles_m;                 // number of boundary triangles

    std::vector<Vector_t> TriNormals_m; // oriented normal vector of triangles
    std::vector<double> TriAreas_m;     // area of triangles

    Vector_t minExtent_m;               // minimum of geometry coordinate.
    Vector_t maxExtent_m;               // maximum of geometry coordinate.

    struct {
        Vector_t minExtent;
        Vector_t maxExtent;
        Vector_t sizeOfVoxel;
        Vektor<int, 3> nr_m;            // number of intervals of geometry in X,Y,Z direction
        std::unordered_map<int,         // map voxel IDs ->
            std::unordered_set<int>> ids; // intersecting triangles

    } voxelMesh_m;

    int debugFlags_m;

    std::vector<Vector_t> partsp_m;     // particle momenta
    std::vector<Vector_t> partsr_m;     // particle positions

    /*
       An additional structure to hold apperture information
       to prevent that particles go past the geometry. The user
       can specify n trippel with the form: (zmin, zmax, r)
    */
    std::vector<double> apert_m;

    SecondaryEmissionPhysics sec_phys_m;

    bool nEmissionMode_m;
    double eInitThreshold_m;

    // Vaughan's model
    double vSeyZero_m;          // energy related to sey_
    double vEzero_m;            // sey_0
    double vSeyMax_m;           // sey max
    double vEmax_m;             // Emax
    double vKenergy_m;          // roughness of surface for impact energy
    double vKtheta_m;           // roughness of surface for impact angle
    double vVThermal_m;         // thermal velocity of Maxwellian distribution of secondaries

    double ppVw_m;              // velocity scalar for Parallel plate benchmark.
    int seBoundaryMatType_m;    // user defined material type for secondary emission model.

    // Fowler-Nordheim model
    double workFunction_m;      // work function
    double fieldEnhancement_m;  // field factor
    size_t maxFNemission_m;     // maximum emitted number per triangle
    double fieldFNthreshold_m;  // lower threshold electric field
    double parameterFNA_m;      // parameter A. Default: \f$1.54\times 10^{-6}\f$.
    double parameterFNB_m;      // parameter B. Default: \f$6.83\times 10^9\f$.
    double parameterFNY_m;      // parameter Y. Default:\f$3.795\times 10^{-5}\f$.
    double parameterFNVYZe_m;   // zero order fit constant for v(y). Default:\f$0.9632\f$.
    double parameterFNVYSe_m;   // second order fit constant for v(y). Default:\f$1.065\f$.

    gsl_rng *randGen_m;         //

    IpplTimings::TimerRef Tinitialize_m; // initialize geometry
    IpplTimings::TimerRef TisInside_m;
    IpplTimings::TimerRef TfastIsInside_m;
    IpplTimings::TimerRef TRayTrace_m;   // ray tracing
    IpplTimings::TimerRef TPartInside_m; // particle inside

    BoundaryGeometry(const BoundaryGeometry&);
    void operator= (const BoundaryGeometry&);

    // Clone constructor.
    BoundaryGeometry(const std::string& name, BoundaryGeometry* parent);

    inline bool hasApperture() {
        return (apert_m.size() != 0);
    }

    inline const Vector_t& getPoint (const int triangle_id, const int vertex_id) {
        assert (1 <= vertex_id && vertex_id <=3);
        return Points_m[Triangles_m[triangle_id][vertex_id]];
    }

    enum INTERSECTION_TESTS {
        SEGMENT,
        RAY,
        LINE
    };

    int intersectLineTriangle (
        const enum INTERSECTION_TESTS kind,
        const Vector_t& P0,
        const Vector_t& P1,
        const int triangle_id,
        Vector_t& I);

    inline int mapVoxelIndices2ID (const int i, const int j, const int k);
    inline Vector_t mapIndices2Voxel (const int, const int, const int);
    inline Vector_t mapPoint2Voxel (const Vector_t&);
    inline void computeMeshVoxelization (void);

    enum {
        FGEOM,    // file holding the geometry
        LENGTH,   // length of elliptic tube or boxcorner
        S,        // start of the geometry
        L1,       // in case of BOXCORNER first part of geometry with hight B
        L2,       // in case of BOXCORNER second part of geometry with hight B-C
        A,        // major semi-axis of elliptic tube
        B,        // minor semi-axis of ellitpic tube
        C,        // in case of BOXCORNER hight of corner
        TOPO,     // BOX, BOXCORNER, ELLIPTIC if FGEOM is selected topo is over-written
        DISTR,    // Add distribution to generate physics model on the surface
        DISTRS,   // Add distribution array to generate physics model on the surface
        ZSHIFT,   // Shift in z direction
        XYZSCALE, // Multiplicative scaling factor for coordinates
        XSCALE,   // Multiplicative scaling factor for x-coordinates
        YSCALE,   // Multiplicative scaling factor for y-coordinates
        ZSCALE,   // Multiplicative scaling factor for z-coordinates
        APERTURE,    // in addition to the geometry
        SIZE
    };
};

inline Inform &operator<< (Inform& os, const BoundaryGeometry& b) {
    return b.printInfo (os);
}
#endif
// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
