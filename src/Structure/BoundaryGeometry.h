//
// Implementation of the BoundaryGeometry class
//
// Copyright (c) 200x - 2020, Achim Gsell,
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

#include "AbstractObjects/Definition.h"
#include "Attributes/Attributes.h"
#include "Utilities/Util.h"
#include "Utility/IpplTimings.h"
#include "Utility/PAssert.h"

#include <gsl/gsl_rng.h>

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>

enum class Topology: unsigned short {
    RECTANGULAR,
    BOXCORNER,
    ELLIPTIC
};

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

    int partInside (
        const Vector_t& r,
        const Vector_t& v,
        const double dt,
        Vector_t& intecoords,
        int& triId);

    Inform& printInfo (
        Inform& os) const;

    void writeGeomToVtk (std::string fn);

    inline std::string getFilename() const {
        return Attributes::getString(itsAttr[FGEOM]);
    }

    inline Topology getTopology() const {
        static const std::unordered_map<std::string, Topology> stringTopology_s = {
            {"RECTANGULAR", Topology::RECTANGULAR},
            {"BOXCORNER",   Topology::BOXCORNER},
            {"ELLIPTIC",    Topology::ELLIPTIC}
        };
        Topology topo = stringTopology_s.at(Attributes::getString(itsAttr[TOPO]));
        return topo;
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

    /**
       Return number of boundary faces.
    */
    inline size_t getNumBFaces () {
        return Triangles_m.size();
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

    inline bool getInsidePoint (Vector_t& pt) {
        if (haveInsidePoint_m == false) {
            return false;
        }
        pt = insidePoint_m;
        return true;
    }

    bool findInsidePoint (void);

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
    bool isInside (
        const Vector_t& P                    // [in] point to test
        );

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
    std::vector<std::array<unsigned int,4>> Triangles_m;   // boundary faces defined via point IDs
                                        // please note: 4 is correct, historical reasons!

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

    bool haveInsidePoint_m;
    Vector_t insidePoint_m;             // attribute INSIDEPOINT

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

    inline const Vector_t& getPoint (const int triangle_id, const int vertex_id) {
        PAssert (1 <= vertex_id && vertex_id <=3);
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
        L1,       // in case of BOXCORNER first part of geometry with height B
        L2,       // in case of BOXCORNER second part of geometry with height B-C
        A,        // major semi-axis of elliptic tube
        B,        // minor semi-axis of ellitpic tube
        C,        // in case of BOXCORNER height of corner
        TOPO,     // RECTANGULAR, BOXCORNER, ELLIPTIC if FGEOM is selected TOPO is over-written
        ZSHIFT,   // Shift in z direction
        XYZSCALE, // Multiplicative scaling factor for coordinates
        XSCALE,   // Multiplicative scaling factor for x-coordinates
        YSCALE,   // Multiplicative scaling factor for y-coordinates
        ZSCALE,   // Multiplicative scaling factor for z-coordinates
        INSIDEPOINT,
        SIZE
    };
};

inline Inform &operator<< (Inform& os, const BoundaryGeometry& b) {
    return b.printInfo (os);
}
#endif
