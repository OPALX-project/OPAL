/*
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _CLASSIC_FIELDS_SECTORMAGNETICFIELDMAP_HH_
#define _CLASSIC_FIELDS_SECTORMAGNETICFIELDMAP_HH_

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Fields/SectorField.h"

namespace interpolation {
    class VectorMap;
    class ThreeDGrid;
}

/** \class[SectorMagneticFieldMap] 
 *
 *  \brief handles field map grids with sector geometry
 *
 *  SectorMagneticFieldMap provides an interface to the 3D interpolator
 *  routines for 3D field maps in a sector geometry. Interpolation is done from
 *  points in (r, y, phi) geometry off of a 3D field map stored in cartesian
 *  coordinates.
 *  - position coordinates are left in polar because it makes the interpolation
 *    easier
 *  - field coordinates are left in cartesian because it means we can avoid
 *    doing a coordinate transformation (tracking is done in cartesian
 *    coordinates)
 *
 *  A constructor is provided which calls the SectorMagneticFieldMapIO class to
 *  generate the field map, and it is expected that this will be the usual way
 *  to generate the map.
 *
 *  The SectorMagneticFieldMap class caches field maps by file name to avoid
 *  multiple loads of the same field map. This cache should be cleared after
 *  running (otherwise we have a memory leak).
 *
 *  Note that the coordinate system used for tracking by Ring (z is up) is a
 *  bit different to that implemented in SectorField (y is up). I perform a
 *  coordinate transformation in the getFieldstrength function to facilitate
 *  this. Note that OpalBeamline has y is up...
 *
 *  SectorMagneticFieldMap enforces that field is hard against the start edge
 *  with no gap. SectorMagneticFieldMap auto-detects any field phi offset and
 *  removes each time it does a field lookup (just an add operation, so fast).
 *  If a phi offset is desired, user must zero-pad field map.
 */
class SectorMagneticFieldMap : public SectorField {
  public:
    /** Generate the field map by calling SectorMagneticFieldMap::IO
     *
     *  \param file_name name of the file to read in
     *  \param symmetry parameter that controls mirror symmetry applied to the
     *         field map - either "Dipole" or "None"
     *  \param length_units multiplier for lengths
     *  \param field_units multiplier for fields
     *  \param polynomial_order order of polynomial fit to the mesh points
     *  \param smoothing_order order of smoothing (should be >= polynomial_order)
     */
    SectorMagneticFieldMap(std::string file_name,
                           std::string symmetry,
                           double length_units,
                           double field_units,
                           int polynomial_order,
                           int smoothing_order);

    /** Copy constructor - makes a deep copy of the field map */
    SectorMagneticFieldMap(const SectorMagneticFieldMap& field);

    /** Destructor - delete allocated memory */
    ~SectorMagneticFieldMap();

    /** Return the field value in polar coordinates
     *
     *  \param R_p position in cylindrical coordinates at which to
     *         evaluate the field. Should be an array of length 3 like
     *           (r, y, phi)
     *  \param E_p reference to an allocated 3-vector.
     *         The function will fill the block with the value of the field in
     *         cylindrical polar coordinates,
     *           (e_r, e_y, e_phi)
     *         Overwrites any existing data
     *  \param B_p reference to an allocated 3-vector.
     *         The function will fill the block with the value of the field in
     *         cylindrical polar coordinates,
     *           (b_r, b_y, b_phi)
     *         Overwrites any existing data
     *  \returns false if R_p is inside the bounding box
     */
    /* bool getFieldstrengthPolar */
    /*               (const Vector_t &R_p, Vector_t &E_p, Vector_t &B_p) const; */

    /** Get the field value in cartesian coordinates
     *
     *  \param R_p array holding the cartesian 3-position at which the field
     *         is to be evaluated, (x, y, z)
     *  \param E_p array of at 3 doubles, to which the field at point
     *         is written (ex, ey, ez) - note this field is magnetostatic so
     *         always return 0 here.
     *  \param B_p array of at 3 doubles, to which the field at point
     *         is written (bx, by, bz)
     *  \returns false if R_c is inside the bounding box
     */
    bool getFieldstrength
                  (const Vector_t &R_c, Vector_t &E_c, Vector_t &B_c) const;

    /** Get a pointer to the interpolator or nullptr if it is not set
     *
     *  Note SectorMagneticFieldMap still owns this memory.
     */
    interpolation::VectorMap* getInterpolator();

    /** Set the interpolator
     *
     *  \param interpolator set the field map interpolator. Note
     *  SectorMagneticFieldMap now owns the memory pointed to by interpolator
     */
    void setInterpolator(interpolation::VectorMap* interpolator);


    /** Get the field map file name */
    std::string getFieldMapFileName() const {return filename_m;}

    /** Get a string corresponding to the field map symmetry
     *
     *  \returns field map symmetry, either "None" or "Dipole"
     */
    std::string getSymmetry() const;

    /** Set the field map symmetry
     *
     *  \param name field map symmetry, either "None" or "Dipole"
     */
    void setSymmetry(std::string name);

    /** Delete cached fields
     */
    static void clearFieldCache();

    /** Not implemented - throws a LogicalError */
    bool getFieldDerivative(const Vector_t &R, Vector_t &E,
                                  Vector_t &B, const DiffDirection &dir) const;

    /** Does nothing */
    void swap();

    /** Print summary information about the field map to Inform */
    void getInfo(Inform* msg);

    /** Print summary information about the field map to out */
    void print(std::ostream& out);

    /** Read in the field map from the file */
    void readMap();

    /** Delete the field map interpolator and set pointer to nullptr */
    void freeMap();

    /** Magnetostatic field map - so returns 0. */
    double getFrequency() const {return 0;}

    /** Magnetostatic field map - so does nothing */
    void setFrequency(double) {}

    /** Get phi offset */
    double getPhiOffset() const {return phiOffset_m;}

    /** Set phi offset */
    void setPhiOffset(double dphi) {phiOffset_m = dphi;}

    /** Get change in azimuthal angle between entrance and exit */
    double getDeltaPhi() const;

    friend class FieldMap;

  private:
    enum symmetry {none, dipole};

    /** Reflect R_temp about y if below bbmin 
     *  \returns true if symmetry transformation was applied
     */
    bool applySymmetry(double* R_temp) const;

    static symmetry StringToSymmetry(std::string name);
    static std::string SymmetryToString(symmetry sym);

    void Rotate(double* value, double angle);

    interpolation::VectorMap* interpolator_m;
    symmetry symmetry_m;
    std::vector<double> units_m;
    const std::string filename_m;
    double phiOffset_m;
    int poly_order_m;
    int smoothing_order_m;

    static const double fractionalBBPhiTolerance_m;
    static std::map<std::string, SectorMagneticFieldMap*> _fields;
    friend class SectorMagneticFieldMapIO;

    // disabled - dont use
    SectorMagneticFieldMap& operator=(const SectorMagneticFieldMap& field);

    class IO;
};


/** \class[SectorMagneticFieldMap::IO]
 *
 *  \brief handles reading sector field maps
 *
 *  SectorMagneticFieldMap::IO provides routines to read a sector field map for
 *  input to tracking.
 */
class SectorMagneticFieldMap::IO {
  public:
    /** Read in the field map
     *
     *  Read in the field map with some specified geometry
     *
     *  \param file_name name of the file containing the field map
     *  \param units units of the file - should be 6-vector
     *  \param sym symmetry of the file - either "none" (full field map) or
     *         "dipole" (field map is reflected about y = 0)
     *  \param poly_order order of the polynomial fit
     *  \param smoothing_order order of the polynomial smoothing
     */
    static interpolation::VectorMap* readMap(std::string file_name,
                                           std::vector<double> units,
                                           SectorMagneticFieldMap::symmetry sym,
                                           int poly_order,
                                           int smoothing_order);

  private:
    static const double floatTolerance_m;
    static const int sortOrder_m[3];

    // read and sort the lines of the map file
    static std::vector< std::vector<double> > readLines
                             (std::string file_name, std::vector<double> units);

    // generate a grid based on the input map file
    static interpolation::ThreeDGrid* generateGrid(
                         const std::vector< std::vector<double> > field_points,
                         SectorMagneticFieldMap::symmetry sym);

    // get the interpolator based on field points and grid information
    static interpolation::VectorMap* getInterpolator(
                         const std::vector< std::vector<double> > field_points,
                         interpolation::ThreeDGrid* grid,
                         SectorMagneticFieldMap::symmetry sym);

    static interpolation::VectorMap* getInterpolatorPolyPatch(
                         const std::vector< std::vector<double> > field_points,
                         interpolation::ThreeDGrid* grid,
                         SectorMagneticFieldMap::symmetry sym,
                         int poly_order,
                         int smoothing_order);

    // Compare two floats are same within tolerance
    static bool floatGreaterEqual(double in1, double in2);

    // comparator for sorting field_points by r, y, phi. Sort order is given by
    // sort_order variable
    static bool comparator
             (std::vector<double> field_item1, std::vector<double> field_item2);

    // private constructor i.e. disabled
    IO();
    // private copy constructor i.e. disabled
    IO(const IO& map);

    // private destructor i.e. disabled
    ~IO();

    static const std::string errMsg1;
};

#endif
