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

#ifndef _CLASSIC_FIELDS_SECTORFIELD_HH_
#define _CLASSIC_FIELDS_SECTORFIELD_HH_

#include <vector>

#include "Fields/Fieldmap.h"

/** \class SectorField
 *
 *  SectorField is an abstraction type for a sector field map i.e. a field map
 *  with a bent geometry as in a ring. Functions are provided for converting
 *  from cartesian coordinate systems to ring coordinate systems and vice versa
 *  and field calculation functions are defined for the field in cartesian
 *  coordinates and additionally ring coordinates.
 *
 *  The sector field map holds two bounding boxes - the standard rectangular
 *  bounding box in cartesian coordinates and additionally a sector bounding box
 *  in polar coordinates. The bounding box should always be given in polar
 *  coordinates - SectorField will do the conversion to cartesian coordinates
 *  (note this is a lossy procedure)
 *
 *  Ring coordinates are typically given as three vectors going like
 *  (radius, y, angle) also denoted by (r, y, phi). The ring coordinate system
 *  is defined with (x, y, z) = (0, 0, 0) the centre of the ring; y is the same
 *  in ring and cartesian coordinate systems. Positive phi is anticlockwise in
 *  the (x, z) plane.
 */
class SectorField {
  public:
    /** Make an empty sector field; set bounding box to max double
     */
    SectorField(const std::string& file_name);

    /** Destructor (does nothing)
     */
    virtual ~SectorField();

    /** Return the field value in polar coordinates
     *
     *  \param R_p position in cylindrical coordinates at which to
     *         evaluate the field. Should be an array of length 4 like
     *           (r, y, phi, t)
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
     *  \returns true if any field value is non-zero
     */
    /* virtual bool getFieldstrengthPolar */
    /*               (const Vector_t &R_p, Vector_t &E_p, Vector_t &B_p) const = 0; */

    /** Return the field value in cartesian coordinates
     *
     *  \param R_c position in cartesian coordinates at which to
     *         evaluate the field. Should be an array of length 4 like
     *           (x, y, z, t)
     *  \param E_p reference to an allocated 3-vector.
     *         The function will fill the block with the value of the field in
     *         cartesian coordinates,
     *           (e_x, e_y, e_z)
     *         Overwrites any existing data
     *  \param B_p reference to an allocated 3-vector.
     *         The function will fill the block with the value of the field in
     *         cartesian coordinates,
     *           (b_x, b_y, b_z)
     *         Overwrites any existing data
     *  \returns true if any field value is non-zero
     */
    virtual bool getFieldstrength
                  (const Vector_t &R_c, Vector_t &E_c, Vector_t &B_c) const = 0;

    /** Convert a position from cartesian to polar coordinates
     *
     *  \param position position in cartesian coordinates to convert to
     *         cylindrical polar coordinates. Input should be an allocated block
     *         of at least 3 doubles containing values (x, y, z). This is over
     *         written with 3 double containing values (r, y, phi).
     */
    static void convertToPolar(double* position);

    /** Convert a 3 vector from cartesian to polar coordinate system
     *
     *  \param position position in polar coordinates at which the value is
     *         valid
     *  \param value pointer to an allocated block of at least 3 doubles. The
     *         function will apply a rotation to the existing data to render it
     *         from cartesian coordinates (a_x, a_y, a_z) to polar coordinates
     *         (a_r, a_y, a_phi) appropriate for the specified position.
     */
    static void convertToPolar(const double* position_polar, double* value);

    /** Convert a position from polar coordinates to cartesian
     *
     *  \param position position in polar coordinates to convert to
     *         cartesian polar coordinates. Input should be an allocated block
     *         of at least 3 doubles containing values (r, y, phi). This is over
     *         written with 3 double containing values (x, y, z).
     */
    static void convertToCartesian(double* position);

    /** Convert a 3 vector from polar to cartesian coordinate system
     *
     *  \param position position in polar coordinates at which the value is
     *         valid
     *  \param value pointer to an allocated block of at least 3 doubles. The
     *         function will apply a rotation to the existing data to render it
     *         from polar coordinates (a_r, a_y, a_phi) to polar coordinates
     *         (a_x, a_y, a_z) appropriate for the specified position.
     */
    static void convertToCartesian(const double* position_polar, double* value);

    /** Get the minimum bounding box in polar coordinates
     *
     *  \returns bounding box minimum as a 3-vector like (r_min, y_min, phi_min) 
     */
    virtual std::vector<double> getPolarBoundingBoxMin() const;

    /** Get the maximum bounding box in polar coordinates
     *
     *  \returns bounding box maximum as a 3-vector like (r_max, y_max, phi_max) 
     */
    virtual std::vector<double> getPolarBoundingBoxMax() const;

    /** Fill inputs with the bounding box in Polar coordinates
     *
     *  \param zBegin lower bound on field length in phi direction (units of
     *         distance)
     *  \param zEnd upper bound on field length in phi direction (units of
     *         distance)
     *  \param rBegin lower bound on field length in radial direction
     *  \param rEnd lower bound on field length in radial direction
     */
    void getFieldDimensions(double &zBegin, double &zEnd,
                            double &rBegin, double &rEnd) const;

    /** Fill inputs with the bounding box in Cartesian coordinates
     *  
     *  \param xIni lower bound on field x-position (horizontal)
     *  \param xFinal upper bound on field x-position (horizontal)
     *  \param yIni lower bound on field y-position (vertical)
     *  \param yFinal upper bound on field y-position (vertical)
     *  \param zIni lower bound on field z-position (longitudinal)
     *  \param zFinal upper bound on field z-position (longitudinal)
     */
    void getFieldDimensions(double &xIni, double &xFinal,
                            double &yIni, double &yFinal,
                            double &zIni, double &zFinal) const;

    /** Return true if polar vector R_p is within polar bounding box
     *
     *  \param R_p polar coordinates (r, y, phi)
     */
    inline bool isInBoundingBox(const double R_p[]) const;

  protected:
    /** Set the bounding boxes from polar coordinates
     *
     *  Sets the bounding boxes, both polar and cartesian, based on a set of
     *  minimum and maximum polar coordinates. Note that MinPhi->MaxPhi are
     *  allowed in the domain -2pi to 2pi. If the difference is >= 2pi, then
     *  this will describe the full ring.
     *
     *  \param bbMinR minimum radius - must be positive
     *  \param bbMinY minimum y value
     *  \param bbMinPhi minimum phi value - must be greater than -2*pi
     *  \param bbMaxR maximum radius - must be greater than minimum radius
     *  \param bbMaxY maximum y value - must be greater than minimum y value
     *  \param bbMaxPhi maximum phi value - must be greater than minimum phi and
     *                  less than 2*pi
     */
    void setPolarBoundingBox(double bbMinR, double bbMinY, double bbMinPhi,
                             double bbMaxR, double bbMaxY, double bbMaxPhi,
                             double bbTolR, double bbTolY, double bbTolPhi);

    /** bounding box minimum as a 3-vector like (x_min, y_min, z_min) */
    std::vector<double> bbMin_m;
    /** bounding box maximum as a 3-vector like (x_max, y_max, z_max) */
    std::vector<double> bbMax_m;
    /** bounding box minimum as a 3-vector like (r_min, y_min, phi_min) */
    std::vector<double> polarBBMin_m;
    /** bounding box maximum as a 3-vector like (r_max, y_max, phi_max) */
    std::vector<double> polarBBMax_m;
    /** Keep the filename */
    std::string Filename_m;
   
  private:
    std::vector< std::vector<double> > getCorners
               (double bbMinR, double bbMinPhi, double bbMaxR, double bbMaxPhi);
};

bool SectorField::isInBoundingBox(const double R_p[]) const {
    bool inside = R_p[0] >= polarBBMin_m[0] &&
                  R_p[0] <= polarBBMax_m[0] &&
                  R_p[1] >= polarBBMin_m[1] &&
                  R_p[1] <= polarBBMax_m[1] &&
                  R_p[2] >= polarBBMin_m[2] &&
                  R_p[2] <= polarBBMax_m[2];
    // std::cerr << "BB " << inside << "\n" << R_p[0] << " " << R_p[1] << " " << R_p[2] << "\n"
    //          << polarBBMin_m[0] << " " << polarBBMin_m[1] << " " << polarBBMin_m[2] << "\n"
    //          << polarBBMax_m[0] << " " << polarBBMax_m[1] << " " << polarBBMax_m[2] << std::endl;
    return inside;
}

#endif

