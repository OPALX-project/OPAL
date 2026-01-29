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

#include "Fields/BMultipoleField.h"
#include "Fields/SectorMagneticFieldMap.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "AbsBeamline/Component.h"

#ifndef ABSBEAMLINE_SBEND3D_H
#define ABSBEAMLINE_SBEND3D_H

/** Sector bending magnet from a 3D field map
 *
 *  SBend3D provides interface between Component type and routines for linear
 *  interpolation from a 3D field map in a sector geometry.
 *
 *  \param planarArcGeometry_m Cell geometry is a PlanarArcGeometry with radius
 *         of curvature equal to the radius at the middle of the field map and
 *         angle given by the field map opening angle.
 *  \param map_m the actual field map class (with associated interpolation and
 *         IO routines).
 *  \param field_m a dummy field - never used, only put in to obey inheritance
 *         requirements
 */

class SBend3D : public Component {
  public:
    /** Construct a new SBend3D
     *
     *  \param name User-defined name of the SBend3D
     */
    explicit SBend3D(const std::string &name);

    /** Copy constructor */
    SBend3D(const SBend3D &right);

    /** Destructor - deletes map */
    ~SBend3D();

    /** Inheritable copy constructor */
    ElementBase* clone() const override;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event; field is calculated at this
     *         position
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map
     */
    bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the bend
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B) override;

     /** Initialise the SBend3D
      *
      *  \param bunch the global bunch object
      *  \param startField not used
      *  \param endField not used
      */
    void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

     /** Finalise the SBend3D - sets bunch to NULL */
    void finalise() override;

    /** Return true - SBend3D always bends the reference particle */
    inline bool bends() const override;

    /** Not implemented */
    void getDimensions(double &/*zBegin*/, double &/*zEnd*/) const override {}

    /** Return the cell geometry */
    BGeometryBase& getGeometry() override;

    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const override;

    /** Return a dummy (0.) field value (what is this for?) */
    EMField &getField() override;

    /** Return a dummy (0.) field value (what is this for?) */
    const EMField &getField() const override;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;

    /** Set the field map file
     *
     *  @param name the file name of the field map
     *
     *  This generates a new field map object with the given name. If name == ""
     *  then sets the field map to NULL.
     */
    void setFieldMapFileName(std::string name);

    /** Get the file name of the field map */
    inline std::string getFieldMapFileName() const;

    /** Get the scale factor */
    double getLengthUnits() const {return lengthUnits_m;}

    /** Set the scale factor */
    void setLengthUnits(double lengthUnits) {lengthUnits_m = lengthUnits;}

    /** Get the scale factor */
    double getFieldUnits() const {return fieldUnits_m;}

    /** Set the scale factor */
    void setFieldUnits(double fieldUnits) {fieldUnits_m = fieldUnits;}

    /** Get the polynomial order */
    int getPolynomialOrder() const {return polyOrder_m;}

    /** Set the polynomial order */
    void setPolynomialOrder(int polyOrder) {polyOrder_m = polyOrder;}

    /** Get the smoothing factor */
    int getSmoothingOrder() const {return smoothOrder_m;}

    /** Set the smoothing factor */
    void setSmoothingOrder(int polyOrder) {smoothOrder_m = polyOrder;}

    /** Get the sector magnetic field map */
    SectorMagneticFieldMap* getSectorMagneticFieldMap() const {return map_m;}

  private:
    SectorMagneticFieldMap* map_m;
    // Geometry
    PlanarArcGeometry planarArcGeometry_m;
    // Units used for field maps
    double fieldUnits_m;
    double lengthUnits_m;
    // Polynomial interpolation
    int polyOrder_m;
    double smoothOrder_m;
    // Not implemented (I think this is something to do with error fields?)
    BMultipoleField dummy;
};

std::string SBend3D::getFieldMapFileName() const {
    return map_m->getFieldMapFileName();
}

#endif
