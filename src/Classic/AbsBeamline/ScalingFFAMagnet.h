//
// Class ScalingFFAMagnet
//   Defines the abstract interface for a sector FFA magnet
//   with radially scaling fringe fields.
//
// Copyright (c) 2017 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef ABSBEAMLINE_ScalingFFAMagnet_H
#define ABSBEAMLINE_ScalingFFAMagnet_H

#include "AbsBeamline/Component.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Algorithms/PartBunch.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"

#include <cmath>

/** Sector bending magnet with an FFA-style field index and spiral end shape
 * 
 *  Note about placement and end field; in order to get a user-defined end field
 *  shape, we do a lookup from user defined end fields, which is pulled from
 *  the EndFieldModel object and done by setupEndField(). Because end fields tie
 *  into the geometry, this has to be done before placements; but it also needs
 *  to be done after parsing has finished. setupEndField() is called by e.g.
 *  ParallelCyclotronTracker just before the object is handed to OpalRing for
 *  placement.
 */

class ScalingFFAMagnet: public Component {

public:
    /** Construct a new ScalingFFAMagnet
     *
     *  \param name User-defined name of the ScalingFFAMagnet
     */
    explicit ScalingFFAMagnet(const std::string& name);

    /** Destructor - deletes map */
    ~ScalingFFAMagnet();

    /** Inheritable copy constructor */
    ScalingFFAMagnet* clone() const override;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event override; field is calculated at this
     *         position
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map
     */
    inline bool apply(const size_t& i, const double& t,
                      Vector_t& E, Vector_t& B) override;

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the bend
     *  \param P not used
     *  \param t not used
     *  \param E not used
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    inline bool apply(const Vector_t& R, const Vector_t& P, const double& t,
               Vector_t& E, Vector_t& B) override;

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValue(const Vector_t& R, Vector_t& B) const;

    /** Calculate the field at some arbitrary position in cylindrical coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cylindrical polar coordinates defined like (r, z, phi)
     *  \param B calculated magnetic field defined like (Br, By, Bphi)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValueCylindrical(const Vector_t& R, Vector_t& B) const;

     /** Initialise the ScalingFFAMagnet
      *
      *  \param bunch the global bunch object
      *  \param startField not used
      *  \param endField not used
      */
      void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) override;

     /** Initialise the ScalingFFAMagnet
      *
      *  Sets up the field expansion and the geometry; call after changing any
      *  field parameters
      */
    void initialise();

     /** Finalise the ScalingFFAMagnet - sets bunch to nullptr */
    void finalise() override;

    /** Return true - ScalingFFAMagnet always bends the reference particle */
    inline bool bends() const override;

    /** Not implemented */
    void getDimensions(double& /*zBegin*/, double& /*zEnd*/) const override {}

    /** Return the cell geometry */
    BGeometryBase& getGeometry() override;

    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const override;

    /** Return a dummy (0.) field value (what is this for?) */
    EMField& getField() override;

    /** Return a dummy (0.) field value (what is this for?) */
    const EMField& getField() const override;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;

    /** Get tan delta - delta is the spiral angle */
    double getTanDelta() const {return tanDelta_m;}

    /** Set tan delta - delta is the spiral angle */
    void setTanDelta(double tanDelta) {tanDelta_m = tanDelta;}

    /** Get the field index k */
    double getFieldIndex() const {return k_m;}

    /** Set the field index k */
    void setFieldIndex(double k) {k_m = k;}

    /** Get the dipole constant B_0 */
    double getDipoleConstant() const {return Bz_m;}

    /** Set the dipole constant B_0 */
    void setDipoleConstant(double Bz) {Bz_m = Bz;}

    /** Get the radius constant R_0 */
    double getR0() const {return r0_m;}

    /** Set the radius constant R_0 */
    void setR0(double r0) {r0_m = r0; r0Sign_m = r0_m / std::abs(r0_m);}

    /** Get the centre of the sector */
    Vector_t getCentre() const {return centre_m;}

    /** Set the centre of the sector */
    void setCentre(Vector_t centre) {centre_m = centre;}

    /** Get the fringe field
     *
     *  Returns the fringe field model; ScalingFFAMagnet retains ownership of the
     *  returned memory.
     */
    endfieldmodel::EndFieldModel* getEndField() const {return endField_m;}

    /** Set the fringe field
      * 
      * - endField: the new fringe field; ScalingFFAMagnet takes ownership of the
      *   memory associated with endField.
      */
    void setEndField(endfieldmodel::EndFieldModel* endField);

    /** Get the maximum power of y modelled in the off-midplane expansion; 
     */
    size_t getMaxOrder() const {return maxOrder_m;}

    /** Set the maximum power of y modelled in the off-midplane expansion;
     */
    void setMaxOrder(size_t maxOrder) {maxOrder_m = maxOrder;}

    /** Get the offset of the magnet centre from the start 
     */
    double getPhiStart() const {return phiStart_m;}

    /** Set the offset of the magnet centre from the start 
     */
    void setPhiStart(double phiStart) {phiStart_m = phiStart;}

    /** Get the offset of the magnet end from the start 
     */
    double getPhiEnd() const {return phiEnd_m;}

    /** Set the offset of the magnet end from the start 
     */
    void setPhiEnd(double phiEnd) {phiEnd_m = phiEnd;}

    /** Get the maximum radius
     */
    double getRMin() const {return rMin_m;}

    /** Set the maximum radius 
     */
    void setRMin(double rMin) {rMin_m = rMin;}

    /** Get the maximum radius
     */
    double getRMax() const {return rMax_m;}

    /** Set the maximum radius 
     */
    void setRMax(double rMax) {rMax_m = rMax;}

    /** Get the maximum azimuthal displacement from \psi=0
     */
    double getAzimuthalExtent() const {return azimuthalExtent_m;}

    /** Set the maximum azimuthal displacement from \psi=0
     */
    void setAzimuthalExtent(double azimuthalExtent) {azimuthalExtent_m = azimuthalExtent;}

    /** Get the maximum vertical displacement from the midplane
     */
    double getVerticalExtent() const {return verticalExtent_m;}

    /** Set the maximum vertical displacement from the midplane
     */
    void setVerticalExtent(double verticalExtent) {verticalExtent_m = verticalExtent;}

    /** Return the calculated df coefficients */
    std::vector<std::vector<double> > getDfCoefficients() {return dfCoefficients_m;}

    /** setupEndField does some end field and geometry set-up
     * 
     *  This is normally called just before the magnet is placed; can only set
     *  up the end field after everything has been parsed from input (otherwise
     *  OPAL may not know about an end field model).
     *
     *  sets PhiStart, PhiEnd, AzimuthalExtent and the end field model itself.
     */
    void setupEndField();

    /** Set the end field name.
     * 
     *  Called during parsing of the input file; OPAL looks for the endFieldName
     *  when setupEndField() is called.
     */
    void setEndFieldName(std::string name) {endFieldName_m = name;}

    /** Return the end field name. */
    std::string getEndFieldName() const {return endFieldName_m;}


private:
    /** Calculate the df coefficients, ready for field generation
     *
     *  Must be called following any update to the the field parameters, in
     *  order for correct field to be calculated.
     */
    void calculateDfCoefficients();

    /** Copy constructor */
    ScalingFFAMagnet(const ScalingFFAMagnet& right);

    ScalingFFAMagnet& operator=(const ScalingFFAMagnet& rhs);
    PlanarArcGeometry planarArcGeometry_m;
    BMultipoleField dummy;

    size_t maxOrder_m = 0; // maximum order used in the calculation
    double tanDelta_m = 0.; // tan(spiral angle)
    double k_m = 0.; // field index
    double Bz_m = 0.; // nominal field on the radius
    double r0_m = 0.; // radius; negative for clockwise geometry, else positive
    double rMin_m = 0.; // minimum radius
    double rMax_m = 0.; // maximum radius
    double phiStart_m = 0.; // offsets this element
    double phiEnd_m = 0.; // used for placement of next element
    double azimuthalExtent_m = 0.; // maximum distance used for field calculation
    double verticalExtent_m = 0.; // maximum allowed distance from the midplane
    double r0Sign_m = 0.0; // -1 for clockwise geometry, else +1
    Vector_t centre_m; // nominal ring centre
    endfieldmodel::EndFieldModel* endField_m = nullptr;
    std::string endFieldName_m = ""; 
    const double fp_tolerance = 1e-18;
    std::vector<std::vector<double> > dfCoefficients_m;
};

bool ScalingFFAMagnet::apply(const Vector_t& R, const Vector_t& /*P*/,
                             const double& /*t*/, Vector_t& /*E*/, Vector_t& B) {
    return getFieldValue(R, B);
}

bool ScalingFFAMagnet::apply(const size_t& i, const double& t,
                             Vector_t& E, Vector_t& B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

#endif
