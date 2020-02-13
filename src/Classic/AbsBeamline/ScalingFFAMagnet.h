/*
 *  Copyright (c) 2017, Chris Rogers
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
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/Component.h"

#ifndef ABSBEAMLINE_ScalingFFAMagnet_H
#define ABSBEAMLINE_ScalingFFAMagnet_H

/** Sector bending magnet with an FFA-style field index and spiral end shape
 */

class ScalingFFAMagnet : public Component {
  public:
    /** Construct a new ScalingFFAMagnet
     *
     *  \param name User-defined name of the ScalingFFAMagnet
     */
    explicit ScalingFFAMagnet(const std::string &name);

    /** Destructor - deletes map */
    ~ScalingFFAMagnet();

    /** Inheritable copy constructor */
    ElementBase* clone() const override;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event override; field is calculated at this
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
     *  \param P not used
     *  \param t not used
     *  \param E not used
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B) override;

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValue(const Vector_t &R, Vector_t &B) const;

    /** Calculate the field at some arbitrary position in cylindrical coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cylindrical polar coordinates defined like (r, y, phi)
     *  \param B calculated magnetic field defined like (Br, By, Bphi)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValueCylindrical(const Vector_t &R, Vector_t &B) const;

     /** Initialise the ScalingFFAMagnet
      *
      *  \param bunch the global bunch object
      *  \param startField not used
      *  \param endField not used
      */
      void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

     /** Initialise the ScalingFFAMagnet
      *
      *  Sets up the field expansion and the geometry; call after changing any
      *  field parameters
      */
    void initialise();

     /** Finalise the ScalingFFAMagnet - sets bunch to NULL */
    void finalise() override;

    /** Return true - ScalingFFAMagnet always bends the reference particle */
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
    void setR0(double r0) {r0_m = r0;}

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

  private:
    /** Calculate the df coefficients, ready for field generation
     *
     *  Must be called following any update to the the field parameters, in
     *  order for correct field to be calculated.
     */
    void calculateDfCoefficients();

    /** Copy constructor */
    ScalingFFAMagnet(const ScalingFFAMagnet &right);

    ScalingFFAMagnet& operator=(const ScalingFFAMagnet& rhs);
    PlanarArcGeometry planarArcGeometry_m;
    BMultipoleField dummy;

    size_t maxOrder_m = 0;
    double tanDelta_m = 0.;
    double k_m = 0.;
    double Bz_m = 0.;
    double r0_m = 0.;
    double rMin_m = 0.; // minimum radius
    double rMax_m = 0.; // maximum radius
    double phiStart_m = 0.; // offsets this element
    double phiEnd_m = 0.; // used for placement of next element
    double azimuthalExtent_m = 0.; // maximum distance used for field calculation
    double verticalExtent_m = 0.; // maximum allowed distance from the midplane
    Vector_t centre_m;
    endfieldmodel::EndFieldModel* endField_m = NULL;
    const double fp_tolerance = 1e-18;
    std::vector<std::vector<double> > dfCoefficients_m;
};

#endif

