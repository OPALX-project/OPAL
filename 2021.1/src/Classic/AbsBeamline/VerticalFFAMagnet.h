//
// Header file for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "Fields/BMultipoleField.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "AbsBeamline/Component.h"
#include "Algorithms/PartBunch.h"

#ifndef ABSBEAMLINE_VerticalFFAMagnet_H
#define ABSBEAMLINE_VerticalFFAMagnet_H

namespace endfieldmodel {
    class EndFieldModel;
}

/** Bending magnet with an exponential dependence on field in the vertical plane
 *
 *  VerticalFFAMagnet makes a rectangular bending magnet with a dipole field
 *  that has a dependence like B0 exp(mz)
 */

class VerticalFFAMagnet : public Component {
  public:
    /** Construct a new VerticalFFAMagnet
     *
     *  \param name User-defined name of the VerticalFFAMagnet
     */
    explicit VerticalFFAMagnet(const std::string &name);

    /** Destructor - deletes the field */
    ~VerticalFFAMagnet();

    /** Inheritable copy constructor */
    ElementBase* clone() const;

    /** Calculate the field at the position of the ith particle
     *
     *  \param i index of the particle event; field is calculated at this
     *         position
     *  \param t time at which the field is to be calculated
     *  \param E calculated electric field - always 0 (no E-field)
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map
     */
    inline bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    /** Calculate the field at some arbitrary position
     *
     *  \param R position in the local coordinate system of the magnet
     *  \param P not used
     *  \param t not used
     *  \param E not used
     *  \param B calculated magnetic field
     *  \returns true if particle is outside the field map, else false
     */
    inline bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B);

    /** Calculate the field at some arbitrary position in cartesian coordinates
     *
     *  \param R position in the local coordinate system of the bend, in
     *           cartesian coordinates defined like (x, y, z)
     *  \param B calculated magnetic field defined like (Bx, By, Bz)
     *  \returns true if particle is outside the field map, else false
     */
    bool getFieldValue(const Vector_t &R, Vector_t &B) const;

     /** Initialise the VerticalFFAMagnet
      *
      *  \param bunch the global bunch object (but not used)
      *  \param startField not used
      *  \param endField not used
      */
     void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

     /** Initialise the VerticalFFAMagnet
      *
      *  Sets up the field expansion and the geometry; call after changing any
      *  field parameters
      */
    void initialise();

     /** Finalise the VerticalFFAMagnet - sets bunch to NULL */
    void finalise();

    /** Return false - VerticalFFAMagnet is a straight magnet 
     *
     *  Nb: the VerticalFFAMagnet geometry is straight even though trajectories
     *      are not
     */
    inline bool bends() const {return false;}

    /** Not implemented */
    void getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {}

    /** Return the cell geometry */
    BGeometryBase& getGeometry();

    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const;

    /** Return a dummy (0.) field value (what is this for?) */
    EMField &getField();

    /** Return a dummy (0.) field value (what is this for?) */
    const EMField &getField() const;

    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const;

    /** Get the fringe field
     *
     *  Returns the fringe field model; VerticalFFAMagnet retains ownership of
     *  the returned memory.
     */
    endfieldmodel::EndFieldModel* getEndField() const {return endField_m.get();}

    /** Set the fringe field
      *
      * - endField: the new fringe field; VerticalFFAMagnet takes ownership of
      *   the memory associated with endField.
      */
    void setEndField(endfieldmodel::EndFieldModel* endField);

    /** Get the maximum power of x used in the off-midplane expansion; 
     */
    size_t getMaxOrder() const {return maxOrder_m;}

    /** Set the maximum power of x used in the off-midplane expansion;
     */
    void setMaxOrder(size_t maxOrder) {maxOrder_m = maxOrder;}

    /** Get the centre field at z=0 */
    double getB0() const {return Bz_m/Tesla;}

    /** Set the centre field at z=0 */
    void setB0(double Bz) {Bz_m = Bz*Tesla;}

    /** Get the field index */
    double getFieldIndex() const {return k_m*mm;} // units are [m^{-1}]

    /** Set the field index */
    void setFieldIndex(double index) {k_m = index/mm;}

    /** Get the maximum extent below z = 0 */
    double getNegativeVerticalExtent() const {return zNegExtent_m/mm;}

    /** Set the maximum extent below z = 0 */
    inline void setNegativeVerticalExtent(double negativeExtent);

    /** Get the maximum extent above z = 0 */
    double getPositiveVerticalExtent() const {return zPosExtent_m/mm;}

    /** set the maximum extent above z = 0 */
    inline void setPositiveVerticalExtent(double positiveExtent);

    /** Get the length of the bounding box (centred on magnet centre) */
    double getBBLength() const {return bbLength_m/mm;}

    /** Set the length of the bounding box (centred on magnet centre) */
    void setBBLength(double bbLength) {bbLength_m = bbLength*mm;}

    /** Get the full width of the bounding box (centred on magnet centre) */
    double getWidth() const {return halfWidth_m/mm*2.;}

    /** Set the full width of the bounding box (centred on magnet centre) */
    void setWidth(double width) {halfWidth_m = width/2*mm;}

    /** Get the coefficients used for the field expansion
     *  
     *  B_y is given by 
     *     sum_n B_0 exp(ky) f_n x^n
     *  where 
     *     f_n = sum_k c_{nk} partial_k f_0
     * 
     *  Returns a vector of vectors, like c[n][k]. The expansion for the other
     *  field elements can be related back to c[n][k] (see elsewhere for details).
     */
    inline std::vector<std::vector<double> > getDfCoefficients() const;
private:
    void calculateDfCoefficients();

    /** Copy constructor */
    VerticalFFAMagnet(const VerticalFFAMagnet &right);

    VerticalFFAMagnet& operator=(const VerticalFFAMagnet& rhs);
    StraightGeometry straightGeometry_m;
    BMultipoleField dummy;

    size_t maxOrder_m = 0;
    double k_m = 0.;
    double Bz_m = 0.;
    double zNegExtent_m = 0.; // extent downwards from the midplane
    double zPosExtent_m = 0.; // extent upwards from the midplane
    double halfWidth_m = 0.; // extent in either +x or -x
    double bbLength_m = 0.;
    std::unique_ptr<endfieldmodel::EndFieldModel> endField_m;
    std::vector<std::vector<double> > dfCoefficients_m;

    const double mm=1000.;
    const double Tesla=10.;
};

void VerticalFFAMagnet::setNegativeVerticalExtent(double negativeExtent) {
    zNegExtent_m = negativeExtent*mm;
}

void VerticalFFAMagnet::setPositiveVerticalExtent(double positiveExtent) {
    zPosExtent_m = positiveExtent*mm;
}

bool VerticalFFAMagnet::apply(const size_t &i, const double &t,
                              Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool VerticalFFAMagnet::apply(const Vector_t &R, const Vector_t &/*P*/,
                              const double &, Vector_t &/*E*/, Vector_t &B) {
    return getFieldValue(R, B);
}

std::vector<std::vector<double> > VerticalFFAMagnet::getDfCoefficients() const {
    return dfCoefficients_m;
}

#endif

