/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
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


#ifndef CLASSIC_MULTIPOLETBASE_H
#define CLASSIC_MULTIPOLETBASE_H

/** ---------------------------------------------------------------------
  *
  * MultipoleTBase is a base class for a straight or curved combined \n
  * function magnet (up to arbitrary multipole component) with fringe fields.\n
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
  *
  * ---------------------------------------------------------------------
  *
  * The field is obtained from the scalar potential \n
  *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!} +
  *     ...  @f] \n
  *     (x,z,s) -> Frenet-Serret local coordinates along the magnet \n
  *     z -> vertical component \n
  *     assume mid-plane symmetry \n 
  *     set field on mid-plane -> @f$ B_z = f_0(x,s) = T(x) \cdot S(s) @f$ \n
  *     T(x) -> transverse profile; this is a polynomial describing
  *             the field expansion on the mid-plane inside the magnet
  *             (not in the fringe field);
  *             1st term is the dipole strength, 2nd term is the 
  *             quadrupole gradient * x, etc. \n
  *          -> when setting the magnet, one gives the multipole
  *             coefficients of this polynomial (i.e. dipole strength,  
  *             quadrupole gradient, etc.) \n
  * \n
  * ------------- example ----------------------------------------------- \n
  *     Setting a combined function magnet with dipole, quadrupole and 
  *     sextupole components: \n
  *     @f$ T(x) = B_0 + B_1 \cdot x + B_2 \cdot x^2 @f$\n
  *     user gives @f$ B_0, B_1, B_2 @f$ \n
  * ------------- example end ------------------------------------------- \n
  * \n
  *     S(s) -> fringe field \n
  *     recursion -> @f$ f_n (x,s) = (-1)^n \cdot \sum_{i=0}^{n} C_n^i 
  *     \cdot T^{(2i)} \cdot S^{(2n-2i)} @f$ \n
  *     for curved magnets the above recursion is more complicated \n
  *     @f$ C_n^i @f$ -> binomial coeff; 
  *     @f$ T^{(n)} @f$ -> n-th derivative
  *
  * ---------------------------------------------------------------------
  */

#include "Algorithms/PartBunch.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "Fields/BMultipoleField.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include <vector>

class MultipoleTBase: public Component {
public:
    /** Default constructor */
    MultipoleTBase();
    /** Constructor
     *  \param name -> User-defined name
     */
    explicit MultipoleTBase(const std::string &name);
    /** Copy constructor */
    MultipoleTBase(const MultipoleTBase &right);
    /** Destructor */
    ~MultipoleTBase();
    /** Return a dummy field value */
    EMField &getField();
    /** Return a dummy field value */
    const EMField &getField() const;
    /** Calculate the field at some arbitrary position \n
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param R -> Position in the lab coordinate system of the multipole
     *  \param P -> Not used
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B);
    /** Calculate the field at the position of the ith particle
     *  \param i -> Index of the particle event; field is calculated at this
     *  position
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);
    /** Initialise the MultipoleT
     *  \param bunch -> Bunch the global bunch object
     *  \param startField -> Not used
     *  \param endField -> Not used
     */
    void initialise(PartBunchBase<double, 3>*,
                    double &startField,
                    double &endField);
    /** Finalise the MultipoleT - sets bunch to NULL */
    void finalise();
    /** Return true if dipole component not zero */
    bool bends() const;
    /** Get the dipole constant B_0 */
    double getDipoleConstant() const;
    /** Set the dipole constant B_0 */
    void setDipoleConstant(const double &B0);
    /** Get the number of terms used in calculation of field components */
    std::size_t getMaxOrder() const;
    /** Set the number of terms used in calculation of field components \n
     *  Maximum power of z in Bz is 2 * maxOrder_m
     *  \param maxOrder -> Number of terms in expansion in z
     */
    virtual void setMaxOrder(const std::size_t &maxOrder);
    /** Get the maximum order in the given transverse profile */
    std::size_t getTransMaxOrder() const;
    /** Set the maximum order in the given transverse profile
     *  \param transMaxOrder -> Highest power of x in field expansion 
     */
    void setTransMaxOrder(const std::size_t &transMaxOrder);
    /** Set transverse profile T(x)
      * T(x) = B_0 + B1 x + B2 x^2 + B3 x^3 + ...
      * \param n -> Order of the term (d^n/dx^n) to be set
      * \param Bn -> Value of transverse profile coefficient
      */
    void setTransProfile(const std::size_t &n, const double &Bn);
    /** Get transverse profile
     *  \param n -> Power of x
     */
    double getTransProfile(const std::size_t &n) const;
    /** Get all terms of transverse profile */
    std::vector<double> getTransProfile() const;
    /** Set fringe field model \n
     *  Tanh model used here \n
     *  @f[ 1/2 * \left [tanh \left( \frac{s + s_0}{\lambda_{left}} \right) 
     *  - tanh \left( \frac{s - s_0}{\lambda_{right}} \right) \right] @f] 
     *  \param s0 -> Centre field length and
     *  \lambda_{left} -> Left end field length
     *  \lambda_{right} -> Right end field length
     */
    bool setFringeField(const double &s0, 
                        const double &lambda_left, 
                        const double &lambda_right);
    /** Return vector of 2 doubles
      * [left fringe length, right fringelength]  
      */
    std::vector<double> getFringeLength() const;
    /** Set the entrance angle
     *  \param entranceAngle -> Entrance angle
     */
    void setEntranceAngle(const double &entranceAngle);
    /** Set the bending angle of the magnet */
    virtual void setBendAngle(const double &angle);
    /** Get the bending angle of the magnet */
    virtual double getBendAngle() const;
    /** Get the entrance angle */
    double getEntranceAngle() const;
    /** Set the length of the magnet
      * If straight-> Actual length
      * If curved -> Arc length
      */
    void setLength(const double &length);
    /** Get the length of the magnet */
    double getLength() const;
    /** Set the aperture dimensions \n
      * This element only supports a rectangular aperture
      * \param vertAp -> Vertical aperture length
      * \param horizAp -> Horisontal aperture length
      */
    void setAperture(const double &vertAp, const double &horizAp);
    /** Get the aperture dimensions
      * Returns a vector of 2 doubles
      */
    std::vector<double> getAperture() const;
    /** Set the angle of rotation of the magnet around its axis \n
     *  To make skew components
     *  \param rot -> Angle of rotation
     */
    void setRotation(const double &rot);
    /** Get the angle of rotation of the magnet around its axis */
    double getRotation() const;
    /** Get distance between centre of magnet and entrance */
    double getBoundingBoxLength() const;
    /** Set distance between centre of magnet and enctrance
     *  \param boundingBoxLength -> Distance between centre and entrance
     */
    void setBoundingBoxLength(const double &boundingBoxLength);
    /** Not implemented */
    virtual void getDimensions(double &zBegin, double &zEnd) const;
protected:
    /** Returns the value of the fringe field n-th derivative at s
     *  \param n -> nth derivative
     *  \param s -> Coordinate s
     */
    double getFringeDeriv(const std::size_t &n, const double &s);
    /** Returns the value of the transverse field n-th derivative at x \n
     *  Transverse field is a polynomial in x, differentiation follows
     *  usual polynomial rules of differentiation
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     */
    double getTransDeriv(const std::size_t &n, const double &x);
private:
    //MultipoleTBase operator=(const MultipoleTBase &rhs);
    // End fields
    endfieldmodel::Tanh fringeField_l; // Left
    endfieldmodel::Tanh fringeField_r; // Right
    /** Number of terms in z expansion used in calculating field components */
    std::size_t maxOrder_m;
    /** Highest power in given mid-plane field */
    std::size_t transMaxOrder_m = 0;
    /** List of transverse profile coefficients */
    std::vector<double> transProfile_m;
    /** Rotate frame for skew elements \n
     *  Consecutive rotations:
     *  1st -> about central axis
     *  2nd -> azimuthal rotation
     *  \param R -> Vector to be rotated
     */
    Vector_t rotateFrame(const Vector_t &R);
    /** Inverse of the 1st rotation in rotateFrame() method \n
     *  Used to rotate B field back to global coordinate system
     */
    Vector_t rotateFrameInverse(Vector_t &B);
    /** Transform to Frenet-Serret coordinates for sector magnets */
    virtual void transformCoords(Vector_t &R) = 0;
    /** Transform B-field from Frenet-Serret coordinates to lab coordinates */
    virtual void transformBField(Vector_t &B, const Vector_t &R) = 0;
    /** Magnet parameters */
    double length_m;
    double entranceAngle_m;
    double rotation_m;
    /** Distance between centre of magnet and entrance */
    double boundingBoxLength_m;
    /** Returns the radial component of the field \n
     *  Returns zero far outside fringe field
     *  @f$ Bx = sum_n z^(2n+1) / (2n+1)! * \partial_x f_n @f$
     */
    virtual double getBx (const Vector_t &R);
    /** Returns the vertical field component \n
     *  Returns zero far outside fringe field
     *  @f$ Bz = sum_n  f_n * z^(2n) / (2n)! @f$
     */
    double getBz (const Vector_t &R);
    /** Returns the component of the field along the central axis \n
     *  Returns zero far outside fringe field
      * @f$ Bs = sum_n z^(2n+1) / (2n+1)! \partial_s f_n / h_s @f$
      */
    virtual double getBs (const Vector_t &R);
    /** Assume rectangular aperture with these dimensions */
    double verticalApert_m;
    double horizontalApert_m;
    /** Not implemented */
    BMultipoleField dummy;
    /** Tests if inside the magnet
     *  \param R -> Coordinate vector
     */
    bool insideAperture(const Vector_t &R);
    /** Radius of curvature
     *  \param s -> Coordinate s
     */
    virtual double getRadius(const double &s) = 0;
    /** Returns the scale factor @f$ h_s = 1 + x / \rho(s) @f$
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getScaleFactor(const double &x, const double &s) = 0;
    /** Calculate partial derivative of fn wrt x using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivX(const std::size_t &n,
                       const double &x, 
                       const double &s);
    /** Calculate partial derivative of fn wrt s using a 5-point
     *  finite difference formula
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivS(const std::size_t &n, 
                       const double &x,
                       const double &s);
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getFn(const std::size_t &n,
                         const double &x,
                         const double &s) = 0;
};

inline
    void MultipoleTBase::finalise() {
        RefPartBunch_m = NULL;
}
inline
    bool MultipoleTBase::apply(const size_t &i, const double &t,
                           Vector_t &E, Vector_t &B) {
        return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}
inline
void MultipoleTBase::setBendAngle(const double &/*angle*/) {}
inline
    double MultipoleTBase::getBendAngle() const {
        return 0.0;
}
inline
    void MultipoleTBase::setEntranceAngle(const double &entranceAngle) {
        entranceAngle_m = entranceAngle;
}
inline
    bool MultipoleTBase::insideAperture(const Vector_t &R) {
        return (fabs(R[1]) <= (verticalApert_m / 2.0) &&
                fabs(R[0]) <= (horizontalApert_m / 2.0));
}
inline
    double MultipoleTBase::getEntranceAngle() const {
        return entranceAngle_m;
}
inline
    double MultipoleTBase::getTransProfile(const std::size_t &n) const {
        return transProfile_m[n];
}
inline
    std::vector<double> MultipoleTBase::getTransProfile() const {
        return transProfile_m;
}
inline
    double MultipoleTBase::getDipoleConstant() const {
         return transProfile_m[0];
} 
inline
    void MultipoleTBase::setMaxOrder(const std::size_t &maxOrder) {
        maxOrder_m = maxOrder;
}
inline
    std::size_t MultipoleTBase::getMaxOrder() const {
         return maxOrder_m;
}
inline
    std::size_t MultipoleTBase::getTransMaxOrder() const {
        return transMaxOrder_m;
}
inline
    void MultipoleTBase::setTransMaxOrder(const std::size_t &transMaxOrder) {
        transMaxOrder_m = transMaxOrder;
        transProfile_m.resize(transMaxOrder + 1, 0.);
}
inline
    double MultipoleTBase::getRotation() const {
         return rotation_m;
}
inline
    void MultipoleTBase::setRotation(const double &rot) {
         rotation_m = rot;
}
inline
    void MultipoleTBase::setLength(const double &length) {
        length_m = std::abs(length);
}
inline
    double MultipoleTBase::getLength() const {
        return length_m;
}
inline
    double MultipoleTBase::getBoundingBoxLength() const {
        return boundingBoxLength_m;
}
inline
    void MultipoleTBase::setBoundingBoxLength(const double &boundingBoxLength) {
        boundingBoxLength_m = boundingBoxLength;
}
inline
    void MultipoleTBase::setTransProfile(const std::size_t &n,
                                         const double &dTn) {
        if (n > transMaxOrder_m) {
            transMaxOrder_m = n;
            transProfile_m.resize(n + 1, 0.0);
        }
        transProfile_m[n] = dTn;
}
inline
    void MultipoleTBase::setDipoleConstant(const double &B0) {
        if (transMaxOrder_m < 1) {
            transProfile_m.resize(1, 0.);
        }
        transProfile_m[0] = B0;
}
inline
    void MultipoleTBase::setAperture(const double &vertAp,
                                     const double &horizAp) {
        verticalApert_m = vertAp;
        horizontalApert_m = horizAp;
}
inline
    std::vector<double> MultipoleTBase::getAperture() const {
        std::vector<double> temp(2, 0.0);
        temp[0] = verticalApert_m;
        temp[1] = horizontalApert_m;
        return temp;
}
inline
    std::vector<double> MultipoleTBase::getFringeLength() const {
        std::vector<double> temp(2, 0.0);
        temp[0] = fringeField_l.getLambda();
        temp[1] = fringeField_r.getLambda();
        return temp;
}
inline
void MultipoleTBase::initialise(PartBunchBase<double, 3>* /*bunch*/,
                                double &/*startField*/,
                                double &/*endField*/) {
}
inline
    bool MultipoleTBase::bends() const {
        return transProfile_m[0] != 0;
}
inline
    EMField &MultipoleTBase::getField() {
        return dummy;
}
inline
    const EMField &MultipoleTBase::getField() const {
        return dummy;
}
inline
    void MultipoleTBase::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {
}

#endif
