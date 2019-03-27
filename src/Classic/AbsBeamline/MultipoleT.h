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


#ifndef CLASSIC_MULTIPOLET_H
#define CLASSIC_MULTIPOLET_H

/** ---------------------------------------------------------------------
  *
  * MultipoleT defines a straight or curved combined function magnet (up 
  * to arbitrary multipole component) with Fringe Fields
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

#include "AbsBeamline/EndFieldModel/Tanh.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/MultipoleTFunctions/RecursionRelation.h"
#include "AbsBeamline/MultipoleTFunctions/RecursionRelationTwo.h"
#include "gsl/gsl_sf.h"
#include <vector>

class MultipoleT: public Component {
public: 
    /** Constructor
     *  \param name -> User-defined name
     */
    explicit MultipoleT(const std::string &name);
    /** Copy constructor */
    MultipoleT(const MultipoleT &right);
    /** Destructor */ 
    ~MultipoleT();
    /** Inheritable copy constructor */
    ElementBase* clone() const override;
    /** Return a dummy field value */
    EMField &getField() override;
    /** Return a dummy field value */
    const EMField &getField() const override;
    /** Not implemented */
    void getDimensions(double &zBegin, double &zEnd) const override;
    /** Calculate the field at some arbitrary position
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param R -> Position in the lab coordinate system of the multipole
     *  \param P -> Not used
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(const Vector_t &R, const Vector_t &P, const double &t,
               Vector_t &E, Vector_t &B) override;
    /** Calculate the field at the position of the ith particle
     *  \param i -> Index of the particle event override; field is calculated at this
     *  position
     *  If particle is outside field map true is returned,
     *  otherwise false is returned
     *  \param t -> Time at which the field is to be calculated
     *  \param E -> Calculated electric field - always 0 (no E-field)
     *  \param B -> Calculated magnetic field
     */
    bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;
    /** Initialise the MultipoleT
     *  \param bunch -> Bunch the global bunch object
     *  \param startField -> Not used
     *  \param endField -> Not used
     */
    void initialise(PartBunchBase<double, 3>*,
                    double &startField,
                    double &endField) override;
    /** Initialises the geometry */
    void initialise();

    /** Finalise the MultipoleT - sets bunch to NULL */
    void finalise() override;
    /** Return true if dipole component not zero */
    bool bends() const override;
    /** Return the cell geometry */
    PlanarArcGeometry& getGeometry() override;
    /** Return the cell geometry */
    const PlanarArcGeometry& getGeometry() const override;
    /** Accept a beamline visitor */
    void accept(BeamlineVisitor& visitor) const override;
    /** Get the dipole constant B_0 */
    double getDipoleConstant() const;
    /** Set the dipole constant B_0 */
    void setDipoleConstant(double B0);
    /** Get the number of terms used in calculation of field components */
    std::size_t getMaxOrder() const;
    /** Set the number of terms used in calculation of field components \n
     *  Maximum power of z in Bz is 2 * maxOrder_m \n
     *  \param maxOrder -> Number of terms in expansion in z
     */
    void setMaxOrder(std::size_t maxOrder);
    /** Get highest power of x in polynomial expansions */
    std::size_t getMaxXOrder() const;
    /** Set the number of terms used in polynomial expansions
     *  \param maxXOrder -> Number of terms in expansion in z
     */
    void setMaxXOrder(std::size_t maxXOrder);
    /** Get the maximum order in the given transverse profile */
    std::size_t getTransMaxOrder() const;
    /** Set the maximum order in the given transverse profile
     *  \param transMaxOrder -> Highest power of x in field expansion 
     */
    void setTransMaxOrder(std::size_t transMaxOrder);
    /** Set transverse profile T(x)
      * T(x) = B_0 + B1 x + B2 x^2 + B3 x^3 + ...
      * \param n -> Order of the term (d^n/dx^n) to be set
      * \param Bn -> Value of transverse profile coefficient
      */
    void setTransProfile(std::size_t n, double Bn);
    /** Get transverse profile
     *  \param n -> Power of x
     */
    double getTransProfile(int n) const;
    /** Get all terms of transverse profile */
    std::vector<double> getTransProfile() const;
    /** Set fringe field model \n
     *  Tanh model used here \n
     *  @f[ 1/2 * \left [tanh \left( \frac{s + s_0}{\lambda_{left}} \right) 
     *  - tanh \left( \frac{s - s_0}{\lambda_{right}} \right) \right] @f] 
     *  \param s0 -> Centre field length
     *  \param \lambda_{left} -> Left end field length
     *  \param \lambda_{right} -> Right end field length
     */
    bool setFringeField(double s0, double lambda_left, double lambda_right);
    /** Return vector of 2 doubles
      * [left fringe length, right fringelength]  
      */
    std::vector<double> getFringeLength() const;
    /** Set the bending angle of the magnet */
    void setBendAngle(double angle);
    /** Get the bending angle of the magnet */
    double getBendAngle() const;
    /** Set the entrance angle
     *  \param entranceAngle -> Entrance angle
     */
    void setEntranceAngle(double entranceAngle);
    /** Get the entrance angle */
    double getEntranceAngle() const;
    /** Get the bending radius \n 
      * Not used, when needed radius is found from length_m / angle_m
      */
    double getBendRadius() const;
    /** Set the length of the magnet \n
      * If straight-> Actual length \n
      * If curved -> Arc length \n
      */
    void setLength(double length);
    /** Get the length of the magnet */
    double getLength() const;
    /** Not used */
    double getChordLength() const;
    /** Set the aperture dimensions
      * This element only supports a rectangular aperture
      * \param vertAp -> Vertical aperture length
      * \param horizAp -> Horisontal aperture length
      */
    void setAperture(double vertAp, double horizAp);
    /** Get the aperture dimensions
      * Returns a vector of 2 doubles
      */
    std::vector<double> getAperture() const;
    /** Set the angle of rotation of the magnet around its axis \n
     *  To make skew components
     *  \param rot -> Angle of rotation
     */
    void setRotation(double rot);
    /** Get the angle of rotation of the magnet around its axis */
    double getRotation() const;
    /** Set variable radius flag to true */
    void setVarRadius();
    /** Get the value of variableRadius_m */
    bool getVarRadius() const;
    /** Get distance between centre of magnet and entrance */
    double getBoundingBoxLength() const;
    /** Set distance between centre of magnet and enctrance
     *  \param boundingBoxLength -> Distance between centre and entrance
     */
    void setBoundingBoxLength(const double &boundingBoxLength);

private:
    MultipoleT operator=(const MultipoleT& rhs);
    // End fields
    endfieldmodel::Tanh fringeField_l; // Left
    endfieldmodel::Tanh fringeField_r; // Right
    /** Field expansion parameters \n
     *  Number of terms in z expansion used in calculating field components
     */
    std::size_t maxOrder_m = 0;
    /** Highest order of polynomial expansions in x */
    std::size_t maxOrderX_m = 0;
    /** Objects for storing differential operator acting on Fn */
    std::vector<polynomial::RecursionRelationTwo> recursion_VarRadius_m;
    std::vector<polynomial::RecursionRelation> recursion_ConstRadius_m;
    /** Highest power in given mid-plane field */
    std::size_t transMaxOrder_m = 0; 
    /** List of transverse profile coefficients */
    std::vector<double> transProfile_m;
    /** Geometry */
    PlanarArcGeometry planarArcGeometry_m;
    /** Rotate frame for skew elements
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
    Vector_t transformCoords(const Vector_t &R);
    /** Magnet parameters */
    double length_m;
    double angle_m;
    double entranceAngle_m;
    double rotation_m;
    /** Variable radius flag */
    bool variableRadius_m;
    /** Distance between centre of magnet and entrance */
    double boundingBoxLength_m;
    /** Get field component methods
     */
    double getBx (const Vector_t &R);
    double getBz (const Vector_t &R);
    double getBs (const Vector_t &R);
    /** Assume rectangular aperture with these dimensions */
    double verticalApert_m;
    double horizApert_m;
    /** Not implemented */
    BMultipoleField dummy;
    /** Returns the value of the fringe field n-th derivative at s
     *  \param n -> nth derivative
     *  \param s -> Coordinate s
     */
    double getFringeDeriv(int n, double s);
    /** Returns the value of the transverse field n-th derivative at x
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     */
    double getTransDeriv(std::size_t n, double x);
    /** Tests if inside the magnet
     *  \param R -> Coordinate vector
     */
    bool insideAperture(const Vector_t &R);
    /** Radius of curvature \n
     *  If radius of curvature is infinite, -1 is returned \n
     *  If radius is constant, then \n
     *  @f$ \rho(s) = length_m / angle_m @f$ \n
     *  If radius is variable, then
     *  @f$ \rho(s) = rho(0) * S(0) / S(s) @f$
     *  where S(s) is the fringe field
     *  \param s -> Coordinate s
     */
    double getRadius(double s);
    /** Returns the scale factor @f$ h_s = 1 + x / \rho(s) @f$
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getScaleFactor(double x, double s);
    /** Calculate partial derivative of fn wrt x using a 5-point
     *  finite difference formula \n
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivX(std::size_t n, double x, double s);
    /** Calculate partial derivative of fn wrt s using a 5-point
     *  finite difference formuln
     *  Error of order stepSize^4
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFnDerivS(std::size_t n, double x, double s);
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    double getFn(std::size_t n, double x, double s);
};

inline
    void MultipoleT::setVarRadius() {
        variableRadius_m = true;
}
inline
    bool MultipoleT::getVarRadius() const {
        return variableRadius_m;
}
inline
    void MultipoleT::setEntranceAngle(double entranceAngle) {
        entranceAngle_m = entranceAngle;
}
inline
    double MultipoleT::getEntranceAngle() const {
        return entranceAngle_m;
}
inline 
    double MultipoleT::getTransProfile(int n) const {
        return transProfile_m[n];
}
inline
    std::vector<double> MultipoleT::getTransProfile() const {
        return transProfile_m;
}
inline 
    double MultipoleT::getDipoleConstant() const {
         return transProfile_m[0];
}
inline
    std::size_t MultipoleT::getMaxOrder() const {
         return maxOrder_m;
}

inline
    std::size_t MultipoleT::getMaxXOrder() const {
        return maxOrderX_m;
}

inline
    void MultipoleT::setMaxXOrder(std::size_t maxOrderX) {
        maxOrderX_m = maxOrderX;
}
inline
    std::size_t MultipoleT::getTransMaxOrder() const {
        return transMaxOrder_m;
}
inline 
    void MultipoleT::setTransMaxOrder(std::size_t transMaxOrder) {
        transMaxOrder_m = transMaxOrder;
	transProfile_m.resize(transMaxOrder + 1, 0.);
}
inline
    double MultipoleT::getRotation() const {
         return rotation_m;
}
inline 
    void MultipoleT::setRotation(double rot) {
         rotation_m = rot;
}
inline
    void MultipoleT::setBendAngle(double angle) {
        angle_m = angle;
}
inline
    double MultipoleT::getBendAngle() const {
        return angle_m;
}
inline
    void MultipoleT::setLength(double length) {
        length_m = std::abs(length);
}
inline
    double MultipoleT::getLength() const {
        return length_m;
}
inline
    double MultipoleT::getBoundingBoxLength() const {
        return boundingBoxLength_m;
}
inline
    void MultipoleT::setBoundingBoxLength(const double &boundingBoxLength) {
        boundingBoxLength_m = boundingBoxLength;
}

#endif
