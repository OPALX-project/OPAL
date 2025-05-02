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

#include "Algorithms/Vektor.h"
class MultipoleT;
class BeamlineVisitor;
class BGeometryBase;

// Base class for the implementation part of the MultipoleT magnet element
class MultipoleTBase {
public:
    /** Constructor */
    explicit MultipoleTBase(MultipoleT* element);
    /** Destructor */
    virtual ~MultipoleTBase() = default;
protected:
    MultipoleT* element_m;
public:
    /** Initialise the element */
    virtual void initialise() = 0;
    /** Return the cell geometry */
    virtual BGeometryBase& getGeometry() = 0;
    /** Return the cell geometry */
    virtual const BGeometryBase& getGeometry() const = 0;
    /** Transform to Frenet-Serret coordinates for sector magnets */
    virtual void transformCoords(Vector_t &R) = 0;
    /** Transform B-field from Frenet-Serret coordinates to lab coordinates */
    virtual void transformBField(Vector_t &B, const Vector_t &R) = 0;
    /** Returns the scale factor @f$ h_s = 1 + x / \rho(s) @f$
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getScaleFactor(double x, double s) = 0;
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
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getFn(size_t n, double x, double s) = 0;
    /** Set the number of terms used in calculation of field components \n
     *  \param orderZ -> Number of terms in expansion in z
     *  \param orderX -> Number of terms in expansion in x
     */
    virtual void setMaxOrder(size_t /*orderZ*/, size_t /*orderX*/) {}
    virtual Vector_t localCartesianToOpalCartesian(const Vector_t& r) { return r; }
    virtual double localCartesianRotation() { return 0.0; }
};

#endif
