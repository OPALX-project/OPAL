/*
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

#ifndef CLASSIC_VarRadiusGeometry_HH
#define CLASSIC_VarRadiusGeometry_HH

/** ---------------------------------------------------------------------
  *
  * VarRadiusGeometry represents a Geometry with variable radius. Assuming \n
  * a Tanh model for fringe fields, the radius of curvature varies inversely \n
  * proportional with the fringe field. Such a magnet will follow the \n
  * trajectory of the reference particle.
  * The origin is defined at the centre and extends from -length / 2 to \n
  * +length / 2.
  * Transformations are calculated using a CoordinateTransformation object, \n
  * which integrates to find the reference trajectory.
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * $Author: Martin Duy Tat, Chris Rogers
  *
  * ---------------------------------------------------------------------
  *
  *
  * ---------------------------------------------------------------------
  */

#include "BeamlineGeometry/Geometry.h"
#include "Utilities/GeneralClassicException.h"
#include <algorithm>

class VarRadiusGeometry: public BGeometryBase {
public:
    /** Build VarRadiusGeometry with given length, centre radius of curvature
     *  and fringe field
     *  \param length -> Length of geometry
     *  \param rho -> Centre radius of curvature of geometry
     *  \param s_0 -> Length of central fringe field
     *  \param lambda_left -> Length of left end fringe field
     *  \param lambda_right -> Length of right end fringe field
     */
    VarRadiusGeometry(double length,
                      double rho,
                      double s_0,
                      double lambda_left,
                      double lambda_right);
    /** Copy constructor */
    VarRadiusGeometry(const VarRadiusGeometry &right);
    /** Destructor */
    virtual ~VarRadiusGeometry();
    /** Assigment operator */
    const VarRadiusGeometry &operator=(const VarRadiusGeometry &right);
    /** Arc length along the design arc */
    virtual double getArcLength() const;
    /** Get element length measured along hte design arc */
    virtual double getElementLength() const;
    /** Set arc length
     *  \param length -> Length of element
     */
    virtual void setElementLength(double length);
    /** Get centre radius of curvature */
    double getRadius() const;
    /** Set centre radius of curvature
     *  \param rho -> Central radius of curvature
     */
    void setRadius(const double &rho);
    /** Get central fringe field length */
    double getS0() const;
    /** Set central fringe field length
     *  \param s_0 -> Central fringe field length
     */
    void setS0(const double &s_0);
    /** Get left end fringe field length */
    double getLambdaLeft() const;
    /** Set left end fringe field length
     *  \param lambda_left -> Left end fringe field length
     */
    void setLambdaLeft(const double &lambda_left);
    /** Get right end fringe field length */
    double getLambdaRight() const;
    /** Set right end fringe field length
     *  \param lambda_right -> Right end fringe field length
     */
    void setLambdaRight(const double &lambda_right);
    /** Transform of the local coordinate system
     *  \param fromS -> Transform from this position
     *  \param toS -> Transform to this position
     */
    virtual Euclid3D getTransform(double fromS, double toS) const;
    /** Transform of the local coordinate system from the origin
     *  to the entrance of the element
     *  Equivalent to getTransform(0.0, getEntrance())
     */
    virtual Euclid3D getEntranceFrame() const;
    /** Transform of the local coordinate system from the origin
     *  to the exit of the element
     *  Equivalent to getTransform(0.0, getExit())
     */
    virtual Euclid3D getExitFrame() const;
private:
    double length_m;
    double rho_m;
    double s_0_m;
    double lambda_left_m;
    double lambda_right_m;
};

// inlined (trivial) member functions

inline
    VarRadiusGeometry::VarRadiusGeometry(double length,
                                         double rho,
                                         double s_0,
                                         double lambda_left,
                                         double lambda_right):
    length_m(length), rho_m(rho), s_0_m(s_0),
    lambda_left_m(lambda_left), lambda_right_m(lambda_right) {
}

inline
    VarRadiusGeometry::VarRadiusGeometry(const VarRadiusGeometry &rhs):
    BGeometryBase(rhs),
    length_m(rhs.length_m), rho_m(rhs.rho_m), s_0_m(rhs.s_0_m),
    lambda_left_m(rhs.lambda_left_m), lambda_right_m(rhs.lambda_right_m) {
}


inline
    const VarRadiusGeometry &VarRadiusGeometry::operator= (
                            const VarRadiusGeometry &rhs) {
    length_m = rhs.length_m;
    rho_m = rhs.rho_m;
    s_0_m = rhs.s_0_m;
    lambda_left_m = rhs.lambda_left_m;
    lambda_right_m = rhs.lambda_right_m;
    return *this;
}
inline
    VarRadiusGeometry::~VarRadiusGeometry() {
}
inline
    double VarRadiusGeometry::getArcLength() const {
        return length_m;
}
inline
    double VarRadiusGeometry::getElementLength() const {
        return length_m;
}
inline
    void VarRadiusGeometry::setElementLength(double length) {
        if (length < 0.0) {
            throw GeneralClassicException("VarRadiusGeometry::setElementLength",
                                          "The length of an element has to be positive");
        }
        length_m = std::max(0.0, length);
}
inline
    double VarRadiusGeometry::getRadius() const {
        return rho_m;
}
inline
    void VarRadiusGeometry::setRadius(const double &rho) {
        rho_m = rho;
}
inline
    double VarRadiusGeometry::getS0() const {
        return s_0_m;
}
inline
    void VarRadiusGeometry::setS0(const double &s_0) {
        s_0_m = s_0;
}
inline
    double VarRadiusGeometry::getLambdaLeft() const {
        return lambda_left_m;
}
inline
    void VarRadiusGeometry::setLambdaLeft(const double &lambda_left) {
        lambda_left_m = lambda_left;
}
inline
    double VarRadiusGeometry::getLambdaRight() const {
        return lambda_right_m;
}
inline
    void VarRadiusGeometry::setLambdaRight(const double &lambda_right) {
        lambda_right_m = lambda_right;
}
inline
    Euclid3D VarRadiusGeometry::getEntranceFrame() const {
        return getTransform(getOrigin(), getEntrance());
}
inline
    Euclid3D VarRadiusGeometry::getExitFrame() const {
        return getTransform(getOrigin(), getExit());
}

#endif