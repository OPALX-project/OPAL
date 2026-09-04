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

#ifndef ENDFIELDMODEL_ENGE_H_
#define ENDFIELDMODEL_ENGE_H_

#include <iostream>
#include <vector>

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace endfieldmodel {

/** Enge class is a symmetric Enge function for analytical field models
 * 
 *  Enge function is
 * 
 *  \f$f(x) = 1/(1+exp(h(x-x0)))+1/(1+exp(h(-x-x0)))-1\f$.
 * 
 *  where h is a polynomial in x/lambda with polynomial coefficients a
 */

class Enge : public EndFieldModel {
  public:
    /** Default constructor */
    Enge() : _a(), _lambda(0.) {setEngeDiffIndices(10);}
    /** Builds Enge function with parameters a_0, a_1, ..., lambda and x0.
     * 
     *  Note that this class is in the inner loop of tracking, so many function
     *  calls are _not_ checked for correct indexing. Call setMaximumDerivative
     *  before use.
     */
    Enge(const std::vector<double> a, double x0, double lambda);

    /** Destructor - no mallocs, so does nothing */
    ~Enge() {}

    /** Inheritable copy constructor - no mallocs, so does nothing */
    Enge* clone() const;

    /** Rescale so Enge(x) -> Enge(scaleFactor*x)
     * 
     *  Sets x0 to scaleFactor*x0 and lambda to scaleFactor*lambda
     */
    void rescale(double scaleFactor);

    /** Return the value of enge(x+x0) + enge(-x-x0) at some point x */
    inline double function(double x, int n) const;

    /** Nominal end length is lambda */
    inline double getEndLength() const;

    /** Nominal centre length is x0/2 */
    inline double getCentreLength() const;

    /** Nominal start position is -x0/2 */
    inline double getStartOffset() const;

    /** Print human-readable version of enge */
    std::ostream& print(std::ostream& out) const;

    /** Returns the enge polynomial coefficients (a_i) */
    std::vector<double> getCoefficients() const {return _a;}

    /** Sets the enge polynomial coefficients (a_i) */
    void setCoefficients(std::vector<double> a) {_a = a;}

    /** Returns the value of lambda */
    inline double getLambda() const {return _lambda;}

    /** Sets the value of lambda */
    inline void setLambda(double lambda) {_lambda = lambda;}

    /** Returns the value of x0 */
    inline double getX0() const {return _x0;}

    /** Sets the value of x0 */
    inline void setX0(double x0) {_x0 = x0;}

    /** Calls setEngeDiffIndices to set the maximum derivative */
    inline void setMaximumDerivative(size_t n);

    /** Returns the value of the Enge function or its \f$n^{th}\f$ derivative.
     *
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    double getEnge(double x, int n) const;

    /** Returns \f$Enge(x-x0) + Enge(-x-x0)-1\f$ and its derivatives */
    inline double getDoubleEnge(double x, int n) const;

    /** Returns \f$h(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$h(x) = a_0 + a_1 x/\lambda + a_2 x^2/lambda^2 + \ldots \f$
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    double hN(double x, int n) const;

    /** Returns \f$g(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$g(x) = 1+exp(h(x))\f$.
     *  Please call setEngeDiffIndices(n) before calling if n > max_index
     */
    double gN(double x, int n) const;

    /** Recursively calculate the indices for Enge and H
     *
     *  This will calculate the indices for Enge and H that are required to
     *  calculate the differential up to order n.
     */
    static void   setEngeDiffIndices(size_t n);

    /** Return the indices for calculating the nth derivative of Enge ito g(x) */
    inline static std::vector< std::vector<int> > getQIndex(int n);

    /** Return the indices for calculating the nth derivative of g(x) ito h(x) */
    inline static std::vector< std::vector<int> > getHIndex(int n);
  private:
    Enge(const Enge& enge);
    Enge& operator=(const Enge& enge);
    std::vector<double> _a;
    double              _lambda, _x0;

    /** Indexes the derivatives of enge in terms of g */
    static std::vector< std::vector< std::vector<int> > > _q;
    /** Indexes the derivatives of g in terms of h */
    static std::vector< std::vector< std::vector<int> > > _h;
};

void Enge::setMaximumDerivative(size_t n) {
    Enge::setEngeDiffIndices(n);
}

double Enge::function(double x, int n) const {
    return getDoubleEnge(x, n);
}

std::vector< std::vector<int> > Enge::getQIndex(int n) {
  return _q[n];
}

std::vector< std::vector<int> > Enge::getHIndex(int n) {
  return _h[n];
}

double Enge::getDoubleEnge(double x, int n) const {
  if (n == 0) {
    return (getEnge(x-_x0, n)+getEnge(-x-_x0, n))-1.;
  } else {
    if (n%2 != 0) return getEnge(x-_x0, n)-getEnge(-x-_x0, n);
    else          return getEnge(x-_x0, n)+getEnge(-x-_x0, n);
  }
}

double Enge::getCentreLength() const {
  return _x0*2.0;
}

double Enge::getEndLength() const {
  return _lambda;
}

double Enge::getStartOffset() const {
    return _x0;
}
}

#endif
