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

class Enge : public EndFieldModel {
  public:
    /** Default constructor */
    Enge() : _a(), _lambda(0.) {SetEngeDiffIndices(10);}
    /** Builds Enge function with parameters a_0, a_1, ..., lambda and x0.
     *
     *  max_index is the maximum derivative that will be used in calculation
     *  if, after setup, you find you need to calculate higher derivatives, you
     *  can just call SetEngeDiffIndices(n) where n is the highest derivative
     *  you think you will need.
     */
    Enge(const std::vector<double> a, double x0, double lambda, int max_index)
        : _a(a), _lambda(lambda), _x0(x0) {SetEngeDiffIndices(max_index);}
    /** Destructor - no mallocs, so does nothing */
    ~Enge() {}

    /** Inheritable copy constructor - no mallocs, so does nothing */
    Enge* Clone() const;

    /** Returns the enge parameters (a_i) */
    std::vector<double> GetEngeParameters() const {return _a;}

    /** Sets the enge parameters (a_i) */
    void                SetEngeParameters(std::vector<double> a) {_a = a;}

    /** Returns the value of lambda */
    inline double       GetLambda() const {return _lambda;}

    /** Sets the value of lambda */
    inline void         SetLambda(double lambda) {_lambda = lambda;}

    /** Returns the value of x0 */
    inline double       GetX0() const {return _x0;}

    /** Sets the value of x0 */
    inline void         SetX0(double x0) {_x0 = x0;}

    /** Returns the value of the Enge function or its \f$n^{th}\f$ derivative.
     *
     *  Please call SetEngeDiffIndices(n) before calling if n > max_index
     */
    double GetEnge(double x, int n) const;

    /** Returns \f$Enge(x-x0) + Enge(-x-x0)-1\f$ and its derivatives */
    inline double GetDoubleEnge(double x, int n) const;

    /** Returns \f$h(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$h(x) = a_0 + a_1 x/\lambda + a_2 x^2/lambda^2 + \ldots \f$
     *  Please call SetEngeDiffIndices(n) before calling if n > max_index
     */
    double HN(double x, int n) const;

    /** Returns \f$g(x)\f$ or its \f$n^{th}\f$ derivative.
     *
     *  Here \f$g(x) = 1+exp(h(x))\f$.
     *  Please call SetEngeDiffIndices(n) before calling if n > max_index
     */
    double GN(double x, int n) const;

    /** Recursively calculate the indices for Enge and H
     *
     *  This will calculate the indices for Enge and H that are required to
     *  calculate the differential up to order n.
     */
    static void   SetEngeDiffIndices(size_t n);

    /** Return the indices for calculating the nth derivative of Enge ito g(x) */
    inline static std::vector< std::vector<int> > GetQIndex(int n);

    /** Return the indices for calculating the nth derivative of g(x) ito h(x) */
    inline static std::vector< std::vector<int> > GetHIndex(int n);
  private:
    std::vector<double> _a;
    double              _lambda, _x0;

    /** Indexes the derivatives of enge in terms of g */
    static std::vector< std::vector< std::vector<int> > > _q;
    /** Indexes the derivatives of g in terms of h */
    static std::vector< std::vector< std::vector<int> > > _h;
};

std::vector< std::vector<int> > Enge::GetQIndex(int n) {
  SetEngeDiffIndices(n);
  return _q[n];
}

std::vector< std::vector<int> > Enge::GetHIndex(int n) {
  SetEngeDiffIndices(n);
  return _h[n];
}

double Enge::GetDoubleEnge(double x, int n) const {
  if (n == 0) {
    return (GetEnge(x-_x0, n)+GetEnge(-x-_x0, n))-1.;
  } else {
    if (n%2 != 0) return GetEnge(x-_x0, n)-GetEnge(-x-_x0, n);
    else          return GetEnge(x-_x0, n)+GetEnge(-x-_x0, n);
  }
}

}

#endif
