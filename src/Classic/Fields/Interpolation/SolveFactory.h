/* 
 *  Copyright (c) 2015, Chris Rogers
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

#include "Fields/Interpolation/MMatrix.h"
#include "Fields/Interpolation/SquarePolynomialVector.h"

#ifndef SolveFactory_hh
#define SolveFactory_hh

namespace interpolation {

/** \class SolveFactory
 * \brief SolveFactory is a factory class for solving a set of linear
 *  equations to generate a polynomial based on nearby points.
 *
 *  See also PPSolveFactory for more details
 *
 */

class SolveFactory {
  public:
    /** Construct a new SolveFactory.
     *
     *  \param polynomial_order Order of the polynomial part of the fit
     *  \param smoothing_order Order of the smoothing part of the fit; should
     *         always be >= polynomial_order
     *  \param point_dim Dimension of the points (i.e. space of x), for
     *                    y = f(x)
     *  \param value_dim Dimension of the values (i.e. space of y), for
     *                    y = f(x)
     *  \param positions Position of the values; should be a vector of length
     *                    polynomial_order^point_dim. Each element should be a
     *                    vector of length point_dim
     *  \param deriv_positions Position of the derivatives; should be a vector
     *                    of length 
     *                    smoothing_order^point_dim - polynomial_order^point_dim.
     *                    Each element should be a vector of length point_dim.
     *  \param deriv_indices Index the derivatives. Should be a vector with
     *                    same length as deriv_positions. Each element should
     *                    be a vector of length point_dim. Each element tells us
     *                    what derivative will be defined, indexing by the order
     *                    of the derivative.
     *
     *  The constructor takes the set of positions; then performs a matrix
     *  inversion. Each call to PolynomialSolve uses the same matrix inversion
     *  in solution of the linear equations (hence when we do this lots of times
     *  as in PolynomialPatch, we get a lot faster). The matrix inversion can
     *  fail for badly formed sets of positions/deriv_positions; caveat emptor!
     */
    SolveFactory(int smoothing_order,
                 int point_dim,
                 int value_dim,
                 std::vector< std::vector<double> > positions,
                 std::vector< std::vector<double> > deriv_positions,
                 std::vector< std::vector<int> >& deriv_indices);

    /** Destructor; does nothing */
    ~SolveFactory() {}

    /** Solve to get a SquarePolynomialVector 
     * 
     *  \param values Set of y values for the solve, giving function values at
     *                 the positions defined in the constructor
     *  \param deriv_values Set of y derivatives for the solve, giving q^th
     *                 derivatives, q defined by deriv_indices, at positions
     *                 defined by deriv_positions (in constructor)
     */
    SquarePolynomialVector* PolynomialSolve(
             const std::vector< std::vector<double> >& values,
             const std::vector< std::vector<double> >& deriv_values);

    /** Convert a position vector to a set of polynomial products
     *
     *  \param position a position vector \f$ x_{b} \f$
     *
     *  Return value is a set of polynomial products \f$ x_{b}^{p} \f$
     */
    std::vector<double> MakeSquareVector(std::vector<double> position);

    /** Convert a position vector to a derivative of a set of polynomial products
     *
     *  \param position a position vector \f$ x_{b} \f$
     *  \param derivIndices indexes the derivative, \f$\vec{q} \f$
     *
     *  Return value is a set of derivatives of polynomial products
     *  \f$ p_j!/(p_j-q_j)! x_{b_j}^{p_j-q_j} \f$
     */
    std::vector<double> MakeSquareDerivVector(std::vector<double> position, std::vector<int> derivIndices);

  private:
    void BuildHInvMatrix(std::vector< std::vector<double> > positions,
                         std::vector< std::vector<double> > deriv_positions,
                         std::vector< std::vector<int> >& deriv_indices);
    int n_poly_coeffs_;
    std::vector< std::vector<int> > square_points_;
    std::vector< std::vector<int> > square_deriv_nearby_points_;

    SquarePolynomialVector square_temp_;
    MMatrix<double> h_inv_;

};
}

#endif // SolveFactory_hh


