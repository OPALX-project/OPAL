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

#include <gsl/gsl_sf_pow_int.h>

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/PPSolveFactory.h"
#include "Fields/Interpolation/SolveFactory.h"

namespace interpolation {

    SolveFactory::SolveFactory(int /*polynomial_order*/,
                           int smoothing_order,
                           int point_dim,
                           int value_dim,
                           std::vector< std::vector<double> > positions,
                           std::vector< std::vector<double> > deriv_positions,
                           std::vector< std::vector<int> >& deriv_indices)
// :FIXME: unused!
// : polynomial_order_(polynomial_order), smoothing_order_(smoothing_order)
{
    n_poly_coeffs_ = SquarePolynomialVector::NumberOfPolynomialCoefficients(point_dim, smoothing_order);
    square_points_ = PPSolveFactory::getNearbyPointsSquares(point_dim, -1, smoothing_order);
    square_deriv_nearby_points_ = PPSolveFactory::getNearbyPointsSquares(point_dim, -1, smoothing_order);
    BuildHInvMatrix(positions, deriv_positions, deriv_indices);
    MMatrix<double> A_temp(value_dim, n_poly_coeffs_, 0.);
    square_temp_ = SquarePolynomialVector(point_dim, A_temp);
}

void SolveFactory::BuildHInvMatrix(
                           std::vector< std::vector<double> > positions,
                           std::vector< std::vector<double> > deriv_positions,
                           std::vector< std::vector<int> >& deriv_indices) {
    int nCoeffs = positions.size();
    h_inv_ = MMatrix<double>(n_poly_coeffs_, n_poly_coeffs_, 0.);
    for (int i = 0; i < nCoeffs; ++i) {
        std::vector<double> poly_vec = MakeSquareVector(positions[i]);
        for (int j = 0; j < n_poly_coeffs_; ++j) {
            h_inv_(i+1, j+1) = poly_vec[j];
        }
    }
    for (size_t i = 0; i < deriv_positions.size(); ++i) {
        std::vector<double> deriv_vec = MakeSquareDerivVector(deriv_positions[i],
                                                              deriv_indices[i]);
        for (int j = 0; j < n_poly_coeffs_; ++j) {
            h_inv_(i+1+nCoeffs, j+1) = deriv_vec[j];
        }
    }
    h_inv_.invert();
}

std::vector<double> SolveFactory::MakeSquareVector(std::vector<double> x) {
    std::vector<double> square_vector(square_points_.size(), 1.);
    for (size_t i = 0; i < square_points_.size(); ++i) {
        std::vector<int>& point = square_points_[i];
        for (size_t j = 0; j < point.size(); ++j)
            square_vector[i] *= gsl_sf_pow_int(x[j], point[j]);
    }
    return square_vector;
}

std::vector<double> SolveFactory::MakeSquareDerivVector(std::vector<double> positions, std::vector<int> deriv_indices) {
    std::vector<double> deriv_vec(square_deriv_nearby_points_.size(), 1.);
    int square_deriv_nearby_points_size = square_deriv_nearby_points_.size();
    int dim = square_deriv_nearby_points_[0].size();
    for (int i = 0; i < square_deriv_nearby_points_size; ++i) {
        for (int j = 0; j < dim; ++j) {
            int power = square_deriv_nearby_points_[i][j] - deriv_indices[j]; // p_j - q_j
            if (power < 0) {
                deriv_vec[i] = 0.;
                break;
            } else {
                // x^(p_j-q_j)
                deriv_vec[i] *= gsl_sf_pow_int(positions[j], power);
            }
            // p_j*(p_j-1)*(p_j-2)*...*(p_j-q_j)
            for (int k = square_deriv_nearby_points_[i][j]; k > power; --k) {
                deriv_vec[i] *= k;
            }
        }
    }
    return deriv_vec;
}

SquarePolynomialVector* SolveFactory::PolynomialSolve(
         const std::vector< std::vector<double> >& values,
         const std::vector< std::vector<double> >& deriv_values) {
    // Algorithm:
    // define G_i = vector of F_i and d^pF/dF^p with values taken from coeffs
    // and derivs respectively
    // define H_ij = vector of f_i and d^pf/df^p)
    // Then use G = H A => A = H^-1 G
    // where A is vector of polynomial coefficients
    // First set up G_i from coeffs and derivs; then calculate H_ij;
    // then do the matrix inversion
    // It is an error to include the same polynomial index more than once in
    // coeffs or derivs or both; any polynomial indices that are missing will be
    // assigned 0 value; polynomial order will be taken as the maximum
    // polynomial order from coeffs and derivs.
    // PointDimension and ValueDimension will be taken from coeffs and derivs;
    // it is an error if these do not all have the same dimensions.
    
    // OPTIMISATION - if we are doing this many times and only changing values,
    // can reuse H
    int nCoeffs = values.size();
    int nDerivs = deriv_values.size();
    if (values.size()+deriv_values.size() != size_t(n_poly_coeffs_)) {
        throw GeneralClassicException(
              "SolveFactory::PolynomialSolve",
              "Values and derivatives over or under constrained"
        );
    }
    int valueDim = 0;
    if (values.size() != 0) {
        valueDim = values[0].size();
    } else if (deriv_values.size() != 0) {
        valueDim = deriv_values[0].size();
    }

    MMatrix<double> A(valueDim, n_poly_coeffs_, 0.);
    for (size_t y_index = 0; y_index < values[0].size(); ++y_index) {
        MVector<double> G(n_poly_coeffs_, 0.);
        // First fill the values
        for (int i = 0; i < nCoeffs && i < n_poly_coeffs_; ++i) {
            G(i+1) = values[i][y_index];
        }
        // Now fill the derivatives
        for (int i = 0; i < nDerivs; ++i) {
            G(i+nCoeffs+1) = deriv_values[i][y_index];
        }
        MVector<double> A_vec = h_inv_*G;
        for (int j = 0; j < n_poly_coeffs_; ++j) {
            A(y_index+1, j+1) = A_vec(j+1);
        }
    }
    SquarePolynomialVector* temp = new SquarePolynomialVector(square_temp_);
    temp->SetCoefficients(A);
    return temp;
}
}

