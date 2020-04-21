//
// File with Boost uBLAS extensions
//   This file provides additional functions for the matrix classes
//   of uBLAS that is part of the Boost library (http://www.boost.org/).
//
// Copyright (c) 2014, Matthias Frey, ETH Zürich
// All rights reserved
//
// Implemented as part of the Semester thesis by Matthias Frey
// "Matched Distributions in Cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef MATRIX_OPERATION_H
#define MATRIX_OPERATION_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include <stdexcept>

/*!
 *  \addtogroup matt_boost
 *  @{
 */

namespace matt_boost {
    namespace ublas = boost::numeric::ublas;

    /// Computes the trace of a square matrix
    template<class V>
        BOOST_UBLAS_INLINE
        V trace(ublas::matrix<V>& e) {
            V tr = 0;
            if (e.size1() == e.size2()) {
                ublas::matrix_vector_range<ublas::matrix<V> > diag(e,ublas::range(0,e.size1()),ublas::range(0,e.size2()));
                tr = sum(diag);
            }
            else
                throw std::length_error("Error in function trace() of matrix_vector_operation.h: Wrong matrix dimensions.");
            
            return tr;
        }

    /// Computes the cross product \f$ v_{1}\times v_{2}\f$ of two vectors in \f$ \mathbb{R}^{3} \f$
    template<class V>
        BOOST_UBLAS_INLINE
        ublas::vector<V> cross_prod(ublas::vector<V>& v1, ublas::vector<V>& v2) {
            ublas::vector<V> v(v1.size());
            if (v1.size() == v2.size() && v1.size() == 3) {
                v(0) = v1(1) * v2(2) - v1(2) * v2(1);
                v(1) = v1(2) * v2(0) - v1(0) * v2(2);
                v(2) = v1(0) * v2(1) - v1(1) * v2(0);
            }
            else
                throw std::length_error("Error in function cross_prod() of matrix_vector_operation.h: Wrong vector dimensions.");
            
            return v;
        }

    /// Computes Taylor-Series of M(s) = exp(F*s)
    template<class V>
        BOOST_UBLAS_INLINE
        ublas::matrix<V> taylor_exp(const ublas::matrix<V>& F, const V ds, const unsigned int order) {
            double fac = 1.0;
            ublas::matrix<V> Fn = ublas::identity_matrix<V>(6);
            ublas::matrix<V> M = Fn;

            for (unsigned int k = 1; k < order; ++k) {
                fac *= ds / V(k);
                Fn = prod(Fn,F);
                M = M + fac * Fn;
            }
            return M;
        }

    /// Generalized matrix-matrix-matrix multiplication \f$ e_{1}\cdot e_{2}\cdot e_{3} \f$
    template<class M, class E1, class E2, class E3>
        BOOST_UBLAS_INLINE
        M gemmm(const ublas::matrix_expression<E1>& e1, const ublas::matrix_expression<E2>& e2, const ublas::matrix_expression<E3>& e3) {
            M tmp = prod(e2,e3);
            return prod(e1,tmp);
        }
}

/*! @} End of Doxygen Groups*/

#endif