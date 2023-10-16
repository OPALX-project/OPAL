//
// This header file is made to facilitate boost's matrix and vector_t operations in OPAL.
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
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

#ifndef OPAL_MATRIX_HH
#define OPAL_MATRIX_HH

#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> matrix_t;

template<class T>
T prod_boost_vector(const boost::numeric::ublas::matrix<double>& rotation, const T& vector) {

    // Create a temporary result vector
    T result;

    // Perform matrix-vector multiplication directly without copying
    for (std::size_t i = 0; i < 3; ++i) {
        result(i) = 0;
        for (std::size_t j = 0; j < 3; ++j) {
            result(i) += rotation(i, j) * vector(j);
        }
    }

    return result;
}

#endif
