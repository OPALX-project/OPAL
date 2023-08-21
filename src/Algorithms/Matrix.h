#ifndef OPAL_MATRIX_HH
#define OPAL_MATRIX_HH

#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::matrix<double> matrix_t;

template<class T>
T prod_boost_vector(boost::numeric::ublas::matrix<double> rotation, const T& vector) {
        boost::numeric::ublas::vector<double> boostVector(3);
        boostVector(0) = vector(0);
        boostVector(1) = vector(1);
        boostVector(2) = vector(2);

        boostVector = boost::numeric::ublas::prod(rotation, boostVector);

        T prodVector;
        prodVector(0) = boostVector(0);
        prodVector(1) = boostVector(1);
        prodVector(2) = boostVector(2);

        return prodVector;
}

#endif
