#ifndef CLASSIC_FMonomial_HH
#define CLASSIC_FMonomial_HH

// ------------------------------------------------------------------------
// $RCSfile: FMonomial.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FMonomial<N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include <algorithm>
#include <numeric>


// Class FMonomial
// ------------------------------------------------------------------------
/// Representation of the exponents for a monomial with fixed dimension.

template <int N>
class FMonomial {

public:

    /// Constructor, defines variable var.
    explicit FMonomial(int var);

    FMonomial();
    FMonomial(const FMonomial &);
    ~FMonomial();
    const FMonomial &operator=(const FMonomial &);

    /// Return reference to exponent.
    int &operator[](int index);

    /// Return value of exponent.
    int operator[](int index) const;

    /// Get exponent set of a product.
    FMonomial operator*(const FMonomial &rhs) const;

    /// Compute the monomial's order.
    int getOrder() const;

    /// Get the monomial's number of variables.
    int getVariables() const;

private:

    // The monomial's exponents.
    int array[N];
};


// Class FMonomial
// ------------------------------------------------------------------------

template <int N>
FMonomial<N>::FMonomial() {
    std::fill(array + 0, array + N, 0u);
}


template <int N>
FMonomial<N>::FMonomial(const FMonomial &rhs) {
    std::copy(rhs.array + 0, rhs.array + N, array + 0);
}


template <int N>
FMonomial<N>::FMonomial(int var) {
    std::fill(array + 0, array + N, 0u);
    if(var < N) array[var] = 1;
}


template <int N>
FMonomial<N>::~FMonomial()
{}


template <int N>
const FMonomial<N> &FMonomial<N>::operator=(const FMonomial &rhs) {
    if(&rhs != this) {
        std::copy(rhs.array + 0, rhs.array + N, array + 0);
    }

    return *this;
}


template <int N>
int &FMonomial<N>::operator[](int index) {
    return array[index];
}


template <int N>
int FMonomial<N>::operator[](int index) const {
    return array[index];
}


template <int N>
FMonomial<N> FMonomial<N>::operator*(const FMonomial &rhs) const {
    FMonomial result;
    for(int i = 0; i < N; ++i) result.array[i] = array[i] + rhs.array[i];
    return result;
}


template <int N>
int FMonomial<N>::getOrder() const {
    return std::accumulate(array + 0, array + N, 0);
}


template <int N>
int FMonomial<N>::getVariables() const {
    return N;
}

#endif // CLASSIC_FMonomial_HH
