#ifndef OPAL_VEKTOR_HH
#define OPAL_VEKTOR_HH

#include "AppTypes/Vektor.h"

#include <complex>

typedef Vektor<double, 3> Vector_t;
class ComplexVector {
    Vector_t real_m;
    Vector_t imag_m;

public:
    ComplexVector():
        real_m(0.0),
        imag_m(0.0)
    { }

    ComplexVector(ComplexVector const& other):
        real_m(other.real_m),
        imag_m(other.imag_m)
    { }

    ComplexVector(Vector_t const& real, Vector_t const& imag):
        real_m(real),
        imag_m(imag)
    { }

    Vector_t& real() {
        return real_m;
    }
    Vector_t const& real() const{
        return real_m;
    }
    Vector_t& imag() {
        return imag_m;
    }
    Vector_t const& imag() const {
        return imag_m;
    }

    ComplexVector& operator+=(ComplexVector const& other) {
        real_m += other.real_m;
        imag_m += other.imag_m;

        return *this;
    }

    ComplexVector& operator-=(ComplexVector const& other) {
        real_m -= other.real_m;
        imag_m -= other.imag_m;

        return *this;
    }

    ComplexVector& operator=(ComplexVector const& other) {
        real_m = other.real_m;
        imag_m = other.imag_m;

        return *this;
    }

    ComplexVector& operator=(double a) {
        real_m = a;
        imag_m = 0;

        return *this;
    }
};

inline
ComplexVector operator+(ComplexVector const& a, ComplexVector const& b) {
    return ComplexVector(a.real() + b.real(), a.imag() + b.imag());
}

inline
ComplexVector operator-(ComplexVector const& a, ComplexVector const& b) {
    return ComplexVector(a.real() - b.real(), a.imag() - b.imag());
}

inline
ComplexVector operator*(std::complex<double> const& lhs, ComplexVector const& rhs) {
    return ComplexVector(lhs.real() * rhs.real() - lhs.imag() * rhs.imag(),
                         lhs.imag() * rhs.real() + lhs.real() * rhs.imag());
}

inline
std::ostream& operator<<(std::ostream &out, ComplexVector const& vec) {
    out << vec.real() << " + i * " << vec.imag();

    return out;
}

typedef ComplexVector ComplexVector_t;

#endif