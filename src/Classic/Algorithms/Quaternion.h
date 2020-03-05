#ifndef OPAL_QUATERNION_H
#define OPAL_QUATERNION_H

#include "AppTypes/Vektor.h"
#include "Algorithms/Vektor.h"

template <class, unsigned>
class Tenzor;

class Quaternion: public Vektor<double, 4> {
public:
    Quaternion();
    Quaternion(const Quaternion &);
    Quaternion(const double &, const double &, const double &, const double &);
    Quaternion(const Vector_t &);
    Quaternion(const double &, const Vector_t &);
    Quaternion(const Tenzor<double, 3> &);

    Quaternion operator*(const double &) const;
    Quaternion operator*(const Quaternion &) const;
    Quaternion& operator=(const Quaternion &) = default;
    Quaternion& operator*=(const Quaternion &);
    Quaternion operator/(const double &) const;

    double Norm() const;
    double length() const;
    Quaternion & normalize();

    bool isUnit() const;
    bool isPure() const;
    bool isPureUnit() const;

    Quaternion inverse() const;
    Quaternion conjugate() const;

    double real() const;
    Vector_t imag() const;

    Vector_t rotate(const Vector_t &) const;

    Tenzor<double, 3> getRotationMatrix() const;
};

typedef Quaternion Quaternion_t;

Quaternion getQuaternion(Vector_t vec, Vector_t reference);


inline
Quaternion::Quaternion():
    Vektor<double, 4>(1.0, 0.0, 0.0, 0.0)
{}

inline
Quaternion::Quaternion(const Quaternion & quat):
    Vektor<double, 4>(quat)
{}

inline
Quaternion::Quaternion(const double & x0, const double & x1, const double & x2, const double & x3):
    Vektor<double, 4>(x0, x1, x2, x3)
{}

inline
Quaternion::Quaternion(const Vector_t & vec):
    Quaternion(0.0, vec(0), vec(1), vec(2))
{}

inline
Quaternion::Quaternion(const double & realPart, const Vector_t & vec):
    Quaternion(realPart, vec(0), vec(1), vec(2))
{}

inline
double Quaternion::Norm() const
{
    return dot(*this, *this);
}

inline
double Quaternion::length() const
{
    return sqrt(this->Norm());
}

inline
bool Quaternion::isUnit() const
{
    return (std::abs(this->Norm() - 1.0) < 1e-12);
}

inline
bool Quaternion::isPure() const
{
    return (std::abs((*this)(0)) < 1e-12);
}

inline
bool Quaternion::isPureUnit() const
{
    return (this->isPure() && this->isUnit());
}

inline
Quaternion Quaternion::conjugate() const
{
    Quaternion quat(this->real(), -this->imag());

    return quat;
}

inline
double Quaternion::real() const
{
    return (*this)(0);
}

inline
Vector_t Quaternion::imag() const
{
    Vector_t vec((*this)(1), (*this)(2), (*this)(3));

    return vec;
}


#endif
