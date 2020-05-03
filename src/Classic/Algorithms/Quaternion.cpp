#include "Algorithms/Quaternion.h"
#include "AppTypes/Tenzor.h"
#include "Physics/Physics.h"
#include "Utility/RandomNumberGen.h"
#include "Utilities/GeneralClassicException.h"

namespace {
    Vector_t normalize(const Vector_t & vec)
    {
        double length = std::sqrt(dot(vec, vec));

#ifndef NOPAssert
        if (length < 1e-12)
            throw GeneralClassicException("normalize()",
                                          "length of vector less than 1e-12");
#endif

        return vec / length;
    }
}

Quaternion::Quaternion(const Tenzor<double, 3> &M):
    Vektor<double, 4>(0.0)
{
    (*this)(0) = std::sqrt(std::max(0.0, 1 + M(0, 0) + M(1, 1) + M(2, 2))) / 2;
    (*this)(1) = std::sqrt(std::max(0.0, 1 + M(0, 0) - M(1, 1) - M(2, 2))) / 2;
    (*this)(2) = std::sqrt(std::max(0.0, 1 - M(0, 0) + M(1, 1) - M(2, 2))) / 2;
    (*this)(3) = std::sqrt(std::max(0.0, 1 - M(0, 0) - M(1, 1) + M(2, 2))) / 2;
    (*this)(1) = std::abs(M(2, 1) - M(1, 2)) > 0? std::copysign((*this)(1), M(2, 1) - M(1, 2)): 0.0;
    (*this)(2) = std::abs(M(0, 2) - M(2, 0)) > 0? std::copysign((*this)(2), M(0, 2) - M(2, 0)): 0.0;
    (*this)(3) = std::abs(M(1, 0) - M(0, 1)) > 0? std::copysign((*this)(3), M(1, 0) - M(0, 1)): 0.0;
}

Quaternion getQuaternion(Vector_t u, Vector_t ref)
{
    const double tol = 1e-12;

    u = normalize(u);
    ref = normalize(ref);

    Vector_t axis = cross(u, ref);
    double normAxis = std::sqrt(dot(axis,axis));

    if (normAxis < tol) {
        if (std::abs(dot(u, ref) - 1.0) < tol) {
            return Quaternion(1.0, Vector_t(0.0));
        }
        // vectors are parallel or antiparallel
        do { // find any vector in plane with ref as normal
            double u = IpplRandom();
            double v = 2 * Physics::pi * IpplRandom();
            axis(0) = std::sqrt(1 - u*u) * std::cos(v);
            axis(1) = std::sqrt(1 - u*u) * std::sin(v);
            axis(2) = u;
        } while(std::abs(dot(axis, ref)) > 0.9);

        axis -= dot(axis, ref) * ref;
        axis = normalize(axis);

        return Quaternion(0, axis);
    }

    axis /= normAxis;

    double cosAngle = sqrt(0.5 * (1 + dot(u, ref)));
    double sinAngle = sqrt(1 - cosAngle * cosAngle);

    return Quaternion(cosAngle, sinAngle * axis);
}

Quaternion Quaternion::operator*(const double & d) const
{
    Quaternion result(d * this->real(), d * this->imag());

    return result;
}

Quaternion Quaternion::operator*(const Quaternion & other) const
{
    Quaternion result(*this);
    return result *= other;
}

Quaternion & Quaternion::operator*=(const Quaternion & other)
{
    Vector_t imagThis = this->imag();
    Vector_t imagOther = other.imag();

    *this = Quaternion((*this)(0) * other(0) - dot(imagThis, imagOther),
                       (*this)(0) * imagOther + other(0) * imagThis + cross(imagThis, imagOther));

    return *this;
}

Quaternion Quaternion::operator/(const double & d) const
{
    Quaternion result(this->real() / d, this->imag() / d);

    return result;
}

Quaternion & Quaternion::normalize()
{
#ifndef NOPAssert
    if (this->Norm() < 1e-12)
        throw GeneralClassicException("Quaternion::normalize()",
                                      "length of quaternion less than 1e-12");
#endif

    (*this) /= this->length();

    return (*this);
}

Quaternion Quaternion::inverse() const
{
    Quaternion returnValue = conjugate();

    return returnValue.normalize();
}

Vector_t Quaternion::rotate(const Vector_t & vec) const
{
#ifndef NOPAssert
    if (!this->isUnit())
        throw GeneralClassicException("Quaternion::rotate()",
                                      "quaternion isn't unit quaternion. Norm: " + std::to_string(this->Norm()));
#endif

    Quaternion quat(vec);

    return ((*this) * (quat * (*this).conjugate())).imag();
}

Tenzor<double, 3> Quaternion::getRotationMatrix() const
{
    Quaternion rot(*this);
    rot.normalize();
    Tenzor<double, 3> mat(1 - 2 * (rot(2) * rot(2) + rot(3) * rot(3)),
                          2 * (-rot(0) * rot(3) + rot(1) * rot(2)),
                          2 * (rot(0) * rot(2) + rot(1) * rot(3)),
                          2 * (rot(0) * rot(3) + rot(1) * rot(2)),
                          1 - 2 * (rot(1) * rot(1) + rot(3) * rot(3)),
                          2 * (-rot(0) * rot(1) + rot(2) * rot(3)),
                          2 * (-rot(0) * rot(2) + rot(1) * rot(3)),
                          2 * (rot(0) * rot(1) + rot(2) * rot(3)),
                          1 - 2 * (rot(1) * rot(1) + rot(2) * rot(2)));

    return mat;
}