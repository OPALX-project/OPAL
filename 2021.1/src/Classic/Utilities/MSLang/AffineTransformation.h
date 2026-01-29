#ifndef MSLANG_AFFINETRANSFORMATION_H
#define MSLANG_AFFINETRANSFORMATION_H

#include "Algorithms/Vektor.h"
#include "AppTypes/Tenzor.h"

#include <iostream>
#include <fstream>

namespace mslang {
    struct AffineTransformation: public Tenzor<double, 3> {
        AffineTransformation(const Vector_t& row0,
                             const Vector_t& row1):
            Tenzor(row0[0], row0[1], row0[2], row1[0], row1[1], row1[2], 0.0, 0.0, 1.0) {
        }

        AffineTransformation():
            Tenzor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) { }

        AffineTransformation getInverse() const {
            AffineTransformation Ret;
            double det = (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);

            Ret(0, 0) = (*this)(1, 1) / det;
            Ret(1, 0) = -(*this)(1, 0) / det;
            Ret(0, 1) = -(*this)(0, 1) / det;
            Ret(1, 1) = (*this)(0, 0) / det;

            Ret(0, 2) = - Ret(0, 0) * (*this)(0, 2) - Ret(0, 1) * (*this)(1, 2);
            Ret(1, 2) = - Ret(1, 0) * (*this)(0, 2) - Ret(1, 1) * (*this)(1, 2);
            Ret(2, 2) = 1.0;

            return Ret;
        }

        Vector_t getOrigin() const {
            return Vector_t(-(*this)(0, 2), -(*this)(1, 2), 0.0);
        }

        double getAngle() const {
            return atan2((*this)(1, 0), (*this)(0, 0));
        }

        Vector_t transformTo(const Vector_t &v) const {
            const Tenzor<double, 3> &A = *static_cast<const Tenzor<double, 3>* >(this);
            Vector_t b(v[0], v[1], 1.0);
            Vector_t w = dot(A, b);

            return Vector_t(w[0], w[1], 0.0);
        }

        Vector_t transformFrom(const Vector_t &v) const {
            AffineTransformation inv = getInverse();
            return inv.transformTo(v);
        }

        AffineTransformation mult(const AffineTransformation &B) {
            AffineTransformation Ret;
            const Tenzor<double, 3> &A = *static_cast<const Tenzor<double, 3> *>(this);
            const Tenzor<double, 3> &BTenz = *static_cast<const Tenzor<double, 3> *>(&B);
            Tenzor<double, 3> &C = *static_cast<Tenzor<double, 3> *>(&Ret);

            C = dot(A, BTenz);

            return Ret;
        }
    };
}

#endif