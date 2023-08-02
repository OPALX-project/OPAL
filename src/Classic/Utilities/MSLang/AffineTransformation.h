#ifndef MSLANG_AFFINETRANSFORMATION_H
#define MSLANG_AFFINETRANSFORMATION_H

#include "Algorithms/Vektor.h"
#include "Algorithms/Matrix.h"

#include <iostream>
#include <fstream>

namespace mslang {
     struct AffineTransformation: public matrix_t {
        AffineTransformation(const Vector_t& row0,
                             const Vector_t& row1):
          matrix_t(3,3){
            (*this)(0, 0) = row0[0];
            (*this)(0, 1) = row0[1];
            (*this)(0, 2) = row0[2];
            (*this)(1, 0) = row1[0];
            (*this)(1, 1) = row1[1];
            (*this)(1, 2) = row1[2];
            (*this)(2, 0) = 0.0;
            (*this)(2, 1) = 0.0;
            (*this)(2, 2) = 1.0;
          }

        AffineTransformation():
          matrix_t(3,3){
            (*this)(0, 0) = 1.0;
            (*this)(0, 1) = 0.0;
            (*this)(0, 2) = 0.0;
            (*this)(1, 0) = 0.0;
            (*this)(1, 1) = 1.0;
            (*this)(1, 2) = 0.0;
            (*this)(2, 0) = 0.0;
            (*this)(2, 1) = 0.0;
            (*this)(2, 2) = 1.0;       
          }

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
            matrix_t b = matrix_t(3, 1); // Create a 3x1 matrix for b
            b(0, 0) = v(0);
            b(1, 0) = v(1);
            b(2, 0) = 1.0;
            matrix_t w = boost::numeric::ublas::prod(*this, b);
            return Vector_t(w(0, 0), w(1, 0), 0.0);
        }

        Vector_t transformFrom(const Vector_t &v) const {
            AffineTransformation inv = getInverse();
            return inv.transformTo(v);
        }

        AffineTransformation mult(const AffineTransformation &B) {
            AffineTransformation Ret;
            matrix_t &A = *this;
            const matrix_t &BTenz = B;
            matrix_t &C = Ret;
            C = boost::numeric::ublas::prod(A, BTenz);
            return Ret;
        }
    };
}

#endif
