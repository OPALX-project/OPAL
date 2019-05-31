/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
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


#include "MultipoleTBase.h"
#include <gsl/gsl_sf_gamma.h>
#include "AbsBeamline/MultipoleTFunctions/tanhDeriv.h"

using namespace endfieldmodel;

MultipoleTBase::MultipoleTBase():
    MultipoleTBase("")
{}

MultipoleTBase::MultipoleTBase(const std::string &name):
    Component(name),
    fringeField_l(endfieldmodel::Tanh()),
    fringeField_r(endfieldmodel::Tanh()),
    maxOrder_m(3),
    transMaxOrder_m(0),
    length_m(1.0),
    entranceAngle_m(0.0),
    rotation_m(0.0),
    boundingBoxLength_m(0.0),
    verticalApert_m(0.5),
    horizontalApert_m(0.5) {
}

MultipoleTBase::MultipoleTBase(const MultipoleTBase &right):
    Component(right),
    fringeField_l(right.fringeField_l),
    fringeField_r(right.fringeField_r),
    maxOrder_m(right.maxOrder_m),
    transMaxOrder_m(right.transMaxOrder_m),
    transProfile_m(right.transProfile_m),
    length_m(right.length_m),
    entranceAngle_m(right.entranceAngle_m),
    rotation_m(right.rotation_m),
    boundingBoxLength_m(right.boundingBoxLength_m),
    verticalApert_m(right.verticalApert_m),
    horizontalApert_m(right.horizontalApert_m),
    dummy() {
    RefPartBunch_m = right.RefPartBunch_m;
}

MultipoleTBase::~MultipoleTBase() {
}

bool MultipoleTBase::apply(const Vector_t &R, const Vector_t &P,
                           const double &t,Vector_t &E, Vector_t &B) {
    /** Rotate coordinates around the central axis of the magnet */
    Vector_t R_prime = rotateFrame(R);
    /** Go to local Frenet-Serret coordinates */
    R_prime[2] *= -1; //OPAL uses a different sign convention...
    transformCoords(R_prime);
    if (insideAperture(R_prime)) {
        /** Calculate B-field in the local Frenet-Serret frame */
        B[0] = getBx(R_prime);
        B[1] = getBz(R_prime);
        B[2] = getBs(R_prime);
        /** Transform B-field from local to lab coordinates */
        transformBField(B, R_prime);
        B[2] *= -1; //OPAL uses a different sign convention...
        return false;
    } else {
        for(int i = 0; i < 3; i++) {
            B[i] = 0.0;
        }
        return true;
    }
}

Vector_t MultipoleTBase::rotateFrame(const Vector_t &R) {
    Vector_t R_prime(3), R_pprime(3);
    /** Apply two 2D rotation matrices to coordinate vector
      * Rotate around central axis => skew fields
      * Rotate azymuthaly => entrance angle
      */
    // 1st rotation
    R_prime[0] = R[0] * cos(rotation_m) + R[1] * sin(rotation_m);
    R_prime[1] = - R[0] * sin(rotation_m) + R[1] * cos(rotation_m);
    R_prime[2] = R[2];
    // 2nd rotation
    R_pprime[0] = R_prime[2] * sin(entranceAngle_m) + 
                  R_prime[0] * cos(entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * cos(entranceAngle_m) - 
                  R_prime[0] * sin(entranceAngle_m);
    return R_pprime;

}

Vector_t MultipoleTBase::rotateFrameInverse(Vector_t &B) {
    /** This function represents the inverse of the rotation 
      * around the central axis performed by rotateFrame() method
      * Used to rotate B field back to global coordinate system
      */
    Vector_t B_prime(3);
    B_prime[0] = B[0] * cos(rotation_m) + B[1] * sin(rotation_m);
    B_prime[1] = - B[0] * sin(rotation_m) + B[1] * cos(rotation_m);
    B_prime[2] = B[2];
    return B_prime;
}

bool MultipoleTBase::setFringeField(const double &s0,
                                    const double &lambda_l,
                                    const double &lambda_r) {
    fringeField_l.Tanh::setLambda(lambda_l);
    fringeField_l.Tanh::setX0(s0);
    fringeField_l.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);
    fringeField_r.Tanh::setLambda(lambda_r);
    fringeField_r.Tanh::setX0(s0);
    fringeField_r.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);
    return true;
}

double MultipoleTBase::getBz(const Vector_t &R) {
    if (fabs(getFringeDeriv(0, R[2])) < 1.0e-12) {
        return 0.0; 
    }
    double Bz = 0.0;
    std::size_t n = getMaxOrder() + 1;
    while (n != 0) {
        n--;
        Bz = Bz * R[1] * R[1] + getFn(n, R[0], R[2])
             / gsl_sf_fact(2 * n);
    }
    return Bz;
} 

double MultipoleTBase::getBx(const Vector_t &R) {
    if (fabs(getFringeDeriv(0, R[2])) < 1.0e-12) {
        return 0.0; 
    }
    double Bx = 0.0;
    std::size_t n = getMaxOrder() + 1;
    while (n != 0) {
        n--;
        Bx = Bx * R[1] * R[1] + getFnDerivX(n, R[0], R[2])
             / gsl_sf_fact(2 * n + 1);
    }
    Bx *= R[1];
    return Bx;
} 

double MultipoleTBase::getBs(const Vector_t &R) {
    if (fabs(getFringeDeriv(0, R[2])) < 1.0e-12) {
        return 0.0; 
    }
    double Bs = 0.0;
    std::size_t n = getMaxOrder() + 1;
    while (n != 0) {
        n--;
        Bs = Bs * R[1] * R[1] + getFnDerivS(n, R[0], R[2])
             / gsl_sf_fact(2 * n + 1);
    }
    Bs *= R[1] / getScaleFactor(R[0], R[2]);
    return Bs;
}

double MultipoleTBase::getFringeDeriv(const std::size_t &n, const double &s) {
    if (n <= 10) {
        return (fringeField_l.getTanh(s, n) - fringeField_r.getNegTanh(s, n))
                / 2;
    } else {
        return tanhderiv::integrate(s,
                                    fringeField_l.Tanh::getX0(),
                                    fringeField_l.Tanh::getLambda(),
                                    fringeField_r.Tanh::getLambda(),
                                    n);
    }
}

double MultipoleTBase::getTransDeriv(const std::size_t &n, const double &x) {
    double func = 0.0;
    std::size_t transMaxOrder = getTransMaxOrder();
    if (n > transMaxOrder) {
        return func;
    }
    std::vector<double> temp = getTransProfile();
    for(std::size_t i = 1; i <= n; i++) {
        for(std::size_t j = 0; j <= transMaxOrder; j++) {
            if (j <= transMaxOrder_m - i) {
                temp[j] = temp[j + 1] * (j + 1);
            } else {
                temp[j] = 0.0;
            }
        }
    }
    std::size_t k = transMaxOrder - n + 1;
    while (k != 0) {
        k--;
        func = func * x + temp[k];
    }
    return func;
}

double MultipoleTBase::getFnDerivX(const std::size_t &n,
                                   const double &x,
                                   const double &s) {
    if (n == 0) {
        return getTransDeriv(1, x) * getFringeDeriv(0, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * getFn(n, x - 2. * stepSize, s);
    deriv += -8. * getFn(n, x - stepSize, s);
    deriv += 8. * getFn(n, x + stepSize, s);
    deriv += -1. * getFn(n, x + 2. * stepSize, s);
    deriv /= 12 * stepSize;
    return deriv;
}

double MultipoleTBase::getFnDerivS(const std::size_t &n,
                                   const double &x, 
                                   const double &s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(1, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * getFn(n, x, s - 2. * stepSize);
    deriv += -8. * getFn(n, x, s - stepSize);
    deriv += 8. * getFn(n, x, s + stepSize);
    deriv += -1. * getFn(n, x, s + 2. * stepSize);
    deriv /= 12 * stepSize;
    return deriv;
}
