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
#include <vector>
#include "gsl/gsl_sf_pow_int.h"
#include "MultipoleTCurvedVarRadius.h"
#include "MultipoleT.h"
#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"

MultipoleTCurvedVarRadius::MultipoleTCurvedVarRadius(MultipoleT* element):
    MultipoleTBase(element),
    varRadiusGeometry_m(1.0, 1.0, 1.0, 1.0, 1.0) {
}

void MultipoleTCurvedVarRadius::initialise() {
    // Record geometry information
    varRadiusGeometry_m.setElementLength(element_m->getLength() +
                                         2.0 * element_m->getEntryOffset());
    varRadiusGeometry_m.setRadius(element_m->getLength() / element_m->getBendAngle());
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    varRadiusGeometry_m.setS0(s0);
    varRadiusGeometry_m.setLambdaLeft(leftFringe);
    varRadiusGeometry_m.setLambdaRight(rightFringe);
    setMaxOrder(element_m->getMaxFOrder(), element_m->getMaxXOrder());
    // Now work out where the entry point will be in the local cartesian coordinate system
    // whose origin is at the center of the magnet.
    localCartesianEntryPoint_ = curvilinearToLocalCartesian(
        Vector_t{0.0, 0.0, element_m->getLength() / 2.0 + element_m->getEntryOffset()});
    // The tangent to the curve at this point forms the z axis of the coordinate system
    // opal addresses us in, so we can calculate the rotation required
    auto secondPoint = curvilinearToLocalCartesian(
        Vector_t{0.0, 0.0,
                 element_m->getLength() / 2.0 + element_m->getEntryOffset() + TangentStep});
    localCartesianRotation_ = -atan2(secondPoint[0] - localCartesianEntryPoint_[0],
                                    secondPoint[2] - localCartesianEntryPoint_[2]);
}

void MultipoleTCurvedVarRadius::transformCoords(Vector_t &R) {
    // Rotate Opal supplied cartesian coordinates around its origin
    auto x_rotated = R[0] * cos(localCartesianRotation_) -
                     R[2] * sin(localCartesianRotation_);
    auto z_rotated = R[0] * sin(localCartesianRotation_) +
                     R[2] * cos(localCartesianRotation_);
    // Offset to the center of the magnet
    R = {x_rotated + localCartesianEntryPoint_[0], R[1],
         z_rotated + localCartesianEntryPoint_[2]};
    // And finally into curvilinear coordinates
    R = localCartesianToCurvilinear(R);
}

Vector_t MultipoleTCurvedVarRadius::localCartesianToOpalCartesian(const Vector_t& r) {
    // Offset to the Opal origin
    auto x_offset = r[0] - localCartesianEntryPoint_[0];
    auto z_offset = r[2] - localCartesianEntryPoint_[2];
    // And rotate
    auto x_rotated = x_offset * cos(-localCartesianRotation_) -
                     z_offset * sin(-localCartesianRotation_);
    auto z_rotated = x_offset * sin(-localCartesianRotation_) +
                     z_offset * cos(-localCartesianRotation_);
    return {x_rotated, r[1], -z_rotated};
}

Vector_t MultipoleTCurvedVarRadius::localCartesianToCurvilinear(const Vector_t& r) {
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    double rho = element_m->getLength() / element_m->getBendAngle();
    coordinatetransform::CoordinateTransform t(r[0], r[1], r[2], s0, leftFringe, rightFringe, rho);
    std::vector<double> result = t.getTransformation();
    return {result[0], result[1], result[2]};
}

void MultipoleTCurvedVarRadius::transformBField(Vector_t &B, const Vector_t &R) {
    auto [s0, leftFringe, rightFringe] = element_m->getFringeField();
    double rho = element_m->getLength() / element_m->getBendAngle();
    double prefactor = rho * (tanh(s0 / leftFringe) + tanh(s0 / rightFringe));
    double theta = leftFringe * log(cosh((R[2] + s0) / leftFringe)) -
                   rightFringe * log(cosh((R[2] - s0) / rightFringe));
    theta /= prefactor;
    double Bx = B[0], Bs = B[2];
    B[0] = Bx * cos(theta) - Bs * sin(theta);
    B[2] = Bx * sin(theta) + Bs * cos(theta);
}

void MultipoleTCurvedVarRadius::setMaxOrder(size_t orderZ, size_t orderX) {
    std::size_t N = recursion_m.size();
    while (orderZ >= N) {
        polynomial::RecursionRelationTwo r(N, 2 * (N + orderX + 1));
        r.resizeX(element_m->getTransMaxOrder());
        r.truncate(orderX);
        recursion_m.push_back(r);
        N = recursion_m.size();
    }
}

double MultipoleTCurvedVarRadius::getScaleFactor(double x, double s) {
    double result;
    if (element_m->getFringeDeriv(0, s) > 1.0e-12 && element_m->getBendAngle() != 0.0) {
        double radius = element_m->getLength() * element_m->getFringeDeriv(0, 0) /
               (element_m->getFringeDeriv(0, s) * element_m->getBendAngle());
        result = 1.0 + x / radius;
    } else {
        result = 1.0;
    }
    return result;
}

double MultipoleTCurvedVarRadius::getFn(size_t n, double x, double s) {
    double result{};
    if (n == 0) {
        result = element_m->getTransDeriv(0, x) * element_m->getFringeDeriv(0, s);
    } else {
        double rho  = element_m->getLength() / element_m->getBendAngle();
        double S_0  = element_m->getFringeDeriv(0, 0);
        double y    = element_m->getFringeDeriv(0, s) / (S_0 * rho);
        std::vector<double> fringeDerivatives;
        for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
            fringeDerivatives.push_back(element_m->getFringeDeriv(j, s) / (S_0 * rho));
        }
        for (std::size_t i = 0; i <= recursion_m.at(n).getMaxXDerivatives(); i++) {
            double temp = 0.0;
            for (std::size_t j = 0; j <= recursion_m.at(n).getMaxSDerivatives(); j++) {
                temp += recursion_m.at(n).evaluatePolynomial(x, y, i, j, fringeDerivatives)
                        * fringeDerivatives.at(j);
            }
            result += temp * element_m->getTransDeriv(i, x);
        }
        result *= gsl_sf_pow_int(-1.0, static_cast<int>(n)) * S_0 * rho;
    }
    return result;
}

double MultipoleTCurvedVarRadius::reverseTransformResidual(const Vector_t& r,
                                                           const Vector_t& target) {
    // Return the distance between the vector r and the target.
    // We only consider the first and last coordinates as the height coordinate
    // is invariant across these transforms.
    auto c = localCartesianToCurvilinear(r);
    double dx = c[0] - target[0];
    double ds = c[2] - target[2];
    return sqrt(dx * dx + ds * ds);
}

Vector_t MultipoleTCurvedVarRadius::curvilinearToLocalCartesian(const Vector_t& r) {
    // This functions uses a minimize loop and coordinate descent with backtracking
    // to implement the inverse coordinate transform from the magnet's curvilinear
    // system to the local cartesian system whose origins are the centre of the magnet.
    // Note that this function is iterative and should therefore only be used occasionally.
    Vector_t result{r};
    double step = 1.0;
    double best_res = reverseTransformResidual(result, r);
    for (size_t iter = 0; iter < ReverseTransformMaxIterations; ++iter) {
        bool improved = false;
        for (int dim = 0; dim < 2; ++dim) {
            for (int dir = -1; dir <= 1; dir += 2) {
                Vector_t trial = result;
                if (dim == 0) {
                    trial[0] += dir * step;
                } else {
                    trial[2] += dir * step;
                }
                double res = reverseTransformResidual(trial, r);
                if (res < best_res) {
                    result = trial;
                    best_res = res;
                    improved = true;
                    break;
                }
            }
            if (improved) {
                break;
            }
        }
        if (!improved) {
            step *= 0.5;
        }
        if (step < ReverseTransformTolerance) {
            break;
        }
    }
    return result;
}