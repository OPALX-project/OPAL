/*
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

#include <cmath>
#include <vector>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_pow_int.h"
#include "Utility/Inform.h" // ippl/src/

#include "CoordinateTransform.h"

extern Inform *gmsg;

namespace coordinatetransform {

const double CoordinateTransform::error = 1e-10;
const int CoordinateTransform::workspaceSize = 1000;
const int CoordinateTransform::algorithm = GSL_INTEG_GAUSS61;

CoordinateTransform::CoordinateTransform(const double &xlab,
                                         const double &ylab,
                                         const double &zlab,
                                         const double &s_0,
                                         const double &lambdaleft,
                                         const double &lambdaright,
                                         const double &rho):
    s_0_m(s_0), lambdaleft_m(lambdaleft),
    lambdaright_m(lambdaright), rho_m(rho) {
    std::vector<double> coordinates;
    coordinates.push_back(xlab);
    coordinates.push_back(zlab);
    //transformFromEntranceCoordinates(coordinates, boundingBoxLength);
    calcSCoordinate(coordinates[0], coordinates[1]);
    calcXCoordinate(coordinates[0], coordinates[1]);
    z_m = ylab;
}

CoordinateTransform::CoordinateTransform(const CoordinateTransform &transform):
    s_0_m(transform.s_0_m), lambdaleft_m(transform.lambdaleft_m),
    lambdaright_m(transform.lambdaright_m),
    rho_m(transform.rho_m), x_m(transform.x_m),
    z_m(transform.z_m), s_m(transform.s_m) {
}

CoordinateTransform::~CoordinateTransform() {
}

CoordinateTransform& CoordinateTransform::operator= (
                     const CoordinateTransform& transform) {
    s_0_m = transform.s_0_m;
    lambdaleft_m = transform.lambdaleft_m;
    lambdaright_m = transform.lambdaright_m;
    x_m = transform.x_m;
    z_m = transform.z_m;
    s_m = transform.s_m;
    rho_m = transform.rho_m;
    return *this;
}

std::vector<double> CoordinateTransform::getTransformation() const {
    std::vector<double> result;
    result.push_back(x_m);
    result.push_back(z_m);
    result.push_back(s_m);
    return result;
}

std::vector<double> CoordinateTransform::getUnitTangentVector(
                                         const double &s) const {
    std::vector<double> result;
    double prefactor = rho_m * (std::tanh(s_0_m / lambdaleft_m) -
                       std::tanh(-s_0_m / lambdaright_m));
    result.push_back(-std::sin((lambdaleft_m * std::log(std::cosh((s + s_0_m) / lambdaleft_m))
                     - lambdaright_m * std::log(std::cosh((s - s_0_m) / lambdaright_m)))
                     / prefactor));
    result.push_back(std::cos((lambdaleft_m * std::log(std::cosh((s + s_0_m) / lambdaleft_m))
                     - lambdaright_m * std::log(std::cosh((s - s_0_m) / lambdaright_m)))
                     / prefactor));
    return result;
}

std::vector<double> CoordinateTransform::calcReferenceTrajectory(
                                         const double &s) const {
    struct myParams params = {s_0_m, lambdaleft_m, lambdaright_m, rho_m};
    gsl_function FX, FY;
    FX.function = &getUnitTangentVectorX;
    FY.function = &getUnitTangentVectorY;
    FX.params = &params;
    FY.params = &params;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspaceSize);
    double resultX, resultY, absErrX, absErrY;
    gsl_error_handler_t* err_default = gsl_set_error_handler_off();
    int errX = gsl_integration_qag(&FX, 0, s, error, error, workspaceSize,
                                   algorithm, w, &resultX, &absErrX);
    int errY = gsl_integration_qag(&FY, 0, s, error, error, workspaceSize,
                                   algorithm, w, &resultY, &absErrY);
    if (errX || errY) {
        *gmsg << "Warning - failed to reach specified error " << error 
              << " in multipoleT coordinateTransform" << endl;
        *gmsg << "  X " << errX << " absErr: " << absErrX << " s: " << s << endl;
        *gmsg << "  Y " << errY << " absErr: " << absErrY << " s: " << s << endl;
    }
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(err_default);
    std::vector<double> result;
    result.push_back(resultX);
    result.push_back(resultY);
    return result;
}

void CoordinateTransform::calcSCoordinate(const double &xlab,
                                          const double &zlab) {
    double s1 = -(5 * lambdaleft_m + s_0_m), s2 = (5 * lambdaright_m + s_0_m);
    std::vector<double> shat1 = getUnitTangentVector(s1);
    std::vector<double> shat2 = getUnitTangentVector(s2);
    std::vector<double> r_1 = calcReferenceTrajectory(s1);
    std::vector<double> r_2 = calcReferenceTrajectory(s2);
    double eqn1 = (r_1[0] - xlab) * shat1[0] + (r_1[1] - zlab) * shat1[1];
    double eqn2 = (r_2[0] - xlab) * shat2[0] + (r_2[1] - zlab) * shat2[1];
    if (eqn1 * eqn2 > 0) {
        s_m = 10 * s_0_m;
        return;
    }
    // if (eqn1 > 0 && eqn2 < 0) {
    //     std::swap(eqn1, eqn2)
    // }
    int n = 0;
    while (n < 10000 && std::abs(s1 - s2) > 1e-12) {
        double stemp = (s2 + s1) / 2;
        std::vector<double> r_temp = calcReferenceTrajectory(stemp);
        std::vector<double> shattemp = getUnitTangentVector(stemp);
        double eqntemp = (r_temp[0] - xlab) * shattemp[0] +
                         (r_temp[1] - zlab) * shattemp[1];
        if (eqntemp > 0) {
            s2 = stemp;
        } else if (eqntemp < 0) {
            s1 = stemp;
        } else {
            s_m = stemp;
            return;
        }
        n++;
    }
    if (n == 10000) {
        s_m = 0.0;
        return;
    }
    s_m = (s2 + s1) / 2;
}

void CoordinateTransform::calcXCoordinate(const double &xlab,
                                          const double &zlab) {
    std::vector<double> shat = getUnitTangentVector(s_m);
    std::vector<double> r_0 = calcReferenceTrajectory(s_m);
    x_m = (xlab * shat[1] - zlab * shat[0] -
           r_0[0] * shat[1] + r_0[1] * shat[0]);
}

void CoordinateTransform::transformFromEntranceCoordinates(
                          std::vector<double> &coordinates,
                          const double &boundingBoxLength) {
    std::vector<double> r_entrance =
                        calcReferenceTrajectory(-boundingBoxLength);
    std::vector<double> shat = getUnitTangentVector(-boundingBoxLength);
    double x = coordinates[0], z = coordinates[1];
    coordinates[0] =  x * shat[1] + z * shat[0] + r_entrance[0];
    coordinates[1] = -x * shat[0] + z * shat[1] + r_entrance[1];
}

double getUnitTangentVectorX(double s, void *p) {
    struct myParams *params = (struct myParams *)p;
    double prefactor = params->rho * (std::tanh(params->s_0 / params->lambdaleft) -
                       std::tanh(-params->s_0 / params->lambdaright));
    return -std::sin((params->lambdaleft * std::log(std::cosh((s + params->s_0)
           / params->lambdaleft)) - params->lambdaright
           * std::log(std::cosh((s - params->s_0) / params->lambdaright))) / prefactor);
}

double getUnitTangentVectorY(double s, void *p) {
    struct myParams *params = (struct myParams *)p;
    double prefactor = params->rho * (std::tanh(params->s_0 / params->lambdaleft) -
                       std::tanh(-params->s_0 / params->lambdaright));
    return std::cos((params->lambdaleft * std::log(std::cosh((s + params->s_0)
           / params->lambdaleft)) - params->lambdaright
           * std::log(std::cosh((s - params->s_0) / params->lambdaright))) / prefactor);
}

}