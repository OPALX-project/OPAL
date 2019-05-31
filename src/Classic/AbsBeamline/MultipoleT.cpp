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


#include "Algorithms/PartBunch.h"
#include "MultipoleT.h"
#include <gsl/gsl_math.h>
#include "MultipoleTFunctions/RecursionRelation.h"
#include "MultipoleTFunctions/RecursionRelationTwo.h"
#include "MultipoleTFunctions/tanhDeriv.h"
#include "MultipoleTFunctions/CoordinateTransform.h"

using namespace endfieldmodel;

MultipoleT::MultipoleT(const std::string &name):
    Component(name),
    fringeField_l(endfieldmodel::Tanh()),
    fringeField_r(endfieldmodel::Tanh()),
    maxOrder_m(5),
    maxOrderX_m(10),
    transMaxOrder_m(0),
    planarArcGeometry_m(1., 1.),
    length_m(1.0),
    angle_m(0.0),
    entranceAngle_m(0.0),
    rotation_m(0.0),
    variableRadius_m(false),
    boundingBoxLength_m(0.0),
    verticalApert_m(0.5),
    horizApert_m(0.5) {
}

MultipoleT::MultipoleT(const MultipoleT &right):
    Component(right),
    fringeField_l(right.fringeField_l),
    fringeField_r(right.fringeField_r),
    maxOrder_m(right.maxOrder_m),
    maxOrderX_m(right.maxOrderX_m),
    recursion_VarRadius_m(right.recursion_VarRadius_m),
    recursion_ConstRadius_m(right.recursion_ConstRadius_m),
    transMaxOrder_m(right.transMaxOrder_m),
    transProfile_m(right.transProfile_m),
    planarArcGeometry_m(right.planarArcGeometry_m),
    length_m(right.length_m),
    angle_m(right.angle_m),
    entranceAngle_m(right.entranceAngle_m),
    rotation_m(right.rotation_m),
    variableRadius_m(right.variableRadius_m),
    boundingBoxLength_m(right.boundingBoxLength_m),
    verticalApert_m(right.verticalApert_m),
    horizApert_m(right.horizApert_m),
    dummy() {
    RefPartBunch_m = right.RefPartBunch_m;
}

MultipoleT::~MultipoleT() {
}

ElementBase* MultipoleT::clone() const {
    MultipoleT* newMultipole = new MultipoleT(*this);
    newMultipole->initialise();
    return newMultipole;
}

void MultipoleT::finalise() {
    RefPartBunch_m = NULL;
}

bool MultipoleT::apply(const Vector_t &R, const Vector_t &P,
               const double &t,Vector_t &E, Vector_t &B) {
    /** Rotate coordinates around the central axis of the magnet */
    Vector_t R_prime = rotateFrame(R);
    /** If magnet is not straight go to local Frenet-Serret coordinates */
    //R_prime[2] *= -1; //OPAL uses different sign convention...
    Vector_t X = R_prime;
    R_prime = transformCoords(X);
    if (insideAperture(R_prime)) {
        /** Transform B-field from local to lab coordinates */
        double theta;
        double prefactor = (length_m / angle_m) * (tanh((fringeField_l.getX0())
                           / fringeField_l.getLambda()) +
                           tanh((fringeField_r.getX0())
                           / fringeField_r.getLambda()));
        if (angle_m == 0.0) {
            theta = 0.0;
        } else if (!variableRadius_m) {
            theta = R_prime[2] * angle_m / length_m;
        } else { // variableRadius_m == true
            theta = fringeField_l.getLambda() * log(cosh((R_prime[2] +
                    fringeField_l.getX0()) / fringeField_l.getLambda())) -
                    fringeField_r.getLambda() * log(cosh((R_prime[2] -
                    fringeField_r.getX0()) / fringeField_r.getLambda()));
            theta /= prefactor;
        }
        double Bx = getBx(R_prime);
        double Bs = getBs(R_prime);
        B[0] = Bx * cos(theta) - Bs * sin(theta);
        B[2] = Bx * sin(theta) + Bs * cos(theta);
        B[1] = getBz(R_prime);
        // B[2] *= -1; //OPAL uses different sign convention
        return false;
    } else {
        for(int i = 0; i < 3; i++) {
            B[i] = 0.0;
        }
        return true;
    }
}

bool MultipoleT::apply(const size_t &i, const double &t,
               Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

Vector_t MultipoleT::rotateFrame(const Vector_t &R) {
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

Vector_t MultipoleT::rotateFrameInverse(Vector_t &B) {
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

bool MultipoleT::insideAperture(const Vector_t &R) {
    if (std::abs(R[1]) <= verticalApert_m / 2. && std::abs(R[0]) <= horizApert_m / 2.) {
        return true;
    }
    else {
        return false;
    }
}

Vector_t MultipoleT::transformCoords(const Vector_t &R) {
    Vector_t X(3);
    // If radius is constant
    if (!variableRadius_m) {
        double radius = length_m / angle_m;
        // Transform from Cartesian to Frenet-Serret along the magnet
        double alpha = atan(R[2] / (R[0] + radius ));
        if (alpha != 0.0 && angle_m != 0.0) {
            X[0] = R[2] / sin(alpha) - radius;
            X[1] = R[1];
            X[2] = radius * alpha;// + boundingBoxLength_m;
        } else {
            X[0] = R[0];
            X[1] = R[1];
            X[2] = R[2];// + boundingBoxLength_m;
        }
    } else {
        // If radius is variable
        coordinatetransform::CoordinateTransform t(R[0], R[1], R[2],
                                                   fringeField_l.getX0(),
                                                   fringeField_l.getLambda(),
                                                   fringeField_r.getLambda(),
                                                   (length_m / angle_m));
        std::vector<double> r = t.getTransformation();
        X[0] = r[0];
        X[1] = r[1];
        X[2] = r[2];
    }
    return X;
}

void MultipoleT::setMaxOrder(std::size_t maxOrder) {
    if (variableRadius_m && angle_m != 0.0) {
        std::size_t N = recursion_VarRadius_m.size();
        while (maxOrder >= N) {
            polynomial::RecursionRelationTwo r(N, 2 * (N + maxOrderX_m + 1));
            r.resizeX(transMaxOrder_m);
            r.truncate(maxOrderX_m);
            recursion_VarRadius_m.push_back(r);
            N = recursion_VarRadius_m.size();
        }
    } else if (!variableRadius_m && angle_m != 0.0) {
        std::size_t N = recursion_ConstRadius_m.size();
        while (maxOrder >= N) {
            polynomial::RecursionRelation r(N, 2 * (N + maxOrderX_m + 1));
            r.resizeX(transMaxOrder_m);
            r.truncate(maxOrderX_m);
            recursion_ConstRadius_m.push_back(r);
            N = recursion_ConstRadius_m.size();
        }
    }
    maxOrder_m = maxOrder; 
}

void MultipoleT::setTransProfile(std::size_t n, double dTn) {
    if (n > transMaxOrder_m) {
      transMaxOrder_m = n;
      transProfile_m.resize(n+1, 0.0);
    }
    transProfile_m[n] = dTn;
}

bool MultipoleT::setFringeField(double s0, double lambda_l, double lambda_r) {
    fringeField_l.Tanh::setLambda(lambda_l);
    fringeField_l.Tanh::setX0(s0);
    fringeField_l.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);

    fringeField_r.Tanh::setLambda(lambda_r);
    fringeField_r.Tanh::setX0(s0);
    fringeField_r.Tanh::setTanhDiffIndices(2 * maxOrder_m + 1);

    return true;
}

double MultipoleT::getBz(const Vector_t &R) {
    /** Returns the vertical field component
      * sum_n  f_n * z^(2n) / (2n)!
      */
    double Bz = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> Use corresponding field expansion directly
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double f_n = 0.0;
            for(std::size_t i = 0; i <= n; i++)  {
                f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, R[0]) * 
                       getFringeDeriv(2 * n - 2 * i, R[2]);
            }
            f_n *= gsl_sf_pow_int(-1.0, n);
            Bz += gsl_sf_pow_int(R[1], 2 * n) / gsl_sf_fact(2 * n) * f_n;
        }
    } else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // Return 0 if end of fringe field is reached
            // This is to avoid functions being called at infinite radius
            return 0.0; 
        }
        // Curved geometry -> Use full machinery of differential operators
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double f_n = getFn(n, R[0], R[2]);
            Bz += gsl_sf_pow_int(R[1], 2 * n) / gsl_sf_fact(2 * n) * f_n;
        }
    }
    return Bz;
}

double MultipoleT::getBx(const Vector_t &R) {
    /** Returns the radial component of the field
      * sum_n z^(2n+1) / (2n+1)! * \partial_x f_n 
      */
    double Bx = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> Use corresponding field expansion directly
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double f_n = 0.0;
            for(std::size_t i = 0; i <= n; i++) {
                f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i + 1, R[0]) * 
                   getFringeDeriv(2 * n - 2 * i, R[2]);
            }
            f_n *= gsl_sf_pow_int(-1.0, n);
            Bx += gsl_sf_pow_int(R[1], 2 * n + 1) / 
                  gsl_sf_fact(2 * n + 1) * f_n;
        }
    } else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // Return 0 if end of fringe field is reached
            // This is to avoid functions being called at infinite radius
            return 0.0; 
        }
        // Curved geometry -> Use full machinery of differential operators
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double partialX_fn = getFnDerivX(n, R[0], R[2]);
            Bx += partialX_fn * gsl_sf_pow_int(R[1], 2 * n + 1)
                  / gsl_sf_fact(2 * n + 1); 
        }
    }
    return Bx;
} 

double MultipoleT::getBs(const Vector_t &R) {
    /** Returns the component of the field along the central axis 
      * 1/h_s * sum_n z^(2n+1) / (2n+1)! \partial_s f_n
      */
    double Bs = 0.0;
    if (angle_m == 0.0) {
        // Straight geometry -> Use corresponding field expansion directly
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double f_n = 0.0;
            for(std::size_t i = 0; i <= n; i++) {
                f_n += gsl_sf_choose(n, i) * getTransDeriv(2 * i, R[0]) *
                       getFringeDeriv(2 * n - 2 * i + 1, R[2]);
            }
            f_n *= gsl_sf_pow_int(-1.0, n);
            Bs += gsl_sf_pow_int(R[1], 2 * n + 1) / 
                      gsl_sf_fact(2 * n + 1) * f_n;
        }
    } else {
        if (variableRadius_m == true and getFringeDeriv(0, R[2]) < 1.0e-12) {
            // Return 0 if end of fringe field is reached
            // This is to avoid functions being called at infinite radius
            return 0.0; 
        }
        // Curved geometry -> Use full machinery of differential operators
        for(std::size_t n = 0; n <= maxOrder_m; n++) {
            double partialS_fn = getFnDerivS(n, R[0], R[2]);
            Bs += partialS_fn * gsl_sf_pow_int(R[1], 2 * n + 1)
                  / gsl_sf_fact(2 * n + 1);
        }
        Bs /= getScaleFactor(R[0], R[2]);
    }
    return Bs;
}

double MultipoleT::getFringeDeriv(int n, double s) {
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

double MultipoleT::getTransDeriv(std::size_t n, double x) {
    /** Sets a vector of the coefficients in the polynomial expansion
      * of transverse profile; shifts them to the left and multiply by
      * corresponding power each time to take derivative once; 
      * repeats until desired derivative is reached
      */   
    double func = 0;
    std::vector<double> temp = transProfile_m;
    if (n <= transMaxOrder_m) {
        if (n != 0) {
            for(std::size_t i = 1; i <= n; i++) {
                for(std::size_t j = 0; j <= transMaxOrder_m; j++) {
                    if (j <= transMaxOrder_m - i) {
                        // Move terms to the left and multiply by power
                        temp[j] = temp[j + 1] * (j + 1);
                        } else {
                        // Put 0 at the end for missing higher powers
                        temp[j] = 0.0;
                    }
                }
            }
        }
        // Now use the vector to calculate value of the function
        for(std::size_t k = 0; k <= transMaxOrder_m; k++) {
            func += temp[k] * gsl_sf_pow_int(x, k);
        }
    }
    return func;
}

void MultipoleT::setDipoleConstant(double B0) {
    if (transMaxOrder_m < 1) {
        transProfile_m.resize(1, 0.);
    }
    transProfile_m[0] = B0;
}

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    visitor.visitMultipoleT(*this);
}

void MultipoleT::getDimensions(double &zBegin, double &zEnd) const {
}

void MultipoleT::setAperture(double vertAp, double horizAp) {
    verticalApert_m = vertAp;
    horizApert_m = horizAp;
}

std::vector<double> MultipoleT::getAperture() const {
    std::vector<double> temp(2, 0.0);   
    temp[0] = verticalApert_m;
    temp[1] = horizApert_m;
    return temp;
}

std::vector<double> MultipoleT::getFringeLength() const {
    std::vector<double> temp(2, 0.0);
    temp[0] = fringeField_l.getLambda();
    temp[1] = fringeField_r.getLambda();
    return temp;
}

double MultipoleT::getRadius(double s) {
    double centralRadius = length_m / angle_m;
    if (!variableRadius_m) {
        return centralRadius;
    }
    if (getFringeDeriv(0, s) != 0.0) {
       return centralRadius * getFringeDeriv(0, 0) / getFringeDeriv(0, s);
    } else {
       return -1; // Return -1 if radius is infinite 
    }
}

double MultipoleT::getScaleFactor(double x, double s) {
    if (!variableRadius_m) {
        return 1 + x * angle_m/ length_m;
    } else {
        return 1 + x / getRadius(s);
    }
}

double MultipoleT::getFnDerivX(std::size_t n, double x, double s) {
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

double MultipoleT::getFnDerivS(std::size_t n, double x, double s) {
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

double MultipoleT::getFn(std::size_t n, double x, double s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(0, s);
    }
    if (!variableRadius_m) {
        double rho = length_m / angle_m;
        double func = 0.0;
        
        for (std::size_t j = 0;
             j <= recursion_ConstRadius_m.at(n).getMaxSDerivatives();
             j++) {
            double FringeDerivj = getFringeDeriv(2 * j, s);
            for (std::size_t i = 0;
                 i <= recursion_ConstRadius_m.at(n).getMaxXDerivatives();
                 i++) {
                if (recursion_ConstRadius_m.at(n).isPolynomialZero(i, j)) {
                    continue;
                }
                func += (recursion_ConstRadius_m.at(n)
                        .evaluatePolynomial(x / rho, i, j) *
                        getTransDeriv(i, x) * FringeDerivj) /
                        gsl_sf_pow_int(rho, 2 * n - i - 2 * j);
            }
        }
        func *= gsl_sf_pow_int(-1.0, n);
        return func;
    } else {
        double rho = length_m / angle_m;
        double S_0 = getFringeDeriv(0, 0);
        double y = getFringeDeriv(0, s) / (S_0 * rho);
        double func = 0.0;
        std::vector<double> fringeDerivatives;
        for (std::size_t j = 0;
             j <= recursion_VarRadius_m.at(n).getMaxSDerivatives();
             j++) {
            fringeDerivatives.push_back(getFringeDeriv(j, s) / (S_0 * rho));
        }
        for (std::size_t i = 0;
             i <= recursion_VarRadius_m.at(n).getMaxXDerivatives();
             i++) {
            double temp = 0.0;
            for (std::size_t j = 0;
                 j <= recursion_VarRadius_m.at(n).getMaxSDerivatives();
                 j++) {
                temp += recursion_VarRadius_m.at(n)
                        .evaluatePolynomial(x, y, i, j, fringeDerivatives)
                        * fringeDerivatives.at(j);
            }
            func += temp * getTransDeriv(i, x);
        }
        func *= gsl_sf_pow_int(-1.0, n) * S_0 * rho;
        return func;
    }
}

void MultipoleT::initialise() {
    planarArcGeometry_m.setElementLength(2 * boundingBoxLength_m);
    planarArcGeometry_m.setCurvature(angle_m / length_m);
}

void MultipoleT::initialise(PartBunchBase<double, 3>* bunch, 
                double &startField, 
                double &endField) {
    RefPartBunch_m = bunch;
    initialise();
}

bool MultipoleT::bends() const {
    return (transProfile_m[0] != 0);
}

PlanarArcGeometry& MultipoleT::getGeometry() {
    return planarArcGeometry_m;
}

const PlanarArcGeometry& MultipoleT::getGeometry() const {
    return planarArcGeometry_m;
}

EMField &MultipoleT::getField() {
    return dummy;
}

const EMField &MultipoleT::getField() const {
    return dummy;
}
