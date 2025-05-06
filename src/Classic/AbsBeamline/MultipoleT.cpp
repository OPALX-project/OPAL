/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
 *  Copyright (c) 2025, Jon Thompson
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

#include "MultipoleT.h"
#include "AbsBeamline/MultipoleTFunctions/tanhDeriv.h"
#include "MultipoleTStraight.h"
#include "MultipoleTCurvedConstRadius.h"
#include "MultipoleTCurvedVarRadius.h"

using namespace endfieldmodel;

MultipoleT::MultipoleT(const std::string& name)
    : Component(name) {
    chooseImplementation();
}

MultipoleT::MultipoleT(const MultipoleT& right)
    : Component(right),
      fringeField_l(right.fringeField_l),
      fringeField_r(right.fringeField_r),
      maxFOrder_m(right.maxFOrder_m),
      maxXOrder_m(right.maxXOrder_m),
      transProfile_m(right.transProfile_m),
      transMaxOrder_m(right.transMaxOrder_m),
      length_m(right.length_m),
      entranceAngle_m(right.entranceAngle_m),
      rotation_m(right.rotation_m),
      bendAngle_m(right.bendAngle_m),
      variableRadius_m(right.variableRadius_m),
      boundingBoxLength_m(right.boundingBoxLength_m),
      verticalApert_m(right.verticalApert_m),
      horizontalApert_m(right.horizontalApert_m) {
    RefPartBunch_m = right.RefPartBunch_m;
    chooseImplementation();
}

ElementBase* MultipoleT::clone() const {
    return new MultipoleT(*this);
}

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    visitor.visitMultipoleT(*this);
}

Vector_t MultipoleT::rotateFrame(const Vector_t &R) const {
    Vector_t R_prime(3), R_pprime(3);
    /** Apply two 2D rotation matrices to coordinate vector
      * Rotate around central axis => skew fields
      * Rotate azimuthally => entrance angle
     */
    // 1st rotation
    R_prime[0] = R[0] * std::cos(rotation_m) + R[1] * std::sin(rotation_m);
    R_prime[1] = - R[0] * std::sin(rotation_m) + R[1] * std::cos(rotation_m);
    R_prime[2] = R[2];
    // 2nd rotation
    R_pprime[0] = R_prime[2] * std::sin(entranceAngle_m) +
                  R_prime[0] * std::cos(entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * std::cos(entranceAngle_m) -
                  R_prime[0] * std::sin(entranceAngle_m);
    return R_pprime;

}

bool MultipoleT::insideAperture(const Vector_t& R) const {
    return std::abs(R[1]) <= verticalApert_m / 2.0 &&
            std::abs(R[0]) <= horizontalApert_m / 2.0;
}

bool MultipoleT::insideBoundingBox(const Vector_t& R) const {
    return boundingBoxLength_m == 0.0 || fabs(R[2]) <= boundingBoxLength_m / 2.0;
}

Vector_t MultipoleT::toMagnetCoords(const Vector_t& R) {
    /** Rotate coordinates around the central axis of the magnet */
    /* TODO:  I'm not sure this is the correct thing to do */
    Vector_t result = rotateFrame(R);
    /** Go to local Frenet-Serret coordinates */
    result[2] *= -1;   // OPAL uses a different sign convention...
    implementation_->transformCoords(result);
    return result;
}

Vector_t MultipoleT::getField(const Vector_t& magnetCoords) {
    Vector_t result;
    /** Calculate B-field in the local Frenet-Serret frame */
    result[0] = implementation_->getBx(magnetCoords);
    result[1] = implementation_->getBz(magnetCoords);
    result[2] = implementation_->getBs(magnetCoords);
    /** Transform B-field from local to lab coordinates */
    implementation_->transformBField(result, magnetCoords);
    result[2] *= -1;  // OPAL uses a different sign convention...
    return result;
}

bool MultipoleT::apply(
    const Vector_t& R, const Vector_t& /*P*/, const double& /*t*/, Vector_t& /*E*/, Vector_t& B) {
    Vector_t R_prime = toMagnetCoords(R);
    bool result;
    if (insideAperture(R_prime)) {
        if(insideBoundingBox(R_prime)) {
            B = getField(R_prime);
        } else {
            B = {0.0, 0.0, 0.0};
        }
        result = false;
    } else {
        B = {0.0, 0.0, 0.0};
        result = getFlagDeleteOnTransverseExit();
    }
    return result;
}

void MultipoleT::setFringeField(const double& s0, const double& lambda_l, const double& lambda_r) {
    fringeField_l.setLambda(lambda_l);
    fringeField_l.setX0(s0);
    fringeField_r.setLambda(lambda_r);
    fringeField_r.setX0(s0);
    Tanh::setTanhDiffIndices(2 * maxFOrder_m + 1);  // Static table builds highest order encountered
    implementation_->initialise();
}

std::tuple<double, double, double> MultipoleT::getFringeField() const {
    return {fringeField_l.getX0(), fringeField_l.getLambda(), fringeField_r.getLambda()};
}

double MultipoleT::getFringeDeriv(const std::size_t& n, const double& s) {
    double result;
    if (n <= 10) {
        result = (fringeField_l.getTanh(s, static_cast<int>(n)) -
                fringeField_r.getNegTanh(s, static_cast<int>(n))) / 2;
    } else {
        result = tanhderiv::integrate(
            s, fringeField_l.getX0(), fringeField_l.getLambda(),
            fringeField_r.getLambda(), static_cast<int>(n));
    }
    return result;
}

double MultipoleT::getTransDeriv(const std::size_t& n, const double& x) const {
    double func               = 0.0;
    std::size_t transMaxOrder = getTransMaxOrder();
    if (n > transMaxOrder) {
        return func;
    }
    std::vector<double> temp = getTransProfile();
    for (std::size_t i = 1; i <= n; i++) {
        for (std::size_t j = 0; j <= transMaxOrder; j++) {
            if (j <= transMaxOrder_m - i) {
                temp[j] = temp[j + 1] * static_cast<double>(j + 1);
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

double MultipoleT::getFnDerivX(const std::size_t& n, const double& x, const double& s) {
    if (n == 0) {
        return getTransDeriv(1, x) * getFringeDeriv(0, s);
    }
    double deriv    = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * implementation_->getFn(n, x - 2. * stepSize, s);
    deriv += -8. * implementation_->getFn(n, x - stepSize, s);
    deriv += 8. * implementation_->getFn(n, x + stepSize, s);
    deriv += -1. * implementation_->getFn(n, x + 2. * stepSize, s);
    deriv /= 12 * stepSize;
    return deriv;
}

double MultipoleT::getFnDerivS(const std::size_t& n, const double& x, const double& s) {
    if (n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(1, s);
    }
    double deriv    = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * implementation_->getFn(n, x, s - 2. * stepSize);
    deriv += -8. * implementation_->getFn(n, x, s - stepSize);
    deriv += 8. * implementation_->getFn(n, x, s + stepSize);
    deriv += -1. * implementation_->getFn(n, x, s + 2. * stepSize);
    deriv /= 12 * stepSize;
    return deriv;
}

void MultipoleT::finalise() {
    RefPartBunch_m = nullptr;
}

bool MultipoleT::apply(const size_t &i, const double &t,Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

void MultipoleT::setElementLength(double length) {
    // Base class first
    Component::setElementLength(length);
    // Then me
    length_m = length;
    implementation_->initialise();
}

void MultipoleT::setBendAngle(double angle, bool variableRadius) {
    // Record information
    bendAngle_m = angle;
    variableRadius_m = variableRadius;
    chooseImplementation();
}

void MultipoleT::chooseImplementation() {
    if(bendAngle_m == 0.0) {
        implementation_ = std::make_unique<MultipoleTStraight>(this);
    } else if(variableRadius_m) {
        implementation_ = std::make_unique<MultipoleTCurvedVarRadius>(this);
    } else {
        implementation_ = std::make_unique<MultipoleTCurvedConstRadius>(this);
    }
    implementation_->initialise();
}

void MultipoleT::setAperture(const double &vertAp, const double &horizAp) {
    verticalApert_m = vertAp;
    horizontalApert_m = horizAp;
}

void MultipoleT::setBoundingBoxLength(double boundingBoxLength) {
    boundingBoxLength_m = boundingBoxLength;
}

void MultipoleT::setTransProfile(const std::vector<double>& profile) {
    transProfile_m = profile;
    if(transProfile_m.empty()) {
        transProfile_m = {0.0};
    }
    transMaxOrder_m = transProfile_m.size() - 1;
}

void MultipoleT::setMaxOrder(size_t orderZ, size_t orderX) {
    maxFOrder_m = orderZ;
    maxXOrder_m = orderX;
    implementation_->setMaxOrder(maxFOrder_m, maxXOrder_m);
}

void MultipoleT::setRotation(double rot) {
    rotation_m = rot;
}

void MultipoleT::setEntranceAngle(double entranceAngle) {
    entranceAngle_m = entranceAngle;
}

void MultipoleT::setEntryOffset(double offset){
    entryOffset_m = offset;
}

bool MultipoleT::bends() const {
    return transProfile_m[0] != 0 || bendAngle_m != 0.0;
}

void MultipoleT::initialise(PartBunchBase<double, 3>* bunch,
                            double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
    implementation_->initialise();
}

BGeometryBase& MultipoleT::getGeometry() {
    return implementation_->getGeometry();
}

const BGeometryBase& MultipoleT::getGeometry() const {
    return implementation_->getGeometry();
}

Vector_t MultipoleT::localCartesianToOpalCartesian(const Vector_t& r) {
    return implementation_->localCartesianToOpalCartesian(r);
}

double MultipoleT::localCartesianRotation() {
    return implementation_->localCartesianRotation();
}
