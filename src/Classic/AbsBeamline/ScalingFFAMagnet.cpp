/*
 *  Copyright (c) 2017, Chris Rogers
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

#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"

ScalingFFAMagnet::ScalingFFAMagnet(const std::string &name)
        : Component(name),
         planarArcGeometry_m(1., 1.), dummy(), endField_m(NULL) {
    setElType(isDrift);
}

ScalingFFAMagnet::ScalingFFAMagnet(const ScalingFFAMagnet &right)
        : Component(right),
          planarArcGeometry_m(right.planarArcGeometry_m),
          dummy(), maxOrder_m(right.maxOrder_m), tanDelta_m(right.tanDelta_m),
          k_m(right.k_m), Bz_m(right.Bz_m), r0_m(right.r0_m),
          rMin_m(right.rMin_m), rMax_m(right.rMax_m), phiStart_m(right.phiStart_m),
          phiEnd_m(right.phiEnd_m), azimuthalExtent_m(right.azimuthalExtent_m),
          verticalExtent_m(right.verticalExtent_m), centre_m(right.centre_m),
          dfCoefficients_m(right.dfCoefficients_m) {
    if (endField_m != NULL) {
        delete endField_m;
    }
    endField_m = right.endField_m->clone();
    RefPartBunch_m = right.RefPartBunch_m;
    setElType(isDrift);
    Bz_m = right.Bz_m;
    r0_m = right.r0_m;
}

ScalingFFAMagnet::~ScalingFFAMagnet() {
    if (endField_m != NULL) {
        delete endField_m;
    }
}

ElementBase* ScalingFFAMagnet::clone() const {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

EMField &ScalingFFAMagnet::getField() {
    return dummy;
}

const EMField &ScalingFFAMagnet::getField() const {
    return dummy;
}

bool ScalingFFAMagnet::apply(const size_t &i, const double &t,
                    Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

void ScalingFFAMagnet::initialise() {
    calculateDfCoefficients();
    planarArcGeometry_m.setElementLength(r0_m*phiEnd_m); // length = phi r
    planarArcGeometry_m.setCurvature(1./r0_m);
}

void ScalingFFAMagnet::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    initialise();
}

void ScalingFFAMagnet::finalise() {
    RefPartBunch_m = NULL;
}

bool ScalingFFAMagnet::bends() const {
    return true;
}

BGeometryBase& ScalingFFAMagnet::getGeometry() {
    return planarArcGeometry_m;
}

const BGeometryBase& ScalingFFAMagnet::getGeometry() const {
    return planarArcGeometry_m;
}

void ScalingFFAMagnet::accept(BeamlineVisitor& visitor) const {
    visitor.visitScalingFFAMagnet(*this);
}


bool ScalingFFAMagnet::getFieldValue(const Vector_t &R, Vector_t &B) const {
    Vector_t pos = R - centre_m;
    double r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    double phi = -atan2(pos[0], pos[2]); // angle between y-axis and position vector in anticlockwise direction
    Vector_t posCyl(r, pos[1], phi);
    Vector_t bCyl(0., 0., 0.); //br bz bphi
    bool outOfBounds = getFieldValueCylindrical(posCyl, bCyl);
    // this is cartesian coordinates
    B[1] += bCyl[1];
    B[0] += -bCyl[2]*cos(phi) -bCyl[0]*sin(phi);
    B[2] += +bCyl[0]*cos(phi) -bCyl[2]*sin(phi);
    return outOfBounds;

}


bool ScalingFFAMagnet::getFieldValueCylindrical(const Vector_t &pos, Vector_t &B) const {

    double r = pos[0];
    double z = pos[1];
    double phi = pos[2];
    if (r < rMin_m || r > rMax_m) {
        return true;
    }

    double normRadius = r/r0_m;
    double g = tanDelta_m*log(normRadius);
    double phiSpiral = phi - g - phiStart_m;
    double h = pow(normRadius, k_m)*Bz_m;
    if (phiSpiral < -azimuthalExtent_m || phiSpiral > azimuthalExtent_m) {
        return true;
    }
    if (z < -verticalExtent_m || z > verticalExtent_m) {
        return true;
    }
    std::vector<double> fringeDerivatives(maxOrder_m+1, 0.);
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = endField_m->function(phiSpiral, i); // d^i_phi f
    }
    for (size_t n = 0; n < dfCoefficients_m.size(); n += 2) {
        double f2n = 0;
        Vector_t deltaB;
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            f2n += dfCoefficients_m[n][i]*fringeDerivatives[i];
        }
        deltaB[1] = f2n*h*pow(z/r, n); // Bz = sum(f_2n * h * (z/r)^2n
        if (maxOrder_m >= n+1) {
            double f2nplus1 = 0;
            for (size_t i = 0; i < dfCoefficients_m[n+1].size() && n+1 < dfCoefficients_m.size(); ++i) {
                f2nplus1 += dfCoefficients_m[n+1][i]*fringeDerivatives[i];
            }
            deltaB[0] = (f2n*(k_m-n)/(n+1) - tanDelta_m*f2nplus1)*h*pow(z/r, n+1); // Br
            deltaB[2] = f2nplus1*h*pow(z/r, n+1); // Bphi = sum(f_2n+1 * h * (z/r)^2n+1
        }
        B += deltaB;
    }
    return false;
}


bool ScalingFFAMagnet::apply(const Vector_t &R, const Vector_t &P,
                    const double &t, Vector_t &E, Vector_t &B) {
    return getFieldValue(R, B);
}

void ScalingFFAMagnet::calculateDfCoefficients() {
    dfCoefficients_m = std::vector<std::vector<double> >(maxOrder_m+1);
    dfCoefficients_m[0] = std::vector<double>(1, 1.); // f_0 = 1.*0th derivative
    for (size_t n = 0; n < maxOrder_m; n += 2) { // n indexes the power in z
        dfCoefficients_m[n+1] = std::vector<double>(dfCoefficients_m[n].size()+1, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+1][i+1] = dfCoefficients_m[n][i]/(n+1);
        }
        if (n+1 == maxOrder_m) {
            break;
        }
        dfCoefficients_m[n+2] = std::vector<double>(dfCoefficients_m[n].size()+2, 0);
        for(size_t i = 0; i < dfCoefficients_m[n].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] = -(k_m-n)*(k_m-n)/(n+1)*dfCoefficients_m[n][i]/(n+2);
        }
        for(size_t i = 0; i < dfCoefficients_m[n+1].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] += 2*(k_m-n)*tanDelta_m*dfCoefficients_m[n+1][i]/(n+2);
            dfCoefficients_m[n+2][i+1] -= (1+tanDelta_m*tanDelta_m)*dfCoefficients_m[n+1][i]/(n+2);
        }
    }

}

void ScalingFFAMagnet::setEndField(endfieldmodel::EndFieldModel* endField) {
    if (endField_m != NULL) {
        delete endField_m;
    }
    endField_m = endField;
}

