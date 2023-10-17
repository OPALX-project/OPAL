//
// Class ScalingFFAMagnet
//   Defines the abstract interface for a sector FFA magnet
//   with radially scaling fringe fields.
//
// Copyright (c) 2017 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "AbsBeamline/ScalingFFAMagnet.h"

#include "AbsBeamline/BeamlineVisitor.h"

ScalingFFAMagnet::ScalingFFAMagnet(const std::string& name):
    Component(name),
    planarArcGeometry_m(1., 1.),
    dummy(),
    endField_m(nullptr) {
}

ScalingFFAMagnet::ScalingFFAMagnet(const ScalingFFAMagnet& right):
    Component(right),
    planarArcGeometry_m(right.planarArcGeometry_m),
    dummy(), maxOrder_m(right.maxOrder_m), tanDelta_m(right.tanDelta_m),
    k_m(right.k_m), Bz_m(right.Bz_m), r0_m(right.r0_m),
    rMin_m(right.rMin_m), rMax_m(right.rMax_m), phiStart_m(right.phiStart_m),
    phiEnd_m(right.phiEnd_m), azimuthalExtent_m(right.azimuthalExtent_m),
    verticalExtent_m(right.verticalExtent_m), centre_m(right.centre_m),
    endField_m(nullptr), endFieldName_m(right.endFieldName_m),
    dfCoefficients_m(right.dfCoefficients_m)
{
    endField_m = right.endField_m->clone();
    RefPartBunch_m = right.RefPartBunch_m;
    Bz_m = right.Bz_m;
    r0_m = right.r0_m;
    r0Sign_m = right.r0Sign_m;
}

ScalingFFAMagnet::~ScalingFFAMagnet() {
    delete endField_m;
}

ScalingFFAMagnet* ScalingFFAMagnet::clone() const {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

EMField& ScalingFFAMagnet::getField() {
    return dummy;
}

const EMField& ScalingFFAMagnet::getField() const {
    return dummy;
}

void ScalingFFAMagnet::initialise() {
    calculateDfCoefficients();
}

void ScalingFFAMagnet::initialise(PartBunchBase<double, 3>* bunch, double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
    initialise();
}

void ScalingFFAMagnet::finalise() {
    RefPartBunch_m = nullptr;
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

bool ScalingFFAMagnet::getFieldValue(const Vector_t& R, Vector_t& B) const {
    double x = r0Sign_m * (r0_m + R[0]);
    double r = std::sqrt(x * x + R[2] * R[2]);
    double phi = std::atan2(R[2], x); // angle between y-axis and position vector in anticlockwise direction
    Vector_t posCyl(r, R[1], phi);
    //Vector_t posCyl(r, pos[1], phi);
    Vector_t bCyl(0., 0., 0.); //br bz bphi
    bool outOfBounds = getFieldValueCylindrical(posCyl, bCyl);

    // this is cartesian coordinates
    B[1] += bCyl[1];
    B[0] += -r0Sign_m * (-bCyl[0] * std::cos(phi) + bCyl[2] * std::sin(phi));
    B[2] += bCyl[0] * std::sin(phi) + bCyl[2] * std::cos(phi);
    return outOfBounds;
}

bool ScalingFFAMagnet::getFieldValueCylindrical(const Vector_t& pos, Vector_t& B) const {
    double r = pos[0];
    double z = pos[1];
    double phi = pos[2];
    if (r < rMin_m || r > rMax_m) {
        return true;
    }
    if (z < -verticalExtent_m || z > verticalExtent_m) {
        return true;
    }
    double normRadius = r/std::abs(r0_m);
    double g = tanDelta_m*std::log(normRadius);
    double phiSpiral = phi - g - phiStart_m;
    if (phiSpiral < -azimuthalExtent_m || phiSpiral > azimuthalExtent_m) {
        return true;
    }

    double h = std::pow(normRadius, k_m)*Bz_m;
    std::vector<double> fringeDerivatives(maxOrder_m+1, 0.);
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = endField_m->function(phiSpiral, i); // d^i_phi f
    }

    double zOverR = r0Sign_m * z/r;
    std::vector<double> zOverRVec(dfCoefficients_m.size()+1); // zOverR^n
    zOverRVec[0] = 1.0;
    for (size_t n = 1; n < zOverRVec.size(); ++n) {
        zOverRVec[n] = zOverRVec[n-1] * zOverR;
    }
    for (size_t n = 0; n < dfCoefficients_m.size(); n += 2) {
        double f2n = 0;
        Vector_t deltaB;
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            f2n += dfCoefficients_m[n][i]*fringeDerivatives[i];
        }
        deltaB[1] = f2n * h * zOverRVec[n]; // Bz = sum(f_2n * h * (z/r)^2n
        if (maxOrder_m >= n+1) {
            double f2nplus1 = 0;
            for (size_t i = 0; i < dfCoefficients_m[n+1].size() && n+1 < dfCoefficients_m.size(); ++i) {
                f2nplus1 += dfCoefficients_m[n+1][i]*fringeDerivatives[i];
            }
            deltaB[0] = r0Sign_m * (f2n * (k_m-n) / (n+1) - tanDelta_m * f2nplus1) * h * zOverRVec[n+1]; // Br
            deltaB[2] = r0Sign_m * f2nplus1 * h * zOverRVec[n+1]; // Bphi = sum(f_2n+1 * h * (z/r)^2n+1
        }
        B += deltaB;
    }
    return false;
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
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] = -(k_m-n)*(k_m-n)/(n+1)*dfCoefficients_m[n][i]/(n+2);
        }
        for (size_t i = 0; i < dfCoefficients_m[n+1].size(); ++i) { // i indexes the derivative
            dfCoefficients_m[n+2][i] += 2*(k_m-n)*tanDelta_m*dfCoefficients_m[n+1][i]/(n+2);
            dfCoefficients_m[n+2][i+1] -= (1+tanDelta_m*tanDelta_m)*dfCoefficients_m[n+1][i]/(n+2);
        }
    }

}

void ScalingFFAMagnet::setEndField(endfieldmodel::EndFieldModel* endField) {
    if (endField_m != nullptr) {
        delete endField_m;
    }
    endField_m = endField;
}

// Note this is tested in OpalScalingFFAMagnetTest.*
void ScalingFFAMagnet::setupEndField() {
    if (endFieldName_m.empty()) { // no end field is defined
        return;
    }

    std::shared_ptr<endfieldmodel::EndFieldModel> efm
            = endfieldmodel::EndFieldModel::getEndFieldModel(endFieldName_m);

    endfieldmodel::EndFieldModel* newEFM = efm->clone();
    newEFM->rescale(1.0/std::abs(r0_m));
    setEndField(newEFM);
    newEFM->setMaximumDerivative(maxOrder_m+2);

    double defaultExtent = (newEFM->getEndLength()*4.+
                            newEFM->getCentreLength());
    if (phiStart_m < 0.0) {
        setPhiStart(defaultExtent/2.0);
    } else {
        setPhiStart(getPhiStart()+newEFM->getCentreLength()*0.5);
    }
    if (phiEnd_m < 0.0) {
        setPhiEnd(defaultExtent);
    }
    if (azimuthalExtent_m < 0.0) {
        setAzimuthalExtent(newEFM->getEndLength()*5.+
                           newEFM->getCentreLength()/2.0);
    }
    planarArcGeometry_m.setElementLength(std::abs(r0_m)*phiEnd_m); // length = phi r
    planarArcGeometry_m.setCurvature(1./r0_m);
}
