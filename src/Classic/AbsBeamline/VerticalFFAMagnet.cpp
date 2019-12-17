//
// Source file for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/EndFieldModel/EndFieldModel.h"

#include "AbsBeamline/VerticalFFAMagnet.h"

VerticalFFAMagnet::VerticalFFAMagnet(const std::string &name)
        : Component(name), straightGeometry_m(1.) {
    setElType(isDrift);
}

VerticalFFAMagnet::VerticalFFAMagnet(const VerticalFFAMagnet &right) 
    : Component(right),
    straightGeometry_m(right.straightGeometry_m),
    dummy(right.dummy),
    maxOrder_m(right.maxOrder_m),
    k_m(right.k_m),
    Bz_m(right.Bz_m),
    zNegExtent_m(right.zNegExtent_m),
    zPosExtent_m(right.zPosExtent_m),
    halfWidth_m(right.halfWidth_m),
    bbLength_m(right.bbLength_m),
    endField_m(right.endField_m->clone()),
    dfCoefficients_m(right.dfCoefficients_m) {
        RefPartBunch_m = right.RefPartBunch_m;
}


VerticalFFAMagnet::~VerticalFFAMagnet() {
}

ElementBase* VerticalFFAMagnet::clone() const {
    VerticalFFAMagnet* magnet = new VerticalFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

EMField &VerticalFFAMagnet::getField() {
    return dummy;
}

const EMField &VerticalFFAMagnet::getField() const {
    return dummy;
}

void VerticalFFAMagnet::initialise() {
    calculateDfCoefficients();
    straightGeometry_m.setElementLength(bbLength_m); // length = phi r
}

void VerticalFFAMagnet::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    initialise();
}

void VerticalFFAMagnet::finalise() {
    RefPartBunch_m = NULL;
}

BGeometryBase& VerticalFFAMagnet::getGeometry() {
    return straightGeometry_m;
}

const BGeometryBase& VerticalFFAMagnet::getGeometry() const {
    return straightGeometry_m;
}

void VerticalFFAMagnet::accept(BeamlineVisitor& visitor) const {
    visitor.visitVerticalFFAMagnet(*this);
}


bool VerticalFFAMagnet::getFieldValue(const Vector_t &R, Vector_t &B) const {
    if (abs(R[0]) > halfWidth_m ||
        R[2] < 0. || R[2] > bbLength_m ||
        R[1] < -zNegExtent_m || R[1] > zPosExtent_m) {
        return true;
    }
    std::vector<double> fringeDerivatives(maxOrder_m+2, 0.);
    double zRel = R[2]-bbLength_m/2.; // z relative to centre of magnet
    for (size_t i = 0; i < fringeDerivatives.size(); ++i) {
        fringeDerivatives[i] = endField_m->function(zRel, i); // d^i_phi f
    }

    std::vector<double> x_n(maxOrder_m+1); // x^n
    x_n[0] = 1.; // x^0
    for (size_t i = 1; i < x_n.size(); ++i) {
        x_n[i] = x_n[i-1]*R[0];
    }

    // note that the last element is always 0, because dfCoefficients_m is 
    // of size maxOrder_m+1. This leads to better Maxwellianness in testing.
    std::vector<double> f_n(maxOrder_m+2, 0.);
    std::vector<double> dz_f_n(maxOrder_m+1, 0.);
    for (size_t n = 0; n < dfCoefficients_m.size(); ++n) {
        const std::vector<double>& coefficients = dfCoefficients_m[n];
        for (size_t i = 0; i < coefficients.size(); ++i) {
            f_n[n] += coefficients[i]*fringeDerivatives[i];
            dz_f_n[n] += coefficients[i]*fringeDerivatives[i+1];
        }
    }
    double bref = Bz_m*exp(k_m*R[1]);
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 0.;
    for (size_t n = 0; n < x_n.size(); ++n) {
        B[0] += bref*f_n[n+1]*(n+1)/k_m*x_n[n];
        B[1] += bref*f_n[n]*x_n[n];
        B[2] += bref*dz_f_n[n]/k_m*x_n[n];
    }
    return false;
}

void VerticalFFAMagnet::calculateDfCoefficients() {
    dfCoefficients_m = std::vector< std::vector<double> >(maxOrder_m+1);
    dfCoefficients_m[0] = std::vector<double>(1, 1.);
    if (maxOrder_m > 0) {
        dfCoefficients_m[1] = std::vector<double>();
    }
    // n indexes like the polynomial order of the midplane expansion
    // e.g. Bz = exp(mz) f_n y^n
    // where y is distance from the midplane and z is height
    for (size_t n = 2; n < dfCoefficients_m.size(); n+=2) {
        const std::vector<double>& oldCoefficients = dfCoefficients_m[n-2];
        std::vector<double> coefficients(oldCoefficients.size()+2, 0);
        // j indexes the derivative of f_0
        for (size_t j = 0; j < oldCoefficients.size(); ++j) {
            coefficients[j] += -1./(n)/(n-1)*k_m*k_m*oldCoefficients[j];
            coefficients[j+2] += -1./(n)/(n-1)*oldCoefficients[j];
        }
        dfCoefficients_m[n] = coefficients;
    }
}

void VerticalFFAMagnet::setEndField(endfieldmodel::EndFieldModel* endField) {
    endField_m.reset(endField);
}
