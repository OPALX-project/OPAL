// ------------------------------------------------------------------------
// $RCSfile: PartData.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PartData
//   PartData represents a set of reference values for use in algorithms.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2003/12/02 23:04:59 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartData.h"
#include "Utilities/LogicalError.h"
#include <cmath>


// Class PartData
// ------------------------------------------------------------------------

PartData::PartData(double q, double m, double momentum) {
    charge = q;
    mass = m;
    setP(momentum);
}


PartData::PartData() {
    charge = 1.0;
    mass = 0.0;
    beta = 1.0;
    gamma = 1.0e10;
}


void PartData::setP(double p) {
    if(mass == 0.0) {
        throw LogicalError("PartData::setP()",
                           "Particle mass must not be zero.");
    }

    if(p == 0.0) {
        throw LogicalError("PartData::setP()",
                           "Particle momentum must not be zero.");
    }

    double e = std::sqrt(p * p + mass * mass);
    beta = p / e;
    gamma = e / mass;
}


void PartData::setE(double energy) {
    if(energy <= mass) {
        throw LogicalError("PartData::setE()", "Energy should be > mass.");
    }

    gamma = energy / mass;
    //beta = std::sqrt(energy*energy - mass*mass) / energy;
    double ginv = 1.0 / gamma;
    beta = std::sqrt((1.0 - ginv) * (1.0 + ginv));
}


void PartData::setBeta(double v) {
    if(v >= 1.0) {
        throw LogicalError("PartData::setBeta()", "Beta should be < 1.");
    }

    beta = v;
    gamma = 1.0 / std::sqrt(1.0 - beta * beta);
}


void PartData::setGamma(double v) {
    if(v <= 1.0) {
        throw LogicalError("PartData::setGamma()", "Gamma should be > 1.");
    }

    gamma = v;
    beta = std::sqrt(gamma * gamma - 1.0) / gamma;
}