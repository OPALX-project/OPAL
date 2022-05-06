//
// Class Astra1DDynamic_fast
//
// This class provides a reader for Astra style field maps. It pre-computes the field
// on a lattice to increase the performance during simulation.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#include "Fields/Astra1DDynamic_fast.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <fstream>
#include <ios>

Astra1DDynamic_fast::Astra1DDynamic_fast(std::string aFilename):
    Astra1D_fast(aFilename)
{
    numHeaderLines_m = 3;

    onAxisField_m = nullptr;

    Type = TAstraDynamic;

    // open field map, parse it and disable element on error
    std::ifstream file(Filename_m.c_str());
    if(!file.good()) {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
        return;
    }

    bool parsing_passed = readFileHeader(file);

    parsing_passed = parsing_passed && determineNumSamplingPoints(file);
    file.close();

    if(!parsing_passed) {
        disableFieldmapWarning();
        zend_m = zbegin_m - 1e-3;
        throw GeneralClassicException("Astra1DDynamic_fast::Astra1DDynamic_fast",
                                      "An error occured when reading the fieldmap '" + Filename_m + "'");
    }

    // conversion from MHz to Hz and from frequency to angular frequency
    frequency_m *= Physics::two_pi * Units::MHz2Hz;
    xlrep_m = frequency_m / Physics::c;

    hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    length_m = 2.0 * num_gridpz_m * hz_m;
}

Astra1DDynamic_fast::~Astra1DDynamic_fast() {
    freeMap();
}

void Astra1DDynamic_fast::readMap() {
    if(onAxisField_m == nullptr) {
        std::ifstream file(Filename_m.c_str());

        onAxisField_m = new double[num_gridpz_m];
        zvals_m = new double[num_gridpz_m];

        int accuracy = stripFileHeader(file);
        double maxEz = readFieldData(file);
        file.close();

        normalizeFieldData(maxEz * Units::Vpm2MVpm);

        std::vector<double> zvals = getEvenlyDistributedSamplingPoints();
        std::vector<double> evenFieldSampling = interpolateFieldData(zvals);
        std::vector<double> fourierCoefs = computeFourierCoefficients(accuracy, evenFieldSampling);

        computeFieldDerivatives(fourierCoefs, accuracy);

        checkMap(accuracy,
                 length_m,
                 zvals,
                 fourierCoefs,
                 onAxisInterpolants_m[0],
                 onAxisAccel_m[0]);

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m + "'", "info") << endl);
    }
}

bool Astra1DDynamic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do fourier interpolation in z-direction
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    double ez, ezp, ezpp, ezppp;
    try {
        ez = gsl_spline_eval(onAxisInterpolants_m[0], R(2) - zbegin_m, onAxisAccel_m[0]);
        ezp = gsl_spline_eval(onAxisInterpolants_m[1], R(2) - zbegin_m, onAxisAccel_m[1]);
        ezpp = gsl_spline_eval(onAxisInterpolants_m[2], R(2) - zbegin_m, onAxisAccel_m[2]);
        ezppp = gsl_spline_eval(onAxisInterpolants_m[3], R(2) - zbegin_m, onAxisAccel_m[3]);
    } catch (OpalException const& e) {
        throw OpalException("Astra1DDynamic_fast::getFieldstrength",
                            "The requested interpolation point, " + std::to_string(R(2)) + " is out of range");
    }
    // expand the field off-axis
    const double f  = -(ezpp  + ez *  xlrep_m * xlrep_m) / 16.;
    const double fp = -(ezppp + ezp * xlrep_m * xlrep_m) / 16.;

    const double EfieldR = -(ezp / 2. + fp * RR2);
    const double BfieldT = (ez / 2. + f * RR2) * xlrep_m / Physics::c;

    E(0) +=  EfieldR * R(0);
    E(1) +=  EfieldR * R(1);
    E(2) +=  ez + 4. * f * RR2;
    B(0) += -BfieldT * R(1);
    B(1) +=  BfieldT * R(0);

    return false;
}

bool Astra1DDynamic_fast::getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    double ezp = gsl_spline_eval(onAxisInterpolants_m[1], R(2) - zbegin_m, onAxisAccel_m[1]);

    E(2) +=  ezp;

    return false;
}

void Astra1DDynamic_fast::getFieldDimensions(double &zBegin, double &zEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
}

void Astra1DDynamic_fast::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void Astra1DDynamic_fast::swap()
{ }

void Astra1DDynamic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double Astra1DDynamic_fast::getFrequency() const {
    return frequency_m;
}

void Astra1DDynamic_fast::setFrequency(double freq) {
    frequency_m = freq;
}

void Astra1DDynamic_fast::getOnaxisEz(std::vector<std::pair<double, double> > & F) {
    F.resize(num_gridpz_m);
    if(onAxisField_m == nullptr) {
        double Ez_max = 0.0;
        double tmpDouble;
        int tmpInt;
        std::string tmpString;

        std::ifstream in(Filename_m.c_str());
        interpretLine<std::string, int>(in, tmpString, tmpInt);
        interpretLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpretLine<double>(in, tmpDouble);
        interpretLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for(int i = 0; i < num_gridpz_m; ++ i) {
            F[i].first = hz_m * i;
            interpretLine<double>(in, F[i].second);
            if(std::abs(F[i].second) > Ez_max) {
                Ez_max = std::abs(F[i].second);
            }
        }
        in.close();

        for(int i = 0; i < num_gridpz_m; ++ i) {
            F[i].second /= Ez_max;
        }
    } else {
        for(int i = 0; i < num_gridpz_m; ++ i) {
            F[i].first = zvals_m[i];
            F[i].second = onAxisField_m[i] / 1e6;
        }
    }
}

bool Astra1DDynamic_fast::readFileHeader(std::ifstream &file) {
    std::string tmpString;
    int tmpInt;
    bool passed = true;

    try {
        passed = interpretLine<std::string, int>(file, tmpString, tmpInt);
    } catch (GeneralClassicException &e) {
        passed = interpretLine<std::string, int, std::string>(file, tmpString, tmpInt, tmpString);

        tmpString = Util::toUpper(tmpString);
        if (tmpString != "TRUE" &&
            tmpString != "FALSE")
            throw GeneralClassicException("Astra1DDynamic_fast::readFileHeader",
                                          "The third string on the first line of 1D field "
                                          "maps has to be either TRUE or FALSE");

        normalize_m = (tmpString == "TRUE");
    }

    passed = passed && interpretLine<double>(file, frequency_m);

    return passed;
}

int Astra1DDynamic_fast::stripFileHeader(std::ifstream &file) {
    std::string tmpString;
    double tmpDouble;
    int accuracy;

    try {
        interpretLine<std::string, int>(file, tmpString, accuracy);
    } catch (GeneralClassicException &e) {
        interpretLine<std::string, int, std::string>(file, tmpString, accuracy, tmpString);
    }

    interpretLine<double>(file, tmpDouble);

    return accuracy;
}
