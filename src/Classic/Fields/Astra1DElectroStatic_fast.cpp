#include "Fields/Astra1DElectroStatic_fast.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

Astra1DElectroStatic_fast::Astra1DElectroStatic_fast(std::string aFilename):
    Astra1D_fast(aFilename) {
    numHeaderLines_m = 2;

    onAxisField_m = NULL;

    Type = TAstraElectroStatic;

    // open field map, parse it and disable element on error
    std::ifstream file(Filename_m.c_str());
    if(!file.good()) {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }

    bool parsing_passed = readFileHeader(file);
    parsing_passed = parsing_passed && determineNumSamplingPoints(file);
    file.close();

    if(!parsing_passed && !file.eof()) {
        disableFieldmapWarning();
        zend_m = zbegin_m - 1e-3;
        throw GeneralClassicException("Astra1DElectroStatic_fast::Astra1DElectroStatic_fast",
                                      "An error occured when reading the fieldmap '" + Filename_m + "'");
    }
    hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    length_m = 2.0 * num_gridpz_m * hz_m;
}

Astra1DElectroStatic_fast::~Astra1DElectroStatic_fast() {
    freeMap();
}

void Astra1DElectroStatic_fast::readMap() {
    if(onAxisField_m == NULL) {
        // declare variables and allocate memory

        onAxisField_m = new double[num_gridpz_m];
        zvals_m = new double[num_gridpz_m];

        std::ifstream file(Filename_m.c_str());
        int accuracy = stripFileHeader(file);
        double maxEz = readFieldData(file);
        file.close();

        normalizeFieldData(maxEz * 1e-6);

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


        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

bool Astra1DElectroStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &/*B*/) const {
    // do fourier interpolation in z-direction
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    double ez = gsl_spline_eval(onAxisInterpolants_m[0], R(2), onAxisAccel_m[0]);
    double ezp = gsl_spline_eval(onAxisInterpolants_m[1], R(2), onAxisAccel_m[1]);
    double ezpp = gsl_spline_eval(onAxisInterpolants_m[2], R(2), onAxisAccel_m[2]);
    double ezppp = gsl_spline_eval(onAxisInterpolants_m[3], R(2), onAxisAccel_m[3]);

    // expand to off-axis
    const double EfieldR = -ezp / 2. + ezppp / 16. * RR2;

    E(0) += EfieldR * R(0);
    E(1) += EfieldR * R(1);
    E(2) += ez - ezpp * RR2 / 4.;
    return false;
}

bool Astra1DElectroStatic_fast::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void Astra1DElectroStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &/*rBegin*/, double &/*rEnd*/) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
}
void Astra1DElectroStatic_fast::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void Astra1DElectroStatic_fast::swap()
{ }

void Astra1DElectroStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D electrostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double Astra1DElectroStatic_fast::getFrequency() const {
    return 0.0;
}

void Astra1DElectroStatic_fast::setFrequency(double /*freq*/)
{ }

bool Astra1DElectroStatic_fast::readFileHeader(std::ifstream &file) {
    std::string tmpString;
    int tmpInt;

    bool passed;
    try {
        passed = interpreteLine<std::string, int>(file, tmpString, tmpInt);
    } catch (GeneralClassicException &e) {
        passed = interpreteLine<std::string, int, std::string>(file, tmpString, tmpInt, tmpString);

        tmpString = Util::toUpper(tmpString);
        if (tmpString != "TRUE" &&
            tmpString != "FALSE")
            throw GeneralClassicException("Astra1DDynamic_fast::readFileHeader",
                                          "The third string on the first line of 1D field "
                                          "maps has to be either TRUE or FALSE");

        normalize_m = (tmpString == "TRUE");
    }

    return passed;
}

int Astra1DElectroStatic_fast::stripFileHeader(std::ifstream &file) {
    std::string tmpString;
    int accuracy;

    try {
        interpreteLine<std::string, int>(file, tmpString, accuracy);
    } catch (GeneralClassicException &e) {
        interpreteLine<std::string, int, std::string>(file, tmpString, accuracy, tmpString);
    }

    return accuracy;
}

