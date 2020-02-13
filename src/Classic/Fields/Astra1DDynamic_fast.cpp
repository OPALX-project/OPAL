#include "Fields/Astra1DDynamic_fast.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

Astra1DDynamic_fast::Astra1DDynamic_fast(std::string aFilename):
    Astra1D_fast(aFilename)
{
    numHeaderLines_m = 3;

    onAxisField_m = NULL;

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
    frequency_m *= Physics::two_pi * 1e6;
    xlrep_m = frequency_m / Physics::c;

    hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    length_m = 2.0 * num_gridpz_m * hz_m;
}

Astra1DDynamic_fast::~Astra1DDynamic_fast() {
    freeMap();
}

void Astra1DDynamic_fast::readMap() {
    if(onAxisField_m == NULL) {
        std::ifstream file(Filename_m.c_str());

        onAxisField_m = new double[num_gridpz_m];
        zvals_m = new double[num_gridpz_m];

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

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m + "'", "info") << endl);
    }
}

bool Astra1DDynamic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do fourier interpolation in z-direction
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    double ez = gsl_spline_eval(onAxisInterpolants_m[0], R(2), onAxisAccel_m[0]);
    double ezp = gsl_spline_eval(onAxisInterpolants_m[1], R(2), onAxisAccel_m[1]);
    double ezpp = gsl_spline_eval(onAxisInterpolants_m[2], R(2), onAxisAccel_m[2]);
    double ezppp = gsl_spline_eval(onAxisInterpolants_m[3], R(2), onAxisAccel_m[3]);

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
    double ezp = gsl_spline_eval(onAxisInterpolants_m[1], R(2), onAxisAccel_m[1]);

    E(2) +=  ezp;

    return false;
}

void Astra1DDynamic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &/*rBegin*/, double &/*rEnd*/) const {
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
    if(onAxisField_m == NULL) {
        double Ez_max = 0.0;
        double tmpDouble;
        int tmpInt;
        std::string tmpString;

        std::ifstream in(Filename_m.c_str());
        interpreteLine<std::string, int>(in, tmpString, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double>(in, tmpDouble);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for(int i = 0; i < num_gridpz_m; ++ i) {
            F[i].first = hz_m * i;
            interpreteLine<double>(in, F[i].second);
            if(fabs(F[i].second) > Ez_max) {
                Ez_max = fabs(F[i].second);
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

    passed = passed && interpreteLine<double>(file, frequency_m);

    return passed;
}

int Astra1DDynamic_fast::stripFileHeader(std::ifstream &file) {
    std::string tmpString;
    double tmpDouble;
    int accuracy;

    try {
        interpreteLine<std::string, int>(file, tmpString, accuracy);
    } catch (GeneralClassicException &e) {
        interpreteLine<std::string, int, std::string>(file, tmpString, accuracy, tmpString);
    }

    interpreteLine<double>(file, tmpDouble);

    return accuracy;
}

