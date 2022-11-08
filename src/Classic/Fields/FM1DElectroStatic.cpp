#include "Fields/FM1DElectroStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include "gsl/gsl_fft_real.h"

#include <fstream>
#include <ios>

_FM1DElectroStatic::_FM1DElectroStatic(const std::string& filename):
    _Fieldmap(filename) {

    Type = T1DElectroStatic;

    std::ifstream fieldFile(Filename_m.c_str());
    if (fieldFile.good()) {

        bool parsingPassed = readFileHeader(fieldFile);
        parsingPassed = checkFileData(fieldFile, parsingPassed);
        fieldFile.close();

        if (!parsingPassed) {
            disableFieldmapWarning();
            zEnd_m = zBegin_m - 1.0e-3;
        } else
            convertHeaderData();

        length_m = (2.0 * numberOfGridPoints_m * (zEnd_m - zBegin_m)
                    / (numberOfGridPoints_m - 1));
    } else {
        noFieldmapWarning();
        zBegin_m = 0.0;
        zEnd_m = -1.0e-3;
    }
}

_FM1DElectroStatic::~_FM1DElectroStatic() {
    freeMap();
}

FM1DElectroStatic _FM1DElectroStatic::create(const std::string& filename)
{
    return FM1DElectroStatic(new _FM1DElectroStatic(filename));
}

void _FM1DElectroStatic::readMap() {

    if (fourierCoefs_m.empty()) {

        std::ifstream fieldFile(Filename_m.c_str());
        stripFileHeader(fieldFile);

        double *fieldData = new double[2 * numberOfGridPoints_m - 1];
        double maxEz = readFileData(fieldFile, fieldData);
        fieldFile.close();
        computeFourierCoefficients(maxEz, fieldData);
        delete [] fieldData;

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")
                << endl);
    }
}

void _FM1DElectroStatic::freeMap() {

    if (!fourierCoefs_m.empty()) {
        fourierCoefs_m.clear();
    }
}

bool _FM1DElectroStatic::getFieldstrength(const Vector_t &R, Vector_t &E,
                                         Vector_t &B) const {

    std::vector<double> fieldComponents;
    computeFieldOnAxis(R(2) - zBegin_m, fieldComponents);
    computeFieldOffAxis(R, E, B, fieldComponents);

    return false;
}

bool _FM1DElectroStatic::getFieldDerivative(const Vector_t &R,
                                           Vector_t &E,
                                           Vector_t &/*B*/,
                                           const DiffDirection &/*dir*/) const {

    double kz = Physics::two_pi * (R(2) - zBegin_m) / length_m + Physics::pi;
    double eZPrime = 0.0;

    int coefIndex = 1;
    for (int n = 1; n < accuracy_m; ++ n) {

        eZPrime += n * Physics::two_pi / length_m
            * (-fourierCoefs_m.at(coefIndex) * sin(kz * n)
               - fourierCoefs_m.at(coefIndex + 1) * cos(kz * n));
        coefIndex += 2;

    }

    E(2) +=  eZPrime;

    return false;
}

void _FM1DElectroStatic::getFieldDimensions(double &zBegin, double &zEnd) const {
    zBegin = zBegin_m;
    zEnd = zEnd_m;
}
void _FM1DElectroStatic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/,
                                           double &/*yIni*/, double &/*yFinal*/,
                                           double &/*zIni*/, double &/*zFinal*/) const {}

void _FM1DElectroStatic::swap()
{ }

void _FM1DElectroStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m
           << " (1D electrostatic); zini= "
           << zBegin_m << " m; zfinal= "
           << zEnd_m << " m;" << endl;
}

double _FM1DElectroStatic::getFrequency() const {
    return 0.0;
}

void _FM1DElectroStatic::setFrequency(double /*freq*/)
{ }

bool _FM1DElectroStatic::checkFileData(std::ifstream &fieldFile,
                                      bool parsingPassed) {

    double tempDouble;
    for (int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex)
        parsingPassed = parsingPassed
            && interpretLine<double>(fieldFile, tempDouble);

    return parsingPassed && interpreteEOF(fieldFile);

}

void _FM1DElectroStatic::computeFieldOffAxis(const Vector_t &R,
                                            Vector_t &E,
                                            Vector_t &/*B*/,
                                            std::vector<double> fieldComponents) const {

    double radiusSq = pow(R(0), 2.0) + pow(R(1), 2.0);
    double transverseEFactor = -fieldComponents.at(1) / 2.0
        + radiusSq * fieldComponents.at(3) / 16.0;

    E(0) += R(0) * transverseEFactor;
    E(1) += R(1) * transverseEFactor;
    E(2) += fieldComponents.at(0) - fieldComponents.at(2) * radiusSq / 4.0;

}

void _FM1DElectroStatic::computeFieldOnAxis(double z,
                                           std::vector<double> &fieldComponents) const {

    double kz = Physics::two_pi * z / length_m + Physics::pi;
    fieldComponents.push_back(fourierCoefs_m.at(0));
    fieldComponents.push_back(0.0);
    fieldComponents.push_back(0.0);
    fieldComponents.push_back(0.0);

    int coefIndex = 1;
    for (int n = 1; n < accuracy_m; ++ n) {

        double kn = n * Physics::two_pi / length_m;
        double coskzn = cos(kz * n);
        double sinkzn = sin(kz * n);

        fieldComponents.at(0) += fourierCoefs_m.at(coefIndex) * coskzn
            - fourierCoefs_m.at(coefIndex + 1) * sinkzn;

        fieldComponents.at(1) += kn * (-fourierCoefs_m.at(coefIndex) * sinkzn
                                       - fourierCoefs_m.at(coefIndex + 1) * coskzn);

        double derivCoeff = pow(kn, 2.0);
        fieldComponents.at(2) += derivCoeff * (-fourierCoefs_m.at(coefIndex) * coskzn
                                               + fourierCoefs_m.at(coefIndex + 1) * sinkzn);
        derivCoeff *= kn;
        fieldComponents.at(3) += derivCoeff * (fourierCoefs_m.at(coefIndex) * sinkzn
                                               + fourierCoefs_m.at(coefIndex + 1) * coskzn);

        coefIndex += 2;
    }
}

void _FM1DElectroStatic::computeFourierCoefficients(double maxEz,
                                                   double fieldData[]) {
    const unsigned int totalSize = 2 * numberOfGridPoints_m - 1;
    gsl_fft_real_wavetable *waveTable = gsl_fft_real_wavetable_alloc(totalSize);
    gsl_fft_real_workspace *workSpace = gsl_fft_real_workspace_alloc(totalSize);

    gsl_fft_real_transform(fieldData, 1, totalSize, waveTable, workSpace);

    /*
     * Normalize the Fourier coefficients such that the max field
     * value is 1 V/m.
     */

    fourierCoefs_m.push_back(fieldData[0] / (totalSize * maxEz * Units::Vpm2MVpm));
    for (int coefIndex = 1; coefIndex < 2 * accuracy_m - 1; ++ coefIndex)
        fourierCoefs_m.push_back(2.0 * fieldData[coefIndex] / (totalSize * maxEz * Units::Vpm2MVpm));

    gsl_fft_real_workspace_free(workSpace);
    gsl_fft_real_wavetable_free(waveTable);

}

void _FM1DElectroStatic::convertHeaderData() {

    // Convert to m.
    rBegin_m *= Units::cm2m;
    rEnd_m *= Units::cm2m;
    zBegin_m *= Units::cm2m;
    zEnd_m *= Units::cm2m;
}

double _FM1DElectroStatic::readFileData(std::ifstream &fieldFile,
                                       double fieldData[]) {

    double maxEz = 0.0;
    for (int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex) {
        interpretLine<double>(fieldFile,
                               fieldData[numberOfGridPoints_m - 1 + dataIndex]);
        if (std::abs(fieldData[numberOfGridPoints_m + dataIndex]) > maxEz)
            maxEz = std::abs(fieldData[numberOfGridPoints_m + dataIndex]);

        /*
         * Mirror the field map about minimum z value to ensure that
         * it is periodic.
         */
        if (dataIndex != 0)
            fieldData[numberOfGridPoints_m - 1 - dataIndex]
                = fieldData[numberOfGridPoints_m + dataIndex];
    }

    if (!normalize_m)
        maxEz = 1.0;

    return maxEz;
}

bool _FM1DElectroStatic::readFileHeader(std::ifstream &fieldFile) {

    std::string tempString;
    int tempInt;

    bool parsingPassed = true;
    try {
        parsingPassed = interpretLine<std::string, int>(fieldFile,
                                                         tempString,
                                                         accuracy_m);
    } catch (GeneralClassicException &e) {
        parsingPassed = interpretLine<std::string, int, std::string>(fieldFile,
                                                                      tempString,
                                                                      accuracy_m,
                                                                      tempString);

        tempString = Util::toUpper(tempString);
        if (tempString != "TRUE" &&
            tempString != "FALSE")
            throw GeneralClassicException("_FM1DElectroStatic::readFileHeader",
                                          "The third string on the first line of 1D field "
                                          "maps has to be either TRUE or FALSE");

        normalize_m = (tempString == "TRUE");
    }

    parsingPassed = parsingPassed &&
        interpretLine<double, double, int>(fieldFile,
                                            zBegin_m,
                                            zEnd_m,
                                            numberOfGridPoints_m);
    parsingPassed = parsingPassed &&
        interpretLine<double, double, int>(fieldFile, rBegin_m,
                                            rEnd_m, tempInt);

    ++ numberOfGridPoints_m;

    if (accuracy_m > numberOfGridPoints_m)
        accuracy_m = numberOfGridPoints_m;

    return parsingPassed;
}

void _FM1DElectroStatic::stripFileHeader(std::ifstream &fieldFile) {

    std::string tempString;

    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
}
