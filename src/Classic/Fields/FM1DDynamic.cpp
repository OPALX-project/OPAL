#include "Fields/FM1DDynamic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include "gsl/gsl_fft_real.h"

#include <fstream>
#include <iostream>

FM1DDynamic::FM1DDynamic(std::string aFilename):
    Fieldmap(aFilename) {

    Type = T1DDynamic;

    std::ifstream fieldFile(Filename_m.c_str());
    if(fieldFile.good()) {

        bool parsingPassed = readFileHeader(fieldFile);
        parsingPassed = checkFileData(fieldFile, parsingPassed);
        fieldFile.close();

        if(!parsingPassed) {
            disableFieldmapWarning();
            zEnd_m = zBegin_m - 1.0e-3;
        } else
            convertHeaderData();

        length_m = (2.0 * numberOfGridPoints_m * (zEnd_m - zBegin_m) /
                    (numberOfGridPoints_m - 1));

    } else {
        noFieldmapWarning();
        zBegin_m = 0.0;
        zEnd_m = -1e-3;
    }
}

FM1DDynamic::~FM1DDynamic() {
    freeMap();
}

void FM1DDynamic::readMap() {

    if(fourierCoefs_m.empty()) {

        std::ifstream fieldFile(Filename_m.c_str());
        stripFileHeader(fieldFile);

        double *fieldData = new double[2 * numberOfGridPoints_m - 1];
        double maxEz = readFileData(fieldFile, fieldData);
        fieldFile.close();
        computeFourierCoefficients(maxEz, fieldData);
        delete [] fieldData;

        INFOMSG(typeset_msg("Read in field map '" + Filename_m + "'", "info") << endl);
    }
}

void FM1DDynamic::freeMap() {

    if(!fourierCoefs_m.empty()) {
        fourierCoefs_m.clear();

        INFOMSG(typeset_msg("Freed field map '" + Filename_m  + "'", "info")
                << endl);
    }
}

bool FM1DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E,
                                   Vector_t &B) const {

    std::vector<double> fieldComponents;
    computeFieldOnAxis(R(2), fieldComponents);
    computeFieldOffAxis(R, E, B, fieldComponents);

    return false;
}

bool FM1DDynamic::getFieldDerivative(const Vector_t &R,
                                     Vector_t &E,
                                     Vector_t &/*B*/,
                                     const DiffDirection &/*dir*/) const {

    double kz = Physics::two_pi * R(2) / length_m + Physics::pi;
    double eZPrime = 0.0;

    int coefIndex = 1;
    for(int n = 1; n < accuracy_m; ++ n) {

        eZPrime += (n * Physics::two_pi / length_m
                    * (-fourierCoefs_m.at(coefIndex) * sin(kz * n)
                       - fourierCoefs_m.at(coefIndex + 1) * cos(kz * n)));
        coefIndex += 2;

    }

    E(2) +=  eZPrime;

    return false;
}

void FM1DDynamic::getFieldDimensions(double &zBegin, double &zEnd,
                                     double &rBegin, double &rEnd) const {
    zBegin = zBegin_m;
    zEnd = zEnd_m;
    rBegin = rBegin_m;
    rEnd = rEnd_m;
}

void FM1DDynamic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/,
                                     double &/*yIni*/, double &/*yFinal*/,
                                     double &/*zIni*/, double &/*zFinal*/) const {
}

void FM1DDynamic::swap()
{ }

void FM1DDynamic::getInfo(Inform *msg) {
    (*msg) << Filename_m
           << " (1D dynamic); zini= "
           << zBegin_m << " m; zfinal= "
           << zEnd_m << " m;" << endl;
}

double FM1DDynamic::getFrequency() const {
    return frequency_m;
}

void FM1DDynamic::setFrequency(double frequency) {
    frequency_m = frequency;
}

void FM1DDynamic::getOnaxisEz(std::vector<std::pair<double, double> > &eZ) {

    eZ.resize(numberOfGridPoints_m);
    std::ifstream fieldFile(Filename_m.c_str());
    stripFileHeader(fieldFile);
    double maxEz = readFileData(fieldFile, eZ);
    fieldFile.close();
    scaleField(maxEz, eZ);

}

bool FM1DDynamic::checkFileData(std::ifstream &fieldFile, bool parsingPassed) {

    double tempDouble;
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex)
        parsingPassed = parsingPassed
            && interpreteLine<double>(fieldFile, tempDouble);

    return parsingPassed && interpreteEOF(fieldFile);

}

void FM1DDynamic::computeFieldOffAxis(const Vector_t &R,
                                      Vector_t &E,
                                      Vector_t &B,
                                      std::vector<double> fieldComponents) const {

    double radiusSq = pow(R(0), 2.0) + pow(R(1), 2.0);
    double transverseEFactor = (fieldComponents.at(1)
                                * (0.5 - radiusSq * twoPiOverLambdaSq_m / 16.0)
                                - radiusSq * fieldComponents.at(3) / 16.0);
    double transverseBFactor = ((fieldComponents.at(0)
                                 * (0.5 - radiusSq * twoPiOverLambdaSq_m / 16.0)
                                 - radiusSq * fieldComponents.at(2) / 16.0)
                                * twoPiOverLambdaSq_m / frequency_m);

    E(0) += - R(0) * transverseEFactor;
    E(1) += - R(1) * transverseEFactor;
    E(2) += (fieldComponents.at(0)
             * (1.0 - radiusSq * twoPiOverLambdaSq_m / 4.0)
             - radiusSq * fieldComponents.at(2) / 4.0);

    B(0) += - R(1) * transverseBFactor;
    B(1) += R(0) * transverseBFactor;

}

void FM1DDynamic::computeFieldOnAxis(double z,
                                     std::vector<double> &fieldComponents) const {

    double kz = Physics::two_pi * z / length_m + Physics::pi;
    fieldComponents.push_back(fourierCoefs_m.at(0));
    fieldComponents.push_back(0.0);
    fieldComponents.push_back(0.0);
    fieldComponents.push_back(0.0);

    int coefIndex = 1;
    for(int n = 1; n < accuracy_m; ++ n) {

        double kn = n * Physics::two_pi / length_m;
        double coskzn = cos(kz * n);
        double sinkzn = sin(kz * n);

        fieldComponents.at(0) += (fourierCoefs_m.at(coefIndex) * coskzn
                                  - fourierCoefs_m.at(coefIndex + 1) * sinkzn);

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

void FM1DDynamic::computeFourierCoefficients(double maxEz, double fieldData[]) {
    const unsigned int totalSize = 2 * numberOfGridPoints_m - 1;
    gsl_fft_real_wavetable *waveTable = gsl_fft_real_wavetable_alloc(totalSize);
    gsl_fft_real_workspace *workSpace = gsl_fft_real_workspace_alloc(totalSize);
    gsl_fft_real_transform(fieldData, 1, totalSize, waveTable, workSpace);

    /*
     * Normalize the Fourier coefficients such that the max field value
     * is 1 V/m.
     */
    maxEz /= 1.0e6;

    fourierCoefs_m.push_back(fieldData[0] / (totalSize * maxEz));
    for(int coefIndex = 1; coefIndex < 2 * accuracy_m - 1; ++ coefIndex)
        fourierCoefs_m.push_back(2.0 * fieldData[coefIndex] / (totalSize * maxEz));

    gsl_fft_real_workspace_free(workSpace);
    gsl_fft_real_wavetable_free(waveTable);

}

void FM1DDynamic::convertHeaderData() {

    // Convert to angular frequency in Hz.
    frequency_m *= Physics::two_pi * 1.0e6;

    // Convert to m.
    rBegin_m /= 100.0;
    rEnd_m /= 100.0;
    zBegin_m /= 100.0;
    zEnd_m /= 100.0;

    twoPiOverLambdaSq_m = pow(frequency_m / Physics::c, 2.0);
}

double FM1DDynamic::readFileData(std::ifstream &fieldFile, double fieldData[]) {

    double maxEz = 0.0;
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex) {
        interpreteLine<double>(fieldFile, fieldData[numberOfGridPoints_m
                                                    - 1 + dataIndex]);
        if(std::abs(fieldData[numberOfGridPoints_m + dataIndex]) > maxEz)
            maxEz = std::abs(fieldData[numberOfGridPoints_m + dataIndex]);

        /*
         * Mirror the field map about minimum z value to ensure that it
         * is periodic.
         */
        if(dataIndex != 0)
            fieldData[numberOfGridPoints_m - 1 - dataIndex]
                = fieldData[numberOfGridPoints_m + dataIndex];
    }

    if (!normalize_m)
        maxEz = 1.0;

    return maxEz;
}

double FM1DDynamic::readFileData(std::ifstream &fieldFile,
                                 std::vector<std::pair<double, double>> &eZ) {

    double maxEz = 0.0;
    double deltaZ = (zEnd_m - zBegin_m) / (numberOfGridPoints_m - 1);
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex) {
        eZ.at(dataIndex).first = deltaZ * dataIndex;
        interpreteLine<double>(fieldFile, eZ.at(dataIndex).second);
        if(std::abs(eZ.at(dataIndex).second) > maxEz)
            maxEz = std::abs(eZ.at(dataIndex).second);
    }

    if (!normalize_m)
        maxEz = 1.0;

    return maxEz;
}

bool FM1DDynamic::readFileHeader(std::ifstream &fieldFile) {

    std::string tempString;
    int tempInt;

    bool parsingPassed = true;
    try {
        parsingPassed = interpreteLine<std::string, int>(fieldFile,
                                                         tempString,
                                                         accuracy_m);
    } catch (GeneralClassicException &e) {
        parsingPassed = interpreteLine<std::string, int, std::string>(fieldFile,
                                                                      tempString,
                                                                      accuracy_m,
                                                                      tempString);

        tempString = Util::toUpper(tempString);
        if (tempString != "TRUE" &&
            tempString != "FALSE")
            throw GeneralClassicException("FM1DDynamic::FM1DDynamic",
                                          "The third string on the first line of 1D field "
                                          "maps has to be either TRUE or FALSE");

        normalize_m = (tempString == "TRUE");
    }

    parsingPassed = parsingPassed &&
        interpreteLine<double, double, int>(fieldFile,
                                            zBegin_m,
                                            zEnd_m,
                                            numberOfGridPoints_m);
    parsingPassed = parsingPassed &&
        interpreteLine<double>(fieldFile, frequency_m);
    parsingPassed = parsingPassed &&
        interpreteLine<double, double, int>(fieldFile,
                                            rBegin_m,
                                            rEnd_m,
                                            tempInt);

    ++ numberOfGridPoints_m;

    if(accuracy_m > numberOfGridPoints_m)
        accuracy_m = numberOfGridPoints_m;

    return parsingPassed;
}

void FM1DDynamic::scaleField(double maxEz, std::vector<std::pair<double, double>> &eZ) {
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex)
        eZ.at(dataIndex).second /= maxEz;
}

void FM1DDynamic::stripFileHeader(std::ifstream &fieldFile) {

    std::string tempString;

    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
}

