#include "Fields/FM1DMagnetoStatic_fast.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/Util.h"
#include "Utilities/GeneralClassicException.h"

#include "gsl/gsl_fft_real.h"

#include <fstream>
#include <ios>

FM1DMagnetoStatic_fast::FM1DMagnetoStatic_fast(std::string aFilename):
    Fieldmap(aFilename) {

    Type = T1DMagnetoStatic;
    onAxisField_m = NULL;

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

        deltaZ_m = (zEnd_m - zBegin_m) / (numberOfGridPoints_m - 1);
        length_m = 2.0 * numberOfGridPoints_m * deltaZ_m;

    } else {
        noFieldmapWarning();
        zBegin_m = 0.0;
        zEnd_m = -1.0e-3;
    }
}

FM1DMagnetoStatic_fast::~FM1DMagnetoStatic_fast() {
    freeMap();
}

void FM1DMagnetoStatic_fast::readMap() {

    if(onAxisField_m == NULL) {

        std::ifstream fieldFile(Filename_m.c_str());
        stripFileHeader(fieldFile);

        onAxisField_m = new double[numberOfGridPoints_m];
        double maxBz = readFileData(fieldFile, onAxisField_m);
        fieldFile.close();

        std::vector<double> fourierCoefs
            = computeFourierCoefficients(onAxisField_m);
        normalizeField(maxBz, fourierCoefs);

        double *onAxisFieldP = new double[numberOfGridPoints_m];
        double *onAxisFieldPP = new double[numberOfGridPoints_m];
        double *onAxisFieldPPP = new double[numberOfGridPoints_m];
        computeFieldDerivatives(fourierCoefs, onAxisFieldP,
                                onAxisFieldPP, onAxisFieldPPP);
        computeInterpolationVectors(onAxisFieldP, onAxisFieldPP,
                                    onAxisFieldPPP);

        prepareForMapCheck(fourierCoefs);

        delete [] onAxisFieldP;
        delete [] onAxisFieldPP;
        delete [] onAxisFieldPPP;

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")
                << endl);
    }
}

void FM1DMagnetoStatic_fast::freeMap() {
    if(onAxisField_m != NULL) {
        delete [] onAxisField_m;
        onAxisField_m = NULL;

        gsl_spline_free(onAxisFieldInterpolants_m);
        gsl_spline_free(onAxisFieldPInterpolants_m);
        gsl_spline_free(onAxisFieldPPInterpolants_m);
        gsl_spline_free(onAxisFieldPPPInterpolants_m);
        gsl_interp_accel_free(onAxisFieldAccel_m);
        gsl_interp_accel_free(onAxisFieldPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPPAccel_m);

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info")
                << endl);
    }
}

bool FM1DMagnetoStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E,
                                              Vector_t &B) const {

    std::vector<double> fieldComponents;
    computeFieldOnAxis(R(2), fieldComponents);
    computeFieldOffAxis(R, E, B, fieldComponents);

    return false;

}

bool FM1DMagnetoStatic_fast::getFieldDerivative(const Vector_t &R,
                                                Vector_t &/*E*/,
                                                Vector_t &B,
                                                const DiffDirection &/*dir*/) const {

    B(2) += gsl_spline_eval(onAxisFieldPInterpolants_m, R(2),
                            onAxisFieldPAccel_m);

    return false;

}

void FM1DMagnetoStatic_fast::getFieldDimensions(double &zBegin, double &zEnd,
                                                double &rBegin, double &rEnd) const {
    zBegin = zBegin_m;
    zEnd = zEnd_m;
    rBegin = rBegin_m;
    rEnd = rEnd_m;
}

void FM1DMagnetoStatic_fast::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/,
                                                double &/*yIni*/, double &/*yFinal*/,
                                                double &/*zIni*/, double &/*zFinal*/) const {}

void FM1DMagnetoStatic_fast::swap()
{ }

void FM1DMagnetoStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m
           << " (1D magnetostatic); zini= "
           << zBegin_m << " m; zfinal= "
           << zEnd_m << " m;" << endl;
}

double FM1DMagnetoStatic_fast::getFrequency() const {
    return 0.0;
}

void FM1DMagnetoStatic_fast::setFrequency(double /*freq*/)
{ }

bool FM1DMagnetoStatic_fast::checkFileData(std::ifstream &fieldFile,
                                           bool parsingPassed) {

    double tempDouble;
    for (unsigned int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex)
        parsingPassed = parsingPassed
            && interpreteLine<double>(fieldFile, tempDouble);

    return parsingPassed && interpreteEOF(fieldFile);

}

void FM1DMagnetoStatic_fast::computeFieldDerivatives(std::vector<double> fourierCoefs,
                                                     double onAxisFieldP[],
                                                     double onAxisFieldPP[],
                                                     double onAxisFieldPPP[]) {

    for (unsigned int zStepIndex = 0; zStepIndex < numberOfGridPoints_m; ++ zStepIndex) {

        double z = deltaZ_m * zStepIndex;
        double kz = Physics::two_pi * z / length_m + Physics::pi;
        onAxisFieldP[zStepIndex] = 0.0;
        onAxisFieldPP[zStepIndex] = 0.0;
        onAxisFieldPPP[zStepIndex] = 0.0;

        int coefIndex = 1;
        for (unsigned int n = 1; n < accuracy_m; ++ n) {

            double kn = n * Physics::two_pi / length_m;
            double coskzn = cos(kz * n);
            double sinkzn = sin(kz * n);

            onAxisFieldP[zStepIndex] += kn * (-fourierCoefs.at(coefIndex) * sinkzn
                                              - fourierCoefs.at(coefIndex + 1) * coskzn);

            double derivCoeff = pow(kn, 2.0);
            onAxisFieldPP[zStepIndex] += derivCoeff * (-fourierCoefs.at(coefIndex) * coskzn
                                                       + fourierCoefs.at(coefIndex + 1) * sinkzn);
            derivCoeff *= kn;
            onAxisFieldPPP[zStepIndex] += derivCoeff * (fourierCoefs.at(coefIndex) * sinkzn
                                                        + fourierCoefs.at(coefIndex + 1) * coskzn);

            coefIndex += 2;
        }
    }
}

void FM1DMagnetoStatic_fast::computeFieldOffAxis(const Vector_t &R,
                                                 Vector_t &/*E*/,
                                                 Vector_t &B,
                                                 std::vector<double> fieldComponents) const {

    double radiusSq = pow(R(0), 2.0) + pow(R(1), 2.0);
    double transverseBFactor = -fieldComponents.at(1) / 2.0
        + radiusSq * fieldComponents.at(3) / 16.0;

    B(0) += R(0) * transverseBFactor;
    B(1) += R(1) * transverseBFactor;
    B(2) += fieldComponents.at(0) - fieldComponents.at(2) * radiusSq / 4.0;

}

void FM1DMagnetoStatic_fast::computeFieldOnAxis(double z,
                                                std::vector<double> &fieldComponents) const {

    fieldComponents.push_back(gsl_spline_eval(onAxisFieldInterpolants_m,
                                              z, onAxisFieldAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPInterpolants_m,
                                              z, onAxisFieldPAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPPInterpolants_m,
                                              z, onAxisFieldPPAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPPPInterpolants_m,
                                              z, onAxisFieldPPPAccel_m));
}

std::vector<double> FM1DMagnetoStatic_fast::computeFourierCoefficients(double fieldData[]) {

    const unsigned int totalSize = 2 * numberOfGridPoints_m - 1;
    gsl_fft_real_wavetable *waveTable = gsl_fft_real_wavetable_alloc(totalSize);
    gsl_fft_real_workspace *workSpace = gsl_fft_real_workspace_alloc(totalSize);

    // Reflect field about minimum z value to ensure that it is periodic.
    double *fieldDataReflected = new double[totalSize];
    for (unsigned int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex) {
        fieldDataReflected[numberOfGridPoints_m - 1 + dataIndex]
            = fieldData[dataIndex];
        if(dataIndex != 0)
            fieldDataReflected[numberOfGridPoints_m - 1 - dataIndex]
                = fieldData[dataIndex];
    }

    gsl_fft_real_transform(fieldDataReflected, 1, totalSize, waveTable, workSpace);

    std::vector<double> fourierCoefs;
    fourierCoefs.push_back(fieldDataReflected[0] / totalSize);
    for (unsigned int coefIndex = 1; coefIndex + 1 < 2 * accuracy_m; ++ coefIndex)
        fourierCoefs.push_back(2.0 * fieldDataReflected[coefIndex] / totalSize);

    delete [] fieldDataReflected;
    gsl_fft_real_workspace_free(workSpace);
    gsl_fft_real_wavetable_free(waveTable);

    return fourierCoefs;

}

void FM1DMagnetoStatic_fast::computeInterpolationVectors(double onAxisFieldP[],
                                                         double onAxisFieldPP[],
                                                         double onAxisFieldPPP[]) {

    onAxisFieldInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                                 numberOfGridPoints_m);
    onAxisFieldPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                                  numberOfGridPoints_m);
    onAxisFieldPPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                                   numberOfGridPoints_m);
    onAxisFieldPPPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                                    numberOfGridPoints_m);

    double *z = new double[numberOfGridPoints_m];
    for (unsigned int zStepIndex = 0; zStepIndex < numberOfGridPoints_m; ++ zStepIndex)
        z[zStepIndex] = deltaZ_m * zStepIndex;

    gsl_spline_init(onAxisFieldInterpolants_m, z,
                    onAxisField_m, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPInterpolants_m, z,
                    onAxisFieldP, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPPInterpolants_m, z,
                    onAxisFieldPP, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPPPInterpolants_m, z,
                    onAxisFieldPPP, numberOfGridPoints_m);

    onAxisFieldAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPPAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPPPAccel_m = gsl_interp_accel_alloc();

    delete [] z;
}

void FM1DMagnetoStatic_fast::convertHeaderData() {

    // Convert to m.
    rBegin_m /= 100.0;
    rEnd_m /= 100.0;
    zBegin_m /= 100.0;
    zEnd_m /= 100.0;
}

void FM1DMagnetoStatic_fast::normalizeField(double maxBz, std::vector<double> &fourierCoefs) {
    if (!normalize_m) return;

    for (unsigned int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex)
        onAxisField_m[dataIndex] /= maxBz;

    for (std::vector<double>::iterator fourierIterator = fourierCoefs.begin();
        fourierIterator < fourierCoefs.end(); ++ fourierIterator)
        *fourierIterator /= maxBz;

}

double FM1DMagnetoStatic_fast::readFileData(std::ifstream &fieldFile,
                                            double fieldData[]) {

    double maxBz = 0.0;
    for (unsigned int dataIndex = 0; dataIndex < numberOfGridPoints_m; ++ dataIndex) {
        interpreteLine<double>(fieldFile, fieldData[dataIndex]);
        if(std::abs(fieldData[dataIndex]) > maxBz)
            maxBz = std::abs(fieldData[dataIndex]);
    }

    return maxBz;
}

bool FM1DMagnetoStatic_fast::readFileHeader(std::ifstream &fieldFile) {

    std::string tempString;
    int tempInt;

    bool parsingPassed = true;

    try {
        parsingPassed = interpreteLine<std::string, unsigned int>(fieldFile,
                                                                  tempString,
                                                                  accuracy_m);
    } catch (GeneralClassicException &e) {
        parsingPassed = interpreteLine<std::string, unsigned int, std::string>(fieldFile,
                                                                               tempString,
                                                                               accuracy_m,
                                                                               tempString);

        tempString = Util::toUpper(tempString);
        if (tempString != "TRUE" &&
            tempString != "FALSE")
            throw GeneralClassicException("FM1DMagnetoStatic_fast::readFileHeader",
                                          "The third string on the first line of 1D field "
                                          "maps has to be either TRUE or FALSE");

        normalize_m = (tempString == "TRUE");
    }

    parsingPassed = parsingPassed &&
        interpreteLine<double, double, unsigned int>(fieldFile,
                                                     zBegin_m,
                                                     zEnd_m,
                                                     numberOfGridPoints_m);
    parsingPassed = parsingPassed &&
        interpreteLine<double, double, int>(fieldFile,
                                            rBegin_m,
                                            rEnd_m,
                                            tempInt);

    ++ numberOfGridPoints_m;

    if (accuracy_m > numberOfGridPoints_m)
        accuracy_m = numberOfGridPoints_m;

    return parsingPassed;
}

void FM1DMagnetoStatic_fast::stripFileHeader(std::ifstream &fieldFile) {

    std::string tempString;

    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
    getLine(fieldFile, tempString);
}

void FM1DMagnetoStatic_fast::prepareForMapCheck(std::vector<double> &fourierCoefs) {
    std::vector<double> zSampling(numberOfGridPoints_m);
    for (unsigned int zStepIndex = 0; zStepIndex < numberOfGridPoints_m; ++ zStepIndex)
        zSampling[zStepIndex] = deltaZ_m * zStepIndex;

    checkMap(accuracy_m,
             length_m,
             zSampling,
             fourierCoefs,
             onAxisFieldInterpolants_m,
             onAxisFieldAccel_m);
}

