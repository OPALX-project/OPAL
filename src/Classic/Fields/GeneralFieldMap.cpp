#include "Fields/GeneralFieldMap.h"

#include "Utilities/GeneralClassicException.h"
#include "Physics/Physics.h"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
// #include <boost/range/combine.hpp>
// #include <boost/foreach.hpp>

#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

#include <map>
#include <list>
#include <regex>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

namespace {
    std::pair<std::vector<double>, std::vector<double>>
    smoothData(unsigned int accuracy,
               std::vector<double> const& EzRealValues,
               std::vector<double> const& EzImagValues) {
        const unsigned int size = EzRealValues.size();
        const unsigned int totalSize = 2 * size - 1;
        gsl_fft_real_wavetable *waveTable = gsl_fft_real_wavetable_alloc(totalSize);
        gsl_fft_real_workspace *workSpace = gsl_fft_real_workspace_alloc(totalSize);

        // Reflect field about minimum z value to ensure that it is periodic.
        std::vector<double> fieldDataReflected(totalSize);

        for (unsigned int dataIndex = 0; dataIndex < size; ++ dataIndex) {
            fieldDataReflected[size - 1 + dataIndex]
                = EzRealValues[dataIndex];
            if (dataIndex > 0) {
                fieldDataReflected[size - 1 - dataIndex]
                    = EzRealValues[dataIndex];
            }
        }

        gsl_fft_real_transform(&(fieldDataReflected[0]), 1,
                               totalSize,
                               waveTable, workSpace);

        std::pair<std::vector<double>, std::vector<double>> fourierCoefs;
        fourierCoefs.first.resize(2 * accuracy - 1);
        fourierCoefs.second.resize(2 * accuracy - 1);

        fourierCoefs.first[0] = fieldDataReflected[0] / totalSize;
        for (unsigned int coefIndex = 1; coefIndex < 2 * accuracy - 1; ++ coefIndex) {
            fourierCoefs.first[coefIndex] = 2.0 * fieldDataReflected[coefIndex] / totalSize;
        }

        for (unsigned int dataIndex = 0; dataIndex < size; ++ dataIndex) {
            fieldDataReflected[size - 1 + dataIndex]
                = EzImagValues[dataIndex];
            if (dataIndex != 0) {
                fieldDataReflected[size - 1 - dataIndex]
                    = EzImagValues[dataIndex];
            }
        }

        gsl_fft_real_transform(&(fieldDataReflected[0]), 1,
                               totalSize,
                               waveTable, workSpace);

        fourierCoefs.second[0] = fieldDataReflected[0] / totalSize;
        for (unsigned int coefIndex = 1; coefIndex < 2 * accuracy - 1; ++ coefIndex) {
            fourierCoefs.second[coefIndex] = 2.0 * fieldDataReflected[coefIndex] / totalSize;
        }
        gsl_fft_real_workspace_free(workSpace);
        gsl_fft_real_wavetable_free(waveTable);

        return fourierCoefs;
    }

    void computeOnAxisFieldAndDerivatives(unsigned int const size, double const dz,
                                          std::vector<double> const& fourierCoefs,
                                          std::vector<double> & onAxisField,
                                          std::vector<double> & onAxisFieldP,
                                          std::vector<double> & onAxisFieldPP,
                                          std::vector<double> & onAxisFieldPPP) {
        onAxisField.resize(size);
        onAxisFieldP.resize(size);
        onAxisFieldPP.resize(size);
        onAxisFieldPPP.resize(size);

        unsigned int accuracy = fourierCoefs.size() / 2;
        double lengthOverall = 2 * size * dz - dz;
        for (unsigned int zStepIndex = 0; zStepIndex < size; ++ zStepIndex) {

            double z = dz * (zStepIndex - 0.5) ;
            double kz = Physics::two_pi * z / lengthOverall + Physics::pi;
            onAxisField[zStepIndex] = fourierCoefs[0];
            onAxisFieldP[zStepIndex] = 0.0;
            onAxisFieldPP[zStepIndex] = 0.0;
            onAxisFieldPPP[zStepIndex] = 0.0;

            int coefIndex = 1;
            for (unsigned int n = 1; n < accuracy; ++ n) {

                double kn = n * Physics::two_pi / lengthOverall;
                double coskzn = cos(kz * n);
                double sinkzn = sin(kz * n);

                double derivCoeff = kn;
                onAxisField[zStepIndex] += (fourierCoefs[coefIndex] * coskzn
                                            - fourierCoefs[coefIndex + 1] * sinkzn);

                onAxisFieldP[zStepIndex] += derivCoeff * (-fourierCoefs[coefIndex] * sinkzn
                                                          - fourierCoefs[coefIndex + 1] * coskzn);

                derivCoeff *= kn;
                onAxisFieldPP[zStepIndex] += derivCoeff * (-fourierCoefs[coefIndex] * coskzn
                                                           + fourierCoefs[coefIndex + 1] * sinkzn);

                derivCoeff *= kn;
                onAxisFieldPPP[zStepIndex] += derivCoeff * (fourierCoefs[coefIndex] * sinkzn
                                                            + fourierCoefs[coefIndex + 1] * coskzn);

                coefIndex += 2;
            }
        }
    }
}

GeneralFieldMap::GeneralFieldMap(std::string const& filename):
    Fieldmap(filename),
    accuracy_m(0)
{
    Type = T1GeneralFieldMap;
    readMap();
}

GeneralFieldMap::~GeneralFieldMap() {
    freeMap();
}

GeneralFieldMap::FieldData GeneralFieldMap::readFieldData() const {
    FieldData result;
    std::ifstream input(Filename_m);
    std::list<std::string> lines;
    std::string tmp;
    char lineBuffer[1024];
    std::regex rexp("//.*");
    std::streampos pos = input.tellg();

    while(!input.eof()) {
        input >> tmp;
        if (input.rdstate()) break;

        input.seekg(pos);
        input.getline(lineBuffer, 1024);
        pos = input.tellg();

        std::string line(lineBuffer);
        line = std::regex_replace(line, rexp, "");

        if (line.length() == 0) continue;

        lines.push_back(line);
    }
    if (lines.size() < 3u) {
        throw GeneralClassicException("GeneralFieldMap::readMap",
                                      "File '" + Filename_m + "' only contains " + std::to_string(lines.size()) + " lines");
    }
    std::string const& columnDescription = *std::next(lines.begin(), 2);
    std::istringstream splitter(columnDescription);
    std::map<std::string, unsigned int> columnNumByName;
    unsigned int numColumns = 0;
    while(!splitter.eof()) {
        std::string columnName;
        splitter >> columnName;
        if (splitter.bad()) break;

        columnNumByName.insert(std::make_pair(columnName, numColumns));
        ++ numColumns;

        if (splitter.eof()) break;
    }

    for (std::string col: {"z", "Ez_r", "Ez_i"}) {
        if (columnNumByName.find(col) == columnNumByName.end()) {
            throw GeneralClassicException("GeneralFieldMap::readMap",
                                          "No column for " + col + " could be found in file '" + Filename_m + "'");
        }
    }

    std::vector<double> lineValues(numColumns);
    const unsigned int zIndex = columnNumByName["z"];
    const unsigned int EzRealIndex = columnNumByName["Ez_r"];
    const unsigned int EzImagIndex = columnNumByName["Ez_i"];

    std::list<std::string>::const_iterator it = lines.begin();
    std::stringstream line(*it);
    line >> tmp >> result.accuracy;
    line.clear();

    line.str(*(++ it));
    line >> result.frequency;
    result.frequency *= 1e6 * Physics::two_pi;

    it = std::next(it, 2);
    while (it != lines.end()) {
        line.clear();
        line.str(*it);

        for (unsigned int i = 0; i < numColumns; ++ i) {
            double value;
            line >> value;
            if (line.bad() || (i + 1 < numColumns && line.eof())) {
                throw GeneralClassicException("GeneralFieldMap::readMap",
                                              "Line num " + std::to_string(result.zValues.size()) +
                                              " in file '" + Filename_m + "' doesn't contain" +
                                              " enough values.\nIts content is '" + *it +
                                              "'.\nFailed at column " + std::to_string(i) + ".");
            }

            lineValues[i] = value;
        }

        result.zValues.push_back(lineValues[zIndex]);
        result.EzRealValues.push_back(lineValues[EzRealIndex]);
        result.EzImagValues.push_back(lineValues[EzImagIndex]);

        ++ it;
        lines.pop_front();
    }
    return result;
}

void GeneralFieldMap::readMap() {
    if (onAxisFieldInterpolants_m.first) {
        return;
    }

    namespace fs = boost::filesystem;
    if (!fs::exists(Filename_m)) {
        throw GeneralClassicException("GeneralFieldMap::readMap",
                                      "File '" + Filename_m + "' doesn't exist");
    }

    FieldData fieldData = readFieldData();
    accuracy_m = fieldData.accuracy;
    frequency_m = fieldData.frequency;
    zBegin_m = fieldData.zValues.front();
    length_m = fieldData.zValues.back() - zBegin_m;

    std::cout << "Found " << fieldData.zValues.size() << " lines of data in the file" << std::endl;
    std::cout << "The frequency is " << frequency_m << " MHz" << std::endl;

    allocateMemory(fieldData.zValues.size());
    evenlyDistributeData(fieldData.zValues, fieldData.EzRealValues, fieldData.EzImagValues);
    std::pair<std::vector<double>, std::vector<double>> data = ::smoothData(accuracy_m,
                                                                            fieldData.EzRealValues,
                                                                            fieldData.EzImagValues);
    initInterpolation(fieldData.zValues, data);
}

void GeneralFieldMap::freeMap() {

}

void GeneralFieldMap::allocateMemory(unsigned int numGridPoints) {
    onAxisFieldInterpolants_m.real.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldInterpolants_m.imag.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPInterpolants_m.real.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPInterpolants_m.imag.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPPInterpolants_m.real.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPPInterpolants_m.imag.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPPPInterpolants_m.real.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));
    onAxisFieldPPPInterpolants_m.imag.reset(gsl_spline_alloc(gsl_interp_steffen, numGridPoints));

    onAxisFieldAccel_m.reset(gsl_interp_accel_alloc());
}

void GeneralFieldMap::evenlyDistributeData(std::vector<double>& zValues,
                                           std::vector<double>& EzRealValues,
                                           std::vector<double>& EzImagValues) const {
    const unsigned int size = zValues.size();
    gsl_spline_init(onAxisFieldInterpolants_m.real.get(), &(zValues[0]), &(EzRealValues[0]), size);
    gsl_spline_init(onAxisFieldInterpolants_m.imag.get(), &(zValues[0]), &(EzImagValues[0]), size);

    double hz = length_m / (size - 1);
    for (unsigned int i = 1; i + 1u < size; ++ i) {
        zValues[i] = zValues.front() + i * hz;
        EzRealValues[i] = gsl_spline_eval(onAxisFieldInterpolants_m.real.get(), zValues[i], onAxisFieldAccel_m.get());
        EzImagValues[i] = gsl_spline_eval(onAxisFieldInterpolants_m.imag.get(), zValues[i], onAxisFieldAccel_m.get());
    }
}

void GeneralFieldMap::initInterpolation(std::vector<double> const& zValues,
                                        std::pair<std::vector<double>, std::vector<double>> const& fourCoefs) {

    const unsigned int size = zValues.size();
    std::vector<double> onAxisField;
    std::vector<double> onAxisFieldP;
    std::vector<double> onAxisFieldPP;
    std::vector<double> onAxisFieldPPP;
    std::vector<double> shiftedZValues(size);

    const double dz = length_m / (size - 1);
    for (unsigned int i = 0; i < size; ++ i) {
        shiftedZValues[i] = i * dz;
    }

    computeOnAxisFieldAndDerivatives(size,
                                     dz,
                                     fourCoefs.first,
                                     onAxisField,
                                     onAxisFieldP,
                                     onAxisFieldPP,
                                     onAxisFieldPPP);
    gsl_spline_init(onAxisFieldInterpolants_m.real.get(), &(shiftedZValues[0]),
                    &(onAxisField[0]), size);
    gsl_spline_init(onAxisFieldPInterpolants_m.real.get(), &(shiftedZValues[0]),
                    &(onAxisFieldP[0]), size);
    gsl_spline_init(onAxisFieldPPInterpolants_m.real.get(), &(shiftedZValues[0]),
                    &(onAxisFieldPP[0]), size);
    gsl_spline_init(onAxisFieldPPPInterpolants_m.real.get(), &(shiftedZValues[0]),
                    &(onAxisFieldPPP[0]), size);

    computeOnAxisFieldAndDerivatives(size,
                                     dz,
                                     fourCoefs.second,
                                     onAxisField,
                                     onAxisFieldP,
                                     onAxisFieldPP,
                                     onAxisFieldPPP);

    gsl_spline_init(onAxisFieldInterpolants_m.imag.get(), &(shiftedZValues[0]),
                    &(onAxisField[0]), size);
    gsl_spline_init(onAxisFieldPInterpolants_m.imag.get(), &(shiftedZValues[0]),
                    &(onAxisFieldP[0]), size);
    gsl_spline_init(onAxisFieldPPInterpolants_m.imag.get(), &(shiftedZValues[0]),
                    &(onAxisFieldPP[0]), size);
    gsl_spline_init(onAxisFieldPPPInterpolants_m.imag.get(), &(shiftedZValues[0]),
                    &(onAxisFieldPPP[0]), size);
}

bool GeneralFieldMap::getFieldstrength(const Vector_t &R, ComplexVector_t &E, ComplexVector_t &B) const {
    double z = R(2);
    if (z < 0.0 || z > length_m) {
        return false;
    }

    std::vector<std::pair<double, double>> fieldComponents;
    computeFieldOnAxis(z, fieldComponents);
    computeFieldOffAxis(R, E, B, fieldComponents);

    return false;
}

void GeneralFieldMap::computeFieldOnAxis(double z,
                                         std::vector<std::pair<double, double>> &fieldComponents) const {

    fieldComponents.push_back(std::make_pair(gsl_spline_eval(onAxisFieldInterpolants_m.real.get(),
                                                             z, onAxisFieldAccel_m.get()),
                                             gsl_spline_eval(onAxisFieldInterpolants_m.imag.get(),
                                                             z, onAxisFieldAccel_m.get())));
    fieldComponents.push_back(std::make_pair(gsl_spline_eval(onAxisFieldPInterpolants_m.real.get(),
                                                             z, onAxisFieldAccel_m.get()),
                                             gsl_spline_eval(onAxisFieldPInterpolants_m.imag.get(),
                                                             z, onAxisFieldAccel_m.get())));
    fieldComponents.push_back(std::make_pair(gsl_spline_eval(onAxisFieldPPInterpolants_m.real.get(),
                                                             z, onAxisFieldAccel_m.get()),
                                             gsl_spline_eval(onAxisFieldPPInterpolants_m.imag.get(),
                                                             z, onAxisFieldAccel_m.get())));
    fieldComponents.push_back(std::make_pair(gsl_spline_eval(onAxisFieldPPPInterpolants_m.real.get(),
                                                             z, onAxisFieldAccel_m.get()),
                                             gsl_spline_eval(onAxisFieldPPPInterpolants_m.imag.get(),
                                                             z, onAxisFieldAccel_m.get())));
}

void GeneralFieldMap::computeFieldOffAxis(const Vector_t &R,
                                          ComplexVector_t &E,
                                          ComplexVector_t &B,
                                          std::vector<std::pair<double, double>> const& fieldComponents) const {

    double radiusSq = pow(R(0), 2.0) + pow(R(1), 2.0);
    auto get = [](unsigned int n, std::pair<double, double> & p) -> double& {return (n == 0? p.first: p.second);};
    auto const_get = [](unsigned int n, std::pair<double, double> const& p) -> double {return (n == 0? p.first: p.second);};

    std::pair<double, double> transverseEFactor, longitudinalEFactor;
    std::pair<double, double> transverseBFactor;

    for (unsigned int i: {0, 1}) {
        get(i, transverseEFactor) = (const_get(i, fieldComponents[0])
                                     * (8.0 - radiusSq * twoPiOverLambdaSq_m)
                                     - radiusSq * const_get(i, fieldComponents[3])) / 16.0;
        get(i, longitudinalEFactor) = (const_get(i, fieldComponents[0])
                                       * (4.0 - radiusSq * twoPiOverLambdaSq_m)
                                       - radiusSq * const_get(i, fieldComponents[2])) / 4.0;
        get(i, transverseBFactor) = ((const_get(i, fieldComponents[0])
                                      * (8.0 - radiusSq * twoPiOverLambdaSq_m)
                                      - radiusSq * const_get(i, fieldComponents[2])) / 16.0
                                     * twoPiOverLambdaSq_m / frequency_m);
    }


    E.real()(0) += - R(0) * transverseEFactor.first;
    E.imag()(0) += - R(0) * transverseEFactor.second;
    E.real()(1) += - R(1) * transverseEFactor.first;
    E.imag()(1) += - R(1) * transverseEFactor.second;

    E.real()(2) += longitudinalEFactor.first;
    E.imag()(2) += longitudinalEFactor.second;

    B.real()(0) += - R(1) * transverseBFactor.first;
    B.real()(1) +=   R(0) * transverseBFactor.first;
    B.imag()(0) += - R(1) * transverseBFactor.second;
    B.imag()(1) +=   R(0) * transverseBFactor.second;
}

bool GeneralFieldMap::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/,
                                         Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void GeneralFieldMap::getFieldDimensions(double &zBegin, double &zEnd) const {
    zBegin = zBegin_m;
    zEnd = zBegin_m + length_m;
}
void GeneralFieldMap::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/,
                                         double &/*yIni*/, double &/*yFinal*/,
                                         double &/*zIni*/, double &/*zFinal*/) const {

}

void GeneralFieldMap::swap() {

}

void GeneralFieldMap::getInfo(Inform *) {

}

double GeneralFieldMap::getFrequency() const {
    return frequency_m;
}

void GeneralFieldMap::setFrequency(double /*freq*/) {

}
