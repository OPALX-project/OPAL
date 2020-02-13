#include "Fields/Astra1D_fast.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

#include <fstream>
#include <ios>

Astra1D_fast::Astra1D_fast(std::string aFilename):
    Fieldmap(aFilename) {
    normalize_m = true;
}

Astra1D_fast::~Astra1D_fast() {
    freeMap();
}

void Astra1D_fast::readMap()
{ }

void Astra1D_fast::freeMap() {
    if(onAxisField_m != NULL) {
        for(int i = 0; i < 4; ++i) {
            gsl_spline_free(onAxisInterpolants_m[i]);
            gsl_interp_accel_free(onAxisAccel_m[i]);
        }

        delete[] onAxisField_m;
        onAxisField_m = NULL;
        delete[] zvals_m;
        zvals_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

void Astra1D_fast::getOnaxisEz(std::vector<std::pair<double, double> > & /*F*/)
{ }

bool Astra1D_fast::determineNumSamplingPoints(std::ifstream &file) {
    double tmpDouble, tmpDouble2;
    unsigned int skippedValues = 0;
    bool passed = interpreteLine<double, double>(file, zbegin_m, tmpDouble);

    tmpDouble2 = zbegin_m;
    while (!file.eof() && passed) {
        passed = interpreteLine<double, double>(file, zend_m, tmpDouble, false);
        if (zend_m - tmpDouble2 > 1e-10) {
            tmpDouble2 = zend_m;
        } else if (passed) {
            ++ skippedValues;
        }
    }

    num_gridpz_m = lines_read_m - numHeaderLines_m - skippedValues;
    lines_read_m = 0;

    return passed || file.eof();
}

double Astra1D_fast::readFieldData(std::ifstream &file) {
    double tmpDouble = zbegin_m - hz_m;
    double maxValue = 0;
    int nsp;

    for (nsp = 0; nsp < num_gridpz_m; ) {
        interpreteLine<double, double>(file, zvals_m[nsp], onAxisField_m[nsp]);

        if(zvals_m[nsp] - tmpDouble > 1e-10) {
            if(std::abs(onAxisField_m[nsp]) > maxValue) {
                maxValue = std::abs(onAxisField_m[nsp]);
            }
            tmpDouble = zvals_m[nsp];
            ++ nsp;
        }
    }
    num_gridpz_m = nsp;
    hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);

    if (!normalize_m) return 1.0;

    return maxValue;
}

void Astra1D_fast::normalizeFieldData(double maxValue) {
    int ii = num_gridpz_m - 1;
    for(int i = 0; i < num_gridpz_m; ++ i, -- ii) {
        onAxisField_m[i] /= maxValue;
        zvals_m[ii] -= zvals_m[0];
    }
}

std::vector<double> Astra1D_fast::getEvenlyDistributedSamplingPoints() {
    std::vector<double> zvals(num_gridpz_m, 0.0);

    for(int i = 1; i < num_gridpz_m - 1; ++ i) {
        zvals[i] = zvals[i - 1] + hz_m;
    }
    zvals[num_gridpz_m - 1] = zvals_m[num_gridpz_m - 1];

    return zvals;
}

std::vector<double> Astra1D_fast::interpolateFieldData(std::vector<double> &samplingPoints) {
    std::vector<double> realValues(num_gridpz_m);

    onAxisInterpolants_m[0] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
    onAxisAccel_m[0] = gsl_interp_accel_alloc();

    gsl_spline_init(onAxisInterpolants_m[0], zvals_m, onAxisField_m, num_gridpz_m);

    for (int i = 0; i < num_gridpz_m - 1; ++ i) {
        realValues[i] = gsl_spline_eval(onAxisInterpolants_m[0], samplingPoints[i], onAxisAccel_m[0]);
    }
    realValues[num_gridpz_m - 1] = onAxisField_m[num_gridpz_m - 1];

    return realValues;
}

std::vector<double> Astra1D_fast::computeFourierCoefficients(int accuracy,
                                                                    std::vector<double> &evenFieldSampling) {
    std::vector<double> fourierCoefficients(2 * accuracy - 1);
    std::vector<double> RealValues(2 * num_gridpz_m, 0.0);
    int ii = num_gridpz_m, iii = num_gridpz_m - 1;

    for(int i = 0; i < num_gridpz_m; ++ i, ++ ii, -- iii) {
        RealValues[ii] = evenFieldSampling[i];
        RealValues[iii] = evenFieldSampling[i];
    }

    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

    gsl_fft_real_transform(&RealValues[0], 1, 2 * num_gridpz_m, real, work);

    gsl_fft_real_workspace_free(work);
    gsl_fft_real_wavetable_free(real);

    // normalize to Ez_max = 1 MV/m
    fourierCoefficients[0] = RealValues[0] / (2 * num_gridpz_m);
    for(int i = 1; i < 2 * accuracy - 1; i++) {
        fourierCoefficients[i] = RealValues[i] / num_gridpz_m;
    }

    return fourierCoefficients;
}

void Astra1D_fast::computeFieldDerivatives(std::vector<double> & fourierComponents, int accuracy) {
    double interiorDerivative, base;
    double coskzl, sinkzl, z = 0.0;
    std::vector<double> zvals(num_gridpz_m);
    std::vector<double> higherDerivatives[3];
    higherDerivatives[0].resize(num_gridpz_m, 0.0);
    higherDerivatives[1].resize(num_gridpz_m, 0.0);
    higherDerivatives[2].resize(num_gridpz_m, 0.0);

    for(int i = 0; i < num_gridpz_m; ++ i, z += hz_m) {
        zvals[i] = z;
        const double kz = Physics::two_pi * (zvals[i] / length_m + 0.5);

        for(int l = 1; l < accuracy; ++l) {
            int coefIndex = 2 * l - 1;
            base = Physics::two_pi / length_m * l;
            interiorDerivative = base;
            coskzl = cos(kz * l);
            sinkzl = sin(kz * l);

            higherDerivatives[0][i] += interiorDerivative * (-fourierComponents[coefIndex] * sinkzl
                                                             - fourierComponents[coefIndex + 1] * coskzl);
            interiorDerivative *= base;
            higherDerivatives[1][i] += interiorDerivative * (-fourierComponents[coefIndex] * coskzl
                                                             + fourierComponents[coefIndex + 1] * sinkzl);
            interiorDerivative *= base;
            higherDerivatives[2][i] += interiorDerivative * (fourierComponents[coefIndex] * sinkzl
                                                             + fourierComponents[coefIndex + 1] * coskzl);
        }
    }

    for(int i = 1; i < 4; ++i) {
        onAxisInterpolants_m[i] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
        onAxisAccel_m[i] = gsl_interp_accel_alloc();
        gsl_spline_init(onAxisInterpolants_m[i], &zvals[0], &higherDerivatives[i - 1][0], num_gridpz_m);
    }
}

