#include "Fields/Astra1DMagnetoStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_fft_real.h"

#include <fstream>
#include <ios>

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

Astra1DMagnetoStatic::Astra1DMagnetoStatic(std::string aFilename)
    : Fieldmap(aFilename),
      FourCoefs_m(NULL) {
    ifstream file;
    int skippedValues = 0;
    std::string tmpString;
    double tmpDouble;
    double tmpDouble2;

    Type = TAstraMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpreteLine<std::string, int>(file, tmpString, accuracy_m);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, int, std::string>(file, tmpString, accuracy_m, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("Astra1DMagnetoStatic::Astra1DMagnetoStatic",
                                              "The third string on the first line of 1D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double>(file, zbegin_m, tmpDouble);

        tmpDouble2 = zbegin_m;
        while(!file.eof() && parsing_passed) {
            parsing_passed = interpreteLine<double, double>(file, zend_m, tmpDouble, false);
            if(zend_m - tmpDouble2 > 1e-10) {
                tmpDouble2 = zend_m;
            } else if(parsing_passed) {
                ++ skippedValues;
            }
        }

        file.close();
        num_gridpz_m = lines_read_m - 2 - skippedValues;
        lines_read_m = 0;

        if(!parsing_passed && !file.eof()) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("Astra1DMagnetoStatic::Astra1DMagnetoStatic",
                                          "An error occured when reading the fieldmap '" + Filename_m + "'");
        }
        length_m = 2.0 * num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

Astra1DMagnetoStatic::~Astra1DMagnetoStatic() {
    freeMap();
}

void Astra1DMagnetoStatic::readMap() {
    if(FourCoefs_m == NULL) {
        // declare variables and allocate memory
        ifstream in;

        bool parsing_passed = true;

        std::string tmpString;

        double Bz_max = 0.0;
        double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
        double tmpDouble = zbegin_m - dz;

        double *RealValues = new double[2 * num_gridpz_m];
        double *zvals = new double[num_gridpz_m];

        gsl_spline *Bz_interpolant = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
        gsl_interp_accel *Bz_accel = gsl_interp_accel_alloc();

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

        FourCoefs_m = new double[2 * accuracy_m - 1];

        // read in and parse field map
        in.open(Filename_m.c_str());
        getLine(in, tmpString);

        for(int i = 0; i < num_gridpz_m && parsing_passed;/* skip increment on i here */) {
            parsing_passed = interpreteLine<double, double>(in, zvals[i], RealValues[i]);
            // the sequence of z-position should be strictly increasing
            // drop sampling points that don't comply to this
            if(zvals[i] - tmpDouble > 1e-10) {
                if(fabs(RealValues[i]) > Bz_max) {
                    Bz_max = fabs(RealValues[i]);
                }
                tmpDouble = zvals[i];
                ++ i; // increment i only if sampling point is accepted
            }
        }
        in.close();

        gsl_spline_init(Bz_interpolant, zvals, RealValues, num_gridpz_m);

        // get equidistant sampling from the, possibly, non-equidistant sampling
        // using cubic spline for this
        int ii = num_gridpz_m;
        for(int i = 0; i < num_gridpz_m - 1; ++ i, ++ ii) {
            double z = zbegin_m + dz * i;
            RealValues[ii] = gsl_spline_eval(Bz_interpolant, z, Bz_accel);
        }
        RealValues[ii ++] = RealValues[num_gridpz_m - 1];
        // prepend mirror sampling points such that field values are periodic for sure
        -- ii; // ii == 2*num_gridpz_m at the moment
        for(int i = 0; i < num_gridpz_m; ++ i, -- ii) {
            RealValues[i] = RealValues[ii];
        }

        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        if (!normalize_m) {
            Bz_max = 1.0;
        }
        // normalize to Bz_max = 1 A/m
        FourCoefs_m[0] = RealValues[0] / (Bz_max * 2. * num_gridpz_m);
        for(int i = 1; i < 2 * accuracy_m - 1; i++) {
            FourCoefs_m[i] = RealValues[i] / (Bz_max * num_gridpz_m);
        }

        gsl_spline_free(Bz_interpolant);
        gsl_interp_accel_free(Bz_accel);

        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] zvals;
        delete[] RealValues;

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

void Astra1DMagnetoStatic::freeMap() {
    if(FourCoefs_m != NULL) {

        delete[] FourCoefs_m;
        FourCoefs_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

bool Astra1DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &/*E*/, Vector_t &B) const {
    // do fourier interpolation in z-direction
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    const double kz = two_pi * R(2) / length_m + Physics::pi;

    double ez = FourCoefs_m[0];
    double ezp = 0.0;
    double ezpp = 0.0;
    double ezppp = 0.0;
    double somefactor_base, somefactor;
    double coskzl;
    double sinkzl;

    int n = 1;
    for(int l = 1; l < accuracy_m ; l++, n += 2) {
        somefactor_base =  two_pi / length_m * l;    // = \frac{d(kz*l)}{dz}
        somefactor = 1.0;
        coskzl = cos(kz * l);
        sinkzl = sin(kz * l);
        ez    += (FourCoefs_m[n] * coskzl - FourCoefs_m[n + 1] * sinkzl);
        somefactor *= somefactor_base;
        ezp   += somefactor * (-FourCoefs_m[n] * sinkzl - FourCoefs_m[n + 1] * coskzl);
        somefactor *= somefactor_base;
        ezpp  += somefactor * (-FourCoefs_m[n] * coskzl + FourCoefs_m[n + 1] * sinkzl);
        somefactor *= somefactor_base;
        ezppp += somefactor * (FourCoefs_m[n] * sinkzl + FourCoefs_m[n + 1] * coskzl);
    }
    // expand the field off-axis
    const double BfieldR = -ezp / 2. + ezppp / 16. * RR2;

    B(0) += BfieldR * R(0);
    B(1) += BfieldR * R(1);
    B(2) += ez - ezpp * RR2 / 4.;

    return false;
}

bool Astra1DMagnetoStatic::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void Astra1DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &/*rBegin*/, double &/*rEnd*/) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
}
void Astra1DMagnetoStatic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void Astra1DMagnetoStatic::swap()
{ }

void Astra1DMagnetoStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double Astra1DMagnetoStatic::getFrequency() const {
    return 0.0;
}

void Astra1DMagnetoStatic::setFrequency(double /*freq*/)
{ }
