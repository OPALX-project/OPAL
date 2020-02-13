#include "Fields/FM3DMagnetoStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

extern Inform *gmsg;

using namespace std;
using Physics::mu_0;

FM3DMagnetoStatic::FM3DMagnetoStatic(std::string aFilename):
    Fieldmap(aFilename),
    FieldstrengthBz_m(NULL),
    FieldstrengthBx_m(NULL),
    FieldstrengthBy_m(NULL) {

    std::string tmpString;
    double tmpDouble;

    Type = T3DMagnetoStatic;
    ifstream file(Filename_m.c_str());

    if(file.good()) {
        bool parsing_passed = true;
        try {
            interpreteLine<std::string>(file, tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, std::string>(file, tmpString, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("FM3DMagnetoStatic::FM3DMagnetoStatic",
                                              "The second string on the first line of 3D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, unsigned int>(file, xbegin_m, xend_m, num_gridpx_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, unsigned int>(file, ybegin_m, yend_m, num_gridpy_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, unsigned int>(file, zbegin_m, zend_m, num_gridpz_m);

        for(unsigned long i = 0; (i < (num_gridpz_m + 1) * (num_gridpx_m + 1) * (num_gridpy_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file,
                                                    tmpDouble,
                                                    tmpDouble,
                                                    tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("FM3DMagnetoStatic::FM3DMagnetoStatic",
                                          "An error occured when reading the fieldmap '" + Filename_m + "'");
        } else {
            xbegin_m /= 100.0;
            xend_m /= 100.0;
            ybegin_m /= 100.0;
            yend_m /= 100.0;
            zbegin_m /= 100.0;
            zend_m /= 100.0;

            hx_m = (xend_m - xbegin_m) / num_gridpx_m;
            hy_m = (yend_m - ybegin_m) / num_gridpy_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            ++ num_gridpx_m;
            ++ num_gridpy_m;
            ++ num_gridpz_m;

        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}


FM3DMagnetoStatic::~FM3DMagnetoStatic() {
    freeMap();
}

void FM3DMagnetoStatic::readMap() {
    if(FieldstrengthBz_m == NULL) {

        ifstream in(Filename_m.c_str());
        std::string tmpString;
        const size_t totalSize = num_gridpx_m * num_gridpy_m * num_gridpz_m;

        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        FieldstrengthBz_m = new double[totalSize];
        FieldstrengthBx_m = new double[totalSize];
        FieldstrengthBy_m = new double[totalSize];

        long ii = 0;
        for(unsigned int i = 0; i < num_gridpx_m; ++ i) {
            for(unsigned int j = 0; j < num_gridpy_m; ++ j) {
                for(unsigned int k = 0; k < num_gridpz_m; ++ k) {
                    interpreteLine<double>(in,
                                           FieldstrengthBx_m[ii],
                                           FieldstrengthBy_m[ii],
                                           FieldstrengthBz_m[ii]);
                    ++ ii;
                }
            }
        }
        in.close();

        if (normalize_m) {
            double Bymax = 0.0;
            // find maximum field
            unsigned int centerX = static_cast<unsigned int>(std::floor(-xbegin_m / hx_m + 0.5));
            unsigned int centerY = static_cast<unsigned int>(std::floor(-ybegin_m / hy_m + 0.5));
            for(unsigned int k = 0; k < num_gridpz_m; ++ k) {
                double By = FieldstrengthBy_m[getIndex(centerX, centerY, k)];
                if(std::abs(By) > Bymax) {
                    Bymax = std::abs(By);
                }
            }

            // normalize field
            for(size_t i = 0; i < totalSize; ++ i) {
                FieldstrengthBx_m[i] /= Bymax;
                FieldstrengthBy_m[i] /= Bymax;
                FieldstrengthBz_m[i] /= Bymax;
            }
        }

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl);
    }
}

void FM3DMagnetoStatic::freeMap() {
    if(FieldstrengthBz_m != NULL) {
        delete[] FieldstrengthBz_m;
        delete[] FieldstrengthBx_m;
        delete[] FieldstrengthBy_m;

        FieldstrengthBz_m = NULL;
        FieldstrengthBx_m = NULL;
        FieldstrengthBy_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
                << endl);
    }
}

bool FM3DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &/*E*/, Vector_t &B) const {
    if (isInside(R))
        B += interpolateTrilinearly(R);

    return false;
}

Vector_t FM3DMagnetoStatic::interpolateTrilinearly(const Vector_t &X) const {
    IndexTriplet idx = getIndex(X);
    Vector_t B(0.0);

    B(0) = (getWeightedData(FieldstrengthBx_m, idx, LX|LY|LZ) +
            getWeightedData(FieldstrengthBx_m, idx, LX|LY|HZ) +
            getWeightedData(FieldstrengthBx_m, idx, LX|HY|LZ) +
            getWeightedData(FieldstrengthBx_m, idx, LX|HY|HZ) +
            getWeightedData(FieldstrengthBx_m, idx, HX|LY|LZ) +
            getWeightedData(FieldstrengthBx_m, idx, HX|LY|HZ) +
            getWeightedData(FieldstrengthBx_m, idx, HX|HY|LZ) +
            getWeightedData(FieldstrengthBx_m, idx, HX|HY|HZ));

    B(1) = (getWeightedData(FieldstrengthBy_m, idx, LX|LY|LZ) +
            getWeightedData(FieldstrengthBy_m, idx, LX|LY|HZ) +
            getWeightedData(FieldstrengthBy_m, idx, LX|HY|LZ) +
            getWeightedData(FieldstrengthBy_m, idx, LX|HY|HZ) +
            getWeightedData(FieldstrengthBy_m, idx, HX|LY|LZ) +
            getWeightedData(FieldstrengthBy_m, idx, HX|LY|HZ) +
            getWeightedData(FieldstrengthBy_m, idx, HX|HY|LZ) +
            getWeightedData(FieldstrengthBy_m, idx, HX|HY|HZ));

    B(2) = (getWeightedData(FieldstrengthBz_m, idx, LX|LY|LZ) +
            getWeightedData(FieldstrengthBz_m, idx, LX|LY|HZ) +
            getWeightedData(FieldstrengthBz_m, idx, LX|HY|LZ) +
            getWeightedData(FieldstrengthBz_m, idx, LX|HY|HZ) +
            getWeightedData(FieldstrengthBz_m, idx, HX|LY|LZ) +
            getWeightedData(FieldstrengthBz_m, idx, HX|LY|HZ) +
            getWeightedData(FieldstrengthBz_m, idx, HX|HY|LZ) +
            getWeightedData(FieldstrengthBz_m, idx, HX|HY|HZ));

    return B;
}

double FM3DMagnetoStatic::getWeightedData(double *data, const IndexTriplet &idx, unsigned short corner) const {
    unsigned short switchX = ((corner & HX) >> 2), switchY = ((corner & HY) >> 1), switchZ = (corner & HZ);
    double factorX = 0.5 + (1 - 2 * switchX) * (0.5 - idx.weight(0));
    double factorY = 0.5 + (1 - 2 * switchY) * (0.5 - idx.weight(1));
    double factorZ = 0.5 + (1 - 2 * switchZ) * (0.5 - idx.weight(2));

    unsigned long i = idx.i + switchX, j = idx.j + switchY, k = idx.k + switchZ;

    return factorX * factorY * factorZ * data[getIndex(i, j, k)];
}

bool FM3DMagnetoStatic::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM3DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}
void FM3DMagnetoStatic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM3DMagnetoStatic::swap() {}

void FM3DMagnetoStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D magnetostatic) "
           << " xini= " << xbegin_m << " xfinal= " << xend_m
           << " yini= " << ybegin_m << " yfinal= " << yend_m
           << " zini= " << zbegin_m << " zfinal= " << zend_m << " (m) " << endl;
}
