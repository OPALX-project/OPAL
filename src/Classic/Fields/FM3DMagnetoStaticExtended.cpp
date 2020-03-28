#include "Fields/FM3DMagnetoStaticExtended.h"
#include "Fields/Fieldmap.hpp"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <cmath>
#include <fstream>
#include <ios>
#include <algorithm>

extern Inform *gmsg;

FM3DMagnetoStaticExtended::FM3DMagnetoStaticExtended(std::string aFilename):
    Fieldmap(aFilename),
    FieldstrengthBx_m(NULL),
    FieldstrengthBy_m(NULL),
    FieldstrengthBz_m(NULL) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpreteLine<std::string>(file, tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, std::string>(file, tmpString, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("FM3DMagnetoStaticExtended::FM3DMagnetoStaticExtended",
                                              "The second string on the first line of 3D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }
        parsing_passed = (parsing_passed &&
                          interpreteLine<double, double, unsigned int>(file, xbegin_m, xend_m, num_gridpx_m));
        parsing_passed = (parsing_passed &&
                          interpreteLine<double, double, unsigned int>(file, ybegin_m, yend_m, num_gridpy_m));
        parsing_passed = (parsing_passed &&
                          interpreteLine<double, double, unsigned int>(file, zbegin_m, zend_m, num_gridpz_m));

        for(unsigned long i = 0; (i < (num_gridpz_m + 1) * (num_gridpx_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double>(file, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("FM3DMagnetoStaticExtended::FM3DMagnetoStaticExtended(std::string)",
                                          "Format of fieldmap '" + Filename_m + "' didn't pass basic test");
        } else {
            // conversion from cm to m
            xbegin_m /= 100.;
            xend_m /= 100.;
            yend_m = std::max(std::abs(ybegin_m), yend_m) / 100.;
            ybegin_m = 0.0;
            zbegin_m /= 100.;
            zend_m /= 100.;
            length_m = zend_m - zbegin_m;

            hx_m = (xend_m - xbegin_m) / num_gridpx_m;
            hy_m = (yend_m - ybegin_m) / num_gridpy_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            // num spacings -> num grid points
            num_gridpx_m++;
            num_gridpy_m++;
            num_gridpz_m++;
        }
    } else {
        throw GeneralClassicException("FM3DMagnetoStaticExtended::FM3DMagnetoStaticExtended(std::string)",
                                      "Couldn't read fieldmap '" + Filename_m + "'");
    }
}

FM3DMagnetoStaticExtended::~FM3DMagnetoStaticExtended() {
    freeMap();
}

void FM3DMagnetoStaticExtended::readMap() {
    if(FieldstrengthBz_m == NULL) {
        // declare variables and allocate memory
        std::ifstream in;
        std::string tmpString;
        const size_t totalSize = num_gridpx_m * num_gridpy_m * num_gridpz_m;
        FieldstrengthBx_m = new double[totalSize];
        FieldstrengthBy_m = new double[totalSize];
        FieldstrengthBz_m = new double[totalSize];

        // read in and parse field map
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        for(unsigned int i = 0; i < num_gridpx_m; i++) {
            for(unsigned int k = 0; k < num_gridpz_m; k++) {
                unsigned long index = getIndex(i,0,k);
                interpreteLine<double>(in, FieldstrengthBy_m[index]);
            }
        }
        in.close();

        std::fill(FieldstrengthBx_m, FieldstrengthBx_m + totalSize, 0.0);
        std::fill(FieldstrengthBz_m, FieldstrengthBz_m + totalSize, 0.0);

        if (normalize_m) {
            double Bymax = 0.0;

            // find maximum field
            unsigned int centerX = static_cast<unsigned int>(std::round(-xbegin_m / hx_m));
            for(unsigned int k = 0; k < num_gridpz_m; ++ k) {
                double By = FieldstrengthBy_m[getIndex(centerX, 0, k)];
                if(std::abs(By) > std::abs(Bymax)) {
                    Bymax = By;
                }
            }

            // normalize field
            for(unsigned int i = 0; i < num_gridpx_m; i ++) {
                for(unsigned int k = 0; k < num_gridpz_m; k ++) {
                    unsigned long index = getIndex(i,0,k);
                    FieldstrengthBy_m[index] /= Bymax;
                }
            }
        }

        smoothData(FieldstrengthBy_m, 0);
        for (unsigned int j = 1; j < num_gridpy_m; j ++) {
            integrateBx(j);
            integrateBz(j);
            integrateBy(j);

            smoothData(FieldstrengthBx_m, j);
            smoothData(FieldstrengthBz_m, j);
            smoothData(FieldstrengthBy_m, j);
        }

        // saveField("data/" + Filename_m + ".ext", num_gridpy_m - 1);

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")  << endl);
    }
}

void FM3DMagnetoStaticExtended::integrateBx(unsigned j) {
    if (j == 1) {
        {
            unsigned int i = 0;
            // treat cells at i = 0, j = 1, k = 0:num_gridpz_m - 1;
            for (unsigned int k = 0; k < num_gridpz_m; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = 0.5 * hy_m / hx_m * (-3 * FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                                                4 * FieldstrengthBy_m[getIndex(i + 1, j - 1, k)] -
                                                                FieldstrengthBy_m[getIndex(i + 2, j - 1, k)]);
            }
        }
        // treat cells at i = 1:num_gridpx_m - 2, j = 1, k = 0:num_gridpz_m - 1;
        for(unsigned int i = 1; i < num_gridpx_m - 1; i ++) {
            for (unsigned int k = 0; k < num_gridpz_m; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = 0.5 * hy_m / hx_m * (FieldstrengthBy_m[getIndex(i + 1, j - 1, k)] -
                                                                FieldstrengthBy_m[getIndex(i - 1, j - 1, k)]);
            }
        }
        {
            unsigned int i = num_gridpx_m - 1;
            // treat cells at i = num_gridpx_m - 1, j = 1, k = 0:num_gridpz_m - 1;
            for (unsigned int k = 0; k < num_gridpz_m; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = 0.5 * hy_m / hx_m * (FieldstrengthBy_m[getIndex(i - 2, j - 1, k)] -
                                                                4 * FieldstrengthBy_m[getIndex(i - 1, j - 1, k)] +
                                                                3 * FieldstrengthBy_m[getIndex(i, j - 1, k)]);
            }
        }

    } else {
        { // treat cells at i = 0, j > 1, k = 0:num_gridpz_m - 1;
            unsigned int i = 0;
            for (unsigned int k = 0; k < num_gridpz_m; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = (4 * FieldstrengthBx_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBx_m[getIndex(i, j - 2, k)] +
                                            hy_m / hx_m * (-3 * FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                                           4 * FieldstrengthBy_m[getIndex(i + 1, j - 1, k)] -
                                                           FieldstrengthBy_m[getIndex(i + 2, j - 1, k)])) / 3;
            }
        }
        // treat cells at i = 1:num_gridpx_m - 2, j > 1, k = 0:num_gridpz_m - 1;
        for(unsigned int i = 1; i < num_gridpx_m - 1; i ++) {
            for (unsigned int k = 0; k < num_gridpz_m; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = (4 * FieldstrengthBx_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBx_m[getIndex(i, j - 2, k)] +
                                            hy_m / hx_m * (FieldstrengthBy_m[getIndex(i + 1, j - 1, k)] -
                                                           FieldstrengthBy_m[getIndex(i - 1, j - 1, k)])) / 3;
            }
        }
        { // treat cells at i = num_gridpx_m - 1, j > 1, k = 0:num_gridpz_m - 1;
            unsigned int i = num_gridpx_m - 1;
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBx_m[index] = (4 * FieldstrengthBx_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBx_m[getIndex(i, j - 2, k)] +
                                            hy_m / hx_m * (FieldstrengthBy_m[getIndex(i - 2, j - 1, k)] -
                                                           4 * FieldstrengthBy_m[getIndex(i - 1, j - 1, k)] +
                                                           3 * FieldstrengthBy_m[getIndex(i, j - 1, k)])) / 3;
            }
        }
    }
}

void FM3DMagnetoStaticExtended::integrateBz(unsigned j) {
    if (j == 1) {
        // treat cells at i = 0:num_gridpx_m - 1, j = 1, k = 1:num_gridpz_m - 2;
        for(unsigned int i = 0; i < num_gridpx_m; i ++) {
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = 0.5 * hy_m / hz_m * (FieldstrengthBy_m[getIndex(i, j - 1, k + 1)] -
                                                                FieldstrengthBy_m[getIndex(i, j - 1, k - 1)]);
            }

            { // treat cells at i = 0:num_gridpx_m - 1, j = 1, k = 0;
                unsigned int k = 0;
                unsigned int index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = 0.5 * hy_m / hz_m * (-3 * FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                                                4 * FieldstrengthBy_m[getIndex(i, j - 1, k + 1)] -
                                                                FieldstrengthBy_m[getIndex(i, j - 1, k + 2)]);
            }
            { // treat cells at i = 0:num_gridpx_m - 1, j = 1, k = num_gridpz_m - 1;
                unsigned int k = num_gridpz_m - 1;
                unsigned int index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = 0.5 * hy_m / hz_m * (FieldstrengthBy_m[getIndex(i, j - 1, k - 2)] -
                                                                4 * FieldstrengthBy_m[getIndex(i, j - 1, k - 1)] +
                                                                3 * FieldstrengthBy_m[getIndex(i,j - 1, k)]);
            }
        }
    } else {
        // treat cells at i = 0:num_gridpx_m - 1, j > 1, k = 1:num_gridpz_m - 2;
        for(unsigned int i = 0; i < num_gridpx_m; i ++) {
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = (4 * FieldstrengthBz_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBz_m[getIndex(i, j - 2, k)] +
                                            hy_m / hz_m * (FieldstrengthBy_m[getIndex(i, j - 1, k + 1)] -
                                                           FieldstrengthBy_m[getIndex(i, j - 1, k - 1)])) / 3;
            }

            { // treat cells at i = 0:num_gridpx_m - 1, j > 1, k = 0;
                unsigned int k = 0;
                unsigned int index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = (4 * FieldstrengthBz_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBz_m[getIndex(i, j - 2, k)] +
                                            hy_m / hz_m * (-3 * FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                                           4 * FieldstrengthBy_m[getIndex(i, j - 1, k + 1)] -
                                                           FieldstrengthBy_m[getIndex(i, j - 1, k + 2)])) / 3;
            }
            { // treat cells at i = 0:num_gridpx_m - 1, j > 1, k = num_gridpz_m - 1;
                unsigned int k = num_gridpz_m - 1;
                unsigned int index = getIndex(i,j,k);
                FieldstrengthBz_m[index] = (4 * FieldstrengthBz_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBz_m[getIndex(i, j - 2, k)] +
                                            hy_m / hz_m * (FieldstrengthBy_m[getIndex(i, j - 1, k - 2)] -
                                                           4 * FieldstrengthBy_m[getIndex(i, j - 1, k - 1)] +
                                                           3 * FieldstrengthBy_m[getIndex(i,j - 1, k)])) / 3;
            }
        }
    }
}

void FM3DMagnetoStaticExtended::integrateBy(unsigned j) {
    if (j == 1) {
        for (unsigned int i = 1; i < num_gridpx_m - 1; i ++) {
            { // treat cells at i = 1:num_gridpx_m - 2, j = 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m));
            }
            // treat cells at i = 1:num_gridpx_m - 2, j = 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m));
            }
            { // treat cells at i = 1:num_gridpx_m - 2, j = 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m));
            }
        }
        {
            unsigned int i = 0;
            { // treat cells at i = 0, j = 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m));
            }
            // treat cells at i = 0, j = 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m));
            }
            { // treat cells at i = 0, j = 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m));
            }
        }
        {
            unsigned int i = num_gridpx_m - 1;
            { // treat cells at i = 0, j = 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m));
            }
            // treat cells at i = 0, j = 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m));
            }
            { // treat cells at i = 0, j = 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (FieldstrengthBy_m[getIndex(i, j - 1, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m));
            }
        }
    } else {
        for (unsigned int i = 1; i < num_gridpx_m - 1; i ++) {
            { // treat cells at i = 1:num_gridpx_m - 2, j > 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m)) / 3;
            }
            // treat cells at i = 1:num_gridpx_m - 2, j > 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m)) / 3;
            }
            { // treat cells at i = 1:num_gridpx_m - 2, j > 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i - 1, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m)) / 3;
            }
        }
        {
            unsigned int i = 0;
            { // treat cells at i = 0, j > 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m)) / 3;
            }
            // treat cells at i = 0, j > 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m)) / 3;
            }
            { // treat cells at i = 0, j > 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((-3 * FieldstrengthBx_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBx_m[getIndex(i + 1, j, k)] -
                                                     FieldstrengthBx_m[getIndex(i + 2, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m)) / 3;
            }
        }
        {
            unsigned int i = num_gridpx_m - 1;
            { // treat cells at i = 0, j > 1, k = 0
                unsigned int k = 0;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (-3 * FieldstrengthBz_m[getIndex(i, j, k)] +
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k + 2)]) / hz_m)) / 3;
            }
            // treat cells at i = 0, j > 1, k = 1:num_gridpz_m - 2
            for (unsigned int k = 1; k < num_gridpz_m - 1; k ++) {
                unsigned long index = getIndex(i,j,k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k + 1)] -
                                                     FieldstrengthBz_m[getIndex(i, j, k - 1)]) / hz_m)) / 3;
            }
            { // treat cells at i = 0, j > 1, k = num_gridpz_m - 1
                unsigned int k = num_gridpz_m - 1;
                unsigned long index = getIndex(i, j, k);
                FieldstrengthBy_m[index] = (4 * FieldstrengthBy_m[getIndex(i, j - 1, k)] -
                                            FieldstrengthBy_m[getIndex(i, j - 2, k)] +
                                            hy_m * ((FieldstrengthBx_m[getIndex(i - 2, j, k)] -
                                                     4 * FieldstrengthBx_m[getIndex(i - 1, j, k)] +
                                                     3 * FieldstrengthBx_m[getIndex(i, j, k)]) / hx_m -
                                                    (FieldstrengthBz_m[getIndex(i, j, k - 2)] -
                                                     4 * FieldstrengthBz_m[getIndex(i, j, k - 1)] +
                                                     3 * FieldstrengthBz_m[getIndex(i, j, k)]) / hz_m)) / 3;
            }
        }
    }
}

void FM3DMagnetoStaticExtended::smoothData(double * data, unsigned j)
{
    const double offWeight = 0.1, sumWeightInv = 1.0 / (1.0 + 4 * (1 + offWeight) * offWeight);
    double *tmp = new double[num_gridpx_m * num_gridpy_m * num_gridpz_m];

    for (unsigned int i = 1; i < num_gridpx_m - 1; ++ i) {
        for (unsigned int k = 1; k < num_gridpz_m - 1; ++ k) {
            double sum = 0.0;
            for (int i2 = -1; i2 < 2; ++ i2) {
                for (int k2 = -1; k2 < 2; ++ k2) {
                    double weight = std::pow(offWeight, std::abs(i2) + std::abs(k2));
                    sum += weight * data[getIndex(i + i2, j, k + k2)];
                }
            }

            tmp[getIndex(i, j, k)] = sum * sumWeightInv;
        }
    }
    {
        const double sumWeightInv = 1.0 / (1.0 + (3 + 2 * offWeight) * offWeight);
        unsigned int i = 0;
        for (unsigned int k = 1; k < num_gridpz_m - 1; ++ k) {
            double sum = 0.0;
            for (int i2 = 0; i2 < 2; ++ i2) {
                for (int k2 = -1; k2 < 2; ++ k2) {
                    double weight = std::pow(offWeight, std::abs(i2) + std::abs(k2));
                    sum += weight * data[getIndex(i + i2, j, k + k2)];
                }
            }

            tmp[getIndex(i, j, k)] = sum * sumWeightInv;
        }
    }
    {
        const double sumWeightInv = 1.0 / (1.0 + (3 + 2 * offWeight) * offWeight);
        unsigned int i = num_gridpx_m - 1;
        for (unsigned int k = 1; k < num_gridpz_m - 1; ++ k) {
            double sum = 0.0;
            for (int i2 = -1; i2 < 1; ++ i2) {
                for (int k2 = -1; k2 < 2; ++ k2) {
                    double weight = std::pow(offWeight, std::abs(i2) + std::abs(k2));
                    sum += weight * data[getIndex(i + i2, j, k + k2)];
                }
            }

            tmp[getIndex(i, j, k)] = sum * sumWeightInv;
        }
    }
    {
        const double sumWeightInv = 1.0 / (1.0 + (3 + 2 * offWeight) * offWeight);
        unsigned int k = 0;
        for (unsigned int i = 1; i < num_gridpx_m - 1; ++ i) {
            double sum = 0.0;
            for (int i2 = -1; i2 < 2; ++ i2) {
                for (int k2 = 0; k2 < 2; ++ k2) {
                    double weight = std::pow(offWeight, std::abs(i2) + std::abs(k2));
                    sum += weight * data[getIndex(i + i2, j, k + k2)];
                }
            }

            tmp[getIndex(i, j, k)] = sum * sumWeightInv;
        }
    }
    {
        const double sumWeightInv = 1.0 / (1.0 + (3 + 2 * offWeight) * offWeight);
        unsigned int k = num_gridpz_m - 1;
        for (unsigned int i = 1; i < num_gridpx_m - 1; ++ i) {
            double sum = 0.0;
            for (int i2 = -1; i2 < 2; ++ i2) {
                for (int k2 = -1; k2 < 1; ++ k2) {
                    double weight = std::pow(offWeight, std::abs(i2) + std::abs(k2));
                    sum += weight * data[getIndex(i + i2, j, k + k2)];
                }
            }

            tmp[getIndex(i, j, k)] = sum * sumWeightInv;
        }
    }
    {
        unsigned long index = getIndex(0, j, 0);
        tmp[index] = data[index];
    }
    {
        unsigned long index = getIndex(num_gridpx_m - 1, j, 0);
        tmp[index] = data[index];
    }
    {
        unsigned long index = getIndex(num_gridpx_m - 1, j, num_gridpz_m - 1);
        tmp[index] = data[index];
    }
    {
        unsigned long index = getIndex(0, j, num_gridpz_m - 1);
        tmp[index] = data[index];
    }

    for (unsigned int i = 0; i < num_gridpx_m; ++ i) {
        for (unsigned int k = 0; k < num_gridpz_m; ++ k) {
            unsigned long index = getIndex(i, j, k);
            data[index] = tmp[index];
        }
    }

    delete[] tmp;
}

void FM3DMagnetoStaticExtended::saveField(const std::string &fname, unsigned int j) const
{
    std::ofstream out(fname);
    out.precision(6);
    double x = xbegin_m, y = j * hy_m;
    for (unsigned int i = 0; i < num_gridpx_m; ++ i, x += hx_m) {
        double z = zbegin_m;
        for (unsigned int k = 0; k < num_gridpz_m; ++ k, z += hz_m) {
            unsigned long index = getIndex(i, j, k);
            out << std::setw(14) << x
                << std::setw(14) << y
                << std::setw(14) << z
                << std::setw(14) << FieldstrengthBx_m[index]
                << std::setw(14) << FieldstrengthBy_m[index]
                << std::setw(14) << FieldstrengthBz_m[index]
                << "\n";
        }
    }

    out.close();
}

void FM3DMagnetoStaticExtended::freeMap() {
    if(FieldstrengthBz_m != NULL) {
        delete[] FieldstrengthBx_m;
        delete[] FieldstrengthBy_m;
        delete[] FieldstrengthBz_m;

        FieldstrengthBx_m = NULL;
        FieldstrengthBy_m = NULL;
        FieldstrengthBz_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

Vector_t FM3DMagnetoStaticExtended::interpolateTrilinearly(const Vector_t &X) const {
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

double FM3DMagnetoStaticExtended::getWeightedData(double *data, const IndexTriplet &idx, unsigned short corner) const {
    unsigned short switchX = ((corner & HX) >> 2), switchY = ((corner & HY) >> 1), switchZ = (corner & HZ);
    double factorX = 0.5 + (1 - 2 * switchX) * (0.5 - idx.weight(0));
    double factorY = 0.5 + (1 - 2 * switchY) * (0.5 - idx.weight(1));
    double factorZ = 0.5 + (1 - 2 * switchZ) * (0.5 - idx.weight(2));

    unsigned long i = idx.i + switchX, j = idx.j + switchY, k = idx.k + switchZ;

    return factorX * factorY * factorZ * data[getIndex(i, j, k)];
}

bool FM3DMagnetoStaticExtended::getFieldstrength(const Vector_t &R, Vector_t &/*E*/, Vector_t &B) const {
    if (isInside(R)) {
        Vector_t suppB = interpolateTrilinearly(R);
        suppB(0) *= copysign(1, R(1));
        suppB(2) *= copysign(1, R(1));

        B += suppB;
    }

    return false;
}

bool FM3DMagnetoStaticExtended::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM3DMagnetoStaticExtended::getFieldDimensions(double &zBegin, double &zEnd, double &xbegin, double &xend) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    xbegin = xbegin_m;
    xend = xend_m;
}

void FM3DMagnetoStaticExtended::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {
    xIni = xbegin_m;
    xFinal = xend_m;
    yIni = -yend_m;
    yFinal = yend_m;
    zIni = zbegin_m;
    zFinal = zend_m;
}

void FM3DMagnetoStaticExtended::swap() {
}

void FM3DMagnetoStaticExtended::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D magnetostatic, extended); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM3DMagnetoStaticExtended::getFrequency() const {
    return 0.0;
}

void FM3DMagnetoStaticExtended::setFrequency(double /*freq*/)
{ ;}
