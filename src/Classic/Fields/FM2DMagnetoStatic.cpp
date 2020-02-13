#include "Fields/FM2DMagnetoStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <fstream>
#include <ios>

extern Inform *gmsg;

FM2DMagnetoStatic::FM2DMagnetoStatic(std::string aFilename):
    Fieldmap(aFilename),
    FieldstrengthBz_m(NULL),
    FieldstrengthBr_m(NULL) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpreteLine<std::string, std::string>(file,
                                                                      tmpString,
                                                                      tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, std::string, std::string>(file,
                                                                                   tmpString,
                                                                                   tmpString,
                                                                                   tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("FM2DMagnetoStatic::FM2DMagnetoStatic",
                                              "The third string on the first line of 2D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        if(tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if(tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            std::cerr << "unknown orientation of 2D magnetostatic fieldmap" << std::endl;
            parsing_passed = false;
        }

        for(long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double>(file, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("FM2DMagnetoStatic::FM2DMagnetoStatic",
                                          "An error occured when reading the fieldmap '" + Filename_m + "'");
        } else {
            // conversion from cm to m
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;

            hr_m = (rend_m - rbegin_m) / num_gridpr_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            // num spacings -> num grid points
            num_gridpr_m++;
            num_gridpz_m++;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM2DMagnetoStatic::~FM2DMagnetoStatic() {
    freeMap();
}

void FM2DMagnetoStatic::readMap() {
    if(FieldstrengthBz_m == NULL) {
        // declare variables and allocate memory
        std::ifstream in;
        std::string tmpString;

        FieldstrengthBz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthBr_m = new double[num_gridpz_m * num_gridpr_m];

        // read in and parse field map
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        if(swap_m) {
            for(int i = 0; i < num_gridpz_m; i++) {
                for(int j = 0; j < num_gridpr_m; j++) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthBr_m[i + j * num_gridpz_m],
                                                   FieldstrengthBz_m[i + j * num_gridpz_m]);
                }
            }
        } else {
            for(int j = 0; j < num_gridpr_m; j++) {
                for(int i = 0; i < num_gridpz_m; i++) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthBz_m[i + j * num_gridpz_m],
                                                   FieldstrengthBr_m[i + j * num_gridpz_m]);
                }
            }
        }
        in.close();

        if (normalize_m) {
            double Bzmax = 0.0;
            // find maximum field
            for(int i = 0; i < num_gridpz_m; ++ i) {
                if(std::abs(FieldstrengthBz_m[i]) > Bzmax) {
                    Bzmax = std::abs(FieldstrengthBz_m[i]);
                }
            }

            // normalize field
            for(int i = 0; i < num_gridpr_m * num_gridpz_m; ++ i) {
                FieldstrengthBz_m[i] /= Bzmax;
                FieldstrengthBr_m[i] /= Bzmax;
            }
        }

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")  << endl);
    }
}

void FM2DMagnetoStatic::freeMap() {
    if(FieldstrengthBz_m != NULL) {
        delete[] FieldstrengthBz_m;
        delete[] FieldstrengthBr_m;

        FieldstrengthBz_m = NULL;
        FieldstrengthBr_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

bool FM2DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &/*E*/, Vector_t &B) const {
    // do bi-linear interpolation
    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    const int indexr = abs((int)floor(RR / hr_m));
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = abs((int)floor((R(2)) / hz_m));
    const double leverz = (R(2) / hz_m) - indexz;

    if((indexz < 0) || (indexz + 2 > num_gridpz_m))
        return false;

    if(indexr + 2 > num_gridpr_m)
        return true;

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;
    const double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index2]
                           + leverz         * leverr         * FieldstrengthBr_m[index2 + 1];

    const double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index2]
                           + leverz         * leverr         * FieldstrengthBz_m[index2 + 1];

    if(RR != 0) {
        B(0) += BfieldR * R(0) / RR;
        B(1) += BfieldR * R(1) / RR;
    }
    B(2) += BfieldZ;
    return false;
}

bool FM2DMagnetoStatic::getFieldDerivative(const Vector_t &R, Vector_t &/*E*/, Vector_t &B, const DiffDirection &dir) const {

    double BfieldR, BfieldZ;

    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = (int)floor((R(2)) / hz_m);
    const double leverz = (R(2) / hz_m) - indexz;

    if((indexz < 0) || (indexz + 2 > num_gridpz_m))
        return false;

    if(indexr + 2 > num_gridpr_m)
        return true;

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;
    if(dir == DZ) {
        if(indexz > 0) {
            if(indexz < num_gridpz_m - 1) {
                BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1 - 1] - 2. * FieldstrengthBr_m[index1] + FieldstrengthBr_m[index1 + 1]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                            + (FieldstrengthBr_m[index1 + 1] - FieldstrengthBr_m[index1 - 1]) / (2.*hz_m))
                          + leverr         * ((FieldstrengthBr_m[index2 - 1] - 2. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index2 + 1]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                              + (FieldstrengthBr_m[index2 + 1] - FieldstrengthBr_m[index2 - 1]) / (2.*hz_m));

                BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1 - 1] - 2. * FieldstrengthBz_m[index1] + FieldstrengthBz_m[index1 + 1]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                            + (FieldstrengthBz_m[index1 + 1] - FieldstrengthBz_m[index1 - 1]) / (2.*hz_m))
                          + leverr         * ((FieldstrengthBz_m[index2 - 1] - 2. * FieldstrengthBz_m[index2] + FieldstrengthBz_m[index2 + 1]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                              + (FieldstrengthBz_m[index2 + 1] - FieldstrengthBz_m[index2 - 1]) / (2.*hz_m));
            } else {
                BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index1 - 1] + FieldstrengthBr_m[index1 - 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                            - (3. * FieldstrengthBr_m[index1] - 4. * FieldstrengthBr_m[index1 - 1] + FieldstrengthBr_m[index1 - 2]) / (2.*hz_m))
                          + leverr         * ((FieldstrengthBr_m[index2] - 2. * FieldstrengthBr_m[index2 - 1] + FieldstrengthBr_m[index2 - 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                              - (3. * FieldstrengthBr_m[index2]  - 4. * FieldstrengthBr_m[index1 - 1] + FieldstrengthBr_m[index2 - 2]) / (2.*hz_m));

                BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index1 - 1] + FieldstrengthBz_m[index1 - 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                            - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index1 - 1] + FieldstrengthBz_m[index1 - 2]) / (2.*hz_m))
                          + leverr         * ((FieldstrengthBz_m[index2] - 2. * FieldstrengthBz_m[index2 - 1] + FieldstrengthBz_m[index2 - 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                              - (3. * FieldstrengthBz_m[index2]  - 4. * FieldstrengthBr_m[index1 - 1] + FieldstrengthBz_m[index2 - 2]) / (2.*hz_m));
            }
        } else {
            BfieldR = (1.0 - leverr) * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBr_m[index1 + 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                        - (3. * FieldstrengthBr_m[index1] - 4. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBr_m[index1 + 2]) / (2.*hz_m))
                      + leverr         * ((FieldstrengthBr_m[index2] - 2. * FieldstrengthBr_m[index2 + 1] + FieldstrengthBr_m[index2 + 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                          - (3. * FieldstrengthBr_m[index2]  - 4. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBr_m[index2 + 2]) / (2.*hz_m));

            BfieldZ = (1.0 - leverr) * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index1 + 1] + FieldstrengthBz_m[index1 + 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                        - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBz_m[index1 + 2]) / (2.*hz_m))
                      + leverr         * ((FieldstrengthBz_m[index2] - 2. * FieldstrengthBz_m[index2 + 1] + FieldstrengthBz_m[index2 + 2]) * (R(2) - indexz * hz_m) / (hz_m * hz_m)
                                          - (3. * FieldstrengthBz_m[index2]  - 4. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBz_m[index2 + 2]) / (2.*hz_m));
        }
    } else {
        if(indexr > 0) {
            const int index_1 = index1 - num_gridpz_m;
            if(indexr < num_gridpr_m - 1) {
                BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index_1 + 1] - 2. * FieldstrengthBr_m[index1 + 1] + FieldstrengthBr_m[index2 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                            + (FieldstrengthBr_m[index2 + 1] - FieldstrengthBr_m[index_1 + 1]) / (2.*hr_m))
                          + leverz         * ((FieldstrengthBr_m[index_1] - 2. * FieldstrengthBr_m[index1] + FieldstrengthBr_m[index2]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                              + (FieldstrengthBr_m[index2] - FieldstrengthBr_m[index_1]) / (2.*hr_m));

                BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index_1 + 1] - 2. * FieldstrengthBz_m[index1 + 1] + FieldstrengthBz_m[index2 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                            + (FieldstrengthBz_m[index2 + 1] - FieldstrengthBz_m[index_1 + 1]) / (2.*hr_m))
                          + leverz         * ((FieldstrengthBz_m[index_1] - 2. * FieldstrengthBz_m[index1] + FieldstrengthBz_m[index2]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                              + (FieldstrengthBz_m[index2] - FieldstrengthBz_m[index_1]) / (2.*hr_m));
            } else {
                const int index_2 = index_1 - num_gridpz_m;
                BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index1 + 1] - 2. * FieldstrengthBr_m[index_1 + 1] + FieldstrengthBr_m[index_2 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                            - (3. * FieldstrengthBr_m[index1 + 1] - 4. * FieldstrengthBr_m[index_1 + 1] + FieldstrengthBr_m[index_2 + 1]) / (2.*hr_m))
                          + leverz         * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index_1] + FieldstrengthBr_m[index_2]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                              - (3. * FieldstrengthBr_m[index1]  - 4. * FieldstrengthBr_m[index_1] + FieldstrengthBr_m[index_2]) / (2.*hr_m));

                BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index1 + 1] - 2. * FieldstrengthBz_m[index_1 + 1] + FieldstrengthBz_m[index_2 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                            - (3. * FieldstrengthBz_m[index1 + 1]  - 4. * FieldstrengthBr_m[index_1 + 1] + FieldstrengthBz_m[index_2 + 1]) / (2.*hr_m))
                          + leverz         * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index_1] + FieldstrengthBz_m[index_2]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                              - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index_1] + FieldstrengthBz_m[index_2]) / (2.*hr_m));
            }
        } else {
            const int index3 = index2 + num_gridpz_m;
            BfieldR = (1.0 - leverz) * ((FieldstrengthBr_m[index1 + 1] - 2. * FieldstrengthBr_m[index2 + 1] + FieldstrengthBr_m[index3 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                        - (3. * FieldstrengthBr_m[index1 + 1] - 4. * FieldstrengthBr_m[index2 + 1] + FieldstrengthBr_m[index3 + 1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBr_m[index1] - 2. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index3]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                          - (3. * FieldstrengthBr_m[index1]  - 4. * FieldstrengthBr_m[index2] + FieldstrengthBr_m[index3]) / (2.*hr_m));

            BfieldZ = (1.0 - leverz) * ((FieldstrengthBz_m[index1 + 1] - 2. * FieldstrengthBz_m[index2 + 1] + FieldstrengthBz_m[index3 + 1]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                        - (3. * FieldstrengthBz_m[index1 + 1]  - 4. * FieldstrengthBr_m[index2 + 1] + FieldstrengthBz_m[index3 + 1]) / (2.*hr_m))
                      + leverz         * ((FieldstrengthBz_m[index1] - 2. * FieldstrengthBz_m[index2] + FieldstrengthBz_m[index3]) * (RR - indexr * hr_m) / (hr_m * hr_m)
                                          - (3. * FieldstrengthBz_m[index1]  - 4. * FieldstrengthBr_m[index2] + FieldstrengthBz_m[index3]) / (2.*hr_m));
        }

    }

    if(RR > 1.e-12) {
        B(0) += BfieldR * R(0) / RR;
        B(1) += BfieldR * R(1) / RR;
    }
    B(2) += BfieldZ;
    return false;
}

void FM2DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM2DMagnetoStatic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM2DMagnetoStatic::swap() {
    if(swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DMagnetoStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (2D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM2DMagnetoStatic::getFrequency() const {
    return 0.0;
}

void FM2DMagnetoStatic::setFrequency(double /*freq*/)
{ ;}
