#include "Fields/FM2DElectroStatic.h"
#include "Fields/Fieldmap.hpp"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <fstream>
#include <ios>

using namespace std;

FM2DElectroStatic::FM2DElectroStatic(std::string aFilename)
    : Fieldmap(aFilename),
      FieldstrengthEz_m(NULL),
      FieldstrengthEr_m(NULL) {
    ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type =  T2DElectroStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if (file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpreteLine<std::string, std::string>(file, tmpString, tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, std::string, std::string>(file,
                                                                                   tmpString,
                                                                                   tmpString,
                                                                                   tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("FM2DElectroStatic::FM2DElectroStatic",
                                              "The third string on the first line of 2D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        if (tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if (tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            cerr << "unknown orientation of 2D electrostatic fieldmap" << endl;
            parsing_passed = false;
        }

        for (long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double>(file, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if (!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("FM2DElectroStatic::FM2DElectroStatic",
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
            ++ num_gridpr_m;
            ++ num_gridpz_m;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM2DElectroStatic::~FM2DElectroStatic() {
    freeMap();
}

void FM2DElectroStatic::readMap() {
    if (FieldstrengthEz_m == NULL) {
        // declare variables and allocate memory
        ifstream in;
        std::string tmpString;
        double Ezmax = 0.0;

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];

        // read in and parse field map
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);


        if (swap_m) {
            for (int i = 0; i < num_gridpz_m; ++ i) {
                for (int j = 0; j < num_gridpr_m; ++ j) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthEr_m[i + j * num_gridpz_m],
                                                   FieldstrengthEz_m[i + j * num_gridpz_m]);
                }
                if (fabs(FieldstrengthEz_m[i]) > Ezmax) Ezmax = fabs(FieldstrengthEz_m[i]);
            }
        } else {
            for (int j = 0; j < num_gridpr_m; ++ j) {
                for (int i = 0; i < num_gridpz_m; ++ i) {
                    interpreteLine<double, double>(in,
                                                   FieldstrengthEz_m[i + j * num_gridpz_m],
                                                   FieldstrengthEr_m[i + j * num_gridpz_m]);
                }
            }

            for (int i = 0; i < num_gridpz_m; ++ i) {
                if (std::abs(FieldstrengthEz_m[i]) > Ezmax) {
                    Ezmax = std::abs(FieldstrengthEz_m[i]);
                }
            }
        }
        in.close();

        if (!normalize_m) 
            Ezmax = 1.0;

        // conversion MV/m to V/m and normalization to Ez_max = 1 MV/m
        for (int i = 0; i < num_gridpr_m * num_gridpz_m; ++ i) {
            FieldstrengthEz_m[i] *= 1e6 / Ezmax;
            FieldstrengthEr_m[i] *= 1e6 / Ezmax;
        }
        INFOMSG(typeset_msg("read in field map '" + Filename_m + "'", "info") << "\n"
                << endl);
    }
}

void FM2DElectroStatic::freeMap() {
    if (FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        FieldstrengthEz_m = NULL;
        delete[] FieldstrengthEr_m;
        FieldstrengthEr_m = NULL;

        INFOMSG(typeset_msg("freed field map '" + Filename_m + "'", "info") << "\n"
                << endl)
    }
}

bool FM2DElectroStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &/*B*/) const {
    // do bi-linear interpolation
    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = (int)floor((R(2)) / hz_m);
    const double leverz = (R(2) / hz_m) - indexz;

    if ((indexz < 0) || (indexz + 2 > num_gridpz_m))
        return false;
    if (indexr + 2 > num_gridpr_m)
        return true;

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;
    const double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
                           + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

    const double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
                           + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

    if (RR > 1e-10) {
        E(0) += EfieldR * R(0) / RR;
        E(1) += EfieldR * R(1) / RR;
    }
    E(2) += EfieldZ;
    return false;
}

bool FM2DElectroStatic::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM2DElectroStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}
void FM2DElectroStatic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM2DElectroStatic::swap() {
    if (swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DElectroStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (2D electrostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM2DElectroStatic::getFrequency() const {
    return 0.0;
}

void FM2DElectroStatic::setFrequency(double /*freq*/)
{ ;}
