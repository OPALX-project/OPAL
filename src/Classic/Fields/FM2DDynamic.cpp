#include "Fields/FM2DDynamic.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <fstream>
#include <ios>

extern Inform *gmsg;

// using namespace std;
using Physics::mu_0;

FM2DDynamic::FM2DDynamic(std::string aFilename)
    : Fieldmap(aFilename),
      FieldstrengthEz_m(NULL),
      FieldstrengthEr_m(NULL),
      FieldstrengthBt_m(NULL) {
    std::ifstream file;
    std::string tmpString;
    double tmpDouble;

    Type = T2DDynamic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = true;
        try {
            parsing_passed = interpreteLine<std::string, std::string>(file, tmpString, tmpString);
        } catch (GeneralClassicException &e) {
            parsing_passed = interpreteLine<std::string, std::string, std::string>(file, tmpString, tmpString, tmpString);

            tmpString = Util::toUpper(tmpString);
            if (tmpString != "TRUE" &&
                tmpString != "FALSE")
                throw GeneralClassicException("FM2DDynamic::FM2DDynamic",
                                              "The third string on the first line of 2D field "
                                              "maps has to be either TRUE or FALSE");

            normalize_m = (tmpString == "TRUE");
        }

        if(tmpString == "ZX") {
            swap_m = true;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, frequency_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        } else if(tmpString == "XZ") {
            swap_m = false;
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, frequency_m);
            parsing_passed = parsing_passed &&
                             interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        } else {
            std::cerr << "unknown orientation of 2D dynamic fieldmap" << std::endl;
            parsing_passed = false;
            zbegin_m = 0.0;
            zend_m = -1e-3;
        }

        for(long i = 0; (i < (num_gridpz_m + 1) * (num_gridpr_m + 1)) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed && interpreteLine<double, double, double, double>(file, tmpDouble, tmpDouble, tmpDouble, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
            throw GeneralClassicException("FM2DDynamic::FM2DDynamic",
                                          "An error occured when reading the fieldmap '" + Filename_m + "'");
        } else {
            // convert MHz to Hz and frequency to angular frequency
            frequency_m *= Physics::two_pi * 1e6;

            // convert cm to m
            rbegin_m /= 100.0;
            rend_m /= 100.0;
            zbegin_m /= 100.0;
            zend_m /= 100.0;

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


FM2DDynamic::~FM2DDynamic() {
    freeMap();
}

void FM2DDynamic::readMap() {
    if(FieldstrengthEz_m == NULL) {
        // declare variables and allocate memory
        std::ifstream in;
        std::string tmpString;
        double tmpDouble, Ezmax = 0.0;

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthBt_m = new double[num_gridpz_m * num_gridpr_m];

        // read in field map and parse it
        in.open(Filename_m.c_str());
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);
        getLine(in, tmpString);

        if(swap_m) {
            for(int i = 0; i < num_gridpz_m; i++) {
                for(int j = 0; j < num_gridpr_m; j++) {
                    interpreteLine<double, double, double, double>(in,
                            FieldstrengthEr_m[i + j * num_gridpz_m],
                            FieldstrengthEz_m[i + j * num_gridpz_m],
                            FieldstrengthBt_m[i + j * num_gridpz_m],
                            tmpDouble);
                }
            }
        } else {
            for(int j = 0; j < num_gridpr_m; j++) {
                for(int i = 0; i < num_gridpz_m; i++) {
                    interpreteLine<double, double, double, double>(in,
                            FieldstrengthEz_m[i + j * num_gridpz_m],
                            FieldstrengthEr_m[i + j * num_gridpz_m],
                            tmpDouble,
                            FieldstrengthBt_m[i + j * num_gridpz_m]);
                }
            }
        }

        in.close();


        if (normalize_m) {
            // find maximum field
            for(int i = 0; i < num_gridpz_m; ++ i) {
                if(std::abs(FieldstrengthEz_m[i]) > Ezmax) {
                    Ezmax = std::abs(FieldstrengthEz_m[i]);
                }
            }
        } else {
            Ezmax = 1.0;
        }

        for(int i = 0; i < num_gridpr_m * num_gridpz_m; i++) {
            FieldstrengthEz_m[i] *= 1.e6 / Ezmax; // conversion MV/m to V/m and normalization
            FieldstrengthEr_m[i] *= 1.e6 / Ezmax;
            FieldstrengthBt_m[i] *= mu_0 / Ezmax; // H -> B
        }

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl);

    }
}

void FM2DDynamic::freeMap() {
    if(FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        FieldstrengthEz_m = NULL;
        delete[] FieldstrengthEr_m;
        FieldstrengthEr_m = NULL;
        delete[] FieldstrengthBt_m;
        FieldstrengthBt_m = NULL;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
                << endl);
    }
}

bool FM2DDynamic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do bi-linear interpolation
    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = RR / hr_m - indexr;

    const int indexz = (int)floor(R(2) / hz_m);
    const double leverz = R(2) / hz_m - indexz;

    if((indexz < 0) || (indexz + 2 > num_gridpz_m))
        return false;
    if(indexr + 2 > num_gridpr_m)
        return true;

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;

    double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
                     + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
                     + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
                     + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

    double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
                     + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
                     + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
                     + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

    double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
                     + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index1 + 1]
                     + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index2]
                     + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];

    if(RR > 1e-10) {
        E(0) += EfieldR * R(0) / RR;
        E(1) += EfieldR * R(1) / RR;
        B(0) -= BfieldT * R(1) / RR;
        B(1) += BfieldT * R(0) / RR;
    }
    E(2) += EfieldZ;

    return false;
}

bool FM2DDynamic::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM2DDynamic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}
void FM2DDynamic::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM2DDynamic::swap() {
    if(swap_m) swap_m = false;
    else swap_m = true;
}

void FM2DDynamic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (2D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM2DDynamic::getFrequency() const {
    return frequency_m;
}

void FM2DDynamic::setFrequency(double freq) {
    frequency_m = freq;
}

void FM2DDynamic::getOnaxisEz(std::vector<std::pair<double, double> > & F) {
    double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    F.resize(num_gridpz_m);

    for(int i = 0; i < num_gridpz_m; ++ i) {
        F[i].first = dz * i;
        F[i].second = FieldstrengthEz_m[i] / 1e6;

    }
}
