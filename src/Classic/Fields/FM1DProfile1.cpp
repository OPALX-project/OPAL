#include "Fields/FM1DProfile1.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

#include <fstream>
#include <ios>

FM1DProfile1::FM1DProfile1(std::string aFilename):
    Fieldmap(aFilename),
    entranceParameter1_m(0.0),
    entranceParameter2_m(0.0),
    entranceParameter3_m(0.0),
    exitParameter1_m(0.0),
    exitParameter2_m(0.0),
    exitParameter3_m(0.0),
    polyOrderEntry_m(5),
    polyOrderExit_m(5),
    gapHeight_m(0.02),
    sBegin_m(0.0),
    sEnd_m(0.0) {

    // Read file header information. Set field type.
    Type = T1DProfile1;

    // Read file. Check if we are using the default field profile.
    if(Filename_m == "1DPROFILE1-DEFAULT") {

        /*
          * Use default field profile coefficients:
          *
          * polyOrderEntry_m = 5
          * polyOrderExit_m  = 5
          * gapHeight_m      = 0.02 m
          *
          */
        polyOrderEntry_m = 5;
        polyOrderExit_m = 5;
        gapHeight_m = 0.02;

        entranceParameter1_m = -0.1;
        entranceParameter2_m = 0.0;
        entranceParameter3_m = 0.1;

        exitParameter1_m = -0.1;
        exitParameter2_m = 0.0;
        exitParameter3_m = 0.1;

    } else {

        std::ifstream inputFile(Filename_m.c_str());

        if(inputFile.good()) {

            int tempInt = 0;
            std::string tempString = "";
            double tempDouble = 0.0;

            bool parsingPassed = interpreteLine<std::string, int, int, double>
                                 (inputFile,
                                  tempString,
                                  polyOrderEntry_m,
                                  polyOrderExit_m,
                                  gapHeight_m);

            parsingPassed = parsingPassed &&
                            interpreteLine<double, double, double, int>
                            (inputFile,
                             entranceParameter1_m,
                             entranceParameter2_m,
                             entranceParameter3_m,
                             tempInt,
                             false);

            parsingPassed = parsingPassed &&
                            interpreteLine<double, double, double, int>
                            (inputFile,
                             exitParameter1_m,
                             exitParameter2_m,
                             exitParameter3_m,
                             tempInt);

            for(int index = 0;
                (index < polyOrderEntry_m + polyOrderExit_m + 2) && parsingPassed;
                index++)
                parsingPassed = parsingPassed &&
                                interpreteLine<double>(inputFile, tempDouble);

            parsingPassed = parsingPassed && interpreteEOF(inputFile);

            inputFile.close();

            if(!parsingPassed) {

                disableFieldmapWarning();
                sBegin_m = 0.0;
                sEnd_m = sBegin_m - 1.0e-3;

            } else {

                // Convert from cm to m. Make sure gap is positive.
                entranceParameter1_m /= 100.0;
                entranceParameter2_m /= 100.0;
                entranceParameter3_m /= 100.0;

                exitParameter1_m /= 100.0;
                exitParameter2_m /= 100.0;
                exitParameter3_m /= 100.0;

                gapHeight_m = std::abs(gapHeight_m / 100.0);

            }

        } else {

            // No field map file or field map file is somehow flawed.
            noFieldmapWarning();
            sBegin_m = 0.0;
            sEnd_m = sBegin_m - 1.0e-3;

        }
    }
}

FM1DProfile1::~FM1DProfile1() {
}

void FM1DProfile1::readMap() {

    if(!engeCoeffsEntry_m.empty())
        engeCoeffsEntry_m.clear();

    if(!engeCoeffsExit_m.empty())
        engeCoeffsExit_m.clear();

    if(Filename_m == "1DPROFILE1-DEFAULT") {

        engeCoeffsEntry_m.push_back(0.478959);
        engeCoeffsEntry_m.push_back(1.911289);
        engeCoeffsEntry_m.push_back(-1.185953);
        engeCoeffsEntry_m.push_back(1.630554);
        engeCoeffsEntry_m.push_back(-1.082657);
        engeCoeffsEntry_m.push_back(0.318111);

        engeCoeffsExit_m.push_back(0.478959);
        engeCoeffsExit_m.push_back(1.911289);
        engeCoeffsExit_m.push_back(-1.185953);
        engeCoeffsExit_m.push_back(1.630554);
        engeCoeffsExit_m.push_back(-1.082657);
        engeCoeffsExit_m.push_back(0.31811);

    } else {

        std::ifstream inputFile(Filename_m.c_str());

        int tempInt;
        std::string tempString;
        double tempDouble;

        interpreteLine<std::string, int, int, double>(inputFile,
                tempString,
                tempInt,
                tempInt,
                tempDouble);
        interpreteLine<double, double, double, int>(inputFile,
                tempDouble,
                tempDouble,
                tempDouble,
                tempInt);
        interpreteLine<double, double, double, int>(inputFile,
                tempDouble,
                tempDouble,
                tempDouble,
                tempInt);

        for(int index = 0; index < polyOrderEntry_m + 1; index++) {
            interpreteLine<double>(inputFile,  tempDouble);
            engeCoeffsEntry_m.push_back(tempDouble);
        }

        for(int index = 0; index < polyOrderExit_m + 1; index++) {
            interpreteLine<double>(inputFile, tempDouble);
            engeCoeffsExit_m.push_back(tempDouble);
        }

        inputFile.close();

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")
                << "\n" << endl);

    }

    if (computeEntranceFringe(entranceParameter1_m) < computeEntranceFringe(entranceParameter3_m)) {
        for (int index = 0; index < polyOrderEntry_m + 1; ++ index) {
            if (index % 2 == 0) continue;
            engeCoeffsEntry_m[index] *= -1;
        }
    }

    if (computeExitFringe(exitParameter1_m) < computeExitFringe(exitParameter3_m)) {
        for (int index = 0; index < polyOrderExit_m + 1; ++ index) {
            if (index % 2 == 0) continue;
            engeCoeffsExit_m[index] *= -1;
        }
    }
}

void FM1DProfile1::freeMap() {
}

bool FM1DProfile1::getFieldstrength(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/) const {

    /*
     * For this type of field map, the elements who use it calculate the field
     * amplitude using the field map parameters.
     */
    return true;

}

bool FM1DProfile1::getFieldDerivative(const Vector_t &/*R*/,
                                      Vector_t &/*E*/,
                                      Vector_t &/*B*/,
                                      const DiffDirection &/*dir*/) const {
    return false;
}

void FM1DProfile1::getFieldDimensions(double &sBegin,
                                      double &sEnd,
                                      double &/*rBegin*/,
                                      double &/*rEnd*/) const {
    sBegin = sBegin_m;
    sEnd = sEnd_m;
}
void FM1DProfile1::getFieldDimensions(double &/*xIni*/,
                                      double &/*xFinal*/,
                                      double &/*yIni*/,
                                      double &/*yFinal*/,
                                      double &/*zIni*/,
                                      double &/*zFinal*/) const {

}

void FM1DProfile1::swap()
{}

void FM1DProfile1::getInfo(Inform *msg) {
    (*msg) << Filename_m
           << " (1D Profile type 1)"
           << endl;
}

double FM1DProfile1::getFrequency() const {
    return 0.0;
}

void FM1DProfile1::setFrequency(double /*freq*/)
{}

void FM1DProfile1::get1DProfile1EngeCoeffs(std::vector<double> &engeCoeffsEntry,
        std::vector<double> &engeCoeffsExit) {
    engeCoeffsEntry = engeCoeffsEntry_m;
    engeCoeffsExit = engeCoeffsExit_m;

}

void FM1DProfile1::get1DProfile1EntranceParam(double &entranceParameter1,
        double &entranceParameter2,
        double &entranceParameter3) {
    entranceParameter1 = entranceParameter1_m;
    entranceParameter2 = entranceParameter2_m;
    entranceParameter3 = entranceParameter3_m;
}

void FM1DProfile1::get1DProfile1ExitParam(double &exitParameter1,
        double &exitParameter2,
        double &exitParameter3) {
    exitParameter1 = exitParameter1_m;
    exitParameter2 = exitParameter2_m;
    exitParameter3 = exitParameter3_m;
}

double FM1DProfile1::getFieldGap() {
    return gapHeight_m;
}
void FM1DProfile1::setFieldGap(double gap) {

    gapHeight_m = gap;

}

double FM1DProfile1::computeEntranceFringe(double z) const {
    return computeFringe(engeCoeffsEntry_m, z / gapHeight_m);
}

double FM1DProfile1::computeExitFringe(double z) const {
    return computeFringe(engeCoeffsExit_m, z / gapHeight_m);
}

double FM1DProfile1::computeFringe(const std::vector<double> &coefs, double z) const {

    const size_t N = coefs.size();
    double expSum = coefs.at(0);

    for (size_t i = 1; i < N; ++ i) {
        expSum += std::pow(z, i) * coefs.at(i);
    }

    return 1.0 / (1.0 + exp(expSum));
}
