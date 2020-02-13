#include "Fields/FM1DProfile2.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DProfile2::FM1DProfile2(std::string aFilename)
    : Fieldmap(aFilename),
      EngeCoefs_entry_m(NULL),
      EngeCoefs_exit_m(NULL),
      exit_slope_m(0.0),
      xExit_m(0.0),
      zExit_m(0.0),
      cosExitRotation_m(1.0),
      sinExitRotation_m(0.0) {

    int tmpInt;
    std::string tmpString;
    double tmpDouble;
    int num_gridpz = -1;

    Type = T1DProfile2;
    ifstream file(Filename_m.c_str());

    if(file.good()) {
        bool parsing_passed =                               \
                interpreteLine<std::string, int, int, double>(file,
                        tmpString,
                        polynomialOrder_entry_m,
                        polynomialOrder_exit_m,
                        gapHeight_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, double, int>(file,
                                 zbegin_entry_m,
                                 polynomialOrigin_entry_m,
                                 zend_entry_m,
                                 num_gridpz);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, double, int>(file,
                                 zbegin_exit_m,
                                 polynomialOrigin_exit_m,
                                 zend_exit_m,
                                 tmpInt);
        for(int i = 0; (i < num_gridpz + 1) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_exit_m = zbegin_entry_m - 1e-3;
            zend_entry_m = zbegin_entry_m - 1e-3;
            zbegin_exit_m = zbegin_entry_m - 1e-3;
        } else {
            // conversion cm to m
            zbegin_entry_m /= 100.;
            zend_entry_m /= 100.;
            polynomialOrigin_entry_m /= 100.;
            zbegin_exit_m /= 100.;
            zend_exit_m /= 100.;
            polynomialOrigin_exit_m /= 100.;
            gapHeight_m /= 100.0;
        }
        length_m = zend_exit_m - zbegin_entry_m;
    } else {
        noFieldmapWarning();
        zbegin_entry_m = 0.0;
        zend_exit_m = zbegin_entry_m - 1e-3;
        zend_entry_m = zbegin_entry_m - 1e-3;
        zbegin_exit_m = zbegin_entry_m - 1e-3;
    }
}

FM1DProfile2::~FM1DProfile2() {
    if(EngeCoefs_entry_m != NULL) {
        delete[] EngeCoefs_entry_m;
        delete[] EngeCoefs_exit_m;
    }
}

void FM1DProfile2::readMap() {
    if(EngeCoefs_entry_m == NULL) {
        double tolerance = 1e-8;

        ifstream in(Filename_m.c_str());

        int tmpInt;
        std::string tmpString;
        double tmpDouble;

        double minValue = 99999999.99, maxValue = -99999999.99;
        int num_gridpz;
        int num_gridp_fringe_entry, num_gridp_fringe_exit;
        int num_gridp_before_fringe_entry, num_gridp_before_fringe_exit;
        double *RealValues;
        double *rightHandSide;
        double *leastSquareMatrix;
        double dZ;

        interpreteLine<std::string, int, int, double>(in, tmpString, tmpInt, tmpInt, tmpDouble);
        interpreteLine<double, double, double, int>(in, tmpDouble, tmpDouble, tmpDouble, num_gridpz);
        interpreteLine<double, double, double, int>(in, tmpDouble, tmpDouble, tmpDouble, tmpInt);

        dZ = length_m / num_gridpz;

        EngeCoefs_entry_m = new double[polynomialOrder_entry_m + 1];
        EngeCoefs_exit_m = new double[polynomialOrder_exit_m + 1];

        RealValues = new double[num_gridpz + 1];

        for(int i = 0; i < num_gridpz + 1; ++i) {
            interpreteLine<double>(in, RealValues[i]);
            if(RealValues[i] > maxValue) maxValue = RealValues[i];
            else if(RealValues[i] < minValue) minValue = RealValues[i];
        }
        in.close();

        // normalise the values //
        for(int i = 0; i < num_gridpz + 1; ++i)
            RealValues[i] = (RealValues[i] - minValue) / (maxValue - minValue);

        // find begin of entry fringe field //
        int i = 0;
        while(i < num_gridpz + 1 && RealValues[i] < tolerance) ++i;
        num_gridp_before_fringe_entry = i - 1;

        // find end of entry fringe field //
        while(i < num_gridpz + 1 && RealValues[i] < 1. - tolerance) ++i;
        num_gridp_fringe_entry = i - 1 - num_gridp_before_fringe_entry;

        // find begin of exit fringe field //
        while(i < num_gridpz + 1 && RealValues[i] >= 1. - tolerance) ++i;
        num_gridp_before_fringe_exit = i - 1;

        while(i < num_gridpz + 1 && RealValues[i] > tolerance) ++i;
        num_gridp_fringe_exit = i - 1 - num_gridp_before_fringe_exit;

        // set the origin of the polynomials

        int num_gridp_fringe = max(num_gridp_fringe_entry, num_gridp_fringe_exit);
        int polynomialOrder = max(polynomialOrder_entry_m, polynomialOrder_exit_m);
        leastSquareMatrix = new double[(polynomialOrder + 1) * num_gridp_fringe];
        rightHandSide = new double[num_gridp_fringe];


        double first = polynomialOrigin_entry_m - zbegin_entry_m - dZ * (num_gridp_before_fringe_entry + 1);
        for(int i = 0; i < num_gridp_fringe_entry; ++i) {
            double powerOfZ = 1.;
            double Z = (first - dZ * i) / gapHeight_m;
            rightHandSide[i] = log(1. / RealValues[num_gridp_before_fringe_entry + i + 1] - 1.);
            for(int j = 0; j < polynomialOrder_entry_m + 1; ++j) {
                leastSquareMatrix[i * (polynomialOrder_entry_m + 1) + j] = powerOfZ;
                powerOfZ *= Z;
            }
        }

        QRDecomposition::solve(leastSquareMatrix, EngeCoefs_entry_m, rightHandSide, num_gridp_fringe_entry, polynomialOrder_entry_m + 1);

        first = polynomialOrigin_exit_m - zbegin_entry_m - dZ * (num_gridp_before_fringe_exit + 1);
        for(int i = 0; i < num_gridp_fringe_exit; ++i) {
            double powerOfZ = 1.;
            double Z = (dZ * i - first) / gapHeight_m;
            rightHandSide[i] = log(1. / RealValues[num_gridp_before_fringe_exit + i + 1] - 1.);
            for(int j = 0; j < polynomialOrder_exit_m + 1; ++j) {
                leastSquareMatrix[i * (polynomialOrder_exit_m + 1) + j] = powerOfZ;
                powerOfZ *= Z;
            }
        }
        QRDecomposition::solve(leastSquareMatrix, EngeCoefs_exit_m, rightHandSide, num_gridp_fringe_exit, polynomialOrder_exit_m + 1);

        zbegin_exit_m = zbegin_entry_m + num_gridp_before_fringe_exit * dZ;
        zend_exit_m = zbegin_exit_m + num_gridp_fringe_exit * dZ;
        zbegin_entry_m += num_gridp_before_fringe_entry * dZ;
        zend_entry_m = zbegin_entry_m + num_gridp_fringe_entry * dZ;

        delete[] RealValues;
        delete[] leastSquareMatrix;
        delete[] rightHandSide;

        INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl);

    }
}

void FM1DProfile2::freeMap() {
    if(EngeCoefs_entry_m != NULL) {

        delete[] EngeCoefs_entry_m;
        delete[] EngeCoefs_exit_m;

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl);

    }
}

bool FM1DProfile2::getFieldstrength(const Vector_t &R, Vector_t &strength, Vector_t &info) const {

    info = Vector_t(0.0);

    // Find coordinates in the entrance frame.
    Vector_t REntrance(R(0), 0.0, R(2) + zbegin_entry_m);

    // Find coordinates in the exit frame.
    Vector_t RExit(0.0, R(1), 0.0);

    RExit(0) = (R(0) - xExit_m) * cosExitRotation_m - (R(2) + zbegin_entry_m - zExit_m) * sinExitRotation_m;
    RExit(2) = (R(0) - xExit_m) * sinExitRotation_m + (R(2) + zbegin_entry_m - zExit_m) * cosExitRotation_m + polynomialOrigin_exit_m;


    if(REntrance(2) >= zend_entry_m && RExit(2) <= zbegin_exit_m) {
        strength = Vector_t(1.0, 0.0, 0.0);
    } else {
        double S, dSdz, d2Sdz2 = 0.0;
        double expS, f, dfdz, d2fdz2;
        double z;
        double *EngeCoefs;
        int polynomialOrder;
        info(0) = 1.0;
        if(REntrance(2) >= zbegin_entry_m && REntrance(2) < zend_entry_m) {
            z = -(REntrance(2) - polynomialOrigin_entry_m) / gapHeight_m;
            EngeCoefs = EngeCoefs_entry_m;
            polynomialOrder = polynomialOrder_entry_m;
        } else if(RExit(2) > zbegin_exit_m && RExit(2) <= zend_exit_m) {
            z = (RExit(2) - polynomialOrigin_exit_m) / gapHeight_m;
            EngeCoefs = EngeCoefs_exit_m;
            polynomialOrder = polynomialOrder_exit_m;
            info(1) = 1.0;
        } else {
            return true;
        }

        S = EngeCoefs[polynomialOrder] * z;
        S += EngeCoefs[polynomialOrder - 1];
        dSdz = polynomialOrder * EngeCoefs[polynomialOrder];

        for(int i = polynomialOrder - 2; i >= 0; i--) {
            S = S * z + EngeCoefs[i];
            dSdz = dSdz * z + (i + 1) * EngeCoefs[i + 1];
            d2Sdz2 = d2Sdz2 * z + (i + 2) * (i + 1) * EngeCoefs[i + 2];
        }
        expS = exp(S);
        f = 1.0 / (1.0 + expS);
        if(f > 1.e-30) {
            // First derivative of Enge function, f.
            dfdz = - f * ((f * expS) * dSdz);

            // Second derivative of Enge functioin, f.
            d2fdz2 = ((-d2Sdz2 - dSdz * dSdz * (1. - 2. * (expS * f))) * (f * expS) * f) / (gapHeight_m * gapHeight_m);

            strength(0) = f;
            strength(1) = dfdz / gapHeight_m;
            strength(2) = d2fdz2;
        } else {
            strength = Vector_t(0.0);
        }

    }
    info(2) = exit_slope_m;

    return true;

    //    info = Vector_t(0.0);
    //    const Vector_t tmpR(R(0), R(1), R(2) + zbegin_entry_m);
    //
    //    if(tmpR(2) >= zend_entry_m && tmpR(2) <= exit_slope_m * tmpR(0) + zbegin_exit_m) {
    //        strength = Vector_t(1.0, 0.0, 0.0);
    //        info(0) = 3.0;
    //    } else {
    //        double S, dSdz, d2Sdz2 = 0.0;
    //        double expS, f, dfdz, d2fdz2;
    //        double z;
    //        double *EngeCoefs;
    //        int polynomialOrder;
    //        if(tmpR(2) >= zbegin_entry_m && tmpR(2) < zend_entry_m) {
    //            z = -(tmpR(2) - polynomialOrigin_entry_m) / gapHeight_m;
    //            EngeCoefs = EngeCoefs_entry_m;
    //            polynomialOrder = polynomialOrder_entry_m;
    //            info(0) = 1.0;
    //        } else if(tmpR(2) > exit_slope_m * tmpR(0) + zbegin_exit_m && tmpR(2) <= exit_slope_m * tmpR(0) + zend_exit_m) {
    //            z = (tmpR(2) - exit_slope_m * tmpR(0) - polynomialOrigin_exit_m) / sqrt(exit_slope_m * exit_slope_m + 1) / gapHeight_m;
    //            EngeCoefs = EngeCoefs_exit_m;
    //            polynomialOrder = polynomialOrder_exit_m;
    //            info(0) = 2.0;
    //        } else {
    //            return true;
    //        }
    //
    //        S = EngeCoefs[polynomialOrder] * z;
    //        S += EngeCoefs[polynomialOrder - 1];
    //        dSdz = polynomialOrder * EngeCoefs[polynomialOrder];
    //
    //        for(int i = polynomialOrder - 2; i >= 0; i--) {
    //            S = S * z + EngeCoefs[i];
    //            dSdz = dSdz * z + (i + 1) * EngeCoefs[i+1];
    //            d2Sdz2 = d2Sdz2 * z + (i + 2) * (i + 1) * EngeCoefs[i+2];
    //        }
    //        expS = exp(S);
    //        f = 1.0 / (1.0 + expS);
    //        if(f > 1.e-30) {
    //            dfdz = - f * ((f * expS) * dSdz); // first derivative of f
    //            d2fdz2 = ((-d2Sdz2 - dSdz * dSdz * (1. - 2. * (expS * f))) * (f * expS) * f) / (gapHeight_m * gapHeight_m);  // second derivative of f
    //
    //            strength(0) = f;
    //            strength(1) = dfdz / gapHeight_m;
    //            strength(2) = d2fdz2;
    //        } else {
    //            strength = Vector_t(0.0);
    //        }
    //
    //    }
    //    info(1) = exit_slope_m;
    //    return true;

}

bool FM1DProfile2::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM1DProfile2::getFieldDimensions(double &zBegin, double &zEnd, double &/*rBegin*/, double &/*rEnd*/) const {
    zBegin = zbegin_entry_m;
    zEnd = zend_exit_m;
}
void FM1DProfile2::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM1DProfile2::swap()
{}

void FM1DProfile2::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D Profile type 2); zini= " << zbegin_entry_m << " m; zfinal= " << zend_exit_m << " m;" << endl;
}

double FM1DProfile2::getFrequency() const {
    return 0.0;
}

void FM1DProfile2::setFrequency(double /*freq*/)
{}

void FM1DProfile2::setExitFaceSlope(const double &m) {
    exit_slope_m = m;
}

void FM1DProfile2::setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle) {

    double deltaZ = polynomialOrigin_exit_m - polynomialOrigin_entry_m;

    zExit_m = polynomialOrigin_entry_m + deltaZ * cos(bendAngle / 2.0);
    xExit_m = -deltaZ * sin(bendAngle / 2.0);

    cosExitRotation_m = cos(-bendAngle + entranceAngle + exitAngle);
    sinExitRotation_m = sin(-bendAngle + entranceAngle + exitAngle);
}

namespace QRDecomposition {
    void solve(double *Matrix, double *Solution, double *rightHandSide, const int &M, const int &N) {
        double sinphi;
        double cosphi;
        double tempValue;
        double len;
        double *R = new double[M * N];
        double *tempVector = new double[M];
        double *residuum = new double[M];

        for(int i = 0; i < M; ++i) {
            for(int j = 0; j < N; ++j)
                R[i * N + j] = Matrix[i * N + j];
            tempVector[i] = rightHandSide[i];
        }

        /* using Givens rotations */
        for(int i = 0; i < N; ++i) {
            for(int j = i + 1; j < M; ++j) {
                len = sqrt(R[j * N + i] * R[j * N + i] + R[i * (N + 1)] * R[i * (N + 1)]);
                sinphi = -R[j * N + i] / len;
                cosphi = R[i * (N + 1)] / len;

                for(int k = 0; k < N; ++k) {
                    tempValue = cosphi * R[ i * N + k] - sinphi * R[ j * N + k];
                    R[j * N + k] = sinphi * R[ i * N + k] + cosphi * R[ j * N + k];
                    R[i * N + k] = tempValue;
                }
            }
        }

        /* one step of iterative refinement */

        //     cout << "A^T*b" << endl;
        for(int i = 0; i < N; ++i) {     /* A^T*b */
            tempValue = 0.0;
            for(int j = 0; j < M; ++j) {
                tempValue += Matrix[j * N + i] * rightHandSide[j];
            }
            Solution[i] = tempValue;
        }
        //     cout << "R^-TA^T*b" << endl;
        for(int i = 0; i < N; ++i) {    /* R^-T*A^T*b */
            tempValue = 0.0;
            for(int j = 0; j < i; ++j)
                tempValue += R[j * N + i] * residuum[j];
            residuum[i] = (Solution[i] - tempValue) / R[i * (N + 1)];
        }
        //     cout << "R^-1R^-TA^T*b" << endl;
        for(int i = N - 1; i >= 0; --i) { /* R^-1*R^-T*A^T*b */
            tempValue = 0.0;
            for(int j = N - 1; j > i; --j)
                tempValue += R[i * N + j] * Solution[j];
            Solution[i] = (residuum[i] - tempValue) / R[i * (N + 1)];
        }
        //     cout << "b - A*x" << endl;
        for(int i = 0; i < M; ++i) {
            tempValue = 0.0;
            for(int j = 0; j < N; ++j)
                tempValue += Matrix[i * N + j] * Solution[j];
            residuum[i] = rightHandSide[i] - tempValue;
        }
        //     cout << "A^T*r" << endl;
        for(int i = 0; i < N; ++i) {
            tempValue = 0.0;
            for(int j = 0; j < M; ++j)
                tempValue += Matrix[j * N + i] * residuum[j];
            tempVector[i] = tempValue;
        }
        //     cout << "R^-TA^T*r" << endl;
        for(int i = 0; i < N; ++i) {
            tempValue = 0.0;
            for(int j = 0; j < i; ++j)
                tempValue += R[j * N + i] * residuum[j];
            residuum[i] = (tempVector[i] - tempValue) / R[i * (N + 1)];
        }
        //     cout << "R^-1R^-TA^T*r" << endl;
        for(int i = N - 1; i >= 0; --i) {
            tempValue = 0.0;
            for(int j = N - 1; j > i; --j)
                tempValue += R[i * N + j] * tempVector[j];
            tempVector[i] = (residuum[i] - tempValue) / R[i * (N + 1)];
            Solution[i] += tempVector[i];
        }

        //     for (int i = 0; i < N; ++i)
        //        cout << Solution[i] << endl;
        //     cout << endl;

        delete[] residuum;
        delete[] tempVector;
        delete[] R;

    }

}
