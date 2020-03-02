#include "gtest/gtest.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTBase.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"

#include "opal_test_utilities/SilenceTest.h"

#include <fstream>
#include <cmath>

using namespace std;

/*vector< vector<double> > partialsDerivB(const Vector_t &R,const Vector_t B, double stepSize, Component* dummyField)
{
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
    {
      // B at the previous and next grid points R_prev,  R_next
      Vector_t R_prev = R, R_next = R;
      R_prev[i] -= stepSize;
      R_next[i] += stepSize;
      Vector_t B_prev, B_next;
      dummyField->apply(R_prev, P, t, E, B_prev);
      dummyField->apply(R_next, P, t, E, B_next);
      for(int j = 0; j < 3; j++)
        allPartials[i][j] = (B_next[j] - B_prev[j]) / (2 * stepSize);
    }
     return allPartials;
}*/

vector< vector<double> > partialsDerivB(const Vector_t &R,const Vector_t /*B*/, double stepSize, Component* dummyField)
{
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
    {
      // B at the previous and next grid points R_prev,  R_next
      Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
      R_pprev(i) -= 2 * stepSize;
      R_nnext(i) += 2 * stepSize;
      R_prev(i) -= stepSize;
      R_next(i) += stepSize;
      Vector_t B_prev, B_next, B_pprev, B_nnext;
      dummyField->apply(R_prev, P, t, E, B_prev);
      dummyField->apply(R_next, P, t, E, B_next);
      dummyField->apply(R_pprev, P, t, E, B_pprev);
      dummyField->apply(R_nnext, P, t, E, B_nnext);
      for(int j = 0; j < 3; j++)
        allPartials[i][j] = (B_pprev[j] - 8 * B_prev[j] + 8 * B_next[j] - B_nnext[j]) / (12 * stepSize);
    }
    return allPartials;
}

double calcDivB(Vector_t &R, Vector_t B, double stepSize, Component* dummyField )
{
    double div = 0;
    vector< vector<double> > partials (3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    for(int i = 0; i < 3; i++)
        div += partials[i][i];
    return div;
}

vector<double> calcCurlB(Vector_t &R, Vector_t B, double stepSize, Component* dummyField)
{
    vector<double> curl(3);
    vector< vector<double> > partials(3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}

TEST(MultipoleTTest, Maxwell) {
    OpalTestUtilities::SilenceTest silencer;

    MultipoleT* myMagnet = new MultipoleT("Quadrupole");
    double centralField = 5;
    double fringeLength = 0.5;
    // the highest differential in the fringe field
    double max_index = 0.5;
    //Set the magnet
    myMagnet->setFringeField(centralField, fringeLength, max_index);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setDipoleConstant(1.0);
    myMagnet->setTransProfile(1, 100.0);
    //highest power in the field is z ^ (2 * maxOrder + 1)
    //      !!!  should be less than max_index / 2 !!!
    myMagnet->setMaxOrder(3);
    //ofstream fout("Quad_CurlB_off");
    Vector_t R(0., 0., 0.), P(3), E(3);
    double t = 0., x, z, s, stepSize= 1e-7;
    for(x = -0.2; x <= 0.2 ; x += 0.1) {
        for(z = 0.0; z <= 0.02 ; z += 0.001) {
            for(s = -10; s <= 10 ; s += 0.5) {
                R[0] = x;
                R[1] = z;
                R[2] = s;
                Vector_t B(0., 0., 0.);
                myMagnet->apply(R, P, t, E, B);
                double div = calcDivB(R, B, stepSize, myMagnet);
                vector<double> curl = calcCurlB(R, B, stepSize, myMagnet);
                EXPECT_NEAR(div, 0.0, 0.01);
                EXPECT_NEAR(curl[0], 0.0, 1e-4);
                EXPECT_NEAR(curl[1], 0.0, 1e-4);
                EXPECT_NEAR(curl[2], 0.0, 1e-4);
            }
        }
    }
    delete myMagnet;
}

TEST(MultipoleTTest, CurvedMagnet) {
    OpalTestUtilities::SilenceTest silencer;

    MultipoleT* myMagnet = new MultipoleT("Combined function");
    myMagnet->setLength(4.4);
    myMagnet->setBoundingBoxLength(0.0);
    myMagnet->setBendAngle(0.628);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setVarRadius();
    myMagnet->setTransMaxOrder(1);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile(0, 1);
    myMagnet->setTransProfile(1, 1);
    myMagnet->setMaxXOrder(3);
    myMagnet->setMaxOrder(3);
    double t = 0.0;
    double stepSize = 1e-3;
    double x[21] = {-1.12, -0.99, -0.86, -0.77, -0.65, -0.53, -0.42, -0.29, -0.19, -0.11, -0.039, 0.00, -0.030, -0.12, -0.26, -0.40, -0.56, -0.72, -0.86, -0.96, -1.12};
    double y[21] = {-2.74, -2.58, -2.30, -2.27, -2.00, -1.83, -1.62, -1.45, -1.13, -0.87, 0.53, 0.00, 0.46, 0.90, 1.36, 1.60, 1.83, 2.17, 2.30, 2.45, 2.77};
    double z = 0.2;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for (int n = 0; n < 21; n++) {
        R[0] = x[n];
        R[1] = z;
        R[2] = y[n];
        Vector_t B(0., 0., 0.);
        myMagnet->apply(R, P, t, E, B);
        double div = calcDivB(R, B, stepSize, myMagnet);
        vector<double> curl = calcCurlB(R, B, stepSize, myMagnet);
        double curlMag = 0.0;
        curlMag += gsl_sf_pow_int(curl[0], 2.0);
        curlMag += gsl_sf_pow_int(curl[1], 2.0);
        curlMag += gsl_sf_pow_int(curl[2], 2.0);
        curlMag = sqrt(curlMag);
        coordinatetransform::CoordinateTransform t(x[n], z, y[n], 2.2, 0.3, 0.3, 4.4 / 0.628);
        std::vector<double> r = t.getTransformation();
        EXPECT_NEAR(div, 0, 2e-2)
                     << "R: " << r[0] << " " << r[1] << " " << r[2] << std::endl
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
        EXPECT_NEAR(curlMag, 0, 1e-9)
                     << "R: " << r[0] << " " << r[1] << " " << r[2] << std::endl
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
    }
    delete myMagnet;
}

TEST(MultipoleTTest, Straight) {
    // failing
    OpalTestUtilities::SilenceTest silencer;

    MultipoleTStraight* myMagnet = new MultipoleTStraight("Combined function");
    myMagnet->setLength(4.4);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile(0, 1);
    myMagnet->setTransProfile(1, 1);
    myMagnet->setMaxOrder(5);
    double t = 0.0;
    double stepSize = 1e-3;
    double z = -0.3;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for (double x = -0.3; x <= 0.300001; x += 0.1) {
        for (double y = -3.0; y <= 3.00001; y += 1.0) {
            R[0] = x;
            R[1] = z;
            R[2] = y;
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet);
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet);
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
        }
    }
    delete myMagnet;
}

TEST(MultipoleTTest, CurvedConstRadius) {
    OpalTestUtilities::SilenceTest silencer;

    MultipoleTCurvedConstRadius* myMagnet = new MultipoleTCurvedConstRadius("Combined function");
    myMagnet->setLength(4.4);
    myMagnet->setBendAngle(0.628);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile(0, 1);
    myMagnet->setTransProfile(1, 1);
    myMagnet->setMaxXOrder(20);
    myMagnet->setMaxOrder(5);
    double t = 0.0;
    double stepSize = 1e-3;
    double radius = 4.4 / 0.628;
    double z = 0.2;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for (double theta = 0; theta <= 0.3001; theta += 0.2) {
        double x = radius * cos(theta) - radius;
        double y = radius * sin(theta);
        for (double delta = -0.3; delta <= 0.3001; delta += 0.02) {
            Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
            R[0] = x + delta * cos(theta);
            R[1] = z;
            R[2] = y + delta * sin(theta);
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet);
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet);
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 5e-6)
                     << "R: " << delta << " " << z << " " << radius * theta << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-9)
                     << "R: " << delta << " " << z << " " << radius * theta << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
        }
    }
    delete myMagnet;
}

TEST(MultipoleTTest, CurvedVarRadius) {
    OpalTestUtilities::SilenceTest silencer;

    MultipoleTCurvedVarRadius* myMagnet = new MultipoleTCurvedVarRadius("Combined function");
    myMagnet->setLength(4.4);
    myMagnet->setBendAngle(0.628);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setTransMaxOrder(1);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile(0, 1);
    myMagnet->setTransProfile(1, 1);
    myMagnet->setMaxXOrder(3);
    myMagnet->setMaxOrder(3);
    double t = 0.0;
    double stepSize = 1e-3;
    double x[21] = {-1.12, -0.99, -0.86, -0.77, -0.65, -0.53, -0.42, -0.29, -0.19, -0.11, -0.039, 0.00, -0.030, -0.12, -0.26, -0.40, -0.56, -0.72, -0.86, -0.96, -1.12};
    double y[21] = {-2.74, -2.58, -2.30, -2.27, -2.00, -1.83, -1.62, -1.45, -1.13, -0.87, 0.53, 0.00, 0.46, 0.90, 1.36, 1.60, 1.83, 2.17, 2.30, 2.45, 2.77};
    double z = 0.2;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for (int n = 0; n < 21; n++) {
        R[0] = x[n];
        R[1] = z;
        R[2] = y[n];
        Vector_t B(0., 0., 0.);
        myMagnet->apply(R, P, t, E, B);
        double div = calcDivB(R, B, stepSize, myMagnet);
        vector<double> curl = calcCurlB(R, B, stepSize, myMagnet);
        double curlMag = 0.0;
        curlMag += gsl_sf_pow_int(curl[0], 2.0);
        curlMag += gsl_sf_pow_int(curl[1], 2.0);
        curlMag += gsl_sf_pow_int(curl[2], 2.0);
        curlMag = sqrt(curlMag);
        EXPECT_NEAR(div, 0, 2e-2)
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
        EXPECT_NEAR(curlMag, 0, 1e-9)
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
    }
    delete myMagnet;
}