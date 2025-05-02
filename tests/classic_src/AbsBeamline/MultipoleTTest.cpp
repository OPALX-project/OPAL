#include "gtest/gtest.h"
#include "AbsBeamline/Component.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"
#include "AbsBeamline/EndFieldModel/Tanh.h"
#include <gsl/gsl_sf_pow_int.h>
#include "opal_test_utilities/SilenceTest.h"

#include <cmath>
#include <sstream>

using namespace std;

std::ofstream* fout = nullptr;

Vector_t rotateBy(const Vector_t& center, const Vector_t& point, double theta) {
    auto x_prime = point[0] - center[0];
    auto z_prime = point[2] - center[2];
    auto x_rotated = x_prime * cos(theta) - z_prime * sin(theta) + center[0];
    auto z_rotated = x_prime * sin(theta) + z_prime * cos(theta) + center[2];
    return {x_rotated, point[1], z_rotated};
}

vector< vector<double> > partialsDerivB(const Vector_t &R,const Vector_t /*B*/, double stepSize,
                                       Component* dummyField, double theta = 0.0)
{
    auto* magnet = dynamic_cast<MultipoleT*>(dummyField);
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++) {
        // B at the previous and next grid points R_prev,  R_next
        Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
        R_pprev(i) -= 2 * stepSize;
        R_nnext(i) += 2 * stepSize;
        R_prev(i) -= stepSize;
        R_next(i) += stepSize;
        // Rotate these points by the angle
        R_pprev = rotateBy(R, R_pprev, theta);
        R_nnext = rotateBy(R, R_nnext, theta);
        R_prev = rotateBy(R, R_prev, theta);
        R_next = rotateBy(R, R_next, theta);
        if(fout) {
            *fout << "R=" << magnet->toMagnetCoords(R)
                  << ", R_pprev=" << magnet->toMagnetCoords(R_pprev)
                  << ", R_nnext=" << magnet->toMagnetCoords(R_nnext)
                  << ", R_prev=" << magnet->toMagnetCoords(R_prev)
                  << ", R_next=" << magnet->toMagnetCoords(R_next) << std::endl;
        }
        // Get the magnetic fields and derivatives
        Vector_t B_prev, B_next, B_pprev, B_nnext;
        dummyField->apply(R_prev, P, t, E, B_prev);
        dummyField->apply(R_next, P, t, E, B_next);
        dummyField->apply(R_pprev, P, t, E, B_pprev);
        dummyField->apply(R_nnext, P, t, E, B_nnext);
        for(int j = 0; j < 3; j++) {
            allPartials[i][j] =
                    (B_pprev[j] - 8 * B_prev[j] + 8 * B_next[j] - B_nnext[j]) / (12 * stepSize);
        }
    }
    return allPartials;
}

double calcDivB(Vector_t &R, Vector_t B, double stepSize, Component* dummyField,
                double theta = 0.0) {
    double div = 0;
    vector< vector<double> > partials (3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField, theta);
    for(int i = 0; i < 3; i++) {
        div += partials[i][i];
    }
    return div;
}

vector<double> calcCurlB(Vector_t &R, Vector_t B, double stepSize, Component* dummyField,
                         double theta = 0.0) {
    vector<double> curl(3);
    vector< vector<double> > partials(3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField, theta);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}

TEST(MultipoleTTest, Maxwell) {
    OpalTestUtilities::SilenceTest silencer;
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    myMagnet->setBendAngle(0.0, false);
    double centralField = 5;
    double fringeLength = 0.5;
    // the highest differential in the fringe field
    double max_index = 0.5;
    //Set the magnet
    myMagnet->setElementLength(1.0);
    myMagnet->setFringeField(centralField, fringeLength, max_index);
    myMagnet->setTransProfile({1.0, 100.0});
    //highest power in the field is z ^ (2 * maxOrder + 1)
    //      !!!  should be less than max_index / 2 !!!
    myMagnet->setMaxOrder(3, 10);
    //ofstream fout("Quad_CurlB_off");
    Vector_t R(0., 0., 0.), P(3), E(3);
    double t = 0., stepSize= 1e-7;
    for(size_t i = 0; i < 5; ++i) {
        double x = -0.2 + static_cast<double>(i) * 0.1;
        for(size_t j = 0; j < 20; ++j) {
            double z = static_cast<double>(j) * 0.001;
            for(size_t k = 0; k < 41; ++k) {
                double s = -10.0 + static_cast<double>(k) * 0.5;
                R[0] = x;
                R[1] = z;
                R[2] = s - 0.5;
                Vector_t B(0., 0., 0.);
                myMagnet->apply(R, P, t, E, B);
                double div = calcDivB(R, B, stepSize, myMagnet.get());
                vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get());
                EXPECT_NEAR(div, 0.0, 0.01);
                EXPECT_NEAR(curl[0], 0.0, 1e-4);
                EXPECT_NEAR(curl[1], 0.0, 1e-4);
                EXPECT_NEAR(curl[2], 0.0, 1e-4);
            }
        }
    }
}

TEST(MultipoleTTest, CurvedMagnet) {
    OpalTestUtilities::SilenceTest silencer;
    // Build the magnet
    double angle = 0.628;
    double length = 4.4;
    double rho = length / angle;
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    myMagnet->setBendAngle(angle, true);
    myMagnet->setElementLength(length);
    myMagnet->setBoundingBoxLength(0.0);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(length / 2.0, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(3, 3);
    // Account for the magnet origin being at the center
    //double rho = 4.4 / 0.628;
    //double dz = rho * sin(0.628 / 2.0);
    //double dx = rho * (1 - cos(0.628 / 2.0));
    // Test the magnet
    double t = 0.0;
    double stepSize = 1e-3;
    std::vector<double> x = {-1.12, -0.99, -0.86, -0.77,  -0.65,  -0.53, -0.42,
                             -0.29, -0.19, -0.11, -0.039, -0.030, -0.12, -0.26,
                             -0.40, -0.56, -0.72, -0.86,  -0.96,  -1.12};
    std::vector<double> y = {-2.74, -2.58, -2.30, -2.27, -2.00, -1.83, -1.62,
                             -1.45, -1.13, -0.87,  0.53,  0.46,  0.90,  1.36,
                              1.60,  1.83,  2.17,  2.30,  2.45,  2.77};
    double z = 0.2;
    Vector_t centerR{-rho, z, 0.0};
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    double localTheta = myMagnet->localCartesianRotation();
    for (size_t n = 0; n < x.size() && n < y.size(); n++) {
        R = myMagnet->localCartesianToOpalCartesian({x[n], z, -y[n]});
        Vector_t B(0., 0., 0.);
        myMagnet->apply(R, P, t, E, B);
        double div = calcDivB(R, B, stepSize, myMagnet.get(), localTheta);
        vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get(), localTheta);
        double curlMag = 0.0;
        curlMag += gsl_sf_pow_int(curl[0], 2.0);
        curlMag += gsl_sf_pow_int(curl[1], 2.0);
        curlMag += gsl_sf_pow_int(curl[2], 2.0);
        curlMag = sqrt(curlMag);
        coordinatetransform::CoordinateTransform ct(x[n], z, y[n], 2.2, 0.3, 0.3, 4.4 / 0.628);
        std::vector<double> r = ct.getTransformation();
        std::stringstream msgStr;
        msgStr << "R: " << r[0] << " " << r[1] << " " << r[2] << std::endl
               << "R: " << x[n] << " " << z << " " << y[n] << std::endl
               << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
               << "Del: " << div << " " << curl[0] << " " << curl[1] << " " << curl[2] << std::endl;
        std::cerr << msgStr.str() << std::endl;
        EXPECT_NEAR(div, 0, 2e-2) << msgStr.str();
        EXPECT_NEAR(curlMag, 0, 1e-9) << msgStr.str();
    }
}

TEST(MultipoleTTest, Straight) {
    // failing
    OpalTestUtilities::SilenceTest silencer;
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    myMagnet->setBendAngle(0.0, false);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(5, 20);
    double t = 0.0;
    double stepSize = 1e-3;
    double z = -0.3;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for(size_t i = 0; i < 6; ++i) {
        double x = -0.3 + static_cast<double>(i) * 0.1;
        for(size_t j = 0; j < 6; ++j) {
            double y = -3.0 + static_cast<double>(j) * 1.0;
            R[0] = x;
            R[1] = z;
            R[2] = y - length / 2.0;
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet.get());
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get());
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1] << " "
                << curl[2] << std::endl;
        }
    }
}

TEST(MultipoleTTest, ClonedStraight) {
    // failing
    OpalTestUtilities::SilenceTest silencer;
    // Build the magnet
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    myMagnet->setBendAngle(0.0, false);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(2.2, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(5, 20);
    // Make the clone
    myMagnet.reset(dynamic_cast<MultipoleT*>(myMagnet->clone()));
    // Test it
    double t = 0.0;
    double stepSize = 1e-3;
    double z = -0.3;
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    for(size_t i = 0; i < 6; ++i) {
        double x = -0.3 + static_cast<double>(i) * 0.1;
        for(size_t j = 0; j < 6; ++j) {
            double y = -3.0 + static_cast<double>(j) * 1.0;
            R[0] = x;
            R[1] = z;
            R[2] = y - length / 2.0;
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet.get());
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get());
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-1)
                << "R: " << x << " " << z << " " << y << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1] << " "
                << curl[2] << std::endl;
        }
    }
}

TEST(MultipoleTTest, CurvedConstRadius) {
    OpalTestUtilities::SilenceTest silencer;
    // Build the magnet
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    double angle = 0.628;
    myMagnet->setBendAngle(angle, false);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(length / 2.0, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(5, 20);
    // Test it
    double t = 0.0;
    double stepSize = 1e-3;
    double radius = length / angle;
    double z = 0.2;
    double dTheta = angle / 2.0;  // Account for the magnet origin originally at the center
    for(size_t i = 0; i < 2; ++i) {
        double theta = 0.0 + static_cast<double>(i) * 0.2;
        double x = radius * cos(theta + dTheta) - radius;
        double y = radius * sin(theta + dTheta);
        for(size_t k = 0; k < 31; ++k) {
            double delta = -0.3 + static_cast<double>(k) * 0.02;
            Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
            R[0] = x + delta * cos(theta + dTheta);
            R[1] = z;
            R[2] = y + delta * sin(theta + dTheta);
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet.get(), dTheta);
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get(), dTheta);
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 5e-6)
                << "R: " << delta << " " << z << " " << radius * theta << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-9)
                << "R: " << delta << " " << z << " " << radius * theta << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
        }
    }
}

TEST(MultipoleTTest, ClonedCurvedConstRadius) {
    OpalTestUtilities::SilenceTest silencer;
    // Build the magnet
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    double angle = 0.628;
    myMagnet->setBendAngle(angle, false);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(length / 2.0, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(5, 20);
    // Make the clone
    myMagnet.reset(dynamic_cast<MultipoleT*>(myMagnet->clone()));
    // Test it
    double t = 0.0;
    double stepSize = 1e-3;
    double radius = length / angle;
    double z = 0.2;
    double dTheta = angle / 2.0;  // Account for the magnet origin originally at the center
    for(size_t i = 0; i < 2; ++i) {
        double theta = static_cast<double>(i) * 0.1;
        double x = radius * cos(theta + dTheta) - radius;
        double y = radius * sin(theta + dTheta);
        for(size_t k = 0; k < 31; ++k) {
            double delta = -0.3 + static_cast<double>(k) * 0.02;
            Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
            R[0] = x + delta * cos(theta + dTheta);
            R[1] = z;
            R[2] = y + delta * sin(theta + dTheta);
            Vector_t B(0., 0., 0.);
            myMagnet->apply(R, P, t, E, B);
            double div = calcDivB(R, B, stepSize, myMagnet.get(), dTheta);
            vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get(), dTheta);
            double curlMag = 0.0;
            curlMag += gsl_sf_pow_int(curl[0], 2.0);
            curlMag += gsl_sf_pow_int(curl[1], 2.0);
            curlMag += gsl_sf_pow_int(curl[2], 2.0);
            curlMag = sqrt(curlMag);
            EXPECT_NEAR(div, 0, 5e-6)
                << "R: " << delta << " " << z << " " << radius * theta << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
            EXPECT_NEAR(curlMag, 0, 1e-9)
                << "R: " << delta << " " << z << " " << radius * theta << std::endl
                << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                << "Del: " << div << " " << curl[0] << " " << curl[1]
                << " " << curl[2] << std::endl;
        }
    }
}

TEST(MultipoleTTest, CurvedVarRadius) {
    fout = new std::ofstream("TestData.txt");
    OpalTestUtilities::SilenceTest silencer;
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    double angle = 0.628;
    double rho = length / angle;
    myMagnet->setBendAngle(angle, true);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(length / 2.0, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(3, 3);
    double t = 0.0;
    double stepSize = 1e-3;
    double x[21] = {-1.12, -0.99, -0.86, -0.77, -0.65, -0.53, -0.42, -0.29, -0.19, -0.11, -0.039,
                       0.00, -0.030, -0.12, -0.26, -0.40, -0.56, -0.72, -0.86, -0.96, -1.12};
    double y[21] = {-2.74, -2.58, -2.30, -2.27, -2.00, -1.83, -1.62, -1.45, -1.13, -0.87, 0.53,
                       0.00, 0.46, 0.90, 1.36, 1.60, 1.83, 2.17, 2.30, 2.45, 2.77};
    double z = 0.2;
    Vector_t centerR{-rho, z, 0.0};
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    double localTheta = myMagnet->localCartesianRotation();
    for (int n = 0; n < 21; n++) {
        R = myMagnet->localCartesianToOpalCartesian({x[n], z, -y[n]});
        Vector_t B(0., 0., 0.);
        myMagnet->apply(R, P, t, E, B);
        double div = calcDivB(R, B, stepSize, myMagnet.get(), localTheta);
        vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get(), localTheta);
        double curlMag = 0.0;
        curlMag += gsl_sf_pow_int(curl[0], 2.0);
        curlMag += gsl_sf_pow_int(curl[1], 2.0);
        curlMag += gsl_sf_pow_int(curl[2], 2.0);
        curlMag = sqrt(curlMag);
        EXPECT_NEAR(div, 0, 4e-2 /*2e-2*/)
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1]
                     << " " << curl[2] << std::endl;
        EXPECT_NEAR(curlMag, 0, 0.1 /*1e-9*/)
                     << "R: " << x[n] << " " << z << " " << y[n] << std::endl
                     << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
                     << "Del: " << div << " " << curl[0] << " " << curl[1] << " "
                     << curl[2] << std::endl;
    }
}

TEST(MultipoleTTest, ClonedCurvedVarRadius) {
    OpalTestUtilities::SilenceTest silencer;
    // Build a magnet
    auto myMagnet = std::make_unique<MultipoleT>("Combined function");
    double length = 4.4;
    double angle = 0.628;
    double rho = length / angle;
    myMagnet->setBendAngle(angle, true);
    myMagnet->setElementLength(length);
    myMagnet->setAperture(3.5, 3.5);
    myMagnet->setFringeField(length / 2.0, 0.3, 0.3);
    myMagnet->setRotation(0.0);
    myMagnet->setEntranceAngle(0.0);
    myMagnet->setTransProfile({1.0, 1.0});
    myMagnet->setMaxOrder(3, 3);
    // Make the clone
    myMagnet.reset(dynamic_cast<MultipoleT*>(myMagnet->clone()));
    // Test it
    double t = 0.0;
    double stepSize = 1e-3;
    double x[21] = {-1.12, -0.99, -0.86, -0.77, -0.65, -0.53, -0.42, -0.29, -0.19, -0.11, -0.039,
                       0.00, -0.030, -0.12, -0.26, -0.40, -0.56, -0.72, -0.86, -0.96, -1.12};
    double y[21] = {-2.74, -2.58, -2.30, -2.27, -2.00, -1.83, -1.62, -1.45, -1.13, -0.87, 0.53,
                       0.00, 0.46, 0.90, 1.36, 1.60, 1.83, 2.17, 2.30, 2.45, 2.77};
    double z = 0.2;
    Vector_t centerR{-rho, z, 0.0};
    Vector_t R(0.0, 0.0, 0.0), P(3), E(3);
    double localTheta = myMagnet->localCartesianRotation();
    for (int n = 0; n < 21; n++) {
        R = myMagnet->localCartesianToOpalCartesian({x[n], z, -y[n]});
        Vector_t B(0., 0., 0.);
        myMagnet->apply(R, P, t, E, B);
        double div = calcDivB(R, B, stepSize, myMagnet.get(), localTheta);
        vector<double> curl = calcCurlB(R, B, stepSize, myMagnet.get(), localTheta);
        double curlMag = 0.0;
        curlMag += gsl_sf_pow_int(curl[0], 2.0);
        curlMag += gsl_sf_pow_int(curl[1], 2.0);
        curlMag += gsl_sf_pow_int(curl[2], 2.0);
        curlMag = sqrt(curlMag);
        EXPECT_NEAR(div, 0, 2e-2)
            << "R: " << x[n] << " " << z << " " << y[n] << std::endl
            << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
            << "Del: " << div << " " << curl[0] << " " << curl[1]
            << " " << curl[2] << std::endl;
        EXPECT_NEAR(curlMag, 0, 1e-9)
            << "R: " << x[n] << " " << z << " " << y[n] << std::endl
            << "B: " << B[0] << " " << B[1] << " " << B[2] << std::endl
            << "Del: " << div << " " << curl[0] << " " << curl[1] << " "
            << curl[2] << std::endl;
    }
}

