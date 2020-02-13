#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Solvers/GreenWakeFunction.hh"

#include "TestLambda.h" // Profile of the bunch

#include <iostream>

TEST(GreenWakeFunctionTest, TestApply)
{
    OpalTestUtilities::SilenceTest silencer;

    // wake types with green wake function
    std::vector<std::string> types      = {"LONG-SHORT-RANGE", "TRANSV-SHORT-RANGE"};
    std::vector<int> acmodes = {1, 2}; // 1: ac, 2: dc
    std::vector<Filter *> filters;
    // not sure if values are realistic, but couldn't find code example
    double Z0     = 1e3;
    int nbin      = 100;
    double radius = 0.1;
    double sigma  = 1;
    double tau    = 1;
    bool const_length = true;
    std::string fname = "";

    std::vector<double> finalWakeValues   = {1.61757e+10,  2.78525e+19};
    std::vector<double> finalEnergyValues = {-0.00125303, -2.15739e+06};
    std::vector<double> relativeErrorWake   = {1e+6, 1e+15};
    std::vector<double> relativeErrorEnergy = {1e-7, 1e+2};

    for (int acmode : acmodes) {
        GreenWakeFunction gwf("opal", nullptr, filters, nbin, Z0, radius, sigma, acmode, tau, 0, const_length, fname);

        double spacing = 1e-6; //IFF: charge in testLambda.h in 1um spacings
        // determine K and charge
        double charge = 0.8e-9; // nC
        double K = 0.20536314319923724e-9; //K normalizes nC data in lambda.h?
        gwf.NBin_m = 294;

        std::cout << "# Z0 = "        << gwf.Z0_m        << std::endl
                  << "# radius = "    << gwf.radius_m    << std::endl
                  << "# sigma = "     << gwf.sigma_m     << std::endl
                  << "# acMode = "    << gwf.acMode_m    << std::endl
                  << "# tau = "       << gwf.tau_m       << std::endl
                  << "# direction = " << gwf.direction_m << std::endl
                  << "# spacing = "   << spacing         << std::endl
                  << "# Lbunch = "    << gwf.NBin_m      << std::endl;

        if(gwf.FftWField_m.empty()) {
            gwf.FftWField_m.resize(2*gwf.NBin_m - 1);
            gwf.CalcWakeFFT(spacing);
        } else if(!gwf.constLength_m) {
            gwf.CalcWakeFFT(spacing);
        }

        std::cout << "# Wake calculated in Opal" << std::endl;
        for(int i = 0; i < gwf.NBin_m; i++) {
            std::cout << i + 1 << "   " << gwf.FftWField_m[i] << std::endl;
        }

        std::vector<double> OutEnergy(gwf.NBin_m);
        gwf.compEnergy(K, charge, testLambda, OutEnergy.data());

        std::cout << "# Energy of the Wake calculated in Opal\n";

        for(int i = 0; i < gwf.NBin_m; i++) {
            std::cout << i + 1 << " " << OutEnergy[i] << std::endl;
        }

        EXPECT_NEAR(gwf.FftWField_m[gwf.NBin_m-1],   finalWakeValues[acmode-1], relativeErrorWake[acmode-1]);
        EXPECT_NEAR(      OutEnergy[gwf.NBin_m-1], finalEnergyValues[acmode-1], relativeErrorEnergy[acmode-1]);
    }
}