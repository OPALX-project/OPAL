//
// Class GreenWakeFunction
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef GREENWAKEFUNCTION_HH
#define GREENWAKEFUNCTION_HH

#include "Filters/Filter.h"
#include "Physics/Physics.h"
#include "Solvers/WakeFunction.h"
#include "Utility/IpplInfo.h"
#include "Utility/PAssert.h"

#include <complex>
#include <map>
#include <string>
#include <vector>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

enum class WakeDirection: unsigned short {
    TRANSVERSAL,
    LONGITUDINAL
};

typedef std::map<std::string, int> FilterOptions;

class GreenWakeFunction: public WakeFunction {
public:
    ~GreenWakeFunction();
    //IFF: changed direction to int (was double)
    //IFF: changed acMode to int (was double)
    GreenWakeFunction(const std::string& name,
                      std::vector<Filter*> filters,
                      int NBIN,
                      double Z0,
                      double radius,
                      double sigma,
                      int acMode,
                      double tau,
                      WakeDirection direction,
                      bool constLength,
                      std::string fname);

    std::pair<int, int> distrIndices(int vectLen);

    void apply(PartBunchBase<double, 3>* bunch) override;
    void setWakeFromFile(int NBin, double spacing);
    virtual WakeType getType() const override;

private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GreenWakeFunctionTest, TestApply);
#endif

    class Wake {

    public:
        Wake(double s, double Z0, double a, double sigma, int acMode, double tau, WakeDirection direction)
            : Z0_(Z0), a_(a), sigma_(sigma), s_(s), acMode_(acMode), tau_(tau), direction_(direction)
        {}

        /**
         * @brief   Used to integrate the function
         *
         * @param[in]   k parameter
         *
         * @return  the function value at position k
         */
        double operator()(double k) {

            std::complex <double> i(0, 1);
            std::complex <double> Z(0, 0);
            double signK = (k > 0 ? 1 : -1);

            //1 == AC
            //2 == DC
            switch(acMode_) {
                case 1:
                    Z = (Z0_ / (2 * Physics::pi * a_)) * 1.0 / (std::sqrt(Z0_ * std::abs(k) / 2) * std::sqrt(sigma_ / (1.0 - i * Physics::c * k * tau_)) * (i + signK) / k - (i * k * a_) / 2.0);
                    break;
                case 2:
                    Z = (Z0_ / (2 * Physics::pi * a_)) * 1.0 / (std::sqrt(sigma_ * Z0_ * std::abs(k) / 2) * (i + signK) / k - (i * k * a_) / 2.0);
                    break;
            }
            switch(direction_) {
                case WakeDirection::LONGITUDINAL:
                    return real(Z) * std::cos(k * s_) * 2.0 * Physics::c / Physics::pi;
                    break;
                case WakeDirection::TRANSVERSAL:
                    return real(Z) * Physics::c / k * std::cos(k * s_) * 2.0 * Physics::c / Physics::pi;
                    break;
            }
            ERRORMSG("We should not be here: " << __FILE__ << " L" << __LINE__ << endl);

            return 0.0;
        }

    private:
        /// impedance
        double Z0_;
        /// radius
        double a_;
        /// material constant
        double sigma_;
        /// distance from the particle
        double s_;
        /// conductivity either 1="AC" or 2="DC"
        int acMode_;
        /// material constant
        double tau_;
        /// direction either 1="Longitudinal" 0= "Transversal"
        WakeDirection direction_;
    };

    /**
     * @brief   Simpson-Integration from the function f from a to b with N steps
     *
     *
     * @param[in]   f the function to integrate
     * @param[in]   a integrate from a
     * @param[in]   b integrate to b
     * @param[in]   N Number of integration points
     * @return  function value of the integration
     *
     */
    template<class F> double simpson(F& f, double a, double b, unsigned int N) {
        PAssert(b > a);
        PAssert(N > 0);

        double result = 0;
        double h = (b - a) / N;

        // boundary values
        result += (f(a) + 4 * f(a + h / 2) + f(b)) / 2.0;

        // values between boundaries
        for(unsigned int i = 1; i < N; ++ i) {
            result += f(a + i * h) + 2 * f(a + (i + 0.5) * h);
        }

        result *= h / 3.0;

        return result;
    }

    /// save the line Density of the particle bunch
    std::vector<double> lineDensity_m;
    /// FFT of the zero padded wakefield
    std::vector<double>  FftWField_m;

    /// divides the particle bunch in NBin slices
    int NBin_m;
    /// impedance
    double Z0_m;
    /// radius
    double radius_m;
    /// material constant
    double sigma_m;
    /// conductivity either 1="AC" or 2="DC"
    int acMode_m;
    /// material constant
    double tau_m;
    /// direction either "Longitudinal" - "Transversal"
    WakeDirection direction_m;
    /// true if the length of the particle bunch is considered as constant
    bool constLength_m;
    /// filename of the wakefield
    std::string filename_m;

    std::vector<Filter*> filters_m;

    static const std::map<WakeDirection, std::string> wakeDirectiontoString_s;

    void compEnergy(const double K, const double charge, const double* lambda, double* OutEnergy);
    void compEnergy(const double K, const double charge, std::vector<double> lambda, double* OutEnergy);
    void CalcWakeFFT(double spacing);
    static std::string getWakeDirectionString(const WakeDirection& direction);
};
#endif //GREENWAKEFUNCTION_HH
