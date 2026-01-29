#ifndef ROOTFINDERFORCSR_H
#define ROOTFINDERFORCSR_H

#include "Utilities/GeneralClassicException.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <complex>
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>

class RootFinderForCSR {
public:
    RootFinderForCSR(double a,
                     double b,
                     double d,
                     double e):
        a_m(a),
        b_m(b),
        d_m(d),
        e_m(e),
        searchRange_m(0.0, 0.0) {
        deltaZero_m.real(std::pow(3. * b, 2));
        deltaZero_m.imag(0.0);
        deltaOne_m.real(2. * std::pow(3. * b, 3) + 27. * std::pow(4. * a, 2) * d);
        deltaOne_m.imag(0.0);
    }

    bool hasPositiveRealRoots() {
        std::complex<double> tmp = sqrt(std::pow(deltaOne_m, 2) - 4.0 * std::pow(deltaZero_m, 3));
        std::complex<double> C1 = std::pow(0.5 * (deltaOne_m + copysign(1.0, deltaOne_m.real()) * tmp), 1.0 / 3.0);

        std::complex<double> x1 = -(3.0 * b_m + C1 + deltaZero_m / C1) / (12. * a_m);
        if (std::abs(x1.imag()) < 1e-9 && x1.real() > 0.0)
            thirdOrderRoots_m.push_back(x1.real());

        std::complex<double> C2 = C1 * std::complex<double>(-0.5, -0.5 * sqrt(3));
        std::complex<double> x2 = -(3. * b_m + C2 + deltaZero_m / C2) / (12. * a_m);
        if (std::abs(x2.imag()) < 1e-9 && x2.real() > 0.0)
            thirdOrderRoots_m.push_back(x2.real());

        std::complex<double> C3 = C1 * std::complex<double>(-0.5, 0.5 * sqrt(3));
        std::complex<double> x3 = -(3. * b_m + C3 + deltaZero_m / C3) / (12. * a_m);
        if (std::abs(x3.imag()) < 1e-9 && x3.real() > 0.0)
            thirdOrderRoots_m.push_back(x3.real());

        if (thirdOrderRoots_m.size() == 0) return false;

        std::sort(thirdOrderRoots_m.begin(),
                  thirdOrderRoots_m.end());

        thirdOrderRoots_m.insert(thirdOrderRoots_m.begin(), 0.0);

        double rangeMax = thirdOrderRoots_m.back() + 0.1;
        while (computeValue(rangeMax) < 0.0) {
            rangeMax += 0.1;
        }
        thirdOrderRoots_m.push_back(rangeMax);

        double oldValue = computeValue(0.0);
        unsigned int size = thirdOrderRoots_m.size();
        unsigned int i;
        for (i = 1; i < size; ++ i) {
            const double &x = thirdOrderRoots_m[i];
            double value = computeValue(x);

            if (oldValue * value < 0.0) {
                searchRange_m.first = thirdOrderRoots_m[i - 1];
                searchRange_m.second = x;

                return true;
            }

            oldValue = value;
        }

        return false;
    }

    template <class T>
    T computeValue(const T &x) const {
        T xcube = std::pow(x, 3);

        return a_m * x * xcube + b_m * xcube + d_m * x + e_m;
    }

        template <class T>
    T computeDerivative(T x) const {
        T xsqr = std::pow(x, 2);

        return 4.0 * a_m * x * xsqr + 3.0 * b_m * xsqr + d_m;
    }

    double searchRoot(const double &tol) {
        int status;
        int iter = 0, max_iter = 100;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *solver;
        double root = 0;
        double x_lo = searchRange_m.first, x_hi = searchRange_m.second;
        gsl_function F;
        struct PolyParams params = {a_m, b_m, d_m, e_m};

        F.function = &computeValueGSL;
        F.params = &params;

        T = gsl_root_fsolver_brent;
        solver = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (solver, &F, x_lo, x_hi);

        do
            {
                iter++;
                status = gsl_root_fsolver_iterate (solver);
                root = gsl_root_fsolver_root (solver);
                x_lo = gsl_root_fsolver_x_lower (solver);
                x_hi = gsl_root_fsolver_x_upper (solver);
                status = gsl_root_test_interval (x_lo, x_hi,
                                                 0, 0.000001);
            }
        while (status == GSL_CONTINUE && iter < max_iter && computeValue(root) > tol);

        gsl_root_fsolver_free (solver);

        return root;
    }

    std::pair<double, double> getSearchRange() const {
        return searchRange_m;
    }

private:

    double a_m; // x^4
    double b_m; // x^3
    double d_m; // x
    double e_m; // 1

    struct PolyParams {
        double a_m;
        double b_m;
        double d_m;
        double e_m;
    };

    static
    double computeValueGSL(double x, void *params) {
        double xcube = std::pow(x, 3);

        struct PolyParams *p
            = (struct PolyParams *) params;

        return p->a_m * x * xcube + p->b_m * xcube + p->d_m * x + p->e_m;
    }

    std::pair<double, double> searchRange_m;

    std::complex<double> deltaZero_m;
    std::complex<double> deltaOne_m;

    std::vector<double> thirdOrderRoots_m;
    // double C_m;

};

#endif