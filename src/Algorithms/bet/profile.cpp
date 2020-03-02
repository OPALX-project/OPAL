/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 07-03-2006

   calculates a functional profile from a mapping
*/

#include <cmath>
#include <string>
#include <cstdio>
#include <algorithm>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "Algorithms/bet/profile.h"
#include "Utilities/OpalException.h"

// global internal functions for integration
// -----------------------------------------
static Profile *cProfile = NULL;

static double f1(double x, void*) {
    return (cProfile ? cProfile->get(x) : 0.0);
}

static double f2(double x, void*) {
    return (cProfile ? std::pow(cProfile->get(x), 2) : 0.0);
}

static double f3(double x, void*) {
    return (cProfile ? std::abs(cProfile->get(x)) : 0.0);
}

Profile::Profile(double v) {
    n = 0;

    sf   = 1.0;
    yMin = v;
    yMax = v;
    acc = nullptr;
    spline = nullptr;
}

Profile::Profile(double *_x, double *_y, int _n) :
    n(_n), x(_x, _x + _n), y(_y, _y + _n) {
    create();
}

Profile::Profile(char *fname, double eps) {
    FILE  *f;
    int   i, i0;
    double a, b, m;

    f = fopen(fname, "r");
    if(!f) {
        std::string s(fname);
        throw OpalException(
            "Profile::Profile()",
            "Cannot load profile mapping \"" + s + "\".");
    }
    i = 0;
    m = 0.0;
    while(fscanf(f, "%lf %lf", &a, &b) == 2) {
        if(std::abs(b) > m) m = std::abs(b);
        ++i;
    }
    fclose(f);

    n = i;
    x.resize(n);
    y.resize(n);

    // read all values
    f = fopen(fname, "r");
    for(i = 0; i < n; i++) {
        int res = fscanf(f, "%lf %lf", &x[i], &y[i]);
        if (res !=0) {
            std::string s(fname);
            throw OpalException(
                "Profile::Profile()",
                "Cannot read profile mapping \"" + s + "\".");
        }
    }
    fclose(f);

    // cut tails (if applicable)
    m = std::abs(m * eps);
    // cut start
    i0 = 0;
    while((i0 < n) && (std::abs(y[i0]) < m)) ++i0;
    if((i0 > 0) && (i0 < n)) {
        for(i = i0; i < n; i++) {
            x[i-i0] = x[i];
            y[i-i0] = y[i];
        }
        n -= i0;
    }
    // cut end
    i0 = n - 1;
    while((i0 >= 0) && (std::abs(y[i0]) < m)) --i0;
    if((i0 < (n - 1)) && (i0 >= 0)) n = i0;

    create();
}

Profile::~Profile() {
    if (spline != nullptr) {
        gsl_spline_free (spline);
    }
    if (acc != nullptr) {
        gsl_interp_accel_free (acc);
    }
}

void Profile::create() {
    int
    i, j, k;
    double
    sx, sy;

    sf = 1.0;

    y2.resize(n);
    yMax = y[0];
    yMin = yMax;

    for(i = 0; i < n; i++) {
        if(y[i] < yMin) yMin = y[i];
        if(y[i] > yMax) yMax = y[i];
    }

    //sort y according to x
    std::vector<std::pair<double,double>> xy(n);
    for (int i=0; i<n; i++) {
        xy[i] = std::make_pair(x[i],y[i]);
    }

    std::sort(xy.begin(), xy.end(),
              [](const std::pair<double,double> &left, const std::pair<double,double> &right)
              {return left.first < right.first;});

    for (int i=0; i<n; i++) {
        x[i] = xy[i].first;
        y[i] = xy[i].second;
    }


    /* remove points with identical x-value
       take the average if the situation does occur */
    j = 0;
    for(i = 0; i + j < n; i++) {
        k  = 1;
        sx = x[i+j];
        sy = y[i+j];
        while((i + j + k < n) && (x[i+j] == x[i+j+k])) {
            sx += x[i+j+k];
            sy += y[i+j];
            ++k;
        }
        x[i] = sx / k;
        y[i] = sy / k;
        //printf("i = %d:\tx = %e, y = %e\n",i,x[i],y[i]);
        j   += (k - 1);
    }
    n = i;

    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, &x[0], &y[0], n);
}

double Profile::get(double xa, Interpol_type /*tp*/) {
    double val = 0.0;

    if (x.empty() || xa < x[0] || xa > x[n-1]) {
        return 0.0;
    }
    /*
      natural cubic spline interpolation.
      The interpolated value is limited between the y-values
      of the adjacent points.
    */
    val = gsl_spline_eval (spline, xa, acc);
    size_t low = 0;
    size_t high = n - 1;
    while (high - low > 1) {
        size_t k = (low + high) >> 1;
        if (x[k] > xa) {
            high = k;
        } else {
            low = k;
        }
    }
    double y_min = std::min (y[low], y[high]);
    double y_max = std::max (y[low], y[high]);
    if (val < y_min) {
        val = y_min;
    } else if (val > y_max) {
        val = y_max;
    }
    return (sf * val);
}

void Profile::normalize() {
    if(yMax > 0.0)
        sf = 1.0 / yMax;
    else if(yMin != 0.0)
        sf = 1.0 / std::abs(yMin);
    else
        sf = 1.0;
}

void Profile::scale(double v) {
    sf *= v;
}

double Profile::set(double f) {
    double v = std::abs(std::abs(yMax) > std::abs(yMin) ? yMax : yMin);

    if(v > 0.0) sf = f / v;
    else sf = 1.0;

    return sf;
}

void Profile::setSF(double value) {
    sf =  value;
}

double Profile::getSF() {
    return sf;
}

void Profile::dump(char fname[], double dx) {
    FILE *f;

    f = fopen(fname, "w");
    if(!f) {
        std::string s(fname);
        throw OpalException(
            "Profile::dump()",
            "Cannot dump profile \"" + s + "\".");
    }
    dump(f, dx);
    fclose(f);
}

void Profile::dump(FILE *f, double dx) {
    int    i;
    double xx, dxx;

    std::fprintf(f, "SDDS1\n");
    std::fprintf(f, "&parameter name=n, type=long, fixed_value=%d &end\n", n);
    std::fprintf(f, "&parameter name=sf, type=double, fixed_value=%20.12le &end\n", sf);
    std::fprintf(f, "&column name=x,    type=double &end\n");
    std::fprintf(f, "&column name=y,    type=double &end\n");
    std::fprintf(f, "&data mode=ascii &end\n");

    std::fprintf(f, "! next page\n");
    std::fprintf(f, "  %d\n", n);
    for(i = 0; i < n; i++) {
        std::fprintf(f, "%20.12le \t %20.12le\n", x[i] + dx, sf * y[i]);
    }
    std::fprintf(f, "! next page\n");
    std::fprintf(f, " %d\n", 10 * n);

    dxx = (x[n-1] - x[0]) / (10 * n - 1);
    for(i = 0; i < 10 * n; i++) {
        xx = x[0] + dxx * i;
        std::fprintf(f, "%20.12le \t %20.12le\n", xx + dx, get(xx));
    }
}

int Profile::getN() {
    return n;
}

double Profile::min() {
    return sf * yMin;
}

double Profile::max() {
    return sf * yMax;
}

double Profile::xMax() {
    return ((x.empty()==false) ? x[n-1] : 0.0);
}

double Profile::xMin() {
    return ((x.empty()==false) ? x[0] : 0.0);
}

double Profile::Leff() {
    double ym;

    ym       = std::abs((std::abs(yMin) > std::abs(yMax)) ? yMin : yMax);
    cProfile = this;
    double result = 0.0;
    if (!x.empty() && (x[n-1] != x[0]) && (ym != 0.0)) {
        gsl_function F = { &f1, NULL};
        gsl_integration_romberg_workspace* w = gsl_integration_romberg_alloc (30);
        size_t neval = 0;
        gsl_integration_romberg (&F, x[0], x[n-1], 1.0e-4, 1000, &result, &neval, w);
        gsl_integration_romberg_free (w);
        result /= ym;
    }
    return result;
}

double Profile::Leff2() {
    double ym;

    ym       = std::pow((std::abs(yMin) > std::abs(yMax)) ? yMin : yMax, 2);
    cProfile = this;
    double result = 0.0;
    if (!x.empty() && (x[n-1] != x[0]) && (ym != 0.0)) {
        gsl_function F = {&f2, NULL};
        gsl_integration_romberg_workspace* w = gsl_integration_romberg_alloc (30);
        size_t neval = 0;
        gsl_integration_romberg (&F, x[0], x[n-1], 1.0e-4, 1000, &result, &neval, w);
        gsl_integration_romberg_free (w);
        result /= ym;
    }
    return result;
}

double Profile::Labs() {
    double ym;

    ym       = std::abs((std::abs(yMin) > std::abs(yMax)) ? yMin : yMax);
    cProfile = this;
    double result = 0.0;
    if (!x.empty() && (x[n-1] != x[0]) && (ym != 0.0)) {
        gsl_function F = {&f3, NULL};
        gsl_integration_romberg_workspace* w = gsl_integration_romberg_alloc (30);
        size_t neval = 0;
        gsl_integration_romberg (&F, x[0], x[n-1], 1.0e-4, 1000, &result, &neval, w);
        gsl_integration_romberg_free (w);
        result /= ym;
    }
    return result;
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End: