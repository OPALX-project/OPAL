/* profile.C
   profile interpolation class

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: profile.C 106 2007-05-08 19:12:24Z bakker $
*/

#include "Ippl.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string.h>

#include "Algorithms/bet/math/sort.h"
#include "Algorithms/bet/math/interpol.h"
#include "Algorithms/bet/math/integrate.h"
#include "Algorithms/bet/profile.h"


// global internal functions for integration
// -----------------------------------------
static Profile *cProfile = NULL;

static double f1(double x) {
    return (cProfile ? cProfile->get(x) : 0.0);
}

static double f2(double x) {
    return (cProfile ? pow(cProfile->get(x), 2) : 0.0);
}

static double f3(double x) {
    return (cProfile ? fabs(cProfile->get(x)) : 0.0);
}

Profile::Profile(double v) {
    n = 0;
    x = NULL;
    y = NULL;
    y2 = NULL;

    sf   = 1.0;
    yMin = v;
    yMax = v;
}

Profile::Profile(double *_x, double *_y, int _n) {
    n  = _n;

    x  = (double *) malloc(sizeof(double) * n);
    y  = (double *) malloc(sizeof(double) * n);
    if((x == NULL) || (y == NULL))
        std::cout << "Profile::Profile: Insuffient memory for Profile(n = " << n << " )" << std::endl;
    memcpy(x, _x, sizeof(double)*n);
    memcpy(y, _y, sizeof(double)*n);

    create();
} /* Profile::Profile() */

Profile::Profile(char *fname, double eps) {
    FILE  *f;
    int   i, i0;
    double a, b, m;

    f = fopen(fname, "r");
    if(!f) {
        std::cout << "Profile::Profile: Cannot load profile mapping " << fname << std::endl;
    }
    i = 0;
    m = 0.0;
    while(fscanf(f, "%lf %lf", &a, &b) == 2) {
        if(fabs(b) > m) m = fabs(b);
        ++i;
    }
    fclose(f);

    n = i;
    x = (double *) malloc(sizeof(double) * n);
    y = (double *) malloc(sizeof(double) * n);

    if((x == NULL) || (y == NULL)) {
        std::cout << "Profile::Profile: Insuffient memory for Profile(n = " << n << " )" << std::endl;
    }
    // read all values
    f = fopen(fname, "r");
    for(i = 0; i < n; i++) {
      int res = fscanf(f, "%lf %lf", &x[i], &y[i]);
      if (res !=0)
          ERRORMSG("fscanf in profile.cpp has res!=0" << endl);
    }
    fclose(f);

    // cut tails (if applicable)
    m = fabs(m * eps);
    // cut start
    i0 = 0;
    while((i0 < n) && (fabs(y[i0]) < m)) ++i0;
    if((i0 > 0) && (i0 < n)) {
        for(i = i0; i < n; i++) {
            x[i-i0] = x[i];
            y[i-i0] = y[i];
        }
        n -= i0;
    }
    // cut end
    i0 = n - 1;
    while((i0 >= 0) && (fabs(y[i0]) < m)) --i0;
    if((i0 < (n - 1)) && (i0 >= 0)) n = i0;

    create();
} /* Profile::Profile() */

Profile::~Profile() {
    if(x)  free(x);
    if(y)  free(y);
    if(y2) free(y2);
} /* Profile::~Profile() */

void Profile::create() {
    int
    i, j, k;
    double
    sx, sy;

    sf = 1.0;

    y2   = (double *) malloc(sizeof(double) * n);
    yMax = y[0];
    yMin = yMax;

    if(y2 == NULL) {
        std::cout << "Profile::Profile: Insuffient memory for Profile y2 (n = " << n << " )" << std::endl;
    }
    for(i = 0; i < n; i++) {
        if(y[i] < yMin) yMin = y[i];
        if(y[i] > yMax) yMax = y[i];
    }
    sort2(x, y, n);

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

    spline(x, y, n, y2);
} /* Profile::create() */

double Profile::get(double xa, Interpol_type tp) {
    double val = 0.0;

    if(x) {
        if(xa < x[0]) val = 0.0;
        else if(xa > x[n-1]) val = 0.0;
        else switch(tp) {
                case itype_lin :
                    lsplint(x, y, y2, n, xa, &val);
                    break;
                default :
                    lsplint(x, y, y2, n, xa, &val);
                    break;
            }
    }
    return (sf * val);
} /* Profile::get() */

void Profile::normalize() {
    if(yMax > 0.0) sf = 1.0 / yMax;
    else if(yMin != 0.0) sf = 1.0 / fabs(yMin);
    else sf = 1.0;
} /* Profile:: normalize */

void Profile::scale(double v) {
    sf *= v;
} /* Profile::scale() */

double Profile::set(double f) {
    double v = fabs(fabs(yMax) > fabs(yMin) ? yMax : yMin);

    if(v > 0.0) sf = f / v;
    else sf = 1.0;

    return sf;
} /* Profile::set() */

void Profile::setSF(double value) {
    sf =  value;
} /* Profile::setSF() */

double Profile::getSF() {
    return sf;
} /* Profile::getSF() */

void Profile::dump(char *fname, double dx) {
    FILE *f;

    f = fopen(fname, "w");
    if(f) {
        dump(f, dx);
        fclose(f);
    } else {
        std::cout << "Profile::Profile: Failed to create profile output:" << fname << std::endl;
    }
}

void Profile::dump(FILE *f, double dx) {
    int    i;
    double xx, dxx;

    fprintf(f, "SDDS1\n");
    fprintf(f, "&parameter name=n, type=long, fixed_value=%d &end\n", n);
    fprintf(f, "&parameter name=sf, type=double, fixed_value=%20.12le &end\n", sf);
    fprintf(f, "&column name=x,    type=double &end\n");
    fprintf(f, "&column name=y,    type=double &end\n");
    fprintf(f, "&data mode=ascii &end\n");

    fprintf(f, "! next page\n");
    fprintf(f, "  %d\n", n);
    for(i = 0; i < n; i++) {
        fprintf(f, "%20.12le \t %20.12le\n", x[i] + dx, sf * y[i]);
    }
    fprintf(f, "! next page\n");
    fprintf(f, " %d\n", 10 * n);

    dxx = (x[n-1] - x[0]) / (10 * n - 1);
    for(i = 0; i < 10 * n; i++) {
        xx = x[0] + dxx * i;
        fprintf(f, "%20.12le \t %20.12le\n", xx + dx, get(xx));
    }
} /* Profile::dump() */

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
    return (x ? x[n-1] : 0.0);
}

double Profile::xMin() {
    return (x ? x[0] : 0.0);
}

double Profile::Leff() {
    double ym;

    ym       = fabs((fabs(yMin) > fabs(yMax)) ? yMin : yMax);
    cProfile = this;
    return (((x == NULL) || (x[n-1] == x[0]) || (ym == 0.0)) ? 0.0 :
            fabs(qromb(f1, x[0], x[n-1]) / ym));
}

double Profile::Leff2() {
    double ym;

    ym       = pow((fabs(yMin) > fabs(yMax)) ? yMin : yMax, 2);
    cProfile = this;
    return (((x == NULL) || (x[n-1] == x[0]) || (ym == 0.0)) ? 0.0 :
            fabs(qromb(f2, x[0], x[n-1]) / ym));
}

double Profile::Labs() {
    double ym;

    ym       = fabs((fabs(yMin) > fabs(yMax)) ? yMin : yMax);
    cProfile = this;
    return (((x == NULL) || (x[n-1] == x[0]) || (ym == 0.0)) ? 0.0 :
            fabs(qromb(f3, x[0], x[n-1]) / ym));
}