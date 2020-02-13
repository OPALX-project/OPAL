#include "Filters/SavitzkyGolay.h"
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"
#include <iostream>

using namespace std;

SavitzkyGolayFilter::SavitzkyGolayFilter(int np, int nl, int nr, int m):
    NumberPoints_m(np),
    NumberPointsLeft_m(nl),
    NumberPointsRight_m(nr),
    PolynomialOrder_m(m),
    Coefs_m(np, 0.0),
    CoefsDeriv_m(np, 0.0) {
    savgol(Coefs_m, NumberPoints_m, NumberPointsLeft_m, NumberPointsRight_m, 0, PolynomialOrder_m);
    savgol(CoefsDeriv_m, NumberPoints_m, NumberPointsLeft_m, NumberPointsRight_m, 1, PolynomialOrder_m);
}

void SavitzkyGolayFilter::apply(vector<double> &LineDensity) {
    vector<double> temp(LineDensity.size(), 0.0);
    convlv(LineDensity, Coefs_m, 1, temp);
    LineDensity.assign(temp.begin(), temp.end());
}

void SavitzkyGolayFilter::calc_derivative(vector<double> &LineDensity, const double &/*h*/) {
    vector<double> temp(LineDensity.size(), 0.0);
    convlv(LineDensity, CoefsDeriv_m, 1, temp);
    LineDensity.assign(temp.begin(), temp.end());
}


void savgol(vector<double> &c, const int &np, const int &nl, const int &nr, const int &ld, const int &m) {
    int j, k, imj, ipj, kk, mm;
    double d, fac, sum;

    if(np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m) {
        cerr << "bad args in savgol" << endl;
        return;
    }
    vector<int> indx(m + 1, 0);
    vector<double> a((m + 1) * (m + 1), 0.0);
    vector<double> b(m + 1, 0.0);

    for(ipj = 0; ipj <= (m << 1); ++ipj) {
        sum = (ipj ? 0.0 : 1.0);
        for(k = 1; k <= nr; ++k)
            sum += (int)pow((double)k, (double)ipj);
        for(k = 1; k <= nl; ++k)
            sum += (int)pow((double) - k, (double)ipj);
        mm = (ipj < 2 * m - ipj ? ipj : 2 * m - ipj);
        for(imj = -mm; imj <= mm; imj += 2)
            a[(ipj + imj) / 2 * (m + 1) + (ipj - imj) / 2] = sum;
    }
    ludcmp(a, indx, d);

    for(j = 0; j < m + 1; ++j)
        b[j] = 0.0;
    b[ld] = 1.0;

    lubksb(a, indx, b);
    for(kk = 0; kk < np; ++kk)
        c[kk] = 0.0;
    for(k = -nl; k <= nr; ++k) {
        sum = b[0];
        fac = 1.0;
        for(mm = 1; mm <= m; ++mm)
            sum += b[mm] * (fac *= k);
        kk = (np - k) % np;
        c[kk] = sum;
    }

}

void convlv(const vector<double> &data, const vector<double> &respns, const int &isign, vector<double> &ans) {
    int n = data.size();
    int m = respns.size();

    double tmp;

    double *tempfd1 = new double[n];
    double *tempfd2 = new double[n];

    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(n);


    for(int i = 0; i < n; ++ i) {
        tempfd1[i] = 0.0;
        tempfd2[i] = data[i];
    }

    tempfd1[0] = respns[0];
    for(int i = 1; i < (m + 1) / 2; ++i) {
        tempfd1[i] = respns[i];
        tempfd1[n - i] = respns[m - i];
    }

    gsl_fft_real_transform(tempfd1, 1, n, real, work);
    gsl_fft_real_transform(tempfd2, 1, n, real, work);

    gsl_fft_real_wavetable_free(real);
    hc = gsl_fft_halfcomplex_wavetable_alloc(n);

    if(isign == 1) {
        for(int i = 1; i < n - 1; i += 2) {
            tmp = tempfd1[i];
            tempfd1[i] = (tempfd1[i] * tempfd2[i] - tempfd1[i + 1] * tempfd2[i + 1]);
            tempfd1[i + 1] = (tempfd1[i + 1] * tempfd2[i] + tmp * tempfd2[i + 1]);
        }
        tempfd1[0] *= tempfd2[0];
        tempfd1[n - 1] *= tempfd2[n - 1];
    }

    gsl_fft_halfcomplex_inverse(tempfd1, 1, n, hc, work);

    for(int i = 0; i < n; ++i) {
        ans[i] = tempfd1[i];
    }

    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);

    delete[] tempfd1;
    delete[] tempfd2;
}

void ludcmp(vector<double> &a, vector<int> &indx,  double &d) {
    const double TINY = 1.0e-20;
    int i, imax = -1, j, k;
    double big, dum, sum, temp;

    int n = indx.size();
    vector<double> vv(n, 0.0);

    d = 1.0;
    for(i = 0; i < n; ++i) {
        big = 0.0;
        for(j = 0; j < n; ++j)
            if((temp = fabs(a[i * n + j])) > big) big = temp;

        if(big == 0.0) {
            cerr << "Singular matrix in routine ludcmp" << endl;
            return;
        }
        vv[i] = 1. / big;
    }

    for(j = 0; j < n; ++j) {
        for(i = 0; i < j; ++i) {
            sum = a[i * n + j];
            for(k = 0; k < i; ++k)
                sum -= a[i * n + k] * a[k * n + j];
            a[i * n + j] = sum;
        }
        big = 0.0;
        for(i = j; i < n; ++i) {
            sum = a[i * n + j];
            for(k = 0; k < j; ++k)
                sum -= a[i * n + k] * a[k * n + j];
            a[i * n + j] = sum;
            if((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }

        if(j != imax) {
            for(k = 0; k < n; ++k) {
                dum = a[imax * n + k];
                a[imax * n + k] = a[j * n + k];
                a[j * n + k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if(a[j * n + j] == 0.0) a[j * n + j] = TINY;
        if(j != n - 1) {
            dum = 1. / a[j * n + j];
            for(i = j + 1; i < n; ++i)
                a[i * n + j] *= dum;
        }
    }
}

void lubksb(vector<double> &a, vector<int> &indx, vector<double> &b) {
    int i, ii = 0, ip, j;
    double sum;
    int n = indx.size();

    for(i = 0; i < n; ++i) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if(ii != 0)
            for(j = ii - 1; j < i; ++j)
                sum -= a[i * n + j] * b[j];
        else if(sum != 0.0)
            ii = i + 1;
        b[i] = sum;
    }
    for(i = n - 1; i >= 0; --i) {
        sum = b[i];
        for(j = i + 1; j < n; ++j)
            sum -= a[i * n + j] * b[j];
        b[i] = sum / a[i * n + i];
    }

}
