#include "Filters/RelativeFFTLowPass.h"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

#include <cmath>

RelativeFFTLowPassFilter::RelativeFFTLowPassFilter(const double &threshold):
    threshold_m(threshold)
{ }

void RelativeFFTLowPassFilter::apply(std::vector<double> &LineDensity) {
    const int M = LineDensity.size();
    double max_four_coef = 0.0;

    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(M);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(M);
    double *LD = new double[M];

    for (int i = 0; i < M; ++ i) {
        LD[i] = LineDensity[i];
    }
    gsl_fft_real_transform(LD, 1, M, real, work);

    gsl_fft_real_wavetable_free(real);

    for (int i = 0; i < M; ++ i) {
        if (std::abs(LD[i]) > max_four_coef) {
            max_four_coef = std::abs(LD[i]);
        }
    }
    max_four_coef *= threshold_m;

    for (int i = 0; i < M; ++ i) {
        if (std::abs(LD[i]) < max_four_coef) {
            LD[i] = 0.0;
        }
    }

    gsl_fft_halfcomplex_wavetable *hc = gsl_fft_halfcomplex_wavetable_alloc(M);

    gsl_fft_halfcomplex_inverse(LD, 1, M, hc, work);

    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);

    for (int i = 0; i < M; ++ i) {
        LineDensity[i] = LD[i];
    }

    delete[] LD;
}

void RelativeFFTLowPassFilter::calc_derivative(std::vector<double> &LineDensity, const double &h) {
    const int M = LineDensity.size();
    const double gff = 2.* Physics::pi / (h * (M - 1));
    double max_four_coef = 0.0;

    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(M);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(M);
    double *LD = new double[M];

    for (int i = 0; i < M; ++ i) {
        LD[i] = LineDensity[i];
    }
    gsl_fft_real_transform(LD, 1, M, real, work);

    gsl_fft_real_wavetable_free(real);

    for (int i = 1; i < M; ++ i) {
        if (std::abs(LD[i]) > max_four_coef) {
            max_four_coef = std::abs(LD[i]);
        }
    }
    max_four_coef *= threshold_m;

    LD[0] = 0.0;
    for (int i = 1; i < M; i += 2) {
        double temp = LD[i];
        if (std::abs(LD[i+1]) > max_four_coef) {
            temp = LD[i];
            LD[i] = -LD[i+1] * gff * (i + 1) / 2;
        } else {
            LD[i] = 0.0;
        }
        if (std::abs(temp) > max_four_coef) {
            LD[i+1] = temp * gff * (i + 1) / 2;
        } else {
            LD[i+1] = 0.0;
        }
    }

    gsl_fft_halfcomplex_wavetable *hc = gsl_fft_halfcomplex_wavetable_alloc(M);

    gsl_fft_halfcomplex_inverse(LD, 1, M, hc, work);

    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);

    for (int i = 0; i < M; ++ i) {
        LineDensity[i] = LD[i];
    }

    delete[] LD;
}