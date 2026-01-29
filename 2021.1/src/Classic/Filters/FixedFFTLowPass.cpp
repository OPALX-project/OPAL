#include "Filters/FixedFFTLowPass.h"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"


FixedFFTLowPassFilter::FixedFFTLowPassFilter(const int &N):
    number_frequencies_m(N)
{ }

void FixedFFTLowPassFilter::apply(std::vector<double> &LineDensity) {
    const int M = LineDensity.size();
    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(M);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(M);
    double *LD = new double[M];
    for (int i = 0; i < M; ++ i) {
        LD[i] = LineDensity[i];
    }
    gsl_fft_real_transform(LD, 1, M, real, work);

    for (int i = 2 * number_frequencies_m + 1; i < M; ++ i) {
        LD[i] = 0.0;
    }

    gsl_fft_real_wavetable_free(real);
    hc = gsl_fft_halfcomplex_wavetable_alloc(M);

    gsl_fft_halfcomplex_inverse(LD, 1, M, hc, work);

    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);

    for (int i = 0; i < M; ++ i) {
        LineDensity[i] = LD[i];
    }
    delete[] LD;
}

void FixedFFTLowPassFilter::calc_derivative(std::vector<double> &LineDensity, const double &h) {
    const int M = LineDensity.size();
    const double gff = 2. * Physics::pi / (h * (M - 1));

    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(M);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(M);
    double *LD = new double[M];

    for (int i = 0; i < M; ++ i) {
        LD[i] = LineDensity[i];
    }
    gsl_fft_real_transform(LD, 1, M, real, work);

    gsl_fft_real_wavetable_free(real);

    LD[0] = 0.0;
    for (int i = 1; i < 2 * number_frequencies_m + 1; i += 2) {
    	double temp = LD[i];
        LD[i] = -LD[i+1] * gff * (i + 1) / 2;
        LD[i+1] = temp * gff * (i + 1) / 2;
    }

    for (int i = 2 * number_frequencies_m + 1; i < M; ++ i) {
        LD[i] = 0.0;
    }

    hc = gsl_fft_halfcomplex_wavetable_alloc(M);

    gsl_fft_halfcomplex_inverse(LD, 1, M, hc, work);

    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);

    for (int i = 0; i < M; ++ i) {
        LineDensity[i] = LD[i];
    }

    delete[] LD;
}

