#include "Filters/Stencil.h"


void StencilFilter::apply(std::vector<double> &LineDensity) {
    int size = LineDensity.size();
    std::vector<double> temp1(size + 8, 0.0);
    std::vector<double> temp2(5, 0.0);

    for(int i = 0; i < size; ++i)
        temp1[i+4] = LineDensity[i];

    for(int i = 0; i < 4; ++ i) {
        temp1[3-i] = temp1[5+i];
        temp1[size+5+i] = temp1[size+3-i];
    }

    for(int i = 0; i < size; ++i) {
        temp2[i%5] = (7. * (temp1[i] + temp1[i+8]) + 24. * (temp1[i+2] + temp1[i+6]) + 34.*temp1[i+4]) / 96.;
        temp1[i] = temp2[(i+1)%5];
    }

    temp1[2] = temp1[6];
    temp1[3] = temp1[5];
    temp1[size+5] = temp1[size+1];
    temp1[size+4] = temp1[size+2];

    for(int i = 0; i < size; ++i)
        LineDensity[i] = (7. * (temp1[i+2] + temp1[i+6]) + 24. * (temp1[i+3] + temp1[i+5]) + 34.*temp1[i+4]) / 96.;
}

void StencilFilter::calc_derivative(std::vector<double> &LineDensity, const double &h) {
    std::vector<double> temp(LineDensity.begin(), LineDensity.end());
    const double N = LineDensity.size();
    LineDensity[0] = (-25.*temp[0] + 48.*temp[1] - 36.*temp[2] + 16.*temp[3] + temp[4]) / (12.*h);
    LineDensity[1] = (-3.*temp[0] - 10.*temp[1] + 18.*temp[2] - 6.*temp[3] + temp[4]) / (12.*h);
    for(int i = 2; i < N - 2; ++ i) {
        LineDensity[i] = (temp[i-2] - 8.*temp[i-1] + 8.*temp[i+1] - temp[i+2]) / (12.*h);
    }
    LineDensity[N-2] = (-temp[N-5] + 6.*temp[N-4] - 18.*temp[N-3] + 10.*temp[N-2] + 3.*temp[N-1]) / (12.*h);
    LineDensity[N-1] = (3.*temp[N-5] - 16.*temp[N-4] + 36.*temp[N-3] - 48.*temp[N-2] + 25.*temp[N-1]) / (12.*h);
    apply(LineDensity);
}