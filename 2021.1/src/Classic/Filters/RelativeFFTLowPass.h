#ifndef CLASSIC_RELATIVE_LOW_PASS_FFT_HH
#define CLASSIC_RELATIVE_LOW_PASS_FFT_HH

#include "Filters/Filter.h"

class RelativeFFTLowPassFilter: public Filter {
public:
    RelativeFFTLowPassFilter(const double &threshold);
    void apply(std::vector<double> &LineDensity);
    void calc_derivative(std::vector<double> &histogram, const double &h);

private:
    double threshold_m;
};

#endif // CLASSIC_RELATIVE_LOW_PASS_FFT_HH
