#ifndef CLASSIC_STENCIL_HH
#define CLASSIC_STENCIL_HH

#include "Filters/Filter.h"

class StencilFilter: public Filter {
public:
    StencilFilter() { ;}
    void apply(std::vector<double> &histogram);
    void calc_derivative(std::vector<double> &histogram, const double &h);
};

#endif // CLASSIC_STENCIL_HH