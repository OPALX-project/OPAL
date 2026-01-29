#ifndef CLASSIC_FILTER_HH
#define CLASSIC_FILTER_HH

//#include <iostream>
#include <vector>
//#include <cmath>

class Filter {

public:
    Filter() { ;}
    virtual ~Filter() { };
    virtual void apply(std::vector<double> &histogram) = 0;
    virtual void calc_derivative(std::vector<double> &histogram, const double &h) = 0;

};

#endif //CLASSIC_FILTER_HH
