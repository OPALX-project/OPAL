#ifndef CLASSIC_SAVITZKY_GOLAY_FILTER_HH
#define CLASSIC_SAVITZKY_GOLAY_FILTER_HH

#include "Filters/Filter.h"

class SavitzkyGolayFilter: public Filter {
public:
    SavitzkyGolayFilter(int np, int nl, int nr, int m);
    ~SavitzkyGolayFilter()
    { ;}

    void apply(std::vector<double> &histogram);
    void calc_derivative(std::vector<double> &histogram, const double &h);

private:
    int NumberPoints_m;
    int NumberPointsLeft_m;
    int NumberPointsRight_m;
    int PolynomialOrder_m;
    std::vector<double> Coefs_m;
    std::vector<double> CoefsDeriv_m;

};

void savgol(std::vector<double> &c, const int &np, const int &nl, const int &nr, const int &ld, const int &m);
void convlv(const std::vector<double> &data, const std::vector<double> &respns, const int &isign, std::vector<double> &ans);
void ludcmp(std::vector<double> &a, std::vector<int> &indx, double &d);
void lubksb(std::vector<double> &a, std::vector<int> &indx, std::vector<double> &b);

#endif // CLASSIC_SAVITZKY_GOLAY_FILTER_HH
