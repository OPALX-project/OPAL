#include "TrimCoils/TrimCoilFit.h"

#include <cmath>

TrimCoilFit::TrimCoilFit(double bmax,
                         double rmin,
                         double rmax,
                         const std::vector<double>& coefnum,
                         const std::vector<double>& coefdenom,
                         const std::vector<double>& coefnumphi,
                         const std::vector<double>& coefdenomphi):

    TrimCoil(bmax, rmin, rmax)
{
    coefs.resize(4);
    coefs[NUM]      = coefnum;
    coefs[DENOM]    = coefdenom;
    coefs[NUMPHI]   = coefnumphi;
    coefs[DENOMPHI] = coefdenomphi;

    // normal polynom if no denominator coefficients (denominator = 1)
    if (coefs[DENOM].empty())
      coefs[DENOM].push_back(1.0);
    if (coefs[DENOMPHI].empty())
      coefs[DENOMPHI].push_back(1.0);
    // default constant nominator
    if (coefs[NUM].empty())
      coefs[NUM].push_back(1.0);
    if (coefs[NUMPHI].empty())
      coefs[NUMPHI].push_back(1.0);
}

void TrimCoilFit::calculateRationalFunction(FunctionType type, double val, double& quot, double& der_quot) const
{
    double num    = 0.0; // numerator
    double dnum   = 0.0; // derivative of numerator
    double powval = 1.0; // power of value

    const std::vector<double>& coefnum   = coefs[type];
    const std::vector<double>& coefdenom = coefs[type+1];

    // add constant
    num += coefnum[0];
    for (std::size_t i = 1; i < coefnum.size(); ++i) {
        dnum   += coefnum[i] * powval * i;
        powval *= val;
        num    += coefnum[i] * powval;
    }
    double denom  = 0.0; // denominator
    double ddenom = 0.0; // derivative of denominator
    powval        = 1.0; // power of value

    // add constant
    denom += coefdenom[0];
    for (std::size_t i = 1; i < coefdenom.size(); ++i) {
        ddenom += coefdenom[i] * powval * i;
        powval *= val;
        denom  += coefdenom[i] * powval;
    }

    quot     = num / denom;
    // derivative with quotient rule
    der_quot = (dnum * denom - ddenom * num) / (denom*denom);
}

void TrimCoilFit::calculateRationalFunction(FunctionType type, double val, double& quot, double& der_quot, double& der2_quot) const
{
    double num     = 0.0; // numerator
    double d_num   = 0.0; // derivative of numerator
    double d2_num  = 0.0; // second derivative of numerator
    double powval  = 1.0; // power of r

    const std::vector<double>& coefnum   = coefs[type];
    const std::vector<double>& coefdenom = coefs[type+1];

    unsigned int order = coefnum.size();

    // add constant and first term
    num += coefnum[0];
    if (order > 1) {
        num   += coefnum[1] * val;
        d_num += coefnum[1];
    }
    for (std::size_t i = 2; i < coefnum.size(); ++i) {
        d2_num += coefnum[i] * powval * i * (i-1);
        powval *= val; // r^(i-1)
        d_num  += coefnum[i] * powval * i;
        num    += coefnum[i] * powval * val;
    }

    double denom    = 0.0; // denominator
    double d_denom  = 0.0; // derivative of denominator
    double d2_denom = 0.0; // derivative of denominator
    powval            = 1.0; // power of r
    order           = coefdenom.size();

    // add constant
    denom += coefdenom[0];
    if (order > 1) {
        denom   += coefdenom[1] * val;
        d_denom += coefdenom[1];
    }
    for (std::size_t i = 2; i < coefdenom.size(); ++i) {
        d2_denom += coefdenom[i] * powval * i * (i-1);
        powval   *= val;
        d_denom  += coefdenom[i] * powval * i;
        denom    += coefdenom[i] * powval * val;
    }

    quot = num / denom;

    // derivative of phase with quotient rule (B - field)
    der_quot = (d_num * denom - d_denom * num) / (denom*denom);

    // second derivitive of phase (dB/dr) with quotient rule
    // (d2_num - 2*(num/denom)' * d_denom - (num/denom) * d2_denom) / denom
    der2_quot = (d2_num - 2*der_quot*d_denom - quot * d2_denom) / denom;
}