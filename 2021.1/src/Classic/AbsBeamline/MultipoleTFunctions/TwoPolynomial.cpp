/*
 *  Copyright (c) 2018, Martin Duy Tat
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "gsl/gsl_fft_complex.h"
#include "TwoPolynomial.h"

namespace polynomial {

TwoPolynomial::TwoPolynomial(): maxXorder_m(0), maxSorder_m(0) {
    std::vector<int> temp;
    temp.push_back(0);
    coefficients_m.push_back(temp);
}

TwoPolynomial::TwoPolynomial(const std::vector<std::vector<int>> &coefficients):
    maxXorder_m(coefficients.size() - 1),
    coefficients_m(coefficients) {
    if (coefficients.size() == 0) {
        maxXorder_m = 0;
        maxSorder_m = 0;
        std::vector<int> temp;
        temp.push_back(0);
        coefficients_m.push_back(temp);
        return;
    }
    maxSorder_m = coefficients[0].size() - 1;
    std::size_t len = maxSorder_m + 1;
    for (std::size_t i = 0; i < coefficients.size(); i++) {
        if (coefficients[i].size() != len) {
            throw std::length_error("2D vector not rectangular");
        }
        if (coefficients[i].size() == 0) {
            maxXorder_m = 0;
            maxSorder_m = 0;
            std::vector<int> temp;
            temp.push_back(0);
            coefficients_m.resize(0);
            coefficients_m.push_back(temp);
            return;
        }
    }
}

TwoPolynomial::TwoPolynomial(const TwoPolynomial &poly):
    maxXorder_m(poly.maxXorder_m), maxSorder_m(poly.maxSorder_m),
    coefficients_m(poly.coefficients_m), dSfactors_m(poly.dSfactors_m) {
}

TwoPolynomial::~TwoPolynomial() {
}

TwoPolynomial& TwoPolynomial::operator= (const TwoPolynomial &poly) {
    maxXorder_m = poly.maxXorder_m;
    maxSorder_m = poly.maxSorder_m;
    coefficients_m = poly.coefficients_m;
    dSfactors_m = poly.dSfactors_m;
    return *this;
}

void TwoPolynomial::differentiateX() {
    if (maxXorder_m == 0 && maxSorder_m == 0 && coefficients_m[0][0] == 0) {
        return;
    }
    if (maxXorder_m == 0) {
        std::vector<std::vector<int>> temp {{0}};
        coefficients_m = temp;
        maxSorder_m = 0;
        return;
    }
    for (std::size_t i = 0; i < maxXorder_m; i++) {
        for (std::size_t j = 0; j <= maxSorder_m; j++) {
            coefficients_m[i][j] = coefficients_m[i + 1][j] * (i + 1);
        }
    }
    coefficients_m.pop_back();
    maxXorder_m--;
}

void TwoPolynomial::differentiateS() {
    if (maxXorder_m == 0 && maxSorder_m == 0 && coefficients_m[0][0] == 0) {
        return;
    }
    if (maxSorder_m == 0) {
        std::vector<std::vector<int>> temp {{0}};
        coefficients_m = temp;
        maxXorder_m = 0;
        return;
    }
    for (std::size_t i = 0; i <= maxXorder_m; i++) {
        for (std::size_t j = 0; j < maxSorder_m; j++) {
            coefficients_m[i][j] = coefficients_m[i][j + 1] * (j + 1);
        }
        coefficients_m[i].pop_back();
    }
    if (dSfactors_m.size() != 0) {
        dSfactors_m[0] = dSfactors_m[0] + 1;
    } else {
        dSfactors_m.push_back(1);
    }
    maxSorder_m--;
}

void TwoPolynomial::multiplyConstant (const int &constant) {
    if (maxXorder_m == 0 && maxSorder_m == 0 && coefficients_m[0][0] == 0) {
        return;
    }
    if (constant == 0) {
        setZero();
        return;
    }
    for (std::size_t i = 0; i <= maxXorder_m; i++) {
        for (std::size_t j = 0; j <= maxSorder_m; j++) {
            coefficients_m[i][j] = constant * coefficients_m[i][j];
        }
    }
}

void TwoPolynomial::multiplydSfactors(const TwoPolynomial &poly) {
    if (poly.dSfactors_m.size() > dSfactors_m.size()) {
        dSfactors_m.resize(poly.dSfactors_m.size(), 0);
    }
    for (std::size_t k = 0; k < poly.dSfactors_m.size(); k++) {
        dSfactors_m[k] = dSfactors_m[k] + poly.dSfactors_m[k];
    }
}

void TwoPolynomial::convert2Dto1Darray(double *vec1, double *vec2,
                                       const vectorLengths &nn,
                                       const TwoPolynomial &poly) const {
    for (std::size_t i = 0; i < nn.nx; i++) {
        for (std::size_t j = 0; j < nn.ny; j++) {
            if (i < nn.nx1 && j < nn.ny1) {
                vec1[2 * (nn.ny * i + j)] = coefficients_m[i][j];
                vec1[2 * (nn.ny * i + j) + 1] = 0;
            } else {
                vec1[2 * (nn.ny * i + j)] = 0;
                vec1[2 * (nn.ny * i + j) + 1] = 0;
            }
            if (i < nn.nx2 && j < nn.ny2) {
                vec2[2 * (nn.ny * i + j)] = poly.coefficients_m[i][j];
                vec2[2 * (nn.ny * i + j) + 1] = 0;
            } else {
                vec2[2 * (nn.ny * i + j)] = 0;
                vec2[2 * (nn.ny * i + j) + 1] = 0;
            }
        }
    }
}

void TwoPolynomial::convolution(double *vec1,
                                double *vec2,
                                const std::size_t &nx,
                                const std::size_t &ny) const {
    /** Do 2D FFT, 1D FFT on each row and column */
    for (std::size_t i = 0; i < nx; i++) {
        double *p1 = &vec1[0] + 2 * i * ny;
        gsl_fft_complex_radix2_forward(p1, 1, ny);
        double *p2 = &vec2[0] + 2 * i * ny;
        gsl_fft_complex_radix2_forward(p2, 1, ny);
    }
    for (std::size_t j = 0; j < ny; j++) {
        double *p1 = &vec1[0] + 2 * j;
        gsl_fft_complex_radix2_forward(p1, ny, nx);
        double *p2 = &vec2[0] + 2 * j;
        gsl_fft_complex_radix2_forward(p2, ny, nx);
    }
    /** Convolution, multiply each element (complex numbers) */
    for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 0; j < ny; j++) {
            double real1 = vec1[2 * (ny * i + j)];
            double real2 = vec2[2 * (ny * i + j)];
            double imag1 = vec1[2 * (ny * i + j) + 1];
            double imag2 = vec2[2 * (ny * i + j) + 1];
            vec1[2 * (ny * i + j)] = real1 * real2 - imag1 * imag2;
            vec1[2 * (ny * i + j) + 1] = real1 * imag2 + imag1 * real2;
        }
    }
    /** Inverse 2D FFT */
    for (std::size_t i = 0; i < nx; i++) {
        double *p1 = &vec1[0] + 2 * i * ny;
        gsl_fft_complex_radix2_inverse(p1, 1, ny);
    }
    for (std::size_t j = 0; j < ny; j++) {
        double *p1 = &vec1[0] + 2 * j;
        gsl_fft_complex_radix2_inverse(p1, ny, nx);
    }
}
    
void TwoPolynomial::convert1Dto2Darray(double *vec1,
                                       const vectorLengths &nn) {
    coefficients_m.resize(nn.nx1 + nn.nx2 - 1);
    for (std::size_t i = 0; i < (nn.nx1 + nn.nx2 - 1); i++) {
        coefficients_m[i].resize(nn.ny1 + nn.ny2 - 1);
        for (std::size_t j = 0; j < (nn.ny1 + nn.ny2 - 1); j++) {
            coefficients_m[i][j] = std::round(vec1[2 * (nn.ny * i + j)]);
        }
    }
}

void TwoPolynomial::multiplyPolynomial(const TwoPolynomial &poly) {
    if (maxXorder_m == 0 && maxSorder_m == 0 && coefficients_m[0][0] == 0) {
        return;
    }
    if (poly.maxXorder_m == 0 &&
        poly.maxSorder_m == 0 &&
        poly.coefficients_m[0][0] == 0) {
        setZero();
        return;
    }
    /** Multiply the S(s)-derivatives */
    multiplydSfactors(poly);
    /** Find vector sizes and round up to nearest power of 2 */
    vectorLengths nn;
    nn.nx = 1;
    nn.ny = 1;
    nn.nx1 = coefficients_m.size();
    nn.nx2 = poly.coefficients_m.size();
    nn.ny1 = coefficients_m[0].size();
    nn.ny2 = poly.coefficients_m[0].size();
    while ((nn.nx1 + nn.nx2 - 1) > nn.nx) {
        nn.nx *= 2;
    }
    while ((nn.ny1 + nn.ny2 - 1) > nn.ny) {
        nn.ny *= 2;
    }
    /** Allocate memory for 1D array */
    double *vec1 = new double[2 * nn.nx * nn.ny];
    double *vec2 = new double[2 * nn.nx * nn.ny];
    /** Convert 2D vectors into 1D vectors with indexing, change to double */
    convert2Dto1Darray(vec1, vec2, nn, poly);
    /** Do convolution */
    convolution(vec1, vec2, nn.nx, nn.ny);
    /** Convert back to 2D vector and change to int by rounding */
    convert1Dto2Darray(vec1, nn);
    delete[] vec1;
    delete[] vec2;
    maxXorder_m = coefficients_m.size() - 1;
    maxSorder_m = coefficients_m[0].size() - 1;
}

void TwoPolynomial::printPolynomial() const {
    if (maxXorder_m == 0 && maxSorder_m == 0 && coefficients_m[0][0] == 0) {
        return;
    }
    std::cout << "(";
    for (std::size_t i = 0; i <= maxXorder_m; i++) {
        for (std::size_t j = 0; j <= maxSorder_m; j++) {
            if (coefficients_m[i][j] > 0) {
                std::cout << " + " << coefficients_m[i][j];
                std::cout << "x^" << i << "s^" << j;
            } else if (coefficients_m[i][j] < 0) {
                std::cout << " - " << -coefficients_m[i][j];
                std::cout << "x^" << i << "s^" << j;
            }
        }
    }
    std::cout << ")";
    for (std::size_t i = 0; i < dSfactors_m.size(); i++) {
        std::cout << "(d^" << i + 1 << "S/ds)^" << dSfactors_m[i];
    }
}

int TwoPolynomial::getCoefficient(const std::size_t &Xorder,
                                  const std::size_t &Sorder) const {
    if (Xorder > maxXorder_m || Sorder > maxSorder_m) {
        return 0;
    }
    return coefficients_m[Xorder][Sorder];
}

void TwoPolynomial::setCoefficient(const int &coefficient,
                                   const std::size_t &Xorder,
                                   const std::size_t &Sorder) {
    if (Xorder > maxXorder_m) {
        coefficients_m.resize(Xorder + 1);
        for (std::size_t i = (maxXorder_m + 1); i <= Xorder; i++) {
            coefficients_m[i].resize(maxSorder_m + 1, 0);
        }
        maxXorder_m = coefficients_m.size() - 1;
    }
    if (Sorder > maxSorder_m) {
        for (std::size_t i = 0; i < maxXorder_m; i++) {
            coefficients_m[i].resize(Sorder + 1, 0);
        }
        maxSorder_m = coefficients_m[0].size() - 1;
    }
    coefficients_m[Xorder][Sorder] = coefficient;
}
    
void TwoPolynomial::setMaxXorder(const std::size_t &maxXorder) {
    coefficients_m.resize(maxXorder + 1);
    for (std::size_t i = (maxXorder_m + 1); i <= maxXorder; i++) {
        coefficients_m[i].resize(maxSorder_m + 1);
    }
    maxXorder_m = maxXorder;
}

void TwoPolynomial::setMaxSorder(const std::size_t &maxSorder) {
    for (std::size_t i = 0; i <= maxXorder_m; i++) {
        coefficients_m[i].resize(maxSorder + 1);
    }
    maxSorder_m = maxSorder;
}

void TwoPolynomial::setZero() {
    maxXorder_m = 0;
    maxSorder_m = 0;
    coefficients_m.resize(1);
    coefficients_m[0].resize(1);
    coefficients_m[0][0] = 0;
}

double TwoPolynomial::evaluatePolynomial(const double &x,
                                         const double &s) const {
    if (maxXorder_m == 0 && maxSorder_m == 0 &&
        coefficients_m[0][0] == 0) {
        return 0;
    }
    double result = 0.0;
    std::size_t i = maxXorder_m + 1;
    while (i != 0) {
        i--;
        std::size_t j = maxSorder_m + 1;
        double temp = 0.0;
        while (j != 0) {
            j--;
            temp = temp * s + coefficients_m[i][j];
        }
        result = result * x + temp;
    }
    return result;
}

void TwoPolynomial::addPolynomial(const TwoPolynomial &poly) {
    if (maxXorder_m < poly.maxXorder_m) {
        setMaxXorder(poly.maxXorder_m);
    }
    if (maxSorder_m < poly.maxSorder_m) {
        setMaxSorder(poly.maxSorder_m);
    }
    for (std::size_t i = 0; i <= poly.maxXorder_m; i++) {
        for (std::size_t j = 0; j <= poly.maxSorder_m; j++) {
            coefficients_m[i][j] += poly.coefficients_m[i][j];
        }
    }
}

bool operator < (const TwoPolynomial &left, const TwoPolynomial &right) {
    std::size_t nleft = left.dSfactors_m.size();
    std::size_t nright = right.dSfactors_m.size();
    if (nleft < nright) {
        return true;
    } else if (nright < nleft) {
        return false;
    } else {
        for (std::size_t n = 0; n < nleft; n++) {
            std::size_t elemleft = left.dSfactors_m[n];
            std::size_t elemright = right.dSfactors_m[n];
            if (elemleft < elemright) {
                return true;
            } else if (elemright < elemleft) {
                return false;
            }
        }
        return false;
    }
}

bool operator == (const TwoPolynomial &left, const TwoPolynomial &right) {
    if (left < right || right < left) {
        return false;
    } else {
        return true;
    }
}

}
