/*
 *  Copyright (c) 2015, Chris Rogers
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


#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "gsl/gsl_sf_gamma.h"

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/MMatrix.h"
#include "Fields/Interpolation/MVector.h"
#include "Fields/Interpolation/SquarePolynomialVector.h"

namespace interpolation {

bool SquarePolynomialVector::_printHeaders=true;
std::vector< std::vector< std::vector<int> > > SquarePolynomialVector::_polyKeyByPower;
std::vector< std::vector< std::vector<int> > > SquarePolynomialVector::_polyKeyByVector;

SquarePolynomialVector::SquarePolynomialVector()
    : _pointDim(0), _polyCoeffs() {
}

SquarePolynomialVector::SquarePolynomialVector (const SquarePolynomialVector& pv)
    : _pointDim(pv._pointDim),
      _polyCoeffs(pv._polyCoeffs.num_row(), pv._polyCoeffs.num_col(), 0.) {
    SetCoefficients(pv._polyCoeffs);
}


SquarePolynomialVector::SquarePolynomialVector(int numberOfInputVariables, MMatrix<double> polynomialCoefficients)
        :  _pointDim(numberOfInputVariables), _polyCoeffs(polynomialCoefficients) {
    SetCoefficients(numberOfInputVariables, polynomialCoefficients);
}

SquarePolynomialVector::SquarePolynomialVector(std::vector<PolynomialCoefficient> coefficients)
        :  _pointDim(0), _polyCoeffs() {
    SetCoefficients(coefficients);
}

void SquarePolynomialVector::SetCoefficients(int pointDim, MMatrix<double> coeff) {
    _pointDim   = pointDim;
    _polyCoeffs = coeff;

    if (int(_polyKeyByVector.size()) < pointDim || _polyKeyByVector[pointDim-1].size() < coeff.num_col())
        IndexByVector(coeff.num_col(), pointDim); // sets _polyKeyByVector and _polyKeyByPower
}

void SquarePolynomialVector::SetCoefficients(std::vector<PolynomialCoefficient> coeff) {
  _pointDim        = 0;
  int maxPolyOrder = 0;
  int valueDim     = 0;
  for(unsigned int i=0; i<coeff.size(); i++) {
    int polyOrder = coeff[i].InVariables().size();
    for(unsigned int j=0; j<coeff[i].InVariables().size(); j++)
      if(coeff[i].InVariables()[j] > _pointDim) _pointDim = coeff[i].InVariables()[j];
    if(coeff[i].OutVariable()    > valueDim)  valueDim  = coeff[i].OutVariable();
    if(polyOrder > maxPolyOrder) maxPolyOrder = polyOrder;
  }
  _pointDim++; //PointDim indexes from 0
  _polyCoeffs = MMatrix<double>(valueDim+1, NumberOfPolynomialCoefficients(_pointDim, maxPolyOrder+1), 0.);
  if (int(_polyKeyByVector.size()) < _pointDim || _polyKeyByVector[_pointDim-1].size() < _polyCoeffs.num_col())
      IndexByVector(_polyCoeffs.num_col(), _pointDim); // sets _polyKeyByVector and _polyKeyByPower

  for(size_t i=0; i<_polyCoeffs.num_col(); i++) {
      for(unsigned int j=0; j<coeff.size(); j++)
          if(coeff[j].InVariables() == _polyKeyByVector[_pointDim-1].back()) {
            int dim = coeff[j].OutVariable();
            _polyCoeffs(dim+1,i+1) = coeff[j].Coefficient();
          }
  }
}

void SquarePolynomialVector::SetCoefficients(MMatrix<double> coeff)
{
  for(size_t i=0; i<coeff.num_row() && i<_polyCoeffs.num_row(); i++)
    for(size_t j=0; j<coeff.num_col() && j<_polyCoeffs.num_col(); j++)
      _polyCoeffs(i+1,j+1) = coeff(i+1,j+1);
}

void  SquarePolynomialVector::F(const double*   point,    double* value)          const
{
    MVector<double> pointV(PointDimension(), 1), valueV(ValueDimension(), 1);
    for(unsigned int i=0; i<PointDimension(); i++) pointV(i+1) = point[i];
    F(pointV, valueV);
    for(unsigned int i=0; i<ValueDimension(); i++) value[i]  = valueV(i+1);
}

void  SquarePolynomialVector::F(const MVector<double>& point,   MVector<double>& value) const
{
    MVector<double> polyVector(_polyCoeffs.num_col(), 1);
    MakePolyVector(point, polyVector);
    MVector<double> answer = _polyCoeffs*polyVector;
    for(unsigned int i=0; i<ValueDimension(); i++) value(i+1) = answer(i+1);
}


MVector<double>& SquarePolynomialVector::MakePolyVector(const MVector<double>& point, MVector<double>& polyVector) const
{
    for(unsigned int i=0; i<_polyCoeffs.num_col(); i++) {
        polyVector(i+1) = 1.;
        for(unsigned int j=0; j<_polyKeyByVector[_pointDim-1][i].size(); j++)
            polyVector(i+1) *= point( _polyKeyByVector[_pointDim-1][i][j]+1 );
    }
    return polyVector;
}

double* SquarePolynomialVector::MakePolyVector(const double* point, double* polyVector) const
{
    for(unsigned int i=0; i<_polyCoeffs.num_col(); i++)
    {
        polyVector[i] = 1.;
        for(unsigned int j=0; j<_polyKeyByVector[i].size(); j++)
            polyVector[i] *= point[_polyKeyByVector[_pointDim-1][i][j] ];
    }
    return polyVector;
}


void SquarePolynomialVector::IndexByPowerRecursive(std::vector<int> check, size_t check_index, size_t poly_power, std::vector<std::vector<int> >& nearby_points) {
    check[check_index] = poly_power;
    nearby_points.push_back(check);
    if (check_index+1 == check.size())
        return;
    for (int i = 1; i < int(poly_power); ++i)
        IndexByPowerRecursive(check, check_index+1, i, nearby_points);
    for (int i = 0; i <= int(poly_power); ++i) {
        check[check_index] = i;
        IndexByPowerRecursive(check, check_index+1, poly_power, nearby_points);
    }
}

std::vector<int> SquarePolynomialVector::IndexByPower(int index, int point_dim) {
    if (point_dim < 1)
        throw(GeneralClassicException(
            "PPSolveFactory::GetNearbyPointsSquares",
            "Point dimension must be > 0"
        ));
    while (int(_polyKeyByPower.size()) < point_dim)
        _polyKeyByPower.push_back(std::vector< std::vector<int> >());
    if (index < int(_polyKeyByPower[point_dim-1].size()))
        return _polyKeyByPower[point_dim-1][index];
    // poly_order 0 means constant term
    std::vector<std::vector<int> > nearby_points(1, std::vector<int>(point_dim, 0));
    // poly_order > 0 means corresponding poly terms (1 is linear, 2 is quadratic...)
    int this_poly_order = 0;
    while (int(nearby_points.size()) < index+1) {
        IndexByPowerRecursive(nearby_points[0], 0, this_poly_order+1, nearby_points);
        this_poly_order += 1;
    }
    _polyKeyByPower[point_dim-1] = nearby_points;
    return _polyKeyByPower[point_dim-1][index];
}

//Turn int index into an index for element of polynomial
// e.g. x_1^4 x_2^3 = {1,1,1,1,2,2,2}
std::vector<int> SquarePolynomialVector::IndexByVector(int index, int point_dim) {
    while (int(_polyKeyByVector.size()) < point_dim)
        _polyKeyByVector.push_back(std::vector< std::vector<int> >());
    // if _polyKeyByVector is big enough, just return the value
    if (index < int(_polyKeyByVector[point_dim-1].size()))
        return _polyKeyByVector[point_dim-1][index];
    // make sure _polyKeyByPower is big enough
    IndexByPower(index, point_dim);
    // update _polyKeyByVector with values from _polyKeyByPower
    for (size_t i = _polyKeyByVector[point_dim-1].size(); i < _polyKeyByPower[point_dim-1].size(); ++i) {
        _polyKeyByVector[point_dim-1].push_back(std::vector<int>());
        for (size_t j = 0; j < _polyKeyByPower[point_dim-1][i].size(); ++j)
            for (int k = 0; k < _polyKeyByPower[point_dim-1][i][j]; ++k)
                _polyKeyByVector[point_dim-1][i].push_back(j);
    }
    return _polyKeyByVector[point_dim-1][index];
}

unsigned int SquarePolynomialVector::NumberOfPolynomialCoefficients(int pointDimension, int order) {
    int number = 1;
    for (int i = 0; i < pointDimension; ++i)
        number *= order+1;
    return number;
}

std::ostream& operator<<(std::ostream& out, const SquarePolynomialVector& spv)
{
    if(SquarePolynomialVector::_printHeaders) spv.PrintHeader(out, '.', ' ', 14, true);
    out << '\n' << spv.GetCoefficientsAsMatrix();
    return out;
}

void SquarePolynomialVector::PrintHeader(std::ostream& out, char int_separator, char str_separator, int length, bool pad_at_start) const
{
  if(_polyKeyByPower[_pointDim-1].size() > 0) PrintContainer< std::vector<int> >(out, _polyKeyByPower[_pointDim-1][0], int_separator, str_separator, length-1, pad_at_start);
  for(unsigned int i=1; i<_polyCoeffs.num_col(); ++i)
    PrintContainer<std::vector<int> >(out, _polyKeyByPower[_pointDim-1][i], int_separator, str_separator, length, pad_at_start);
}


template <class Container >
void SquarePolynomialVector::PrintContainer(std::ostream& out, const Container& container, char T_separator, char str_separator, int length, bool pad_at_start)
{
//  class Container::iterator it;
  std::stringstream strstr1("");
  std::stringstream strstr2("");
  typename Container::const_iterator it1 = container.begin();
  typename Container::const_iterator it2 = it1;
  while(it1!=container.end())
  {
    it2++;
    if(it1 != container.end() && it2 != container.end())
         strstr1 << (*it1) << T_separator;
    else strstr1 << (*it1);
    it1++;
  }

  if(int(strstr1.str().size()) > length && length > 0) strstr2 << strstr1.str().substr(1, length);
  else                                                 strstr2 << strstr1.str();

  if(!pad_at_start) out << strstr2.str();
  for(int i=0; i<length - int(strstr2.str().size()); i++) out << str_separator;
  if(pad_at_start) out << strstr2.str();
}

unsigned int SquarePolynomialVector::PolynomialOrder() const {
    std::vector<int> index = IndexByPower(_polyCoeffs.num_col(), _pointDim);
    size_t maxPower = *std::max_element(index.begin(), index.end());
    return maxPower;
}


}