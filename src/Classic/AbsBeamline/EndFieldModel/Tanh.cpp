/*
 *  Copyright (c) 2017, Chris Rogers
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

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_pow_int.h"

#include "AbsBeamline/EndFieldModel/Tanh.h"

namespace endfieldmodel {

std::vector< std::vector< std::vector<int> > > Tanh::_tdi;

Tanh::Tanh(double x0, double lambda, int max_index) : _x0(x0), _lambda(lambda) {
  setTanhDiffIndices(max_index);
}

Tanh::~Tanh() {}

EndFieldModel* Tanh::clone() const {
    return new Tanh(*this);
}


double Tanh::getTanh(double x, int n) const {
  if (n == 0) return tanh((x+_x0)/_lambda);
  double t = 0;
  double lam_n = gsl_sf_pow_int(_lambda, n);
  double tanh_x = tanh((x+_x0)/_lambda);
  for (size_t i = 0; i < _tdi[n].size(); i++)
    t += 1./lam_n*static_cast<double>(_tdi[n][i][0])
            *gsl_sf_pow_int(tanh_x, _tdi[n][i][1]);
  return t;
}

double Tanh::getNegTanh(double x, int n) const {
  if (n == 0) return tanh((x-_x0)/_lambda);
  double t = 0;
  double lam_n = gsl_sf_pow_int(_lambda, n);
  double tanh_x = tanh((x-_x0)/_lambda);
  for (size_t i = 0; i < _tdi[n].size(); i++)
    t += 1./lam_n*static_cast<double>(_tdi[n][i][0])
            *gsl_sf_pow_int(tanh_x, _tdi[n][i][1]);
  return t;
}

double Tanh::function(double x, int n) const {
  return (getTanh(x, n)-getNegTanh(x, n))/2.;
}

void Tanh::setTanhDiffIndices(size_t n) {
  if (_tdi.size() == 0) {
    _tdi.push_back(std::vector< std::vector<int> >(1, std::vector<int>(2)));
    _tdi[0][0][0] = 1;  // 1*tanh(x) - third index is redundant
    _tdi[0][0][1] = 1;
  }
  for (size_t i = _tdi.size(); i < n+1; ++i) {
    _tdi.push_back(std::vector< std::vector<int> >());
    for (size_t j = 0; j < _tdi[i-1].size(); ++j) {
      if (_tdi[i-1][j][1] != 0) {
        std::vector<int> new_vec(_tdi[i-1][j]);
        new_vec[0] *= _tdi[i-1][j][1];
        new_vec[1] -= 1;
        _tdi[i].push_back(new_vec);
      }
      if (_tdi[i-1][j][1] != 0) {
        std::vector<int> new_vec(_tdi[i-1][j]);
        new_vec[0] *= -_tdi[i-1][j][1];
        new_vec[1] += 1;
        _tdi[i].push_back(new_vec);
      }
    }
    _tdi[i] = CompactVector(_tdi[i]);
  }
}

std::vector< std::vector<int> > Tanh::getTanhDiffIndices(size_t n) {
  setTanhDiffIndices(n);
  return _tdi[n];
}


std::ostream& Tanh::print(std::ostream& out) const {
    out << "Tanh model with centre length: " << _x0 
        << " end length: " << _lambda;
    return out;
}
}

