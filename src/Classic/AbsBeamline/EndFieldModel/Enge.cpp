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

#include "AbsBeamline/EndFieldModel/Enge.h"

namespace endfieldmodel {

// Use
// d^n E/dx^n = a_n1m1 F(n1) g(m1) + a_n2m1m2 F(n2) g(m1)g(m2)+...
// where 
double Enge::GetEnge(double x, int n) const {
  std::vector< std::vector<int> > qt = GetQIndex(n);
  std::vector<double> g;
  double e(0.);
  for (size_t i = 0; i < qt.size(); ++i) {
    double ei(qt[i][0]);
    for (size_t j = 1; j < qt[i].size(); ++j) {
      if (j > g.size()) g.push_back(GN(x, j-1));
      ei *= gsl_sf_pow_int(g[j-1], qt[i][j]);
    }
    if (ei != ei) ei = 0;  // div 0, usually g == 0 and index < 0
    e += ei;
  }
  return e;
}


// h     = a_0+a_1 (x/w)+a_2 (x/w)^2+a_3 (x/w)^3+...+a_m (x/w)^m
// h^(n) = d^nh/dx^n = sum^m_{i=n} a_i x^{i-n}/w^i i!/n!
double Enge::HN(double x, int n) const {
  double hn   = 0;
  // optimise by precalculating factor
  for (unsigned int i = n; i <_a.size(); i++)
    hn += _a[i]/gsl_sf_pow_int(_lambda, i)*gsl_sf_pow_int(x, i-n)*
          static_cast<double>(gsl_sf_fact(i))
         /static_cast<double>(gsl_sf_fact(i-n));
  return hn;
}

// g     = 1+exp(h)
// g^(n) = d^ng/dx^n
double Enge::GN(double x, int n) const {
  if (n == 0) return 1+exp(HN(x, 0));  // special case
  std::vector<double> hn(n+1);
  for (int i = 0; i <= n; i++) hn[i] = HN(x, i);
  double exp_h0 = exp(hn[0]);
  double gn = 0;
  for (size_t i = 0; i < _h[n].size(); ++i) {
    double gnj = _h[n][i][0]*exp_h0;
    for (size_t j = 1; j < _h[n][i].size(); ++j)
      gnj *= gsl_sf_pow_int(hn[j], _h[n][i][j]);
    gn += gnj;
  }
  return gn;
}

// _q[i][j][k]; urk, 3d vector
//              i indexes the derivative of f;
//              j indexes the element in f derivative
//              k indexes the derivative of g
// this will quickly become grotesque
std::vector< std::vector< std::vector<int> > > Enge::_q;
std::vector< std::vector< std::vector<int> > > Enge::_h;
void Enge::SetEngeDiffIndices(size_t n) {
  if (_q.size() == 0) {
    _q.push_back(std::vector< std::vector<int> >(1, std::vector<int>(3)) );
    _q[0][0][0] = +1;  // f_0 = 1*g^(-1)
    _q[0][0][1] = -1;
    _q[0][0][2] =  0;
  }

  for (size_t i = _q.size(); i < n+1; ++i) {
    _q.push_back(std::vector< std::vector<int> >() );
    for (size_t j = 0; j < _q[i-1].size(); ++j) {
      size_t k_max = _q[i-1][j].size();
      std::vector<int> new_vec(_q[i-1][j]);
      // derivative of g^-n0 = -n0*g^(-n0-1)*g(1)
      new_vec[0] *= new_vec[1];  //  alpha *= g(0) power
      new_vec[1] -= 1;  // g(0) power -= 1
      new_vec[2] += 1;  // g(1) power += 1
      _q[i].push_back(new_vec);
      for (size_t k = 2; k < k_max; ++k) {  //  0 is alpha; 1 is g(0)
        // derivative of g(k)^nk = nk g(k+1) g(k)^(nk-1)
        if (_q[i-1][j][k] > 0) {
          std::vector<int> new_vec(_q[i-1][j]);
          if ( k == k_max-1 ) new_vec.push_back(0);  // need enough coefficients
          new_vec[0]   *= new_vec[k];
          new_vec[k]   -= 1;
          new_vec[k+1] += 1;
          _q[i].push_back(new_vec);
        }
      }
    }
  }

  if (_h.size() == 0) {
     // first one is special case (1+e^h dealt with explicitly)
    _h.push_back(std::vector< std::vector<int> >());
     // second is (1*e^h h'^1)
    _h.push_back(std::vector< std::vector<int> >());
    _h[1].push_back(std::vector<int>(2, 1));
  }
  for (size_t i = _h.size(); i < n+1; ++i) {
    _h.push_back(std::vector< std::vector<int> >());
    for (size_t j = 0; j < _h[i-1].size(); ++j) {
      // d/dx k0 e^g g(1)^k1 ... g(n)^kn ... = k0 e^g g(1)^(k1+1) ... g(n)^kn
      //                              + SUM_n k0 kn e^g ... g(n)^(kn-1) g(n-1)
      std::vector<int> new_vec(_h[i-1][j]);
      new_vec[1] += 1;
      _h[i].push_back(new_vec);
      for (size_t k = 1; k <_h[i-1][j].size(); ++k) {
        if (_h[i-1][j][k] > 0) {
          std::vector<int> new_vec(_h[i-1][j]);
          if (k == _h[i-1][j].size()-1) new_vec.push_back(0);
          new_vec[0]   *= new_vec[k];
          new_vec[k]   -= 1;
          new_vec[k+1] += 1;
          _h[i].push_back(new_vec);
        }
      }
    }
    _h[i] = CompactVector(_h[i]);
  }
}

}

