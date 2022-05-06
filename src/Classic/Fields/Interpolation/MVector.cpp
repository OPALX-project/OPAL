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

#include "Fields/Interpolation/MVector.h"
#include "Fields/Interpolation/MMatrix.h"

namespace interpolation {
///////////////////////// m_complex ////////////////////////////

std::ostream& operator<<(std::ostream& out, m_complex  c)
{
  out << re(c) << " r " << im(c) << " i";
  return out;
}

std::istream& operator>>(std::istream& in,  m_complex& c)
{
  std::string dummy;
  in >> re(c) >> dummy >> im(c) >> dummy;
  return in;
}


///////////////////////// MVector //////////////////////////////

template MVector<double>   ::MVector(size_t i);
template MVector<m_complex>::MVector(size_t i);

template MVector<double>   ::MVector(const MVector<double>&);
template MVector<m_complex>::MVector(const MVector<m_complex>&);
template MVector<double>   ::MVector(size_t i, double    value);
template MVector<m_complex>::MVector(size_t i, m_complex value);

template void MVector<double>::build_vector( const double* data_begin, const double* data_end );
template void MVector<m_complex>::build_vector( const m_complex* data_begin, const m_complex* data_end );

template std::ostream& operator<<(std::ostream& out, MVector<double>    v);
template std::ostream& operator<<(std::ostream& out, MVector<m_complex> v);
template std::istream& operator>>(std::istream& out, MVector<double>&    v);
template std::istream& operator>>(std::istream& out, MVector<m_complex>& v);
template MVector<double>    MVector<double>   ::sub(size_t n1, size_t n2) const;
template MVector<m_complex> MVector<m_complex>::sub(size_t n1, size_t n2) const;


template <typename Tmplt>
MVector<Tmplt>::MVector( size_t i ) : _vector(nullptr)
{
  build_vector(i);
}


template <typename Tmplt>
MVector<Tmplt>::MVector( const MVector<Tmplt>& mv) : _vector(nullptr)
{ *this = mv; }


template <typename Tmplt>
MVector<Tmplt>::MVector( size_t size, Tmplt  value ) : _vector(nullptr)
{
  build_vector(size);
  for(size_t i=0; i<size; i++) operator()(i+1) = value;
}


template <>
void MVector<double>::build_vector   ( size_t size )
{
  if(_vector != nullptr) gsl_vector_free((gsl_vector*)_vector);
  _vector = gsl_vector_alloc(size);
}

template <>
void MVector<m_complex>::build_vector( size_t size )
{
  if(_vector != nullptr) gsl_vector_complex_free((gsl_vector_complex*)_vector);
  _vector = gsl_vector_complex_alloc(size);
}

template <class Tmplt>
void MVector<Tmplt>::build_vector   ( const Tmplt* data_begin, const Tmplt* data_end )
{
  build_vector(data_end - data_begin);
  for(size_t i=0; i<num_row(); i++) operator()(i+1) = data_begin[i];
}


template <class Tmplt>
size_t MVector<Tmplt>::num_row() const
{ if(_vector != nullptr) return ((gsl_vector*)_vector)->size; else return 0;}
template size_t MVector<double>   ::num_row() const;
template size_t MVector<m_complex>::num_row() const;

template <>
const double& MVector<double>::operator()(const size_t i) const
{return *gsl_vector_ptr( (gsl_vector*)_vector, i-1);}
template <>
const m_complex& MVector<m_complex>::operator()(const size_t i) const
{return *gsl_vector_complex_ptr((gsl_vector_complex*)_vector, i-1);}
template <>
double& MVector<double>::operator()(const size_t i)
{return *gsl_vector_ptr( (gsl_vector*)_vector, i-1);}
template <>
m_complex& MVector<m_complex>::operator()(const size_t i)
{return *gsl_vector_complex_ptr((gsl_vector_complex*)_vector, i-1);}


template <class Tmplt>
MMatrix<Tmplt> MVector<Tmplt>::T() const
{
  MMatrix<Tmplt> mat = MMatrix<Tmplt>(1,num_row());
  for(size_t i=1; i<=num_row(); i++)
    mat(1, i) = this->operator()(i);
  return mat;
}
template MMatrix<double>    MVector<double>   ::T() const;
template MMatrix<m_complex> MVector<m_complex>::T() const;

template <class Tmplt> bool operator==(const MVector<Tmplt>& c1, const MVector<Tmplt>& c2)
{
  if(c1.num_row() != c2.num_row())
    return false;
  for(size_t i=0; i<c1.num_row(); i++)
    if(c1(i+1) != c2(i+1)) return false;
  return true;
}
template bool operator==(const MVector<double>&    c1, const MVector<double>& c2);
template bool operator==(const MVector<m_complex>& c1, const MVector<m_complex>& c2);

template <>
MVector<double>& MVector<double>::operator= (const MVector<double>& mv)
{
  if (&mv == this) return *this;
  delete_vector();
  if(!mv._vector) { _vector = nullptr; return *this; }
  _vector = gsl_vector_alloc( mv.num_row() );
  gsl_vector_memcpy((gsl_vector*)_vector, (const gsl_vector*)mv._vector);
  return *this;
}

template <>
MVector<m_complex>& MVector<m_complex>::operator= (const MVector<m_complex>& mv)
{
  if (&mv == this) return *this;
  delete_vector();
  if(!mv._vector) { _vector = nullptr; return *this; }
  _vector = gsl_vector_complex_alloc( mv.num_row() );
  gsl_vector_complex_memcpy((gsl_vector_complex*)_vector, (const gsl_vector_complex*)mv._vector);
  return *this;
}

template <class Tmplt> std::ostream& operator<<(std::ostream& out, MVector<Tmplt> v)
{
  out << v.num_row() << "\n";
  for(size_t i=0; i<v.num_row(); i++) out << "  " << v(i+1) << "\n";
  return out;
}

template <class Tmplt> std::istream& operator>>(std::istream& in, MVector<Tmplt>& v)
{
  size_t n;
  in >> n;
  v = MVector<Tmplt>(n);
  for(size_t i=1; i<=v.num_row(); i++) in >> v(i);
  return in;
}

template <class Tmplt>
MVector<Tmplt> MVector<Tmplt>::sub(size_t n1, size_t n2) const
{
  MVector<Tmplt> temp(n2-n1+1);
  for(size_t i=n1; i<=n2; i++) temp(i-n1+1) = operator()(i);
  return temp;
}

MVector<m_complex> complex(MVector<double> real)
{
  MVector<m_complex> c(real.num_row());
  for(size_t i=1; i<=real.num_row(); i++) c(i) = m_complex_build(real(i));
  return c;
}

MVector<m_complex> complex(MVector<double> real, MVector<double> im)
{
  MVector<m_complex> c(real.num_row());
  for(size_t i=1; i<=real.num_row(); i++) c(i) = m_complex_build(real(i), im(i));
  return c;
}

MVector<double>    re     (MVector<m_complex> c)
{
  MVector<double> d(c.num_row());
  for(size_t i=1; i<=c.num_row(); i++) d(i) = re( c(i) );
  return d;
}

MVector<double>    im     (MVector<m_complex> c)
{
  MVector<double> d(c.num_row());
  for(size_t i=1; i<=c.num_row(); i++) d(i) = im( c(i) );
  return d;
}

}
