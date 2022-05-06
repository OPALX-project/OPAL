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

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"

#include "Fields/Interpolation/MMatrix.h"

///////////////////////////////// MMATRIX ///////////////////////////////

//template class MMatrix<double>;

/////////// CONSTRUCTORS, DESTRUCTORS
namespace interpolation {

template <class Tmplt>
MMatrix<Tmplt>::MMatrix() : _matrix(nullptr)
{}
template MMatrix<double>   ::MMatrix();
template MMatrix<m_complex>::MMatrix();

template MMatrix<double>    MMatrix<double>   ::Diagonal(size_t i, double    diag, double    off_diag);
template MMatrix<m_complex> MMatrix<m_complex>::Diagonal(size_t i, m_complex diag, m_complex off_diag);

template std::istream& operator>>(std::istream& in, MMatrix<double>&    mat);
template std::istream& operator>>(std::istream& in, MMatrix<m_complex>& mat);

template MMatrix<double>   ::MMatrix(size_t i, size_t j, double*    data_beg );
template MMatrix<m_complex>::MMatrix(size_t i, size_t j, m_complex* data_beg );
template MMatrix<double>   ::MMatrix(size_t i, size_t j, double    value);
template MMatrix<m_complex>::MMatrix(size_t i, size_t j, m_complex value);
template MMatrix<double>   ::MMatrix(size_t i, size_t j );
template MMatrix<m_complex>::MMatrix(size_t i, size_t j );

template MMatrix<double>    MMatrix<double>   ::inverse() const;
template MMatrix<m_complex> MMatrix<m_complex>::inverse() const;

template <>
void MMatrix<double>::delete_matrix()
{
  if(_matrix != nullptr) gsl_matrix_free( (gsl_matrix*)_matrix );
  _matrix = nullptr;
}

template <>
void MMatrix<m_complex>::delete_matrix()
{
  if(_matrix != nullptr) gsl_matrix_complex_free( (gsl_matrix_complex*)_matrix );
  _matrix = nullptr;
}


template <>
MMatrix<double>& MMatrix<double>::operator= (const MMatrix<double>& mm)
{
  if (&mm == this) return *this;
  delete_matrix();
  if(!mm._matrix) { _matrix = nullptr; return *this; }
  _matrix = gsl_matrix_alloc( mm.num_row(), mm.num_col() );
  gsl_matrix_memcpy((gsl_matrix*)_matrix, (const gsl_matrix*)mm._matrix);
  return *this;
}

template <>
MMatrix<m_complex>& MMatrix<m_complex>::operator= (const MMatrix<m_complex>& mm)
{
  if (&mm == this) return *this;
  delete_matrix();
  if(!mm._matrix) { _matrix = nullptr; return *this; }
  _matrix = gsl_matrix_complex_alloc( mm.num_row(), mm.num_col() );
  gsl_matrix_complex_memcpy((gsl_matrix_complex*)_matrix, (const gsl_matrix_complex*)mm._matrix);
  return *this;
}

template <>
MMatrix<double>::MMatrix( const MMatrix<double>& mm ) : _matrix(nullptr)
{
  if(mm._matrix)
  {
    _matrix = gsl_matrix_alloc(mm.num_row(), mm.num_col() );
    gsl_matrix_memcpy( (gsl_matrix*)_matrix, (gsl_matrix*)mm._matrix );
  }
}

template <>
MMatrix<m_complex>::MMatrix( const MMatrix<m_complex>& mm ) : _matrix(nullptr)
{
  if(mm._matrix)
  {
    _matrix = gsl_matrix_complex_alloc( mm.num_row(), mm.num_col() );
    gsl_matrix_complex_memcpy( (gsl_matrix_complex*)_matrix, (gsl_matrix_complex*)mm._matrix );
  }
}

template <class Tmplt>
MMatrix<Tmplt>::MMatrix(size_t i, size_t j, Tmplt* data_beg ) : _matrix(nullptr)
{
  build_matrix(i, j, data_beg);
}

template <class Tmplt>
MMatrix<Tmplt>::MMatrix(size_t i, size_t j, Tmplt  value    )
{
    build_matrix(i, j);
    for(size_t a=1; a<=i; a++)
        for(size_t b=1; b<=j; b++)
            operator()(a,b) = value;
}

template <class Tmplt>
MMatrix<Tmplt>::MMatrix(size_t i, size_t j )
{
  build_matrix(i, j);
}

template <class Tmplt>
MMatrix<Tmplt> MMatrix<Tmplt>::Diagonal(size_t i, Tmplt diag_value, Tmplt off_diag_value)
{
  MMatrix<Tmplt> mm(i,i);
  for(size_t a=1; a<=i; a++)
  {
    for(size_t b=1;   b< a; b++) mm(a, b) = off_diag_value;
    for(size_t b=a+1; b<=i; b++) mm(a, b) = off_diag_value;
    mm(a,a) = diag_value;
  }
  return mm;
}


template <class Tmplt>
MMatrix<Tmplt>::~MMatrix()
{
  delete_matrix();
}
template MMatrix<double>::~MMatrix();
template MMatrix<m_complex>::~MMatrix();

template <>
void MMatrix<double>::build_matrix(size_t i, size_t j)
{
  _matrix = gsl_matrix_alloc(i, j);
}

template <>
void MMatrix<m_complex>::build_matrix(size_t i, size_t j)
{
  _matrix = gsl_matrix_complex_alloc(i, j);
}

template <>
void MMatrix<double>::build_matrix(size_t i, size_t j, double* data)
{
  _matrix = gsl_matrix_alloc(i, j);
  for(size_t a=0; a<i; a++)
    for(size_t b=0; b<j; b++)
      operator()(a+1,b+1) = data[b*i + a];
}

template <>
void MMatrix<m_complex>::build_matrix(size_t i, size_t j, m_complex* data)
{
  _matrix = gsl_matrix_complex_alloc(i, j);
  for(size_t a=0; a<i; a++)
    for(size_t b=0; b<j; b++)
     operator()(a+1,b+1) = data[b*i + a];
}


////////////////// HIGH LEVEL FUNCTIONS

template <>
m_complex MMatrix<m_complex>::determinant() const
{
  int signum = 1;

  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix::determinant()", "Attempt to get determinant of non-square matrix"));
  gsl_permutation * p = gsl_permutation_alloc (num_row());
  MMatrix<m_complex> copy(*this);
  gsl_linalg_complex_LU_decomp ((gsl_matrix_complex*)copy._matrix, p, &signum);
  gsl_permutation_free(p);
  return gsl_linalg_complex_LU_det((gsl_matrix_complex*)copy._matrix, signum);
}

template <>
double MMatrix<double>::determinant() const
{
  int signum = 1;
  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix::determinant()", "Attempt to get determinant of non-square matrix"));
  gsl_permutation * p = gsl_permutation_alloc (num_row());
  MMatrix<double> copy(*this);
  gsl_linalg_LU_decomp ((gsl_matrix*)copy._matrix, p, &signum);
  double d = gsl_linalg_LU_det((gsl_matrix*)copy._matrix,  signum);
  gsl_permutation_free(p);
  return d;
}

template <class Tmplt>
MMatrix<Tmplt> MMatrix<Tmplt>::inverse()     const
{
  MMatrix<Tmplt> copy(*this);
  copy.invert();
  return copy;
}


template <>
void MMatrix<m_complex>::invert()
{
  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix::invert()", "Attempt to get inverse of non-square matrix"));
  gsl_permutation* perm = gsl_permutation_alloc(num_row());
  MMatrix<m_complex>  lu(*this); //hold LU decomposition
  int s=0;
  gsl_linalg_complex_LU_decomp ( (gsl_matrix_complex*)lu._matrix, perm, &s);
  gsl_linalg_complex_LU_invert ( (gsl_matrix_complex*)lu._matrix, perm, (gsl_matrix_complex*)_matrix);
  gsl_permutation_free( perm );
  for(size_t i=1; i<=num_row(); i++)
    for(size_t j=1; j<=num_row(); j++)
      if(operator()(i,j) != operator()(i,j)) throw(GeneralClassicException("MMatrix::invert()", "Failed to invert matrix - singular?"));
}

template <>
void MMatrix<double>::invert()
{
  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix::invert()", "Attempt to get inverse of non-square matrix"));
  gsl_permutation* perm = gsl_permutation_alloc(num_row());
  MMatrix<double>  lu(*this); //hold LU decomposition
  int s=0;
  gsl_linalg_LU_decomp ( (gsl_matrix*)lu._matrix, perm, &s);
  gsl_linalg_LU_invert ( (gsl_matrix*)lu._matrix, perm, (gsl_matrix*)_matrix);
  gsl_permutation_free( perm );
  for(size_t i=1; i<=num_row(); i++)
    for(size_t j=1; j<=num_row(); j++)
      if(operator()(i,j) != operator()(i,j)) throw(GeneralClassicException("MMatrix::invert()", "Failed to invert matrix - singular?"));
}

template <>
MMatrix<double> MMatrix<double>::T()     const
{
  MMatrix<double>  in(num_col(), num_row());
  gsl_matrix_transpose_memcpy ((gsl_matrix*)in._matrix, (gsl_matrix*)_matrix);
  return in;
}

template <>
MMatrix<m_complex> MMatrix<m_complex>::T()     const
{
  MMatrix<m_complex>  in(num_col(), num_row());
  gsl_matrix_complex_transpose_memcpy ((gsl_matrix_complex*)in._matrix, (gsl_matrix_complex*)_matrix);
  return in;
}

template <class Tmplt>
MMatrix<Tmplt>  MMatrix<Tmplt>::sub(size_t min_row, size_t max_row, size_t min_col, size_t max_col) const
{
  MMatrix<Tmplt> sub_matrix(max_row-min_row+1, max_col-min_col+1);
  for(size_t i=min_row; i<=max_row; i++)
    for(size_t j=min_col; j<=max_col; j++)
      sub_matrix(i-min_row+1, j-min_col+1) = operator()(i,j);
  return sub_matrix;
}
template MMatrix<double>    MMatrix<double>   ::sub(size_t min_row, size_t max_row, size_t min_col, size_t max_col) const;
template MMatrix<m_complex> MMatrix<m_complex>::sub(size_t min_row, size_t max_row, size_t min_col, size_t max_col) const;


template <class Tmplt>
Tmplt MMatrix<Tmplt>::trace() const
{
  Tmplt t = operator()(1,1);
  for(size_t i=2; i<=num_row() && i<=num_col(); i++) t+= operator()(i,i);
  return t;
}
template double    MMatrix<double>::trace()    const;
template m_complex MMatrix<m_complex>::trace() const;


template <class Tmplt>
MVector<m_complex> MMatrix<Tmplt>::eigenvalues() const
{
  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix<double>::eigenvalues", "Attempt to get eigenvectors of non-square matrix") );
  MMatrix<Tmplt> temp = *this;
  MVector<m_complex> evals(num_row(), m_complex_build(2.,-1.));
  gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc(num_row());
  gsl_eigen_nonsymm_params(0, 1, w);
  int ierr = gsl_eigen_nonsymm(get_matrix(temp), evals.get_vector(evals), w);
  gsl_eigen_nonsymm_free(w);
  if(ierr != 0) throw(GeneralClassicException("MMatrix<Tmplt>::eigenvalues", "Failed to calculate eigenvalue"));
  return evals;
}
template MVector<m_complex> MMatrix<double>   ::eigenvalues() const;

template <class Tmplt>
std::pair<MVector<m_complex>, MMatrix<m_complex> > MMatrix<Tmplt>::eigenvectors() const
{
  if(num_row() != num_col()) throw(GeneralClassicException("MMatrix<double>::eigenvectors", "Attempt to get eigenvectors of non-square matrix") );
  MMatrix<Tmplt>     temp = *this;
  MVector<m_complex> evals(num_row());
  MMatrix<m_complex> evecs(num_row(), num_row());
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(num_row());
  int ierr = gsl_eigen_nonsymmv(get_matrix(temp), evals.get_vector(evals), get_matrix(evecs), w);
  gsl_eigen_nonsymmv_free(w);
  if(ierr != 0) throw(GeneralClassicException("MMatrix<Tmplt>::eigenvectors", "Failed to calculate eigenvalue"));
  return std::pair<MVector<m_complex>, MMatrix<m_complex> > (evals, evecs) ;
}
template std::pair<MVector<m_complex>, MMatrix<m_complex> > MMatrix<double>::eigenvectors() const;

///////////// OPERATORS

MMatrix<double>& operator *= (MMatrix<double>& m1,  MMatrix<double>  m2)
{
  MMatrix<double> out(m1.num_row(), m2.num_col());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., (gsl_matrix*)m1._matrix, (gsl_matrix*)m2._matrix, 0., (gsl_matrix*)out._matrix);
  m1 = out;
  return m1;
}

MMatrix<m_complex> & operator *= (MMatrix<m_complex>& m1,  MMatrix<m_complex>  m2)
{
  MMatrix<m_complex> out(m1.num_row(), m2.num_col());
  gsl_blas_zgemm( CblasNoTrans, CblasNoTrans, m_complex_build(1.), (gsl_matrix_complex*)m1._matrix,
                  (gsl_matrix_complex*)m2._matrix, m_complex_build(0.), (gsl_matrix_complex*)out._matrix);
  m1 = out;
  return m1;
}

template <class Tmplt> std::ostream& operator<<(std::ostream& out, MMatrix<Tmplt> mat)
{
  out << mat.num_row() << " " << mat.num_col() << "\n";
  for(size_t i=1; i<=mat.num_row(); i++)
  {
    out << "  ";
    for(size_t j=1; j<mat.num_col(); j++)
      out << mat(i,j) << "   ";
    out << mat(i,mat.num_col()) << "\n";
  }
  return out;
}
template std::ostream& operator<<(std::ostream& out, MMatrix<double>    mat);
template std::ostream& operator<<(std::ostream& out, MMatrix<m_complex> mat);

template <class Tmplt> std::istream& operator>>(std::istream& in, MMatrix<Tmplt>& mat)
{
  size_t nr, nc;
  in >> nr >> nc;
  mat = MMatrix<Tmplt>(nr, nc);
  for(size_t i=1; i<=nr; i++)
    for(size_t j=1; j<=nc; j++)
      in >> mat(i,j);
  return in;
}


///////////////// INTERFACES

MMatrix<double>    re(MMatrix<m_complex> mc)
{
  MMatrix<double> md(mc.num_row(), mc.num_col());
  for(size_t i=1; i<=mc.num_row(); i++)
    for(size_t j=1; j<=mc.num_col(); j++)
      md(i,j) = re(mc(i,j));
  return md;
}

MMatrix<double>    im(MMatrix<m_complex> mc)
{
  MMatrix<double> md(mc.num_row(), mc.num_col());
  for(size_t i=1; i<=mc.num_row(); i++)
    for(size_t j=1; j<=mc.num_col(); j++)
      md(i,j) = im(mc(i,j));
  return md;
}

MMatrix<m_complex> complex(MMatrix<double> real)
{
  MMatrix<m_complex> mc(real.num_row(), real.num_col());
  for(size_t i=1; i<=mc.num_row(); i++)
    for(size_t j=1; j<=mc.num_col(); j++)
      mc(i,j) = m_complex_build(real(i,j));
  return mc;
}

MMatrix<m_complex> complex(MMatrix<double> real, MMatrix<double> imaginary)
{
  if(real.num_row() != imaginary.num_row() || real.num_col() != imaginary.num_col())
    throw(GeneralClassicException("MMatrix<m_complex>::complex", "Attempt to build complex matrix when real and imaginary matrix don't have the same size"));
  MMatrix<m_complex> mc(real.num_row(), real.num_col());
  for(size_t i=1; i<=mc.num_row(); i++)
    for(size_t j=1; j<=mc.num_col(); j++)
      mc(i,j) = m_complex_build(real(i,j), imaginary(i,j));
  return mc;
}

template <class Tmplt>
MVector<Tmplt> MMatrix<Tmplt>::get_mvector(size_t column) const
{
  MVector<Tmplt> temp(num_row());
  for(size_t i=1; i<=num_row(); i++) temp(i) = operator()(i,column);
  return temp;
}
template MVector<double> MMatrix<double>::get_mvector(size_t column) const;
template MVector<m_complex> MMatrix<m_complex>::get_mvector(size_t column) const;

}
