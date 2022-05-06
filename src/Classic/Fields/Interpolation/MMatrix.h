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

#ifndef _SRC_COMMON_CPP_MATHS_MMATRIX_HH_
#define _SRC_COMMON_CPP_MATHS_MMATRIX_HH_

#include <iostream>
#include <vector>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/MVector.h"

namespace interpolation {

///////////////// MMatrix     /////////////////

/** \class MMatrix
 * 
 *  C++ wrapper for GSL matrix
 *  MMatrix class handles matrix algebra, maths operators and some
 *  higher level calculation like matrix inversion, eigenvector
 *  analysis etc
 *
 *  Use template to define two types:\n
 *  (i) MMatrix\<double\> is a matrix of doubles\n
 *  (ii) MMatrix\<m_complex\> is a matrix of m_complex\n
 *
 *  Maths operators and a few others are defined, but maths operators
 *  don't allow operations between types - i.e. you can't multiply
 *  a complex matrix by a double matrix. Instead use interface methods
 *  like real() and complex() to convert between types first
 */
template <class Tmplt>
class MMatrix
{
public:
  /** @brief default constructor makes an empty MMatrix of size (0,0)
   */
  MMatrix();

  /** @brief Copy constructor makes a deep copy of mv
   */
  MMatrix(const MMatrix<Tmplt>& mv );

  /** @brief Construct a matrix and fill with data from memory data_beginning
   *
   *  @param nrows number of rows
   *  @param ncols number of columns
   *  @param data_beginning pointer to the start of a memory block of size
   *          nrows*ncols with data laid out \<row\> \<row\> \<row\>. Note
   *          MMatrix does not take ownership of memory in data_beginning.
   */
  MMatrix(size_t nrows, size_t ncols, Tmplt* data_beginning );

  /** @brief Construct a matrix and fill with identical data
   *
   *  @param nrows number of rows
   *  @param ncols number of columns
   *  @param value variable to be copied into all items in the matrix
   */
  MMatrix(size_t nrows, size_t ncols, Tmplt  value);

  /** @brief Construct a matrix and fill all fields with 0
   *
   *  @param nrows number of rows
   *  @param ncols number of columns
   */
  MMatrix(size_t nrows, size_t ncols);

  /** @brief Construct a square matrix filling on and off diagonal values
   *
   *  @param i number of rows and number of columns
   *  @param diag_value fill values with row == column (i.e. on the diagonal)
   *          with this value
   *  @param off_diag_value fill values with row != column (i.e. off the
   *          diagonal) with this value
   */
  static MMatrix<Tmplt>  Diagonal(size_t i, Tmplt diag_value, Tmplt off_diag_value);

  /** @brief destructor
   */
  ~MMatrix();
 
  /** @brief returns number of rows in the matrix
   */
  size_t           num_row()     const;

  /** @brief returns number of columns in the matrix
   */
  size_t           num_col()     const;

  /** @brief returns sum of terms with row == column, even if matrix is not
   *         square
   */
  Tmplt            trace()       const;

  /** @brief returns matrix determinant; throws an exception if matrix is not
   *         square
   */
  Tmplt            determinant() const;

  /** @brief returns matrix inverse leaving this matrix unchanged
   */
  MMatrix<Tmplt>   inverse()     const;

  /** @brief turns this matrix into its inverse
   */
  void             invert();

  /** @brief returns matrix transpose T (such that M(i,j) = T(j,i))
   */
  MMatrix<Tmplt>     T()            const;

  /** @brief returns a vector of eigenvalues. Throws an exception if either this
   *         matrix is not square or the eigenvalues could not be found (e.g.
   *         singular matrix or whatever).
   */
  MVector<m_complex> eigenvalues() const; //return vector of eigenvalues
  std::pair<MVector<m_complex>, MMatrix<m_complex> > 
                     eigenvectors() const; //return pair of eigenvalues, eigenvectors (access using my_pair.first or my_pair.second)

  //return a submatrix, with data from min_row to max_row and min_row to max_row inclusive; indexing starts at 1
  MMatrix<Tmplt>        sub(size_t min_row, size_t max_row, size_t min_col, size_t max_col) const;
  //create a column vector from the column^th column
  MVector<Tmplt>        get_mvector(size_t column) const;
  //return a reference to the datum; indexing starts at 1, goes to num_row (inclusive)
  const Tmplt&          operator()(size_t row, size_t column) const;
  Tmplt&                operator()(size_t row, size_t column);

  MMatrix<Tmplt>& operator= (const MMatrix<Tmplt>& mm);

  //TODO - implement iterator
  //class iterator
  //{
  //}


  friend MMatrix<m_complex>& operator *=(MMatrix<m_complex>& m,  m_complex          c);
  friend MMatrix<double>&    operator *=(MMatrix<double>&    m,  double             d);
  friend MMatrix<m_complex>& operator *=(MMatrix<m_complex>& m1, MMatrix<m_complex> m2);
  friend MMatrix<double>&    operator *=(MMatrix<double>&    m1, MMatrix<double>    m2);
  friend MVector<m_complex>  operator * (MMatrix<m_complex>  m,  MVector<m_complex> v);
  friend MVector<double>     operator * (MMatrix<double>     m,  MVector<double>    v);
  friend MMatrix<m_complex>& operator +=(MMatrix<m_complex>& m1, const MMatrix<m_complex>& m2);
  friend MMatrix<double>&    operator +=(MMatrix<double>&    m1, const MMatrix<double>&    m2);
  template <class Tmplt2> friend MMatrix<Tmplt2>      operator + (MMatrix<Tmplt2>      m1, const MMatrix<Tmplt2>      m2);

  friend class MMatrix<double>; //To do the eigenvector problem, MMatrix<double> needs to see MMatrix<complex>'s _matrix


private:
  //build the matrix with size i,j, data not initialised
  void build_matrix(size_t i, size_t j);  
  //build the matrix with size i,j, data initialised to data in temp
  void build_matrix(size_t i, size_t j, Tmplt* temp);  
  //delete the matrix and set it to null  
  void delete_matrix();  
  //return the matrix overloaded to the correct type or throw
  //if _matrix == nullptr
  static gsl_matrix*         get_matrix(const MMatrix<double>&    m);
  static gsl_matrix_complex* get_matrix(const MMatrix<m_complex>& m);
  //gsl object
  void*              _matrix;
};

//Operators do what you would expect - i.e. normal mathematical functions
MMatrix<m_complex>& operator *=(MMatrix<m_complex>& m,  m_complex          c);
MMatrix<double>&    operator *=(MMatrix<double>&    m,  double             d);
MMatrix<m_complex>& operator *=(MMatrix<m_complex>& m1, MMatrix<m_complex> m2); //M1 *= M2 returns M1 = M1*M2
MMatrix<double>&    operator *=(MMatrix<double>&    m1, MMatrix<double>    m2); //M1 *= M2 returns M1 = M1*M2
MVector<m_complex>  operator * (MMatrix<m_complex>  m,  MVector<m_complex> v);
MVector<double>     operator * (MMatrix<double>     m,  MVector<double>    v);
MMatrix<m_complex>& operator +=(MMatrix<m_complex>& m1, const MMatrix<m_complex>& m2);
MMatrix<double>&    operator +=(MMatrix<double>&    m1, const MMatrix<double>&    m2);

template <class Tmplt> MMatrix<Tmplt>& operator -=(MMatrix<double>&    m1, MMatrix<double>    m2);

template <class Tmplt> MMatrix<Tmplt>  operator*(MMatrix<Tmplt>, MMatrix<Tmplt>);
template <class Tmplt> MMatrix<Tmplt>  operator*(MMatrix<Tmplt>, MVector<Tmplt>);
template <class Tmplt> MMatrix<Tmplt>  operator*(MMatrix<Tmplt>, Tmplt);
template <class Tmplt> MMatrix<Tmplt>  operator/(MMatrix<Tmplt>, Tmplt);

template <class Tmplt> MMatrix<Tmplt>  operator+(MMatrix<Tmplt>, MMatrix<Tmplt>);
template <class Tmplt> MMatrix<Tmplt>  operator-(MMatrix<Tmplt>, MMatrix<Tmplt>);
template <class Tmplt> MMatrix<Tmplt>  operator-(MMatrix<Tmplt>); //unary minus returns matrix*(-1)

template <class Tmplt> bool operator==(const MMatrix<Tmplt>& c1, const MMatrix<Tmplt>& c2);
template <class Tmplt> bool operator!=(const MMatrix<Tmplt>& c1, const MMatrix<Tmplt>& c2);

//format is 
// num_row  num_col
// m11 m12 m13 ...
// m21 m22 m23 ...
// ...
template <class Tmplt> std::ostream& operator<<(std::ostream& out, MMatrix<Tmplt>  mat);
template <class Tmplt> std::istream& operator>>(std::istream& in,  MMatrix<Tmplt>& mat);

//return matrix of doubles filled with either real or imaginary part of m
MMatrix<double>    re(MMatrix<m_complex> m);
MMatrix<double>    im(MMatrix<m_complex> m);
//return matrix of m_complex filled with real and imaginary parts
MMatrix<m_complex> complex(MMatrix<double> real);
MMatrix<m_complex> complex(MMatrix<double> real, MMatrix<double> imaginary);

//////////////////////////// MMatrix declaration end ///////////////



//////////////////////////// MMatrix Inlined Functions //////////////

template <>
inline size_t MMatrix<double>::num_row()     const
{ 
  if(_matrix) return ((gsl_matrix*)_matrix)->size1; 
  return 0;
}
template <>
inline size_t MMatrix<m_complex>::num_row()     const
{ 
  if(_matrix) return ((gsl_matrix_complex*)_matrix)->size1; 
  return 0;
}
template <>
inline size_t MMatrix<double>::num_col()     const
{ 
  if(_matrix) return ((gsl_matrix*)_matrix)->size2; 
  return 0;
}

template <>
inline size_t MMatrix<m_complex>::num_col()     const
{
  if(_matrix) return ((gsl_matrix_complex*)_matrix)->size2; 
  return 0;
}

MMatrix<m_complex> inline & operator *=(MMatrix<m_complex>& m,  m_complex c)
{gsl_matrix_complex_scale((gsl_matrix_complex*)m._matrix,  c); return m;}

MMatrix<double>    inline & operator *=(MMatrix<double>& m,  double d)
{gsl_matrix_scale((gsl_matrix*)m._matrix,  d); return m;}

MVector<m_complex> inline     operator * (MMatrix<m_complex> m,  MVector<m_complex>    v)
{
  MVector<m_complex> v0(m.num_row()); 
  gsl_blas_zgemv(CblasNoTrans, m_complex_build(1.), (gsl_matrix_complex*)m._matrix, (gsl_vector_complex*)v._vector, m_complex_build(0.), (gsl_vector_complex*)v0._vector); 
  return v0;
}

MVector<double> inline operator * (MMatrix<double>     m,  MVector<double>    v)
{
  MVector<double> v0(m.num_row());
  gsl_blas_dgemv(CblasNoTrans, 1., (gsl_matrix*)m._matrix, (gsl_vector*)v._vector, 0., (gsl_vector*)v0._vector);
  return v0;
}

MMatrix<m_complex> inline & operator +=(MMatrix<m_complex>& m1, const MMatrix<m_complex>& m2)
{gsl_matrix_complex_add((gsl_matrix_complex*)m1._matrix, (gsl_matrix_complex*)m2._matrix); return m1;}

MMatrix<double> inline &    operator +=(MMatrix<double>& m1, const MMatrix<double>&    m2)
{gsl_matrix_add((gsl_matrix*)m1._matrix, (gsl_matrix*)m2._matrix); return m1;}

template <class Tmplt> MMatrix<Tmplt> inline operator - (MMatrix<Tmplt>    m1)
{for(size_t i=1; i<=m1.num_row(); i++) for(size_t j=1; j<=m1.num_col(); j++) m1(i,j) = -1.*m1(i,j); return m1; }

template <class Tmplt> MMatrix<Tmplt> inline operator*(MMatrix<Tmplt> m1, MMatrix<Tmplt> m2) 
{return m1*=m2;}

template <class Tmplt> MMatrix<Tmplt> inline operator*(MMatrix<Tmplt> m, MVector<Tmplt> v)
{return m*=v;}

template <class Tmplt> MMatrix<Tmplt> inline operator*(MMatrix<Tmplt> m, Tmplt t)
{return m*=t;}

template <class Tmplt> MMatrix<Tmplt> inline & operator/=(MMatrix<Tmplt>& m, Tmplt t)
{return m*=1./t;}
template <class Tmplt> MMatrix<Tmplt> inline operator/(MMatrix<Tmplt>  m, Tmplt t)
{return m/=t;}

template <class Tmplt> MMatrix<Tmplt> inline operator+(MMatrix<Tmplt> m1, MMatrix<Tmplt> m2)
{ MMatrix<Tmplt> m3 = m1+=m2; return m3; }

template <class Tmplt> MMatrix<Tmplt> inline & operator-=(MMatrix<Tmplt>& m1, MMatrix<Tmplt> m2)
{return m1 += -m2;}
template <class Tmplt> MMatrix<Tmplt> inline operator-(MMatrix<Tmplt>  m1, MMatrix<Tmplt> m2)
{return m1 -=  m2;}



template <>
const double inline & MMatrix<double>::operator()(size_t i, size_t j) const
{ return *gsl_matrix_ptr((gsl_matrix*)_matrix, i-1, j-1); }
template <>
const m_complex inline & MMatrix<m_complex>::operator()(size_t i, size_t j) const
{ return *gsl_matrix_complex_ptr((gsl_matrix_complex*)_matrix, i-1, j-1); }

template <>
double inline & MMatrix<double>::operator()(size_t i, size_t j)
{ return *gsl_matrix_ptr((gsl_matrix*)_matrix, i-1, j-1); }
template <>
m_complex inline & MMatrix<m_complex>::operator()(size_t i, size_t j)
{ return *gsl_matrix_complex_ptr((gsl_matrix_complex*)_matrix, i-1, j-1); }

template <class Tmplt>
bool operator==(const MMatrix<Tmplt>& lhs, const MMatrix<Tmplt>& rhs)
{
  if( lhs.num_row() != rhs.num_row() || lhs.num_col() != rhs.num_col() ) return false;
  for(size_t i=1; i<=lhs.num_row(); i++)
    for(size_t j=1; j<=lhs.num_col(); j++)
      if(lhs(i,j) != rhs(i,j)) return false;
  return true;
}

template <class Tmplt>
bool inline operator!=(const MMatrix<Tmplt>& lhs, const MMatrix<Tmplt>& rhs) {return ! (lhs == rhs);}

//iostream
template <class Tmplt> std::ostream& operator<<(std::ostream& out, MMatrix<Tmplt>);
//template <class Tmplt> std::istream& operator>>(std::istream& in,  MMatrix<Tmplt>);

template <class Tmplt>
gsl_matrix         inline* MMatrix<Tmplt>::get_matrix(const MMatrix<double>&    m)
{
  if(m._matrix == nullptr)
    throw(GeneralClassicException("MMatrix::get_matrix", "Attempt to access uninitialised matrix"));
  return (gsl_matrix*)m._matrix;
}

template <class Tmplt>
gsl_matrix_complex inline* MMatrix<Tmplt>::get_matrix(const MMatrix<m_complex>& m)
{
  if(m._matrix == nullptr)
    throw(GeneralClassicException("MMatrix::get_matrix", "Attempt to access uninitialised matrix"));
  return (gsl_matrix_complex*)m._matrix;

}

////////////////////////// MMatrix inlined functions end //////////////////////
}

#endif
