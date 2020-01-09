// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/


///////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//    IpplTypeComputations.h
//
// CREATED
//    July 11, 1997
//
// DESCRIPTION
//    PETE: Portable Expression Template Engine.
//
//    This header file contains IPPL-specific type computations.
//
///////////////////////////////////////////////////////////////////////////

#ifndef IPPL_TYPE_COMPUTATIONS_H
#define IPPL_TYPE_COMPUTATIONS_H


// include files
#include "PETE/TypeComputations.h"
#include "AppTypes/dcomplex.h"


// forward declarations
template<class T, unsigned D> class Vektor;
template<class T, unsigned D> class Tenzor;
template<class T, unsigned D> class AntiSymTenzor;
template<class T, unsigned D> class SymTenzor;
class RNGBitReverse;
template <class T> class RNGLattice;
class RNGSimple;
class RNGRand;
class RNGXDiv;
class RNGXCI;


// definition of global sign function
template<class T>
inline int sign(T a) { return ((a > 0) ? 1 : (a == 0 ? 0 : -1)); }


///////////////////////////////////////////////////////////////////////////
//
// PETE_Type2Index FOR USER TYPES
//
///////////////////////////////////////////////////////////////////////////

// Complex numbers.

template<> struct PETE_Type2Index<dcomplex> {
  enum { val = 8 };
};

#ifdef IPPL_USE_SINGLE_PRECISION
 #define CTYPE float
#else
 #define CTYPE double
#endif

// Return types for scalar ops with RNGs.

#define _SCALAR_RNG_OP_RETURNS_(GEN,SCA,OP)                             \
template <>                                                             \
struct PETEBinaryReturn<GEN,SCA,OP> {                                   \
  typedef PETEBinaryReturn<CTYPE,SCA,OP>::type type;                   \
};                                                                      \
template <>                                                             \
struct PETEBinaryReturn<SCA,GEN,OP> {                                   \
  typedef PETEBinaryReturn<SCA,CTYPE,OP>::type type;                   \
};

#define _SCALAR_RNG_RETURNS_(GEN,SCA)                                   \
_SCALAR_RNG_OP_RETURNS_(GEN,SCA,OpAdd)                                  \
_SCALAR_RNG_OP_RETURNS_(GEN,SCA,OpSubtract)                             \
_SCALAR_RNG_OP_RETURNS_(GEN,SCA,OpMultipply)                             \
_SCALAR_RNG_OP_RETURNS_(GEN,SCA,OpDivide)

#define _PETE_RNG_RETURNS_(GEN)                                         \
                                                                        \
template <> struct PETE_Type2Index< GEN > {                             \
  enum { val = PETE_Type2Index<CTYPE>::val };                          \
};                                                                      \
                                                                        \
_SCALAR_RNG_RETURNS_(GEN,short)                                         \
_SCALAR_RNG_RETURNS_(GEN,int)                                           \
_SCALAR_RNG_RETURNS_(GEN,long)                                          \
_SCALAR_RNG_RETURNS_(GEN,float)                                         \
_SCALAR_RNG_RETURNS_(GEN,double)                                        \
_SCALAR_RNG_RETURNS_(GEN,dcomplex)

_PETE_RNG_RETURNS_(RNGBitReverse)
_PETE_RNG_RETURNS_(RNGLattice<float>)
_PETE_RNG_RETURNS_(RNGLattice<double>)
_PETE_RNG_RETURNS_(RNGSimple)
_PETE_RNG_RETURNS_(RNGRand)
_PETE_RNG_RETURNS_(RNGXDiv)
_PETE_RNG_RETURNS_(RNGXCI)

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

// Life is way easier with this feature.

template<class T, unsigned Dim>
struct PETE_Type2Index< Vektor<T, Dim> > {
  enum { val = 20 + 10 * Dim + PETE_Type2Index<T>::val };
};

template<class T, unsigned Dim>
struct PETE_Type2Index< SymTenzor<T, Dim> > {
  enum { val = 120 + 10 * Dim + PETE_Type2Index<T>::val };
};

template<class T, unsigned Dim>
struct PETE_Type2Index< Tenzor<T, Dim> > {
  enum { val = 220 + 10 * Dim + PETE_Type2Index<T>::val };
};

template<class T, unsigned Dim>
struct PETE_Type2Index< AntiSymTenzor<T, Dim> > {
  enum { val = 320 + 10 * Dim + PETE_Type2Index<T>::val };
};

#else

// Without partial specialization, we must give all the cases. Yuuuck!

// Vectors of size 1.

template<> struct PETE_Type2Index< Vektor<short, 1U> > {
  enum { val = 30 };
};

template<> struct PETE_Type2Index< Vektor<int, 1U> > {
  enum { val = 31 };
};

template<> struct PETE_Type2Index< Vektor<long, 1U> > {
  enum { val = 32 };
};

template<> struct PETE_Type2Index< Vektor<float, 1U> > {
  enum { val = 33 };
};

template<> struct PETE_Type2Index< Vektor<double, 1U> > {
  enum { val = 34 };
};

template<> struct PETE_Type2Index< Vektor<dcomplex, 1U> > {
  enum { val = 35 };
};

// Vectors of size 2.

template<> struct PETE_Type2Index< Vektor<short, 2U> > {
  enum { val = 40 };
};

template<> struct PETE_Type2Index< Vektor<int, 2U> > {
  enum { val = 41 };
};

template<> struct PETE_Type2Index< Vektor<long, 2U> > {
  enum { val = 42 };
};

template<> struct PETE_Type2Index< Vektor<float, 2U> > {
  enum { val = 43 };
};

template<> struct PETE_Type2Index< Vektor<double, 2U> > {
  enum { val = 44 };
};

template<> struct PETE_Type2Index< Vektor<dcomplex, 2U> > {
  enum { val = 45 };
};

// Vectors of size 3.

template<> struct PETE_Type2Index< Vektor<short, 3U> > {
  enum { val = 50 };
};

template<> struct PETE_Type2Index< Vektor<int, 3U> > {
  enum { val = 51 };
};

template<> struct PETE_Type2Index< Vektor<long, 3U> > {
  enum { val = 52 };
};

template<> struct PETE_Type2Index< Vektor<float, 3U> > {
  enum { val = 53 };
};

template<> struct PETE_Type2Index< Vektor<double, 3U> > {
  enum { val = 54 };
};

template<> struct PETE_Type2Index< Vektor<dcomplex, 3U> > {
  enum { val = 55 };
};

// SymTensors of size 1.

template<> struct PETE_Type2Index< SymTenzor<short, 1U> > {
  enum { val = 130 };
};

template<> struct PETE_Type2Index< SymTenzor<int, 1U> > {
  enum { val = 131 };
};

template<> struct PETE_Type2Index< SymTenzor<long, 1U> > {
  enum { val = 132 };
};

template<> struct PETE_Type2Index< SymTenzor<float, 1U> > {
  enum { val = 133 };
};

template<> struct PETE_Type2Index< SymTenzor<double, 1U> > {
  enum { val = 134 };
};

template<> struct PETE_Type2Index< SymTenzor<dcomplex, 1U> > {
  enum { val = 135 };
};

// SymTensors of size 2.

template<> struct PETE_Type2Index< SymTenzor<short, 2U> > {
  enum { val = 140 };
};

template<> struct PETE_Type2Index< SymTenzor<int, 2U> > {
  enum { val = 141 };
};

template<> struct PETE_Type2Index< SymTenzor<long, 2U> > {
  enum { val = 142 };
};

template<> struct PETE_Type2Index< SymTenzor<float, 2U> > {
  enum { val = 143 };
};

template<> struct PETE_Type2Index< SymTenzor<double, 2U> > {
  enum { val = 144 };
};

template<> struct PETE_Type2Index< SymTenzor<dcomplex, 2U> > {
  enum { val = 145 };
};

// SymTensors of size 3.

template<> struct PETE_Type2Index< SymTenzor<short, 3U> > {
  enum { val = 150 };
};

template<> struct PETE_Type2Index< SymTenzor<int, 3U> > {
  enum { val = 151 };
};

template<> struct PETE_Type2Index< SymTenzor<long, 3U> > {
  enum { val = 152 };
};

template<> struct PETE_Type2Index< SymTenzor<float, 3U> > {
  enum { val = 153 };
};

template<> struct PETE_Type2Index< SymTenzor<double, 3U> > {
  enum { val = 154 };
};

template<> struct PETE_Type2Index< SymTenzor<dcomplex, 3U> > {
  enum { val = 155 };
};

// Tensors of size 1.

template<> struct PETE_Type2Index< Tenzor<short, 1U> > {
  enum { val = 230 };
};

template<> struct PETE_Type2Index< Tenzor<int, 1U> > {
  enum { val = 231 };
};

template<> struct PETE_Type2Index< Tenzor<long, 1U> > {
  enum { val = 232 };
};

template<> struct PETE_Type2Index< Tenzor<float, 1U> > {
  enum { val = 233 };
};

template<> struct PETE_Type2Index< Tenzor<double, 1U> > {
  enum { val = 234 };
};

template<> struct PETE_Type2Index< Tenzor<dcomplex, 1U> > {
  enum { val = 235 };
};

// Tensors of size 2.

template<> struct PETE_Type2Index< Tenzor<short, 2U> > {
  enum { val = 240 };
};

template<> struct PETE_Type2Index< Tenzor<int, 2U> > {
  enum { val = 241 };
};

template<> struct PETE_Type2Index< Tenzor<long, 2U> > {
  enum { val = 242 };
};

template<> struct PETE_Type2Index< Tenzor<float, 2U> > {
  enum { val = 243 };
};

template<> struct PETE_Type2Index< Tenzor<double, 2U> > {
  enum { val = 244 };
};

template<> struct PETE_Type2Index< Tenzor<dcomplex, 2U> > {
  enum { val = 245 };
};

// Tensors of size 3.

template<> struct PETE_Type2Index< Tenzor<short, 3U> > {
  enum { val = 250 };
};

template<> struct PETE_Type2Index< Tenzor<int, 3U> > {
  enum { val = 251 };
};

template<> struct PETE_Type2Index< Tenzor<long, 3U> > {
  enum { val = 252 };
};

template<> struct PETE_Type2Index< Tenzor<float, 3U> > {
  enum { val = 253 };
};

template<> struct PETE_Type2Index< Tenzor<double, 3U> > {
  enum { val = 254 };
};

template<> struct PETE_Type2Index< Tenzor<dcomplex, 3U> > {
  enum { val = 255 };
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// PETE_Index2Type FOR USER TYPES
//
///////////////////////////////////////////////////////////////////////////

#if !defined(IPPL_USE_PARTIAL_SPECIALIZATION)

// Complex numbers.

template<> struct PETE_Index2Type<8> {
  typedef dcomplex type;
};

// Vectors of length 1.

template<> struct PETE_Index2Type<30> {
  typedef Vektor<short, 1U> type;
};

template<> struct PETE_Index2Type<31> {
  typedef Vektor<int, 1U> type;
};

template<> struct PETE_Index2Type<32> {
  typedef Vektor<long, 1U> type;
};

template<> struct PETE_Index2Type<33> {
  typedef Vektor<float, 1U> type;
};

template<> struct PETE_Index2Type<34> {
  typedef Vektor<double, 1U> type;
};

template<> struct PETE_Index2Type<35> {
  typedef Vektor<dcomplex, 1U> type;
};

// Vectors of length 2.

template<> struct PETE_Index2Type<40> {
  typedef Vektor<short, 2U> type;
};

template<> struct PETE_Index2Type<41> {
  typedef Vektor<int, 2U> type;
};

template<> struct PETE_Index2Type<42> {
  typedef Vektor<long, 2U> type;
};

template<> struct PETE_Index2Type<43> {
  typedef Vektor<float, 2U> type;
};

template<> struct PETE_Index2Type<44> {
  typedef Vektor<double, 2U> type;
};

template<> struct PETE_Index2Type<45> {
  typedef Vektor<dcomplex, 2U> type;
};

// Vectors of length 3.

template<> struct PETE_Index2Type<50> {
  typedef Vektor<short, 3U> type;
};

template<> struct PETE_Index2Type<51> {
  typedef Vektor<int, 3U> type;
};

template<> struct PETE_Index2Type<52> {
  typedef Vektor<long, 3U> type;
};

template<> struct PETE_Index2Type<53> {
  typedef Vektor<float, 3U> type;
};

template<> struct PETE_Index2Type<54> {
  typedef Vektor<double, 3U> type;
};

template<> struct PETE_Index2Type<55> {
  typedef Vektor<dcomplex, 3U> type;
};

// SymTensors of length 1.

template<> struct PETE_Index2Type<130> {
  typedef SymTenzor<short, 1U> type;
};

template<> struct PETE_Index2Type<131> {
  typedef SymTenzor<int, 1U> type;
};

template<> struct PETE_Index2Type<132> {
  typedef SymTenzor<long, 1U> type;
};

template<> struct PETE_Index2Type<133> {
  typedef SymTenzor<float, 1U> type;
};

template<> struct PETE_Index2Type<134> {
  typedef SymTenzor<double, 1U> type;
};

template<> struct PETE_Index2Type<135> {
  typedef SymTenzor<dcomplex, 1U> type;
};

// SymTensors of length 2.

template<> struct PETE_Index2Type<140> {
  typedef SymTenzor<short, 2U> type;
};

template<> struct PETE_Index2Type<141> {
  typedef SymTenzor<int, 2U> type;
};

template<> struct PETE_Index2Type<142> {
  typedef SymTenzor<long, 2U> type;
};

template<> struct PETE_Index2Type<143> {
  typedef SymTenzor<float, 2U> type;
};

template<> struct PETE_Index2Type<144> {
  typedef SymTenzor<double, 2U> type;
};

template<> struct PETE_Index2Type<145> {
  typedef SymTenzor<dcomplex, 2U> type;
};

// SymTensors of length 3.

template<> struct PETE_Index2Type<150> {
  typedef SymTenzor<short, 3U> type;
};

template<> struct PETE_Index2Type<151> {
  typedef SymTenzor<int, 3U> type;
};

template<> struct PETE_Index2Type<152> {
  typedef SymTenzor<long, 3U> type;
};

template<> struct PETE_Index2Type<153> {
  typedef SymTenzor<float, 3U> type;
};

template<> struct PETE_Index2Type<154> {
  typedef SymTenzor<double, 3U> type;
};

template<> struct PETE_Index2Type<155> {
  typedef SymTenzor<dcomplex, 3U> type;
};

// Tensors of length 1.

template<> struct PETE_Index2Type<230> {
  typedef Tenzor<short, 1U> type;
};

template<> struct PETE_Index2Type<231> {
  typedef Tenzor<int, 1U> type;
};

template<> struct PETE_Index2Type<232> {
  typedef Tenzor<long, 1U> type;
};

template<> struct PETE_Index2Type<233> {
  typedef Tenzor<float, 1U> type;
};

template<> struct PETE_Index2Type<234> {
  typedef Tenzor<double, 1U> type;
};

template<> struct PETE_Index2Type<235> {
  typedef Tenzor<dcomplex, 1U> type;
};

// Tensors of length 2.

template<> struct PETE_Index2Type<240> {
  typedef Tenzor<short, 2U> type;
};

template<> struct PETE_Index2Type<241> {
  typedef Tenzor<int, 2U> type;
};

template<> struct PETE_Index2Type<242> {
  typedef Tenzor<long, 2U> type;
};

template<> struct PETE_Index2Type<243> {
  typedef Tenzor<float, 2U> type;
};

template<> struct PETE_Index2Type<244> {
  typedef Tenzor<double, 2U> type;
};

template<> struct PETE_Index2Type<245> {
  typedef Tenzor<dcomplex, 2U> type;
};

// Tensors of length 3.

template<> struct PETE_Index2Type<250> {
  typedef Tenzor<short, 3U> type;
};

template<> struct PETE_Index2Type<251> {
  typedef Tenzor<int, 3U> type;
};

template<> struct PETE_Index2Type<252> {
  typedef Tenzor<long, 3U> type;
};

template<> struct PETE_Index2Type<253> {
  typedef Tenzor<float, 3U> type;
};

template<> struct PETE_Index2Type<254> {
  typedef Tenzor<double, 3U> type;
};

template<> struct PETE_Index2Type<255> {
  typedef Tenzor<dcomplex, 3U> type;
};

#endif


///////////////////////////////////////////////////////////////////////////
//
// SPECIAL CASES FOR UNARY FUNCTIONS
//
///////////////////////////////////////////////////////////////////////////

// Abs function: special return for complex numbers.

struct FnAbs {
  enum { tag = PETE_UnaryPassThruTag };
};

template<> struct PETEUnaryReturn<dcomplex, FnAbs> {
  typedef double type;
};

// The conj, norm, arg, real, and imag functions for complex numbers.

struct FnConj {
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnNorm {
  typedef double type;
  enum { tag = PETE_Type2Index<double>::val };
};

template<> struct PETEUnaryReturn<dcomplex, FnNorm> {
  typedef double type;
};

struct FnArg {
  typedef double type;
  enum { tag = PETE_Type2Index<double>::val };
};

template<> struct PETEUnaryReturn<dcomplex, FnArg> {
  typedef double type;
};

struct FnReal {
  typedef double type;
  enum { tag = PETE_Type2Index<double>::val };
};

template<> struct PETEUnaryReturn<dcomplex, FnReal> {
  typedef double type;
};

struct FnImag {
  typedef double type;
  enum { tag = PETE_Type2Index<double>::val };
};

template<> struct PETEUnaryReturn<dcomplex, FnImag> {
  typedef double type;
};

// The sign function.

struct FnSign {
  typedef int type;
  enum { tag = PETE_Type2Index<int>::val };
};

template<class TP>
struct OpParens
{
  enum { tag = PETE_UnaryPassThruTag };
  TP Arg;
  OpParens() { Arg = TP(); }
  OpParens(const TP& a) : Arg(a) {}
};

// Tensor functions: trace, det (determinant),  and transpose

struct FnTrace {
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnDet {
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnTranspose {
  enum { tag = PETE_UnaryPassThruTag };
};

struct FnCofactors {
  enum { tag = PETE_UnaryPassThruTag };
};

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

// Life is pretty simple if we have partial specialization.

template<class T, unsigned Dim>
struct PETEUnaryReturn<Tenzor<T,Dim>, FnTrace> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<SymTenzor<T,Dim>, FnTrace> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<AntiSymTenzor<T,Dim>, FnTrace> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<Tenzor<T,Dim>, FnDet> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<SymTenzor<T,Dim>, FnDet> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<AntiSymTenzor<T,Dim>, FnDet> {
  typedef T type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<Tenzor<T,Dim>, FnTranspose> {
  typedef Tenzor<T,Dim> type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<SymTenzor<T,Dim>, FnTranspose> {
  typedef SymTenzor<T,Dim> type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<AntiSymTenzor<T,Dim>, FnTranspose> {
  typedef AntiSymTenzor<T,Dim> type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<Tenzor<T,Dim>, FnCofactors> {
  typedef Tenzor<T,Dim> type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<SymTenzor<T,Dim>, FnCofactors> {
  typedef SymTenzor<T,Dim> type;
};

template<class T, unsigned Dim>
struct PETEUnaryReturn<AntiSymTenzor<T,Dim>, FnCofactors> {
  typedef AntiSymTenzor<T,Dim> type;
};

#else

// Life is pitiful and barely worth living without partial specialization.
// Sigh, I guess we'll limp along...

#define _UNARY_TENSOR_RETURNS_(T, D)
template<> struct PETEUnaryReturn<Tenzor<T,D>, FnTrace>
{ typedef T type; };
template<> struct PETEUnaryReturn<SymTenzor<T,D>, FnTrace>
{ typedef T type; };
template<> struct PETEUnaryReturn<AntiSymTenzor<T,D>, FnTrace>
{ typedef T type; };
template<> struct PETEUnaryReturn<Tenzor<T,D>, FnDet>
{ typedef T type; };
template<> struct PETEUnaryReturn<SymTenzor<T,D>, FnDet>
{ typedef T type; };
template<> struct PETEUnaryReturn<AntiSymTenzor<T,D>, FnDet>
{ typedef T type; };
template<> struct PETEUnaryReturn<Tenzor<T,D>, FnTranspose>
{ typedef Tenzor<T,D> type; };
template<> struct PETEUnaryReturn<SymTenzor<T,D>, FnTranspose>
{ typedef SymTenzor<T,D> type; };
template<> struct PETEUnaryReturn<AntiSymTenzor<T,D>, FnTranspose>
{ typedef AntiSymTenzor<T,D> type; };
template<> struct PETEUnaryReturn<Tenzor<T,D>, FnCofactors>
{ typedef Tenzor<T,D> type; };
template<> struct PETEUnaryReturn<SymTenzor<T,D>, FnCofactors>
{ typedef Tenzor<T,D> type; };
template<> struct PETEUnaryReturn<AntiSymTenzor<T,D>, FnCofactors>
{ typedef SymTenzor<T,D> type; };

_UNARY_TENSOR_RETURNS_(short,1U)
_UNARY_TENSOR_RETURNS_(int,1U)
_UNARY_TENSOR_RETURNS_(long,1U)
_UNARY_TENSOR_RETURNS_(float,1U)
_UNARY_TENSOR_RETURNS_(double,1U)
_UNARY_TENSOR_RETURNS_(dcomplex,1U)

_UNARY_TENSOR_RETURNS_(short,2U)
_UNARY_TENSOR_RETURNS_(int,2U)
_UNARY_TENSOR_RETURNS_(long,2U)
_UNARY_TENSOR_RETURNS_(float,2U)
_UNARY_TENSOR_RETURNS_(double,2U)
_UNARY_TENSOR_RETURNS_(dcomplex,2U)

_UNARY_TENSOR_RETURNS_(short,3U)
_UNARY_TENSOR_RETURNS_(int,3U)
_UNARY_TENSOR_RETURNS_(long,3U)
_UNARY_TENSOR_RETURNS_(float,3U)
_UNARY_TENSOR_RETURNS_(double,3U)
_UNARY_TENSOR_RETURNS_(dcomplex,3U)

#undef _UNARY_TENSOR_RETURNS_

#endif


///////////////////////////////////////////////////////////////////////////
//
// SPECIAL CASES FOR BINARY FUNCTIONS
//
///////////////////////////////////////////////////////////////////////////

// Min and Max functions.

struct FnMin {
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnMax {
  enum { tag = PETE_BinaryPromoteTag };
};

// Dot, dot-dot, and outerProduct functions.

struct FnDot {
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnDotDot {
  enum { tag = PETE_BinaryPromoteTag };
};

struct FnOuterProduct {
  enum { tag = PETE_BinaryPromoteTag };
};

// Cross-product:

struct FnCross {
  enum { tag = PETE_BinaryPromoteTag };
};

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

// Life is pretty simple if we have partial specialization.


// Involving Vektors:

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>, FnCross> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>, FnOuterProduct> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>, FnDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

// Involving Tenzors, but no combination with SymTenzors or AntiSymTenzors:

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,Tenzor<T2,Dim>,FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,Tenzor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,Vektor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,Tenzor<T2,Dim>,FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

// Involving SymTenzors, possibly combined with Tenzors:

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> 
    type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,SymTenzor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Vektor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpMultipply> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

// Involving AntiSymTenzors, possibly combined with Tenzors or SymTenzors:

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> 
    type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot>
{
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Vektor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Vektor<T2,Dim>, FnDot> {
  typedef Vektor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpMultipply> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpMultipply> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,OpAdd> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,OpSubtract> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,FnDot> {
  typedef Tenzor<typename PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type;
};

template<class T1, class T2, unsigned Dim>
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,FnDotDot> {
  typedef typename PETEBinaryReturn<T1,T2,OpMultipply>::type type;
};

// Need to specify scalar operations directly.

#define _SCALAR_VST_RETURNS_(Sca)                                           \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<Vektor<T1,Dim>,Sca,OpMultipply> {                    \
  typedef Vektor<typename PETEBinaryReturn<T1,Sca,OpMultipply>::type,Dim>    \
    type;                                                                   \
};                                                                          \
template<class T2, unsigned Dim>                                            \
struct PETEBinaryReturn<Sca,Vektor<T2,Dim>,OpMultipply> {                    \
  typedef Vektor<typename PETEBinaryReturn<Sca,T2,OpMultipply>::type,Dim>    \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<Vektor<T1,Dim>,Sca,OpDivide> {                      \
  typedef Vektor<typename PETEBinaryReturn<T1,Sca,OpDivide>::type,Dim>      \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<Tenzor<T1,Dim>,Sca,OpMultipply> {                    \
  typedef Tenzor<typename PETEBinaryReturn<T1,Sca,OpMultipply>::type,Dim>    \
    type;                                                                   \
};                                                                          \
template<class T2, unsigned Dim>                                            \
struct PETEBinaryReturn<Sca,Tenzor<T2,Dim>,OpMultipply> {                    \
  typedef Tenzor<typename PETEBinaryReturn<Sca,T2,OpMultipply>::type,Dim>    \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<Tenzor<T1,Dim>,Sca,OpDivide> {                      \
  typedef Tenzor<typename PETEBinaryReturn<T1,Sca,OpDivide>::type,Dim>      \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Sca,OpMultipply> {                 \
  typedef SymTenzor<typename PETEBinaryReturn<T1,Sca,OpMultipply>::type,Dim> \
    type;                                                                   \
};                                                                          \
template<class T2, unsigned Dim>                                            \
struct PETEBinaryReturn<Sca,SymTenzor<T2,Dim>,OpMultipply> {                 \
  typedef SymTenzor<typename PETEBinaryReturn<Sca,T2,OpMultipply>::type,Dim> \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<SymTenzor<T1,Dim>,Sca,OpDivide> {                   \
  typedef SymTenzor<typename PETEBinaryReturn<T1,Sca,OpDivide>::type,Dim>   \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Sca,OpMultipply> {             \
  typedef                                                                   \
  AntiSymTenzor<typename PETEBinaryReturn<T1,Sca,OpMultipply>::type,Dim>     \
    type;                                                                   \
};                                                                          \
template<class T2, unsigned Dim>                                            \
struct PETEBinaryReturn<Sca,AntiSymTenzor<T2,Dim>,OpMultipply> {             \
  typedef                                                                   \
  AntiSymTenzor<typename PETEBinaryReturn<Sca,T2,OpMultipply>::type,Dim>     \
    type;                                                                   \
};                                                                          \
template<class T1, unsigned Dim>                                            \
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Sca,OpDivide> {               \
  typedef                                                                   \
  AntiSymTenzor<typename PETEBinaryReturn<T1,Sca,OpDivide>::type,Dim>       \
    type;                                                                   \
};

_SCALAR_VST_RETURNS_(short)
_SCALAR_VST_RETURNS_(int)
_SCALAR_VST_RETURNS_(long)
_SCALAR_VST_RETURNS_(float)
_SCALAR_VST_RETURNS_(double)
_SCALAR_VST_RETURNS_(dcomplex)

#undef _SCALAR_VST_RETURNS_

#else

// Without partial specialization, we must list all the possibilities.
// Life is too short to type all this stuff out, so we use macros.

#define _VEKTOR_RETURNS_(T1,T2,Dim)				            \
template<> struct                                                           \
  PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>, FnOuterProduct>           \
  { typedef Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim>            \
    type; };                                                                \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>,FnDot>	    \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };	            \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,T2,OpMultipply>	    \
  { typedef Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim>            \
    type; };                                                                \
template<> struct PETEBinaryReturn<T2,Vektor<T1,Dim>,OpMultipply>	    \
  { typedef Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim>            \
    type; };                                                                \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,T2,OpDivide>		    \
  { typedef Vektor<PETEBinaryReturn<T1,T2,OpDivide>::type,Dim>              \
    type; };                                                                \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,Vektor<T2,Dim>,FnCross>   \
  { typedef Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };

_VEKTOR_RETURNS_(short,short,1U)
_VEKTOR_RETURNS_(int,int,1U)
_VEKTOR_RETURNS_(long,long,1U)
_VEKTOR_RETURNS_(float,float,1U)
_VEKTOR_RETURNS_(double,double,1U)
_VEKTOR_RETURNS_(dcomplex,dcomplex,1U)
_VEKTOR_RETURNS_(int,short,1U)
_VEKTOR_RETURNS_(long,short,1U)
_VEKTOR_RETURNS_(long,int,1U)
_VEKTOR_RETURNS_(float,short,1U)
_VEKTOR_RETURNS_(float,int,1U)
_VEKTOR_RETURNS_(float,long,1U)
_VEKTOR_RETURNS_(double,short,1U)
_VEKTOR_RETURNS_(double,int,1U)
_VEKTOR_RETURNS_(double,long,1U)
_VEKTOR_RETURNS_(double,float,1U)
_VEKTOR_RETURNS_(dcomplex,short,1U)
_VEKTOR_RETURNS_(dcomplex,int,1U)
_VEKTOR_RETURNS_(dcomplex,long,1U)
_VEKTOR_RETURNS_(dcomplex,float,1U)
_VEKTOR_RETURNS_(dcomplex,double,1U)

_VEKTOR_RETURNS_(short,short,2U)
_VEKTOR_RETURNS_(int,int,2U)
_VEKTOR_RETURNS_(long,long,2U)
_VEKTOR_RETURNS_(float,float,2U)
_VEKTOR_RETURNS_(double,double,2U)
_VEKTOR_RETURNS_(dcomplex,dcomplex,2U)
_VEKTOR_RETURNS_(int,short,2U)
_VEKTOR_RETURNS_(long,short,2U)
_VEKTOR_RETURNS_(long,int,2U)
_VEKTOR_RETURNS_(float,short,2U)
_VEKTOR_RETURNS_(float,int,2U)
_VEKTOR_RETURNS_(float,long,2U)
_VEKTOR_RETURNS_(double,short,2U)
_VEKTOR_RETURNS_(double,int,2U)
_VEKTOR_RETURNS_(double,long,2U)
_VEKTOR_RETURNS_(double,float,2U)
_VEKTOR_RETURNS_(dcomplex,short,2U)
_VEKTOR_RETURNS_(dcomplex,int,2U)
_VEKTOR_RETURNS_(dcomplex,long,2U)
_VEKTOR_RETURNS_(dcomplex,float,2U)
_VEKTOR_RETURNS_(dcomplex,double,2U)

_VEKTOR_RETURNS_(short,short,3U)
_VEKTOR_RETURNS_(int,int,3U)
_VEKTOR_RETURNS_(long,long,3U)
_VEKTOR_RETURNS_(float,float,3U)
_VEKTOR_RETURNS_(double,double,3U)
_VEKTOR_RETURNS_(dcomplex,dcomplex,3U)
_VEKTOR_RETURNS_(int,short,3U)
_VEKTOR_RETURNS_(long,short,3U)
_VEKTOR_RETURNS_(long,int,3U)
_VEKTOR_RETURNS_(float,short,3U)
_VEKTOR_RETURNS_(float,int,3U)
_VEKTOR_RETURNS_(float,long,3U)
_VEKTOR_RETURNS_(double,short,3U)
_VEKTOR_RETURNS_(double,int,3U)
_VEKTOR_RETURNS_(double,long,3U)
_VEKTOR_RETURNS_(double,float,3U)
_VEKTOR_RETURNS_(dcomplex,short,3U)
_VEKTOR_RETURNS_(dcomplex,int,3U)
_VEKTOR_RETURNS_(dcomplex,long,3U)
_VEKTOR_RETURNS_(dcomplex,float,3U)
_VEKTOR_RETURNS_(dcomplex,double,3U)

#undef _VEKTOR_RETURNS_

#define _TENZOR_RETURNS_(T1,T2,Dim)			                    \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,Tenzor<T2,Dim>, FnDot>    \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,Tenzor<T2,Dim>, FnDotDot> \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };     	    \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,Tenzor<T2,Dim>, FnDot>    \
  { typedef                                                                 \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,Vektor<T2,Dim>, FnDot>    \
  { typedef                                                                 \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<T1,Tenzor<T2,Dim>,OpMultipply>	    \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,T2,OpMultipply>	    \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,T2,OpDivide>		    \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpDivide>::type,Dim> type; };

_TENZOR_RETURNS_(short,short,1U)
_TENZOR_RETURNS_(int,int,1U)
_TENZOR_RETURNS_(long,long,1U)
_TENZOR_RETURNS_(float,float,1U)
_TENZOR_RETURNS_(double,double,1U)
_TENZOR_RETURNS_(dcomplex,dcomplex,1U)
_TENZOR_RETURNS_(int,short,1U)
_TENZOR_RETURNS_(long,short,1U)
_TENZOR_RETURNS_(long,int,1U)
_TENZOR_RETURNS_(float,short,1U)
_TENZOR_RETURNS_(float,int,1U)
_TENZOR_RETURNS_(float,long,1U)
_TENZOR_RETURNS_(double,short,1U)
_TENZOR_RETURNS_(double,int,1U)
_TENZOR_RETURNS_(double,long,1U)
_TENZOR_RETURNS_(double,float,1U)
_TENZOR_RETURNS_(dcomplex,short,1U)
_TENZOR_RETURNS_(dcomplex,int,1U)
_TENZOR_RETURNS_(dcomplex,long,1U)
_TENZOR_RETURNS_(dcomplex,float,1U)
_TENZOR_RETURNS_(dcomplex,double,1U)

_TENZOR_RETURNS_(short,short,2U)
_TENZOR_RETURNS_(int,int,2U)
_TENZOR_RETURNS_(long,long,2U)
_TENZOR_RETURNS_(float,float,2U)
_TENZOR_RETURNS_(double,double,2U)
_TENZOR_RETURNS_(dcomplex,dcomplex,2U)
_TENZOR_RETURNS_(int,short,2U)
_TENZOR_RETURNS_(long,short,2U)
_TENZOR_RETURNS_(long,int,2U)
_TENZOR_RETURNS_(float,short,2U)
_TENZOR_RETURNS_(float,int,2U)
_TENZOR_RETURNS_(float,long,2U)
_TENZOR_RETURNS_(double,short,2U)
_TENZOR_RETURNS_(double,int,2U)
_TENZOR_RETURNS_(double,long,2U)
_TENZOR_RETURNS_(double,float,2U)
_TENZOR_RETURNS_(dcomplex,short,2U)
_TENZOR_RETURNS_(dcomplex,int,2U)
_TENZOR_RETURNS_(dcomplex,long,2U)
_TENZOR_RETURNS_(dcomplex,float,2U)
_TENZOR_RETURNS_(dcomplex,double,2U)

_TENZOR_RETURNS_(short,short,3U)
_TENZOR_RETURNS_(int,int,3U)
_TENZOR_RETURNS_(long,long,3U)
_TENZOR_RETURNS_(float,float,3U)
_TENZOR_RETURNS_(double,double,3U)
_TENZOR_RETURNS_(dcomplex,dcomplex,3U)
_TENZOR_RETURNS_(int,short,3U)
_TENZOR_RETURNS_(long,short,3U)
_TENZOR_RETURNS_(long,int,3U)
_TENZOR_RETURNS_(float,short,3U)
_TENZOR_RETURNS_(float,int,3U)
_TENZOR_RETURNS_(float,long,3U)
_TENZOR_RETURNS_(double,short,3U)
_TENZOR_RETURNS_(double,int,3U)
_TENZOR_RETURNS_(double,long,3U)
_TENZOR_RETURNS_(double,float,3U)
_TENZOR_RETURNS_(dcomplex,short,3U)
_TENZOR_RETURNS_(dcomplex,int,3U)
_TENZOR_RETURNS_(dcomplex,long,3U)
_TENZOR_RETURNS_(dcomplex,float,3U)
_TENZOR_RETURNS_(dcomplex,double,3U)

#undef _TENZOR_RETURNS_


#define _SYMTENZOR_RETURNS_(T1,T2,Dim)				            \
template<> struct                                                           \
  PETEBinaryReturn<SymTenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDot>    	    \
  { typedef                                                                 \
      Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };         \
template<> struct                                                           \
  PETEBinaryReturn<SymTenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDotDot>	    \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };    	    \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,SymTenzor<T2,Dim>, FnDot> \
  { typedef                                                                 \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<SymTenzor<T1,Dim>,Vektor<T2,Dim>, FnDot> \
  { typedef                                                                 \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpAdd>  \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	            \
template<> struct                                                           \
  PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpSubtract>	            \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };           \
template<> struct                                                           \
  PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>,OpMultipply>	            \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDot> \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct                                                           \
  PETEBinaryReturn<Tenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDotDot>    	    \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };	            \
template<> struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpAdd>  \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	            \
template<> struct                                                           \
  PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpSubtract>	            \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };           \
template<> struct PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>, FnDot> \
  { typedef                                                                 \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct                                                           \
  PETEBinaryReturn<SymTenzor<T1,Dim>,Tenzor<T2,Dim>, FnDotDot>    	    \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };    	    \
template<> struct PETEBinaryReturn<SymTenzor<T1,Dim>,T2,OpMultipply>	    \
  { typedef                                                                 \
      SymTenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };      \
template<> struct PETEBinaryReturn<T1,SymTenzor<T2,Dim>,OpMultipply>	    \
  { typedef                                                                 \
      SymTenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim>               \
      type; };                                                              \
template<> struct PETEBinaryReturn<SymTenzor<T1,Dim>,T2,OpDivide>           \
  { typedef SymTenzor<PETEBinaryReturn<T1,T2,OpDivide>::type,Dim> type; };

_SYMTENZOR_RETURNS_(short,short,1U)
_SYMTENZOR_RETURNS_(int,int,1U)
_SYMTENZOR_RETURNS_(long,long,1U)
_SYMTENZOR_RETURNS_(float,float,1U)
_SYMTENZOR_RETURNS_(double,double,1U)
_SYMTENZOR_RETURNS_(dcomplex,dcomplex,1U)
_SYMTENZOR_RETURNS_(int,short,1U)
_SYMTENZOR_RETURNS_(long,short,1U)
_SYMTENZOR_RETURNS_(long,int,1U)
_SYMTENZOR_RETURNS_(float,short,1U)
_SYMTENZOR_RETURNS_(float,int,1U)
_SYMTENZOR_RETURNS_(float,long,1U)
_SYMTENZOR_RETURNS_(double,short,1U)
_SYMTENZOR_RETURNS_(double,int,1U)
_SYMTENZOR_RETURNS_(double,long,1U)
_SYMTENZOR_RETURNS_(double,float,1U)
_SYMTENZOR_RETURNS_(dcomplex,short,1U)
_SYMTENZOR_RETURNS_(dcomplex,int,1U)
_SYMTENZOR_RETURNS_(dcomplex,long,1U)
_SYMTENZOR_RETURNS_(dcomplex,float,1U)
_SYMTENZOR_RETURNS_(dcomplex,double,1U)

_SYMTENZOR_RETURNS_(short,short,2U)
_SYMTENZOR_RETURNS_(int,int,2U)
_SYMTENZOR_RETURNS_(long,long,2U)
_SYMTENZOR_RETURNS_(float,float,2U)
_SYMTENZOR_RETURNS_(double,double,2U)
_SYMTENZOR_RETURNS_(dcomplex,dcomplex,2U)
_SYMTENZOR_RETURNS_(int,short,2U)
_SYMTENZOR_RETURNS_(long,short,2U)
_SYMTENZOR_RETURNS_(long,int,2U)
_SYMTENZOR_RETURNS_(float,short,2U)
_SYMTENZOR_RETURNS_(float,int,2U)
_SYMTENZOR_RETURNS_(float,long,2U)
_SYMTENZOR_RETURNS_(double,short,2U)
_SYMTENZOR_RETURNS_(double,int,2U)
_SYMTENZOR_RETURNS_(double,long,2U)
_SYMTENZOR_RETURNS_(double,float,2U)
_SYMTENZOR_RETURNS_(dcomplex,short,2U)
_SYMTENZOR_RETURNS_(dcomplex,int,2U)
_SYMTENZOR_RETURNS_(dcomplex,long,2U)
_SYMTENZOR_RETURNS_(dcomplex,float,2U)
_SYMTENZOR_RETURNS_(dcomplex,double,2U)

_SYMTENZOR_RETURNS_(short,short,3U)
_SYMTENZOR_RETURNS_(int,int,3U)
_SYMTENZOR_RETURNS_(long,long,3U)
_SYMTENZOR_RETURNS_(float,float,3U)
_SYMTENZOR_RETURNS_(double,double,3U)
_SYMTENZOR_RETURNS_(dcomplex,dcomplex,3U)
_SYMTENZOR_RETURNS_(int,short,3U)
_SYMTENZOR_RETURNS_(long,short,3U)
_SYMTENZOR_RETURNS_(long,int,3U)
_SYMTENZOR_RETURNS_(float,short,3U)
_SYMTENZOR_RETURNS_(float,int,3U)
_SYMTENZOR_RETURNS_(float,long,3U)
_SYMTENZOR_RETURNS_(double,short,3U)
_SYMTENZOR_RETURNS_(double,int,3U)
_SYMTENZOR_RETURNS_(double,long,3U)
_SYMTENZOR_RETURNS_(double,float,3U)
_SYMTENZOR_RETURNS_(dcomplex,short,3U)
_SYMTENZOR_RETURNS_(dcomplex,int,3U)
_SYMTENZOR_RETURNS_(dcomplex,long,3U)
_SYMTENZOR_RETURNS_(dcomplex,float,3U)
_SYMTENZOR_RETURNS_(dcomplex,double,3U)

#undef _SYMTENZOR_RETURNS_


#define _ANTISYMTENZOR_RETURNS_(T1,T2,Dim)				      \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDot>        \
  { typedef                                                                   \
      Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };           \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot>     \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };    	      \
template<> struct PETEBinaryReturn<Vektor<T1,Dim>,AntiSymTenzor<T2,Dim>,FnDot>\
  { typedef                                                                   \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Vektor<T2,Dim>,FnDot>\
  { typedef                                                                   \
    Vektor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpAdd>\
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	              \
template<> struct                                                             \
  PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpSubtract>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpMultipply>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,FnDot>\
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<Tenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot>    	      \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };	              \
template<> struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpAdd>\
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	              \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,OpSubtract>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };             \
template<> struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>,FnDot>\
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,Tenzor<T2,Dim>, FnDotDot>    	      \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };    	      \
template<>                                                                    \
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpAdd>        \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	              \
template<> struct                                                             \
  PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpSubtract>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,OpMultipply>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<>                                                                    \
struct PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>,FnDot>        \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<SymTenzor<T1,Dim>,AntiSymTenzor<T2,Dim>, FnDotDot>         \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };	              \
template<>                                                                    \
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,OpAdd>        \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpAdd>::type,Dim> type; };	              \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,OpSubtract>	      \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpSubtract>::type,Dim> type; };             \
template<>                                                                    \
struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>,FnDot>        \
  { typedef                                                                   \
    Tenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };             \
template<> struct                                                             \
  PETEBinaryReturn<AntiSymTenzor<T1,Dim>,SymTenzor<T2,Dim>, FnDotDot>         \
  { typedef PETEBinaryReturn<T1,T2,OpMultipply>::type type; };    	      \
template<> struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,T2,OpMultipply>	      \
  { typedef                                                                   \
      AntiSymTenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim> type; };    \
template<> struct PETEBinaryReturn<T1,AntiSymTenzor<T2,Dim>,OpMultipply>	      \
  { typedef                                                                   \
      AntiSymTenzor<PETEBinaryReturn<T1,T2,OpMultipply>::type,Dim>             \
      type; };                                                                \
template<> struct PETEBinaryReturn<AntiSymTenzor<T1,Dim>,T2,OpDivide>         \
  { typedef AntiSymTenzor<PETEBinaryReturn<T1,T2,OpDivide>::type,Dim> type; };

_ANTISYMTENZOR_RETURNS_(short,short,1U)
_ANTISYMTENZOR_RETURNS_(int,int,1U)
_ANTISYMTENZOR_RETURNS_(long,long,1U)
_ANTISYMTENZOR_RETURNS_(float,float,1U)
_ANTISYMTENZOR_RETURNS_(double,double,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,dcomplex,1U)
_ANTISYMTENZOR_RETURNS_(int,short,1U)
_ANTISYMTENZOR_RETURNS_(long,short,1U)
_ANTISYMTENZOR_RETURNS_(long,int,1U)
_ANTISYMTENZOR_RETURNS_(float,short,1U)
_ANTISYMTENZOR_RETURNS_(float,int,1U)
_ANTISYMTENZOR_RETURNS_(float,long,1U)
_ANTISYMTENZOR_RETURNS_(double,short,1U)
_ANTISYMTENZOR_RETURNS_(double,int,1U)
_ANTISYMTENZOR_RETURNS_(double,long,1U)
_ANTISYMTENZOR_RETURNS_(double,float,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,short,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,int,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,long,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,float,1U)
_ANTISYMTENZOR_RETURNS_(dcomplex,double,1U)

_ANTISYMTENZOR_RETURNS_(short,short,2U)
_ANTISYMTENZOR_RETURNS_(int,int,2U)
_ANTISYMTENZOR_RETURNS_(long,long,2U)
_ANTISYMTENZOR_RETURNS_(float,float,2U)
_ANTISYMTENZOR_RETURNS_(double,double,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,dcomplex,2U)
_ANTISYMTENZOR_RETURNS_(int,short,2U)
_ANTISYMTENZOR_RETURNS_(long,short,2U)
_ANTISYMTENZOR_RETURNS_(long,int,2U)
_ANTISYMTENZOR_RETURNS_(float,short,2U)
_ANTISYMTENZOR_RETURNS_(float,int,2U)
_ANTISYMTENZOR_RETURNS_(float,long,2U)
_ANTISYMTENZOR_RETURNS_(double,short,2U)
_ANTISYMTENZOR_RETURNS_(double,int,2U)
_ANTISYMTENZOR_RETURNS_(double,long,2U)
_ANTISYMTENZOR_RETURNS_(double,float,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,short,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,int,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,long,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,float,2U)
_ANTISYMTENZOR_RETURNS_(dcomplex,double,2U)

_ANTISYMTENZOR_RETURNS_(short,short,3U)
_ANTISYMTENZOR_RETURNS_(int,int,3U)
_ANTISYMTENZOR_RETURNS_(long,long,3U)
_ANTISYMTENZOR_RETURNS_(float,float,3U)
_ANTISYMTENZOR_RETURNS_(double,double,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,dcomplex,3U)
_ANTISYMTENZOR_RETURNS_(int,short,3U)
_ANTISYMTENZOR_RETURNS_(long,short,3U)
_ANTISYMTENZOR_RETURNS_(long,int,3U)
_ANTISYMTENZOR_RETURNS_(float,short,3U)
_ANTISYMTENZOR_RETURNS_(float,int,3U)
_ANTISYMTENZOR_RETURNS_(float,long,3U)
_ANTISYMTENZOR_RETURNS_(double,short,3U)
_ANTISYMTENZOR_RETURNS_(double,int,3U)
_ANTISYMTENZOR_RETURNS_(double,long,3U)
_ANTISYMTENZOR_RETURNS_(double,float,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,short,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,int,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,long,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,float,3U)
_ANTISYMTENZOR_RETURNS_(dcomplex,double,3U)

#undef _SYMTENZOR_RETURNS_


#endif


///////////////////////////////////////////////////////////////////////////
//
// ASSIGNMENT OPERATORS: min=, max=, &&=, ||=
//
///////////////////////////////////////////////////////////////////////////

struct OpMinAssign {
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpMaxAssign {
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpAndAssign {
  enum { tag = PETE_BinaryUseLeftTag };
};

struct OpOrAssign {
  enum { tag = PETE_BinaryUseLeftTag };
};


///////////////////////////////////////////////////////////////////////////
//
// OPERATOR()
//
///////////////////////////////////////////////////////////////////////////

#if defined(IPPL_USE_PARTIAL_SPECIALIZATION)

template<class T, class TP, unsigned Dim>
struct PETEUnaryReturn< Vektor<T, Dim>, OpParens<TP> > {
  typedef T type;
};

template<class T, class TP, unsigned Dim>
struct PETEUnaryReturn< AntiSymTenzor<T, Dim>, OpParens<TP> > {
  typedef T type;
};

template<class T, class TP, unsigned Dim>
struct PETEUnaryReturn< SymTenzor<T, Dim>, OpParens<TP> > {
  typedef T type;
};

template<class T, class TP, unsigned Dim>
struct PETEUnaryReturn< Tenzor<T, Dim>, OpParens<TP> > {
  typedef T type;
};

#else

// More Yucko non-partial-specialization stuff...

#define _PETE_DEFINE_PARENS_RETURNS_(T,D)                                   \
template<> struct PETEUnaryReturn<class Vektor<T,D>,OpParens<int> >         \
  { typedef T type; };                                                      \
template<> struct                                                           \
  PETEUnaryReturn<class Tenzor<T,D>,OpParens<std::pair<int,int> > >         \
  { typedef T type; };                                                      \
template<> struct                                                           \
  PETEUnaryReturn<class AntiSymTenzor<T,D>,OpParens<std::pair<int,int> > >  \
  { typedef T type; };                                                      \
template<> struct                                                           \
  PETEUnaryReturn<class SymTenzor<T,D>,OpParens<std::pair<int,int> > >      \
  { typedef T type; };

_PETE_DEFINE_PARENS_RETURNS_(short,1)
_PETE_DEFINE_PARENS_RETURNS_(short,2)
_PETE_DEFINE_PARENS_RETURNS_(short,3)
_PETE_DEFINE_PARENS_RETURNS_(int,1)
_PETE_DEFINE_PARENS_RETURNS_(int,2)
_PETE_DEFINE_PARENS_RETURNS_(int,3)
_PETE_DEFINE_PARENS_RETURNS_(long,1)
_PETE_DEFINE_PARENS_RETURNS_(long,2)
_PETE_DEFINE_PARENS_RETURNS_(long,3)
_PETE_DEFINE_PARENS_RETURNS_(float,1)
_PETE_DEFINE_PARENS_RETURNS_(float,2)
_PETE_DEFINE_PARENS_RETURNS_(float,3)
_PETE_DEFINE_PARENS_RETURNS_(double,1)
_PETE_DEFINE_PARENS_RETURNS_(double,2)
_PETE_DEFINE_PARENS_RETURNS_(double,3)
_PETE_DEFINE_PARENS_RETURNS_(dcomplex,1)
_PETE_DEFINE_PARENS_RETURNS_(dcomplex,2)
_PETE_DEFINE_PARENS_RETURNS_(dcomplex,3)

#undef _PETE_DEFINE_PARENS_RETURNS_

#endif


///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif // IPPL_TYPE_COMPUTATIONS_H

/***************************************************************************
 * $RCSfile: IpplTypeComputations.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:28 $
 * IPPL_VERSION_ID: $Id: IpplTypeComputations.h,v 1.1.1.1 2003/01/23 07:40:28 adelmann Exp $ 
 ***************************************************************************/
