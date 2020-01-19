// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 ***************************************************************************/

#ifndef RNG_LATTICE_H
#define RNG_LATTICE_H

/***********************************************************************
 * 
 * class RNGLattice
 * class RNGLatticeSequence : public SequenceGen<RNGLattice>
 *
 * Lattice class that implements random number generator by just
 * providing an evenly spaced set of points between two endpoints.
 * 
 * Note that this only works properly automatically on 1 PE; with multipple PE's
 * you'll get a redundant set of points, the same on every PE, unless you take
 * care to provide appropriate different seeds for each PE and call SetSeed
 * before using it.
 * Use RNGLattice as a scalar or container element, and use
 * RNGLatticeSequence when you need a sequence of numbers to fill a container.
 *
 ***********************************************************************/

#include "Utility/SequenceGen.h"


template<class T>
class RNGLattice {

public:
  // return type
  typedef T Return_t;

  // constructor
  RNGLattice(T minval, T maxval, unsigned long numpoints,
             bool includeEndpoints = true)
    : MinVal(minval), MaxVal(maxval), CurrPoint(0), NumPoints(numpoints),
      IncludeEndpoints(includeEndpoints)
  {
    if (IncludeEndpoints) {
      if (NumPoints == 0 || NumPoints == 1) {
	Spacing = 0;
      }
      else {
	Spacing = (MaxVal - MinVal)/(NumPoints - 1);
      }
    }
    else {
      if (NumPoints == 0) {
	Spacing = 0;
      }
      else {
	Spacing = (MaxVal - MinVal)/NumPoints;
      }
    }
  }

  // copy constructor
  RNGLattice(const RNGLattice<T>& rng)
    : MinVal(rng.MinVal), MaxVal(rng.MaxVal), CurrPoint(rng.CurrPoint),
      Spacing(rng.Spacing), NumPoints(rng.NumPoints),
      IncludeEndpoints(rng.IncludeEndpoints) {}

  // destructor
  ~RNGLattice(void) {}

  // advance current point by n, modulo NumPoints
  inline void AdvanceSeed(unsigned long n = 0) {
    CurrPoint = (CurrPoint + n) % NumPoints;
    return;
  }

  // set current point to nth point, modulo NumPoints
  inline void SetSeed(unsigned long n = 0) {
    CurrPoint = n % NumPoints;
  }

  // get current point 
  inline unsigned long GetSeed(void) const { return CurrPoint; }

  // return the next pseudo-random number (from MinVal ... MaxVal)
  inline Return_t GetRandom(void) const {
    Return_t currVal;
    if (IncludeEndpoints)
      currVal = MinVal + CurrPoint * Spacing;
    else
      currVal = MinVal + (CurrPoint + 0.5) * Spacing;
    CurrPoint++;
    if (CurrPoint == NumPoints) CurrPoint = 0;
    return currVal;
  }

  // pseudonym for GetRandom()
  inline Return_t operator()(void) const { return GetRandom(); }

  // conversion to Return_t, same as GetRandom()
  inline operator Return_t() const { return GetRandom(); }

  // return the period of the RNG
  unsigned long GetRandMax(void) const { return NumPoints; }

private:
  T MinVal, MaxVal, Spacing;
  mutable unsigned long CurrPoint;
  unsigned long NumPoints;
  bool IncludeEndpoints;
};

RNG_BASIC_MATH(RNGLattice<double>)

// A version of RNGLattice with extra constructors to make using this
// class easier.  This is the version that people should use to fill
// containers with a random number sequence in an expression.  This
// class is PETE-aware via its inheritance from SequenceGen.

template<class T>
class RNGLatticeSequence : public SequenceGen< RNGLattice<T> > {

public:
  // constructor
  RNGLatticeSequence(T minval, T maxval, unsigned long numpoints,
                     bool includeEndpoints = true)
    : SequenceGen< RNGLattice<T> >(RNGLattice<T>(minval,maxval,numpoints,
                                                 includeEndpoints)) {}

  // copy constructor
  RNGLatticeSequence(const RNGLatticeSequence<T>& rngseq)
    : SequenceGen< RNGLattice<T> >(rngseq.getGenerator()) {}

  // destructor
  ~RNGLatticeSequence(void) {}

  // wrappers around RNG generator functions
  inline void     AdvanceSeed(unsigned long adv = 0) {
      this->getGenerator().AdvanceSeed(adv);
  }
  inline void     SetSeed(unsigned long seed) { this->getGenerator().SetSeed(seed); }
  inline unsigned long GetSeed(void) const { return this->getGenerator().GetSeed(); }
  inline
  typename RNGLattice<T>::Return_t
  GetRandom(void) { return this->getGenerator().GetRandom(); }
  inline
  typename RNGLattice<T>::Return_t
  operator()(void) { return this->getGenerator().GetRandom(); }
  unsigned long GetRandMax(void) const { return this->getGenerator().GetRandMax(); }
};


#endif // RNG_LATTICE_H

