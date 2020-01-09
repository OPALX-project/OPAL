// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef VEC_H
#define VEC_H

// include files
#include "Message/Message.h"
#include "Utility/PAssert.h"

//////////////////////////////////////////////////////////////////////

template<class T, unsigned Length>
class vec
{
public:
  vec() {}
  vec(T v0);
  vec(T v0, T v1);
  vec(T v0, T v1, T v2);
#ifdef IPPL_USE_FUNKY_VEC_COPIES
  vec( const vec<T,Length>& );
  vec<T,Length>& operator=( const vec<T,Length>& );
#endif
  T& operator[](unsigned d) { return Ptr[d]; }
  const T& operator[](unsigned d) const { return Ptr[d]; }
  Message& putMessage(Message &m) {
#ifdef IPPL_USE_MEMBER_TEMPLATES
    m.put(Ptr, Ptr + Length);
#else
    putMessage(m,Ptr,Ptr+Length);
#endif
    return m;
  }
  Message& getMessage(Message &m) {
#ifdef IPPL_USE_MEMBER_TEMPLATES
    m.get_iter(Ptr);
#else
    getMessage_iter(m,Ptr);
#endif
    return m;
  }

  static T dot(const T*,const T*);
private:
  T Ptr[Length];
};

//////////////////////////////////////////////////////////////////////

#ifdef IPPL_USE_FUNKY_VEC_COPIES
//
// Implementation of copy constructor and op=.
//
// These use a template trick to unroll the loop over elements of the vec.
//

//
// A tag class to let us unroll this loop.
// It isn't entirely clear this needs to be done, but
// it is an interesting demonstration.
//
// It is just a record of how many elements to copy,
// along with a flag.  The flag has the values:
// 4: If Length>=4.  That means the interval should be split in half.
// 0..3: Call an explicit function for this length.
// 

template<unsigned Length, int Flag>
class DivideVecCopyTag {
};

//
// Copy a vector of length 4 or greater.
// Split it in half and copy each half.
// Do it this way so that the depth of inlining is O(log(Length)) 
// instead of O(Length).
//
template<class T1, class T2, unsigned L>
inline void
divide_vec_copy(T1 *p1, T2 *p2, DivideVecCopyTag<L,4> )
{
  divide_vec_copy(p1, p2, DivideVecCopyTag< L/2 , (L>=8 ? 4 : L/2)>());
  divide_vec_copy(p1+(L/2), p2+(L/2),DivideVecCopyTag<L-L/2,(L>=8?4:L-L/2)>());
}

//
// Copy a vector of length 0, 1, 2 or 3.  Just move it.
//
template<class T1, class T2>
inline void
divide_vec_copy(T1 *p1, T2 *p2, DivideVecCopyTag<3,3> )
{
  p1[0] = p2[0];
  p1[1] = p2[1];
  p1[2] = p2[2];
}

template<class T1, class T2>
inline void
divide_vec_copy(T1 *p1, T2 *p2, DivideVecCopyTag<2,2> )
{
  p1[0] = p2[0];
  p1[1] = p2[1];
}

template<class T1, class T2>
inline void
divide_vec_copy(T1 *p1, T2 *p2, DivideVecCopyTag<1,1> )
{
  *p1 = *p2;
}

template<class T1, class T2>
inline void
divide_vec_copy(T1 *,T2 *, DivideVecCopyTag<0,0> )
{
}

//
// The copy ctor and op= just call divide_vec_copy to do the copy
//

template<class T, unsigned L>
inline 
vec<T,L>::vec( const vec<T,L>& v )
{
  divide_vec_copy( Ptr , v.Ptr , DivideVecCopyTag<L,( L>=4 ? 4 : L)>() );
}

template<class T, unsigned L>
inline vec<T,L>&
vec<T,L>::operator=( const vec<T,L>& v )
{
  if ( this != &v )
    divide_vec_copy( Ptr , v.Ptr , DivideVecCopyTag<L,( L>=4 ? 4 : L)>() );
  return *this;
}

#endif // IPPL_USE_FUNKY_VEC_COPIES

//////////////////////////////////////////////////////////////////////

template<class T, unsigned Length>
inline 
vec<T,Length>::vec(T v0)
{
  CTAssert(Length==1);
  Ptr[0] = v0;
}

template<class T, unsigned Length>
inline 
vec<T,Length>::vec(T v0, T v1)
{
  CTAssert(Length==2);
  Ptr[0] = v0;
  Ptr[1] = v1;
}

template<class T, unsigned Length>
inline 
vec<T,Length>::vec(T v0, T v1, T v2)
{
  CTAssert(Length==3);
  Ptr[0] = v0;
  Ptr[1] = v1;
  Ptr[2] = v2;
}

//////////////////////////////////////////////////////////////////////

//
// Define a global function for taking the dot product between two 
// short arrays of objects of type T.
//
template<class T, unsigned Length>
inline T
vec<T,Length>::dot(const T* l, const T* r)
{
  T ret = l[0]*r[0];
  for (int i=1; i<Length; ++i)
    ret += l[i]*r[i];
  return ret;
}

//////////////////////////////////////////////////////////////////////

#endif // VEC_H

/***************************************************************************
 * $RCSfile: Vec.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:34 $
 * IPPL_VERSION_ID: $Id: Vec.h,v 1.1.1.1 2003/01/23 07:40:34 adelmann Exp $ 
 ***************************************************************************/
