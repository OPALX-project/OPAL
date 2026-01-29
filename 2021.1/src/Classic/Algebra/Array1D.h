#ifndef CLASSIC_Array1D_HH
#define CLASSIC_Array1D_HH

// ------------------------------------------------------------------------
// $RCSfile: Array1D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Array1D

// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:05 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Utilities/CLRangeError.h"
#include <algorithm>
#include <ostream>


// Template class Array1D<T>
// ------------------------------------------------------------------------
/// One-dimensional array.
//  A templated representation for one-dimensional arrays.
//  This class implements storage management and component access,
//  but contains no arithmetic operations.

template<class T>
class Array1D {

public:

    /// The value type of this array.
    typedef T value_type;

    /// The iterator type for the array.
    typedef T *iterator;

    /// The iterator type for constant array.
    typedef const T *const_iterator;

    /// Default constructor.
    //  Constructs an array of zero length.
    Array1D();

    /// Constructor.
    //  Reserves space for [b]n[/b] elements and leaves them undefined.
    explicit Array1D(int n);

    /// Constructor.
    //  Reserves space for [b]n[/b] elements and store [b]n[/b] copies of
    //  [b]t[/b].
    Array1D(int n, const T &t);

    Array1D(const Array1D<T> &);
    ~Array1D();
    Array1D<T> &operator=(const Array1D<T> &);

    /// Get reference to element.
    //  Return a reference to the element in position [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    T &operator()(int n);

    /// Get value of element.
    //  Return the value of the element in position [b]n[/b].
    //  Throw RangeError, if [b]n[/b] is out of range.
    const T &operator()(int n) const;

    /// Get reference to element.
    //  Return a reference to the element in position [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    T &operator[](int);

    /// Get value of element.
    //  Return the value of the element in position [b]n[/b].
    //  Result is undefined, if [b]n[/b] is out of range.
    const T &operator[](int) const;

    /// Get beginning of data.
    //  Return pointer to beginning of array.
    //  Version for non-constant objects.
    iterator begin();

    /// Get end of data.
    //  Return pointer past end of array.
    //  Version for non-constant objects.
    iterator end();

    /// Get beginning of data.
    //  Return pointer to beginning of array.
    //  Version for constant objects.
    const_iterator begin() const;

    /// Get end of data.
    //  Return pointer past end of array.
    //  Version for constant objects.
    const_iterator end() const;

    /// Get array size.
    int size() const;

    /// Change array size.
    //  Elements added are left undefined.
    void resize(int size);

protected:

    // The array size.
    int len;

    // The array data.
    T *data;
};


// Template function implementation.
// ------------------------------------------------------------------------

template<class T>
Array1D<T>::Array1D():
    len(0), data(0)
{}


template<class T>
Array1D<T>::Array1D(const Array1D<T> &array):
    len(array.len), data(new T[len]) {
    std::copy(array.begin(), array.end(), begin());
}


template<class T>
Array1D<T>::Array1D(int n):
    len(n), data(new T[len])
{}


template<class T>
Array1D<T>::Array1D(int n, const T &val):
    len(n), data(new T[len]) {
    std::fill(begin(), end(), val);
}


template<class T>
Array1D<T>::~Array1D() {
    delete [] data;
}


template<class T>
Array1D<T> &Array1D<T>::operator=(const Array1D<T> &rhs) {
    if(this != &rhs) {
        if(rhs.len > len) {
            delete [] data;
            data = new T[rhs.len];
        }

        len = rhs.len;
        std::copy(rhs.begin(), rhs.end(), begin());
    }
    return *this;
}


template<class T> inline
T &Array1D<T>::operator[](int i) {
    return data[i];
}


template<class T> inline
const T &Array1D<T>::operator[](int i) const {
    return data[i];
}


template<class T>
T &Array1D<T>::operator()(int i) {
    if(i < 0  ||  i >= size()) {
        throw CLRangeError("Array1D::operator()", "Index out of range.");
    }
    return data[i];
}


template<class T>
const T &Array1D<T>::operator()(int i) const {
    if(i < 0  ||  i >= size()) {
        throw CLRangeError("Array1D::operator()", "Index out of range.");
    }
    return data[i];
}


template<class T>
typename Array1D<T>::iterator Array1D<T>::begin() {
    return data;
}


template<class T>
typename Array1D<T>::iterator Array1D<T>::end() {
    return data + len;
}


template<class T>
typename Array1D<T>::const_iterator Array1D<T>::begin() const {
    return data;
}


template<class T>
typename Array1D<T>::const_iterator Array1D<T>::end() const {
    return data + len;
}


template<class T>
int Array1D<T>::size() const {
    return len;
}


template<class T>
void Array1D<T>::resize(int n) {
    if(len < n) {
        T *old = data;
        data = new T[n];
        std::copy(old, old + len, begin());
        delete [] old;
    }
    len = n;
}


template <class T>
std::ostream &operator<<(std::ostream &os, const Array1D<T> &v) {
    for(int i = 0; i < v.size(); ++i) {
        os << v[i] << " ";
    }

    os << std::endl;
    return os;
}

#endif // CLASSIC_Array1D_HH
