#ifndef CLASSIC_SliceIterator_HH
#define CLASSIC_SliceIterator_HH

// ------------------------------------------------------------------------
// $RCSfile: SliceIterator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: SliceIterator
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include <iterator>

// Class SliceIterator.
// ------------------------------------------------------------------------
/// Iterator for array slice.
//  An iterator permitting to iterate on a non-constant array
//  with a stride different from 1.

template<class T>
class SliceIterator {

public:

    /// The iterator tag, taken from the standard template library.
    typedef std::random_access_iterator_tag iterator_category;

    /// The value type.
    typedef T value_type;

    /// The pointer difference type.
    typedef std::ptrdiff_t difference_type;

    /// The pointer type.
    typedef T *pointer;

    /// The reference type.
    typedef T &reference;

    /// Constructor.
    //  Construct iterator for [b]array[/b], given a [b]stride[/b].
    SliceIterator(T *array, std::ptrdiff_t stride);

    SliceIterator(const SliceIterator<T> &);
    SliceIterator<T> &operator=(const SliceIterator<T> &);
    bool operator==(const SliceIterator<T> &rhs) const;
    bool operator!=(const SliceIterator<T> &rhs) const;

    /// Increment (iterate forward).
    SliceIterator<T> &operator++();

    /// Increment (iterate forward).
    SliceIterator<T>  operator++(int);

    /// Decrement (iterate backward).
    SliceIterator<T> &operator--();

    /// Decrement (iterate backward).
    SliceIterator<T>  operator--(int);

    /// Increment by multiple stride.
    SliceIterator<T> &operator+=(std::ptrdiff_t);

    /// Decrement by multiple stride.
    SliceIterator<T> &operator-=(std::ptrdiff_t);

    /// Add multiple stride.
    SliceIterator<T> operator+(std::ptrdiff_t);

    /// Subtract multiple stride.
    SliceIterator<T> operator-(std::ptrdiff_t);

    /// Difference.
    difference_type operator-(const SliceIterator<T> &) const;

    /// Dereference.
    T &operator*() const;

    /// Delegate.
    T &operator[](int) const;

private:

    T *cursor;
    std::ptrdiff_t stride;
};


// Class ConstSliceIterator.
// ------------------------------------------------------------------------
/// Iterator for array slice.
//  An iterator permitting to iterate on a constant array
//  with a stride different from 1.

template<class T>
class ConstSliceIterator {

public:

    /// The iterator tag, taken from the standard template library.
    typedef std::random_access_iterator_tag iterator_category;

    /// The value type.
    typedef const T value_type;

    /// The pointer difference type.
    typedef std::ptrdiff_t difference_type;

    /// The pointer type.
    typedef const T *pointer;

    /// The reference type.
    typedef const T &reference;

    /// Constructor.
    //  Construct iterator for [b]array[/b], given a [b]stride[/b].
    ConstSliceIterator(const T *array, std::ptrdiff_t stride);

    ConstSliceIterator(const ConstSliceIterator<T> &rhs);
    ConstSliceIterator<T> &operator=(const ConstSliceIterator<T>&);
    bool operator==(const ConstSliceIterator<T> &rhs) const;
    bool operator!=(const ConstSliceIterator<T> &rhs) const;

    /// Increment (iterate forward).
    ConstSliceIterator<T> &operator++();

    /// Increment (iterate forward).
    ConstSliceIterator<T>  operator++(int);

    /// Decrement (iterate backward).
    ConstSliceIterator<T> &operator--();

    /// Decrement (iterate backward).
    ConstSliceIterator<T>  operator--(int);

    /// Add multiple stride and assign.
    ConstSliceIterator<T> &operator+=(std::ptrdiff_t);

    /// Subtract multiple stride and assign.
    ConstSliceIterator<T> &operator-=(std::ptrdiff_t);

    /// Add multiple stride.
    ConstSliceIterator<T> operator+(std::ptrdiff_t) const;

    /// Subtract multiple stride.
    ConstSliceIterator<T> operator-(std::ptrdiff_t) const;

    /// Get pointer difference.
    difference_type operator-(const ConstSliceIterator<T> &) const;

    /// Dereference.
    const T &operator*() const;

    /// Delegate.
    const T &operator[](int) const;

private:

    mutable const T *cursor;
    mutable std::ptrdiff_t stride;
};


// Ancillary functions.
// ------------------------------------------------------------------------

template<class T>
inline std::random_access_iterator_tag
iterator_category(const SliceIterator<T> &) {
    return std::random_access_iterator_tag();
}

template<class T>
inline std::random_access_iterator_tag
iterator_category(const ConstSliceIterator<T> &) {
    return std::random_access_iterator_tag();
}

template<class T>
inline T *value_type(const SliceIterator<T> &) {
    return (T *)(0);
}

template<class T>
inline T *value_type(const ConstSliceIterator<T> &) {
    return (T *)(0);
}


// Template function implementation for SliceIterator.
// ------------------------------------------------------------------------

template<class T>
inline SliceIterator<T>::SliceIterator(T *array, std::ptrdiff_t n):
    cursor(array), stride(n)
{}

template<class T>
inline SliceIterator<T>::SliceIterator(const SliceIterator<T> &rhs):
    cursor(rhs.cursor), stride(rhs.stride)
{}

template<class T>
inline SliceIterator<T> &SliceIterator<T>::operator=
(const SliceIterator<T> &rhs) {
    cursor = rhs.cursor;
    stride = rhs.stride;
    return *this;
}

template<class T>
inline bool SliceIterator<T>::operator==(const SliceIterator<T> &rhs) const {
    return cursor == rhs.cursor;
}

template<class T>
inline bool SliceIterator<T>::operator!=(const SliceIterator<T> &rhs) const {
    return cursor != rhs.cursor;
}

template<class T>
inline SliceIterator<T> &SliceIterator<T>::operator++() {
    cursor += stride;
    return *this;
}

template<class T>
inline SliceIterator<T> SliceIterator<T>::operator++(int) {
    SliceIterator<T> tmp(*this);
    cursor += stride;
    return tmp;
}

template<class T>
inline SliceIterator<T> &SliceIterator<T>::operator--() {
    cursor -= stride;
    return *this;
}

template<class T>
inline SliceIterator<T> SliceIterator<T>::operator--(int) {
    SliceIterator<T> tmp(*this);
    cursor -= stride;
    return tmp;
}

template<class T>
inline SliceIterator<T> &SliceIterator<T>::operator+=(std::ptrdiff_t n) {
    cursor += n * stride;
    return *this;
}

template<class T>
inline SliceIterator<T> &SliceIterator<T>::operator-=(std::ptrdiff_t n) {
    cursor -= n * stride;
    return *this;
}

template<class T>
inline SliceIterator<T> SliceIterator<T>::operator+(std::ptrdiff_t n) {
    SliceIterator<T> tmp(*this);
    return tmp += n;
}

template<class T>
inline SliceIterator<T> SliceIterator<T>::operator-(std::ptrdiff_t n) {
    SliceIterator<T> tmp(*this);
    return tmp -= n;
}

template<class T>
inline typename SliceIterator<T>::difference_type
SliceIterator<T>::operator-(const SliceIterator<T> &rhs) const {
    return (cursor - rhs.cursor) / stride;
}

template<class T>
inline T &SliceIterator<T>::operator*() const {
    return *cursor;
}

template<class T>
inline T &SliceIterator<T>::operator[](int n) const {
    return cursor[n*stride];
}


// Template function implementation for ConstSliceIterator.
// ------------------------------------------------------------------------

template<class T>
inline ConstSliceIterator<T>::ConstSliceIterator
(const T *array, std::ptrdiff_t n):
    cursor(array), stride(n)
{}

template<class T>
inline ConstSliceIterator<T>::ConstSliceIterator
(const ConstSliceIterator<T> &rhs):
    cursor(rhs.cursor), stride(rhs.stride)
{}

template<class T>
inline ConstSliceIterator<T> &ConstSliceIterator<T>::operator=
(const ConstSliceIterator<T> &rhs) {
    cursor = rhs.cursor;
    stride = rhs.stride;
    return *this;
}

template<class T>
inline bool ConstSliceIterator<T>::operator==
(const ConstSliceIterator<T> &rhs) const {
    return cursor == rhs.cursor;
}

template<class T>
inline bool ConstSliceIterator<T>::operator!=
(const ConstSliceIterator<T> &rhs) const {
    return cursor != rhs.cursor;
}

template<class T>
inline ConstSliceIterator<T> &ConstSliceIterator<T>::operator++() {
    cursor += stride;
    return *this;
}

template<class T>
inline ConstSliceIterator<T> ConstSliceIterator<T>::operator++(int) {
    ConstSliceIterator<T> tmp(*this);
    cursor += stride;
    return tmp;
}

template<class T>
inline ConstSliceIterator<T> &ConstSliceIterator<T>::operator--() {
    cursor -= stride;
    return *this;
}

template<class T>
inline ConstSliceIterator<T> ConstSliceIterator<T>::operator--(int) {
    ConstSliceIterator<T> tmp(*this);
    cursor -= stride;
    return tmp;
}

template<class T>
inline ConstSliceIterator<T> &ConstSliceIterator<T>::operator+=
(std::ptrdiff_t n) {
    cursor += n * stride;
    return *this;
}

template<class T>
inline ConstSliceIterator<T> &ConstSliceIterator<T>::operator-=
(std::ptrdiff_t n) {
    cursor -= n * stride;
    return *this;
}

template<class T> inline
ConstSliceIterator<T> ConstSliceIterator<T>::operator+(std::ptrdiff_t n) const {
    ConstSliceIterator<T> tmp(*this);
    return tmp += n;
}

template<class T> inline
ConstSliceIterator<T> ConstSliceIterator<T>::operator-(std::ptrdiff_t n) const {
    ConstSliceIterator<T> tmp(*this);
    return tmp -= n;
}

template<class T>
inline typename ConstSliceIterator<T>::difference_type
ConstSliceIterator<T>::operator-(const ConstSliceIterator<T> &rhs) const {
    return (cursor - rhs.cursor) / stride;
}

template<class T>
inline const T &ConstSliceIterator<T>::operator*() const {
    return *cursor;
}

template<class T>
inline const T &ConstSliceIterator<T>::operator[](int n) const {
    return cursor[n*stride];
}

#endif // CLASSIC_SliceIterator_HH

