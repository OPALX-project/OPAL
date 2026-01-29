#ifndef CLASSIC_FSlice_HH
#define CLASSIC_FSlice_HH

// ------------------------------------------------------------------------
// $RCSfile: FSlice.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FSlice
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <iterator>

// Class FSlice.
// ------------------------------------------------------------------------
/// An iterator permitting to iterate with a stride different from 1.

template<class T, int S>
class FSlice {

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

    /// Constructor for array.
    explicit FSlice(T *);

    FSlice(const FSlice &);
    FSlice &operator=(const FSlice &);
    bool operator==(const FSlice &rhs) const;
    bool operator!=(const FSlice &rhs) const;

    /// Increment (iterate forward).
    FSlice &operator++();

    /// Increment (iterate forward).
    FSlice  operator++(int);

    /// Decrement (iterate backward).
    FSlice &operator--();

    /// Decrement (iterate backward).
    FSlice  operator--(int);

    /// Increment by multiple stride.
    FSlice &operator+=(std::ptrdiff_t);

    /// Decrement by multiple stride.
    FSlice &operator-=(std::ptrdiff_t);

    /// Add multiple stride.
    FSlice operator+(std::ptrdiff_t);

    /// Subtract multiple stride.
    FSlice operator-(std::ptrdiff_t);

    /// Difference.
    difference_type operator-(const FSlice &) const;

    /// Dereference.
    T &operator*() const;

    /// Delegate.
    T &operator[](int) const;

private:

    T *cursor;
};


// Class FConstSlice.
// ------------------------------------------------------------------------
/// Constant version of FSlice

template<class T, int S>
class FConstSlice {

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

    /// Constructor from array and stride.
    FConstSlice(const T *);

    FConstSlice(const FConstSlice &rhs);
    FConstSlice &operator=(const FConstSlice &);
    bool operator==(const FConstSlice &rhs) const;
    bool operator!=(const FConstSlice &rhs) const;

    /// Increment (iterate forward).
    FConstSlice &operator++();

    /// Increment (iterate forward).
    FConstSlice  operator++(int);

    /// Decrement (iterate backward).
    FConstSlice &operator--();

    /// Decrement (iterate backward).
    FConstSlice  operator--(int);

    /// Add multiple stride and assign.
    FConstSlice &operator+=(std::ptrdiff_t);

    /// Subtract multiple stride and assign.
    FConstSlice &operator-=(std::ptrdiff_t);

    /// Add multiple stride.
    FConstSlice operator+(std::ptrdiff_t) const;

    /// Subtract multiple stride.
    FConstSlice operator-(std::ptrdiff_t) const;

    /// Get pointer difference.
    difference_type operator-(const FConstSlice &) const;

    /// Dereference.
    const T &operator*() const;

    /// Delegate.
    const T &operator[](int) const;

private:

    mutable const T *cursor;
};


// Ancillary functions.
// ------------------------------------------------------------------------

template<class T, int S>
inline std::random_access_iterator_tag
iterator_category(const FSlice<T, S> &) {
    return std::random_access_iterator_tag();
}

template<class T, int S>
inline std::random_access_iterator_tag
iterator_category(const FConstSlice<T, S> &) {
    return std::random_access_iterator_tag();
}

template<class T, int S>
inline T *value_type(const FSlice<T, S> &) {
    return (T *)(0);
}

template<class T, int S>
inline T *value_type(const FConstSlice<T, S> &) {
    return (T *)(0);
}


// Template function implementation for FSlice.
// ------------------------------------------------------------------------

template<class T, int S>
inline FSlice<T, S>::FSlice(T *array):
    cursor(array)
{}

template<class T, int S>
inline FSlice<T, S>::FSlice(const FSlice &rhs):
    cursor(rhs.cursor)
{}

template<class T, int S>
inline FSlice<T, S> &FSlice<T, S>::operator=(const FSlice &rhs) {
    cursor = rhs.cursor;
    return *this;
}

template<class T, int S>
inline bool FSlice<T, S>::operator==(const FSlice &rhs) const {
    return cursor == rhs.cursor;
}

template<class T, int S>
inline bool FSlice<T, S>::operator!=(const FSlice &rhs) const {
    return cursor != rhs.cursor;
}

template<class T, int S>
inline FSlice<T, S> &FSlice<T, S>::operator++() {
    cursor += S;
    return *this;
}

template<class T, int S>
inline FSlice<T, S> FSlice<T, S>::operator++(int) {
    FSlice tmp(*this);
    cursor += S;
    return tmp;
}

template<class T, int S>
inline FSlice<T, S> &FSlice<T, S>::operator--() {
    cursor -= S;
    return *this;
}

template<class T, int S>
inline FSlice<T, S> FSlice<T, S>::operator--(int) {
    FSlice tmp(*this);
    cursor -= S;
    return tmp;
}

template<class T, int S>
inline FSlice<T, S> &FSlice<T, S>::operator+=(std::ptrdiff_t n) {
    cursor += n * S;
    return *this;
}

template<class T, int S>
inline FSlice<T, S> &FSlice<T, S>::operator-=(std::ptrdiff_t n) {
    cursor -= n * S;
    return *this;
}

template<class T, int S>
inline FSlice<T, S> FSlice<T, S>::operator+(std::ptrdiff_t n) {
    FSlice tmp(*this);
    return tmp += n;
}

template<class T, int S>
inline FSlice<T, S> FSlice<T, S>::operator-(std::ptrdiff_t n) {
    FSlice tmp(*this);
    return tmp -= n;
}

template<class T, int S>
inline typename FSlice<T, S>::difference_type
FSlice<T, S>::operator-(const FSlice &rhs) const {
    return (cursor - rhs.cursor) / S;
}

template<class T, int S>
inline T &FSlice<T, S>::operator*() const {
    return *cursor;
}

template<class T, int S>
inline T &FSlice<T, S>::operator[](int n) const {
    return cursor[n*S];
}


// Template function implementation for FConstSlice.
// ------------------------------------------------------------------------

template<class T, int S>
inline FConstSlice<T, S>::FConstSlice(const T *array):
    cursor(array)
{}

template<class T, int S>
inline FConstSlice<T, S>::FConstSlice(const FConstSlice &rhs):
    cursor(rhs.cursor)
{}

template<class T, int S>
inline FConstSlice<T, S> &FConstSlice<T, S>::operator=(const FConstSlice &rhs) {
    cursor = rhs.cursor;
    return *this;
}

template<class T, int S>
inline bool FConstSlice<T, S>::operator==(const FConstSlice &rhs) const {
    return cursor == rhs.cursor;
}

template<class T, int S>
inline bool FConstSlice<T, S>::operator!=(const FConstSlice &rhs) const {
    return cursor != rhs.cursor;
}

template<class T, int S>
inline FConstSlice<T, S> &FConstSlice<T, S>::operator++() {
    cursor += S;
    return *this;
}

template<class T, int S>
inline FConstSlice<T, S> FConstSlice<T, S>::operator++(int) {
    FConstSlice tmp(*this);
    cursor += S;
    return tmp;
}

template<class T, int S>
inline FConstSlice<T, S> &FConstSlice<T, S>::operator--() {
    cursor -= S;
    return *this;
}

template<class T, int S>
inline FConstSlice<T, S> FConstSlice<T, S>::operator--(int) {
    FConstSlice tmp(*this);
    cursor -= S;
    return tmp;
}

template<class T, int S>
inline FConstSlice<T, S> &FConstSlice<T, S>::operator+=(std::ptrdiff_t n) {
    cursor += n * S;
    return *this;
}

template<class T, int S>
inline FConstSlice<T, S> &FConstSlice<T, S>::operator-=(std::ptrdiff_t n) {
    cursor -= n * S;
    return *this;
}

template<class T, int S> inline
FConstSlice<T, S> FConstSlice<T, S>::operator+(std::ptrdiff_t n) const {
    FConstSlice tmp(*this);
    return tmp += n;
}

template<class T, int S> inline
FConstSlice<T, S> FConstSlice<T, S>::operator-(std::ptrdiff_t n) const {
    FConstSlice tmp(*this);
    return tmp -= n;
}

template<class T, int S>
inline typename FConstSlice<T, S>::difference_type
FConstSlice<T, S>::operator-(const FConstSlice<T, S> &rhs) const {
    return (cursor - rhs.cursor) / S;
}

template<class T, int S>
inline const T &FConstSlice<T, S>::operator*() const {
    return *cursor;
}

template<class T, int S>
inline const T &FConstSlice<T, S>::operator[](int n) const {
    return cursor[n*S];
}

#endif // CLASSIC_FSlice_HH

