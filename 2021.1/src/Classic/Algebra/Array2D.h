#ifndef CLASSIC_Array2D_HH
#define CLASSIC_Array2D_HH

// ------------------------------------------------------------------------
// $RCSfile: Array2D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Array2D
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "Algebra/SliceIterator.h"
#include "Utilities/CLRangeError.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <iosfwd>


// Template class Array2D<T>
// ------------------------------------------------------------------------
/// Two-dimensional array.
//  A templated representation for 2-dimensional arrays.
//  This class implements storage management and component access,
//  as well as access to rows and columns, but contains no arithmetic
//  operations.

template<class T>
class Array2D {

public:

    /// The value type of the array.
    typedef T value_type;

    /// The iterator type for sequential access to all elements.
    typedef T *iterator;

    /// The iterator type for sequential access to all elements
    //  of a constant array.
    typedef const T *const_iterator;

    /// The iterator type for access by rows.
    typedef T *row_iterator;

    /// The iterator type for access by rows for a constant array.
    typedef const T *const_row_iterator;

    /// The iterator type for access by columns.
    typedef SliceIterator<T> col_iterator;

    /// The iterator type for access by columns for a constant array.
    typedef ConstSliceIterator<T> const_col_iterator;

    /// Default constructor.
    //  Constructs array of zero rows and zero columns.
    Array2D();

    /// Constructor.
    //  Reserves [b]rows[/b] times [b]cols[/b] elements and leaves
    //  them undefined.
    Array2D(int rows, int cols);

    /// Constructor.
    //  Reserves [b]rows[/b] x [b]cols[/b] elements and sets them to [b]t[/b].
    Array2D(int rows, int cols, const T &t);

    Array2D(const Array2D &);
    ~Array2D();
    Array2D<T>& operator=(const Array2D<T> &);

    /// Get element reference.
    //  Return a reference to element in row [b]r[/b] and column [b]c[/b].
    T &operator()(int r, int c);

    /// Get element value.
    //  Return the value of element in row [b]r[/b] and column [b]c[/b].
    const T &operator()(int r, int c) const;

    /// Get number of rows.
    int nrows() const;

    /// Get number of columns.
    int ncols() const;

    /// Get total size (rows times columns).
    int size() const;


    /// Get pointer to beginning of data.
    //  Treat the array as one-dimensional.
    //  Version for non-constant objects.
    iterator begin();

    /// Get pointer past end of data.
    //  Treat the array as one-dimensional.
    //  Version for non-constant objects.
    iterator end();

    /// Get pointer to beginning of data.
    //  Treat the array as one-dimensional.
    //  Version for constant objects.
    const_iterator begin() const;

    /// Get pointer past end of data.
    //  Treat the array as one-dimensional.
    //  Version for constant objects.
    const_iterator end() const;


    /// Get row iterator.
    //  Return pointer to beginning of row r.
    //  Throw RangeError, if row is out of range.
    //  Version for non-constant objects.
    row_iterator row_begin(int r);

    /// Get row iterator.
    //  Return pointer past end of row r.
    //  Throw RangeError, if row is out of range.
    //  Version for non-constant objects.
    row_iterator row_end(int r);

    /// Get row iterator.
    //  Return pointer to beginning of constant row r.
    //  Throws RangeError, if row is out of range.
    //  Version for constant objects.
    const_row_iterator row_begin(int r) const;

    /// Get row iterator.
    //  Return pointer past end of constant row r.
    //  Throws RangeError, if row is out of range.
    //  Version for constant objects.
    const_row_iterator row_end(int r) const;


    /// Get row iterator.
    //  Return pointer to beginning of row r.
    //  Result is undefined, if r is out of range.
    //  Version for non-constant objects.
    row_iterator operator[](int r);

    /// Get row iterator.
    //  Return pointer to beginning of row r.
    //  Result is undefined, if r is out of range.
    //  Version for constant objects.
    const_row_iterator operator[](int r) const;


    /// Get column iterator.
    //  Return pointer to beginning of column c.
    //  Throw RangeError if c is out of range.
    //  Version for non-constant objects.
    col_iterator col_begin(int c);

    /// Get column iterator.
    //  Return pointer past end of column c.
    //  Throw RangeError if c is out of range.
    //  Version for non-constant objects.
    col_iterator col_end(int c);

    /// Get column iterator.
    //  Return pointer to beginning of column c.
    //  Throw RangeError if c is out of range.
    //  Version for constant objects.
    const_col_iterator col_begin(int c) const;

    /// Get column iterator.
    //  Return pointer past end of column c.
    //  Throw RangeError if c is out of range.
    //  Version for constant objects.
    const_col_iterator col_end(int c) const;


    /// Fetch column.
    //  Store the column [b]c[/b] into the array [b]toArray[/b].
    //  Throw RangeError if c is out of range.
    void getColumn(Array1D<T> &toArray, int c) const;

    /// Fetch row.
    //  Store the row [b]r[/b] into the array [b]toArray[/b].
    //  Throw RangeError if r is out of range.
    void getRow(Array1D<T> &toArray, int r) const;

    /// Store column.
    //  Fill the column [b]c[/b] from the array [b]fromArray[/b].
    //  Throw RangeError if c is out of range.
    void putColumn(const Array1D<T> &fromArray, int c);

    /// Store row.
    //  Fill the row [b]r[/b] from the array [b]fromArray[/b].
    //  Throw RangeError if r is out of range.
    void putRow(const Array1D<T> &fromArray, int r);

    /// Exchange columns.
    //  Exchange the columns [b]c1[/b] and [b]c2[/b].


    //  Throw RangeError, if either index is out of range.
    void swapColumns(int c1, int c2);

    /// Exchange rows.
    //  Exchange the rows [b]r1[/b] and [b]r2[/b].
    //  Throw RangeError, if either index is out of range.
    void swapRows(int row1, int row2);

protected:

    // Array dimensions.
    int rows;
    int cols;
    int len;  // rows * cols

    // Array data.
    T *data;
};


// Template implementation.
// ------------------------------------------------------------------------

template <class T>
inline Array2D<T>::Array2D():
    rows(0), cols(0), len(0), data(0)
{}


template <class T>
inline Array2D<T>::Array2D(const Array2D<T>& array):
    rows(array.rows), cols(array.cols), len(array.len), data(new T[len]) {
    std::copy(array.begin(), array.end(), begin());
}


template <class T>
inline Array2D<T>::Array2D(int r, int c):
    rows(r), cols(c), len(r *c), data(new T[len])
{}


template <class T>
inline Array2D<T>::Array2D(int r, int c, const T &val):
    rows(r), cols(c), len(r *c), data(new T[len]) {
    std::fill(begin(), end(), val);
}


template <class T>
inline Array2D<T>::~Array2D() {
    delete [] data;
}


template <class T>
inline Array2D<T>& Array2D<T>::operator=(const Array2D<T>& rhs)

{
    if(this != &rhs) {
        if(rhs.len > len) {
            delete [] data;
            data = new T[rhs.len];
        }
        rows = rhs.rows;
        cols = rhs.cols;
        len  = rhs.len;
        std::copy(rhs.begin(), rhs.end(), begin());
    }
    return *this;
}


template <class T>
inline T &Array2D<T>::operator()(int r, int c) {
    if(r < 0  ||  r >= rows  ||  c < 0  ||  c >= ncols()) {
        throw CLRangeError("Array2d::operator()", "Index out of range.");
    }
    return data[cols*r+c];
}


template <class T>
inline const T &Array2D<T>::operator()(int r, int c) const {
    if(r < 0  ||  r >= rows  ||  c < 0  ||  c >= ncols()) {
        throw CLRangeError("Array2d::operator()", "Index out of range.");
    }
    return data[cols*r+c];
}


template <class T>
inline int Array2D<T>::nrows() const {
    return rows;
}


template <class T>
inline int Array2D<T>::ncols() const {
    return cols;
}


template <class T>
inline int Array2D<T>::size() const {
    return len;
}


template <class T>
inline typename Array2D<T>::iterator Array2D<T>::begin() {
    return data;
}


template <class T>
inline typename Array2D<T>::iterator Array2D<T>::end() {
    return data + len;
}


template <class T>
inline typename Array2D<T>::const_iterator Array2D<T>::begin() const {
    return data;
}


template <class T>
inline typename Array2D<T>::const_iterator Array2D<T>::end() const {
    return data + len;
}


template <class T>
typename Array2D<T>::row_iterator Array2D<T>::row_begin(int r) {
    if(r >= rows) {
        throw CLRangeError("Array2D::row_begin()", "Row index out of range.");
    }
    return &data[cols*r];
}


template <class T>
typename Array2D<T>::row_iterator Array2D<T>::row_end(int r) {
    if(r >= rows) {
        throw CLRangeError("Array2D::row_end()", "Row index out of range.");
    }
    return &data[cols*(r+1)];
}


template <class T>
typename Array2D<T>::const_row_iterator Array2D<T>::row_begin(int r) const {
    if(r >= rows) {
        throw CLRangeError("Array2D::row_begin()", "Row index out of range.");
    }
    return &data[cols*r];
}


template <class T>
typename Array2D<T>::const_row_iterator Array2D<T>::row_end(int r) const {
    if(r >= rows) {
        throw CLRangeError("Array2D::row_end()", "Row index out of range.");
    }
    return &data[cols*(r+1)];
}


template <class T>
typename Array2D<T>::row_iterator Array2D<T>::operator[](int r) {
    return &data[cols*r];
}


template <class T>
typename Array2D<T>::const_row_iterator Array2D<T>::operator[](int r) const {
    return &data[cols*r];
}


template <class T>
typename Array2D<T>::col_iterator Array2D<T>::col_begin(int c) {
    if(c >= cols) {
        throw CLRangeError("Array2D::col_begin()", "Column index out of range.");
    }
    return col_iterator(data + c, cols);
}


template <class T>
typename Array2D<T>::col_iterator Array2D<T>::col_end(int c) {
    if(c >= cols) {
        throw CLRangeError("Array2D::col_end()", "Column index out of range.");
    }
    return col_iterator(data + len + c, cols);
}


template <class T>
typename Array2D<T>::const_col_iterator Array2D<T>::col_begin(int c) const {
    if(c >= cols) {
        throw CLRangeError("Array2D::col_begin()", "Column index out of range.");
    }
    return const_col_iterator(data + c, cols);
}


template <class T>
typename Array2D<T>::const_col_iterator Array2D<T>::col_end(int c) const {
    if(c >= cols) {
        throw CLRangeError("Array2D::col_end()", "Column index out of range.");
    }
    return const_col_iterator(data + len + c, cols);
}


template <class T>
void Array2D<T>::getColumn(Array1D<T> &toArray, int c) const {
    if(toArray.size() != rows) toArray = Array1D<T>(rows);
    std::copy(col_begin(c), col_end(c), toArray.begin());
}


template <class T>
void Array2D<T>::getRow(Array1D<T> &toArray, int r) const {
    if(toArray.size() != cols) toArray = Array1D<T>(cols);
    std::copy(row_begin(r), row_end(r), toArray.begin());
}


template <class T>
void Array2D<T>::putColumn(const Array1D<T> &fromArray, int c) {
    if(fromArray.size() != rows) {
        throw SizeError("Array2D::putColumn()", "Dimensions inconsistent.");
    }
    std::copy(fromArray.begin(), fromArray.end(), col_begin(c));
}


template <class T>
void Array2D<T>::putRow(const Array1D<T> &fromArray, int r) {
    if(fromArray.size() != cols) {
        throw SizeError("Array2D::putRow()", "Dimensions inconsistent.");
    }
    std::copy(fromArray.begin(), fromArray.end(), row_begin(r));
}


template <class T>
void Array2D<T>::swapColumns(int col1, int col2) {
    std::swap_ranges(col_begin(col1), col_end(col1), col_begin(col2));
}


template <class T>
void Array2D<T>::swapRows(int row1, int row2) {
    std::swap_ranges(row_begin(row1), row_end(row1), row_begin(row2));
}


template <class T>
std::ostream &operator<<(std::ostream &os, const Array2D<T> &v) {
    for(int i = 0; i < v.nrows(); ++i) {
        for(int j = 0; j < v.ncols(); ++j) {
            os << v(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}


#endif // CLASSIC_Array2D_HH
