#ifndef CLASSIC_FTps_HH
#define CLASSIC_FTps_HH

// ------------------------------------------------------------------------
// $RCSfile: FTps.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.4.2.9 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FTps<T,V>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:06:20 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FTpsData.h"
#include <iosfwd>
#include <list>
#include <string>

struct TpsSubstitution;

template <class T> class Array1D;
template <class T, int M, int N> class FMatrix;
template <int N> class FMonomial;
template <class T, int N> class FTpsRep;
template <class T, int N> class FVector;
template <class T, int N> class FVps;


// Template class FTps<T,N>.
// ------------------------------------------------------------------------
/// Truncated power series in N variables of type T.
//  All divide operations throw DivideError
//  if the constant part of the divisor is zero.

template <class T, int N>
class FTps {

public:

    /// Default constructor.
    //  Constructs zero value.
    FTps();

    /// Copy constructor.
    FTps(const FTps &x);

    /// Constructor.
    //  Define the internal orders.
    FTps(int minOrder, int maxOrder, int trcOrder);

    /// Conversion.
    FTps(const T &);

    /// Conversion.
    FTps(int);

    /// Destructor.
    ~FTps();

    /// Assign.
    FTps &operator=(const FTps &y);

    /// Convert and assign.
    FTps &operator=(const T &y);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the Giorgilli index.
    //  Return zero if index is out of range.
    const T getCoefficient(int index) const;

    /// Set coefficient.
    //  Assign value of the coefficient denoted by the Giorgilli index.
    //  Ignore if index is out of range.
    void setCoefficient(int index, const T &value);

    /// Get coefficient.
    //  Return value of the coefficient denoted by monomial exponents.
    //  Return zero if monomial order lies outside range [minOrd,maxOrd].
    const T getCoefficient(const FMonomial<N> &monomial) const;

    /// Set coefficient.
    //  Assign value of the coeffient denoted by monomial exponents.
    //  Ignore if monomial order exceeds truncation order.
    void setCoefficient(const FMonomial<N> &monomial, const T &value);

    /// Get coefficient.
    //  Return value of the coefficient denoted by the Giorgilli index.
    //  Result undefined for index not in range [orderStart(minOrd),orderEnd(maxOrd)).
    inline const T operator[](int index) const;

    /// Get (Set) coefficient.
    //  Return a reference to the coefficient denoted by the Giorgilli index.
    //  Result undefined for index not in range [orderStart(minOrd),orderEnd(maxOrd)).
    //  Makes FTps unique(), so use as an lvalue is allowed.
    inline T &operator[](int index);

    /// Get coefficient.
    //  Return value of the coefficient denoted by monomial exponents.
    //  Result undefined if monomial order lies outside range [minOrd,maxOrd].
    const T operator[](const FMonomial<N> &monomial) const;

    /// Get (Set) coefficient.
    //  Return a reference to the coefficient denoted by monomial exponents.
    //  Result undefined for monomial order not in range [minOrd,maxOrd].
    //  Makes FTps unique(), so use as an lvalue is allowed.
    T &operator[](const FMonomial<N> &monomial);

    /// Return beginning of monomial array.
    T *begin() const { return itsRep->begin(); }

    /// Return end of monomial array.
    T *end() const { return itsRep->end(); }

    /// Return beginning of coefficient storage for given order.
    T *begin(int order) const
    { return itsRep->data + FTpsData<N>::orderStart(order); }

    /// Return end of coefficient storage for given order.
    T *end(int order) const
    { return itsRep->data + FTpsData<N>::orderEnd(order); }

    /// Get number of variables.
    static int getVariables()
    { return N; }

    /// Get total number of coefficients.
    int getSize() const
    { return itsRep->len; }

    /// Get number of coefficients of degree [b]order[/b] or less.
    static int getSize(int order);

    /// Get index at which [b]order[/b] starts.
    static int orderStart(int order)
    { return FTpsData<N>::orderStart(order); }

    /// Get one plus index at which [b]order[/b] ends.
    static int orderEnd(int order)
    { return FTpsData<N>::orderEnd(order); }

    // Number of coefficients of degree "order".
    static int orderLength(int order)
    { return FTpsData<N>::orderLength(order); }

    // Number of coefficients such that "orderL" <= degree <= "orderH".
    static int orderLength(int orderL, int orderH)
    { return FTpsData<N>::orderLength(orderL, orderH); }

    /// Get exponents for given index.
    static const FMonomial<N> &getExponents(int index);

    /// Get Giorgilli index for monomial.
    static int getIndex(const FMonomial<N> &mono);

    /// Get minimum order.
    int getMinOrder() const
    { return itsRep->minOrd; }

    /// Set minimum order.
    //  If necessary, this function will insert zeroes or modify the
    //  maximum order.  It will not modify the truncation order.
    void setMinOrder(int order);

    /// Get maximum order.
    int getMaxOrder() const
    { return itsRep->maxOrd; }

    /// Set maximum order.
    //  If necessary, this function will insert zeroes or modify the minimum
    //  or maximum orders.  It will not modify the truncation order.
    void setMaxOrder(int order);

    /// Get truncation order.
    int getTruncOrder() const
    { return itsRep->trcOrd; }

    /// Set truncation order.
    //  This is the only function that may increase the truncation order of an FTps.
    void setTruncOrder(int order);

    /// Return the global truncation order.
    static int getGlobalTruncOrder()
    { return globalTruncOrder; }

    /// Set the global truncation order.
    static void setGlobalTruncOrder(int order);

    /// Index array for products of monomial "index".
    static inline const Array1D<int> &getProductArray(int index);

    /// List of variables contained in monomial "index".
    static inline const Array1D<int> &getVariableList(int index);

    /// Return the substitution table.
    static inline const Array1D<TpsSubstitution> &getSubTable();

    /// Extract given range of orders, with truncation.
    //  If the argument range [minOrder,maxOrder] does not overlap the range of
    //  min and max orders of the calling FTps, then filter() returns a zero FTps
    //  with the requested min and max orders.
    FTps filter(int minOrder, int maxOrder, int trcOrder = EXACT) const;

    /// Truncate.
    //  Limit the truncation order to [b]trunc[/b].
    FTps truncate(int trunc);

    /// Make variable.
    //  Construct the variable identified by the index [b]var[/b].
    static FTps makeVariable(int var);

    /// Make power.
    //  Construct [b]power[/b] of variable [b]var[/b].
    static FTps makeVarPower(int var, int power);

    /// Make monomial.
    //  Construct the monomial with Giorgilli index [b]index[/b] and coefficient [b]t[/b].
    static FTps makeMonomial(int index, const T &t);

    /// Make monomial.
    //  Construct the monomial with exponents [b]m[/b] and coefficient [b]t[/b].
    static FTps makeMonomial(const FMonomial<N> &m, const T &t);

    /// Unary plus.
    FTps operator+() const;

    /// Unary minus.
    FTps operator-() const;

    /// Add and assign.
    FTps &operator+=(const FTps &y);

    /// Subtract and assign.
    FTps &operator-=(const FTps &y);

    /// Multiply and assign.
    FTps &operator*=(const FTps &y);

    /// Divide and assign.
    FTps &operator/=(const FTps &y);

    /// Add constant and assign.
    FTps &operator+=(const T &y);

    /// Subtract constant and assign.
    FTps &operator-=(const T &y);

    /// Multiply by constant and assign.
    FTps &operator*=(const T &y);

    /// Divide by constant and assign.
    FTps &operator/=(const T &y);

    /// Scale monomial coefficients by coefficients in [b]y[/b].
    FTps scaleMonomials(const FTps &y) const;

    /// Multiply by variable [b]var[/b].
    FTps multiplyVariable(int var, int trunc = EXACT) const;

    /// Multiplication.
    FTps multiply(const FTps &y, int trunc = EXACT) const;

    /// Reciprocal, 1/(*this).
    FTps inverse(int trunc = EXACT) const;

    /// Division.
    //  Throws exception if handed a pair of EXACT Tps's with trunc > globalTruncOrder.
    FTps divide(const FTps &y, int trunc = EXACT) const;

    /// Equality operator.
    //  Returns "true" if polynomials agree
    //  (up to the smaller truncation order).
    bool operator==(const FTps &y) const;

    /// Equality with constant.
    bool operator==(const T &y) const;

    /// Inequality operator.
    bool operator!=(const FTps &y) const;

    /// Inequality with constant.
    bool operator!=(const T &y) const;

    /// Evaluate monomials at point.
    // NB: this method uses static local memory.
    static Array1D<T> evalMonoms(const FVector<T, N> &, int);

    /// Evaluate FTps at point.
    // NB: this method uses static local memory.
    T evaluate(const FVector<T, N> &) const;

    /// Return orders {min, max, trc} of f(rhs(z)).
    // NB: This routine is NOT guaranteed to return min <= max <= trc.
    // If min exceeds max, then f(rhs(z)) = 0 + O(z^{trunc+1}).
    // This possibility must be checked for.  Also, if both TPSs are
    // EXACT and of sufficiently low order, then the result can be EXACT
    // even if trunc is finite.
    Array1D<int> getSubstOrders(const FVps<T, N> &rhs, int trunc = EXACT) const;

    /// Return orders {min, max, trc} of f(rhs(z)),
    //  with ordersL = orders(f), and ordersR = orders(rhs).
    // NB: This routine is NOT guaranteed to return min <= max <= trc.
    // If min exceeds max, then f(rhs(z)) = 0 + O(z^{trunc+1}).
    // This possibility must be checked for.  Also, if both TPSs are
    // EXACT and of sufficiently low order, then the result can be EXACT
    // even if trunc is finite.
    static Array1D<int> getSubstOrders(Array1D<int> &ordersL, Array1D<int> &ordersR,
                                       int trunc = EXACT);

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b] to transform the
    //  order n part of an FTps; leave remaining parts unchanged.  Return a new FTps.
    // NB: This method uses static local memory.  Note that instead of static local
    //     memory, one could use the data storage of an FTps---it's the correct size,
    //     it's reference counted, and it's retained on a stack.
    FTps substitute(const FMatrix<T, N, N> &M, int n) const;

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b] to transform the
    //  order nl through nh parts of an FTps; leave remaining parts unchanged.
    //  Return a new FTps.
    // NB: This method uses static local memory.
    FTps substitute(const FMatrix<T, N, N> &M, int nl, int nh) const;

    /// Substitute.
    //  Use the linear map with matrix representation [b]M[/b]
    //  to transform an FTps.  Return a new FTps.
    // NB: This method (indirectly) uses static local memory.
    FTps substitute(const FMatrix<T, N, N> &M) const;

    /// Substitute.
    //  Use the (nonlinear) map [b]m[/b] to transform an FTps.  Return a new FTps.
    // NB: this method uses static local memory.
    FTps substitute(const FVps<T, N> &m, int trunc = EXACT) const;

    /// Partial derivative.
    //  Return partial derivative with respect to variable [b]var[/b].
    //  Return zero for a constant.
    FTps derivative(int var) const;

    /// Gradient.
    FVps<T, N> gradient() const;

    /// Partial integral.
    //  Return partial integral with respect to variable [b]var[/b].
    FTps integral(int var, int trunc = EXACT) const;

    /// Taylor series.
    //  Return a Taylor series with given order and coefficients.
    //  f.taylor(c,n) ==> c[0] + c[1] * (f-f[0]) + ... + c[n] * (f-f[0])^n
    FTps taylor(const Array1D<T> &series, int order) const;

    /// Make representation unique.
    inline void unique();
    
    /// Get a list containing the indexes of non-zero coefficients of a FTps
    // Returns a STL list containing the indexes
    std::list<int> getListOfNonzeroCoefficients() const;
    
    /// Extract exponents of coefficient
    // Retuns a 1D Array containing the exponents to index [b]index[/b].
    FArray1D<int, N> extractExponents(int index) const;
    
    /// Multiply FTps with itself
    // Return the power of the truncated power series
    FTps<T, N> makePower(int power) const;

    /// Read FTps on the stream [b]is[/b].
    std::istream &get(std::istream &is);

    /// Write FTps on the stream [b]os[/b].
    std::ostream &put(std::ostream &os) const;

    /// Representation of infinite precision.
    static const int EXACT;// = INT_MAX;

private:

    // Memory allocator.  Return an unfilled representation.
    static FTpsRep<T, N> *allocate(int minOrder, int maxOrder, int trcOrder);

    // Memory deallocator.  Return the representation to the memory pool.
    static void deallocate(FTpsRep<T, N> *);

    // Grow the representation, if necessary, so as to accomodate a higher order.
    void grow(int maxOrder, int trcOrder);

    // Report orders of underlying FTpsRep.  (Useful for debugging.)
    Array1D<int> getRepOrders() const;

    // Check that min-, max-, and trcOrder's have the correct relationships.
    // If not, complain with a message that names the given "method".
    static void checkOrders(const std::string &method, int minOrder, int maxOrder, int &trcOrder);

    // Pointer to representation.
    // This is the only non-static data member of class FTps<T,N>.
    FTpsRep<T, N> *itsRep;

    // The free storage list: freeList[ord] points to the first entry
    // in a linked list of FTpsRep's allocated for order [b]ord[/b].
    static FTpsRep<T, N> *freeList[100];

    // The global truncation order.
    static int globalTruncOrder;
};


// Global functions.
// ------------------------------------------------------------------------

/// Add.
template <class T, int N>
FTps<T, N> operator+(const FTps<T, N> &, const FTps<T, N> &);

/// Subtract.
template <class T, int N>
FTps<T, N> operator-(const FTps<T, N> &, const FTps<T, N> &);

/// Add.
template <class T, int N>
FTps<T, N> operator+(const FTps<T, N> &, const T &);

/// Subtract.
template <class T, int N>
FTps<T, N> operator-(const FTps<T, N> &, const T &);

/// Add.
template <class T, int N>
FTps<T, N> operator+(const T &, const FTps<T, N> &);

/// Subtract.
template <class T, int N>
FTps<T, N> operator-(const T &, const FTps<T, N> &);

/// Multiply.
template <class T, int N>
FTps<T, N> operator*(const FTps<T, N> &, const FTps<T, N> &);

/// Divide.
template <class T, int N>
FTps<T, N> operator/(const FTps<T, N> &, const FTps<T, N> &);

/// Multiply.
template <class T, int N>
FTps<T, N> operator*(const FTps<T, N> &, const T &);

/// Divide.
template <class T, int N>
FTps<T, N> operator/(const FTps<T, N> &, const T &);

/// Multiply.
template <class T, int N>
FTps<T, N> operator*(const T &, const FTps<T, N> &);

/// Divide.
template <class T, int N>
FTps<T, N> operator/(const T &, const FTps<T, N> &);

/// Equality.
template <class T, int N>
bool operator==(const T &, const FTps<T, N> &);

/// Inequality.
template <class T, int N>
bool operator!=(const T &, const FTps<T, N> &);

/// Build the exponential series.
//  Return the series exp(:H:) z,
//  the Lie transform exp(:H:) acting on the identity map.
template <class T, int N>
FVps<T, N> ExpMap(const FTps<T, N> &H, int trunc = FTps<T, N>::EXACT);

/// Build the exponential series.
//  Return the series exp(:H:) f,
//  the Lie transform exp(:H:) acting on the function [b]f[/b].
template <class T, int N>
FTps<T, N> ExpMap(const FTps<T, N> &H, const FTps<T, N> &f, int trunc = FTps<T, N>::EXACT);

/// Poisson bracket.
template <class T, int N>
FTps<T, N> PoissonBracket(const FTps<T, N> &f, const FTps<T, N> &g, int trunc = FTps<T, N>::EXACT);

/// Extract FTps from stream [b]is[/b].
template <class T, int N>
std::istream &operator>>(std::istream &is, FTps<T, N> &);

/// Insert FTps into stream [b]os[/b].
template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FTps<T, N> &);


// Implementation.
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTps.cpp"

#endif // CLASSIC_FTps_HH
