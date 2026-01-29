#ifndef CLASSIC_FTpsData_H
#define CLASSIC_FTpsData_H

// ------------------------------------------------------------------------
// $RCSfile: FTpsData.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.10 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FTpsData
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:54 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "Algebra/TpsSubstitution.h"
#include "FixedAlgebra/FArray1D.h"
#include "FixedAlgebra/FMonomial.h"

// Template class FTpsData<N>
// ------------------------------------------------------------------------
/// Internal utility class for FTps<T,N> class.

template <int N>
class FTpsData {

public:

    FTpsData();
    ~FTpsData();

    // Return exponents for the monomial with given Giorgilli index.
    static inline const FMonomial<N> &getExponents(int index);

    // Return the Giorgilli index for a given monomial.
    static int getIndex(const FMonomial<N> &);

    // Order of monomial which has Giorgilli index "index".
    static inline int getOrder(int index);

    // Number of terms in an FTps<T,N> object of order "order".
    static inline int getSize(int order);

    // Index at which order "order" begins; this equals the number
    // of monomials of degree <= (order - 1) in an FTps<T,N> object.
    static inline int orderStart(int order);

    // Index at which order "order" begins for polynomials in nv (<= N) variables.
    static inline int orderStart(int order, int nv);

    // One plus the index at which order "order" ends.
    static inline int orderEnd(int order);

    // One plus index at which order "order" ends for polynomials in nv (<= N) variables.
    static inline int orderEnd(int order, int nv);

    // Number of monomials of degree "order".
    static inline int orderLength(int order);

    // Number of monomials such that "orderL" <= degree <= "orderH".
    static inline int orderLength(int orderL, int orderH);

    // Return the product index array for monomial "index".
    static inline const Array1D<int> &getProductArray(int index);

    // Return the variable list for monomial "index".
    static inline const Array1D<int> &getVariableList(int index);

    // Return the substitution table.
    static inline const Array1D<TpsSubstitution> &getSubTable();

    // To ensure proper table setup, this method must be
    // called by FTps whenever a new object is created.
    static void setup(int order);

private:

    // Not implemented.
    FTpsData(const FTpsData &);
    void operator=(const FTpsData &);

    // Build tables for given order.
    void build(int order);

    // Initialise Substitution table.
    void fillSubst(int var, int order, FMonomial<N> &pow, int &next);

    // Binomial coefficients: binom[i-1][k] = (i+k)!/(i! k!); used in the indexing
    // methods getIndex(), getSize(), orderStart(), orderEnd, and orderLength().
    Array1D < FArray1D < int, N + 1 > > binom;

    // Boundary indices for each order.
    Array1D<int> bin;

    // Exponents: expon[index] = set of exponents in monomial "index".
    Array1D < FMonomial<N> > expon;

    // Indexing information, used for Tps multiplication.
    // prod[i][j] = index of monomial that results from
    //              multiplying monomials i and j
    Array1D < Array1D<int> > prod;

    // Indexing information, used for Tps multiplication
    // using the look-back method.  These two arrays are
    // used to identify which monomial pairs contribute
    // to a given monomial in the result.
    Array1D <int> lookBack;
    Array1D <int> lbBound;

    // Indexing information, used for Tps multiplication
    // using the compact index method.
    Array1D<int> Giorgilli2ExponIndex;

    // Array of variable lists; used for Tps substitution.
    Array1D < Array1D<int> > vrblList;

    // The substitution table.
    Array1D < TpsSubstitution > subTable;

    // The current maximum order for which the tables exist.
    static int topOrder;

    // The current maximum size of any FTps object.
    static int topSize;

    // The singleton object for local book-keeping.
    static FTpsData<N> *theBook;
};


// Class FTpsData, public methods.
// ------------------------------------------------------------------------

// The singleton object for local book-keeping.
template <int N>
FTpsData<N> *FTpsData<N>::
theBook = new FTpsData<N>();

// The current maximum order for which the tables exist.
template <int N>
int FTpsData<N>::topOrder;

// The current maximum size of any FTps object.
template <int N>
int FTpsData<N>::topSize;


template <int N> inline
const FMonomial<N> &FTpsData<N>::
getExponents(int index) {
    return theBook->expon[index];
}


template <int N>
int FTpsData<N>::
getIndex(const FMonomial<N> &pows) {
    int order = 0;
    int index = 0;

    for(int var = N; var-- > 0;) {
        order += pows[var];
        index += theBook->binom[order][var];
    }

    return index;
}


template <int N> inline
int FTpsData<N>::
getOrder(int index) {
    return theBook->expon[index].getOrder();
}


template <int N> inline
int FTpsData<N>::
getSize(int order) {
    return theBook->bin[order+1];
}


template <int N> inline
int FTpsData<N>::
orderStart(int order) {
    return theBook->bin[order];
}


template <int N> inline
int FTpsData<N>::
orderStart(int order, int nv) {
    return theBook->binom[order][N-nv];
}


template <int N> inline
int FTpsData<N>::
orderEnd(int order) {
    return theBook->bin[order+1];
}


template <int N> inline
int FTpsData<N>::
orderEnd(int order, int nv) {
    return theBook->binom[order+1][N-nv];
}


template <int N> inline
int FTpsData<N>::
orderLength(int order) {
    //  return theBook->bin[order+1] - theBook->bin[order];
    return theBook->binom[order+1][1];
}


template <int N> inline
int FTpsData<N>::
orderLength(int orderL, int orderH) {
    return theBook->bin[orderH+1] - theBook->bin[orderL];
}


template <int N> inline
const Array1D<int> &FTpsData<N>::
getProductArray(int index) {
    return theBook->prod[index];
}


template <int N> inline
const Array1D<int> &FTpsData<N>::
getVariableList(int index) {
    return theBook->vrblList[index];
}


//template <int N> inline
//const Array1D<int> &FTpsData<N>::
//getVariableList(int index)
//{
//  FMonomial<N> jlist = theBook->expon[index];
//  Array1D<int> *result = new Array1D<int> (jlist.getOrder(),0);
//
//  int k = 0;
//  for (int v = 0; v < N; ++v) {
//    int jv = jlist[v];
//    while (jv-- > 0) (*result)[k++] = v;
//  }
////std::cerr << " vrblList(" << index << "[" << jlist << "]" << ") = " << *result;
//  return *result;
//}


template <int N> inline
const Array1D<TpsSubstitution> &FTpsData<N>::
getSubTable() {
    return theBook->subTable;
}


template <int N>
void FTpsData<N>::
setup(int order) {
    theBook->build(order);
}


// Class FTpsData, protected methods.
// ------------------------------------------------------------------------

template <int N>
FTpsData<N>::
FTpsData() {
    // Initialize with order 4.  May expand later.
    topOrder = 0;
    build(4);
}


template <int N>
FTpsData<N>::
~FTpsData()
{}


// Class FTpsData, private methods.
// ------------------------------------------------------------------------

template <int N>
void FTpsData<N>::
build(int order) {
    if(topOrder < order) {
        topOrder = order;
        int prodSize = 0;

        // Build array containing binomial coefficients:
        // binom[0 ... topOrder+1][0 ... N].
        binom = Array1D < FArray1D < int, N + 1 > > (topOrder + 2);
        // Here var runs from N to 0.
        for(int var = N + 1; var-- > 0;) binom[0][var] = 0;
        for(int ord = 1; ord <= topOrder + 1; ord++) {
            binom[ord][N] = 1;
            // Here var runs from N-1 to 0.
            for(int var = N; var-- > 0;)
                binom[ord][var] = binom[ord-1][var] + binom[ord][var+1];
        }

        // Build array of indices that mark order boundaries.
        bin = Array1D<int>(topOrder + 2);
        for(int i = 0; i < topOrder + 2; i++) bin[i] = binom[i][0];

        // Build table of monomial exponents for orders zero through topOrder.
        // (algorithm due to Liam Healey (Marylie 3.0)).
        {
            topSize = bin[topOrder+1];
            expon = Array1D < FMonomial<N> >(topSize);
            FMonomial<N> power;

            if(N == 1)
                for(int index = 1; index < topSize; index++) expon[index][0] = index;
            else
                for(int index = 1; index < topSize; index++) {
                    int carry = power[N-1];
                    power[N-1] = 0;
                    int lastnz = N - 2;
                    while(power[lastnz] == 0 && lastnz-- > 0) {}

                    if(lastnz == -1) power[0] = 1 + carry;
                    else {
                        power[lastnz]--;
                        power[lastnz+1] += 1 + carry;
                    }

                    expon[index] = power;
                }
        }

        // Build table of products for orders zero through topOrder.
        {
            FMonomial<N> power;
            prod = Array1D < Array1D<int> >(topSize);
            for(int xord = 0; xord <= topOrder; xord++) {
                int yord = topOrder - xord;
                int ysize = bin[yord+1];
                for(int i = bin[xord]; i < bin[xord+1]; i++) {
                    prod[i] = Array1D<int>(ysize);
                    // use symmetry for LL half of prod array
                    for(int j = 0; j < std::min(i, ysize); j++) {
                        prod[i][j] = prod[j][i];
                        ++prodSize;
                    }
                    for(int j = i; j < ysize; j++) {
                        power = expon[i] * expon[j];
                        int ord = 0;
                        int ind = 0;
                        for(int vv = N; vv-- > 0;) {
                            ord += power[vv];
                            ind += binom[ord][vv];
                        }
                        prod[i][j] = ind;
                        ++prodSize;
                    }
                }
            }
        }

        // Build table of variable lists.
        {
            // Allocate variable list pointers and working array.
            vrblList = Array1D< Array1D<int> >(topSize);
            int *vars = new int[topOrder];
            // Initialise counter.
            int j = 1, N1 = N - 1;
            // Loop over the orders.
            for(int ord = 1; ord <= topOrder; ord++) {
                // Allocate first variable list; initialise both it and vars.
                vrblList[j] = Array1D<int>(ord, 0);
                std::fill(vars, vars + ord, int(0));
                // Define end points.
                int last_j = bin[ord+1];
                int *vlast = vars + ord;
                // Build remaining variable lists at this order.
                while(++j < last_j) {
                    // From last variable list (in vars), construct next one.
                    int *vi = vlast;
                    while(*--vi == N1) {};
                    int k = *vi + 1;
                    std::fill(vi, vlast, k);
                    // Allocate j-th variable list and copy from vars.
                    vrblList[j] = Array1D<int>(ord);
                    std::copy(vars, vars + ord, vrblList[j].begin());
                }
            }
            delete [] vars;
        }

        // Build the substitution table.
        {
            subTable = Array1D<TpsSubstitution>(topSize);
            FMonomial<N> power;
            int next = 1;
            fillSubst(0, 1, power, next);
        }
    }
}


template <int N>
void FTpsData<N>::
fillSubst(int var, int order, FMonomial<N> &power, int &next) {
    for(int v = var; v < N; v++) {
        TpsSubstitution &s = subTable[next];
        power[v]++;
        int ord = 0;
        int ind = 0;
        for(int vv = N; vv-- > 0;) {
            ord += power[vv];
            ind += binom[ord][vv];
        }
        s.index = ind;
        s.order = order;
        s.variable = v;
        next++;
        if(order < topOrder) fillSubst(v, order + 1, power, next);
        s.skip = next;
        power[v]--;
    }
}

#endif // CLASSIC_FTpsData_H
