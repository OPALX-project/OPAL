#ifndef CLASSIC_DragtFinnMap_HH
#define CLASSIC_DragtFinnMap_HH

// ------------------------------------------------------------------------
// $RCSfile: DragtFinnMap.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: DragtFinnMap<N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2002/03/25 20:44:16 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FTps.h"

#include <complex>

template <class T, int N> class FLieGenerator;
template <class T, int N> class FVector;
template <class T, int N> class FVps;
template <class T, int N> class LinearMap;


// Template class DragtFinnMap<N>.
// ------------------------------------------------------------------------
/// A Lie algebraic map, factored according to Dragt and Finn.
//  The type of coefficients is double, and the number of variables is N.
//  It should be noted that some algorithms, notably concatenation of maps,
//  are limited to maps with generators up to order 6.

template <int N>
class DragtFinnMap {

public:

    DragtFinnMap(const DragtFinnMap &);
    ~DragtFinnMap();
    const DragtFinnMap &operator=(const DragtFinnMap &);

    /// Construct identity map.
    DragtFinnMap();

    /// Construct identity map with space for given order.
    explicit DragtFinnMap(int);

    /// Build factored representation from a Taylor series map.
    explicit DragtFinnMap(const FVps<double, 2 * N> &map);

    /// Build factored representation from a pseudo-Hamiltonian.
    explicit DragtFinnMap(const FTps<double, 2 * N> &map);

    /// Build factored representation from a linear map.
    explicit DragtFinnMap(const LinearMap<double, 2 * N> &map);

    /// Convert to Taylor series map.
    operator FVps<double, 2 * N>() const;

    /// Convert to Taylor series map.
    operator LinearMap<double, 2 * N>() const;


    /// Assign the matrix representing the linear transform.
    void assign(const FMatrix<double, 2 * N, 2 * N> &);

    /// Assign the complete set of generators.
    void assign(const FTps<double, 2 * N> &);

    /// Assign the generator for a selected order.
    void assign(const FLieGenerator<double, N> &);


    // Factorise a map given by a Hamiltonian,
    // using the method by M. Berz, E. Forest and J. Irwin
    // (Particle Accelerators, vol. 24, pp. 91 / 107).
    // This method works best for Hamiltonians whose H_2 have no
    // zero eigenvalues.
    static DragtFinnMap factorBerzForestIrwin(const FTps<double, 2 * N> &H);

    // Factorise a map given by a Hamiltonian,
    // using the method by D. Douglas  (ML University thesis 1982).
    // This method works best for Hamiltonians with q-independent H_2.
    static DragtFinnMap factorDouglas(const FTps<double, 2 * N> &H);

    // Factorise a map given by a Hamiltonian,
    // which is either p-free or q-free.
    // The Lie generators are then simply the homogeneous parts of H.
    static DragtFinnMap factorSimple(const FTps<double, 2 * N> &H);


    /// Return the matrix representing the linear transform.
    const FMatrix<double, 2 * N, 2 * N> &getMatrix() const;

    /// Return the matrix representing the linear transform.
    FMatrix<double, 2 * N, 2 * N> &getMatrix();

    /// Return the complete set of generators.
    const FTps<double, 2 * N> &getGenerators() const;

    /// Return the complete set of generators.
    FTps<double, 2 * N> &getGenerators();

    /// Return the generator for a selected order.
    const FLieGenerator<double, N> getGenerator(int) const;

    /// Return the order of the generators.
    int getOrder() const;

    /// Compute inverse map.
    //  The first-order generator must be zero.
    DragtFinnMap inverse() const;

    /// Compute inverse factorisation of map.
    //  The first-order generator must be zero.
    DragtFinnMap reverse() const;

    /// Add to set of generators.
    void operator+=(const FTps<double, 2 * N> &);

    /// Subtract from set of generators.
    void operator-=(const FTps<double, 2 * N> &);

    /// Add to generator of given order.
    void operator+=(const FLieGenerator<double, N> &);

    /// Subtract from generator of given order.
    void operator-=(const FLieGenerator<double, N> &);


    /// Substitute (concatenate) with another map in beam order.
    DragtFinnMap catenate(const DragtFinnMap &) const;

    /// Conjugate with another map.
    DragtFinnMap conjugate(const DragtFinnMap &) const;

    /// Compute dynamic fixed point.
    //  First argument is set to the fixed point vector,
    //  Second argument is set to map around the fixed point.
    void dynamicFixedPoint(FVector<double, 2 * N> &fp, DragtFinnMap &map);

    /// Compute static fixed point.
    //  First argument is set to the fixed point vector,
    //  Second argument is set to map around the fixed point.
    //  Assumes the longitudinal values are in the last pair.
    void staticFixedPoint(FVector<double, 2 * N> &fp, DragtFinnMap &map);

    /// Split map into dispersion map and non-dispersive part.
    //  First argument is set to dispersion map,
    //  Second argument is set to map without dispersion.
    void removeDispersion(DragtFinnMap &dm, DragtFinnMap &map);

    /// Track an orbit through the map.
    //  The first argument is the input/output orbit,
    //  the second argument is set to the map around the orbit.
    void trackOrbit(FVector<double, 2 * N> &orbit, DragtFinnMap &map);

    /// Build the terminating exponential series.
    //  Return the series  exp(:f:) g, the Lie transform exp(:f:)
    //  acting on the generators g of this map.
    DragtFinnMap transform(const FLieGenerator<double, N> &g, int topOrder);

private:

    // Get beginning of coefficient storage for generator of given order.
    //  Version for constant object.
    inline const double *begin(int) const;

    // Get beginning of coefficient storage for generator of given order.
    // Version for non-constant object.
    inline double *begin(int);

    // Get ending of coefficient storage for generator of given order.
    // Version for constant object.
    inline const double *end(int) const;

    // Get ending of coefficient storage for generator of given order.
    // Version for non-constant object.
    inline double *end(int);


    // Substitute (concatenate) with another map, ignoring linear terms.
    DragtFinnMap
    catenateZero(const DragtFinnMap &) const;

    // Exponentiate a matrix whose norm is small.
    static FMatrix<double, 2 * N, 2 * N>
    exponentiate(const FMatrix<double, 2 * N, 2 * N> &);

    // Factorize a map given by a Hamiltonian ignoring the linear terms.
    static DragtFinnMap
    factorize(const FTps<double, 2 * N> &H);

    // Make the matrix J*S representing the Lie operator :f_2:.
    static FMatrix<double, 2 * N, 2 * N>
    makeMatrix(const FLieGenerator<double, N> &f_2);

    // Move the g_1 term over the terms of f.
    static void
    move_g_1(DragtFinnMap &f, DragtFinnMap &g);

    // Order the modes according to "stable", "unstable", and "coasting".
    // Return the number of non-coasting modes.
    // The eigenvalues in this case are the logarithms of the usual values.
    static int
    orderModes(FMatrix<double, 2 * N, 2 * N> &, FVector<std::complex<double>, 2 * N> &);


    // The matrix of the repesentation.
    FMatrix<double, 2 * N, 2 * N> itsMatrix;

    // The Lie generators are stored in an FTps<double,2*N> object.
    FTps<double, 2 * N> itsGenerators;
};


template <int N>
std::ostream &operator<<(std::ostream &, const DragtFinnMap<N> &);

#include "FixedAlgebra/FDoubleEigen.h"
#include "FixedAlgebra/FLieGenerator.h"
#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsData.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/LinearMap.h"
#include "FixedAlgebra/Taylor.h"
#include "Utilities/LogicalError.h"
#include <complex>
#include <iosfwd>
//#include <algorithm>
//#include <iostream>


// Template class DragtFinnMap.
// ------------------------------------------------------------------------

template <int N>
DragtFinnMap<N>::
DragtFinnMap(const DragtFinnMap &rhs):
    itsMatrix(rhs.itsMatrix), itsGenerators(rhs.itsGenerators)
{}


template <int N>
DragtFinnMap<N>::
DragtFinnMap():
    itsMatrix(), itsGenerators() {
    for(int i = 0; i < 2 * N; ++i) itsMatrix(i, i) = 1.0;
}


template <int N>
DragtFinnMap<N>::
DragtFinnMap(int order):
    itsMatrix(), itsGenerators(0, order, order) {
    for(int i = 0; i < 2 * N; ++i) itsMatrix(i, i) = 1.0;
}


template <int N>
DragtFinnMap<N>::
~DragtFinnMap()
{}


template <int N>
const DragtFinnMap<N> &DragtFinnMap<N>::
operator=(const DragtFinnMap &rhs) {
    itsMatrix = rhs.itsMatrix;
    itsGenerators = rhs.itsGenerators;
    return *this;
}


template <int N>
DragtFinnMap<N>::
DragtFinnMap(const FVps<double, 2 * N> &rhs):
    itsMatrix(rhs.linearTerms()), itsGenerators() {
    for(int i = 0; i < N; ++i) {
        itsGenerators[2*i+1] =   rhs[2*i+1][0];
        itsGenerators[2*i+2] = - rhs[2*i+0][0];
    }
    // ***** MISSING: higher order terms *****
}


template <int N>
DragtFinnMap<N>::
DragtFinnMap(const LinearMap<double, 2 * N> &rhs):
    itsMatrix(rhs.linearTerms()), itsGenerators() {
    for(int i = 0; i < N; ++i) {
        itsGenerators[2*i+1] =   rhs[2*i+1][0];
        itsGenerators[2*i+2] = - rhs[2*i+0][0];
    }
}


template <int N>
DragtFinnMap<N>::
operator FVps<double, 2 * N>() const {
    int maxOrder = getOrder();
    FTps<double, 2 * N> gen = getGenerators();
    FVps<double, 2 * N> theMap;

    for(int i = maxOrder; i >= 3; --i) {
        FTps<double, 2 * N> tps = gen.filter(i, i);
        theMap = ExpMap(tps, theMap);
    }

    theMap = theMap.substitute(getMatrix());

    FVps<double, 2 * N> xlate;
    for(int i = 0; i < N; ++i) {
        xlate[2*i][0] = -gen[2*i+2];
        xlate[2*i+1][0] = gen[2*i+1];
    }
    theMap = theMap.substitute(xlate);

    return theMap.filter(0, maxOrder < 2 ? 1 : maxOrder - 1);
}


template <int N>
DragtFinnMap<N>::
operator LinearMap<double, 2 * N>() const {
    FLieGenerator<double, N> f_1(getGenerators(), 1);
    f_1 = f_1.transform(getMatrix());
    LinearMap<double, 2 * N> theMap(getMatrix());

    for(int i = 0; i < N; ++i) {
        theMap[2*i+0][0] = - f_1[2*i+2];
        theMap[2*i+1][0] =   f_1[2*i+1];
    }

    return theMap;
}


template <int N>
DragtFinnMap<N>::
DragtFinnMap(const FTps<double, 2 * N> &H) {
    if(H.filter(1, 1) == 0.0) {
        // No linear terms in Hamiltonian; no displacement map required.
        *this = factorize(H);
    } else {
        // Apply the method by L. Healy to remove linear terms.
        // Factorise terms without first order.
        int order = H.getMaxOrder();

        // The h map contains the first-order terms.
        DragtFinnMap<N> h;
        h.assign(-H.filter(1, 1));

        // The g map is the current factor.
        DragtFinnMap<N> g = factorize(H.filter(2, order));

        // The "this" map accumulates the result.
        *this = g;

        // Extract and transform first order terms.
        // The epsilon rank of the remainder terms.
        for(int epsilonRank = 1; epsilonRank <= order - 3; ++epsilonRank) {
            int top = order - epsilonRank;

            // Put first-order terms into h and transform with g map.
            h.assign(h.getGenerators().filter(1, 1).substitute(g.getMatrix()));
            for(int ord = 3; ord <= top; ++ord) {
                h = h.transform(g.getGenerator(ord), top);
            }

            // Compute new g map and catenate with f.
            g = factorize(h.getGenerators().filter(2, top));
            *this = g.catenateZero(*this);
        }

        // Transmit the first-order terms.
        assign(h.getGenerator(1));
    }
}


template <int N>
void DragtFinnMap<N>::
assign(const FMatrix<double, 2 * N, 2 * N> &mat) {
    itsMatrix = mat;
}


template <int N>
void DragtFinnMap<N>::
assign(const FTps<double, 2 * N> &gen) {
    itsGenerators = gen;
}


template <int N>
void DragtFinnMap<N>::
assign(const FLieGenerator<double, N> &gen) {
    int order = gen.getOrder();

    if(itsGenerators.getMaxOrder() < order) {
        FTps<double, 2 * N> temp(0, order, FTps<double, 2 * N>::getGlobalTruncOrder());
        std::copy(itsGenerators.begin(), itsGenerators.end(), temp.begin());
        itsGenerators = temp;
    } else {
        itsGenerators.unique();
    }

    std::copy(gen.begin(), gen.end(), begin(order));
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
factorBerzForestIrwin(const FTps<double, 2 * N> &HH) {
    static const double tol = 1.0e-12;
    int modes = N;

    // The matrix M represents the operator :f_2:.
    FTps<double, 2 * N> H = - HH;
    FLieGenerator<double, N> H_2(H, 2);
    FMatrix<double, 2 * N, 2 * N> M = makeMatrix(H_2);

    // Eigenvalues and eigenvectors of M.
    FVector<std::complex<double>, 2 * N> mu; // eigenvalues
    FMatrix<double, 2 * N, 2 * N> V;   // eigenvectors
    FMatrix<double, 2 * N, 2 * N> V_inv; // inverse of V
    {
        // Find the eigenvalues and eigenvectors of M.
        FDoubleEigen<2 * N> eigen(M, true);
        mu = eigen.eigenValues();
        V  = eigen.packedEigenVectors();

        // Order eigenvectors and normalise.
        modes = orderModes(V, mu);
        FLUMatrix<double, 2 * N> lu(V);
        V_inv = lu.inverse();
    }

    // Set up the auxiliary matrices.
    FMatrix<double, 2 * N, 2 * N> Rot, R_dir, R_inv, I_dir;

    for(int i = 0; i < 2 * N; i += 2) {
        if(std::abs(mu[i]) < tol) {
            // "coasting" eigenvalue pair.
            R_dir(i,  i)   = 1.0;
            R_dir(i + 1, i + 1) = 1.0;
            R_inv(i,  i)   = 1.0;
            R_inv(i + 1, i + 1) = 1.0;
            I_dir(i,  i)   = 1.0;
            I_dir(i + 1, i + 1) = 1.0;
            Rot(i,  i)     = 1.0;
            Rot(i,  i + 1)   = M(i, i + 1);
            Rot(i + 1, i)     = M(i + 1, i);
            Rot(i + 1, i + 1)   = 1.0;
        } else {
            R_dir(i,  i)   = + 0.5;
            R_dir(i,  i + 1) = + 0.5;
            R_dir(i + 1, i)   = + 0.5;
            R_dir(i + 1, i + 1) = - 0.5;
            R_inv(i,  i)   = + 1.0;
            R_inv(i,  i + 1) = + 1.0;
            R_inv(i + 1, i)   = + 1.0;
            R_inv(i + 1, i + 1) = - 1.0;

            if(std::abs(std::imag(mu[i])) > tol) {
                // "stable" eigenvalue pair.
                double c       = std::cos(std::imag(mu[i]));
                double s       = std::sin(std::imag(mu[i]));
                Rot(i,  i)     = + c;
                Rot(i,  i + 1)   = + s;
                Rot(i + 1, i)     = - s;
                Rot(i + 1, i + 1)   = + c;
                I_dir(i,  i + 1) = 1.0;
                I_dir(i + 1, i)   = 1.0;
            } else {
                // "unstable" eigenvalue pair.
                double ch      = std::cosh(std::real(mu[i]));
                double sh      = std::sinh(std::real(mu[i]));
                Rot(i,  i)     = + ch;
                Rot(i,  i + 1)   = - sh;
                Rot(i + 1, i)     = - sh;
                Rot(i + 1, i + 1)   = + ch;
                I_dir(i,  i)   = 1.0;
                I_dir(i + 1, i + 1) = 1.0;
            }
        }
    }

    // Declare the fake normal form and normalising maps.
    int maxOrder = H.getMaxOrder();
    DragtFinnMap<N> N_scr;
    DragtFinnMap<N> A_scr;
    A_scr.assign(V_inv);
    N_scr.assign(Rot);
    N_scr.assign(H.filter(2, maxOrder).substitute(V));

    // Remove non-resonant terms order by order.
    for(int omega = 3; omega <= maxOrder; omega++) {
        // Compute the terms to be removed and store in f.
        FLieGenerator<double, N> f(N_scr.getGenerators(), omega);
        FLieGenerator<double, N> a(omega);
        FLieGenerator<double, N> b(omega);
        FLieGenerator<double, N> pi(omega);

        for(int m = f.getBottomIndex(); m < f.getTopIndex(); m++) {
            const FMonomial<2 * N> &index = FTpsData<2 * N>::getExponents(m);
            std::complex<double> factor = 0.0;
            int count = 0;

            for(int j = 0; j < 2 * modes; j += 2) {
                if(std::abs(std::imag(mu[j])) > tol) count += index[j+1];
                factor += double(index[j]) * mu[j] + double(index[j+1]) * mu[j+1];
            }

            if(std::abs(factor) > tol) {
                // Term can be removed.
                factor = 1.0 / factor;
                a[m] = std::real(factor);
                b[m] = std::imag(factor);
            }

            pi[m] = std::pow(- 1.0, (count + 1) / 2);
        }

        // Compute cal_T^(-1) * f.
        FLieGenerator<double, N> f1 = pi.scale(f).transform(R_dir);
        FLieGenerator<double, N> f2 = f1.transform(I_dir);
        FLieGenerator<double, N> F_omega =
            (a.scale(f1) + b.scale(f2)).transform(R_inv).scale(pi);

        N_scr = N_scr.transform(- F_omega, maxOrder);
        A_scr.assign(F_omega);
    }

    // Conjugate the N-map with the A-map.
    N_scr = N_scr.conjugate(A_scr);
    return N_scr;
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
factorDouglas(const FTps<double, 2 * N> &H) {
    int maxOrder = H.getMaxOrder();
    DragtFinnMap<N> theMap;

    // Build first-order matrix.
    // K2 = - H2 is the second-order generator of the map.
    FLieGenerator<double, N> H2(H, 2);
    FMatrix<double, 2 * N, 2 * N> M = makeMatrix(- H2);
    M = M + 1.0;
    theMap.assign(M);

    // Assign terms of higher order.
    if(maxOrder >= 3) {
        // Build the transformed - H3(exp(:- K2: z).
        Taylor < FLieGenerator<double, N> > H3(3);
        H3[0] = FLieGenerator<double, N>(H, 3);
        H3[1] = PoissonBracket(H2, H3[0]);
        H3[2] = PoissonBracket(H2, H3[1]) / 2.0;
        H3[3] = PoissonBracket(H2, H3[2]) / 3.0;

        // Extract K3.
        Taylor < FLieGenerator<double, N> > K3 = - H3.integrate();
        theMap += K3.sum();

        if(maxOrder >= 4) {
            // Build the transformed H4 = H4(exp(:- K2: z).
            Taylor < FLieGenerator<double, N> > H4(4);
            H4[0] = FLieGenerator<double, N>(H, 4);
            H4[1] = PoissonBracket(H2, H4[0]);
            H4[2] = PoissonBracket(H2, H4[1]) / 2.0;
            H4[3] = PoissonBracket(H2, H4[2]) / 3.0;
            H4[4] = PoissonBracket(H2, H4[3]) / 4.0;

            // Change due to ((exp(:K3:) - 1) / :K3: - 1) H3.
            Taylor < FLieGenerator<double, N> > p33 = PoissonBracket(K3, H3);
            H4 -= p33 / 2.0;

            // Extract K4.
            Taylor < FLieGenerator<double, N> > K4 = - H4.integrate();
            theMap += K4.sum();

            if(maxOrder >= 5) {
                // Build the transformed H5 = H5(exp(:- K2: z).
                Taylor < FLieGenerator<double, N> > H5(5);
                H5[0] = FLieGenerator<double, N>(H, 5);
                H5[1] = PoissonBracket(H2, H5[0]);
                H5[2] = PoissonBracket(H2, H5[1]) / 2.0;
                H5[3] = PoissonBracket(H2, H5[2]) / 3.0;
                H5[4] = PoissonBracket(H2, H5[3]) / 4.0;
                H5[5] = PoissonBracket(H2, H5[4]) / 5.0;

                // Contribution of ((exp(:K3:) - 1) / :K3: - 1) H3
                //   H5 -= PB(K3, PB(K3, H3)) / 6.
                Taylor < FLieGenerator<double, N> > p333 = PoissonBracket(K3, p33);
                H5 -= p333 / 6.0;

                // Contribution of exp(:- K3:) H4
                //   H5 -= PB(K3, H4).
                Taylor < FLieGenerator<double, N> > p34 = PoissonBracket(K3, H4);
                H5 -= p34;

                Taylor < FLieGenerator<double, N> > K_5 = - H5.integrate();
                theMap += K_5.sum();

                if(maxOrder >= 6) {
                    // Build the transformed H6 = H6(exp(:- K2: z).
                    Taylor < FLieGenerator<double, N> > H6(16);
                    H6[0] = FLieGenerator<double, N>(H, 6);
                    H6[1] = PoissonBracket(H2, H6[0]);
                    H6[2] = PoissonBracket(H2, H6[1]) / 2.0;
                    H6[3] = PoissonBracket(H2, H6[2]) / 3.0;
                    H6[4] = PoissonBracket(H2, H6[3]) / 4.0;
                    H6[5] = PoissonBracket(H2, H6[4]) / 5.0;
                    H6[6] = PoissonBracket(H2, H6[5]) / 6.0;

                    // Contibution of ((exp(:K3:) - 1) / :K3: - 1) H_3
                    //   H6 -= PB(K3, PB(K3, PB(K3, H3))).
                    // Contribution of exp(:- K3:) (H4 + H5)
                    //   H6 = H6 + PB(K3, PB(K3, H4)) / 2 - PB(K3, H5)
                    // Note that H5 has already been modified by p34.
                    H6 -= PoissonBracket(K3, p333 + H5 + p34 / 2.0) / 24.0;
                    // Contibution of ((exp(:K4:) - 1) / :K4: - 1) H_4
                    //   H6 -= PB(K4, H4) / 2.
                    H6 -= PoissonBracket(K4, H4) / 2.0;
                    Taylor < FLieGenerator<double, N> > K6 = - H6.integrate();
                    theMap += K6.sum();
                }
            }
        }
    }

    return theMap;
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
factorSimple(const FTps<double, 2 * N> &HH) {
    FTps<double, 2 * N> H = - HH;
    DragtFinnMap<N> theMap;
    theMap.assign(H.filter(3, H.getMaxOrder()));

    // Build first-order matrix.
    const FLieGenerator<double, N> H_2(H, 2);
    const FMatrix<double, 2 * N, 2 * N> M = makeMatrix(H_2);
    theMap.assign(M + 1.0);

    // Assign first-order terms, transformed by inverse matrix.
    const FLieGenerator<double, N> H_1(H, 1);
    if(! H_1.isZero()) theMap.assign(H_1.transform(- M + 1.0));
    return theMap;
}


template <int N>
const FMatrix<double, 2 * N, 2 * N> &DragtFinnMap<N>::
getMatrix() const {
    return itsMatrix;
}


template <int N>
FMatrix<double, 2 * N, 2 * N> &DragtFinnMap<N>::
getMatrix() {
    return itsMatrix;
}


template <int N>
const FTps<double, 2 * N> &DragtFinnMap<N>::
getGenerators() const {
    return itsGenerators;
}


template <int N>
FTps<double, 2 * N> &DragtFinnMap<N>::
getGenerators() {
    return itsGenerators;
}


template <int N>
const FLieGenerator<double, N> DragtFinnMap<N>::
getGenerator(int order) const {
    FLieGenerator<double, N> gen(order);

    if(itsGenerators.getMaxOrder() >= order) {
        std::copy(begin(order), end(order), gen.begin());
    } else {
        std::fill(gen.begin(), gen.end(), 0.0);
    }

    return gen;
}


template <int N>
int DragtFinnMap<N>::
getOrder() const {
    return itsGenerators.getMaxOrder();
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
inverse() const {
    DragtFinnMap h(reverse());
    FMatrix<double, 2 * N, 2 * N> f_mat = getMatrix();
    FLUMatrix<double, 2 * N> lu(f_mat);
    h.assign(lu.inverse());
    h.assign(- h.getGenerators());
    return h;
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
reverse() const {
    int order = getOrder();
    FMatrix<double, 2 * N, 2 * N> f_mat = getMatrix();
    FTps<double, 2 * N> f_gen = getGenerators().substitute(f_mat);
    DragtFinnMap h;
    h.assign(f_mat);
    h.assign(f_gen);

    if(order >= 5) {
        FLieGenerator<double, N> f_3(f_gen, 3);
        FLieGenerator<double, N> f_4(f_gen, 4);
        FLieGenerator<double, N> pb  = PoissonBracket(f_3, f_4);
        h += pb;

        if(order >= 6) {
            FLieGenerator<double, N> f_5(f_gen, 5);
            h += PoissonBracket(f_3, f_5 + 0.5 * pb);
        }
    }

    return h;
}


template <int N>
void DragtFinnMap<N>::
operator+=(const FTps<double, 2 * N> &gen) {
    itsGenerators += gen;
}


template <int N>
void DragtFinnMap<N>::
operator-=(const FTps<double, 2 * N> &gen) {
    itsGenerators -= gen;
}


template <int N>
void DragtFinnMap<N>::
operator+=(const FLieGenerator<double, N> &gen) {
    int order = gen.getOrder();

    if(order > 0) {
        if(itsGenerators.getMaxOrder() < order) {
            FTps<double, 2 * N> temp(0, order, FTps<double, 2 * N>::getGlobalTruncOrder());
            std::copy(itsGenerators.begin(), itsGenerators.end(), temp.begin());
            itsGenerators = temp;
        } else {
            itsGenerators.unique();
        }

        std::transform(begin(order), end(order), gen.begin(),
                       begin(order), std::plus<double>());
    }
}


template <int N>
void DragtFinnMap<N>::
operator-=(const FLieGenerator<double, N> &gen) {
    int order = gen.getOrder();

    if(order > 0) {
        if(itsGenerators.getMaxOrder() < order) {
            FTps<double, 2 * N> temp(0, order, FTps<double, 2 * N>::getGlobalTruncOrder());
            std::copy(itsGenerators.begin(), itsGenerators.end(), temp.begin());
            itsGenerators = temp;
        } else {
            itsGenerators.unique();
        }

        std::transform(begin(order), end(order), gen.begin(),
                       begin(order), std::minus<double>());
    }
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
catenate(const DragtFinnMap &rhs) const {
    if(rhs.itsGenerators.filter(1, 1) == 0.0) {
        // No linear terms in right-hand side; use normal catenation.
        return catenateZero(rhs);
    } else {
        DragtFinnMap f(*this);
        DragtFinnMap g(rhs);
        move_g_1(f, g);
        return f.catenateZero(g);
    }
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
conjugate(const DragtFinnMap &rhs) const {
    return rhs.catenate(catenate(rhs.inverse()));
}


template <int N>
void DragtFinnMap<N>::
dynamicFixedPoint(FVector<double, 2 * N> &fp, DragtFinnMap &map) {
    static const int itmax = 10;
    static const double tol = 1.0e-16;
    fp = FVector<double, 2 * N>();

    for(int iter = 1; iter <= itmax; ++iter) {
        FVector<double, 2 * N> fp1(fp);
        trackOrbit(fp1, map);
        FVector<double, 2 * N> error = fp1 - fp;

        if(euclidean_norm(error) < tol) break;

        FMatrix<double, 2 * N, 2 * N> M = map.getMatrix();
        M = M - 1.0;
        FLUMatrix<double, 2 * N> lu(M);
        lu.backSubstitute(error);
        fp -= error;
    }
}


template <int N>
void DragtFinnMap<N>::
staticFixedPoint(FVector<double, 2 * N> &fp, DragtFinnMap &map) {
    static const int itmax = 10;
    static const double tol = 1.0e-16;
    fp = FVector<double, 2 * N>();

    for(int iter = 1; iter <= itmax; ++iter) {
        FVector<double, 2 * N> fp1(fp);
        trackOrbit(fp1, map);
        FVector<double, 2 * N> error = fp1 - fp;

        if(euclidean_norm(error) < tol) break;

        FMatrix<double, 2 * N, 2 * N> M = map.getMatrix();
        M = M - 1.0;

        FMatrix < double, 2 * N - 2, 2 * N - 2 > MM;
        FVector < double, 2 * N - 2 > X;
        for(int i = 0; i < 2 * N - 2; ++i) {
            for(int j = 0; j < 2 * N - 2; ++j) {
                MM(i, j) = M(i, j);
            }
            X[i] = error[i];
        }

        FLUMatrix < double, 2 * N - 2 > lu(MM);
        lu.backSubstitute(X);
        for(int i = 0; i < 2 * N - 2; ++i) {
            fp[i] -= X[i];
        }
    }
}

template <int N>
void DragtFinnMap<N>::
removeDispersion(DragtFinnMap &dm, DragtFinnMap &map) {
    // compute an1.
    FMatrix<double, 2 * N,  2 * N>   A = getMatrix();
    FMatrix < double, 2 * N - 2, 2 * N - 2 > B;
    FVector < double, 2 * N - 2 >       C;

    // Compute first-order dispersion.
    for(int i = 0; i < 2 * N - 2; ++i) {
        for(int j = 0; j < 2 * N - 2; ++j) {
            B(i, j) = A(i, j);
        }

        B(i, i) -= 1.0;
        C[i] = - A(i, 2 * N - 1);
    }

    FLUMatrix < double, 2 * N - 2 > lu(B);
    lu.backSubstitute(C);
    FMatrix<double, 2 * N, 2 * N> T;
    T = T + 1.0;

    for(int i = 0; i < N - 1; ++i) {
        T(2 * i,  2 * N - 1) =   C[2*i];
        T(2 * i + 1, 2 * N - 1) =   C[2*i+1];
        T(2 * N - 2, 2 * i)   =   C[2*i+1];
        T(2 * N - 2, 2 * i + 1) = - C[2*i];
    }

    dm = DragtFinnMap<N>(getOrder());
    dm.assign(T);

    // Compute t*m*t^(-1) where m is initial map.
    map = conjugate(dm);

    // Compute higher dispersions.
    for(int order = 2; order <= getOrder(); ++order) {
        A = FLUMatrix<double, 2 * N>(map.getMatrix()).inverse();

        for(int i = 0; i < 2 * N - 2; ++i) {
            for(int j = 0; j < 2 * N - 2; ++j) {
                B(i, j) = A(i, j);
            }

            B(i, i) -= 1.0;
            FMonomial<2 * N> mono;
            mono[i] = 1;
            mono[2*N-1] = order;
            C[i] = - map.getGenerators()[mono];
        }

        FLUMatrix < double, 2 * N - 2 > lu(B);
        lu.backSubstitute(C);

        // compute t2.
        DragtFinnMap<N> t1(order + 1);

        for(int i = 0; i < 2 * N - 2; ++i) {
            FMonomial<2 * N> mono;
            mono[i] = 1;
            mono[2*N-1] = order;
            t1.getGenerators()[mono] = C[i];
        }

        // Compute and store t1 o t.
        dm = dm.catenateZero(t1);

        // Compute t1 o map o t2^(-1).
        map = map.conjugate(t1);
    }
}


template <int N>
void DragtFinnMap<N>::
trackOrbit(FVector<double, 2 * N> &orbit, DragtFinnMap<N> &map) {
    FTps<double, 2 * N> fp = getGenerators();
    FMatrix<double, 2 * N, 2 * N> fM = getMatrix();
    FLieGenerator<double, N> f_1(fp, 1);

    // The generator to be moved along.
    FLieGenerator<double, N> g_1(1);
    for(int i = 0; i < N; ++i) {
        g_1[2*i+1] = f_1[2*i+1] + orbit[2*i+1];
        g_1[2*i+2] = f_1[2*i+2] - orbit[2*i];
    }

    // Initialise the map around the orbit.
    DragtFinnMap h;

    // Move orbit across the matrix.
    FMatrix<double, 2 * N, 2 * N> M = FLUMatrix<double, 2 * N>(fM).inverse();
    FLieGenerator<double, N> g_1t = g_1.transform(M);
    h += g_1t;

    // Terms of total rank 3.
    if(getOrder() >= 3) {
        FLieGenerator<double, N> f_33(fp, 3);
        FLieGenerator<double, N> f_32 = PoissonBracket(g_1t, f_33);
        FLieGenerator<double, N> f_31 = PoissonBracket(g_1t, f_32) / 2.0;
        h += f_33;
        h += f_32;
        h += f_31;

        // Terms of total rank 4.
        if(getOrder() >= 4) {
            FLieGenerator<double, N> f_44(fp, 4);
            FLieGenerator<double, N> f_43 = PoissonBracket(g_1t, f_44);
            FLieGenerator<double, N> f_42 = PoissonBracket(g_1t, f_43) / 2.0;
            FLieGenerator<double, N> f_41 = PoissonBracket(g_1t, f_42) / 3.0;
            h += f_44;
            h += f_43 + PoissonBracket(f_33, f_32) / 2.0;
            h += f_42 + PoissonBracket(f_31, f_33) / 2.0;
            h += f_41 + PoissonBracket(f_31, f_32) / 2.0;
        }
    }

    // Calculate the matrix part of the factored exponential.
    FTps<double, 2 * N> hp = h.getGenerators();
    FLieGenerator<double, N> h_2(hp, 2);
    FMatrix<double, 2 * N, 2 * N> hM = makeMatrix(h_2);
    hM = (((hM * (1. / 4.) + 1.) * hM * (1. / 3.) + 1.) * hM * (1. / 2.) + 1.) * hM + 1.;

    // Convert moved g_1 to orbit.
    FLieGenerator<double, N> h_1(hp, 1);
    for(int i = 0; i < N; ++i) {
        orbit[2*i] = - h_1[2*i+2];
        orbit[2*i+1] = h_1[2*i+1];
    }

    // Leave the higher generators in map.
    map.assign(hM * fM);
    map.assign(hp.filter(3, 4));
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
transform(const FLieGenerator<double, N> &g, int topOrder) {
    DragtFinnMap<N> H_trans(*this);
    int gOrder = g.getOrder();

    for(int order = 1; order <= topOrder; ++order) {
        FLieGenerator<double, N> pb(itsGenerators, order);
        int pbOrder = order;

        for(int power = 1; power <= topOrder; ++power) {
            pbOrder += gOrder - 2;
            if(pbOrder <= 0  ||  pbOrder > topOrder) break;
            pb = PoissonBracket(g, pb) / double(power);
            H_trans += pb;
        }
    }

    return H_trans;
}


template <int N>
double *DragtFinnMap<N>::
begin(int order) {
    return itsGenerators.begin() + FTpsData<2 * N>::getSize(order - 1);
}


template <int N>
const double *DragtFinnMap<N>::
begin(int order) const {
    return itsGenerators.begin() + FTpsData<2 * N>::getSize(order - 1);
}


template <int N>
double *DragtFinnMap<N>::
end(int order) {
    return itsGenerators.begin() + FTpsData<2 * N>::getSize(order);
}


template <int N>
const double *DragtFinnMap<N>::
end(int order) const {
    return itsGenerators.begin() + FTpsData<2 * N>::getSize(order);
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
catenateZero(const DragtFinnMap<N> &g) const {
    FMatrix<double, 2 * N, 2 * N> g_mat = g.getMatrix();
    FMatrix<double, 2 * N, 2 * N> f_mat = getMatrix();

    // The combined matrix.
    DragtFinnMap<N> h;
    h.assign(g_mat * f_mat);

    // Extract the generators.
    int order = std::max(getOrder(), g.getOrder());
    FMatrix<double, 2 * N, 2 * N> g_inv = FLUMatrix<double, 2 * N>(g_mat).inverse();
    FTps<double, 2 * N> f_gen =
        itsGenerators.filter(2, itsGenerators.getMaxOrder()).substitute(g_inv)
        + itsGenerators.filter(1, 1);
    FTps<double, 2 * N> g_gen = g.getGenerators();

    // Assign sum of the generators, the terms are now correct to order 3.
    h.assign(f_gen + g_gen);

    if(order >= 4) {
        FLieGenerator<double, N> f_3(f_gen, 3);
        FLieGenerator<double, N> g_3(g_gen, 3);
        FLieGenerator<double, N> pb_f3_g3 = PoissonBracket(f_3, g_3);
        h += pb_f3_g3 / 2.0;

        if(order >= 5) {
            FLieGenerator<double, N> f_4(f_gen, 4);
            FLieGenerator<double, N> pb_f3_f3_g3 = PoissonBracket(f_3, pb_f3_g3);
            FLieGenerator<double, N> pb_g3_f3_g3 = PoissonBracket(g_3, pb_f3_g3);
            FLieGenerator<double, N> pb_f4_g3    = PoissonBracket(f_4, g_3);
            h += pb_f4_g3 - pb_f3_f3_g3 / 6.0 - pb_g3_f3_g3 / 3.0;

            if(order >= 6) {
                FLieGenerator<double, N> f_5(f_gen, 5);
                FLieGenerator<double, N> g_4(g_gen, 4);
                h += (PoissonBracket(f_5, g_3)
                      + PoissonBracket(f_4, g_4) / 2.0
                      - PoissonBracket(g_3, pb_f4_g3) / 2.0
                      - PoissonBracket(f_4 + g_4, pb_f3_g3) / 4.0
                      + PoissonBracket(f_3, pb_f3_f3_g3) / 24.0
                      + PoissonBracket(f_3 + g_3, pb_g3_f3_g3) / 8.0);
            }
        }
    }

    return h;
}


template <int N>
DragtFinnMap<N> DragtFinnMap<N>::
factorize(const FTps<double, 2 * N> &HH) {
    FLieGenerator<double, N> H_2(HH, 2);

    for(int i = 0; i < N; ++i) {
        if(! H_2.derivative(2 * i).isZero()  &&
           ! H_2.derivative(2 * i + 1).isZero()) {
            return factorBerzForestIrwin(HH);
        }
    }

    return factorDouglas(HH);
}


template <int N>
FMatrix<double, 2 * N, 2 * N> DragtFinnMap<N>::
exponentiate(const FMatrix<double, 2 * N, 2 * N> &A) {
    FMatrix<double, 2 * N, 2 * N> B = A * A;
    FMatrix<double, 2 * N, 2 * N> C = A * B;

    for(int j = 0; j < 2 * N; ++j) {
        for(int i = 0; i < 2 * N; ++i) {
            B(i, j) = (A(i, j) - C(i, j) / 12.0) / 2.0;
            C(i, j) = - B(i, j);
        }

        B(j, j) = B(j, j) + 1.0;
        C(j, j) = C(j, j) + 1.0;
    }

    FLUMatrix<double, 2 * N> LU(C);
    LU.backSubstitute(B);
    return B;
}

template <int N>
FMatrix<double, 2 * N, 2 * N> DragtFinnMap<N>::
makeMatrix(const FLieGenerator<double, N> &f_2) {
    FMatrix<double, 2 * N, 2 * N> M;

    // OK, we need to know the Giorgilli numbering scheme for this
    // to work right.  I think I know it...

    unsigned gior_min = 2 * N + 1;
    for(unsigned i = 0, ig = gior_min ; i < N ; ++i) {
        M(2 * i + 1, 2 * i) = 2.0 * f_2[ig++];
        for(unsigned j = 2 * i + 1 ; j < 2 * N ; ++j) {
            M(2 * i + 1, j) = f_2[ig];
            M(j + 1 - 2 * (j % 2), 2 * i) = (1 - 2 * (static_cast<int>(j) % 2)) * f_2[ig++];
        }
        M(2 * i, 2 * i + 1) = -2.0 * f_2[ig++];
        for(unsigned j = 2 * i + 2 ; j < 2 * N ; ++j) {
            M(2 * i, j) = -f_2[ig];
            M(j + 1 - 2 * (j % 2), 2 * i + 1) = (1 - 2 * (static_cast<int>(j) % 2)) * f_2[ig++];
        }
    }

    return M;
}

template <>
FMatrix<double, 6, 6> DragtFinnMap<3>::
makeMatrix(const FLieGenerator<double, 3> &f_2);

template <int N>
void DragtFinnMap<N>::
move_g_1(DragtFinnMap<N> &f, DragtFinnMap<N> &g) {
    // Use method of L. Healy for moving the g_1 term over the f map.
    int order = f.getOrder();

    // Split first-order terms f_1 and g_1 from maps.
    FLieGenerator<double, N> f_1 = f.getGenerator(1);
    FMatrix<double, 2 * N, 2 * N> f_mat = f.getMatrix();
    FLieGenerator<double, N> g_1 = g.getGenerator(1);
    g.assign(g.getGenerators().filter(3, g.getOrder()));

    // Switch on highest order in left-hand side.
    DragtFinnMap<N> h;
    FLieGenerator<double, N> g_neg = - g_1;

    /*
      :FIXME: Is there a better solution for C++ < C++17
      C++17 provides a fallthrough attribute
     */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
    switch(order) {

        case 6:
            // Move g_1 accross f_6.
        {
            FLieGenerator<double, N> j = f.getGenerator(6);
            h += j;
            for(int i = 1; i <= 5; ++i) {
                j = PoissonBracket(g_neg, j) / double(i);
                h += j;
            }
        }
        // Fall through to order = 5.

        case 5:
            // Move g_1 accross f_5.
        {
            FLieGenerator<double, N> j = f.getGenerator(5);
            h += j;
            for(int i = 1; i <= 4; ++i) {
                j = PoissonBracket(g_neg, j) / double(i);
                h += j;
            }
        }
        // Fall through to order = 4.

        case 4:
            // Move g_1 accross f_4.
        {
            FLieGenerator<double, N> j_4 = f.getGenerator(4);
            FLieGenerator<double, N> j_3 = PoissonBracket(g_neg, j_4);
            FLieGenerator<double, N> j_2 = PoissonBracket(g_neg, j_3) / 2.0;
            FLieGenerator<double, N> j_1 = PoissonBracket(g_neg, j_2) / 3.0;
            h += j_4;
            h += j_3;
            h += j_2;
            h += j_1;

            if(order >= 6) {
                h += PoissonBracket(j_4, j_3);
                h += PoissonBracket(j_4, j_2);
                h += PoissonBracket(j_4, j_1);
                h += PoissonBracket(j_3, j_2);
                h += PoissonBracket(j_3, j_1);
                h += PoissonBracket(j_2, j_1);
            }
        }
        // Fall through to order = 3.

        case 3:

            // Move g_1 across f_3.
        {
            // factorize j_(1...3), result is k_(1...5).
            FLieGenerator<double, N> j_3 = f.getGenerator(3);
            FLieGenerator<double, N> j_2 = PoissonBracket(g_neg, j_3);
            FLieGenerator<double, N> j_1 = PoissonBracket(g_neg, j_2) / 2.0;

            FLieGenerator<double, N> k_4, k_5;
            FLieGenerator<double, N> k_3 = j_3;
            FLieGenerator<double, N> k_2 = j_2;
            FLieGenerator<double, N> k_1 = j_1;

            if(order >= 4) {
                FLieGenerator<double, N> pb_21 = PoissonBracket(j_2, j_1);
                FLieGenerator<double, N> pb_31 = PoissonBracket(j_3, j_1);
                FLieGenerator<double, N> pb_32 = PoissonBracket(j_3, j_2);
                k_1 += pb_21 / 2.0;
                k_2 += pb_31 / 2.0;
                k_3 += pb_32 / 2.0;

                if(order >= 5) {
                    FLieGenerator<double, N> pb_131 = PoissonBracket(j_1, pb_31);
                    FLieGenerator<double, N> pb_221 = PoissonBracket(j_2, pb_21);
                    FLieGenerator<double, N> pb_231 = PoissonBracket(j_2, pb_31);
                    FLieGenerator<double, N> pb_232 = PoissonBracket(j_2, pb_32);
                    FLieGenerator<double, N> pb_321 = PoissonBracket(j_3, pb_21);
                    FLieGenerator<double, N> pb_331 = PoissonBracket(j_3, pb_31);
                    FLieGenerator<double, N> pb_332 = PoissonBracket(j_3, pb_32);

                    k_1 -= (pb_131 - pb_221) / 6.0;
                    k_2 -= (pb_231 - 2.0 * pb_321) / 12.0;
                    k_3 -= (pb_232 - pb_331) / 6.0;
                    k_4 -= pb_332 / 12.0;

                    if(order >= 6) {
                        k_1 -= (3.0 * PoissonBracket(j_1, pb_321)
                                + PoissonBracket(j_2, pb_131 - pb_221)) / 24.0;
                        k_2 -= (PoissonBracket(j_2, pb_321)
                                - PoissonBracket(j_3, pb_131 - pb_221)) / 24.0;
                        k_3 += (PoissonBracket(j_2, pb_232 - 3.0 * pb_331)
                                + PoissonBracket(j_3, pb_231 + pb_321)) / 24.0;
                        k_4 += PoissonBracket(j_3, pb_232 - pb_331) / 24.0;
                        k_5 = PoissonBracket(j_3, pb_332) / 24.0;
                    }
                }
            }

            // Concatenate f_k with h = f_4 ... f_6.
            // Third and higher orders.
            h += k_1;
            h += k_2;
            h += k_3;

            if(order >= 4) {
                h += k_4;

                if(order >= 5) {
                    FLieGenerator<double, N> f_2 = h.getGenerator(2);
                    FLieGenerator<double, N> f_3 = h.getGenerator(3);
                    FLieGenerator<double, N> pb_22 = PoissonBracket(k_2, f_2);
                    FLieGenerator<double, N> pb_33 = PoissonBracket(k_3, f_3);

                    h += pb_22 / 2.0;
                    h += PoissonBracket(k_3, f_2);
                    h += pb_33 / 2.0;
                    h += k_5;

                    if(order >= 6) {
                        h += PoissonBracket(f_2, pb_22) / 12.0;
                        h -= PoissonBracket(f_3, pb_33) / 6.0;
                    }
                }
            }

            // Modify the first-order terms.
            g_1 += h.getGenerator(1);

            // Make matrix from h_2;
            FLieGenerator<double, N> h_2 = h.getGenerator(2);
            f.assign(exponentiate(makeMatrix(h_2)) * f.getMatrix());

            // Assign the modified terms to f.
            f.assign(h.getGenerators().filter(3, order));
        }

        default:
            // for order < 3 do nothing.
            // for order > 6 the method is not implemented.
            break;
    }
#pragma GCC diagnostic pop
    // Move first-order terms in g_1 across h_2 and combine with f_1.
    f.assign(f_1 + g_1.transform(f_mat));
}


template <int N>
int DragtFinnMap<N>::
orderModes(FMatrix<double, 2 * N, 2 * N> &V, FVector<std::complex<double>, 2 * N> &mu) {
    // Static constant.
    static const double tol = 1.0e-12;

    FMatrix<double, 2 * N, 2 * N>      tmat(V);
    FVector<std::complex<double>, 2 * N> tmu(mu);
    int nDim = 2 * N;
    int n_c = 0;
    int n_r = 0;

    for(int i = 0; i < 2 * N;) {
        if(std::abs(tmu[i]) < tol) {
            // Collect "coasting" modes in upper indices of V.
            nDim--;
            mu[nDim] = 0.0;
            for(int j = 0; j < 2 * N; ++j) V(j, nDim) = 0.0;
            V(nDim, nDim) = 1.0;
            i++;
        } else if(std::abs(std::imag(tmu[i])) < tol) {
            // Collect "unstable" modes in lower indices of tmat.
            if(n_r != i) {
                tmu[n_r] = tmu[i];
                for(int j = 0; j < 2 * N; ++j) tmat(j, n_r) = tmat(j, i);
            }
            n_r++;
            i++;
        } else if (i + 1 < 2 * N) {
            // Collect "stable" modes in lower indices of V.
            mu[n_c]   = tmu[i];
            mu[n_c+1] = tmu[i+1];

            // Normalise.
            double pb = 0.0;
            for(int j = 0; j < 2 * N; j += 2) {
                pb += tmat(j, i) * tmat(j + 1, i + 1) - tmat(j + 1, i) * tmat(j, i + 1);
            }
            double fact = 1.0 / sqrt(std::abs(pb));

            for(int j = 0; j < 2 * N; j++) {
                V(j, n_c)   = tmat(j, i)   * fact;
                V(j, n_c + 1) = tmat(j, i + 1) * fact;
            }

            i += 2;
            n_c += 2;
        }
    }

    // Order and copy "unstable" modes.
    for(int i = 0; i < n_r;) {
        int m = i + 1;
        for(; m < n_r && m < 2 * N; ++m) {
            if(std::abs(tmu[i] + tmu[m]) < tol) break;
        }

        if(m >= n_r) {
            throw LogicalError("FNormalForm::orderModes()",
                               "Cannot find pair of real eigenvalues.");
        }

        // Swap values to make pair.
        if(m != i + 1) {
            std::swap(tmu[m], tmu[i+1]);
            tmat.swapColumns(m, i + 1);
        }

        // Take positive eigenvalue first.
        int i1 = i;
        int i2 = i;
        if(std::real(tmu[i]) > 0.0) {
            ++i2;
        } else {
            ++i1;
        }
        mu[n_c]   = tmu[i2];
        mu[n_c+1] = tmu[i1];

        // Normalise this real pair.
        double pb = 0.0;
        for(int j = 0; j < 2 * N; j += 2) {
            pb += tmat(j, i1) * tmat(j + 1, i2) - tmat(j + 1, i1) * tmat(j, i2);
        }

        // Take factors such as to make the resulting pb = - 1.0;
        double fact1 = - 1.0 / sqrt(std::abs(2.0 * pb));
        double fact2 = (pb > 0.0) ? (- fact1) : fact1;
        for(int j = 0; j < 2 * N; ++j) {
            V(j, n_c)   = tmat(j, i1) * fact1 + tmat(j, i2) * fact2;
            V(j, n_c + 1) = tmat(j, i1) * fact1 - tmat(j, i2) * fact2;
        }

        n_c += 2;
        i += 2;
    }

    // Re-order eigenvector pairs by their main components.
    int modes = nDim / 2;
    for(int i = 0; i < nDim; i += 2) {
        // Find eigenvector pair with largest component in direction i.
        double big = 0.0;
        int k = i;
        for(int j = i; j < 2 * N; j += 2) {
            double c = std::abs(V(i, j) * V(i + 1, j + 1) - V(i, j + 1) * V(i + 1, j));
            if(c > big) {
                big = c;
                k   = j;
            }
        }

        if(k != i) {
            // Move eigenvector pair to its place.
            std::swap(mu[i],   mu[k]);
            std::swap(mu[i+1], mu[k+1]);
            V.swapColumns(i,     k);
            V.swapColumns(i + 1, k + 1);
        }
    }

    return modes;
}


template <int N>
std::ostream &operator<<(std::ostream &os, const DragtFinnMap<N> &map) {
    // next two lines require changes when g++ library becomes standard.
    std::streamsize old_prec = os.precision(14);
    os.setf(std::ios::floatfield, std::ios::scientific);

    os << "Map" << std::endl
       << map.getMatrix() << std::endl
       << map.getGenerators() << std::endl;

    os.precision(old_prec);
    // next line requires change when g++ library becomes standard.
    os.setf(std::ios::floatfield, std::ios::fixed);
    return os;
}

#endif // CLASSIC_DragtFinnMap_HH