#ifndef MAD_DragtFinnNormalForm_HH
#define MAD_DragtFinnNormalForm_HH

// ------------------------------------------------------------------------
// $RCSfile: DragtFinnNormalForm.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DragtFinnNormalForm
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FVector.h"
#include <algorithm>
#include <complex>

template <int> class DragtFinnMap;
template <class T, int> class FLieGenerator;
template <class T, int, int> class FMatrix;


// Tolerance for accepting an eigenvalue.
namespace {
    const double tol = 1.0e-8;
}


// Class DragtFinnNormalForm
// ------------------------------------------------------------------------
/// Normal form of a truncated Taylor series map.
//  Compute from representation of a map which can be a nil-potent,
//  static or dynamic symplectic map. Implementation of an algorithm
//  described in
//  [center]
//  M. Berz, E. Forest and J. Irwin,[br]
//  Particle Accelerators, 1989, Vol. 24, pp. 91-107.
//  [/center]

template <int N>
class DragtFinnNormalForm {

public:

    /// Constructor.
    //  Perform normal-form analysis of [b]map[/b].
    explicit DragtFinnNormalForm(const DragtFinnMap<N> &map);

    DragtFinnNormalForm();
    DragtFinnNormalForm(const DragtFinnNormalForm &);
    ~DragtFinnNormalForm();

    /// Get number of stable degrees of freedom.
    int degreesOfFreedom() const;

    /// Get normal-form map as a Lie transform.
    const DragtFinnMap<N> &normalForm() const;

    /// Get normalising map as a Lie transform.
    const DragtFinnMap<N> &normalisingMap() const;

    /// Get eigenvalues of the linear part as a complex vector.
    const FVector<std::complex<double>, 2 * N> &eigenValues() const;

    /// Get eigenvectors of the linear part in packed form.
    const FMatrix<double, 2 * N, 2 * N> &eigenVectors() const;

    /// Get anharmonicities as a matrix.
    FMatrix<double, N, N> anharmonicity() const;

    /// Get invariant polynomial for mode i.
    FLieGenerator<double, N> invariant(int i) const;

protected:

    // Order the modes of the map and associate them to the planes.
    void orderModes(FVector<std::complex<double>, 2 * N>, FMatrix<double, 2 * N, 2 * N>);

private:

    // Not implemented.
    void operator=(const DragtFinnNormalForm &);

    // Representation of the normal form analysis results.
    // ----------------------------------------------------------------------
    // Number of degrees of freedom.
    int freedom;

    // The factorised normalising map.
    DragtFinnMap<N> A_scr;

    // The factorised normal form map/
    DragtFinnMap<N> N_scr;

    // The vector of eigenvalues.
    FVector<std::complex<double>, 2 * N> lambda;

    // The matrix of eigenvectors.
    FMatrix<double, 2 * N, 2 * N> V;
};


#include "FixedAlgebra/FDoubleEigen.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMonomial.h"
#include "FixedAlgebra/DragtFinnMap.h"
#include "Physics/Physics.h"
#include "Utilities/LogicalError.h"
#include <functional>


// Class DragtFinnNormalForm
// ------------------------------------------------------------------------

template <int N>
DragtFinnNormalForm<N>::DragtFinnNormalForm():
    freedom(0), A_scr(), N_scr(), lambda(), V()
{}


template <int N>
DragtFinnNormalForm<N>::DragtFinnNormalForm(const DragtFinnNormalForm &form):
    freedom(form.freedom), A_scr(form.A_scr), N_scr(form.N_scr),
    lambda(form.lambda), V(form.V)
{}


template <int N>
DragtFinnNormalForm<N>::DragtFinnNormalForm(const DragtFinnMap<N> &map):
    freedom(N), A_scr(), N_scr(), lambda(), V() {
    // Establish Jacobian matrix of map,
    // and find its eigenvalues and eigenvectors.
    // Order eigenvectors, attach modes to planes and find inverse of V.
    DragtFinnMap<N> M_scr(map);
    FMatrix<double, 2 * N, 2 * N> V_inv;
    {
        FMatrix<double, 2 * N, 2 * N> M = M_scr.getMatrix();
        FDoubleEigen<2 * N> eigen(M, true);
        orderModes(eigen.eigenValues(), eigen.packedEigenVectors());
        FLUMatrix<double, 2 * N> lu(V);
        V_inv = lu.inverse();
    }

    // Initialise.
    FMatrix<double, 2 * N, 2 * N> Rot, R_dir, I_dir, R_inv;

    for(int i = 0; i < N; ++i) {
        Rot(i, i) = R_dir(i, i) = I_dir(i, i) = R_inv(i, i) = 1.0;
    }

    FMonomial<2 * N> pows;

    for(int i = 0; i < 2 * freedom; i += 2) {
        // Set up the auxiliary matrices.
        R_dir(i, i) = R_dir(i, i + 1) = R_dir(i + 1, i) = 0.5;
        R_dir(i + 1, i + 1) = -0.5;
        R_inv(i, i) = R_inv(i, i + 1) = R_inv(i + 1, i) = 1.0;
        R_inv(i + 1, i + 1) = -1.0;

        if(std::abs(std::imag(lambda[i])) > tol) {
            // Complex eigenvalue pair.
            Rot(i, i)   = Rot(i + 1, i + 1) = std::real(lambda[i]);
            Rot(i, i + 1) = std::imag(lambda[i]);
            Rot(i + 1, i) = - std::imag(lambda[i]);

            I_dir(i, i) = I_dir(i + 1, i + 1) = 0.0;
            I_dir(i, i + 1) = I_dir(i + 1, i) = 1.0;
        } else {
            // Real eigenvalue pair.
            Rot(i, i) = Rot(i + 1, i + 1) = std::real(lambda[i] + lambda[i+1]) * 0.5;
            Rot(i, i + 1) = Rot(i + 1, i) = std::real(lambda[i] - lambda[i+1]) * 0.5;
        }
    }

    // Initialise the factored map.
    DragtFinnMap<N> N_scr;
    N_scr.assign(V_inv * M_scr.getMatrix() * V);
    N_scr.assign(M_scr.getGenerators().substitute(V));

    // Remove non-resonant terms order by order.
    int maxOrder = M_scr.getOrder();

    for(int omega = 2; omega < maxOrder; omega++) {
        // Compute the terms to be removed and store in f.
        FLieGenerator<double, N> f = N_scr.getGenerator(omega + 1);

        // Set up the "phi" transformations.
        FLieGenerator<double, N> a(omega + 1);
        FLieGenerator<double, N> b(omega + 1);
        FLieGenerator<double, N> pi(omega + 1);
        FLieGenerator<double, N> t(omega + 1);

        for(int m = f.getBottomIndex(); m < f.getTopIndex(); ++m) {
            const FMonomial<2 * N> &index = FTpsData<2 * N>::getExponents(m);
            std::complex<double> factor = 1.0;
            int count = 0;

            for(int j = 0; j < 2 * freedom; j += 2) {
                if(std::abs(std::imag(lambda[j])) > tol) count += index[j+1];
                factor *= (std::pow(lambda[j], index[j]) * std::pow(lambda[j+1], index[j+1]));
            }

            if(std::abs(1.0 - factor) < tol) {
                // Term cannot be removed.
                t[m] = 0.5;
            } else {
                // Term can be removed.
                factor = 1.0 / (1.0 - factor);
                a[m] = std::real(factor);
                b[m] = std::imag(factor);
            }

            pi[m] = std::pow(-1.0, int(count + 1) / 2);
        }

        // Compute cal_T^(-1) * f and T_omega.
        FLieGenerator<double, N> f1 = pi.scale(f).transform(R_dir);
        FLieGenerator<double, N> f2 = f1.transform(I_dir);

        // Remove terms of order "omega" and accumulate map script(A).
        FLieGenerator<double, N> F_omega =
            (a.scale(f1) + b.scale(f2)).transform(R_inv).scale(pi);
        A_scr.assign(F_omega);

        // Accumulate map script(N).
        N_scr = N_scr.transform(- F_omega, maxOrder);
    }
}


template <int N>
DragtFinnNormalForm<N>::~DragtFinnNormalForm()
{}


template <int N>
int DragtFinnNormalForm<N>::degreesOfFreedom() const {
    return freedom;
}


template <int N>
FMatrix<double, N, N> DragtFinnNormalForm<N>::anharmonicity() const {
    FTps<double, 2 * N> N_tps = N_scr.getGenerators();
    FMatrix<double, N, N> QQ;
    using Physics::pi;

    for(int mode1 = 0; mode1 < freedom; mode1++) {
        {
            FMonomial<2 * N> power1, power2, power3;
            power1[2*mode1] = power3[2*mode1+1] = 4;
            power2[2*mode1] = power2[2*mode1+1] = 2;
            QQ(mode1, mode1) = (- 3.0 / (4.0 * pi)) *
                               (N_tps[power1] + N_tps[power3] + N_tps[power2]);
        }

        for(int mode2 = mode1 + 1; mode2 < freedom; mode2++) {
            FMonomial<2 * N> power1, power2, power3, power4;
            power1[2*mode1]   = power1[2*mode2]   = 2;
            power2[2*mode1+1] = power2[2*mode2]   = 2;
            power3[2*mode1]   = power3[2*mode2+1] = 2;
            power4[2*mode1+1] = power4[2*mode2+1] = 2;
            QQ(mode1, mode2) = QQ(mode2, mode1) = (- 1.0 / (4.0 * pi)) *
                                                  (N_tps[power1] + N_tps[power2] + N_tps[power3] + N_tps[power4]);
        }
    }

    return QQ;
}


template <int N>
FLieGenerator<double, N> DragtFinnNormalForm<N>::invariant(int mode) const {
    FLieGenerator<double, N> a(1, 1);
    FLieGenerator<double, N> b(1, 1);

    for(int j = 0; j < 2 * N; j += 2) {
        a[j+1] = - V(j + 1, 2 * mode);
        a[j+2] =   V(j,  2 * mode);
        b[j+1] = - V(j + 1, 2 * mode + 1);
        b[j+2] =   V(j,  2 * mode + 1);
    }

    if(std::abs(std::imag(lambda[2*mode])) > tol) {
        return (a * a + b * b);
    } else {
        return (a * a - b * b);
    }
}


template <int N>
const DragtFinnMap<N> &DragtFinnNormalForm<N>::normalForm() const {
    return N_scr;
}


template <int N>
const DragtFinnMap<N> &DragtFinnNormalForm<N>::normalisingMap() const {
    return A_scr;
}


template <int N>
const FVector<std::complex<double>, 2 * N> &DragtFinnNormalForm<N>::eigenValues() const {
    return lambda;
}


template <int N>
const FMatrix<double, 2 * N, 2 * N> &DragtFinnNormalForm<N>::eigenVectors() const {
    return V;
}


template <int N>
void DragtFinnNormalForm<N>::orderModes
(FVector<std::complex<double>, 2 * N> tlam, FMatrix<double, 2 * N, 2 * N> tmat) {
    // Static constant.
    int nDim = 2 * N;
    int n_c = 0;
    int n_r = 0;

    for(int i = 0; i < 2 * N;) {
        if(std::abs(tlam[i] - 1.0) < tol) {
            // Collect "coasting" modes in upper indices of V.
            nDim--;
            lambda[nDim] = 1.0;
            std::fill(V.col_begin(nDim), V.col_end(nDim), 0.0);
            V(nDim, nDim) = 1.0;
            i++;
        } else if(std::abs(std::imag(tlam[i])) < tol) {
            // Collect "unstable" modes in lower indices of tmat.
            if(n_r != i) {
                tlam[n_r] = tlam[i];
                std::copy(tmat.col_begin(i), tmat.col_end(i), tmat.col_begin(n_r));
            }
            n_r++;
            i++;
        } else {
            // Collect "stable" modes in lower indices of V.
            double pb = 0.0;
            for(int j = 0; j < 2 * N; j += 2) {
                pb += tmat(j, i) * tmat(j + 1, i + 1) - tmat(j + 1, i) * tmat(j, i + 1);
            }

            int i1 = i;
            int i2 = i;
            if(pb > 0.0) {
                ++i2;
            } else {
                ++i1;
            }

            lambda[n_c]   = tlam[i1];
            lambda[n_c+1] = tlam[i2];

            double fact = sqrt(std::abs(pb));

            for(int j = 0; j < 2 * N; j++) {
                V(j, n_c)   = tmat(j, i1) / fact;
                V(j, n_c + 1) = tmat(j, i2) / fact;
            }

            i += 2;
            n_c += 2;
        }
    }

    // Order and copy "unstable" modes.
    for(int i = 0; i < n_r;) {
        int m;
        for(m = i + 1; m < n_r; ++m) {
            if(std::abs(tlam[i] * tlam[m] - 1.0) < tol) break;
        }

        if(m >= n_r) {
            throw LogicalError("DragtFinnNormalForm::orderModes()",
                               "Cannot find pair of real eigenvalues.");
        }

        // Swap values to make pair.
        if(m != i + 1) {
            std::swap(tlam[m], tlam[i+1]);
            tmat.swapColumns(m, i + 1);
        }

        // Normalise this real pair and convert it to Sinh/Cosh form.
        double pb = 0.0;
        for(int j = 0; j < 2 * N; j += 2) {
            pb += tmat(j, i) * tmat(j + 1, i + 1) - tmat(j + 1, i) * tmat(j, i + 1);
        }

        int i1 = i;
        int i2 = i;
        if(pb > 0.0) {
            ++i2;
        } else {
            ++i1;
        }

        lambda[n_c]   = tlam[i1];
        lambda[n_c+1] = tlam[i2];

        double fact = sqrt(std::abs(2.0 * pb));

        for(int j = 0; j < 2 * N; ++j) {
            V(j, n_c)   = (tmat(j, i1) + tmat(j, i2)) / fact;
            V(j, n_c + 1) = (tmat(j, i1) - tmat(j, i2)) / fact;
        }

        n_c += 2;
        i += 2;
    }

    freedom = nDim / 2;
    if(nDim != 2 * freedom) {
        throw LogicalError("DragtFinnNormalForm::orderModes()",
                           "Map dimension is odd.");
    }

    // Re-order eigenvector pairs by their main components.
    for(int i = 0; i < nDim; i += 2) {
        // Find eigenvector pair with larges component in direction i/2.
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
            std::swap(lambda[i],   lambda[k]);
            std::swap(lambda[i+1], lambda[k+1]);
            V.swapColumns(i, k);
            V.swapColumns(i + 1, k + 1);
        }

        if(std::abs(std::imag(lambda[i])) > tol) {
            // Rotate complex eigenvectors to make their main components real.
            double re = V(i, i)     / sqrt(V(i, i) * V(i, i) + V(i, i + 1) * V(i, i + 1));
            double im = V(i, i + 1) / sqrt(V(i, i) * V(i, i) + V(i, i + 1) * V(i, i + 1));

            for(int j = 0; j < 2 * N; j++) {
                double real_part = re * V(j, i)     + im * V(j, i + 1);
                double imag_part = re * V(j, i + 1) - im * V(j, i);
                V(j, i)     = real_part;
                V(j, i + 1) = imag_part;
            }
        }
    }
}

#endif // MAD_DragtFinnNormalForm_HH
