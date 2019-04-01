#ifndef MAD_FNormalForm_HH
#define MAD_FNormalForm_HH

// ------------------------------------------------------------------------
// $RCSfile: FNormalForm.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FNormalForm
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2002/03/25 20:44:17 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <algorithm>
#include <complex>

template <class T, int, int> class FMatrix;
template <class T, int> class FTps;
template <class T, int> class FVps;


// Tolerance for accepting an eigenvalue.
namespace {
    const double tol = 1.0e-8;
}


// Class FNormalForm
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
class FNormalForm {

public:

    /// Constructor.
    //  Perform normal-form analysis of [b]map[/b].
    explicit FNormalForm(const FVps<double, N> &map);

    FNormalForm();
    FNormalForm(const FNormalForm &);
    ~FNormalForm();

    /// Get number of stable degrees of freedom.
    int degreesOfFreedom() const;

    /// Get normal-form map as a Lie transform.
    const FTps<double, N> &normalForm() const;

    /// Get normalising map as a Lie transform.
    const FTps<double, N> &normalisingMap() const;

    /// Get eigenvalues of the linear part as a complex vector.
    const FVector<std::complex<double>, N> &eigenValues() const;

    /// Get eigenvectors of the linear part in packed form.
    const FMatrix<double, N, N> &eigenVectors() const;

    /// Get anharmonicities as a matrix.
    FMatrix < double, N / 2, N / 2 > anharmonicity() const;

    /// Get invariant polynomial for mode i.
    FTps<double, N> invariant(int i) const;

protected:

    // Order the modes of the map and associate them to the planes.
    void orderModes(FVector<std::complex<double>, N>, FMatrix<double, N, N>);

private:

    // Not implemented.
    void operator=(const FNormalForm &);

    // Representation of the normal form analysis results.
    // ----------------------------------------------------------------------
    // Number of degrees of freedom.
    int freedom;

    // The factorised normalising map.
    FTps<double, N> A_Lie;

    // The factorised normal form map/
    FTps<double, N> N_Lie;

    // The vector of eigenvalues.
    FVector<std::complex<double>, N> lambda;

    // The matrix of eigenvectors.
    FMatrix<double, N, N> V;
};


#include "FixedAlgebra/FDoubleEigen.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMonomial.h"
#include "FixedAlgebra/FVps.h"
#include "Physics/Physics.h"
#include "Utilities/LogicalError.h"
#include <functional>




// Class FNormalForm
// ------------------------------------------------------------------------

template <int N>
FNormalForm<N>::FNormalForm():
    freedom(0), A_Lie(), N_Lie(), lambda(), V()
{}


template <int N>
FNormalForm<N>::FNormalForm(const FNormalForm &form):
    freedom(form.freedom), A_Lie(form.A_Lie), N_Lie(form.N_Lie),
    lambda(form.lambda), V(form.V)
{}


template <int N>
FNormalForm<N>::FNormalForm(const FVps<double, N> &M_scr):
    freedom(N / 2), A_Lie(0.0),
    N_Lie(0, 2, M_scr.getTopOrder() + 1), // Reserve space for second-order.
    lambda(), V() {
    if(N % 2 != 0) {
        throw LogicalError("FNormalForm::FNormalForm()", "Map dimension is odd.");
    }

    // Establish Jacobian matrix of map,
    // and find its eigenvalues and eigenvectors.
    // Order eigenvectors, attach modes to planes and find inverse of V.
    FMatrix<double, N, N> V_inv;
    {
        FMatrix<double, N, N> M = M_scr.linearTerms();
        FDoubleEigen<N> eigen(M, true);
        orderModes(eigen.eigenValues(), eigen.packedEigenVectors());
        FLUMatrix<double, N> lu(V);
        V_inv = lu.inverse();
    }

    // Initialise.
    FMatrix<double, N, N> Rot, R_dir, I_dir, R_inv;

    for(int i = 0; i < N; ++i) {
        Rot(i, i) = R_dir(i, i) = I_dir(i, i) = R_inv(i, i) = 1.0;
    }

    FMonomial<N> pows;

    for(int i = 0; i < 2 * freedom; i += 2) {
        // Store linear part of normal form.
        double c1, c2;

        if(std::abs(std::imag(lambda[i])) > tol) {
            c1 = c2 = std::arg(lambda[i]) / 2.0;
        } else {
            c1 = log(std::abs(lambda[i])) / 2.0;
            c2 = - c1;
        }

        pows[i] = 2;
        N_Lie.setCoefficient(pows, c1);
        pows[i] = 0;

        pows[i+1] = 2;
        N_Lie.setCoefficient(pows, c2);
        pows[i+1] = 0;

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

    // Fill in the values for the coasting modes.
    // Note that the matrices already have ones on the main diagonal.
    for(int i = 2 * freedom; i < N; i += 2) {
        // Store linear part of normal form.
        pows[i+1] = 2;
        N_Lie.setCoefficient(pows, 0.5);
        pows[i+1] = 0;
    }

    // Initialise the factored map.
    FVps<double, N> N_scr = V_inv * M_scr.substitute(V * Rot);
    FVps<double, N> N_acc;

    // Remove non-resonant terms order by order.
    int maxOrder = M_scr.getTopOrder();

    for(int omega = 2; omega < maxOrder; omega++) {
        // Compute the terms to be removed and store in f.
        FVps<double, N> temp = N_scr.substitute(N_acc).filter(omega, omega);
        FTps<double, N> f(0, omega + 1, omega + 1);

        for(int i = 0; i < N; i += 2) {
            f += temp[i+1].multiplyVariable(i) - temp[i].multiplyVariable(i + 1);
        }

        f /= double(omega + 1);

        // Set up the "phi" transformations.
        FTps<double, N> a(0, omega + 1, omega + 1);
        FTps<double, N> b(0, omega + 1, omega + 1);
        FTps<double, N> pi(0, omega + 1, omega + 1);
        FTps<double, N> t(0, omega + 1, omega + 1);

        for(int m = FTps<double, N>::getSize(omega);
            m < FTps<double, N>::getSize(omega + 1); m++) {
            const FMonomial<N> &index = FTpsData<N>::getExponents(m);
            std::complex<double> factor = 1.0;
            int count = 0;

            for(int j = 0; j < 2 * freedom; j += 2) {
                if(std::abs(std::imag(lambda[j])) > tol) count += index[j+1];
                factor *= (std::pow(lambda[j], index[j]) * std::pow(lambda[j+1], index[j+1]));
            }

            if(std::abs(1.0 - factor) < tol) {
                // Term cannot be removed.
                t.setCoefficient(m, 0.5);
            } else {
                // Term can be removed.
                factor = 1.0 / (1.0 - factor);
                a.setCoefficient(m, std::real(factor));
                b.setCoefficient(m, std::imag(factor));
            }

            pi.setCoefficient(m, std::pow(-1.0, int(count + 1) / 2));
        }

        // Compute cal_T^(-1) * f and T_omega.
        FTps<double, N> f1 = pi.scaleMonomials(f).substitute(R_dir);
        FTps<double, N> f2 = f1.substitute(I_dir);

        // Remove terms of order "omega" and accumulate map script(A).
        FTps<double, N> F_omega =
            (a.scaleMonomials(f1) + b.scaleMonomials(f2)).
            substitute(R_inv).scaleMonomials(pi);

        if(F_omega != 0.0) {
            FVps<double, N> A_omega = ExpMap(F_omega.substitute(Rot));
            FVps<double, N> A_invert = ExpMap(- F_omega);
            N_scr = A_invert.substitute(N_scr).substitute(A_omega);
            A_Lie -= F_omega;
        }

        // Accumulate map script(N).
        FTps<double, N> T_omega = t.scaleMonomials(f1 + f2).
                                  substitute(R_inv).scaleMonomials(pi);

        if(T_omega != 0.0) {
            N_acc = ExpMap(- T_omega, N_acc);
            N_Lie -= T_omega;
        }
    }
}


template <int N>
FNormalForm<N>::~FNormalForm()
{}


template <int N>
int FNormalForm<N>::degreesOfFreedom() const {
    return freedom;
}


template <int N>
FMatrix < double, N / 2, N / 2 > FNormalForm<N>::anharmonicity() const {
    FMatrix < double, N / 2, N / 2 > QQ;
    using Physics::pi;

    for(int mode1 = 0; mode1 < freedom; mode1++) {
        {
            FMonomial<N> power1, power2, power3;
            power1[2*mode1] = power3[2*mode1+1] = 4;
            power2[2*mode1] = power2[2*mode1+1] = 2;
            QQ(mode1, mode1) =
                - (3.0 * (N_Lie.getCoefficient(power1) +
                          N_Lie.getCoefficient(power3)) +
                   N_Lie.getCoefficient(power2)) / (4.0 * pi);
        }

        for(int mode2 = mode1 + 1; mode2 < freedom; mode2++) {
            FMonomial<N> power1, power2, power3, power4;
            power1[2*mode1]   = power1[2*mode2]   = 2;
            power2[2*mode1+1] = power2[2*mode2]   = 2;
            power3[2*mode1]   = power3[2*mode2+1] = 2;
            power4[2*mode1+1] = power4[2*mode2+1] = 2;
            QQ(mode1, mode2) = QQ(mode2, mode1) =
                                   - (N_Lie.getCoefficient(power1) +
                                      N_Lie.getCoefficient(power2) +
                                      N_Lie.getCoefficient(power3) +
                                      N_Lie.getCoefficient(power4)) / (4.0 * pi);
        }
    }

    return QQ;
}


template <int N>
FTps<double, N> FNormalForm<N>::invariant(int mode) const {
    FTps<double, N> a(1, 1);
    FTps<double, N> b(1, 1);

    for(int j = 0; j < N; j += 2) {
        a.setCoefficient(j + 1, - V(j + 1, 2 * mode));
        a.setCoefficient(j + 2,   V(j,  2 * mode));
        b.setCoefficient(j + 1, - V(j + 1, 2 * mode + 1));
        b.setCoefficient(j + 2,   V(j,  2 * mode + 1));
    }

    if(std::abs(std::imag(lambda[2*mode])) > tol) {
        return (a * a + b * b);
    } else {
        return (a * a - b * b);
    }
}


template <int N>
const FTps<double, N> &FNormalForm<N>::normalForm() const {
    return N_Lie;
}


template <int N>
const FTps<double, N> &FNormalForm<N>::normalisingMap() const {
    return A_Lie;
}


template <int N>
const FVector<std::complex<double>, N> &FNormalForm<N>::eigenValues() const {
    return lambda;
}


template <int N>
const FMatrix<double, N, N> &FNormalForm<N>::eigenVectors() const {
    return V;
}


template <int N>
void FNormalForm<N>::orderModes(FVector<std::complex<double>, N> tlam,
                                FMatrix<double, N, N> tmat) {
    // Static constant.
    int nDim = N;
    int n_c = 0;
    int n_r = 0;

    for(int i = 0; i < N;) {
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
            for(int j = 0; j < N; j += 2) {
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

            for(int j = 0; j < N; j++) {
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
            throw LogicalError("FNormalForm::orderModes()",
                               "Cannot find pair of real eigenvalues.");
        }

        // Swap values to make pair.
        if(m != i + 1) {
            std::swap(tlam[m], tlam[i+1]);
            tmat.swapColumns(m, i + 1);
        }

        // Normalise this real pair and convert it to Sinh/Cosh form.
        double pb = 0.0;
        for(int j = 0; j < N; j += 2) {
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

        for(int j = 0; j < N; ++j) {
            V(j, n_c)   = (tmat(j, i1) + tmat(j, i2)) / fact;
            V(j, n_c + 1) = (tmat(j, i1) - tmat(j, i2)) / fact;
        }

        n_c += 2;
        i += 2;
    }

    freedom = nDim / 2;
    if(nDim != 2 * freedom) {
        throw LogicalError("FNormalForm::orderModes()", "Map dimension is odd.");
    }

    // Re-order eigenvector pairs by their main components.
    for(int i = 0; i < nDim; i += 2) {
        // Find eigenvector pair with larges component in direction i/2.
        double big = 0.0;
        int k = i;
        for(int j = i; j < N; j += 2) {
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

            for(int j = 0; j < N; j++) {
                double real_part = re * V(j, i)     + im * V(j, i + 1);
                double imag_part = re * V(j, i + 1) - im * V(j, i);
                V(j, i)     = real_part;
                V(j, i + 1) = imag_part;
            }
        }
    }
}

#endif // MAD_FNormalForm_HH
