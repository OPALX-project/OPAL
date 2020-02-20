// ------------------------------------------------------------------------
// $RCSfile: Migrad.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Migrad
//   This class encapsulates a minimisation according to the MIGRAD method,
//   Minimization by a gradient method due to Davidon/Fletcher/Powell.
//   (Computer Journal 13, 317 (1970).
//   Algorithm taken from the MINUIT package.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/Migrad.h"
#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Attributes/Attributes.h"
#include "Match/Match.h"
#include "Utilities/Round.h"
#include <algorithm>
#include <cmath>
#include <cfloat>

using std::max;
using std::min;


// Class Migrad
// ------------------------------------------------------------------------

// The attributes of class Migrad.
namespace {
    enum {
        TOLERANCE,   // The desired tolerance.
        CALLS,       // The maximum number of calls to the matching functions.
        SIZE
    };
}


Migrad::Migrad():
    Action(SIZE, "MIGRAD",
           "The \"MIGRAD\" sub-command adjusts parameters according to the "
           "\"MIGRAD\" method taken from MINUIT.") {
    itsAttr[TOLERANCE] = Attributes::makeReal
                         ("TOLERANCE", "The desired tolerance", 1.0e-6);
    itsAttr[CALLS] = Attributes::makeReal
                     ("CALLS", "Maximum number of calls to the matching functions", 1000.);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


Migrad::Migrad(const std::string &name, Migrad *parent):
    Action(name, parent)
{}


Migrad::~Migrad()
{}


Migrad *Migrad::clone(const std::string &name) {
    return new Migrad(name, this);
}


void Migrad::execute() {
    // Built-in constants.
    static const char method[] = "MIGRAD";

    // fetch command attributes.
    double tol = max(Attributes::getReal(itsAttr[TOLERANCE]), DBL_EPSILON);
    int nfcnmax = int(Round(Attributes::getReal(itsAttr[CALLS])));

    // Working vectors.
    int n = Match::block->countVariables();
    Vector<double> x;
    Vector<double> f;
    double fnorm;
    MatchState state = START;

    // Initial function evaluation.
    Match::block->getVariables(x);
    if(! Match::block->evaluate(x, f)) {
        state = FAILED;
    } else if((fnorm = euclidean_norm(f)) <= tol) {
        state = CONVERGED;
    } else {
        // Initialise algorithm.
        fmin = f * f;
        Vector<double> g(n);
        Vector<double> g2(n);
        Matrix<double> W(n, n);
        int restarts = 0;
        int strategy = 2;
        int icovar = 0;
        bool posdef = false;
        state = START;

        // Outer loop: start/restart the algorithm.
        do {
            // Initialise covariance matrix.
            if(strategy == 2  || ( strategy > 2 && icovar < 2 )) {
                hessenberg(x, f, W, g, g2);
                posdef = false;
            } else {
                derivatives(x, f, g, g2);
                if(icovar < 2) {
                    for(int i = 0; i < n; ++i) {
                        for(int j = 0; j < n; ++j) W(i, j) = 0.0;
                        if(g2(i) == 0.0) g2(i) = 1.0;
                        W(i, i) = 1.0 / g2(i);
                    }
                }

                posdef = true;
            }

            // Initialize for first iteration.
            int improv = 0;
            double edm = min(0.5 * g * (W * g), fmin);

            // Print after initialization.
            Match::block->print(method, state);

            // Inner loop: actual iteration.
            while(Match::block->getCallCount() < nfcnmax) {
                // first derivatives all zero?
                state = CHECK;
                if(g * g < DBL_EPSILON * fnorm) break;

                // find step size according to Newton's method.
                Vector<double> gsave(g);
                Vector<double> s = - (W * g);

                // If s * g >= 0 matrix is not positive definite.
                if(s *g >= 0.0) {
                    if(posdef) {
                        // If already forced positive definite, restart algorithm.
                        restarts++;
                        state = FAILED;
                        if(restarts > strategy) break;
                        state = RESTART;
                        break;
                    } else {
                        // Otherwise force covariance matrix positive definite.
                        invertSymmetric(W);
                        forcePositiveDefinite(W);
                        invertSymmetric(W);
                        posdef = true;
                        continue;
                    }
                }

                // Search for minimum along predicted line.
                state = lineSearch(x, s, f, tol);
                if(state != PROGRESS) break;

                // find gradient in new point.
                derivatives(x, f, g, g2);
                posdef = false;

                // Inner loop.
                double sy;
                double yWy;
                Vector<double> y;
                Vector<double> Wy;

                // Inner loop, try estimating the error.
                while(true) {
                    // Estimated distance to minimum.
                    edm = min(0.5 * g * (W * g), fmin);
                    y = g - gsave;
                    sy = s * y;
                    Wy = W * y;
                    yWy = Wy * y;

                    // Test for convergence and print-out.
                    state = CONVERGED;
                    if(edm >= 0.0  &&  edm < 1.0e-4 * tol) break;
                    Match::block->print(method, PROGRESS);

                    // If covariance matrix is not positive definite,
                    // force it positive definite and retry.
                    state = PROGRESS;
                    if(edm >= 0.0  &&  yWy > 0.0) break;
                    icovar = 0;
                    state = CONVERGED;
                    if(posdef) break;
                    invertSymmetric(W);
                    forcePositiveDefinite(W);
                    invertSymmetric(W);
                    posdef = true;
                }

                // If iteration has converged, break out of outer loop.
                if(state == CONVERGED) break;

                // Update covariance matrix.
                for(int i = 0; i < n; ++i) {
                    for(int j = 0; j < n; ++j) {
                        W(i, j) += s(i) * s(j) / sy - Wy(i) * Wy(j) / yWy;
                    }
                }

                if(sy > yWy) {
                    Vector<double> flnu = s / sy - Wy / yWy;

                    for(int i = 0; i < n; ++i) {
                        for(int j = 0; j < n; ++j) {
                            W(i, j) += 2.0 * yWy * flnu(i) * flnu(j);
                        }
                    }
                    improv++;
                    if(improv >= n) icovar = 3;
                }

                // End of outer loop.
            }

            // End of algorithm, check for special cases.
            if(state == CHECK) {
                state = CONVERGED;
                if(strategy >= 2  || (strategy == 1 && icovar < 3)) {
                    hessenberg(x, f, W, g, g2);
                    posdef = false;
                    if(edm > 1.0e-3 * tol) state = RESTART;
                }
            }
        } while(state == RESTART);

        Match::block->setVariables(x);
    }

    // Common exit point; final print-out.
    Match::block->print(method, state);
}



void Migrad::hessenberg(Vector<double> &x, Vector<double> &f,
                        Matrix<double> &W, Vector<double> &g,
                        Vector<double> &g2) {
    double eps = sqrt(DBL_EPSILON);
    int n = Match::block->countVariables();
    int m = Match::block->countFunctions();
    Vector<double> xstep(n), x1(x), x2(x), f1(m), f2(m), fu(n), fl(n);

    for(int i = 0; i < n; ++i) {
        xstep(i) = eps * std::max(std::abs(x(i)), 1.0);

        for(int icycle = 1; icycle <= 10; ++icycle) {
            x1(i) = x(i) - xstep(i);
            x2(i) = x(i) + xstep(i);

            if(Match::block->evaluate(x2, f2)  &&  Match::block->evaluate(x1, f1)) {
                fl(i) = f1 * f1;
                fu(i) = f2 * f2;
                break;
            } else {
                xstep(i) = 0.5 * xstep(i);
            }
        }

        g(i)  = (fu(i) - fl(i)) / (2.0 * xstep(i));
        g2(i) = (fu(i) - 2.0 * fmin + fl(i)) / (xstep(i) * xstep(i));
        if(g2(i) == 0.0) g2(i) = 1.0;
        W(i, i) = g2(i);

        // Off-diagonal elements.
        x1(i) = x(i);
        x2(i) = x(i) + xstep(i);
        for(int j = 0; j < i; ++j) {
            x2(j) = x(j) + xstep(j);
            if(Match::block->evaluate(x2, f)) {
                W(j, i) = W(i, j) =
                              (f * f + fmin - fu(i) - fu(j)) / (xstep(i) * xstep(j));
            } else {
                W(j, i) = W(i, j) = 0.0;
            }
            x2(j) = x(j);
        }
        x2(i) = x(i);
    }

    // Restore original point.
    Match::block->setVariables(x);

    // ensure positive definiteness and invert.
    forcePositiveDefinite(W);
    invertSymmetric(W);
}


void Migrad::derivatives(Vector<double> &x, Vector<double> &/*f*/,
                         Vector<double> &g, Vector<double> &g2) {
    int n = Match::block->countVariables();
    int m = Match::block->countFunctions();
    Vector<double> xstep(n), x1(x), x2(x);
    double eps = sqrt(DBL_EPSILON);

    for(int i = 0; i < n; ++i) {
        xstep(i) = eps * std::max(std::abs(x(i)), 1.0);

        for(int icycle = 1; icycle <= 10; ++icycle) {
            Vector<double> f1(m), f2(m);
            x1(i) = x(i) - xstep(i);
            x2(i) = x(i) + xstep(i);

            if(Match::block->evaluate(x1, f1)  &&  Match::block->evaluate(x2, f2)) {
                double fl = f1 * f1;
                double fu = f2 * f2;
                g(i)  = (fu - fl) / (2.0 * xstep(i));
                g2(i) = (fu - 2.0 * fmin + fl) / (xstep(i) * xstep(i));
                if(g2(i) == 0.0) g2(i) = 1.0;
                break;
            }

            xstep(i) = 0.5 * xstep(i);
            x1(i) = x2(i) = x(i);
        }
    }

    Match::block->setVariables(x);
}


MatchState Migrad::lineSearch(Vector<double> &x, Vector<double> &s,
                              Vector<double> &f, double tol) {
    // Built-in constants.
    static const double alpha  = 2.0;
    static const double small = 0.05;
    static const int maxpoints = 12;
    int n = Match::block->countVariables();

    // find limits for the parameter lambda.
    double lambda_min = 0.0;
    for(int i = 0; i < n; ++i) {
        if(s(i) != 0.0) {
            double ratio = std::abs(x(i) / s(i));
            if(lambda_min == 0.0  ||  ratio < lambda_min) lambda_min = ratio;
        }
    }
    if(lambda_min == 0.0) lambda_min = DBL_EPSILON;
    double lambda_max = lambda_min = lambda_min * DBL_EPSILON;

    // Store the initial point x for lambda = 0;
    // Declare the optimum point (lambda_opt, funct_opt).
    double lambda_val[3], funct_val[3];
    lambda_val[0] = 0.0;
    funct_val[0]  = fmin;
    double lambda_opt = 0.0;
    double funct_opt  = fmin;
    Vector<double> f_opt(f);

    // Other parameters.
    double over_all = 1000.0;
    double under_all = - 100.0;

    // find a valid point in direction lambda * s.
    double lambda = 1.0;
    while(! Match::block->evaluate(x + lambda * s, f)) {
        lambda *= 0.5;
        if(lambda <= lambda_min) break;
    }

    lambda_val[1] = lambda;
    funct_val[1] = f * f;

    if(funct_val[1] < funct_opt) {
        lambda_opt = lambda_val[1];
        funct_opt = funct_val[1];
        f_opt = f;
    }

    // If lambda had to be reduced, use the optimal point.
    // Then compute function for move by lambda * s / 2.
    // If this point is acceptable, continue.
    if(lambda >= 1.0  &&
       Match::block->evaluate(x + 0.5 * lambda * s, f)) {

        // Store the second value.
        lambda_val[2] = 0.5 * lambda;
        funct_val[2] = f * f;
        if(funct_val[2] < funct_opt) {
            lambda_opt = lambda_val[2];
            funct_opt = funct_val[2];
            f_opt = f;
        }

        // Main iteration loop.
        for(int npts = 2; npts <= maxpoints;) {

            // Begin iteration.
            lambda_max = max(lambda_max, alpha * std::abs(lambda_opt));

            // Quadratic interpolation using three points.
            double s21 = lambda_val[1] - lambda_val[0];
            double s32 = lambda_val[2] - lambda_val[1];
            double s13 = lambda_val[0] - lambda_val[2];
            double den = s21 * s32 * s13;
            double c2 =
                (s32 * funct_val[0] + s13 * funct_val[1] + s21 * funct_val[2]) / den;
            double c1 =
                ((lambda_val[2] + lambda_val[1]) * s32 * funct_val[0] +
                 (lambda_val[0] + lambda_val[2]) * s13 * funct_val[1] +
                 (lambda_val[1] + lambda_val[0]) * s21 * funct_val[2]) / den;

            if(c2 >= 0.0) {
                lambda = (c1 > 2.0 * c2 * lambda_opt) ?
                         (lambda_opt + lambda_max) : (lambda_opt - lambda_max);
            } else {
                lambda = c1 / (2.0 * c2);
                lambda = (lambda > lambda_opt + lambda_max) ?
                         (lambda_opt + lambda_max) : (lambda_opt - lambda_max);
            }

            if(lambda > 0.0) {
                if(lambda > over_all) lambda = over_all;
            } else {
                if(lambda < under_all) lambda = under_all;
            }

            // Inner loop: cut interval, if points are no good.
            bool FAILED;
            int nvmax = 0;
            double f3 = 0.0;
            while(true) {

                // Reject new point, if it coincides with a previous one.
                FAILED = true;
                double tol9 = small * max(1.0, lambda);
                for(int i = 0; i < 3; ++i) {
                    if(std::abs(lambda - lambda_val[i]) < tol9) break;
                }

                // Compute function for interpolated point.
                // Reject, if too many points computed, or if point is invalid.
                if(++npts > maxpoints  ||
                   ! Match::block->evaluate(x + lambda * s, f)) break;
                f3 = f * f;

                // find worst point of previous three.
                nvmax = 0;
                if(funct_val[1] > funct_val[nvmax]) nvmax = 1;
                if(funct_val[2] > funct_val[nvmax]) nvmax = 2;

                // Accept new point, if it is an improvement.
                FAILED = false;
                if(f3 < funct_val[nvmax]) break;

                // Otherwise cut the interval.
                if(lambda > lambda_opt) over_all = min(over_all, lambda - small);
                if(lambda < lambda_opt) under_all = max(under_all, lambda + small);
                lambda = 0.5 * (lambda + lambda_opt);
            } // end of inner loop.

            // Break out of outer loop, if inner loop FAILED.
            if(FAILED) break;

            // Otherwise new point replaces the previous worst point.
            lambda_val[nvmax] = lambda;
            funct_val[nvmax] = f3;
            if(f3 < funct_opt) {
                lambda_opt = lambda;
                funct_opt = f3;
                f_opt = f;
            } else {
                if(lambda > lambda_opt) over_all = min(over_all, lambda - small);
                if(lambda < lambda_opt) under_all = max(under_all, lambda + small);
            }

            // End of: for (int npts = 2; ... )
        }

        // End of: if (lambda >= 1.0)
    }

    // Common exit point: return best point and step used.
    fmin = funct_opt;
    s *= lambda_opt;
    x += s;
    Match::block->evaluate(x, f);

    // Return failure indication.
    if(fmin < tol) {
        return CONVERGED;
    } else if(lambda_opt == 0) {
        return FAILED;
    } else if(s * s > DBL_EPSILON * (x * x)) {
        return PROGRESS;
    } else {
        return ACCURACY_LIMIT;
    }
}


int Migrad::symmetricEigen(const Matrix<double> &V, Vector<double> &eigen) {
    // Built-in constants.
    static const double eps = DBL_EPSILON;
    static const int itmax = 15;
    int n = Match::block->countVariables();

    // Working vectors.
    eigen = Vector<double>(n);

    // Special cases.
    if(n <= 0) {
        return 0;
    } else if(n == 1) {
        eigen(0) = V(0, 0);
        return 1;
    } else if(n == 2) {
        double f = V(0, 0) + V(1, 1);
        double g = sqrt(pow(V(0, 0) - V(1, 1), 2) + 4.0 * pow(V(1, 0), 2));
        eigen(0) = (f - g) / 2.0;
        eigen(1) = (f + g) / 2.0;
        return 2;
    }

    // Dimension is at least 3, reduce to tridiagonal form.
    Vector<double> work(n);
    Matrix<double> A(V);
    for(int i = n; i-- >= 2;) {
        double g = 0.0;
        for(int j = 0; j < i - 1; ++j) g += A(j, i) * A(j, i);
        eigen(i) = A(i, i);

        if(g == 0.0) {
            work(i) = A(i, i - 1);
        } else {
            double h = g + pow(A(i, i - 1), 2);
            work(i) = (A(i, i - 1) >= 0.0) ? sqrt(h) : (- sqrt(h));
            h += A(i, i - 1) * work(i);
            A(i, i - 1) += work(i);
            double f = 0.0;

            for(int j = 0; j < i; ++j) {
                g = 0.0;
                for(int k = 0; k < i - 1; ++k) g += A(k, i) * A(k, j);
                work(j) = g / h;
                f += work(j) * A(i, j);
            }

            for(int j = 0; j < i; ++j) {
                work(j) -= (f / (h + h)) * A(i, j);
                for(int k = 0; k <= j; ++k) {
                    A(j, k) -= A(i, j) * work(k) - work(j) * A(i, k);
                }
            }
        }
    }

    work(1) = A(1, 0);
    work(0) = 0.0;
    eigen(1) = A(1, 1);
    eigen(0) = A(0, 0);

    // Iterate on tridiagonal matrix.
    for(int i = 1; i < n; ++i) work(i - 1) = work(i);

    work(n - 1) = 0.0;
    double f = 0.0;
    double b = 0.0;

    for(int l = 0; l < n; ++l) {
        b = max(eps * (std::abs(eigen(l)) + std::abs(work(l))), b);
        int m;
        for(m = 0; m < n - 1; ++m) if(std::abs(work(m)) <= b) break;

        if(m != l) {
            for(int iter = 1; iter <= itmax; ++iter) {
                double p = (eigen(l + 1) - eigen(l)) / (2.0 * work(l));
                double r = (std::abs(p) > 1.0e10) ? std::abs(p) : sqrt(p * p + 1.0);
                double h = eigen(l) - work(l) / ((p > 0) ? (p + r) : (p - r));
                for(int i = l; i < n; ++i) eigen(i) -= h;
                f += h;
                p = eigen(m);
                double c = 1.0;
                double s = 0.0;
                for(int i = m; i-- > l;) {
                    double g = c * work(i);
                    h = c * p;
                    r = sqrt(work(i) * work(i) + p * p);
                    work(i + 1) = s * r;
                    s = work(i) / r;
                    c = p / r;
                    p = c * eigen(i) - s * g;
                    eigen(i + 1) = h + s * (c * g + s * eigen(i));
                }
                work(l) = s * p;
                eigen(l) = c * p;
                if(std::abs(work(l)) <= b) break;
            }
            return l;
        }

        double p = eigen(l) + f;
        for(int i = l; i >= 2; --i) {
            if(p >= eigen(i - 1)) {
                eigen(0) = p;
                return n;
            }
            eigen(i) = eigen(i - 1);
        }
    }

    return n;
}


bool Migrad::invertSymmetric(Matrix<double> &A) {
    // allocate working space.
    int n = Match::block->countVariables();
    Vector<double> scale(n);
    Vector<double> pivot(n);
    Vector<double> s(n);

    // Scale upper triangle.
    for(int i = 0; i < n; ++i) {
        double si = A(i, i);
        if(si <= 0.0) return false;
        scale(i) = 1.0 / sqrt(si);
    }

    for(int i = 0; i < n; ++i) {
        for(int j = i; j < n; ++j) A(i, j) *= scale(i) * scale(j);
    }

    // Invert upper triangle.
    for(int i = 0; i < n; ++i) {
        if(A(i, i) == 0.0) return false;
        pivot(i) = 1.0;
        s(i) = 1.0 / A(i, i);
        A(i, i) = 0.0;

        for(int j = 0; j < n; ++j) {
            if(j < i) {
                pivot(j) = A(j, i);
                s(j) = pivot(j) * s(i);
                A(j, i) = 0.0;
            } else if(j > i) {
                pivot(j) = A(i, j);
                s(j) = - pivot(j) * s(i);
                A(i, j) = 0.0;
            }
        }

        for(int j = 0; j < n; ++j) {
            for(int k = j; k < n; ++k) A(j, k) += pivot(j) * s(k);
        }
    }

    // Rescale upper triangle and symmetrize.
    for(int i = 0; i < n; ++i) {
        for(int j = i; j < n; ++j) {
            A(j, i) = A(i, j) *= scale(i) * scale(j);
        }
    }

    return true;
}


void Migrad::forcePositiveDefinite(Matrix<double> &V) {
    // Built-in constant.
    static const double eps = 1.0e-3;

    // find eigenvalues.
    int n = Match::block->countVariables();
    Vector<double> eigen(n);
    symmetricEigen(V, eigen);

    // Enforce positive definiteness.
    double pmin = eigen(0);
    double pmax = eigen(0);

    for(int i = 0; i < n; ++i) {
        if(eigen(i) < pmin) pmin = eigen(i);
        if(eigen(i) > pmax) pmax = eigen(i);
    }

    pmax = std::max(std::abs(pmax), 1.0);
    if(pmin <= DBL_EPSILON * pmax) V += eps * pmax - pmin;
}