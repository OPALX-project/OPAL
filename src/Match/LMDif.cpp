// ------------------------------------------------------------------------
// $RCSfile: LMDif.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LMDif
//   This class encapsulates a minimisation according to the method used
//   by Levenberg and Marquardt.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/LMDif.h"
#include "Algebra/Matrix.h"
#include "Algebra/QRSolver.h"
#include "Algebra/Vector.h"
#include "Attributes/Attributes.h"
#include "Match/Match.h"
#include "Utilities/Round.h"
#include <algorithm>
#include <cmath>
#include <cfloat>

using std::max;
using std::min;


// Class LMDif
// ------------------------------------------------------------------------

// The attributes of class LMDif.
namespace {
    enum {
        TOLERANCE,   // The desired tolerance.
        CALLS,       // The maximum number of calls to the matching functions.
        SIZE
    };
}


LMDif::LMDif():
    Action(SIZE, "LMDIF",
           "The \"LMDIF\" sub-command adjusts parameters according to the "
           "Levenberg-Marquardt method.") {
    itsAttr[TOLERANCE] = Attributes::makeReal
                         ("TOLERANCE", "The desired tolerance", 1.0e-6);
    itsAttr[CALLS] = Attributes::makeReal
                     ("CALLS", "Maximum number of calls to the matching functions", 1000.);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


LMDif::LMDif(const std::string &name, LMDif *parent):
    Action(name, parent) {
    fraction = 1.0e-6;
    delta = 10.0;
}


LMDif::~LMDif()
{}


LMDif *LMDif::clone(const std::string &name) {
    return new LMDif(name, this);
}


void LMDif::execute() {
    // Built-in constants.
    static const char method[] = "LMDIF";

    // Fetch command attributes.
    double tol = max(DBL_EPSILON, Attributes::getReal(itsAttr[TOLERANCE]));
    int nfcnmax = int(Round(Attributes::getReal(itsAttr[CALLS])));

    // Initialize.
    double fnorm = 0.0;
    MatchState state = START;

    // Working storage.
    Vector<double> x, f;
    Array1D<double> D;
    int m = Match::block->countFunctions();
    int n = Match::block->countVariables();
    jacobian = Matrix<double>(m, n);

    // Initial function evaluation.
    Match::block->getVariables(x);
    if(! Match::block->evaluate(x, f)) {
        state = FAILED;
    } else if((fnorm = euclidean_norm(f)) <= tol) {
        state = CONVERGED;
    } else {
        // Initialize iterations.
        mu = 0.0;
        double dxnorm = 0.0;
        double ratio = 1.0;
        int iter = 1;
        state = START;

        // Outer iteration.
        do {
            // Initial print-out.
            Match::block->print(method, state);

            // Recalculate the Jacobian matrix.
            if(! findJacobian(x, f)) {
                state = FAILED;
                break;
            }

            // Compute the QR factorization of the Jacobian.
            // Pivoting is enforced.
            QRSolver solver(jacobian, f, true);
            solver.getColNorm(D);

            // On the first iteration calculate the norm of the scaled x
            // and initialize the step bound delta.
            if(iter == 1) {
                for(int j = 0; j < D.size(); ++j) if(D(j) == 0.0) D(j) = 1.0;
                double dxnorm = scaled_norm(D, x);
                delta = dxnorm ? dxnorm : 1.0;
            }

            // Inner iteration.
            do {
                // Inner iteration: Apply Levenberg-Marquardt parameter.
                state = INTERNAL;
                Vector<double> h;
                lmpar(solver, D, f, h);
                double dhnorm = scaled_norm(D, h);

                // On the first iteration, adjust the initial step bound.
                if(iter == 1) delta = min(delta, dhnorm);

                // Evaluate the function at x + h and calculate its norm.
                Vector<double> f1;
                double fnorm1;

                if(Match::block->evaluate(x - h, f1)) {
                    fnorm1 = euclidean_norm(f1);
                } else {
                    fnorm1 = 2.0 * fnorm;
                }

                // Compute the scaled actual reduction.
                double actred = - 1.0;
                if(0.1 * fnorm1 < fnorm) {
                    actred = 1.0 - (fnorm1 * fnorm1) / (fnorm * fnorm);
                }

                // Compute the scaled predicted reduction
                // and the scaled directional derivative.
                double gnorm = euclidean_norm(jacobian * h);
                double temp1 = gnorm / fnorm;
                double temp2 = sqrt(mu) * dhnorm / fnorm;
                double prered = temp1 * temp1 + 2.0 * temp2 * temp2;
                double dirder = - (temp1 * temp1 + temp2 * temp2);

                // Compute the ratio of the actual to the predicted reduction.
                ratio = prered ? actred / prered : 0.0;

                // Update the step bound.
                if(ratio <= 0.25) {
                    double temp = (actred >= 0.0) ? 0.5 : dirder / (2.0 * dirder + actred);
                    if(0.1 * fnorm1 >= fnorm  ||  temp < 0.1) temp = 0.1;
                    delta = temp * min(delta, 10.0 * dhnorm);
                    mu /= temp;
                } else if(mu == 0.0  ||  ratio >= 0.75) {
                    delta = 2.0 * dhnorm;
                    mu = 0.5 * mu;
                }

                // Successive tests for successful iteration,
                // for convergence, and for abnormal termination.
                if(ratio >= 0.0001) {
                    x = x - h;
                    dxnorm = scaled_norm(D, x);
                    f = f1;
                    fnorm = fnorm1;
                    iter++;
                    state = PROGRESS;
                }

                if((std::abs(actred) <= tol && prered <= tol && ratio <= 2.0) ||
                   (delta <= tol * dxnorm) ||
                   (fnorm *fnorm <= tol)) {
                    state = CONVERGED;
                } else if(Match::block->getCallCount() >= nfcnmax) {
                    state = CALL_LIMIT;
                } else if((std::abs(actred) <= DBL_EPSILON && prered <= DBL_EPSILON &&
                           ratio <= 2.0) ||
                          (delta <= DBL_EPSILON * dxnorm) ||
                          (fnorm *fnorm <= DBL_EPSILON) ||
                          (gnorm <= DBL_EPSILON * fnorm)) {
                    state = ACCURACY_LIMIT;
                }
            } while(state < CONVERGED && ratio < 0.0001);
        } while(state < CONVERGED);
    }

    // Final evaluation and printout.
    Match::block->evaluate(x, f);
    Match::block->print(method, state);
}


void LMDif::lmpar(QRSolver &solver, const Array1D<double> &D,
                  const Vector<double> &/*f*/, Vector<double> &h) {
    // The iteration limit.
    static const int limit = 10;

    // Evaluate and test for acceptance.
    solver.solveR(h);
    double dhnorm = scaled_norm(D, h);
    double lm = dhnorm - delta;

    // Accept the estimate, if the function is small enough.
    if(lm <= 0.1 * delta) {
        mu = 0.0;
        return;
    }

    // Start iterations.
    int n = D.size();
    for(int iter = 0; iter < limit; ++iter) {
        if(mu <= 0.0) mu = max(DBL_EPSILON, 0.001 * delta);
        solver.solveS(D, mu, h);
        dhnorm = scaled_norm(D, h);
        lm = dhnorm - delta;

        // If the goal function lm is negative or small enough, return.
        if(lm <= 0.1 * delta) break;

        // Compute new parameter mu.
        Vector<double> temp(n);
        for(int i = 0; i < n; ++i) temp[i] = D[i] * D[i] * h[i] / dhnorm;

        solver.solveST(temp);
        double dmu = lm / (delta * (temp * temp));
        mu += dmu;

        if(mu < 0.0) {
            mu = 0.0;
            break;
        } else if(std::abs(dmu) < 1.0e-6 * mu) {
            break;
        }
    }
}


bool LMDif::findJacobian(Vector<double> &x, Vector<double> &f) {
    Vector<double> temp_x(x);
    Vector<double> temp_f(f);
    int n = x.size();
    int m = f.size();
    bool flag = true;

    for(int j = 0; j < n; ++j) {
        double step = (temp_x(j) == 0) ? fraction : fraction * std::abs(temp_x(j));
        temp_x(j) += step;
        flag = Match::block->evaluate(temp_x, temp_f);
        temp_x(j) = x(j);

        if(! flag) break;
        for(int i = 0; i < m; ++i) {
            jacobian(i, j) = (temp_f(i) - f(i)) / step;
        }
    }

    return flag;
}