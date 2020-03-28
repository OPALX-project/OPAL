// ------------------------------------------------------------------------
// $RCSfile: Simplex.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Simplex
//   This class encapsulates a minimisation according to the SIMPLEX method
//   taken from the MINUIT package.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:10:01 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Match/Simplex.h"
#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Attributes/Attributes.h"
#include "Match/Match.h"
#include "Utilities/OpalException.h"
#include <cfloat>
#include <cmath>
#include <numeric>


// Class Simplex
// ------------------------------------------------------------------------

namespace {
    // The attributes of class Simplex.
    enum {
        TOLERANCE,   // The desired tolerance.
        CALLS,       // The maximum number of calls to the matching functions.
        SIZE
    };

    // Auxiliary functions for returning the minimum or maximum.
    inline double minfun(double x, double y)
    { return x < y ? x : y; }

    inline double maxfun(double x, double y)
    { return x > y ? x : y; }
}


Simplex::Simplex():
    Action(SIZE, "SIMPLEX",
           "The \"SIMPLEX\" sub-command adjusts parameters according to the "
           "\"SIMPLEX\" method taken from MINUIT.") {
    itsAttr[TOLERANCE] = Attributes::makeReal
                         ("TOLERANCE", "The desired tolerance", 1.0e-6);
    itsAttr[CALLS] = Attributes::makeReal
                     ("CALLS", "Maximum number of calls to the matching functions", 1000.);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


Simplex::Simplex(const std::string &name, Simplex *parent):
    Action(name, parent)
{}


Simplex::~Simplex()
{}


Simplex *Simplex::clone(const std::string &name) {
    return new Simplex(name, this);
}


void Simplex::execute() {
    // Built-in constants.
    static const char method[] = "SIMPLEX";
    static const double alpha  = 1.0;
    static const double beta   = 0.5;
    static const double gamma  = 2.0;
    static const double rhomin = 4.0;
    static const double rhomax = 8.0;
    static const double rho1   = 1.0 + alpha;
    static const double rho2   = rho1 + alpha * gamma;

    // Working storage.

    /*
    ada 15-6-2000

    double fnorm;

    else if ((fnorm = euclidean_norm(F)) <= tol) {

    */



    int n = Match::block->countVariables();
    int m = Match::block->countFunctions();
    Vector<double> X(n);
    Vector<double> F(m);
    double tol = std::max(Attributes::getReal(itsAttr[TOLERANCE]), DBL_EPSILON);
    MatchState state = START;

    // Initial function evaluation.
    Match::block->getVariables(X);
    if(! Match::block->evaluate(X, F)) {
        state = FAILED;
    } else if(euclidean_norm(F) <= tol) {
        state = CONVERGED;
    } else {
        // Fetch command attributes.
        int nfcnmax = int(std::round(Attributes::getReal(itsAttr[CALLS])));

        // Set up the simplex.
        Fsim = Array1D<double>(n + 1);
        Xsim = Array1D< Vector<double> >(n + 1);

        Vector<double> dX(0.001 * X);
        double Fmin = F * F;
        state = START;

        // Outer loop: start/restart algorithm.
        do {
            // Keep initial point in Xbar/Fbar.
            Vector<double> Xbar = X;
            double Fbar = Fmin;

            for(int i = 0; i < n; ++i) {
                double Xsave = X(i);
                double Fsave = Fmin;
                double step  = dX(i);

                // Find proper initial direction and step.
                int idir;
                for(idir = 1; idir <= 12; ++idir) {
                    X(i) = Xsave + step;
                    if(Match::block->evaluate(X, F)  && (Fmin = F * F) <= Fsave) break;
                    if(idir % 2 == 0) step *= 0.1;
                    step = - step;
                }

                if(idir <= 12) {
                    // Improvement found; try increasing steps.
                    for(int ns = 1; ns <= 3; ++ns) {
                        Xsave = X(i);
                        Fsave = Fmin;
                        step *= 3.0;
                        X(i) += step;

                        // If new point is not an improvement,
                        // restore previous best point and quit.
                        if(! Match::block->evaluate(X, F)  || (Fmin = F * F) > Fsave) {
                            X(i) = Xsave;
                            Fmin = Fsave;
                            break;
                        }
                    }
                }

                // Store local minimum in i'th direction as vertex i.
                Fsim(i) = Fmin;
                Xsim(i) = X;
            }

            // Store initial point as vertex n.
            jh = jl = n;
            razzia(Fbar, Xbar);

            // Extract best point of simplex.
            X = Xsim(jl);
            Fmin = Fsim(jl);
            Match::block->setVariables(X);

            // Inner loop: perform algorithm on simplex.
            bool restart = false;
            while(Fmin > tol) {
                // Iteration print-out.
                Match::block->print(method, state);

                // Test for call limit.
                if(Match::block->getCallCount() >= nfcnmax) {
                    state = CALL_LIMIT;
                } else {
                    // Calculate Xbar and X*.
                    Xbar = std::accumulate
                           (Xsim.begin(), Xsim.end(), - Xsim(jh)) / double(n);
                    Vector<double> Xstar = Xbar + alpha * (Xbar - Xsim(jh));
                    double Fstar = 2.0 * Fsim(jh);
                    int jhold = jh;

                    if(Match::block->evaluate(Xstar, F)) {
                        if((Fstar = F * F) < Fsim(jl)) {
                            // Point X* is better than previous best point;
                            // try expanded point X**.
                            Vector<double> Xstst = Xbar + gamma * (Xstar - Xbar);
                            double Fstst = 0.0;
                            double rho = 0.0;

                            if(Match::block->evaluate(Xstst, F)) {
                                // Fit a parabola through Fsim(jh), f*, f**; minimum = rho.
                                Fstst = F * F;
                                double F1 = (Fstar - Fsim(jh)) * rho2;
                                double F2 = (Fstst - Fsim(jh)) * rho1;
                                rho = 0.5 * (rho2 * F1 - rho1 * F2) / (F1 - F2);
                            }

                            if(rho < rhomin) {
                                // Minimum inexistent or too close to pbar;
                                // Use X** if it gives improvement; otherwise use X*.
                                if(Fstst < Fsim(jl)) {
                                    razzia(Fstst, Xstst);
                                } else {
                                    razzia(Fstar, Xstar);
                                }
                            } else {
                                // Usable minimum found.
                                if(rho > rhomax) rho = rhomax;

                                Vector<double> Xrho = Xsim(jh) + rho * (Xbar - Xsim(jh));
                                double Frho;

                                // Select farthest point which gives decent improvement.
                                if(Match::block->evaluate(Xrho, F)  &&
                                   (Frho = F * F) < Fsim(jl)  &&
                                   Frho < Fstst) {
                                    razzia(Frho, Xrho);
                                } else if(Fstst < Fsim(jl)) {
                                    razzia(Fstst, Xstst);
                                } else {
                                    razzia(Fstar, Xstar);
                                }
                            }
                        } else {
                            // F* is higher than Fsim(jl).
                            if(Fstar < Fsim(jh)) razzia(Fstar, Xstar);

                            // If F* is still highest value, try contraction,
                            // giving point X**.
                            if(jhold == jh) {
                                Vector<double> Xstst = Xbar + beta * (Xsim(jh) - Xbar);
                                double Fstst;

                                if(Match::block->evaluate(Xstst, F)  &&
                                   (Fstst = F * F) < Fsim(jh)) {
                                    razzia(Fstst, Xstst);
                                } else {
                                    restart = true;
                                }
                            }
                        }
                    }

                    // Extract best point as new minimum.
                    Fmin = Fsim(jl);
                    X    = Xsim(jl);
                    Match::block->setVariables(X);
                }

                if(restart) break;
            }

            // Test for convergence.
            if(restart) {
                // Main loop must be restarted.
                state = (state == RESTART) ? FAILED : RESTART;
                restart = false;
            } else {
                // Recompute step sizes.
                if(Fmin >= 2.0 * (tol + DBL_EPSILON)) {
                    Vector<double> Xmin(Xsim(0));
                    Vector<double> Xmax(Xsim(0));

                    for(int j = 1; j < n; ++j) {
                        std::transform(Xmin.begin(), Xmin.end(), Xsim(j).begin(),
                                       Xmin.begin(), minfun);
                        std::transform(Xmax.begin(), Xmax.end(), Xsim(j).begin(),
                                       Xmax.begin(), maxfun);
                    }

                    dX = Xmax - Xmin;
                } else {
                    state = CONVERGED;
                }
            }
        } while(state < CONVERGED);
    }

    // Final print-out.
    Match::block->print(method, state);
}


void Simplex::razzia(double Fnew, Vector<double> &Xnew) {
    // Replace vertex with highest function value.
    Xsim(jh) = Xnew;
    Fsim(jh) = Fnew;

    // Find indices of lowest and highest function value.
    jl = 0;
    jh = 0;
    int n = Match::block->countVariables();
    for(int j = 1; j <= n; ++j) {
        if(Fsim(j) < Fsim(jl)) jl = j;
        if(Fsim(j) > Fsim(jh)) jh = j;
    }
}