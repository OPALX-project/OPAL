// ------------------------------------------------------------------------
// $RCSfile: Micado.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Micado
//   The abstract class Micado implements the interface for a table buffer
//   holding lattice function.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/Micado.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/RangeRep.h"
#include "Algebra/Matrix.h"
#include "Algebra/Vector.h"
#include "Algorithms/ThickMapper.h"
#include "Algorithms/LinearMapper.h"
#include "Algorithms/ThinMapper.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Tables/Flatten.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


// Class Micado
// ------------------------------------------------------------------------

Micado::Micado():
    CorrectionBase(SIZE, "MICADO",
                   "The \"MICADO\" command makes a correction "
                   "using the \"MICADO\" algorithm.") {
    itsAttr[METHOD] = Attributes::makeString
                      ("METHOD", "The method used for orbit tracking", "LINEAR");
    itsAttr[TOL] = Attributes::makeReal
                   ("TOL", "The tolerance for the r.m.s closed orbit.");
    itsAttr[ITERATIONS] = Attributes::makeReal
                          ("ITERATIONS", "The number of iterations");
    itsAttr[CORRECTORS] = Attributes::makeReal
                          ("CORRECTORS", "The number of correctors to be used");
    itsAttr[PLANE] = Attributes::makeString
                     ("PLANE", "The plane(s) for which a Correction is desired", "BOTH");
    itsAttr[LISTC1] = Attributes::makeBool
                      ("LISTC1", "List the correctors before correction");
    itsAttr[LISTC] = Attributes::makeBool
                     ("LISTC", "List the correctors during correction");
    itsAttr[LISTC2] = Attributes::makeBool
                      ("LISTC2", "List the correctors after correction");
    itsAttr[LISTM1] = Attributes::makeBool
                      ("LISTM1", "List the monitors before correction");
    itsAttr[LISTM2] = Attributes::makeBool
                      ("LISTM2", "List the monitors after correction");

    registerOwnership(AttributeHandler::COMMAND);
}


Micado::Micado(const std::string &name, Micado *parent):
    CorrectionBase(name, parent)
{}


Micado::~Micado()
{}


Micado *Micado::clone(const std::string &name) {
    return new Micado(name, this);
}


void Micado::execute() {
    // Find Table definition.
    const std::string &lineName = Attributes::getString(itsAttr[LINE]);
    BeamSequence *use = BeamSequence::find(lineName);

    // Find Beam data.
    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    Beam *beam = Beam::find(beamName);
    reference = beam->getReference();

    // Get the data for correction.
    RangeRep range = Attributes::getRange(itsAttr[RANGE]);
    int iterations = int(std::round(Attributes::getReal(itsAttr[ITERATIONS])));
    bool listm1 = Attributes::getBool(itsAttr[LISTM1]);
    bool listm2 = Attributes::getBool(itsAttr[LISTM2]);
    bool listc1 = Attributes::getBool(itsAttr[LISTC1]);
    bool listc2 = Attributes::getBool(itsAttr[LISTC2]);

    // Make sure all is up-to-date.
    OpalData::getInstance()->update();

    // Create flat line for correction.
    Flatten<Row> flattener(*use->fetchLine(), itsLine, range);
    flattener.execute();

    if(itsLine.empty() && Options::warn) {
        std::cerr << "\n### Warning ### \"MICADO\" table \""
                  << lineName << "\" contains no elements.\n" << std::endl;
    }

    // Decode the plane name.
    const std::string &planeName = Attributes::getString(itsAttr[PLANE]);
    bool planes[2] = { false, false };
    if(planeName == "BOTH") {
        planes[0] = planes[1] = true;
    } else if(planeName == "X") {
        planes[0] = true;
    } else if(planeName == "Y") {
        planes[1] = true;
    } else {
        throw OpalException("Micado::execute()",
                            "Plane name \"" + planeName + "\" is unknown.");
    }

    // Decode the method name.
    const std::string &methodName = Attributes::getString(itsAttr[METHOD]);
    // Create the mapper to be used.
    FTps<double, 6>::setGlobalTruncOrder(2);
    if(methodName == "THICK") {
        itsMapper = new ThickMapper(itsLine, reference, false, false);
    } else if(methodName == "THIN") {
        itsMapper = new ThinMapper(itsLine, reference, false, false);
    } else if(methodName == "LINEAR") {
        itsMapper = new LinearMapper(itsLine, reference, false, false);
    } else {
        throw OpalException("Micado::execute()",
                            "Method name \"" + methodName + "\" is unknown.");
    }

    // Set up the tables of correctors and monitors.
    setupTables();

    // Perform correction.
    for(int iteration = 1; iteration <= iterations; ++iteration) {
        std::cout << "\n\"MICADO\" iteration " << iteration << ":\n";

        for(int plane = 0; plane < 2; ++plane) {
            if(planes[plane]) {
                // Store the current closed orbit.
                findClosedOrbit();
                listCorrectors(listc1, plane);
                listMonitors(listm1, plane);

                // The A matrix.
                // Each row corresponds to one monitor reading,
                // Each column corresponds to one corrector setting.
                int numberMonitors   = monitorTable[plane].size();
                int numberCorrectors = correctorTable[plane].size();
                Matrix<double> A(numberMonitors, numberCorrectors);
                setupInfluence(plane, A);

                // The monitor readings.
                Vector<double> B(numberMonitors);
                setupReadings(plane, B);

                // Solve the system of equations.
                solve(plane, A, B);
            }
        }
    }

    if(listc2 || listm2) {
        std::cout << "\"MICADO\" finished:\n";
        for(int plane = 0; plane < 2; ++plane) {
            listCorrectors(listc2, plane);
            listMonitors(listm2, plane);
        }
    }
}


void Micado::findClosedOrbit() {
    static const int iteration_limit = 20;
    static const double itsTolerance = 1.0e-8;

    for(int count = 0; count < iteration_limit; ++count) {
        // Initial guess for closed orbit.
        LinearMap<double, 6> identity;
        LinearMap<double, 6> currentMap;

        // Compute the one-turn map around the closed orbit.
        itsMapper->setMap(identity + orbitGuess);
        for(TLine::iterator iter = itsLine.begin();
            iter != itsLine.end(); ++iter) {
            iter->accept(*itsMapper);
            itsMapper->getMap(currentMap);
            iter->orbit = currentMap.constantTerm();
        }

        // Get system of equations for fixed point.
        FMatrix<double, 6, 6> A   = currentMap.linearTerms();
        FVector<double, 6> Error = currentMap.constantTerm() - orbitGuess;
        double error = 0.0;

        if(currentMap[5] == LinearFun<double, 6>::makeVariable(5)) {
            // Finding static fixed point.
            for(int i = 0; i < 4; i++) {
                A(i, i) -= 1.0;
                if(std::abs(Error(i)) > error) error = std::abs(Error(i));
            }

            for(int i = 4; i < 6; i++) {
                for(int j = 0; j < 6; j++) A(i, j) = A(j, i) = 0.0;
                A(i, i) = 1.0;
                Error(i) = 0.0;
            }
        } else {
            // Finding dynamic fixed point.
            for(int i = 0; i < 6; i++) {
                A(i, i) -= 1.0;
                if(std::abs(Error(i)) > error) error = std::abs(Error(i));
            }
        }

        // Correction for fixed point.
        FLUMatrix<double, 6> lu(A);
        lu.backSubstitute(Error);
        orbitGuess -= Error;
        if(error < itsTolerance) break;
    }
}


void Micado::setupInfluence(int plane, Matrix<double> &A) {
    // Prepare one turn matrix.
    FMatrix<double, 6, 6> tmat = itsLine.back().matrix;

    // Invert the 4 by 4 block for transverse motion.
    FMatrix<double, 4, 4> R;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            R(i, j) = tmat(i, j);
        }
    }
    FMatrix<double, 4, 4> R1(R - 1.0);
    R1 = - FLUMatrix<double, 4>(R1).inverse();
    FMatrix<double, 6, 6> tmat1;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            tmat1(i, j) = R1(i, j);
        }
    }
    tmat1(4, 4) = tmat1(5, 5) = 1.0;

    const int q = 2 * plane;

    // Loop over the correctors.
    LocalIter cBegin = correctorTable[plane].begin();
    LocalIter cEnd   = correctorTable[plane].end();
    LocalIter mBegin = monitorTable[plane].begin();
    LocalIter mEnd   = monitorTable[plane].end();

    int corr = 0;
    for(LocalIter cIter = cBegin; cIter != cEnd; ++cIter) {
        // Transform kick to begin of line.
        FVector<double, 6> initialOrbit;
        FMatrix<double, 6, 6> cmat = (*cIter)->matrix;

        for(int j = 0; j < 6; j += 2) {
            initialOrbit[j]   = - cmat(q, j + 1);
            initialOrbit[j+1] =   cmat(q, j);
        }

        // Loop over the monitors.
        int moni = 0;
        for(LocalIter mIter = mBegin; mIter != mEnd; ++mIter) {
            FVector<double, 6> orbit;

            if((*mIter)->arc >= (*cIter)->arc) {
                orbit = (*mIter)->matrix * (tmat1 * initialOrbit);
            } else {
                orbit = (*mIter)->matrix * (tmat1 * (tmat * initialOrbit));
            }

            A(moni, corr) = orbit[q];
            ++moni;
        }

        ++corr;
    }
}


void Micado::setupReadings(int plane, Vector<double> &B) {
    const int q = 2 * plane;
    int moni = 0;

    for(LocalIter mIter = monitorTable[plane].begin();
        mIter != monitorTable[plane].end(); ++mIter) {
        FVector<double, 6> orbit = (*mIter)->orbit;
        B(moni) = orbit[q];
        ++moni;
    }
}


void Micado::solve(int plane, Matrix<double> &A, Vector<double> &B) {
    // Set r.m.s. to a huge value in case of failure.
    bool list = Attributes::getBool(itsAttr[LISTC]);
    int m = A.nrows();
    int n = A.ncols();
    Vector<double> x(n);
    Vector<double> r(m);
    Vector<double> sqr(n);
    Vector<double> dot(n);
    std::vector<Row *>
        index(correctorTable[plane].begin(), correctorTable[plane].end());

    // Find scalar products sqr[k] = A[k].A[k] and dot[k] = A[k].b.
    double sum = 0.0;
    for(int k = 0; k < n; ++k) {
        double hh = 0.0;
        double gg = 0.0;
        for(int i = 0; i < m; ++i) {
            hh += A(i, k) * A(i, k);
            gg += A(i, k) * B[i];
        }
        sum += sum;
        sqr[k] = hh;
        dot[k] = gg;
    }
    double sqrmin = 1.e-8 * sum / double(n);

    // Begin of iteration loop.
    double tol = Attributes::getReal(itsAttr[TOL]);
    int correctors = int(std::round(Attributes::getReal(itsAttr[CORRECTORS])));
    if(correctors > n  ||  correctors == 0)  correctors = n;
    int corr = 0;
    while(true) {
        int k = corr;

        // Search the columns not yet used for largest scaled change vector.
        double maxChange  = 0.0;
        int changeIndex = -1;
        for(int j = k; j < n; ++j) {
            if(sqr[j] > sqrmin) {
                double change = dot[j] * dot[j] / sqr[j];
                if(change > maxChange) {
                    changeIndex = j;
                    maxChange  = change;
                }
            }
        }

        // Stop iterations, if no suitable column found.
        if(changeIndex < 0) break;

        // Move the column just found to next position.
        if(changeIndex > k) {
            std::swap(sqr[k], sqr[changeIndex]);
            std::swap(dot[k], dot[changeIndex]);
            std::swap(index[k], index[changeIndex]);
            A.swapColumns(k, changeIndex);
        }

        // Find beta, sigma, and vector u[k].
        double hh = 0.0;
        for(int i = k; i < m; ++i) hh += A(i, k) * A(i, k);
        double sigma = (A(k, k) > 0.0) ? sqrt(hh) : (- sqrt(hh));
        sqr[k] = - sigma;
        A(k, k) = A(k, k) + sigma;
        double beta = 1.0 / (A(k, k) * sigma);

        // Transform remaining columns of A.
        for(int j = k + 1; j < n; ++j) {
            double hh = 0.0;
            for(int i = k; i < m; ++i) hh += A(i, k) * A(i, j);
            hh *= beta;
            for(int i = k; i < m; ++i) A(i, j) -= A(i, k) * hh;
        }

        // Transform vector b.
        hh = 0.0;
        for(int i = k; i < m; ++i) hh += A(i, k) * B[i];
        hh *= beta;
        for(int i = k; i < m; ++i) B[i] -= A(i, k) * hh;

        // Update scalar products sqr[j]=A[j]*A[j] and dot[j]=A[j]*b.
        for(int j = k + 1; j < n; ++j) {
            sqr[j] -= A(k, j) * A(k, j);
            dot[j] -= A(k, j) * B[k];
        }

        // Recalculate solution vector x.
        x[k] = B[k] / sqr[k];
        for(int i = k - 1; i >= 0; --i) {
            x[i] = B[i];
            for(int j = i + 1; j < k; ++j) x[i] -= A(i, j) * x[j];
            x[i] /= sqr[i];
        }

        // Find original residual vector by backward transformation.
        for(int i = 0; i < m; ++i) r[i] = B[i];
        for(int j = k; j >= 0; --j) {
            r[j] = 0.0;
            double hh = 0.0;
            for(int i = j; i < m; ++i) hh += A(i, j) * r[i];
            hh /= sqr[j] * A(j, j);
            for(int i = j; i < m; ++i) {
                r[i] += A(i, j) * hh;
            }
        }

        // Check for convergence.
        hh = r(0) * r(0);
        for(int i = 1; i < m; ++i) hh += r[i] * r[i];
        hh = sqrt(hh / double(std::max(m, 1)));
        ++corr;

        // Give listing of correctors for this iteration.
        if(list) {
            std::cout << "\nCorrector              increment              r.m.s.";
            for(int k = 0; k < corr; ++k) {
                const std::string &name = index[k]->getElement()->getName();
                std::cout << '\n'   << name << std::string(20 - name.length(), ' ')
                          << std::setw(12) << (- x[k]);
            }
            std::cout << std::setw(20) << hh << std::endl;
        }

        // Check for convergence.
        if(hh < tol) break;
        if(corr > correctors) break;
    }

    // End of iteration loop.  Assign corrector strengths.
    for(int k = 0; k < corr; ++k) {
        addKick(plane, *index[k], - x[k]);
    }
}