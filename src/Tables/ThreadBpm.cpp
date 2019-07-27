// ------------------------------------------------------------------------
// $RCSfile: ThreadBpm.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThreadBpm
//   The abstract class ThreadBpm implements the interface for a table buffer
//   holding lattice function.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/ThreadBpm.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/RangeRep.h"
#include "Algorithms/ThickMapper.h"
#include "Algorithms/LinearMapper.h"
#include "Algorithms/ThinMapper.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Tables/Flatten.h"
#include "Utilities/DomainError.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"
#include <iomanip>
#include <iostream>


// Class ThreadBpm
// ------------------------------------------------------------------------

ThreadBpm::ThreadBpm():
    CorrectionBase(SIZE, "THREADBPM",
                   "The \"THREADBPM\" command threads the closed orbit "
                   "using the positions at all orbit monitors.") {
    itsAttr[TOL] = Attributes::makeReal
                   ("TOL", "The tolerance for the closed orbit positions.");
    itsAttr[LISTC] = Attributes::makeBool
                     ("LISTC", "List the correctors after correction");
    itsAttr[LISTM] = Attributes::makeBool
                     ("LISTM", "List the monitors after correction");

    registerOwnership(AttributeHandler::COMMAND);
}


ThreadBpm::ThreadBpm(const std::string &name, ThreadBpm *parent):
    CorrectionBase(name, parent)
{}


ThreadBpm::~ThreadBpm()
{}


ThreadBpm *ThreadBpm::clone(const std::string &name) {
    return new ThreadBpm(name, this);
}


void ThreadBpm::execute() {
    // Find Table definition.
    const std::string &lineName = Attributes::getString(itsAttr[LINE]);
    BeamSequence *use = BeamSequence::find(lineName);

    // Find Beam data.
    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    Beam *beam = Beam::find(beamName);
    reference = beam->getReference();

    // Get the data for threading.
    RangeRep range = Attributes::getRange(itsAttr[RANGE]);
    double tol = Attributes::getReal(itsAttr[TOL]);
    bool listc = Attributes::getBool(itsAttr[LISTC]);
    bool listm = Attributes::getBool(itsAttr[LISTM]);

    // Make sure all is up-to-date.
    OpalData::getInstance()->update();

    // Create flat line for correction.
    Flatten<Row> flattener(*use->fetchLine(), itsLine, range);
    flattener.execute();

    if(itsLine.empty() && Options::warn) {
        std::cerr << "\n### Warning ### Lattice function table \""
                  << lineName << "\" contains no elements.\n" << std::endl;
    }

    // Starting point for threading.
    const int iteration_limit = 100;
    int count = 0;

    // Set up the tables of correctors and monitors.
    setupTables();

    // Orbit at start position.
    itsTracker = new OrbitTracker(itsLine, reference, false, false);
    itsTracker->setOrbit(orbitGuess);
    TLine::iterator iter = itsLine.begin();

    while(iter != itsLine.end()) {
        // Store orbit before current element.
        FVector<double, 6> orbit = itsTracker->getOrbit();
        iter->orbit = orbit;

        // Advance through current element and determine element type.
        try {
            iter->accept(*itsTracker);
            orbit = itsTracker->getOrbit();
            test(iter->getElement());

            if(isMoni[0] && !iter->isUsed[0] && std::abs(orbit[0]) > tol) {
                // Overflow detected in x.
                if(++count > iteration_limit) break;
                iter->isUsed[0] = true;
                correct(0, iter);
            } else if(isMoni[1] && !iter->isUsed[1] && std::abs(orbit[2]) > tol) {
                // Overflow detected in y.
                if(++count > iteration_limit) break;
                iter->isUsed[1] = true;
                correct(1, iter);
            } else {
                // Go to next element.
                ++iter;
            }
        } catch(const DomainError &) {
            break;
        }
    }

    std::cout << "After threading with \"THREADBPM\":\n";
    for(int plane = 0; plane < 2; ++plane) {
        listCorrectors(listc, plane);
        listMonitors(listm, plane);
    }
}


void ThreadBpm::correct(int plane, TLine::iterator &iter) {
    // Common information.
    int p1 = 2 * plane;
    int p2 = 2 * plane + 1;

    // Extract orbit and transfer matrix at current position.
    FVector<double, 6>   orbit = iter->orbit;
    FMatrix<double, 6, 6> omat  = iter->matrix;

    // Backtrack to find two correctors for the given plane.
    TLine::iterator corr[2];
    FMatrix<double, 6, 6> cmat[2];
    int count = 2;

    while(iter != itsLine.begin()) {
        // Backtrack and test for corrector.
        --iter;
        test(iter->getElement());

        if(isCorr[plane]) {
            // Extract transfer matrix to this corrector.
            if(count > 0) --count;
            corr[count] = iter;
            cmat[count] = iter->matrix;

            if(count == 0) {
                // We have found two correctors.
                double r1 = orbit(p1);
                double r2 = orbit(p2);

                double a11 =
                    omat(p1, 1) * cmat[0](p1, 0) - omat(p1, 0) * cmat[0](p1, 1) +
                    omat(p1, 3) * cmat[0](p1, 2) - omat(p1, 2) * cmat[0](p1, 3) +
                    omat(p1, 5) * cmat[0](p1, 4) - omat(p1, 4) * cmat[0](p1, 5);
                double a12 =
                    omat(p1, 1) * cmat[1](p1, 0) - omat(p1, 0) * cmat[1](p1, 1) +
                    omat(p1, 3) * cmat[1](p1, 2) - omat(p1, 2) * cmat[1](p1, 3) +
                    omat(p1, 5) * cmat[1](p1, 4) - omat(p1, 4) * cmat[1](p1, 5);
                double a21 =
                    omat(p2, 1) * cmat[0](p1, 0) - omat(p2, 0) * cmat[0](p1, 1) +
                    omat(p2, 3) * cmat[0](p1, 2) - omat(p2, 2) * cmat[0](p1, 3) +
                    omat(p2, 5) * cmat[0](p1, 4) - omat(p2, 4) * cmat[0](p1, 5);
                double a22 =
                    omat(p2, 1) * cmat[1](p1, 0) - omat(p2, 0) * cmat[1](p1, 1) +
                    omat(p2, 3) * cmat[1](p1, 2) - omat(p2, 2) * cmat[1](p1, 3) +
                    omat(p2, 5) * cmat[1](p1, 4) - omat(p2, 4) * cmat[1](p1, 5);
                double det = a11 * a22 - a12 * a21;

                if(std::abs(det) > 0.01) {
                    double k1 = - (a22 * r1 - a12 * r2) / det;
                    double k2 = - (a11 * r2 - a21 * r1) / det;

                    // Store the kicks.
                    addKick(plane, *corr[0], k1);
                    addKick(plane, *corr[1], k2);

                    // Reset to start before the first kicker.
                    --iter;
                    itsTracker->setOrbit(iter->orbit);
                    return;
                }
            }
        }
    }

    // No corrector pair found; change initial conditions.
    FLUMatrix<double, 6>(omat).backSubstitute(orbit);
    orbitGuess -= orbit;
    itsTracker->setOrbit(orbitGuess);
}