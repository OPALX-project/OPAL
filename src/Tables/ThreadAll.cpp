// ------------------------------------------------------------------------
// $RCSfile: ThreadAll.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ThreadAll
//   The abstract class ThreadAll implements the interface for a table buffer
//   holding lattice function.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/ThreadAll.h"
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
#include <iomanip>
#include <iostream>


// Class ThreadAll
// ------------------------------------------------------------------------

ThreadAll::ThreadAll():
    CorrectionBase(SIZE, "THREADALL",
                   "The \"THREADALL\" command threads the closed orbit "
                   "using the position and angle at all positions.") {
    itsAttr[TOLQ] = Attributes::makeReal
                    ("TOLQ", "The tolerance for the closed orbit positions.");
    itsAttr[TOLP] = Attributes::makeReal
                    ("TOLP", "The tolerance for the closed orbit angles.");
    itsAttr[LISTM] = Attributes::makeBool
                     ("LISTM", "List the monitors after correction");
    itsAttr[LISTC] = Attributes::makeBool
                     ("LISTC", "List the correctors after correction");

    registerOwnership(AttributeHandler::COMMAND);
}


ThreadAll::ThreadAll(const std::string &name, ThreadAll *parent):
    CorrectionBase(name, parent)
{}


ThreadAll::~ThreadAll()
{}


ThreadAll *ThreadAll::clone(const std::string &name) {
    return new ThreadAll(name, this);
}


void ThreadAll::execute() {
    // Find Table definition.
    const std::string &lineName = Attributes::getString(itsAttr[LINE]);
    BeamSequence *use = BeamSequence::find(lineName);

    // Find Beam data.
    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    Beam *beam = Beam::find(beamName);
    reference = beam->getReference();

    // Get the data for correction.
    RangeRep range = Attributes::getRange(itsAttr[RANGE]);
    double tolp = Attributes::getReal(itsAttr[TOLP]);
    double tolq = Attributes::getReal(itsAttr[TOLQ]);
    bool listc = Attributes::getBool(itsAttr[LISTC]);
    bool listm = Attributes::getBool(itsAttr[LISTM]);

    // Make sure all is up-to-date.
    OpalData::getInstance()->update();

    // Create flat list with space for data storage.
    Flatten<Row> flattener(*use->fetchLine(), itsLine, range);
    flattener.execute();

    if(itsLine.empty() && Options::warn) {
        std::cerr << "\n### Warning ### Lattice function table \""
                  << lineName << "\" contains no elements.\n" << std::endl;
        return;
    }

    // Set up the tables of correctors and monitors.
    setupTables();

    // Start point for threading.
    itsTracker = new OrbitTracker(itsLine, reference, false, false);
    itsTracker->setOrbit(orbitGuess);
    static const int iteration_limit = 100;
    int count = 0;
    TLine::iterator iter = itsLine.begin();

    while(iter != itsLine.end()) {
        // Advance through element.
        try {
            iter->accept(*itsTracker);

            // Test for overflow in current element.
            FVector<double, 6> orbit = itsTracker->getOrbit();
            iter->orbit = orbit;

            if(! iter->isUsed[0] &&
               (std::abs(orbit[0]) > tolq || std::abs(orbit[1]) > tolp)) {
                // Overflow detected in x.
                if(count++ > iteration_limit) break;
                correct(0, iter);
            } else if(! iter->isUsed[1] &&
                      (std::abs(orbit[2]) > tolq || std::abs(orbit[3]) > tolp)) {
                // Overflow detected in y.
                if(count++ > iteration_limit) break;
                correct(1, iter);
            } else {
                ++iter;
            }
        } catch(const DomainError &) {
            // Domain error detected;
            // attempt to correct in the plane which has the larger deviation.
            FVector<double, 6> orbit = itsTracker->getOrbit();
            iter->orbit = orbit;

            if(std::abs(orbit[1]) >= std::abs(orbit[3])) {
                correct(0, iter);
            } else {
                correct(1, iter);
            }
        }
    }

    std::cout << "After threading with \"THREADALL\":\n";
    for(int plane = 0; plane < 2; ++plane) {
        listCorrectors(listc, plane);
        listMonitors(listm, plane);
    }
}


void ThreadAll::correct(int plane, TLine::iterator &iter) {
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
        iter->isUsed[plane] = true;
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