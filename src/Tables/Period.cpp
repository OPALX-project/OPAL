// ------------------------------------------------------------------------
// $RCSfile: Period.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Period
//   The class for the OPAL TWISS command.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/Period.h"
#include "Algorithms/IdealMapper.h"
#include "Algorithms/Mapper.h"
#include "Attributes/Attributes.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "FixedAlgebra/FNormalForm.h"
#include "FixedAlgebra/FDynamicFP.h"
#include "FixedAlgebra/FStaticFP.h"
#include "FixedAlgebra/LinearFun.h"
#include "FixedAlgebra/LinearMap.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace std;
using std::setw;


// Class Period
// ------------------------------------------------------------------------

Period::Period():
    Twiss(SIZE, "TWISS",
          "The \"TWISS\" command defines a table of lattice functions\n"
          "which can be matched or tabulated over a periodic line.") {
    itsAttr[MICADO] = Attributes::makeReal
                      ("MICADO", "Number of iterations for MICADO algorithm", 0.0);
    itsAttr[CORRECTORS] = Attributes::makeReal
                          ("CORRECTORS", "Number of correctors for MICADO algorithm", 0.0);
    itsAttr[THREAD] = Attributes::makeString
                      ("THREAD", "Name of method for closed orbit threader");
    itsAttr[TOLQ] = Attributes::makeReal
                    ("TOLQ", "Tolerance for positions in closed orbit threader", 1.0e-3);
    itsAttr[TOLP] = Attributes::makeReal
                    ("TOLP", "Tolerance for momenta in closed orbit threader", 1.0e-3);

    // READ ONLY.
    itsAttr[CIRCUM] = Attributes::makeReal
                      ("CIRCUM", "Circumference in m");
    itsAttr[CIRCUM].setReadOnly(true);

    itsAttr[FREQ] = Attributes::makeReal
                    ("FREQ0", "Revolution frequency in Hz");
    itsAttr[FREQ].setReadOnly(true);

    itsAttr[Q1] = Attributes::makeReal
                  ("Q1", "Tune for mode 1");
    itsAttr[Q1].setReadOnly(true);

    itsAttr[Q2] = Attributes::makeReal
                  ("Q2", "Tune for mode 2");
    itsAttr[Q2].setReadOnly(true);

    itsAttr[Q3] = Attributes::makeReal
                  ("Q3", "Tune for mode 3");
    itsAttr[Q3].setReadOnly(true);

    itsAttr[U0] = Attributes::makeReal
                  ("U0", "Energy loss per turn in MeV");
    itsAttr[U0].setReadOnly(true);

    itsAttr[J1] = Attributes::makeReal
                  ("J1", "Damping partition number for mode 1");
    itsAttr[J1].setReadOnly(true);

    itsAttr[J2] = Attributes::makeReal
                  ("J2", "Damping partition number for mode 2");
    itsAttr[J2].setReadOnly(true);

    itsAttr[J3] = Attributes::makeReal
                  ("J3", "Damping partition number for mode 3");
    itsAttr[J3].setReadOnly(true);

    itsAttr[DELTAP] = Attributes::makeReal
                      ("DELTAP", "Differential momentum variation");
    itsAttr[DELTAP].setReadOnly(true);

    registerOwnership(AttributeHandler::COMMAND);
}


Period::Period(const std::string &name, Period *parent):
    Twiss(name, parent)
{}


Period::~Period()
{}


Period *Period::clone(const std::string &name) {
    return new Period(name, this);
}


void Period::fill() {
    //std::cerr << "==> In Period::fill()..." << std::endl;
    if(refill) {
        // Set truncation order.
        //DTA:    FTps<double,6>::setGlobalTruncOrder(order);

        // Search for closed orbit.
        findClosedOrbit();
        //std::cerr << "Closed orbit \"found\"" << std::endl;

        // Now map is the map around the fixed point.
        FVps<double, 6> map;
        itsMapper->getMap(map);
        orbit = map.constantTerm();
        //std::cerr << "  co = " << orbit << std::endl;
        //std::cerr << " map = " <<  map  << std::endl;
        for(int i = 0; i < 6; ++i) map[i][0] = 0.0;
        FTps<double, 6> A_lie;
        FVps<double, 6> A_scr;
        if(map[5] == FTps<double, 6>::makeVariable(5)) {
            // Map is static.
            //std::cerr << "  static map = " << map << std::endl;
            FStaticFP<6> fix(map);
            //std::cerr << "  fix(map) constructed" << std::endl;
            map = fix.getFixedPointMap();
            //std::cerr << "  got fixedPt map --> " << map << std::endl;
            const FNormalForm<6> normal(map);
            A_lie = normal.normalisingMap();
            //std::cerr << "  got Lie form of normalising map --> " << A_lie << std::endl;
            A_scr = FVps<double, 6>(normal.eigenVectors());
            //std::cerr << "  got matrix form of normalising map --> " << A_scr << std::endl;
            for(int i = order; i >= 3; --i) {
                A_scr = ExpMap(A_lie.filter(i, i), A_scr);
                //std::cerr << "  A_scr(" << i << ") --> " << A_scr << std::endl;
            }
            A_scr = fix.getFixedPoint().substitute(A_scr) + fixPoint;
            //std::cerr << "  A_scr(-1) --> " << A_scr << std::endl;
        } else {
            // Map is dynamic; fixed point is already known.
            const FNormalForm<6> normal(map);
            A_lie = normal.normalisingMap();
            A_scr = FVps<double, 6>(normal.eigenVectors());
            for(int i = order; i >= 3; --i) {
                A_scr = ExpMap(A_lie.filter(i, i), A_scr);
            }
            A_scr += fixPoint;

            if(Options::warn) {
                std::cerr << "\n### Warning ### Momentum is not constant, "
                          << "Twiss is three-dimensional.\n" << std::endl;
            }
        }

        // Initialise mapper, clear local data, and fill table.
        curly_A = A_scr.linearTerms();
        itsMapper->setMap(LinearMap<double, 6>() + orbit);
        put();

        // Fill in the read-only data.
        const Row &row = itsTable->back();
        double arc = getS(row);
        Attributes::setReal(itsAttr[CIRCUM], arc);
        Attributes::setReal(itsAttr[FREQ], Physics::c * reference->getBeta() / arc);
        Attributes::setReal(itsAttr[Q1], getMUi(row, 0));
        Attributes::setReal(itsAttr[Q2], getMUi(row, 1));
        Attributes::setReal(itsAttr[Q3], getMUi(row, 2));

        // Fill is complete.
        refill = false;
    }
    //std::cerr << "==> Leaving Period::fill()" << std::endl;
}


void Period::printTable(std::ostream &os, const CellArray &cells) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Print table header.
    printTableTitle(os, "Periodic lattice functions");

    // Print table body.
    printTableBody(os, cells);

    // Write table specific summary.
    const Row &row = itsTable->back();
    os << "Period length =  " << setw(16)
       << Attributes::getReal(itsAttr[CIRCUM])
       << "    Qx =         " << setw(16) << getMUi(row, 0)
       << "    Qy =         " << setw(16) << getMUi(row, 1)
       << '\n'
       << "DeltaP =         " << setw(16)
       << Attributes::getReal(itsAttr[DELTAP])
       << "    BetaX(max) = " << setw(16)
       << Attributes::getReal(itsAttr[BETXMAX])
       << "    BetaY(max) = " << setw(16)
       << Attributes::getReal(itsAttr[BETYMAX]) << '\n'
       << "                                 "
       << "    x(max) =     " << setw(16)
       << Attributes::getReal(itsAttr[XCMAX])
       << "    y(max) =     " << setw(16)
       << Attributes::getReal(itsAttr[YCMAX]) << '\n'
       << "                                 "
       << "    x(rms) =     " << setw(16)
       << Attributes::getReal(itsAttr[XCRMS])
       << "    y(rms) =     " << setw(16)
       << Attributes::getReal(itsAttr[YCRMS]) << '\n'
       << "                                 "
       << "    Dx(max) =    " << setw(16)
       << Attributes::getReal(itsAttr[DXMAX])
       << "    Dy(max) =    " << setw(16)
       << Attributes::getReal(itsAttr[DYMAX]) << '\n'
       << "                                 "
       << "    Dx(rms) =    " << setw(16)
       << Attributes::getReal(itsAttr[DXRMS])
       << "    Dy(rms) =    " << setw(16)
       << Attributes::getReal(itsAttr[DYRMS]) << '\n';

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


void Period::findClosedOrbit() {
    //std::cerr << "==> In Period::findClosedOrbit()" << std::endl;
    static const int iteration_limit = 20;
    static const double itsTolerance = 1.0e-8;
    static const LinearFun<double, 6> nrgy = LinearFun<double, 6>::makeVariable(5);
    double deltap = Attributes::getReal(itsAttr[DELTAP]);

    // Initialize fixPoint.
    fixPoint[0] = 0.e0;
    fixPoint[1] = 0.e0;
    fixPoint[2] = 0.e0;
    fixPoint[3] = 0.e0;
    fixPoint[4] = 0.e0;
    fixPoint[5] = deltap;
    //fixPoint[5]=0.e0;
    //std::cerr << " [findCO:] fixPoint =" << fixPoint << std::endl;

    double error = 0.0;
    for(int count = 0; count < iteration_limit; ++count) {
        //std::cerr << " [findCO:] count = " << count << std::endl;
        // Initial guess for closed orbit.
        LinearMap<double, 6> identity;
        LinearMap<double, 6> mapAtEnd;

        // Compute the one-turn map around the closed orbit.
        //std::cerr << " [findCO:] computing map about fixPoint = " << fixPoint << std::endl;
        itsMapper->setMap(identity + fixPoint);
        itsMapper->execute();
        itsMapper->getMap(mapAtEnd);
        //std::cerr << " [findCO:] have map about fixPoint" << std::endl;

        // Get system of equations for fixed point.
        FMatrix<double, 6, 6> A   = mapAtEnd.linearTerms();
        FVector<double, 6> Error = mapAtEnd.constantTerm() - fixPoint;
        double errold = error;
        error = 0.0;

        //std::cerr << " [findCO:] count " << count << ": Error =\n {";
        //for(int i=0;i<6;i++) std::cerr << " " << Error(i);
        //std::cerr << " }" << std::endl;
        //std::cerr << " [findCO:] A =\n" << A << std::endl;

        //std::cerr << " [findCO:] finding fixed point ";
        if(mapAtEnd[5] == nrgy + deltap) {
            // Finding static fixed point.
            //std::cerr << "for static map ...\n" << std::endl;
            for(int i = 0; i < 4; i++) {
                A(i, i) -= 1.0;
                if(abs(Error(i)) > error) error = abs(Error(i));
            }
            for(int i = 4; i < 6; i++) {
                for(int j = 0; j < 6; j++) A(i, j) = A(j, i) = 0.0;
                A(i, i) = 1.0;
                Error(i) = 0.0;
            }
        } else {
            // Finding dynamic fixed point.
            //std::cerr << "for dynamic map ...\n" << std::endl;
            for(int i = 0; i < 6; i++) {
                A(i, i) -= 1.0;
                if(abs(Error(i)) > error) error = abs(Error(i));
            }
        }

        //std::cerr << " [findCO:] count " << count << ": Error =\n {";
        //for(int i=0;i<6;i++) std::cerr << " " << Error(i);
        //std::cerr << " }" << std::endl;
        //std::cerr << " [findCO:] A =\n" << A << std::endl;

        // Correction for fixed point.
        FLUMatrix<double, 6> lu(A);
        //std::cerr << " [findCO:] have lu(A)" << std::endl;
        lu.backSubstitute(Error);
        //std::cerr << " [findCO:] backSub done" << std::endl;
        fixPoint -= Error;
        //std::cerr << " [findCO:] count errold error" << count << errold << error << std::endl;
        // return if the error vanishes or the error has fallen below some crude
        // tolerance (machineEps^(1/2)) and bounced (because of round-off error)
        if(count && (error == 0.0 || (error < itsTolerance && error >= errold))) break;
        //if (error<itsTolerance) break;
        //std::cerr << " [findCO:] again" << std::endl;
    }
    //std::cerr << " [findCO:]  fixPoint = " << fixPoint << std::endl;
    //std::cerr << "==> Leaving Period::findClosedOrbit()" << std::endl;
}