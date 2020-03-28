// ------------------------------------------------------------------------
// $RCSfile: Dynamic.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Dynamic
//   The class for the OPAL DYNAMIC command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:19:44 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "PhysicsActions/Dynamic.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/LieMapper.h"
#include "Attributes/Attributes.h"
#include "FixedAlgebra/DragtFinnMap.h"
#include "FixedAlgebra/DragtFinnNormalForm.h"
#include "FixedAlgebra/FDynamicFP.h"
#include "FixedAlgebra/FTps.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;
using std::complex;


// Class Dynamic
// ------------------------------------------------------------------------

// The attributes of class Dynamic.
namespace {
    enum {
        LINE,         // The lattice to be used.
        BEAM,         // The beam to be used.
        FNAME,        // The name of the output file.
        ORDER,        // The desired order.
        SIZE
    };
}


Dynamic::Dynamic():
    Action(SIZE, "DYNAMIC",
           "The \"DYNAMIC\" command analyses a dynamic transfer map.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of the lattice to be analysed");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "Name of the beam to be used");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of the file to be written", "DYNAMIC");
    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of the analysis; must be at least 2", 6.0);

    registerOwnership(AttributeHandler::COMMAND);
}


Dynamic::Dynamic(const std::string &name, Dynamic *parent):
    Action(name, parent)
{}


Dynamic::~Dynamic()
{}


Dynamic *Dynamic::clone(const std::string &name) {
    return new Dynamic(name, this);
}


void Dynamic::execute() {
    // Find BeamSequence and Beam definitions.
    BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    // Open output file.
    std::string file = Attributes::getString(itsAttr[FNAME]);
    std::ofstream os(file.c_str());
    if(os.bad()) {
        throw OpalException("Dynamic::execute()",
                            "Unable to open file \"" + file + "\".");
    }

    // Reference for beam.
    const PartData &reference = beam->getReference();

    // Compute transfer map.
    int order = std::max(int(std::round(Attributes::getReal(itsAttr[ORDER]))), 1);
    FTps<double, 6>::setGlobalTruncOrder(order);
    LieMapper mapper(*use->fetchLine(), reference, false, false, order);
    mapper.execute();
    DragtFinnMap<3> map;
    mapper.getMap(map);
    os << "Accumulated map:\n" << map;

    // Find fixed point.
    FVector<double, 6> fixedPoint;
    DragtFinnMap<3> closedOrbitMap;
    map.dynamicFixedPoint(fixedPoint, closedOrbitMap);
    std::streamsize old_prec = os.precision(12);
    os << "Fixed point:\n";
    for(int i = 0; i < 6; i++) {
        os << std::setw(24) << fixedPoint[i];
    }
    os << "\nMap around fixed point:\n";
    os << closedOrbitMap;

    // Print normal form.
    // const DragtFinnNormalForm<3> normal(map);
    // const FVector<complex<double>,6> &lambda = normal.eigenValues();

    // 21-06-2000 remove const in order to hope KCC likes it ....!
    DragtFinnNormalForm<3> normal(map);
    const FVector<complex<double>, 6> &lambda = normal.eigenValues();

    os << "\nEigenvalues of the linear map:\n";

    for(int i = 0; i < 6; i++) {
        os << lambda[i] << '\n';
    }

    os << "\nFactorised Lie transformation script(N) for normal form map:\n"
       << normal.normalForm();

    os << "\nFactorised Lie Transformation script(A) for normalising map:\n"
       << normal.normalisingMap();

    // Print fractional tunes.
    int freedom = normal.degreesOfFreedom();
    os << "\nFractional tunes:\n";
    for(int mode = 0; mode < freedom; mode++) {
        os << std::setw(24) << std::arg(lambda[2*mode]) / Physics::two_pi;
    }
    os << '\n';

    // Print logarithmic excitation constants.
    os << "\nExcitation constants:\n";
    for(int mode = 0; mode < freedom; mode++) {
        os << std::setw(24) << log(std::norm(lambda[2*mode])) / 2.0;
    }
    os << '\n';

    // Normalised anharmonicities.
    {
        FMatrix<double, 3, 3> QQ = normal.anharmonicity();
        os << "\nNormalised anharmonicites:\n";

        for(int mode1 = 1; mode1 <= freedom; mode1++) {
            for(int mode2 = mode1; mode2 <= freedom; mode2++) {
                os << "dQ" << mode1 << "/dE" << mode2 << " = "
                   << std::setw(24) << QQ(mode1 - 1, mode2 - 1) << '\n';
            }
        }
    }

    // Invariant polynomials.
    for(int mode = 1; mode <= freedom; mode++) {
        os << "\nInvariant polynomial for mode " << mode << ":\n"
           << normal.invariant(mode - 1);
    }

    os << std::flush;
    os.precision(old_prec);
}