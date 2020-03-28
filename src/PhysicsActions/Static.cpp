// ------------------------------------------------------------------------
// $RCSfile: Static.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Static
//   The class for the OPAL STATIC command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:22:04 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "PhysicsActions/Static.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/LieMapper.h"
#include "Attributes/Attributes.h"
#include "FixedAlgebra/DragtFinnMap.h"
#include "FixedAlgebra/DragtFinnNormalForm.h"
#include "FixedAlgebra/FStaticFP.h"
#include "FixedAlgebra/FTps.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Utilities/OpalException.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


// Class Static
// ------------------------------------------------------------------------

// The attributes of class Static.
namespace {
    enum {
        LINE,         // The lattice to be used.
        BEAM,         // The beam to be used.
        FNAME,        // The name of the output file.
        ORDER,        // The desired order.
        SIZE
    };
}


Static::Static():
    Action(SIZE, "STATIC",
           "The \"STATIC\" command analyses a static transfer map.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of the lattice to be analysed");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "Name of the beam to be used");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of the file to be written", "STATIC");
    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "Order of the analysis; must be at least 2", 2.0);

    registerOwnership(AttributeHandler::COMMAND);
}


Static::Static(const std::string &name, Static *parent):
    Action(name, parent)
{}


Static::~Static()
{}


Static *Static::clone(const std::string &name) {
    return new Static(name, this);
}


void Static::execute() {
    // Find BeamSequence and Beam definitions.
    BeamSequence *use = BeamSequence::find(Attributes::getString(itsAttr[LINE]));
    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    // Open output file.
    std::string file = Attributes::getString(itsAttr[FNAME]);
    std::ofstream os(file.c_str());
    if(os.bad()) {
        throw OpalException("Static::execute()",
                            "Unable to open file \"" + file + "\".");
    }

    // Reference attribute.
    const PartData &reference = beam->getReference();

    // Compute transfer map.
    // JMJ comment: minimum sensible order is 2
    int order = std::max(int(std::round(Attributes::getReal(itsAttr[ORDER]))), 2);
    LieMapper mapper(*use->fetchLine(), reference, false, false, order);
    FTps<double, 6>::setGlobalTruncOrder(order);
    std::cout << "mapper.execute(); start " << std::endl;
    mapper.execute();
    std::cout << "mapper.execute(); done " << std::endl;
    DragtFinnMap<3> map;
    mapper.getMap(map);
    os << "Accumulated map:\n" << map;
    std::cout << "Accumulated map done " << std::endl;

    // Print fixed point and map around it.
    std::streamsize old_prec = os.precision(12);
    FVector<double, 6> fixedPoint;
    DragtFinnMap<3> dispersionMap;
    {
        DragtFinnMap<3> temp;
        map.staticFixedPoint(fixedPoint, temp);
        temp.removeDispersion(dispersionMap, map);
    }
    os << "Fixed point:\n";
    for(int i = 0; i < 4; i++) {
        os << std::setw(24) << fixedPoint[i];
    }
    os << "\nDispersion map:\n" << dispersionMap
       << "\nMap with dispersion removed:\n" << map;

    // Print normal form.
    // const DragtFinnNormalForm<3> normal(map);
    // const FVector<std::complex<double>,6> &lambda = normal.eigenValues();

    // 21-06-2000, remove const
    DragtFinnNormalForm<3> normal(map);
    const FVector<std::complex<double>, 6> &lambda = normal.eigenValues();
    os << "\nEigenvalues of the linear map:\n";

    for(int i = 0; i < 4; i++) {
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

        // Invariant polynomials.
        for(int mode = 1; mode <= freedom; mode++) {
            os << "\nInvariant polynomial for mode " << mode << ":\n"
               << normal.invariant(mode - 1);
        }

        os << std::flush;
    }

    os.precision(old_prec);
}