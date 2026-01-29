//
// Class OpalWake
//   The class for the OPAL WAKE command.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "Structure/OpalWake.h"
#include "Solvers/GreenWakeFunction.h"
#include "Solvers/CSRWakeFunction.h"
#include "Solvers/CSRIGFWakeFunction.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Utilities/OpalFilter.h"

extern Inform *gmsg;

// The attributes of class OpalWake.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
        TYPE,       // The type of the wake
        NBIN,       // Number of bins for the line density
        CONST_LENGTH,// True if the length of the Bunch is considered as constant
        CONDUCT,    // Conductivity, either AC or DC
        Z0,     //
        RADIUS, // Radius of the tube
        SIGMA,
        TAU,
        FILTERS, // List of filters to apply on line density
        FNAME,
        SIZE
    };
}

OpalWake::OpalWake():
    Definition(SIZE, "WAKE",
               "The \"WAKE\" statement defines data for the wakefuction "
               "on an element."),
    wf_m(0) {
    itsAttr[TYPE] = Attributes::makePredefinedString
        ("TYPE", "Specifies the wake function.", {"1D-CSR", "1D-CSR-IGF", "LONG-SHORT-RANGE", "TRANSV-SHORT-RANGE"});

    itsAttr[NBIN] = Attributes::makeReal
        ("NBIN", "Number of bins for the line density calculation");

    itsAttr[CONST_LENGTH] = Attributes::makeBool
        ("CONST_LENGTH", "True if the length of the Bunch is considered as constant");

    itsAttr[CONDUCT] = Attributes::makePredefinedString
        ("CONDUCT", "Conductivity.", {"DC", "AC"});

    itsAttr[Z0] = Attributes::makeReal
        ("Z0", "Impedance of the beam pipe ");

    itsAttr[RADIUS] = Attributes::makeReal
        ("RADIUS", "The radius of the beam pipe [m]");

    itsAttr[SIGMA] = Attributes::makeReal
        ("SIGMA", "Material constant dependant on the beam pipe material");

    itsAttr[TAU] = Attributes::makeReal
        ("TAU", "Material constant dependant on the beam pipe material");

    itsAttr[FILTERS] = Attributes::makeStringArray
        ("FILTERS", "List of filters to apply on line density");

    itsAttr[FNAME] = Attributes::makeStringArray
        ("FNAME", "Filename of the wakefield file");

    OpalWake *defWake = clone("UNNAMED_WAKE");
    defWake->builtin = true;

    try {
        defWake->update();
        OpalData::getInstance()->define(defWake);
    } catch(...) {
        delete defWake;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


OpalWake::OpalWake(const std::string &name, OpalWake *parent):
    Definition(name, parent),
    wf_m(parent->wf_m)
{}


OpalWake::~OpalWake() {
    delete wf_m;
}


bool OpalWake::canReplaceBy(Object *object) {
    // Can replace only by another WAKE.
    return dynamic_cast<OpalWake *>(object) != 0;
}


OpalWake *OpalWake::clone(const std::string &name) {
    return new OpalWake(name, this);
}


void OpalWake::execute() {
    update();
}


OpalWake *OpalWake::find(const std::string &name) {
    OpalWake *wake = dynamic_cast<OpalWake *>(OpalData::getInstance()->find(name));

    if (wake == 0) {
        throw OpalException("OpalWake::find()", "Wake \"" + name + "\" not found.");
    }
    return wake;
}


int OpalWake::getNumberOfBins() {
    return (int)Attributes::getReal(itsAttr[NBIN]);
}


void OpalWake::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_WAKE");
}


void OpalWake::initWakefunction(const ElementBase &element) {
    *gmsg << "* ************* W A K E ************************************************************\n";
    *gmsg << "OpalWake::initWakefunction ";
    *gmsg << "for element " << element.getName() << "\n";
    *gmsg << "* **********************************************************************************" << endl;


    std::vector<std::string> filters_str = Attributes::getStringArray(itsAttr[FILTERS]);
    std::vector<Filter *> filters;

    for(std::vector<std::string>::const_iterator fit = filters_str.begin(); fit != filters_str.end(); ++ fit) {
        OpalFilter *f = OpalFilter::find(*fit);

        if (f) {
            f->initOpalFilter();
            filters.push_back(f->filter_m);
        }
    }
    std::string type = Attributes::getString(itsAttr[TYPE]);
    if (type == "1D-CSR") {

        if (filters.size() == 0 && Attributes::getReal(itsAttr[NBIN]) <= 7) {
            throw OpalException("OpalWake::initWakeFunction",
                                "At least 8 bins have to be used, ideally far more");
        }

        wf_m = new CSRWakeFunction(getOpalName(),
                                   filters,
                                   (int)(Attributes::getReal(itsAttr[NBIN])));

    } else if (type == "1D-CSR-IGF") {

        if (filters.size() == 0 && Attributes::getReal(itsAttr[NBIN]) <= 7) {
            throw OpalException("OpalWake::initWakeFunction",
                                "At least 8 bins have to be used, ideally far more");
        }

        wf_m = new CSRIGFWakeFunction(getOpalName(),
                                      filters,
                                      (int)(Attributes::getReal(itsAttr[NBIN])));

    } else if (type == "LONG-SHORT-RANGE") {
        int acMode = Attributes::getString(itsAttr[CONDUCT]) == "DC"? 2: 1;

        wf_m = new GreenWakeFunction(getOpalName(),
                                     filters,
                                     (int)(Attributes::getReal(itsAttr[NBIN])),
                                     Attributes::getReal(itsAttr[Z0]),
                                     Attributes::getReal(itsAttr[RADIUS]),
                                     Attributes::getReal(itsAttr[SIGMA]),
                                     acMode,
                                     Attributes::getReal(itsAttr[TAU]),
                                     1,
                                     Attributes::getBool(itsAttr[CONST_LENGTH]),
                                     Attributes::getString(itsAttr[FNAME]));

    } else if (type == "TRANSV-SHORT-RANGE") {
        int acMode = Attributes::getString(itsAttr[CONDUCT]) == "DC" ? 2: 1;

        wf_m = new GreenWakeFunction(getOpalName(),
                                     filters,
                                     (int)(Attributes::getReal(itsAttr[NBIN])),
                                     Attributes::getReal(itsAttr[Z0]),
                                     Attributes::getReal(itsAttr[RADIUS]),
                                     Attributes::getReal(itsAttr[SIGMA]),
                                     acMode,
                                     Attributes::getReal(itsAttr[TAU]),
                                     0,
                                     Attributes::getBool(itsAttr[CONST_LENGTH]),
                                     Attributes::getString(itsAttr[FNAME]));
    }
}

void OpalWake::print(std::ostream &os) const {
    os << "* ************* W A K E ************************************************************ " << std::endl;
    os << "* WAKE         " << getOpalName() << '\n'
       << "* BINS         " << Attributes::getReal(itsAttr[NBIN]) << '\n'
       << "* CONST_LENGTH " << Attributes::getReal(itsAttr[CONST_LENGTH]) << '\n'
       << "* CONDUCT      " << Attributes::getReal(itsAttr[CONDUCT]) << '\n'
       << "* Z0           " << Attributes::getReal(itsAttr[Z0]) << '\n'
       << "* RADIUS       " << Attributes::getReal(itsAttr[RADIUS]) << '\n'
       << "* SIGMA        " << Attributes::getReal(itsAttr[SIGMA]) << '\n'
       << "* TAU          " << Attributes::getReal(itsAttr[TAU]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}