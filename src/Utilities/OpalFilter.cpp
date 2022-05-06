//
// Class OpalFilter
//   The class for the Filter Object.
//   A FILTER definition is used to define a filter which can be applied
//   to a 1D histogram in order to get rid of noise.
//
// Copyright (c) 2008 - 2022, Christof Metzger-Kraus
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
#include "Utilities/OpalFilter.h"

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Filters/Filters.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include "Utility/IpplInfo.h"

#include <cmath>
#include <unordered_map>

#define NPOINTS_DEFAULT 129
#define NLEFT_DEFAULT 64
#define NRIGHT_DEFAULT 64
#define POLYORDER_DEFAULT 1

extern Inform *gmsg;

// The attributes of class OpalFilter.
namespace {
    enum {
        TYPE,        // The type of filter
        NFREQ,       // Number of frequencies in fixedFFTLowPass filter
        THRESHOLD,   // Relative threshold for amplitude of frequency in relativeFFTLowPass
        NPOINTS,     // Number of points in Savitzky-Golay filter
        NLEFT,       // Number of points to the left in S-G filter
        NRIGHT,      // Number of points to the right in S-G filter
        POLYORDER,   // Polynomial order in S-G Filter
        SIZE
    };
}

OpalFilter::OpalFilter():
    Definition(SIZE, "FILTER",
               "The \"FILTER\" statement defines a 1 dimensional filter to be "
               "applied on histogram."),
    filter_m(0) {
    itsAttr[TYPE] = Attributes::makePredefinedString
        ("TYPE", "Specifies the type of filter.",
         {"SAVITZKY-GOLAY", "FIXEDFFTLOWPASS", "RELATIVEFFTLOWPASS", "STENCIL"});

    itsAttr[NFREQ] = Attributes::makeReal
        ("NFREQ", "Number of frequencies to use in fixedFFTLowPass filter", 9.);

    itsAttr[THRESHOLD] = Attributes::makeReal
        ("THRESHOLD", "Relative threshold for amplitude of frequencies in relativeFFTLowPass filter", 1.e-6);

    itsAttr[NPOINTS] = Attributes::makeReal
        ("NPOINTS", "Number of points in Savitzky-Golay filter", NPOINTS_DEFAULT);

    itsAttr[NLEFT] = Attributes::makeReal
        ("NLEFT", "Number of points to the left in Savitzky-Golay filter", NLEFT_DEFAULT);

    itsAttr[NRIGHT] = Attributes::makeReal
        ("NRIGHT", "Number of points to the right in Savitzky-Golay filter", NRIGHT_DEFAULT);

    itsAttr[POLYORDER] = Attributes::makeReal
        ("POLYORDER", "Polynomial order for local fit-function in Savitzky-Golay filter", POLYORDER_DEFAULT);

    registerOwnership(AttributeHandler::STATEMENT);

    OpalFilter* defFilter = clone("UNNAMED_FILTER");
    defFilter->builtin = true;

    try {
        defFilter->update();
        OpalData::getInstance()->define(defFilter);
    } catch(...) {
        delete defFilter;
    }
}


OpalFilter::OpalFilter(const std::string& name, OpalFilter* parent):
    Definition(name, parent),
    filter_m(0)
{}


OpalFilter::~OpalFilter() {
    if (filter_m)
        delete filter_m;
}


bool OpalFilter::canReplaceBy(Object* object) {
    // Can replace only by another WAKE.
    return dynamic_cast<OpalFilter*>(object) != 0;
}


OpalFilter *OpalFilter::clone(const std::string& name) {
    return new OpalFilter(name, this);
}


void OpalFilter::execute() {
    update();
}


OpalFilter* OpalFilter::find(const std::string& name) {
    OpalFilter* filter = dynamic_cast<OpalFilter*>(OpalData::getInstance()->find(name));

    if (filter == 0) {
        throw OpalException("OpalFilter::find()",
                            "OpalFilter \"" + name + "\" not found.");
    }
    return filter;
}


void OpalFilter::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_FILTER");
}


void OpalFilter::initOpalFilter() {
    if (filter_m == 0) {
        *gmsg << "* ************* F I L T E R ************************************************************" << endl;
        *gmsg << "OpalFilter::initOpalFilterfunction " << endl;
        *gmsg << "* **********************************************************************************" << endl;

        static const std::unordered_map<std::string, FilterType> stringFilterType_s = {
            {"SAVITZKY-GOLAY",     FilterType::SAVITZKYGOLAY},
            {"FIXEDFFTLOWPASS",    FilterType::FIXEDFFTLOWPASS},
            {"RELATIVEFFTLOWPASS", FilterType::RELATIVEFFTLOWPASS},
            {"STENCIL",            FilterType::STENCIL}
        };

        FilterType type = stringFilterType_s.at(Attributes::getString(itsAttr[TYPE]));
        switch (type) {
            case FilterType::SAVITZKYGOLAY: {
                int num_points = (int)(Attributes::getReal(itsAttr[NPOINTS]));
                int num_points_left = (int)(Attributes::getReal(itsAttr[NLEFT]));
                int num_points_right = (int)(Attributes::getReal(itsAttr[NRIGHT]));
                int polynomial_order = std::abs((int)(Attributes::getReal(itsAttr[POLYORDER])));

                Inform svg("Savitzky-Golay: ");
                if (num_points_left < 0) {
                    svg << "Number of points to the left negative; using default (" << NLEFT_DEFAULT << ");" << endl;
                    num_points_left = NLEFT_DEFAULT;
                }
                if (num_points_right < 0) {
                    svg << "Number of points to the right negative; using default (" << NRIGHT_DEFAULT << ");" << endl;
                    num_points_right = NRIGHT_DEFAULT;
                }
                if (num_points < num_points_left + num_points_right) {
                    svg << "Total number of points small than sum of the ones to the left and to the right plus 1; using default (NLEFT + NRIGHT + 1);" << endl;
                    num_points = num_points_left + num_points_right + 1;
                }
                if (polynomial_order > num_points_left + num_points_right) {
                    svg << "Polynomial order bigger than sum of points to the left and to the right; using default (NLEFT + NRIGHT);" << endl;
                    polynomial_order = num_points_left + num_points_right;
                }

                filter_m = new SavitzkyGolayFilter(num_points, num_points_left, num_points_right, polynomial_order);
                break;
            }
            case FilterType::FIXEDFFTLOWPASS: {
                filter_m = new FixedFFTLowPassFilter(std::abs((int)(Attributes::getReal(itsAttr[NFREQ]))));
                break;
            }
            case FilterType::RELATIVEFFTLOWPASS: {
                filter_m = new RelativeFFTLowPassFilter(std::abs(Attributes::getReal(itsAttr[THRESHOLD])));
                break;
            }
            case FilterType::STENCIL: {
                filter_m = new StencilFilter();
                break;
            }
            default: {
                filter_m = 0;
                INFOMSG("No filter attached" << endl);
                break;
            }
        }
    }
}

void OpalFilter::print(std::ostream& os) const {
    os << "* ************* F I L T E R ********************************************************\n"
       << "* FILTER         " << getOpalName() << '\n'
       << "* TYPE           " << Attributes::getString(itsAttr[TYPE]) << '\n'
       << "* NFREQ          " << Attributes::getReal(itsAttr[NFREQ]) << '\n'
       << "* THRESHOLD      " << Attributes::getReal(itsAttr[THRESHOLD]) << '\n'
       << "* NPOINTS        " << Attributes::getReal(itsAttr[NPOINTS]) << '\n'
       << "* NLEFT          " << Attributes::getReal(itsAttr[NLEFT]) << '\n'
       << "* NRIGHT         " << Attributes::getReal(itsAttr[NRIGHT]) << '\n'
       << "* POLYORDER      " << Attributes::getReal(itsAttr[POLYORDER]) << '\n'
       << "* ********************************************************************************** " << std::endl;
}
