// ------------------------------------------------------------------------
// $RCSfile: OpalFilter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Filter
//   The class for the Filter Object.
//
// $Date: 2008/10/14 22:09:00 $
// $Author: C. Kraus $
//
// ------------------------------------------------------------------------

#include "Utilities/OpalFilter.h"

#include "Filters/Filters.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"

#include "Utility/IpplInfo.h"
#include <cmath>

#define NPOINTS_DEFAULT 129
#define NLEFT_DEFAULT 64
#define NRIGHT_DEFAULT 64
#define POLYORDER_DEFAULT 1

extern Inform *gmsg;

using namespace Physics;


// Class OpalFilter
// ------------------------------------------------------------------------

// The attributes of class OpalFilter.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
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
    itsAttr[TYPE] = Attributes::makeUpperCaseString
                    ("TYPE", "Specifies the type of filter: SavitzkyGolay, fixedFFTLowPass, relativeFFTLowPass, Stencil");

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

    OpalFilter *defFilter = clone("UNNAMED_FILTER");
    defFilter->builtin = true;

    try {
        defFilter->update();
        OpalData::getInstance()->define(defFilter);
    } catch(...) {
        delete defFilter;
    }
}


OpalFilter::OpalFilter(const std::string &name, OpalFilter *parent):
    Definition(name, parent),
    filter_m(0)
{}


OpalFilter::~OpalFilter() {
    if (filter_m)
        delete filter_m;
}


bool OpalFilter::canReplaceBy(Object *object) {
    // Can replace only by another WAKE.
    return dynamic_cast<OpalFilter *>(object) != 0;
}


OpalFilter *OpalFilter::clone(const std::string &name) {
    return new OpalFilter(name, this);
}


void OpalFilter::execute() {
    update();
}


OpalFilter *OpalFilter::find(const std::string &name) {
    OpalFilter *filter = dynamic_cast<OpalFilter *>(OpalData::getInstance()->find(name));

    if (filter == 0) {
        throw OpalException("OpalFilter::find()", "OpalFilter \"" + name + "\" not found.");
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

        std::string type = Attributes::getString(itsAttr[TYPE]);
        if (type == "SAVITZKY-GOLAY") {
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
        } else if (type == "FIXEDFFTLOWPASS") {
            filter_m = new FixedFFTLowPassFilter(std::abs((int)(Attributes::getReal(itsAttr[NFREQ]))));
        } else if (type == "RELATIVEFFTLOWPASS") {
            filter_m = new RelativeFFTLowPassFilter(std::abs(Attributes::getReal(itsAttr[THRESHOLD])));
        } else if (type == "STENCIL") {
            filter_m = new StencilFilter();
        } else {
            filter_m = 0;
            INFOMSG("no filter attached" << endl);
        }
    }
}

void OpalFilter::print(std::ostream &os) const {
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