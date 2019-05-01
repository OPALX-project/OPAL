// ------------------------------------------------------------------------
// $RCSfile: Envelope.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Envelope
//   The class for OPAL ENVELOPE commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Envelope.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Tables/Twiss.h"
#include "Utilities/OpalException.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

// Class Envelope
// ------------------------------------------------------------------------

// The attributes of class Envelope.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


Envelope::Envelope():
    Action(SIZE, "ENVELOPE",
           "The \"ENVELOPE\" statement lists the beam envelope for a named "
           "\"TWISS\" table.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "ENVELOPE");

    registerOwnership(AttributeHandler::STATEMENT);
}


Envelope::Envelope(const std::string &name, Envelope *parent):
    Action(name, parent)
{}


Envelope::~Envelope()
{}


Envelope *Envelope::clone(const std::string &name) {
    return new Envelope(name, this);
}


void Envelope::execute() {
    std::string tableName = Attributes::getString(itsAttr[TABLE]);
    Twiss *table = dynamic_cast<Twiss *>(OpalData::getInstance()->find(tableName));

    if(table) {
        std::string fileName = Attributes::getString(itsAttr[FNAME]);
        if(fileName == "TERM") {
            format(std::cout, table);
        } else {
            std::ofstream os(fileName.c_str());

            if(os.good()) {
                format(os, table);
            } else {
                throw OpalException("Envelope::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("Envelope::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void Envelope::format(std::ostream &os, const Twiss *table) {
        formatPrint(os, table);
}


void Envelope::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Print table header.
    table->printTableTitle(os, "Beam envelope");
    os << std::string(118, '-') << '\n';
    os << "Element" << std::string(24, ' ') << "S       Orbit |"
       << std::string(25, ' ') << "S i g m a   M a t r i x\n";
    os << std::string(118, '-') << '\n';

    // Print table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
            os << '\n';
            std::string name = row->getElement()->getName();
            if(int occur = row->getCounter()) {
                std::ostringstream tos;
                tos << name << '[' << occur << ']' << std::ends;
                name = tos.str();
            }

            if(name.length() > 16) {
                // Truncate the element name.
                os << std::string(name, 0, 13) << ".. ";
            } else {
                // Left adjust the element name.
                os << name << std::string(16 - name.length(), ' ');
            }
            os << std::setw(16) << table->getS(*row);

            FVector<double, 6> orbit = table->getOrbit(*row);
            FMatrix<double, 6, 6> sigma = table->getSigma(*row);
            for(int i = 0; i < 6; ++i) {
                if(i != 0) os << std::string(32, ' ');
                os << std::setw(12) << orbit[i] << " |";
                for(int j = 0; j <= i; ++j) {
                    os << std::setw(12) << sigma[i][j];
                }
                os << '\n';
            }
        }
    }

    os << std::string(118, '-') << std::endl;

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}