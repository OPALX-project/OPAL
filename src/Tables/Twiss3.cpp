// ------------------------------------------------------------------------
// $RCSfile: Twiss3.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Twiss3
//   The class for OPAL TWISS3 commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Twiss3.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Tables/Twiss.h"
#include "Utilities/OpalException.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

using std::setw;


// Class Twiss3
// ------------------------------------------------------------------------

// The attributes of class Twiss3.
namespace {
    enum {
        TABLE,       // The name of the table to be listed.
        FNAME,       // The name of the file to be written.
        SIZE
    };
}


Twiss3::Twiss3():
    Action(SIZE, "TWISS3",
           "The \"TWISS3\" statement lists a named \"TWISS\" table in "
           "Mais-Ripken representation.") {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of table to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "TWISS3");

    registerOwnership(AttributeHandler::STATEMENT);
}


Twiss3::Twiss3(const std::string &name, Twiss3 *parent):
    Action(name, parent)
{}


Twiss3::~Twiss3()
{}


Twiss3 *Twiss3::clone(const std::string &name) {
    return new Twiss3(name, this);
}


void Twiss3::execute() {
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
                throw OpalException("Twiss3::execute()",
                                    "Unable to open output stream \"" +
                                    fileName + "\".");
            }
        }
    } else {
        throw OpalException("Twiss3::execute()",
                            "Twiss table \"" + tableName + "\" not found.");
    }
}


void Twiss3::format(std::ostream &os, const Twiss *table) {
        formatPrint(os, table);
}


void Twiss3::formatPrint(std::ostream &os, const Twiss *table) const {
    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table specific header.
    table->printTableTitle(os, "Mais-Ripken lattice functions");
    os << std::string(124, '-') << '\n';
    os << "Element" << std::string(24, ' ') << "S"
       << std::string(10, ' ') << "XC" << std::string(9, ' ') << "PXC"
       << std::string(10, ' ') << "YC" << std::string(9, ' ') << "PYC"
       << std::string(10, ' ') << "TC" << std::string(9, ' ') << "PTC\n";
    os << "Mode          MU"
       << "        BETX        GAMX        ALFX"
       << "        BETY        GAMY        ALFY"
       << "        BETT        GAMY        ALFT\n";
    os << std::string(124, '-') << '\n';
    // Jumbled function names affecting PRINT style output
    // fixed at 15:26:46 on 9 Aug 2000 by JMJ

    // Write table body.
    for(Twiss::TLine::const_iterator row = table->begin();
        row != table->end(); ++row) {
        if(row->getSelectionFlag()) {
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
            os << setw(16) << table->getS(*row);

            FVector<double, 6> orbit = table->getOrbit(*row);
            for(int i = 0; i < 6; ++i) {
                os << setw(12) << orbit[i];
            }
            os << '\n';

            for(int mode = 0; mode < 3; ++mode) {
                os << setw(4) << (mode + 1) << setw(12)
                   << table->getMUi(*row, mode);
                for(int plane = 0; plane < 3; ++plane) {
                    os << setw(12) << table->getBETik(*row, plane, mode)
                       << setw(12) << table->getGAMik(*row, plane, mode)
                       << setw(12) << table->getALFik(*row, plane, mode);
                }
                os << '\n';
            }
        }
    }

    os << std::string(124, '-') << std::endl;

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}