//
// Class AttList
//   The ATTLIST command.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Tables/AttList.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Table.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Elements/AttCell.h"
#include "Elements/OpalElement.h"
#include "Tables/AttWriter.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include <fstream>
#include <iostream>
#include <vector>

using std::vector;


namespace {

    // The attributes of class AttList.
    enum {
        LINE,        // The name of the line to be listed.
        FNAME,       // The name of the file to be written.
        ALL,         // If true, list all columns.
        COLUMN,      // The columns to be written.
        SIZE
    };
}


AttList::AttList():
    Action(SIZE, "ATTLIST",
           "The \"ATTLIST\" statement lists the element strengths "
           "in a beam line or a sequence.") {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of line to be listed");
    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to receive output", "ATTLIST");
    itsAttr[ALL] = Attributes::makeBool
                   ("ALL", "Are all columns desired?");
    itsAttr[COLUMN] = Attributes::makeStringArray
                      ("COLUMN", "The columns to be written");

    registerOwnership(AttributeHandler::STATEMENT);
}


AttList::AttList(const std::string &name, AttList *parent):
    Action(name, parent)
{}


AttList::~AttList()
{}


AttList *AttList::clone(const std::string &name) {
    return new AttList(name, this);
}


void AttList::execute() {
    // Find beam sequence  or table definition.
    const std::string name = Attributes::getString(itsAttr[LINE]);
    const Beamline *line = 0;

    if(Object *obj = OpalData::getInstance()->find(name)) {
        if(BeamSequence *beamLine = dynamic_cast<BeamSequence *>(obj)) {
            line = beamLine->fetchLine();
        } else if(Table *table = dynamic_cast<Table *>(obj)) {
            line = table->getLine();
        } else {
            throw OpalException("Select::execute()",
                                "You cannot do an \"ATTLIST\" on \"" + name +
                                "\", it is neither a line nor a table.");
        }
    } else {
        throw OpalException("Select::execute()",
                            "Object \"" + name + "\" not found.");
    }

    // Select the file to be used.
    const std::string &fileName = Attributes::getString(itsAttr[FNAME]);
    if(fileName == "TERM") {
        writeTable(*line, std::cout);
    } else {
        std::ofstream os(fileName.c_str());

        if(os.good()) {
            writeTable(*line, os);
        } else {
            throw OpalException("AttList::execute()",
                                "Unable to open output stream \"" +
                                fileName + "\".");
        }
    }
}


void AttList::writeTable(const Beamline &line, std::ostream &os) {
    // Construct column access table.
    // This may throw, if a column is unknown.
    vector<std::string> header = Attributes::getStringArray(itsAttr[COLUMN]);
    vector<std::string>::size_type n = header.size();
    vector<AttCell *> buffer(n);
    for(vector<std::string>::size_type i = 0; i < n; ++i) {
        buffer[i] = OpalElement::findRegisteredAttribute(header[i]);
    }

    // Write table descriptors.
    OPALTimer::Timer timer;
    os << "@ TYPE     %s  ATTRIBUTE\n"
       << "@ LINE     %s  " << line.getName() << "\n"
       << "@ DATE     %s  " << timer.date() << "\n"
       << "@ TIME     %s  " << timer.time() << "\n"
       << "@ ORIGIN   %s  OPAL_9.5/4\n"
       << "@ COMMENT  %s  \"\n";
    OpalData::getInstance()->printTitle(os);
    os << "\"\n";

    // Write column header names.
    os << '*';
    for(vector<std::string>::size_type i = 0; i < n; ++i) {
        os << ' ' << header[i];
    }
    os << '\n';

    // Write column header formats.
    os << '$';
    for(vector<std::string>::size_type i = 0; i < n; ++i) {
        os << ' ';
        buffer[i]->printFormat(os);
        buffer[i]->clearValue();
    }
    os << '\n';

    // List the table body.
    AttWriter writer(line, os, buffer);
    writer.execute();
}