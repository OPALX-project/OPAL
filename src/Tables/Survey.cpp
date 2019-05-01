// ------------------------------------------------------------------------
// $RCSfile: Survey.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Survey
//   The Survey class implements a buffer for survey tables.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/Survey.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "AbstractObjects/TableRowRep.h"
#include "Attributes/Attributes.h"
#include "Algorithms/Surveyor.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "Physics/Physics.h"
#include "Tables/Flatten.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>


// Local structures.
// ------------------------------------------------------------------------

namespace {

    // Describe a column with name and function to get it.
    struct ColDesc {

        // The column name.
        const char *colName;

        // Method to return the selected column from the given row.
        double(Survey::*get)(const Survey::Row &, int, int) const;

        int printWidth;
        int printPrecision;

        int ind_1;
        int ind_2;
    };

    // The column entries table.
    const ColDesc allColumns[] = {
        { "S",     &Survey::getS,     14,  6, 0, 0 },

        { "X",     &Survey::getX,     10,  6, 0, 0 },
        { "Y",     &Survey::getY,     10,  6, 0, 0 },
        { "Z",     &Survey::getZ,     10,  6, 0, 0 },

        { "THETA", &Survey::getTheta, 14, 10, 0, 0 },
        { "PHI",   &Survey::getPhi,   14, 10, 0, 0 },
        { "PSI",   &Survey::getPsi,   14, 10, 0, 0 },

        { "W11",   &Survey::getW,     10,  6, 1, 1 },
        { "W21",   &Survey::getW,     10,  6, 2, 1 },
        { "W31",   &Survey::getW,     10,  6, 3, 1 },
        { "W12",   &Survey::getW,     10,  6, 1, 2 },
        { "W22",   &Survey::getW,     10,  6, 2, 2 },
        { "W32",   &Survey::getW,     10,  6, 3, 2 },
        { "W13",   &Survey::getW,     10,  6, 1, 3 },
        { "W23",   &Survey::getW,     10,  6, 2, 3 },
        { "W33",   &Survey::getW,     10,  6, 3, 3 },

        { 0,       0,                  0,  0, 0, 0 }
    };


    // The print column entries table.
    const ColDesc printColumns[] = {
        { "S",     &Survey::getS,     14,  6, 0, 0 },

        { "X",     &Survey::getX,     10,  6, 0, 0 },
        { "Y",     &Survey::getY,     10,  6, 0, 0 },
        { "Z",     &Survey::getZ,     10,  6, 0, 0 },

        { "THETA", &Survey::getTheta, 14, 10, 0, 0 },
        { "PHI",   &Survey::getPhi,   14, 10, 0, 0 },
        { "PSI",   &Survey::getPsi,   14, 10, 0, 0 },

        { 0,       0,                  0,  0, 0, 0 }
    };


    // Find a column by name.
    const ColDesc *findCol(const Survey &table, const std::string &colName) {
        for(const ColDesc *col = allColumns; col->colName; ++col) {
            if(colName == col->colName) {
                return col;
            }
        }

        throw OpalException("Survey::findCol()",
                            "Survey table \"" + table.getOpalName() +
                            "\" has no column named \"" + colName + "\".");
    }


    // Local class Column.
    //   Returns the value for a given column in the current row of a given
    //   table.
    class Column: public Expressions::Scalar<double> {

    public:

        //: Constructor.
        //  Identify the table by its name [b]tab[/b], and the column by its
        //  name [b]col[/b] and the function [b]col[/b].
        //  The row is specified as the ``current'' row of the table.
        Column(const Survey &tab, const std::string &scol, const ColDesc &col);

        Column(const Column &);
        virtual ~Column();

        //: Make clone.
        virtual Expressions::Scalar<double> *clone() const;

        //: Evaluate.
        virtual double evaluate() const;

        //: Print expression.
        virtual void print(std::ostream &os, int precedence = 99) const;

    private:

        // Not implemented.
        Column();
        const Column &operator=(const Column &);

        // The Table referred.
        const Survey &itsTable;

        // Column name.
        std::string colName;

        // The function returning the column value from the given row.
        double(Survey::*get)(const Survey::Row &, int, int) const;

        // Indices for calling get().
        int ind_1, ind_2;
    };


    // Implementation.
    // ------------------------------------------------------------------------

    Column::Column(const Survey &tab,
                   const std::string &name,
                   const ColDesc &desc):
        itsTable(tab),
        colName(name),
        get(desc.get),
        ind_1(desc.ind_1),
        ind_2(desc.ind_2)
    {}


    Column::Column(const Column &rhs):
        Scalar<double>(rhs),
        itsTable(rhs.itsTable),
        colName(rhs.colName),
        get(rhs.get),
        ind_1(rhs.ind_1),
        ind_2(rhs.ind_2)
    {}


    Column::~Column()
    {}


    Expressions::Scalar<double> *Column::clone() const {
        return new Column(*this);
    }


    double Column::evaluate() const {
        // Apply the function to the current row of the table.
        return (itsTable.*get)(itsTable.getCurrent(), ind_1, ind_2);
    }


    void Column::print(std::ostream &os, int) const {
        os << colName;
    }
};


// Class Survey::Row
// ------------------------------------------------------------------------

Survey::Row::Row(ElementBase *elem, int occur):
    FlaggedElmPtr(elem) {
    setCounter(occur);
}


Survey::Row::Row(const FlaggedElmPtr &rhs):
    FlaggedElmPtr(rhs)
{}


Survey::Row::~Row()
{}


double Survey::Row::getS() const {
    return s;
}


const Euclid3D &Survey::Row::getMap() const {
    return euclid;
}


// Class Survey
// ------------------------------------------------------------------------

// The attributes of class Survey.
namespace {
    enum {
        LINE,        // The lattice to be used.
        RANGE,       // The range in the lattice.
        STATIC,      // The flag for suppressing recalculation.
        INIT,        // If given, take initial values from another table.
        X0,          // The initial position  X.
        Y0,          // The initial position  Y.
        Z0,          // The initial position  Z.
        THETA0,      // The initial azimuth   THETA.
        PHI0,        // The initial elevation PHI.
        PSI0,        // The initial roll      PSI.
        L,           // Sum of design lengths.
        SIZE
    };
}


Survey::Survey():
    Table(SIZE, "SURVEY",
          "The \"SURVEY\" command defines a table of survey data which "
          "can be matched or tabulated."),
    itsTable(0), itsVisitor(0) {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "Name of the lattice to be surveyed");
    itsAttr[RANGE] = Attributes::makeRange
                     ("RANGE", "The range in the lattice");
    itsAttr[STATIC] = Attributes::makeBool
                      ("STATIC", "If true, the table is not recalculated at each iteration");
    itsAttr[INIT] = Attributes::makeTableRow
                    ("INIT", "If given, this table position is used to initialise");
    itsAttr[X0] = Attributes::makeReal
                  ("X0", "Initial X position in m");
    itsAttr[Y0] = Attributes::makeReal
                  ("Y0", "Initial Y position in m");
    itsAttr[Z0] = Attributes::makeReal
                  ("Z0", "Initial Z position in m");
    itsAttr[THETA0] = Attributes::makeReal
                      ("THETA0", "Initial azimuth angle in rad");
    itsAttr[PHI0] = Attributes::makeReal
                    ("PHI0", "Initial pitch angle in rad");
    itsAttr[PSI0] = Attributes::makeReal
                    ("PSI0", "Initial roll angle in rad");

    // READ ONLY.
    itsAttr[L] = Attributes::makeReal
                 ("L", "Sum of design lengths in m");
    itsAttr[L].setReadOnly(true);

    registerOwnership(AttributeHandler::COMMAND);
}


Survey::Survey(const std::string &name, Survey *parent):
    Table(name, parent), itsTable(new TLine(name)), itsVisitor(0)
{}


Survey::~Survey() {
    delete itsTable;
    delete itsVisitor;
}


Survey *Survey::clone(const std::string &name) {
    return new Survey(name, this);
}


void Survey::execute() {
    // Find BeamSequence definition.
    itsLine = Attributes::getString(itsAttr[LINE]);
    BeamSequence *use = BeamSequence::find(itsLine);

    // Make sure all is up-to-date.
    OpalData::getInstance()->update();

    // Create flat list with space for data storage.
    RangeRep range = Attributes::getRange(itsAttr[RANGE]);
    Flatten<Row> flattener(*use->fetchLine(), *itsTable, range);
    flattener.execute();

    if(itsTable->empty() && Options::warn) {
        std::cerr << "\n### Warning ### Survey table \"" << getOpalName()
                  << "\" contains no elements.\n" << std::endl;
    } else {
        itsTable->front().setSelectionFlag(true);
        itsTable->back().setSelectionFlag(true);
    }

    // Fill the table.
    itsVisitor = new Surveyor(*itsTable, false);
    fill();

    // Set static flag.
    if(Attributes::getBool(itsAttr[STATIC])) dynamic = false;
    printTable(std::cout, getDefault());
}


void Survey::fill() {
    if(refill) {
        // Reset the surveyor.
        double x, y, z, phi, theta, psi;
        if(itsAttr[INIT]) {
            TableRowRep rowrep = Attributes::getTableRow(itsAttr[INIT]);
            Table *init = Table::find(rowrep.getTabName());

            if(matches(init)) {
                PlaceRep pinit = rowrep.getPosition();
                Survey  *sinit = dynamic_cast<Survey *>(init);
                sinit->fill();
                const Row &row = sinit->findRow(pinit);
                itsVisitor->setMap(sinit->getMap(row));
            } else {
                throw OpalException("Survey::fill()", "Table \"" + rowrep.getTabName() +
                                    "\" is not suitable for initialising survey \"" +
                                    getOpalName() + "\".");
            }
        } else {
            x     = Attributes::getReal(itsAttr[X0]);
            y     = Attributes::getReal(itsAttr[Y0]);
            z     = Attributes::getReal(itsAttr[Z0]);
            phi   = Attributes::getReal(itsAttr[PHI0]);
            theta = Attributes::getReal(itsAttr[THETA0]);
            psi   = Attributes::getReal(itsAttr[PSI0]);

            // Use the OPAL-8 conventions.
            Rotation3D rot = Rotation3D::YRotation(theta) *
                             Rotation3D::XRotation(- phi) * Rotation3D::ZRotation(psi);
            itsVisitor->setMap(Euclid3D(Vector3D(x, y, z), rot));
        }

        s = 0.0;

        for(TLine::iterator i = itsTable->begin(); i != itsTable->end(); ++i) {
            // Advance through element.
            i->accept(*itsVisitor);

            // Update accumulated length.
            ElementBase &elem = *i->getElement();
            if(! dynamic_cast<Beamline *>(&elem)) {
                s += elem.getElementLength();
            }

            // Copy row to table line.
            itsVisitor->getMap(i->euclid);
            i->s = s;
        }

        Euclid3D map;
        itsVisitor->getMap(map);
        map.getAll(x, y, z, phi, theta, psi);
        Attributes::setReal(itsAttr[L], s);
        refill = false;
    }
}


const Survey::Row &Survey::findRow(const PlaceRep &place) {
    PlaceRep row(place);
    row.initialize();
    for(TLine::const_iterator i = itsTable->begin();
        i != itsTable->end(); ++i) {
        row.enter(*i);
        if(row.isActive()) return *i;
        row.leave(*i);
    }

    std::ostringstream os;
    os << row << std::ends;
    throw OpalException("Survey::setRow()", "Row \"" + os.str() +
                        "\" not found in survey table \"" + getOpalName() + "\".");
}


double Survey::getCell(const PlaceRep &place, const std::string &name) {
    const Row &row = findRow(place);
    const ColDesc *col = findCol(*this, name);
    return (this->*(col->get))(row, col->ind_1, col->ind_2);
}


Table::CellArray Survey::getDefault() const {
    CellArray columns;
    for(const ColDesc *col = printColumns; col->colName; ++col) {
        Expressions::PtrToScalar<double> expr =
            new Column(*this, col->colName, *col);
        columns.push_back(Cell(expr, col->printWidth, col->printPrecision));
    }
    return columns;
}


std::vector<double> Survey::getColumn(const RangeRep &rng, const std::string &name) {
    // Find proper column function.
    const ColDesc *col = findCol(*this, name);
    RangeRep range(rng);
    range.initialize();
    std::vector<double> column;

    for(current = itsTable->begin(); current != itsTable->end(); ++current) {
        range.enter(*current);
        if(range.isActive()) {
            column.push_back((this->*(col->get))(*current, col->ind_1, col->ind_2));
        }
        range.leave(*current);
    }

    return column;
}


const Survey::Row &Survey::getCurrent() const {
    return *current;
}


double Survey::getLength() {
    return itsTable->getElementLength();
}


const Beamline *Survey::getLine() const {
    return itsTable;
}


std::vector<double>
Survey::getRow(const PlaceRep &pos, const std::vector<std::string> &columns) {
    const Row &row = findRow(pos);
    std::vector<double> result;

    if(columns.empty()) {
        // Standard column selection.
        for(const ColDesc *col = printColumns; col->colName != 0; ++col) {
            result.push_back((this->*(col->get))(row, col->ind_1, col->ind_2));
        }
    } else {
        // User column selection.
        for(std::vector<std::string>::const_iterator name = columns.begin();
            name != columns.end(); ++name) {
            const ColDesc *col = findCol(*this, *name);
            result.push_back((this->*(col->get))(row, col->ind_1, col->ind_2));
        }
    }

    return result;
}


bool Survey::isDependent(const std::string &name) const {
    // Test if name refers to LINE attribute.
    if(itsLine == name) return true;

    // Test if name occurs in table.
    for(TLine::const_iterator i = itsTable->begin(); i != itsTable->end(); ++i) {
        if(i->getElement()->getName() == name) return true;
    }

    // Otherwise replacement is not required.
    return false;
}


Expressions::PtrToScalar<double>
Survey::makeColumnExpression(const std::string &colName) const {
    return new Column(*this, colName, *findCol(*this, colName));
}


bool Survey::matches(Table *rhs) const {
    return dynamic_cast<Survey *>(rhs) != 0;
}


void Survey::printTable(std::ostream &os, const CellArray &cells) const {
    // Write table specific header.
    OpalData::getInstance()->printTitle(os);
    os << '\n';
    os << "Geometric layout for line: " << itsAttr[LINE]
       << ", range: " << itsAttr[RANGE] << "." << '\n';

    // Find line length.
    int lineLength = 16;
    for(CellArray::const_iterator cell = cells.begin();
        cell < cells.end(); ++cell) {
        lineLength += cell->printWidth;
    }
    os << std::string(lineLength, '-') << '\n';

    // Print table header.
    os << "Element         ";
    for(CellArray::const_iterator cell = cells.begin();
        cell < cells.end(); ++cell) {
        std::ostringstream ss;
        cell->itsExpr->print(ss, 0);
        ss << std::ends;
        std::string image = ss.str();

        if(int(image.length()) < cell->printWidth) {
            // Right adjust the column header.
            os << std::string(cell->printWidth - image.length(), ' ') << image;
        } else {
            // Truncate the column header.
            os << ' ' << std::string(image, 0, cell->printWidth - 3) << "..";
        }
    }
    os << '\n';
    os << std::string(lineLength, '-') << '\n';

    // Save the formatting flags.
    std::streamsize old_prec = os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);

    // Write table body.
    for(current = itsTable->begin(); current != itsTable->end(); ++current) {
        if(current->getSelectionFlag()) {
            std::string name = current->getElement()->getName();
            if(int occur = current->getCounter()) {
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

            for(CellArray::const_iterator cell = cells.begin();
                cell != cells.end(); ++cell) {
                os << std::setw(cell->printWidth)
                   << std::setprecision(cell->printPrecision)
                   << cell->itsExpr->evaluate();
            }
            os << '\n';
        }
    }

    // Write table specific summary.
    os << std::string(lineLength, '-') << '\n';

    // Restore the formatting flags.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
}


double Survey::getS(const Survey::Row &row, int, int) const {
    return current->getS();
}


const Euclid3D &Survey::getMap(const Survey::Row &row) const {
    return row.getMap();
}


double Survey::getX(const Survey::Row &row, int, int) const {
    return row.getMap().getX();
}


double Survey::getY(const Survey::Row &row, int, int) const {
    return row.getMap().getY();
}


double Survey::getZ(const Survey::Row &row, int, int) const {
    return row.getMap().getZ();
}


double Survey::getPhi(const Survey::Row &row, int, int) const {
    const Rotation3D &rot = row.getMap().getRotation();
    return atan2(rot(1, 2), sqrt(rot(1, 0) * rot(1, 0) + rot(1, 1) * rot(1, 1)));
}


double Survey::getTheta(const Survey::Row &row, int, int) const {
    const Rotation3D &rot = row.getMap().getRotation();
    double arg = std::abs(rot(0, 2)) + std::abs(rot(2, 2));
    return (arg > 1.0e-10) ? atan2(rot(0, 2), rot(2, 2)) : 0.0;
}


double Survey::getPsi(const Survey::Row &row, int, int) const {
    const Rotation3D &rot = row.getMap().getRotation();
    double arg = std::abs(rot(1, 0)) + std::abs(rot(1, 1));
    return (arg > 1.0e-10) ? atan2(rot(1, 0), rot(1, 1)) : 0.0;
}


double Survey::getW(const Survey::Row &row, int ind_1, int ind_2) const {
    return row.getMap().M(ind_1, ind_2);
}