// ------------------------------------------------------------------------
// $RCSfile: Twiss.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4.2.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Twiss
//   The abstract class Twiss implements the interface for a table buffer
//   holding lattice function.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/Twiss.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/RangeRep.h"
#include "Algorithms/ThickMapper.h"
#include "Algorithms/LinearMapper.h"
#include "Algorithms/ThinMapper.h"
#include "Attributes/Attributes.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Physics/Physics.h"
#include "Structure/Beam.h"
#include "Tables/Flatten.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Round.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

// Local structures.
// ------------------------------------------------------------------------

namespace {

    struct ColDesc {

        // Column name.
        const char *colName;

        // Pointer to the row method which returns the column value.
        double(Twiss::*get)(const Twiss::Row &, int, int) const;

        // Default field width and precision.
        int printWidth, printPrecision;

        // Indices to be given to get().
        int ind_1, ind_2;
    };


    // The complete column entries table.
    const ColDesc allColumns[] = {

        // Arc length.
        { "S",     &Twiss::getS,     14, 6, 0, 0 },

        // "Naive" Twiss.
        { "ALFX",  &Twiss::getALFi,  10, 6, 0, 0 },
        { "ALFY",  &Twiss::getALFi,  10, 6, 1, 0 },
        { "ALFT",  &Twiss::getALFi,  10, 6, 2, 0 },

        { "BETX",  &Twiss::getBETi,  10, 3, 0, 0 },
        { "BETY",  &Twiss::getBETi,  10, 3, 1, 0 },
        { "BETT",  &Twiss::getBETi,  10, 3, 2, 0 },

        { "MUX",   &Twiss::getMUi,   10, 6, 0, 0 },
        { "MUY",   &Twiss::getMUi,   10, 6, 1, 0 },
        { "MUT",   &Twiss::getMUi,   10, 6, 2, 0 },

        // Mais-Ripken functions.
        { "ALFX1", &Twiss::getALFik,  8, 3, 0, 0 },
        { "ALFX2", &Twiss::getALFik,  8, 3, 0, 1 },
        { "ALFX3", &Twiss::getALFik,  8, 3, 0, 2 },
        { "ALFY1", &Twiss::getALFik,  8, 3, 1, 0 },
        { "ALFY2", &Twiss::getALFik,  8, 3, 1, 1 },
        { "ALFY3", &Twiss::getALFik,  8, 3, 1, 2 },
        { "ALFT1", &Twiss::getALFik,  8, 3, 2, 0 },
        { "ALFT2", &Twiss::getALFik,  8, 3, 2, 1 },
        { "ALFT3", &Twiss::getALFik,  8, 3, 2, 2 },

        { "BETX1", &Twiss::getBETik, 10, 3, 0, 0 },
        { "BETX2", &Twiss::getBETik, 10, 3, 0, 1 },
        { "BETX3", &Twiss::getBETik, 10, 3, 0, 2 },
        { "BETY1", &Twiss::getBETik, 10, 3, 1, 0 },
        { "BETY2", &Twiss::getBETik, 10, 3, 1, 1 },
        { "BETY3", &Twiss::getBETik, 10, 3, 1, 2 },
        { "BETT1", &Twiss::getBETik, 10, 3, 2, 0 },
        { "BETT2", &Twiss::getBETik, 10, 3, 2, 1 },
        { "BETT3", &Twiss::getBETik, 10, 3, 2, 2 },

        { "GAMX1", &Twiss::getGAMik, 10, 3, 0, 0 },
        { "GAMX2", &Twiss::getGAMik, 10, 3, 0, 1 },
        { "GAMX3", &Twiss::getGAMik, 10, 3, 0, 2 },
        { "GAMY1", &Twiss::getGAMik, 10, 3, 1, 0 },
        { "GAMY2", &Twiss::getGAMik, 10, 3, 1, 1 },
        { "GAMY3", &Twiss::getGAMik, 10, 3, 1, 2 },
        { "GAMT1", &Twiss::getGAMik, 10, 3, 2, 0 },
        { "GAMT2", &Twiss::getGAMik, 10, 3, 2, 1 },
        { "GAMT3", &Twiss::getGAMik, 10, 3, 2, 2 },

        // Closed orbit.
        { "XC",    &Twiss::getCO,    10, 6, 0, 0 },
        { "PXC",   &Twiss::getCO,    10, 6, 1, 0 },
        { "YC",    &Twiss::getCO,    10, 6, 2, 0 },
        { "PYC",   &Twiss::getCO,    10, 6, 3, 0 },
        { "TC",    &Twiss::getCO,    10, 6, 4, 0 },
        { "PTC",   &Twiss::getCO,    10, 6, 5, 0 },

        // Dispersion.
        { "DX",    &Twiss::getDisp,  10, 6, 0, 0 },
        { "DPX",   &Twiss::getDisp,  10, 6, 1, 0 },
        { "DY",    &Twiss::getDisp,  10, 6, 2, 0 },
        { "DPY",   &Twiss::getDisp,  10, 6, 3, 0 },
        { "DT",    &Twiss::getDisp,  10, 6, 4, 0 },
        { "DPT",   &Twiss::getDisp,  10, 6, 5, 0 },

        // Eigenvectors.
        { "E11",   &Twiss::getEigen, 10, 6, 0, 0 },
        { "E21",   &Twiss::getEigen, 10, 6, 1, 0 },
        { "E31",   &Twiss::getEigen, 10, 6, 2, 0 },
        { "E41",   &Twiss::getEigen, 10, 6, 3, 0 },
        { "E51",   &Twiss::getEigen, 10, 6, 4, 0 },
        { "E61",   &Twiss::getEigen, 10, 6, 5, 0 },
        { "E12",   &Twiss::getEigen, 10, 6, 0, 1 },
        { "E22",   &Twiss::getEigen, 10, 6, 1, 1 },
        { "E32",   &Twiss::getEigen, 10, 6, 2, 1 },
        { "E42",   &Twiss::getEigen, 10, 6, 3, 1 },
        { "E52",   &Twiss::getEigen, 10, 6, 4, 1 },
        { "E62",   &Twiss::getEigen, 10, 6, 5, 1 },
        { "E13",   &Twiss::getEigen, 10, 6, 0, 2 },
        { "E23",   &Twiss::getEigen, 10, 6, 1, 2 },
        { "E33",   &Twiss::getEigen, 10, 6, 2, 2 },
        { "E43",   &Twiss::getEigen, 10, 6, 3, 2 },
        { "E53",   &Twiss::getEigen, 10, 6, 4, 2 },
        { "E63",   &Twiss::getEigen, 10, 6, 5, 2 },
        { "E14",   &Twiss::getEigen, 10, 6, 0, 3 },
        { "E24",   &Twiss::getEigen, 10, 6, 1, 3 },
        { "E34",   &Twiss::getEigen, 10, 6, 2, 3 },
        { "E44",   &Twiss::getEigen, 10, 6, 3, 3 },
        { "E54",   &Twiss::getEigen, 10, 6, 4, 3 },
        { "E64",   &Twiss::getEigen, 10, 6, 5, 3 },
        { "E15",   &Twiss::getEigen, 10, 6, 0, 4 },
        { "E25",   &Twiss::getEigen, 10, 6, 1, 4 },
        { "E35",   &Twiss::getEigen, 10, 6, 2, 4 },
        { "E45",   &Twiss::getEigen, 10, 6, 3, 4 },
        { "E55",   &Twiss::getEigen, 10, 6, 4, 4 },
        { "E65",   &Twiss::getEigen, 10, 6, 5, 4 },
        { "E16",   &Twiss::getEigen, 10, 6, 0, 5 },
        { "E26",   &Twiss::getEigen, 10, 6, 1, 5 },
        { "E36",   &Twiss::getEigen, 10, 6, 2, 5 },
        { "E46",   &Twiss::getEigen, 10, 6, 3, 5 },
        { "E56",   &Twiss::getEigen, 10, 6, 4, 5 },
        { "E66",   &Twiss::getEigen, 10, 6, 5, 5 },

        // Sigma matrix.
        { "S11",   &Twiss::getSigma, 10, 6, 0, 0 },
        { "S21",   &Twiss::getSigma, 10, 6, 1, 0 },
        { "S31",   &Twiss::getSigma, 10, 6, 2, 0 },
        { "S41",   &Twiss::getSigma, 10, 6, 3, 0 },
        { "S51",   &Twiss::getSigma, 10, 6, 4, 0 },
        { "S61",   &Twiss::getSigma, 10, 6, 5, 0 },
        { "S12",   &Twiss::getSigma, 10, 6, 0, 1 },
        { "S22",   &Twiss::getSigma, 10, 6, 1, 1 },
        { "S32",   &Twiss::getSigma, 10, 6, 2, 1 },
        { "S42",   &Twiss::getSigma, 10, 6, 3, 1 },
        { "S52",   &Twiss::getSigma, 10, 6, 4, 1 },
        { "S62",   &Twiss::getSigma, 10, 6, 5, 1 },
        { "S13",   &Twiss::getSigma, 10, 6, 0, 2 },
        { "S23",   &Twiss::getSigma, 10, 6, 1, 2 },
        { "S33",   &Twiss::getSigma, 10, 6, 2, 2 },
        { "S43",   &Twiss::getSigma, 10, 6, 3, 2 },
        { "S53",   &Twiss::getSigma, 10, 6, 4, 2 },
        { "S63",   &Twiss::getSigma, 10, 6, 5, 2 },
        { "S14",   &Twiss::getSigma, 10, 6, 0, 3 },
        { "S24",   &Twiss::getSigma, 10, 6, 1, 3 },
        { "S34",   &Twiss::getSigma, 10, 6, 2, 3 },
        { "S44",   &Twiss::getSigma, 10, 6, 3, 3 },
        { "S54",   &Twiss::getSigma, 10, 6, 4, 3 },
        { "S64",   &Twiss::getSigma, 10, 6, 5, 3 },
        { "S15",   &Twiss::getSigma, 10, 6, 0, 4 },
        { "S25",   &Twiss::getSigma, 10, 6, 1, 4 },
        { "S35",   &Twiss::getSigma, 10, 6, 2, 4 },
        { "S45",   &Twiss::getSigma, 10, 6, 3, 4 },
        { "S55",   &Twiss::getSigma, 10, 6, 4, 4 },
        { "S65",   &Twiss::getSigma, 10, 6, 5, 4 },
        { "S16",   &Twiss::getSigma, 10, 6, 0, 5 },
        { "S26",   &Twiss::getSigma, 10, 6, 1, 5 },
        { "S36",   &Twiss::getSigma, 10, 6, 2, 5 },
        { "S46",   &Twiss::getSigma, 10, 6, 3, 5 },
        { "S56",   &Twiss::getSigma, 10, 6, 4, 5 },
        { "S66",   &Twiss::getSigma, 10, 6, 5, 5 },

        // Transfer matrix.
        { "R11",   &Twiss::getMatrix, 10, 6, 0, 0 },
        { "R21",   &Twiss::getMatrix, 10, 6, 1, 0 },
        { "R31",   &Twiss::getMatrix, 10, 6, 2, 0 },
        { "R41",   &Twiss::getMatrix, 10, 6, 3, 0 },
        { "R51",   &Twiss::getMatrix, 10, 6, 4, 0 },
        { "R61",   &Twiss::getMatrix, 10, 6, 5, 0 },
        { "R12",   &Twiss::getMatrix, 10, 6, 0, 1 },
        { "R22",   &Twiss::getMatrix, 10, 6, 1, 1 },
        { "R32",   &Twiss::getMatrix, 10, 6, 2, 1 },
        { "R42",   &Twiss::getMatrix, 10, 6, 3, 1 },
        { "R52",   &Twiss::getMatrix, 10, 6, 4, 1 },
        { "R62",   &Twiss::getMatrix, 10, 6, 5, 1 },
        { "R13",   &Twiss::getMatrix, 10, 6, 0, 2 },
        { "R23",   &Twiss::getMatrix, 10, 6, 1, 2 },
        { "R33",   &Twiss::getMatrix, 10, 6, 2, 2 },
        { "R43",   &Twiss::getMatrix, 10, 6, 3, 2 },
        { "R53",   &Twiss::getMatrix, 10, 6, 4, 2 },
        { "R63",   &Twiss::getMatrix, 10, 6, 5, 2 },
        { "R14",   &Twiss::getMatrix, 10, 6, 0, 3 },
        { "R24",   &Twiss::getMatrix, 10, 6, 1, 3 },
        { "R34",   &Twiss::getMatrix, 10, 6, 2, 3 },
        { "R44",   &Twiss::getMatrix, 10, 6, 3, 3 },
        { "R54",   &Twiss::getMatrix, 10, 6, 4, 3 },
        { "R64",   &Twiss::getMatrix, 10, 6, 5, 3 },
        { "R15",   &Twiss::getMatrix, 10, 6, 0, 4 },
        { "R25",   &Twiss::getMatrix, 10, 6, 1, 4 },
        { "R35",   &Twiss::getMatrix, 10, 6, 2, 4 },
        { "R45",   &Twiss::getMatrix, 10, 6, 3, 4 },
        { "R55",   &Twiss::getMatrix, 10, 6, 4, 4 },
        { "R65",   &Twiss::getMatrix, 10, 6, 5, 4 },
        { "R16",   &Twiss::getMatrix, 10, 6, 0, 5 },
        { "R26",   &Twiss::getMatrix, 10, 6, 1, 5 },
        { "R36",   &Twiss::getMatrix, 10, 6, 2, 5 },
        { "R46",   &Twiss::getMatrix, 10, 6, 3, 5 },
        { "R56",   &Twiss::getMatrix, 10, 6, 4, 5 },
        { "R66",   &Twiss::getMatrix, 10, 6, 5, 5 },

        { 0,       0,                 0, 0, 0, 0 }
    };


    // The default column entries table
    const ColDesc defaultColumns[] = {
        { "S",     &Twiss::getS,    14, 6, 0, 0 },

        { "MUX",   &Twiss::getMUi,  10, 6, 0, 0 },
        { "BETX",  &Twiss::getBETi, 10, 3, 0, 0 },
        { "ALFX",  &Twiss::getALFi, 10, 6, 0, 0 },

        { "MUY",   &Twiss::getMUi,  10, 6, 1, 0 },
        { "BETY",  &Twiss::getBETi, 10, 3, 1, 0 },
        { "ALFY",  &Twiss::getALFi, 10, 6, 1, 0 },

        { "XC",    &Twiss::getCO,   10, 6, 0, 0 },
        { "YC",    &Twiss::getCO,   10, 6, 2, 0 },
        { "DX",    &Twiss::getDisp, 10, 6, 0, 0 },
        { "DPX",   &Twiss::getDisp, 10, 6, 1, 0 },

        { 0,       0,                0, 0, 0, 0 }
    };


    const ColDesc *findCol(const Twiss &table, const std::string &colName) {
        for(const ColDesc *col = allColumns; col->colName; ++col) {
            if(colName == col->colName) {
                return col;
            }
        }

        throw OpalException("Twiss::findCol()",
                            "Twiss table \"" + table.getOpalName() +
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
        Column(const Twiss &tab, const std::string &colName, const ColDesc &desc);

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
        const Twiss &itsTable;

        // Column name.
        std::string colName;

        // The function returning the column value.
        double(Twiss::*get)(const Twiss::Row &, int, int) const;

        // The indices to be transmitted to get().
        int ind_1, ind_2;
    };


    // Implementation.
    // ------------------------------------------------------------------------

    Column::Column(const Twiss &tab, const std::string &colName, const ColDesc &desc):
        itsTable(tab), colName(colName),
        get(desc.get), ind_1(desc.ind_1), ind_2(desc.ind_2)
    {}


    Column::Column(const Column &rhs):
        Scalar<double>(rhs),
        itsTable(rhs.itsTable), colName(rhs.colName),
        get(rhs.get), ind_1(rhs.ind_1), ind_2(rhs.ind_2)
    {}


    Column::~Column()
    {}


    Expressions::Scalar<double> *Column::clone() const {
        return new Column(*this);
    }


    double Column::evaluate() const {
        return (itsTable.*get)(itsTable.getCurrent(), ind_1, ind_2);
    }


    void Column::print(std::ostream &os, int) const {
        os << colName;
    }

};


// Class Twiss::Row
// ------------------------------------------------------------------------

Twiss::Row::Row(ElementBase *elem, int occur):
    FlaggedElmPtr(elem) {
    setCounter(occur);
}


Twiss::Row::Row(const FlaggedElmPtr &rhs):
    FlaggedElmPtr(rhs)
{}


Twiss::Row::~Row()
{}


const FVector<double, 6> &Twiss::Row::getCO() const {
    return orbit;
}


const FMatrix<double, 6, 6> &Twiss::Row::getMatrix() const {
    return matrix;
}


double Twiss::Row::getS() const {
    return arc;
}


double Twiss::Row::getMUi(int i) const {
    return mu[i];
}


// Class Twiss
// ------------------------------------------------------------------------

Twiss::Twiss(int size, const char *name, const char *help):
    Table(size, name, help), itsTable(0), itsMapper(0) {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "The beam line use for filling");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "The beam to be used", "UNNAMED_BEAM");
    itsAttr[RANGE] = Attributes::makeRange
                     ("RANGE", "The range in the lattice");
    itsAttr[STATIC] = Attributes::makeBool
                      ("STATIC", "If true, the table is not recalculated at each iteration");
    itsAttr[ORDER] = Attributes::makeReal
                     ("ORDER", "The order of the calculation", 2);
    itsAttr[METHOD] = Attributes::makeString
                      ("METHOD", "the algorithm used for filling:\n"
                       "\t\t\t\"LINEAR\", \"THIN\", \"THICK\", or \"TRANSPORT\"\n"
                       "\t\t\tDefault value is \"LINEAR\"", "LINEAR");
    itsAttr[REVBEAM] = Attributes::makeBool
                       ("REVBEAM", "Set true to run beam backwards through lattice");
    itsAttr[REVTRACK] = Attributes::makeBool
                        ("REVTRACK", "Set true to track against the beam");

    // READ ONLY.
    itsAttr[BETXMAX] = Attributes::makeReal
                       ("BETXMAX", "Maximum horizontal beta in m");
    itsAttr[BETXMAX].setReadOnly(true);

    itsAttr[BETYMAX] = Attributes::makeReal
                       ("BETYMAX", "Maximum vertical beta in m");
    itsAttr[BETYMAX].setReadOnly(true);

    itsAttr[XCMAX] = Attributes::makeReal
                     ("XCMAX", "Maximum horizontal closed orbit in m");
    itsAttr[XCMAX].setReadOnly(true);

    itsAttr[YCMAX] = Attributes::makeReal
                     ("YCMAX", "Maximum vertical closed orbit in m");
    itsAttr[YCMAX].setReadOnly(true);

    itsAttr[XCRMS] = Attributes::makeReal
                     ("XCRMS", "R.m.s. horizontal closed orbit in m");
    itsAttr[XCRMS].setReadOnly(true);

    itsAttr[YCRMS] = Attributes::makeReal
                     ("YCRMS", "R.m.s. vertical closed orbit in m");
    itsAttr[YCRMS].setReadOnly(true);

    itsAttr[DXMAX] = Attributes::makeReal
                     ("DXMAX", "Maximum horizontal dispersion in m");
    itsAttr[DXMAX].setReadOnly(true);

    itsAttr[DYMAX] = Attributes::makeReal
                     ("DYMAX", "Maximum vertical dispersion in m");
    itsAttr[DYMAX].setReadOnly(true);

    itsAttr[DXRMS] =  Attributes::makeReal
                      ("DXRMS", "R.m.s. horizontal dispersion in m");
    itsAttr[DXRMS].setReadOnly(true);

    itsAttr[DYRMS] = Attributes::makeReal
                     ("DYRMS", "R.m.s. vertical dispersion in m");
    itsAttr[DYRMS].setReadOnly(true);
}


Twiss::Twiss(const std::string &name, Twiss *parent):
    Table(name, parent), itsTable(new TLine(name)), itsMapper(0)
{}


Twiss::~Twiss() {
    delete itsTable;
    delete itsMapper;
}

/*
void Twiss::doomPut(DoomWriter &writer) const {
    // Write the definition for this table.
    Object::doomPut(writer);

    // Make sure the table is up to date.
    // Cast away const, to allow logically constant table to update.
    const_cast<Twiss *>(this)->fill();

    // Make the table header.
    const std::string &tableName = getOpalName();
    {
        static const int headKeyList[] = { 1, TABLE_HEAD };
        DoomWriter headWriter(tableName, headKeyList);
        headWriter.setParentName(itsLine);
        headWriter.setTypeName("TABLE_HEADER");

        static const char *columnName[numColumns] = {
            "ALFX", "ALFY", "BETX", "BETY", "DPX", "DPY", "DX", "DY", "MUX", "MUY", "S"
        };
        static const int columnNumber[numColumns] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        };
        for(int i = 0; i < numColumns; ++i) {
            headWriter.putInt(i, columnNumber[i]);
            headWriter.putString(i, columnName[i]);
        }

        int numRows = itsTable->size();
        headWriter.putInt(numColumns + 0, 1);               // 1 delta value.
        headWriter.putInt(numColumns + 1, numRows);         // number of rows.
        headWriter.putInt(numColumns + 3, 3);               // position at end.
        // Header object will be written when headWriter goes out of scope.
    }

    // Write the table row by row.
    {
        static const int rowKeyList[] = { 2, TABLE_BODY, 0 };
        DoomWriter rowWriter(tableName, rowKeyList);
        rowWriter.setParentName(itsLine);
        rowWriter.setTypeName("SUB_TABLE");
        int lineCount = 0;
        int realCount = 0;
        for(TLine::const_iterator row = begin(); row != end(); ++row) {
            rowWriter.putInt(lineCount,       row->getCounter());
            rowWriter.putString(lineCount,    row->getElement()->getName());
            rowWriter.putReal(realCount +  0, getALFi(*row, 0, 0));
            rowWriter.putReal(realCount +  1, getALFi(*row, 1, 0));
            rowWriter.putReal(realCount +  2, getBETi(*row, 0, 0));
            rowWriter.putReal(realCount +  3, getBETi(*row, 1, 0));
            rowWriter.putReal(realCount +  4, getDisp(*row, 1, 0));
            rowWriter.putReal(realCount +  5, getDisp(*row, 3, 0));
            rowWriter.putReal(realCount +  6, getDisp(*row, 0, 0));
            rowWriter.putReal(realCount +  7, getDisp(*row, 2, 0));
            rowWriter.putReal(realCount +  8, getMUi(*row, 0, 0));
            rowWriter.putReal(realCount +  9, getMUi(*row, 1, 0));
            rowWriter.putReal(realCount + 10, getS(*row, 0, 0));
            lineCount++;
            realCount += numColumns;
        }
        // Body object will be written when rowWriter goes out of scope.
    }
}
*/

Twiss::TLine::iterator Twiss::begin() {
    return itsTable->begin();
}


Twiss::TLine::const_iterator Twiss::begin() const {
    return itsTable->begin();
}


Twiss::TLine::iterator Twiss::end() {
    return itsTable->end();
}


Twiss::TLine::const_iterator Twiss::end() const {
    return itsTable->end();
}


void Twiss::execute() {
    //std::cerr << "In Twiss::execute()" << std::endl;
    // Find Table definition.
    itsLine = Attributes::getString(itsAttr[LINE]);
    BeamSequence *use = BeamSequence::find(itsLine);

    // Find Beam data.
    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    beam = Beam::find(beamName);
    reference = &beam->getReference();

    // Get the algorithm name.
    std::string method = Attributes::getString(itsAttr[METHOD]);

    // Make sure all is up-to-date.
    OpalData::getInstance()->update();

    // Create flat list with space for data storage.
    RangeRep range = Attributes::getRange(itsAttr[RANGE]);
    Flatten<Row> flattener(*use->fetchLine(), *itsTable, range);
    flattener.execute();

    if(itsTable->empty() && Options::warn) {
        std::cerr << "\n### Warning ### Lattice function table \""
                  << getOpalName() << "\" contains no elements.\n" << std::endl;
    } else {
        itsTable->front().setSelectionFlag(true);
        itsTable->back().setSelectionFlag(true);
    }

    // Assign the correct mapper.
    revBeam  = Attributes::getBool(itsAttr[REVBEAM]);
    revTrack = Attributes::getBool(itsAttr[REVTRACK]);
    revPath  = ( revBeam && !revTrack ) || ( !revBeam && revTrack );
    order = int(Round(Attributes::getReal(itsAttr[ORDER])));

    if(method == "THICK") {
        //std::cerr << "  method == \"THICK\"" << std::endl;
        itsMapper = new ThickMapper(*itsTable, *reference, revBeam, revTrack);
    } else if(method == "THIN") {
        //std::cerr << "  method == \"THIN\"" << std::endl;
        itsMapper = new ThinMapper(*itsTable, *reference, revBeam, revTrack);
    } else if(method == "LINEAR") {
        //std::cerr << "  method == \"LINEAR\"" << std::endl;
        itsMapper = new LinearMapper(*itsTable, *reference, revBeam, revTrack);
    } else {
        throw OpalException("Twiss::execute()",
                            "Method name \"" + method + "\" is unknown.");
    }

    // Fill the table.
    fill();

    // Set static flag.
    if(Attributes::getBool(itsAttr[STATIC])) dynamic = false;
    printTable(std::cout, getDefault());
    //std::cerr << "Leaving Twiss::execute()" << std::endl;
}


double Twiss::getCell(const PlaceRep &place, const std::string &colName) {
    Row &row = findRow(place);
    const ColDesc *col = findCol(*this, colName);
    return (this->*(col->get))(row, col->ind_1, col->ind_2);
}


Table::CellArray Twiss::getDefault() const {
    CellArray columns;
    for(const ColDesc *col = defaultColumns; col->colName; ++col) {
        Expressions::PtrToScalar<double> expr =
            new Column(*this, col->colName, *col);
        columns.push_back(Cell(expr, col->printWidth, col->printPrecision));
    }
    return columns;
}


std::vector<double>
Twiss::getColumn(const RangeRep &rng, const std::string &colName) {
    // Find proper column function.
    const ColDesc *col = findCol(*this, colName);
    RangeRep range(rng);
    range.initialize();
    std::vector<double> column;

    for(TLine::const_iterator row = begin(); row != end(); ++row) {
        range.enter(*row);
        if(range.isActive()) {
            column.push_back((this->*(col->get))(*row, col->ind_1, col->ind_2));
        }
        range.leave(*row);
    }

    return column;
}


const Twiss::Row &Twiss::getCurrent() const {
    return *current;
}


FMatrix<double, 6, 6> Twiss::getCurlyA() const {
    return curly_A;
}


FMatrix<double, 6, 6> Twiss::getCurlyA(const Row &row) const {
    return row.matrix * curly_A;
}


double Twiss::getEX() const {
    return beam->getEX();
}


double Twiss::getEY() const {
    return beam->getEY();
}


double Twiss::getET() const {
    return beam->getET();
}


double Twiss::getLength() {
    return itsTable->getElementLength();
}


const Beamline *Twiss::getLine() const {
    return itsTable;
}


FMatrix<double, 6, 6> Twiss::getMatrix(const Row &row) const {
    return row.matrix;
}


FVector<double, 6> Twiss::getOrbit(const Row &row) const {
    return row.orbit;
}


FVector<double, 6> Twiss::getOrbit() const {
    return orbit;
}


FMatrix<double, 6, 6> Twiss::getSigma(const Row &row) const {
    double E1 = getEX();
    double E2 = getEY();
    double E3 = getET();
    const FMatrix<double, 6, 6> &eigen = getCurlyA(row);
    FMatrix<double, 6, 6> sigma;
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j <= i; ++j) {
            sigma[i][j] = sigma[j][i] =
                              E1 * (eigen[i][0] * eigen[j][0] + eigen[i][1] * eigen[j][1]) +
                              E2 * (eigen[i][2] * eigen[j][2] + eigen[i][3] * eigen[j][3]) +
                              E3 * (eigen[i][4] * eigen[j][4] + eigen[i][5] * eigen[j][5]);
        }
    }
    return sigma;
}


std::vector<double>
Twiss::getRow(const PlaceRep &pos, const std::vector<std::string> &cols) {
    Row &row = findRow(pos);
    std::vector<double> result;

    if(cols.empty()) {
        // Standard column selection.
        for(const ColDesc *col = defaultColumns; col->colName != 0; ++col) {
            result.push_back((this->*col->get)(row, col->ind_1, col->ind_2));
        }
    } else {
        // User column selection.
        for(std::vector<std::string>::const_iterator iter = cols.begin();
            iter != cols.end(); ++iter) {
            const ColDesc *col = findCol(*this, *iter);
            result.push_back((this->*(col->get))(row, col->ind_1, col->ind_2));
        }
    }

    return result;
}


bool Twiss::isDependent(const std::string &name) const {
    // Test if name refers to USE attribute.
    if(itsLine == name) return true;

    // Test if name occurs in table.
    for(TLine::const_iterator row = begin(); row != end(); ++row) {
        if(row->getElement()->getName() == name) return true;
    }

    // Otherwise replacement is not required.
    return false;
}


Expressions::PtrToScalar<double>
Twiss::makeColumnExpression(const std::string &colName) const {
    const ColDesc *col = findCol(*this, colName);
    return new Column(*this, colName, *col);
}


bool Twiss::matches(Table *rhs) const {
    return dynamic_cast<Twiss *>(rhs) != 0;
}


void Twiss::printTableBody(std::ostream &os, const CellArray &cells) const {
    // Find line length.
    int lineLength = 16;
    for(CellArray::const_iterator cell = cells.begin();
        cell < cells.end(); ++cell) {
        lineLength += cell->printWidth;
    }

    // Write table header.
    os << std::string(lineLength, '-') << '\n';
    os << "Element        ";
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

    // Write table body.
    for(current = begin(); current != end(); ++current) {
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

    os << std::string(lineLength, '-') << '\n';
}


void Twiss::printTableTitle(std::ostream &os, const char *title) const {
    OpalData::getInstance()->printTitle(os);
    os << '\n' << title
       << ", LINE: " << itsAttr[LINE]
       << ", BEAM: " << itsAttr[BEAM]
       << ", RANGE: " << itsAttr[RANGE]
       << ", METHOD: " << itsAttr[METHOD]
       << ", ORDER: " << order << ".\n";
}


// Protected methods.
// ------------------------------------------------------------------------

Twiss::Row &Twiss::findRow(const PlaceRep &place) {
    PlaceRep row(place);
    row.initialize();

    for(TLine::iterator i = begin(); i != end(); ++i) {
        row.enter(*i);
        if(row.isActive()) return *i;
        row.leave(*i);
    }

    std::ostringstream os;
    os << row << std::ends;
    throw OpalException("Twiss::findRow()", "Row \"" + os.str() +
                        "\" not found in twiss table \"" + getOpalName() + "\".");
}


void Twiss::put() {
    //std::cerr << "==> In Twiss::put()..." << std::endl;
    // Initial (final) phase angles.
    static const double epsilon = 1.0e-8;
    static const LinearFun<double, 6> nrgy = LinearFun<double, 6>::makeVariable(5);
    double arc = 0.0;
    LinearMap<double, 6> map;
    itsMapper->getMap(map);
    bool staticQ = (map[5] == map[5][0] + nrgy);

    //std::cerr << " [Twiss::put] map =\n" << map << std::endl;
    //std::cerr << " [Twiss::put] curly_A =\n" << curly_A << std::endl;
    int imux = 0;
    int imuy = 0;
    int imut = 0;
    double mux = atan2(curly_A[0][1], curly_A[0][0]) / Physics::two_pi;
    double muy = atan2(curly_A[2][3], curly_A[2][2]) / Physics::two_pi;
    double mut = atan2(curly_A[4][5], curly_A[4][4]) / Physics::two_pi;
    //std::cerr << " tunes = (" << mux << ", " << muy << ", " << mut << ")" << std::endl;

    if(revPath) {
        //std::cerr << "Twiss::put(): doing reverse path ..." << std::endl;
        for(TLine::reverse_iterator row = itsTable->rbegin();
            row != itsTable->rend(); ++row) {
            // Store values at end of element to the table.
            row->orbit = map.constantTerm();
            row->matrix = map.linearTerms();
            row->arc = arc;
            row->mu[0] = mux + double(imux);
            row->mu[1] = muy + double(imuy);
            row->mu[2] = mut + double(imut);

            // Traverse element.
            row->accept(*itsMapper);

            // Update values for beginning of element.
            ElementBase &elem = *row->getElement();
            if(! dynamic_cast<Beamline *>(&elem)) {
                arc -= elem.getElementLength();
            }

            itsMapper->getMap(map);
            double A11 = 0.0;
            double A12 = 0.0;
            double A33 = 0.0;
            double A34 = 0.0;
            double A55 = 0.0;
            double A56 = 0.0;
            FMatrix<double, 6, 6> matrix = map.linearTerms();
            const double *row1 = matrix[0];
            const double *row3 = matrix[2];
            const double *row5 = matrix[4];

            for(int i = 0; i < 6; ++i) {
                A11 += row1[i] * curly_A[i][0];
                A12 += row1[i] * curly_A[i][1];
                A33 += row3[i] * curly_A[i][2];
                A34 += row3[i] * curly_A[i][3];
                A55 += row5[i] * curly_A[i][4];
                A56 += row5[i] * curly_A[i][5];
            }

            double mu = atan2(A12, A11) / Physics::two_pi;
            if(mu - epsilon > mux) --imux;
            mux = mu;
            mu = atan2(A34, A33) / Physics::two_pi;
            if(mu - epsilon > muy) --imuy;
            muy = mu;
            mu = atan2(A56, A55) / Physics::two_pi;
            if(mu - epsilon > mut) --imut;
            mut = mu;

            //        FMatrix<double,6,6> A_scr = map.linearTerms() * curly_A;
            //        double mu;
            //        mu = atan2(A_scr[0][1], A_scr[0][0]) / Physics::two_pi;
            //        if (mu - epsilon > mux) --imux;
            //        mux = mu;
            //        mu = atan2(A_scr[2][3], A_scr[2][2]) / Physics::two_pi;
            //        if (mu - epsilon > muy) --imuy;
            //        muy = mu;
            //        mu = atan2(A_scr[4][5], A_scr[4][4]) / Physics::two_pi;
            //        if (mu - epsilon > mut) --imut;
            //        mut = mu;
        }

        double arc1 = arc;
        arc = - arc;
        double mux1 = mux + double(imux);
        double muy1 = muy + double(imuy);
        double mut1 = mut + double(imut);

        // Change arc length and phases to their values from the origin.
        for(TLine::iterator row = begin(); row != end(); ++row) {
            row->arc -= arc1;
            row->mu[0] -= mux1;
            row->mu[1] -= muy1;
            row->mu[2] -= mut1;
        }
    } else {
        //std::cerr << "Twiss::put(): doing forward path ..." << std::endl;
        double arc1 = 0.0;
        double mux1 = mux;
        double muy1 = muy;
        double mut1 = mut;

        if(staticQ) {  // reset 'T' component of map
            //std::cerr << "Twiss::put(): resetting 'T' component of map ..." << std::endl;
            itsMapper->getMap(map);
            map[4][0] = 0.0;
            itsMapper->setMap(map);
        }

        for(TLine::iterator row = begin(); row != end(); ++row) {
            // Traverse element.
            row->accept(*itsMapper);

            // Update values for end of element.
            ElementBase &elem = *row->getElement();
            if(! dynamic_cast<Beamline *>(&elem)) {
                arc += elem.getElementLength();
            }

            // Store values at end of element to the table.
            itsMapper->getMap(map);
            row->orbit = map.constantTerm();
            //std::cerr << "Twiss::put(): row->orbit =\n" << row->orbit << std::endl;
            row->matrix = map.linearTerms();
            double A11 = 0.0;
            double A12 = 0.0;
            double A33 = 0.0;
            double A34 = 0.0;
            double A55 = 0.0;
            double A56 = 0.0;
            FMatrix<double, 6, 6> &matrix = row->matrix;
            const double *row1 = matrix[0];
            const double *row3 = matrix[2];
            const double *row5 = matrix[4];

            //std::cerr << " matrix =\n" << matrix << std::endl;
            //std::cerr << " curly_A =\n" << curly_A << std::endl;

            for(int i = 0; i < 6; ++i) {
                A11 += row1[i] * curly_A[i][0];
                A12 += row1[i] * curly_A[i][1];
                A33 += row3[i] * curly_A[i][2];
                A34 += row3[i] * curly_A[i][3];
                A55 += row5[i] * curly_A[i][4];
                A56 += row5[i] * curly_A[i][5];
            }

            double mu = atan2(A12, A11) / Physics::two_pi;
            if(mu + epsilon < mux) ++imux;
            mux = mu;
            mu = atan2(A34, A33) / Physics::two_pi;
            if(mu + epsilon < muy) ++imuy;
            muy = mu;
            mu = atan2(A56, A55) / Physics::two_pi;
            if(mu + epsilon < mut) ++imut;
            mut = mu;

            //        FMatrix<double,6,6> A_scr = row->matrix * curly_A;
            //        double mu;
            //        mu = atan2(A_scr[0][1], A_scr[0][0]) / Physics::two_pi;
            //        if (mu + epsilon < mux) ++imux;
            //        mux = mu;
            //        mu = atan2(A_scr[2][3], A_scr[2][2]) / Physics::two_pi;
            //        if (mu + epsilon < muy) ++imuy;
            //        muy = mu;
            //        mu = atan2(A_scr[4][5], A_scr[4][4]) / Physics::two_pi;
            //        if (mu + epsilon < mut) ++imut;
            //        mut = mu;

            row->arc = arc - arc1;
            row->mu[0] = mux + double(imux) - mux1;
            row->mu[1] = muy + double(imuy) - muy1;
            row->mu[2] = mut + double(imut) - mut1;
        }
    }

    // Compute statistic quantities.
    double betxmax = 0.0;
    double betymax = 0.0;
    double xmax = 0.0;
    double ymax = 0.0;
    double xrms = 0.0;
    double yrms = 0.0;
    double dxmax = 0.0;
    double dymax = 0.0;
    double dxrms = 0.0;
    double dyrms = 0.0;

    for(TLine::iterator row = begin(); row != end(); ++row) {
        betxmax = std::max(std::abs(getBETi(*row, 0, 0)), betxmax);
        betymax = std::max(std::abs(getBETi(*row, 1, 0)), betymax);
        double x = getCO(*row, 0, 0);
        double y = getCO(*row, 2, 0);
        double dx = getDisp(*row, 0, 0);
        double dy = getDisp(*row, 2, 0);
        xmax = std::max(std::abs(x), xmax);
        ymax = std::max(std::abs(y), ymax);
        xrms += x * x;
        yrms += y * y;
        dxmax = std::max(std::abs(dx), dxmax);
        dymax = std::max(std::abs(dy), dymax);
        dxrms += dx * dx;
        dyrms += dy * dy;
    }

    double size = itsTable->size();
    Attributes::setReal(itsAttr[BETXMAX], betxmax);
    Attributes::setReal(itsAttr[BETYMAX], betymax);
    Attributes::setReal(itsAttr[XCMAX], xmax);
    Attributes::setReal(itsAttr[YCMAX], ymax);
    Attributes::setReal(itsAttr[DXMAX], dxmax);
    Attributes::setReal(itsAttr[DYMAX], dymax);
    Attributes::setReal(itsAttr[XCRMS], sqrt(xrms / size));
    Attributes::setReal(itsAttr[YCRMS], sqrt(yrms / size));
    Attributes::setReal(itsAttr[DXRMS], sqrt(dxrms / size));
    Attributes::setReal(itsAttr[DYRMS], sqrt(dyrms / size));
    //  std::cerr << "==> Leaving Twiss::put()" << std::endl;
}


double Twiss::getS(const Twiss::Row &row, int, int) const {
    return row.getS();
}


double Twiss::getMUi(const Twiss::Row &row, int i1, int) const {
    return row.getMUi(i1);
}


double Twiss::getBETi(const Twiss::Row &row, int i1, int) const {
    const double *row1 = row.matrix[2*i1];
    const double *row2 = row.matrix[2*i1+1];
    double r11 = 0.0;
    double r12 = 0.0;
    double r21 = 0.0;
    double r22 = 0.0;

    for(int i = 0; i < 6; ++i) {
        double t1 = row1[i];
        r11 += t1 * curly_A[i][2*i1];
        r12 += t1 * curly_A[i][2*i1+1];
        double t2 = row2[i];
        r21 += t2 * curly_A[i][2*i1];
        r22 += t2 * curly_A[i][2*i1+1];
    }

    //    const FMatrix<double,6,6> eigen = row.matrix * curly_A;
    //    const double r11 = eigen[2*i1][2*i1];
    //    const double r12 = eigen[2*i1][2*i1+1];
    //    const double r21 = eigen[2*i1+1][2*i1];
    //    const double r22 = eigen[2*i1+1][2*i1+1];
    return (r11 * r11 + r12 * r12) / (r11 * r22 - r12 * r21);
}


double Twiss::getALFi(const Twiss::Row &row, int i1, int) const {
    const double *row1 = row.matrix[2*i1];
    const double *row2 = row.matrix[2*i1+1];
    double r11 = 0.0;
    double r12 = 0.0;
    double r21 = 0.0;
    double r22 = 0.0;

    for(int i = 0; i < 6; ++i) {
        double t1 = row1[i];
        r11 += t1 * curly_A[i][2*i1];
        r12 += t1 * curly_A[i][2*i1+1];
        double t2 = row2[i];
        r21 += t2 * curly_A[i][2*i1];
        r22 += t2 * curly_A[i][2*i1+1];
    }

    // ada Mon Mar 27 23:45:26 CEST 2000
    //    const FMatrix<double,6,6> eigen = row.matrix * curly_A;
    //    const double r11 = eigen[2*i1][2*i1];
    //    const double r12 = eigen[2*i1][2*i1+1];
    //    const double r21 = eigen[2*i1+1][2*i1];
    //    const double r22 = eigen[2*i1+1][2*i1+1];
    return - (r11 * r21 + r12 * r22) / (r11 * r22 - r12 * r21);
}


double Twiss::getBETik(const Twiss::Row &row, int i1, int i2) const {
    const double *row1 = row.matrix[2*i1];
    double r11 = 0.0;
    double r12 = 0.0;

    for(int i = 0; i < 6; ++i) {
        double t1 = row1[i];
        r11 += t1 * curly_A[i][2*i2];
        r12 += t1 * curly_A[i][2*i2+1];
    }

    // ada Mon Mar 27 23:45:26 CEST 2000
    //    const FMatrix<double,6,6> eigen = row.matrix * curly_A;
    //    const double r11 = eigen[2*i1][2*i2];
    //    const double r12 = eigen[2*i1][2*i2+1];
    return (r11 * r11 + r12 * r12);
}


double Twiss::getALFik(const Twiss::Row &row, int i1, int i2) const {
    const double *row1 = row.matrix[2*i1];
    const double *row2 = row.matrix[2*i1+1];
    double r11 = 0.0;
    double r12 = 0.0;
    double r21 = 0.0;
    double r22 = 0.0;

    for(int i = 0; i < 6; ++i) {
        double t1 = row1[i];
        r11 += t1 * curly_A[i][2*i2];
        r12 += t1 * curly_A[i][2*i2+1];
        double t2 = row2[i];
        r21 += t2 * curly_A[i][2*i2];
        r22 += t2 * curly_A[i][2*i2+1];
    }

    // ada Mon Mar 27 23:45:26 CEST 2000
    //  const FMatrix<double,6,6> eigen = row.matrix * curly_A;
    //    const double r11 = eigen[2*i1][2*i2];
    //    const double r12 = eigen[2*i1][2*i2+1];
    //    const double r21 = eigen[2*i1+1][2*i2];
    //    const double r22 = eigen[2*i1+1][2*i2+1];
    return (r11 * r21 + r12 * r22);
}


double Twiss::getGAMik(const Twiss::Row &row, int i1, int i2) const {
    const double *row1 = row.matrix[2*i1+1];
    double r21 = 0.0;
    double r22 = 0.0;

    for(int i = 0; i < 6; ++i) {
        double t1 = row1[i];
        r21 += t1 * curly_A[i][2*i2];
        r22 += t1 * curly_A[i][2*i2+1];
    }

    // ada Mon Mar 27 23:45:26 CEST 2000
    //  const FMatrix<double,6,6> eigen = row.matrix * curly_A;
    //  const double r21 = eigen[2*i1+1][2*i2];
    //  const double r22 = eigen[2*i1+1][2*i2+1];
    return (r21 * r21 + r22 * r22);
}


double Twiss::getCO(const Twiss::Row &row, int i1, int) const {
    return row.getCO()[i1];
}


double Twiss::getDisp(const Twiss::Row &row, int i1, int) const {
    //  FMatrix<double,6,6> matrix = row.matrix;
    const double *row1 = row.matrix[i1];
    double result = 0.0;

    for(int i = 0; i < 6; ++i) {
        //    result += matrix[i1][i] * curly_A[i][5];
        result += row1[i] * curly_A[i][5];
    }

    return result;
}


double Twiss::getEigen(const Twiss::Row &row, int i1, int i2) const {
    //  FMatrix<double,6,6> matrix = row.matrix;
    const double *row1 = row.matrix[i1];
    double result = 0.0;

    for(int i = 0; i < 6; ++i) {
        //    result += matrix[i1][i] * curly_A[i][i2];
        result += row1[i] * curly_A[i][i2];
    }

    return result;
}

/*
1334    double Twiss::getEigen(const Twiss::Row &row, int i1, int i2) const
1335    {
1336      FMatrix<double,6,6> matrix = row.matrix;
1337      const double *row1 = row.matrix[i1];
1338      double result = 0.0;
1339
1340      for (int i = 0; i < 6; ++i) {
1341         result += matrix[i1][i] * curly_A[i][i2];
1342        //result += row1[i] * curly_A[i][i2];
1343      }
1344
1345      return result;
1346    }
*/




double Twiss::getSigma(const Twiss::Row &row, int i1, int i2) const {
    const FMatrix<double, 6, 6> eigen = row.matrix * curly_A;
    const double E1 = getEX();
    const double E2 = getEY();
    const double E3 = getET();
    return
        (E1 * (eigen[i1][0] * eigen[i2][0] + eigen[i1][1] * eigen[i2][1]) +
         E2 * (eigen[i1][2] * eigen[i2][2] + eigen[i1][3] * eigen[i2][3]) +
         E3 * (eigen[i1][4] * eigen[i2][4] + eigen[i1][5] * eigen[i2][5]));
}


double Twiss::getMatrix(const Twiss::Row &row, int i1, int i2) const {
    return row.matrix[i1][i2];
}