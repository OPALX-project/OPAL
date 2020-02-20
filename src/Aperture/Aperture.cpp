#include "Aperture/Aperture.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Fields/BMultipoleField.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/Expressions.h"
#include "Utilities/OpalException.h"
#include "Elements/OpalElement.h"
#include "AbsBeamline/Multipole.h"
#include "Physics/Physics.h"
#include "AbstractObjects/RangeRep.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Cyclotron.h"
#include "Utilities/Timer.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;

// Local structures.
// ------------------------------------------------------------------------

namespace {

    //return the inverse of the "alpha,beta,gamma" transformation matrix

    FMatrix<double, 3, 3> invers(double c, double c_prim, double s, double s_prim) {
        FMatrix<double, 3, 3> M;
        double det = -pow(c_prim, 3) * pow(s, 3) + 3 * c * s_prim * pow(c_prim, 2)
            * pow(s, 2) - 3 * pow(c, 2) * c_prim * s * pow(s_prim, 2) + pow(c, 3) * pow(s_prim, 3);

        M(0, 0) = (-c_prim * s * pow(s_prim, 2) + c * pow(s_prim, 3)) / det;
        M(0, 1) = (-2 * c_prim * pow(s, 2) * s_prim + 2 * c * s * pow(s_prim, 2)) / det;
        M(0, 2) = (-c_prim * pow(s, 3) + c * pow(s, 2) * s_prim) / det;
        M(1, 0) = (-pow(c_prim, 2) * s * s_prim + c * c_prim * pow(s_prim, 2)) / det;
        M(1, 1) = (-pow(c_prim, 2) * pow(s, 2) + pow(c, 2) * pow(s_prim, 2)) / det;
        M(1, 2) = (-c * c_prim * pow(s, 2) + pow(c, 2) * s * s_prim) / det;
        M(2, 0) = (-pow(c_prim, 3) * s + c * pow(c_prim, 2) * s_prim) / det;
        M(2, 1) = (-2 * c * pow(c_prim, 2) * s + 2 * pow(c, 2) * c_prim * s_prim) / det;
        M(2, 2) = (-pow(c, 2) * c_prim * s + pow(c, 3) * s_prim) / det;

        return M;
    }

    struct ColDesc {

        // Column name.
        const char *colName;

        // Pointer to the row method which returns the column value.
        double(Aperture::*get)(const Aperture::A_row &, int, int) const;

        // Default field width and precision.
        int printWidth, printPrecision;

        // Indices to be given to get().
        int ind_1, ind_2;
    };


    // The complete column entries table.
    const ColDesc allColumns[] = {

        // "Naive" maximum Twiss.
        { "BETXMAX",  &Aperture::getBETXMAX,  10, 6, 0, 0 },
        { "BETYMAX",  &Aperture::getBETYMAX,  10, 6, 0, 0 },
        { "APERTMIN", &Aperture::getAPERTMIN, 10, 6, 0, 0},
        //{ "ALFY",  &Aperture::getALFi,  10, 6, 1, 0 },
        //{ "ALFT",  &Aperture::getALFi,  10, 6, 2, 0 },
        { 0,       0,                 0, 0, 0, 0 }
    };


    // The default column entries table.
    const ColDesc defaultColumns[] = {
        // "Naive" maximum Twiss.
        { "BETXMAX",  &Aperture::getBETXMAX,  10, 6, 0, 0 },
        { "BETYMAX",  &Aperture::getBETYMAX,  10, 6, 0, 0 },
        { "APERTMIN", &Aperture::getAPERTMIN, 10, 6, 0, 0},
        { 0,       0,                0, 0, 0, 0 }
    };

    const ColDesc *findCol(const Aperture &table, const std::string &colName) {
        for(const ColDesc *col = allColumns; col->colName; ++col) {
            if(colName == col->colName) {
                return col;
            }
        }

        throw OpalException("Aperture::findCol()",
                            "Aperture table \"" + table.getOpalName() +
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
        Column(const Aperture &tab, const std::string &colName, const ColDesc &desc);

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
        const Aperture &itsTable;

        // Column name.
        std::string colName;

        // The function returning the column value.
        double(Aperture::*get)(const Aperture::A_row &, int, int) const;

        // The indices to be transmitted to get().
        int ind_1, ind_2;
    };


    // Implementation.
    // ------------------------------------------------------------------------

    Column::Column(const Aperture &tab, const std::string &colName, const ColDesc &desc):
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

Aperture::Aperture():
    DefaultVisitor(itsTable, false, false),
    Table(SIZE, "APERTURE", "help"), itsTable() {
    itsAttr[TABLE] = Attributes::makeString
                     ("TABLE", "Name of the TWISS table to be used.");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "The beam to be used", "UNNAMED_BEAM");
    itsAttr[NSLICE] = Attributes::makeReal
                      ("NSLICE", "the number of slices of the interpolation inside the optic element", 10);
    itsAttr[STATIC] = Attributes::makeBool
                      ("STATIC", "recalculation if static equal false", true);
    itsAttr[DATA] = Attributes::makeRealArray
                    ("DATA", "The data needed for aperture calculation(co,deltap,BBeat,HALO)");
    itsAttr[DEFAULTAPERTURE] = Attributes::makeRealArray
                               ("DEFAULTAPERTURE", "The default beam screen for markers and drift generated in sequences");
    itsAttr[FILE] = Attributes::makeString
                    ("FILE", "Name of file to receive APERTURE output", "APERTURE.dat");

    registerOwnership(AttributeHandler::STATEMENT);
}
Aperture::Aperture(const std::string &name, Aperture *parent):
    DefaultVisitor(itsTable, false, false),
    Table(name, parent), itsTable(name)
{}

Aperture::A_row::A_row(const A_row &a): FlaggedElmPtr(a), Type_elm(a.Type_elm),
    Orb(a.Orb), Interpol(a.Interpol) {}


inline double Aperture::A_row::getBeta_x(int ind) {
    return Interpol[ind].Beta_x;
}
inline double Aperture::A_row::getBeta_y(int ind) {
    return Interpol[ind].Beta_y;
}

inline double Aperture::A_row::getDisp_x(int ind) {
    return Interpol[ind].Disp_x;
}
inline double Aperture::A_row::getDisp_x_prim(int ind) {
    return Interpol[ind].Disp_x_prim;
}
inline double Aperture::A_row::getDisp_y(int ind) {
    return Interpol[ind].Disp_y;
}
inline double Aperture::A_row::getDisp_y_prim(int ind) {
    return Interpol[ind].Disp_y_prim;
}
inline double Aperture::A_row::getApert(int ind) {
    return Interpol[ind].apert;
}
inline std::string Aperture::A_row::getType_elm() {
    return Type_elm;
}
inline double Aperture::A_row::getOrb() {
    return Orb;
}

//iterator used for the displacement in the list of optic elments

Aperture::A_Tline::iterator Aperture::begin() {
    return itsTable.begin();
}
Aperture::A_Tline::const_iterator Aperture::begin() const {
    return itsTable.begin();
}
Aperture::A_Tline::iterator Aperture::end() {
    return itsTable.end();
}
Aperture::A_Tline::const_iterator Aperture::end() const {
    return itsTable.end();
}

double Aperture::getBETXMAX(const A_row &row, int /*i1*/, int /*i2*/) const {
    double max;
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    max = row.Interpol[0].Beta_x;
    for(int i = 1; i < nslice; ++i) {
        if(max < row.Interpol[i].Beta_x) max = row.Interpol[i].Beta_x;
    }
    return max;
}

double Aperture::getBETYMAX(const A_row &row, int /*i1*/, int /*i2*/) const {
    double max;
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    max = row.Interpol[0].Beta_y;
    for(int i = 1; i < nslice; ++i) {
        if(max < row.Interpol[i].Beta_y) max = row.Interpol[i].Beta_y;

    }
    return max;
}

double Aperture::getAPERTMIN(const A_row &row, int /*i1*/, int /*i2*/) const {
    double min;
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    min = row.Interpol[0].apert;
    for(int i = 1; i < nslice; ++i) {
        if(min > row.Interpol[i].apert) min = row.Interpol[i].apert;

    }
    return min;
}


void Aperture::visitMultipole(const Multipole &m) {
    BMultipoleField f = m.getField();
    double K = f.normal(2) * 0.299792458;

    const PartData &p = beam->getReference();

    data.kq = K / (p.getP() / 1e9);
    data.k0s = 0;
    data.l = m.getElementLength();
    data.phi = 0;

    data.courb = 0;
    data.type = "multi";
}

void Aperture::visitRBend(const RBend &rb) {
    const RBendGeometry &geom = rb.getGeometry();
    data.l = geom.getElementLength();
    if(data.l != 0) {
        data.phi = geom.getBendAngle();
        data.e1 = rb.getEntryFaceRotation();
        data.e2 = rb.getExitFaceRotation();
        data.type = "rbend";
        data.courb = 2 * sin(data.phi / 2) / data.l;
        data.l = data.phi / data.courb;

        BMultipoleField f = rb.getField();
        double K = f.normal(2) * 0.299792458;

        const PartData &p = beam->getReference();


        data.kq = K / (p.getP() / 1e9);
        data.k0s = f.skew(1) * 0.299792458 / (p.getP() / 1e9);
    }
}

void Aperture::visitCyclotron(const Cyclotron &/*cy*/) {
    ERRORMSG("MSplit::visitCyclotron(const Cyclotron &cy) not implemented");
}


void Aperture::visitSBend(const SBend &sb) {
    const PlanarArcGeometry &geom = sb.getGeometry();
    data.l = geom.getElementLength();
    if(data.l != 0) {
        data.phi = geom.getBendAngle();
        data.e1 = sb.getEntryFaceRotation();
        data.e2 = sb.getExitFaceRotation();
        data.type = "sbend";
        data.courb = data.phi / data.l;

        BMultipoleField f = sb.getField();
        double K = f.normal(2) * 0.299792458;
        const PartData &p = beam->getReference();

        data.kq = K / (p.getP() / 1e9);
        data.k0s = f.skew(1) * 0.299792458 / (p.getP() / 1e9);

    }

}

void Aperture::applyDefault(const ElementBase &eb) {
    data.type = "drift";
    data.l = eb.getElementLength();
    data.courb = 0;
    data.phi = 0;
    data.kq = 0;
    data.k0s = 0;
}


void Aperture::calcul(Twiss::TLine::iterator i, A_row &a, int nslice, Twiss *tp) {
    double Kx, Ky;

    if(data.l != 0) {
        FMatrix<double, 3, 3> Mx;
        FMatrix<double, 3, 3> My;

        Kx = data.kq + pow(data.courb, 2);
        Ky = -data.kq + pow(data.k0s, 2);

        a.Interpol[nslice-1].Beta_x = tp->getBETi(*i, 0, 0);
        a.Interpol[nslice-1].Beta_y = tp->getBETi(*i, 1, 0);
        a.Interpol[nslice-1].Disp_x = tp->getDisp(*i, 0, 0);
        a.Interpol[nslice-1].Disp_y = tp->getDisp(*i, 2, 0);
        a.Interpol[nslice-1].Disp_x_prim = tp->getDisp(*i, 1, 0);
        a.Interpol[nslice-1].Disp_y_prim = tp->getDisp(*i, 3, 0);

        Euler_x(0) = tp->getBETi(*i, 0, 0);
        Euler_x(1) = tp->getALFi(*i, 0, 0);
        Euler_x(2) = (1 + pow(Euler_x(1), 2)) / Euler_x(0);
        Euler_y(0) = tp->getBETi(*i, 1, 0);
        Euler_y(1) = tp->getALFi(*i, 1, 0);
        Euler_y(2) = (1 + pow(Euler_y(1), 2)) / Euler_y(0);

        Dispx(0) = tp->getDisp(*i, 0, 0);
        Dispx(1) = tp->getDisp(*i, 1, 0);
        Dispy(0) = tp->getDisp(*i, 2, 0);
        Dispy(1) = tp->getDisp(*i, 3, 0);

        if(data.type.compare("rbend") == 0) {
            //propagation of x twiss function through the exit face
            Transf_mat(0, 0) = 1;
            Transf_mat(0, 1) = 0;
            Transf_mat(1, 0) = data.courb * tan(data.e2 + data.phi / 2);
            Transf_mat(1, 1) = 1;
            Transf_mat(0, 2) = 0;
            Transf_mat(1, 2) = 0;

            Mx = invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

            Euler_x(0) = Mx(0, 0) * tp->getBETi(*i, 0, 0) + Mx(0, 1) * tp->getALFi(*i, 0, 0) + Mx(0, 2) *
                         (1 + pow(tp->getALFi(*i, 0, 0), 2)) / tp->getBETi(*i, 0, 0);
            Euler_x(1) = Mx(1, 0) * tp->getBETi(*i, 0, 0) + Mx(1, 1) * tp->getALFi(*i, 0, 0) + Mx(1, 2) *
                         (1 + pow(tp->getALFi(*i, 0, 0), 2)) / tp->getBETi(*i, 0, 0);
            Euler_x(2) = (1 + pow(Euler_x(1), 2)) / Euler_x(0);

            Dispx(0) = (Transf_mat(1, 1) * (tp->getDisp(*i, 0, 0) - Transf_mat(0, 2)) +
                        Transf_mat(0, 1) * (Transf_mat(1, 2) - tp->getDisp(*i, 1, 0))) /
                       (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
            Dispx(1) = (Transf_mat(0, 0) * (tp->getDisp(*i, 1, 0) - Transf_mat(1, 2)) +
                        Transf_mat(1, 0) * (Transf_mat(0, 2) - tp->getDisp(*i, 0, 0))) /
                       (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));

            //propagation of y twiss function through the exit face

            Transf_mat(3, 3) = 1;
            Transf_mat(3, 4) = 0;
            Transf_mat(4, 3) = -data.courb * tan(data.e2 + data.phi / 2);
            Transf_mat(4, 4) = 1;
            Transf_mat(3, 5) = 0;
            Transf_mat(4, 5) = 0;

            My = invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));

            Euler_y(0) = My(0, 0) * tp->getBETi(*i, 1, 0) + My(0, 1) * tp->getALFi(*i, 1, 0) + My(0, 2) *
                         (1 + pow(tp->getALFi(*i, 1, 0), 2)) / tp->getBETi(*i, 1, 0);
            Euler_y(1) = My(1, 0) * tp->getBETi(*i, 1, 0) + My(1, 1) * tp->getALFi(*i, 1, 0) + My(1, 2) *
                         (1 + pow(tp->getALFi(*i, 1, 0), 2)) / tp->getBETi(*i, 1, 0);
            Euler_y(2) = (1 + pow(Euler_y(1), 2)) / Euler_y(0);

            Dispy(0) = (Transf_mat(4, 4) * (tp->getDisp(*i, 2, 0) - Transf_mat(3, 5)) +
                        Transf_mat(3, 4) * (Transf_mat(4, 5) - tp->getDisp(*i, 3, 0))) /
                       (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
            Dispy(1) = (Transf_mat(3, 3) * (tp->getDisp(*i, 3, 0) - Transf_mat(4, 5)) +
                        Transf_mat(4, 3) * (Transf_mat(3, 5) - tp->getDisp(*i, 2, 0))) /
                       (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));

        } else if(data.type.compare("sbend") == 0) {
            //propagation of x twiss function through the exit face
            Transf_mat(0, 0) = 1;
            Transf_mat(0, 1) = 0;
            Transf_mat(1, 0) = data.courb * tan(data.e2);
            Transf_mat(1, 1) = 1;
            Transf_mat(0, 2) = 0;
            Transf_mat(1, 2) = 0;

            Mx = invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

            Euler_x(0) = Mx(0, 0) * tp->getBETi(*i, 0, 0) + Mx(0, 1) * tp->getALFi(*i, 0, 0) + Mx(0, 2) *
                         (1 + pow(tp->getALFi(*i, 0, 0), 2)) / tp->getBETi(*i, 0, 0);
            Euler_x(1) = Mx(1, 0) * tp->getBETi(*i, 0, 0) + Mx(1, 1) * tp->getALFi(*i, 0, 0) + Mx(1, 2) *
                         (1 + pow(tp->getALFi(*i, 0, 0), 2)) / tp->getBETi(*i, 0, 0);
            Euler_x(2) = (1 + pow(Euler_x(1), 2)) / Euler_x(0);

            Dispx(0) = (Transf_mat(1, 1) * (tp->getDisp(*i, 0, 0) - Transf_mat(0, 2)) +
                        Transf_mat(0, 1) * (Transf_mat(1, 2) - tp->getDisp(*i, 1, 0))) /
                       (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
            Dispx(1) = (Transf_mat(0, 0) * (tp->getDisp(*i, 1, 0) - Transf_mat(1, 2)) +
                        Transf_mat(1, 0) * (Transf_mat(0, 2) - tp->getDisp(*i, 0, 0))) /
                       (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));

            //propagation of y twiss function through the exit face

            Transf_mat(3, 3) = 1;
            Transf_mat(3, 4) = 0;
            Transf_mat(4, 3) = -data.courb * tan(data.e2);
            Transf_mat(4, 4) = 1;
            Transf_mat(3, 5) = 0;
            Transf_mat(4, 5) = 0;

            My = invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));

            Euler_y(0) = My(0, 0) * tp->getBETi(*i, 1, 0) + My(0, 1) * tp->getALFi(*i, 1, 0) + My(0, 2) *
                         (1 + pow(tp->getALFi(*i, 1, 0), 2)) / tp->getBETi(*i, 1, 0);
            Euler_y(1) = My(1, 0) * tp->getBETi(*i, 1, 0) + My(1, 1) * tp->getALFi(*i, 1, 0) + My(1, 2) *
                         (1 + pow(tp->getALFi(*i, 1, 0), 2)) / tp->getBETi(*i, 1, 0);
            Euler_y(2) = (1 + pow(Euler_y(1), 2)) / Euler_y(0);

            Dispy(0) = (Transf_mat(4, 4) * (tp->getDisp(*i, 2, 0) - Transf_mat(3, 5)) +
                        Transf_mat(3, 4) * (Transf_mat(4, 5) - tp->getDisp(*i, 3, 0))) /
                       (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
            Dispy(1) = (Transf_mat(3, 3) * (tp->getDisp(*i, 3, 0) - Transf_mat(4, 5)) +
                        Transf_mat(4, 3) * (Transf_mat(3, 5) - tp->getDisp(*i, 2, 0))) /
                       (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));

        }

        if(abs(Kx) < 1e-6) {
            for(int j = 1; j <= nslice - 1; ++j) {

                //expansion of the transformation matrix coefficients at third order
                Transf_mat(0, 0) = 1 - Kx * pow(data.l * j / nslice, 2) / 2;
                Transf_mat(0, 1) = data.l * j / nslice * (1 - Kx / 6 * pow(data.l * j / nslice, 2));
                Transf_mat(1, 0) = -Kx * Transf_mat(0, 1);
                Transf_mat(1, 1) = Transf_mat(0, 0);
                Transf_mat(0, 2) = data.courb / 2 * pow(data.l * j / nslice, 2);
                Transf_mat(1, 2) = data.courb * data.l * j / nslice * (1 - Kx / 6 * pow(data.l * j / nslice, 2));

                Mx = invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

                //Filling of the vectors containing Beta et Disp inside the element

                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);

                a.Interpol[nslice-j-1].Disp_x = (Transf_mat(1, 1) * (Dispx(0) - Transf_mat(0, 2)) +
                                                 Transf_mat(0, 1) * (Transf_mat(1, 2) - Dispx(1))) /
                                                (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
                a.Interpol[nslice-j-1].Disp_x_prim = (Transf_mat(0, 0) * (Dispx(1) - Transf_mat(1, 2)) +
                                                      Transf_mat(1, 0) * (Transf_mat(0, 2) - Dispx(0))) /
                                                     (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));

            }
        } else if(Kx > 0) {
            for(int j = 1; j <= nslice - 1; ++j) {
                Transf_mat(0, 0) = cos(data.l / nslice * j * sqrt(Kx));
                Transf_mat(0, 1) = sin(data.l / nslice * j * sqrt(Kx)) / sqrt(Kx);
                Transf_mat(1, 0) = -Kx * Transf_mat(0, 1);
                Transf_mat(1, 1) = Transf_mat(0, 0);
                Transf_mat(0, 2) = data.courb * (1 - cos(data.l / nslice * j * sqrt(Kx))) / Kx;
                Transf_mat(1, 2) = data.courb * sin(data.l / nslice * j * sqrt(Kx)) / sqrt(Kx);

                Mx = invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));
                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);

                a.Interpol[nslice-j-1].Disp_x = (Transf_mat(1, 1) * (Dispx(0) - Transf_mat(0, 2)) +
                                                 Transf_mat(0, 1) * (Transf_mat(1, 2) - Dispx(1))) /
                                                (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
                a.Interpol[nslice-j-1].Disp_x_prim = (Transf_mat(0, 0) * (Dispx(1) - Transf_mat(1, 2)) +
                                                      Transf_mat(1, 0) * (Transf_mat(0, 2) - Dispx(0))) /
                                                     (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
            }
        }

        else { //if(Kx<0)
            for(int j = 1; j <= nslice - 1; ++j) {
                Transf_mat(0, 0) = cosh(data.l / nslice * j * sqrt(-Kx));
                Transf_mat(0, 1) = sinh(data.l / nslice * j * sqrt(-Kx)) / sqrt(-Kx);
                Transf_mat(1, 0) = -Kx * Transf_mat(0, 1);
                Transf_mat(1, 1) = Transf_mat(0, 0);
                Transf_mat(0, 2) = data.courb * (1 - cosh(data.l / nslice * j * sqrt(-Kx))) / -Kx;
                Transf_mat(1, 2) = data.courb * sinh(data.l / nslice * j * sqrt(-Kx)) / sqrt(-Kx);

                Mx = invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));
                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);

                a.Interpol[nslice-j-1].Disp_x = (Transf_mat(1, 1) * (Dispx(0) - Transf_mat(0, 2)) +
                                                 Transf_mat(0, 1) * (Transf_mat(1, 2) - Dispx(1))) /
                                                (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
                a.Interpol[nslice-j-1].Disp_x_prim = (Transf_mat(0, 0) * (Dispx(1) - Transf_mat(1, 2)) +
                                                      Transf_mat(1, 0) * (Transf_mat(0, 2) - Dispx(0))) /
                                                     (Transf_mat(1, 1) * Transf_mat(0, 0) - Transf_mat(1, 0) * Transf_mat(0, 1));
            }
        }
        if(abs(Ky) < 1e-6) {
            for(int j = 1; j <= nslice - 1; ++j) {
                //expansion of the transformation matrix coefficients at third order

                Transf_mat(3, 3) = 1 - Ky * pow(data.l * j / nslice, 2) / 2;
                Transf_mat(3, 4) = data.l * j / nslice * (1 - Ky / 6 * pow(data.l * j / nslice, 2));
                Transf_mat(4, 3) = -Ky * Transf_mat(3, 4);
                Transf_mat(4, 4) = Transf_mat(3, 3);
                Transf_mat(3, 5) = data.k0s / 2 * pow(data.l * j / nslice, 2);
                Transf_mat(4, 5) = data.k0s * data.l * j / nslice * (1 - Ky / 6 * pow(data.l * j / nslice, 2));

                My = invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);

                a.Interpol[nslice-j-1].Disp_y = (Transf_mat(4, 4) * (Dispy(0) - Transf_mat(3, 5)) +
                                                 Transf_mat(3, 4) * (Transf_mat(4, 5) - Dispy(1))) /
                                                (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
                a.Interpol[nslice-j-1].Disp_y_prim = (Transf_mat(3, 3) * (Dispy(1) - Transf_mat(4, 5)) +
                                                      Transf_mat(4, 3) * (Transf_mat(3, 5) - Dispy(0))) /
                                                     (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
            }

        } else if(Ky > 0) {
            for(int j = 1; j <= nslice - 1; ++j) {
                Transf_mat(3, 3) = cos(data.l / nslice * j * sqrt(Ky));
                Transf_mat(3, 4) = sin(data.l / nslice * j * sqrt(Ky)) / sqrt(Ky);
                Transf_mat(4, 3) = -Ky * Transf_mat(3, 4);
                Transf_mat(4, 4) = Transf_mat(3, 3);
                Transf_mat(3, 5) = data.k0s * (1 - cos(data.l / nslice * j * sqrt(Ky))) / Ky;
                Transf_mat(4, 5) = data.k0s * sin(data.l / nslice * j * sqrt(Ky)) / sqrt(Ky);

                My = invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);

                a.Interpol[nslice-j-1].Disp_y = (Transf_mat(4, 4) * (Dispy(0) - Transf_mat(3, 5)) +
                                                 Transf_mat(3, 4) * (Transf_mat(4, 5) - Dispy(1))) /
                                                (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
                a.Interpol[nslice-j-1].Disp_y_prim = (Transf_mat(3, 3) * (Dispy(1) - Transf_mat(4, 5)) +
                                                      Transf_mat(4, 3) * (Transf_mat(3, 5) - Dispy(0))) /
                                                     (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
            }
        }

        else { //if(Ky<0)
            for(int j = 1; j <= nslice - 1; ++j) {
                Transf_mat(3, 3) = cosh(data.l / nslice * j * sqrt(-Ky));
                Transf_mat(3, 4) = sinh(data.l / nslice * j * sqrt(-Ky)) / sqrt(-Ky);
                Transf_mat(4, 3) = -Ky * Transf_mat(3, 4);
                Transf_mat(4, 4) = Transf_mat(3, 3);
                Transf_mat(3, 5) = data.k0s * (1 - cosh(data.l / nslice * j * sqrt(-Ky))) / -Ky;
                Transf_mat(4, 5) = data.k0s * sinh(data.l / nslice * j * sqrt(-Ky)) / sqrt(-Ky);

                My = invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);

                a.Interpol[nslice-j-1].Disp_y = (Transf_mat(4, 4) * (Dispy(0) - Transf_mat(3, 5)) +
                                                 Transf_mat(3, 4) * (Transf_mat(4, 5) - Dispy(1))) /
                                                (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
                a.Interpol[nslice-j-1].Disp_y_prim = (Transf_mat(3, 3) * (Dispy(1) - Transf_mat(4, 5)) +
                                                      Transf_mat(4, 3) * (Transf_mat(3, 5) - Dispy(0))) /
                                                     (Transf_mat(4, 4) * Transf_mat(3, 3) - Transf_mat(4, 3) * Transf_mat(3, 4));
            }
        }

    } else {
        for(int j = 0; j < nslice; ++j) {
            a.Interpol[j].Beta_x = tp->getBETi(*i, 0, 0);
            a.Interpol[j].Beta_y = tp->getBETi(*i, 1, 0);
            a.Interpol[j].Disp_x = tp->getDisp(*i, 0, 0);
            a.Interpol[j].Disp_y = tp->getDisp(*i, 2, 0);
            a.Interpol[j].Disp_x_prim = tp->getDisp(*i, 1, 0);
            a.Interpol[j].Disp_y_prim = tp->getDisp(*i, 3, 0);
        }
    }

    const std::string &nam1 = a.getElement()->getName();


    if(nam1 == "[DRIFT]") {
        a.Type_elm = "DRIFT";
        std::vector<double> dat = Attributes::getRealArray(itsAttr[DEFAULTAPERTURE]);

        tolerance_x = dat[0];
        tolerance_y = dat[1];
        offsetx = dat[2];
        offsety = dat[3];
        dat.erase(dat.begin(), dat.begin() + 4);
        ShapeBeamScreen = getShape(dat);
    } else {
        OpalElement &elem = dynamic_cast<OpalElement &>(*Element::find(nam1));
        a.Type_elm = elem.getBaseObject()->getOpalName();
        if(a.Type_elm == "MARKER") {
            std::vector<double> dat = Attributes::getRealArray(itsAttr[DEFAULTAPERTURE]);
            tolerance_x = dat[0];
            tolerance_y = dat[1];
            offsetx = dat[2];
            offsety = dat[3];
            dat.erase(dat.begin(), dat.begin() + 4);
            ShapeBeamScreen = getShape(dat);
        } else {
            auto vec = elem.getApert();

            tolerance_x = vec.second[0];
            tolerance_y = vec.second[1];
            offsetx = vec.second[2];
            offsety = vec.second[3];
            vec.second.erase(vec.second.begin(), vec.second.begin() + 4);
            ShapeBeamScreen = getShape(vec.second);

        }
    }

    if(data.l != 0) {
        for(int j = 0; j < nslice; ++j) {
            a.Interpol[j].delta_x = K_beta * delta_p * a.Interpol[j].Disp_x + tolerance_x + offsetx;
            a.Interpol[j].delta_y = K_beta * delta_p * a.Interpol[j].Disp_y + tolerance_y + offsety;
            calcul_Apert(a, j, tp);
        }
    } else {
        a.Interpol[0].delta_x = K_beta * delta_p * a.Interpol[0].Disp_x + tolerance_x + offsetx;
        a.Interpol[0].delta_y = K_beta * delta_p * a.Interpol[0].Disp_y + tolerance_y + offsety;
        calcul_Apert(a, 0, tp);
        for(int j = 1; j < nslice; ++j) {
            a.Interpol[j].delta_x = a.Interpol[0].delta_x;
            a.Interpol[j].delta_y = a.Interpol[0].delta_y;
            a.Interpol[j].apert = a.Interpol[0].apert;
        }
    }
    a.Orb = tp->getS(*i, 0, 0);
}

void Aperture::calcul_Apert(A_row &a, int slice, Twiss *tp) {
    int pt_halo = ShapeHalo.size() - 1;

    //non normalised edge of the halo induced by the the two collimators
    double *Pn = new double[pt_halo], *Qn = new double[pt_halo];
    double x = 0.0;

    int pt_delta(5);
    double *Apert = new double[pt_halo*pt_delta];

    double deltax, deltay;

    for(int k = 0; k < pt_delta; ++k) {
        deltax = a.Interpol[slice].delta_x + co * cos(k * Physics::pi / 8);
        deltay = a.Interpol[slice].delta_y + co * sin(k * Physics::pi / 8);

        const PartData &p = beam->getReference();

        for(int i = 0; i < pt_halo; ++i) {
            Pn[i] = sqrt(p.getM() * tp->getEX() / p.getP() * a.Interpol[slice].Beta_x) * ShapeHalo[i].x;
            Qn[i] = sqrt(p.getM() * tp->getEY() / p.getP() * a.Interpol[slice].Beta_y) * ShapeHalo[i].y;
            double a_d = Qn[i] / Pn[i];
            double b_d = deltay - a_d * deltax;

            int size = ShapeBeamScreen.size();
            for(int j = 0; j < size - 1; ++j) {

                double a_arc = (ShapeBeamScreen[j+1].y - ShapeBeamScreen[j].y) / (ShapeBeamScreen[j+1].x - ShapeBeamScreen[j].x);
                double b_arc = (ShapeBeamScreen[j].y * ShapeBeamScreen[j+1].x - ShapeBeamScreen[j+1].y * ShapeBeamScreen[j].x) /
                               (ShapeBeamScreen[j+1].x - ShapeBeamScreen[j].x);

                x = (b_d - b_arc) / (a_arc - a_d);
                double y = a_d * x + b_d;

                double ang = atan2(y, x);
                double ang1 = atan2(ShapeBeamScreen[j].y, ShapeBeamScreen[j].x);
                double ang2 = atan2(ShapeBeamScreen[j+1].y, ShapeBeamScreen[j+1].x);

                if((ang > ang1) && (ang <= ang2)) break;
            }
            Apert[i+k*pt_halo] = n1 * (x - deltax) / Pn[i];
        }
    }
    a.Interpol[slice].apert = Apert[0];
    for(int i = 1; i < pt_halo * pt_delta; ++i) {
        if(a.Interpol[slice].apert > Apert[i]) a.Interpol[slice].apert = Apert[i];
    }

    delete [] Apert;
    delete [] Qn;
    delete [] Pn;
}

vector<Aperture::coord> Aperture::getShape(vector<double> vec) {
    vector<coord> S;
    coord pt;
    double nb_pt(11);

    double r = vec[0];
    double h = vec[1];
    double v = vec[2];

    if((r > -1.5) && (r < -0.5)) {
        double ang = atan2(v, h);
        for(int i = 0; i < nb_pt; ++i) {
            if(i *Physics::pi / 20 <= ang) {
                pt.x = h;
                pt.y = h * tan(i * Physics::pi / 20);
                S.push_back(pt);
            } else {
                pt.x = v / tan(i * Physics::pi / 20);
                pt.y = v;
                S.push_back(pt);
            }
        }
    } else if((r > -10.5) && (r < -9.5)) {
        vector<double>::iterator iter = vec.begin();
        ++iter;
        while(iter != vec.end()) {
            pt.x = *iter;
            pt.y = *(iter + 1);
            S.push_back(pt);
            iter = iter + 2;
        }
    } else {
        if((r != 0) && (h == 0) && (v == 0)) {
            for(int i = 0; i < nb_pt; ++i) {
                pt.x = r * cos(i * Physics::pi / 20);
                pt.y = r * sin(i * Physics::pi / 20);
                S.push_back(pt);
            }
        } else if((r != 0) && (h != 0) && (v == 0)) {
            double ang = acos(h / r);
            for(int i = 0; i < nb_pt; ++i) {
                if(i *Physics::pi / 20 <= ang) {
                    pt.x = h;
                    pt.y = h * tan(i * Physics::pi / 20);
                    S.push_back(pt);
                } else {
                    pt.x = r * cos(i * Physics::pi / 20);
                    pt.y = r * sin(i * Physics::pi / 20);
                    S.push_back(pt);
                }
            }
        } else if((r != 0) && (h == 0) && (v != 0)) {
            double ang = asin(v / r);
            for(int i = 0; i < nb_pt; ++i) {
                if(i *Physics::pi / 20 <= ang) {
                    pt.x = r * cos(i * Physics::pi / 20);
                    pt.y = r * sin(i * Physics::pi / 20);
                    S.push_back(pt);
                } else {
                    pt.x = r / tan(i * Physics::pi / 20);
                    pt.y = v;
                    S.push_back(pt);
                }
            }
        } else if((r != 0) && (h != 0) && (v != 0)) {
            double ang1 = acos(h / r);
            double ang2 = asin(v / r);
            for(int i = 0; i < nb_pt; ++i) {
                if(i *Physics::pi / 20 <= ang1) {
                    pt.x = h;
                    pt.y = h * tan(i * Physics::pi / 20);
                    S.push_back(pt);
                } else if((i * Physics::pi / 20 > ang1) && (i *Physics::pi / 20 <= ang2)) {
                    pt.x = r * cos(i * Physics::pi / 20);
                    pt.y = r * sin(i * Physics::pi / 20);
                    S.push_back(pt);
                } else { //i*Physics::pi/20>ang2
                    pt.x = r / tan(i * Physics::pi / 20);
                    pt.y = v;
                    S.push_back(pt);
                }
            }
        } else { //(r=0,h!=0,v!=0) case(000)(0x0)(00x) uninteresting
            for(int i = 0; i < nb_pt; ++i) {
                pt.x = h * cos(i * Physics::pi / 20);
                pt.y = v * sin(i * Physics::pi / 20);
                S.push_back(pt);
            }
        }
    }

    return S;
}

void Aperture::execute() {

    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    beam = Beam::find(beamName);

    std::vector<double> dat = Attributes::getRealArray(itsAttr[DATA]);
    co = dat[0];
    delta_p = dat[1];
    K_beta = dat[2];
    n1 = dat[3];
    dat.erase(dat.begin(), dat.begin() + 4);
    ShapeHalo = getShape(dat);

    //fill the table
    fill();
    if(Attributes::getBool(itsAttr[STATIC])) dynamic = false;

    //output of the Aperture command
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    std::string file_out = Attributes::getString(itsAttr[FILE]);
    ofstream outFile(file_out.c_str());
    if(!outFile) {
        cerr << "Aperture: unable to open output file:" << file_out << endl;
        exit(-1);
    }

    outFile.setf(std::ios::scientific, std::ios::floatfield);

    OPALTimer::Timer timer;
    outFile << "@ NAME    %s " << getOpalName() << endl;
    outFile << "@ DATE    %s " << timer.date() << endl;
    outFile << "@ TIME    %s " << timer.time() << endl;
    outFile << "@ ORIGIN  %s  OPAL_9.5/7\n";
    outFile << "@ COMMENT %s " << endl;
    outFile << "@ NSLICE %e " << nslice << endl;

    outFile << "*" << setw(15) << "NAME" << setw(15) << "Type" << setw(15) << "S" << setw(15) << "L" << setw(15) << "Apert" << endl;
    outFile << "$" << setw(15) << "%23s" << setw(15) << "%12s" << setw(15) << "%e" << setw(15) << "%e" << setw(15) << "%e" << endl;

}

void Aperture::fill() {
    if(refill) {
        run();
        refill = false;
    }
}

void Aperture::run() {
    itsTable.erase(itsTable.begin(), itsTable.end());
    double nslice = Attributes::getReal(itsAttr[NSLICE]);

    itsLine = Attributes::getString(itsAttr[TABLE]);
    Table *t = Table::find(itsLine);
    tp = dynamic_cast<Twiss *>(t);

    if(tp == 0) {
        throw OpalException("Aperture::execute()", "Unknown Twiss Table\"" + itsLine + "\".");
    }
    tp->fill();

    vector<double> dat = Attributes::getRealArray(itsAttr[DEFAULTAPERTURE]);
    if(dat.size() != 0) {
        i = tp->begin();
        while(i != tp->end()) {
            const std::string nam = i->getElement()->getName();
            if(nam == "[DRIFT]") {
                A_row row(*i, static_cast<int>(nslice));
                i->accept(*this);
                calcul(i, row, static_cast<int>(nslice), tp);
                itsTable.append(row);
                i++;
            } else {
                OpalElement &elem = dynamic_cast<OpalElement &>(*Element::find(nam));
                std::string Typ = elem.getBaseObject()->getOpalName();
                if((elem.getApert().second.size() == 0) && (Typ != "MARKER")) i++;
                else {
                    A_row row(*i, static_cast<int>(nslice));
                    i->accept(*this);
                    calcul(i, row, static_cast<int>(nslice), tp);
                    itsTable.append(row);
                    i++;
                }
            }
        }

    } else {
        i = tp->begin();
        while(i != tp->end()) {
            const std::string nam = i->getElement()->getName();
            if(nam == "[DRIFT]") i++;
            else {
                OpalElement &elem = dynamic_cast<OpalElement &>(*Element::find(nam));
                std::string Typ = elem.getBaseObject()->getOpalName();
                if(Typ == "MARKER") i++;
                else if(elem.getApert().second.size() == 0) i++;
                else {
                    A_row row(*i, static_cast<int>(nslice));
                    i->accept(*this);
                    calcul(i, row, static_cast<int>(nslice), tp);
                    itsTable.append(row);
                    i++;
                }
            }
        }
    }

}

const Aperture::A_row &Aperture::getCurrent() const {
    return *current;
}

double Aperture::getLength() {
    return itsTable.getElementLength();
}

double Aperture::getCell(const PlaceRep &place, const std::string &colName) {
    A_row &row = findRow(place);
    const ColDesc *col = findCol(*this, colName);
    return (this->*(col->get))(row, col->ind_1, col->ind_2);
}
const Beamline *Aperture::getLine() const {
    return &itsTable;
}
Table::CellArray Aperture::getDefault() const {
    CellArray columns;
    for(const ColDesc *col = defaultColumns; col->colName; ++col) {
        Expressions::PtrToScalar<double> expr =
            new Column(*this, col->colName, *col);
        columns.push_back(Cell(expr, col->printWidth, col->printPrecision));
    }
    return columns;
}
std::vector<double> Aperture::getColumn(const RangeRep &rng, const std::string
                                        &colName) {
    const ColDesc *col = findCol(*this, colName);
    RangeRep range(rng);
    range.initialize();
    std::vector<double> column;

    for(A_Tline::const_iterator row = begin(); row != end(); ++row) {
        range.enter(*row);
        if(range.isActive()) {
            column.push_back((this->*(col->get))(*row, col->ind_1, col->ind_2));
        }
        range.leave(*row);
    }

    return column;
}
std::vector<double> Aperture::getRow(const PlaceRep &pos, const
                                     std::vector<std::string> &cols) {
    A_row &row = findRow(pos);
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
bool Aperture::isDependent(const std::string &name) const {
    // Test if name refers to USE attribute.
    if(itsLine == name) return true;

    // Test if name occurs in table.
    for(A_Tline::const_iterator row = begin(); row != end(); ++row) {
        if(row->getElement()->getName() == name) return true;
    }

    // Otherwise replacement is not required.
    return false;
}
bool Aperture::matches(Table */*rhs*/) const { return false; }
Expressions::PtrToScalar<double> Aperture::makeColumnExpression(const std::string &colname) const {
    const ColDesc *col = findCol(*this, colname);
    return new Column(*this, colname, *col);
}
Object *Aperture::clone(const std::string &name) {
    return new Aperture(name, this);
}

void Aperture::printTable(std::ostream &, const CellArray &)const {};

Aperture::A_row &Aperture::findRow(const PlaceRep &place) {
    PlaceRep row(place);
    row.initialize();

    for(A_Tline::iterator i = begin(); i != end(); ++i) {
        row.enter(*i);
        if(row.isActive()) return *i;
        row.leave(*i);
    }

    std::ostringstream os;
    os << row << std::ends;
    throw OpalException("Aperture::findRow()", "A_row \"" + os.str() +
                        "\" not found in twiss table \"" + getOpalName() + "\".");
}