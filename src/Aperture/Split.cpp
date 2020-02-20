#include "Aperture/Split.h"
#include "Attributes/Attributes.h"
#include "Structure/Beam.h"
#include "Fields/BMultipoleField.h"
#include "AbstractObjects/Expressions.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/Multipole.h"
#include "Physics/Physics.h"
#include "AbstractObjects/RangeRep.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Cyclotron.h"
#include "Utilities/Timer.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;

// Local structures.
// ------------------------------------------------------------------------

namespace {

    //return the inverse of the "alpha,beta,gamma" transformation matrix

    FMatrix<double, 3, 3> sp_invers(double c, double c_prim, double s, double s_prim) {
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
        double(MSplit::*get)(const MSplit::A_row &, int, int) const;

        // Default field width and precision.
        int printWidth, printPrecision;

        // Indices to be given to get().
        int ind_1, ind_2;
    };


    // The complete column entries table.
    const ColDesc allColumns[] = {

        // "Naive" maximum Twiss.
        { "BETXMAX",  &MSplit::getBETXMAX,  10, 6, 0, 0 },
        { "BETYMAX",  &MSplit::getBETYMAX,  10, 6, 0, 0 },
        { 0,       0,                 0, 0, 0, 0 }
    };


    // The default column entries table.
    const ColDesc defaultColumns[] = {
        // "Naive" maximum Twiss.
        { "BETXMAX",  &MSplit::getBETXMAX,  10, 6, 0, 0 },
        { "BETYMAX",  &MSplit::getBETYMAX,  10, 6, 0, 0 },
        { 0,       0,                0, 0, 0, 0 }
    };

    const ColDesc *findCol(const MSplit &table, const std::string &colName) {
        for(const ColDesc *col = allColumns; col->colName; ++col) {
            if(colName == col->colName) {
                return col;
            }
        }

        throw OpalException("Split::findCol()",
                            "Split table \"" + table.getOpalName() +
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
        Column(const MSplit &tab, const std::string &colName, const ColDesc &desc);

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
        const MSplit &itsTable;

        // Column name.
        std::string colName;

        // The function returning the column value.
        double(MSplit::*get)(const MSplit::A_row &, int, int) const;

        // The indices to be transmitted to get().
        int ind_1, ind_2;
    };


    // Implementation.
    // ------------------------------------------------------------------------

    Column::Column(const MSplit &tab, const std::string &colName, const ColDesc &desc):
        itsTable(tab), colName(colName),
        get(desc.get), ind_1(desc.ind_1), ind_2(desc.ind_2)
    {}


    Column::Column(const Column &rhs):
        Expressions::Scalar<double>(rhs),
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

MSplit::MSplit():
    DefaultVisitor(itsTable, false, false),
    Table(SIZE, "SPLIT", "help"), itsTable() {
    itsAttr[LINE] = Attributes::makeString
                    ("LINE", "The beam line use to filling");
    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "The beam to be used", "UNNAMED_BEAM");
    itsAttr[NSLICE] = Attributes::makeReal
                      ("NSLICE", "the number of slices of the interpolation inside the optic element", 10);
    itsAttr[STATIC] = Attributes::makeBool
                      ("STATIC", "recalculation if static equal false", true);
    itsAttr[FILE] = Attributes::makeString
                    ("FILE", "Name of file to receive SPLIT output", "SPLIT.dat");

    registerOwnership(AttributeHandler::STATEMENT);
}
MSplit::MSplit(const std::string &name, MSplit *parent):
    DefaultVisitor(itsTable, false, false),
    Table(name, parent), itsTable(name)
{}

MSplit::A_row::A_row(const A_row &a): FlaggedElmPtr(a), Interpol(a.Interpol) {}

inline double MSplit::A_row::getBeta_x(int ind) {
    return Interpol[ind].Beta_x;
}
inline double MSplit::A_row::getBeta_y(int ind) {
    return Interpol[ind].Beta_y;
}
inline double MSplit::A_row::getAlpha_x(int ind) {
    return Interpol[ind].Alpha_x;
}
inline double MSplit::A_row::getAlpha_y(int ind) {
    return Interpol[ind].Alpha_y;
}
inline double MSplit::A_row::getDisp_x(int ind) {
    return Interpol[ind].Disp_x;
}
inline double MSplit::A_row::getDisp_x_prim(int ind) {
    return Interpol[ind].Disp_x_prim;
}
inline double MSplit::A_row::getDisp_y(int ind) {
    return Interpol[ind].Disp_y;
}
inline double MSplit::A_row::getDisp_y_prim(int ind) {
    return Interpol[ind].Disp_y_prim;
}

//iterator used for the displacement in the list of optic elments
MSplit::A_Tline::iterator MSplit::begin() {
    return itsTable.begin();
}
MSplit::A_Tline::const_iterator MSplit::begin() const {
    return itsTable.begin();
}
MSplit::A_Tline::iterator MSplit::end() {
    return itsTable.end();
}
MSplit::A_Tline::const_iterator MSplit::end() const {
    return itsTable.end();
}

double MSplit::getBETXMAX(const A_row &row, int /*i1*/, int /*i2*/) const {
    double max;
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    max = row.Interpol[0].Beta_x;
    for(int i = 1; i < nslice; ++i) {
        if(max < row.Interpol[i].Beta_x) max = row.Interpol[i].Beta_x;
    }
    return max;
}

double MSplit::getBETYMAX(const A_row &row, int /*i1*/, int /*i2*/) const {
    double max;
    double nslice = Attributes::getReal(itsAttr[NSLICE]);
    max = row.Interpol[0].Beta_y;
    for(int i = 1; i < nslice; ++i) {
        if(max < row.Interpol[i].Beta_y) max = row.Interpol[i].Beta_y;

    }
    return max;
}

void MSplit::visitMultipole(const Multipole &m) {
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

void MSplit::visitRBend(const RBend &rb) {
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

void MSplit::visitCyclotron(const Cyclotron &/*cy*/) {
    ERRORMSG("MSplit::visitCyclotron(const Cyclotron &cy) not implemented");
}

void MSplit::visitSBend(const SBend &sb) {
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

void MSplit::applyDefault(const ElementBase &eb) {
    data.type = "drift";
    data.l = eb.getElementLength();
    data.courb = 0;
    data.phi = 0;
    data.kq = 0;
    data.k0s = 0;
}


void MSplit::calcul(Twiss::TLine::iterator i, A_row &a, int nslice, Twiss *tp) {
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
        a.Interpol[nslice-1].Alpha_x = tp->getALFi(*i, 0, 0);
        a.Interpol[nslice-1].Alpha_y = tp->getALFi(*i, 1, 0);
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

            Mx = sp_invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

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

            My = sp_invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));

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
            //propagetion of x twiss function through the exit face
            Transf_mat(0, 0) = 1;
            Transf_mat(0, 1) = 0;
            Transf_mat(1, 0) = data.courb * tan(data.e2);
            Transf_mat(1, 1) = 1;
            Transf_mat(0, 2) = 0;
            Transf_mat(1, 2) = 0;

            Mx = sp_invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

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

            //propagetion of y twiss function through the exit face

            Transf_mat(3, 3) = 1;
            Transf_mat(3, 4) = 0;
            Transf_mat(4, 3) = -data.courb * tan(data.e2);
            Transf_mat(4, 4) = 1;
            Transf_mat(3, 5) = 0;
            Transf_mat(4, 5) = 0;

            My = sp_invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));

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

                Mx = sp_invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));

                //Filling of the vectors containing Beta et Disp inside the element

                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);
                a.Interpol[nslice-j-1].Alpha_x = Mx(1, 0) * Euler_x(0) + Mx(1, 1) * Euler_x(1) + Mx(1, 2) * Euler_x(2);

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

                Mx = sp_invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));
                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);
                a.Interpol[nslice-j-1].Alpha_x = Mx(1, 0) * Euler_x(0) + Mx(1, 1) * Euler_x(1) + Mx(1, 2) * Euler_x(2);

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

                Mx = sp_invers(Transf_mat(0, 0), Transf_mat(1, 0), Transf_mat(0, 1), Transf_mat(1, 1));
                a.Interpol[nslice-j-1].Beta_x = Mx(0, 0) * Euler_x(0) + Mx(0, 1) * Euler_x(1) + Mx(0, 2) * Euler_x(2);
                a.Interpol[nslice-j-1].Alpha_x = Mx(1, 0) * Euler_x(0) + Mx(1, 1) * Euler_x(1) + Mx(1, 2) * Euler_x(2);

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

                My = sp_invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);
                a.Interpol[nslice-j-1].Alpha_y = My(1, 0) * Euler_y(0) + My(1, 1) * Euler_y(1) + My(1, 2) * Euler_y(2);

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

                My = sp_invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);
                a.Interpol[nslice-j-1].Alpha_y = My(1, 0) * Euler_y(0) + My(1, 1) * Euler_y(1) + My(1, 2) * Euler_y(2);

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

                My = sp_invers(Transf_mat(3, 3), Transf_mat(4, 3), Transf_mat(3, 4), Transf_mat(4, 4));
                a.Interpol[nslice-j-1].Beta_y = My(0, 0) * Euler_y(0) + My(0, 1) * Euler_y(1) + My(0, 2) * Euler_y(2);
                a.Interpol[nslice-j-1].Alpha_y = My(1, 0) * Euler_y(0) + My(1, 1) * Euler_y(1) + My(1, 2) * Euler_y(2);

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
            a.Interpol[j].Alpha_x = tp->getALFi(*i, 0, 0);
            a.Interpol[j].Alpha_y = tp->getALFi(*i, 1, 0);
            a.Interpol[j].Disp_x_prim = tp->getDisp(*i, 1, 0);
            a.Interpol[j].Disp_y_prim = tp->getDisp(*i, 3, 0);
        }
    }
}

void MSplit::execute() {
    const std::string &beamName = Attributes::getString(itsAttr[BEAM]);
    beam = Beam::find(beamName);
    fill();
    if(Attributes::getBool(itsAttr[STATIC])) dynamic = false;
}

void MSplit::fill() {
    if(refill) run();
    refill = false;
}

void MSplit::run() {
    itsTable.erase(itsTable.begin(), itsTable.end());

    int nslice = static_cast<int>(Attributes::getReal(itsAttr[NSLICE]));
    itsLine = Attributes::getString(itsAttr[LINE]);
    Table *t = Table::find(itsLine);
    tp = dynamic_cast<Twiss *>(t);

    if(tp == 0) {
        throw OpalException("MSplit::execute()", "Unknown Twiss Table\"" + itsLine + "\".");
    }
    tp->fill();
    std::string file = Attributes::getString(itsAttr[FILE]);
    ofstream outFile(file.c_str());
    if(!outFile) {
        cerr << "MSplit: Cannot open output file."
             << endl;
        exit(-1);
    }

    outFile.setf(std::ios::scientific, std::ios::floatfield);

    OPALTimer::Timer timer;
    outFile << "@ NAME    %s " << getOpalName() << endl;
    outFile << "@ DATE    %s " << timer.date() << endl;
    outFile << "@ TIME    %s " << timer.time() << endl;
    outFile <<  "@ ORIGIN  %s  OPAL_9.5/7\n";
    outFile << "@ COMMENT %s " << endl;

    outFile << "*" << setw(15) << "NAME" << setw(15) << "S"
            << setw(15) << "L" << setw(15) << "angle" << setw(15) << "K1"
            << setw(15) << "BETX" << setw(15) << "BETY" << setw(15)
            << "ALFX" << setw(15) << "ALFY" << setw(15) << "DX"
            << setw(15) << "DPX"
            << setw(15) << "DY" << setw(15) << "DPY" << endl;

    outFile << "$" << setw(15) << "%23s" << setw(15) << "%e" << setw(15) << "%e"
            << setw(15) << "%e" << setw(15) << "%e" << setw(15) << "%e"
            << setw(15) << "%e" << setw(15) << "%e" << setw(15)
            << "%e" << setw(15) << "%e" << setw(15) << "%e" << setw(15) << "%e"
            << setw(15) << "%e" << endl;


    for(i = tp->begin(); i != tp->end(); ++i) {
        A_row row(*i, nslice);
        i->accept(*this);
        calcul(i, row, nslice, tp);
        itsTable.append(row);
    }

}

const MSplit::A_row &MSplit::getCurrent() const {
    return *current;
}

double MSplit::getLength() {
    return itsTable.getElementLength();
}

double MSplit::getCell(const PlaceRep &place, const std::string &colName) {
    A_row &row = findRow(place);
    const ColDesc *col = findCol(*this, colName);
    return (this->*(col->get))(row, col->ind_1, col->ind_2);
}
const Beamline *MSplit::getLine() const {
    return &itsTable;
}
Table::CellArray MSplit::getDefault() const {
    CellArray columns;
    for(const ColDesc *col = defaultColumns; col->colName; ++col) {
        Expressions::PtrToScalar<double> expr =
            new Column(*this, col->colName, *col);
        columns.push_back(Cell(expr, col->printWidth, col->printPrecision));
    }
    return columns;
}
std::vector<double> MSplit::getColumn(const RangeRep &rng, const std::string
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
std::vector<double> MSplit::getRow(const PlaceRep &pos, const
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
bool MSplit::isDependent(const std::string &name) const {
    // Test if name refers to USE attribute.
    if(itsLine == name) return true;

    // Test if name occurs in table.
    for(A_Tline::const_iterator row = begin(); row != end(); ++row) {
        if(row->getElement()->getName() == name) return true;
    }

    // Otherwise replacement is not required.
    return false;
}
bool MSplit::matches(Table */*rhs*/) const { return false; }
Expressions::PtrToScalar<double> MSplit::makeColumnExpression(const std::string &colname) const {
    const ColDesc *col = findCol(*this, colname);
    return new Column(*this, colname, *col);
}
Object *MSplit::clone(const std::string &name) {
    return new MSplit(name, this);
}

void MSplit::printTable(std::ostream &, const CellArray &)const {};

MSplit::A_row &MSplit::findRow(const PlaceRep &place) {
    PlaceRep row(place);
    row.initialize();

    for(A_Tline::iterator i = begin(); i != end(); ++i) {
        row.enter(*i);
        if(row.isActive()) return *i;
        row.leave(*i);
    }

    std::ostringstream os;
    os <<  row <<  std::ends;
    throw OpalException("MSplit::findRow()", "A_row \"" + os.str() +
                        "\" not found in twiss table \"" + getOpalName() + "\".");
}