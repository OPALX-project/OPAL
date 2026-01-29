//
// Helpers for unit testing field elements; does numerical derivatives to check
// Maxwell's equations
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//


#include "AbsBeamline/Component.h"
#include "opal_test_utilities/Maxwell.h"

MaxwellTest::MaxwellTest(Vector_t dR, double /*dt*/, Component* field) :
    field_m(field),
    dR_m(dR) /*, dt_m(dt)*/ {
}


Vector_t MaxwellTest::getB(const Vector_t &R, double t) const {
    Vector_t P, E, B;
    field_m->apply(R, P, t, E, B);
    return B;
}

std::vector< std::vector<double> > MaxwellTest::partialsDerivB(const Vector_t &R, double t) const {
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    std::vector< std::vector<double> > allPartials(3, std::vector<double>(3));
    Vector_t P, E;
    for(int i = 0; i < 3; i++) {
        // B at the previous and next grid points R_prev,  R_next
        Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
        R_pprev(i) -= 2 * dR_m(i);
        R_nnext(i) += 2 * dR_m(i);
        R_prev(i) -= dR_m(i);
        R_next(i) += dR_m(i);
        Vector_t B_prev, B_next, B_pprev, B_nnext;
        field_m->apply(R_prev, P, t, E, B_prev);
        field_m->apply(R_next, P, t, E, B_next);
        field_m->apply(R_pprev, P, t, E, B_pprev);
        field_m->apply(R_nnext, P, t, E, B_nnext);
        for(int j = 0; j < 3; j++) {
              allPartials[i][j] =
                   (B_pprev[j]-8*B_prev[j]+8*B_next[j]-B_nnext[j])/(12*dR_m(i));
        }
    }
    return allPartials;
}

std::vector< std::vector<double> > MaxwellTest::partialsDerivA(const Vector_t &R, double t) const {
    // builds a matrix of all partial derivatives of A -> dx_i A_j
    std::vector< std::vector<double> > allPartials(3, std::vector<double>(3));
    double phi;
    for(int i = 0; i < 3; i++) {
        // A at the previous and next grid points R_prev,  R_next
        Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
        R_pprev(i) -= 2 * dR_m(i);
        R_nnext(i) += 2 * dR_m(i);
        R_prev(i) -= dR_m(i);
        R_next(i) += dR_m(i);
        Vector_t A_prev, A_next, A_pprev, A_nnext;
        field_m->getPotential(R_prev, t, A_prev, phi);
        field_m->getPotential(R_next, t, A_next, phi); // next
        field_m->getPotential(R_pprev, t, A_pprev, phi); //pprev
        field_m->getPotential(R_nnext, t, A_nnext, phi); // nnext
        for(int j = 0; j < 3; j++) {
              allPartials[i][j] = 
                   (A_pprev[j]-8*A_prev[j]+8*A_next[j]-A_nnext[j])/(12*dR_m(i));
        }
    }
    return allPartials;
}


double MaxwellTest::divB(const Vector_t &R, double t) const {
    double div = 0;
    std::vector< std::vector<double> > partials = partialsDerivB(R, t);
    for(int i = 0; i < 3; i++)
        div += partials[i][i];
    return div;
}

Vector_t MaxwellTest::curlB(const Vector_t &R, double t) const {
    Vector_t curl;
    std::vector< std::vector<double> > partials = partialsDerivB(R, t);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}

Vector_t MaxwellTest::curlA(const Vector_t &R, double t) const {
    Vector_t curl;
    std::vector< std::vector<double> > partials = partialsDerivA(R, t);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}

std::ostream& MaxwellTest::printHeading(std::ostream& out) const {
    return out;
}

std::ostream& MaxwellTest::printLine(std::ostream& out, const Vector_t& R, double t) const {
    bool printPartials = false;
    Vector_t P, bfield, efield;
    field_m->apply(R, P, t, efield, bfield);
    std::vector< std::vector<double> > partials = partialsDerivB(R, t);
    out << R << " " << t << " ** " << bfield << " " << efield << " ** ";
    if (printPartials) {
        out << partials[0][0] << " " << partials[1][1] << " " << partials[2][2] << " ";
    }
    out << divB(R, t) << " ** " << curlB(R, t) << std::endl;
    return out;
}

