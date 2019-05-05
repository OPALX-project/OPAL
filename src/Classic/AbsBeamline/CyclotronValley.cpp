// ------------------------------------------------------------------------
// $RCSfile: CyclotronValley.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronValley
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/CyclotronValley.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class CyclotronValley
// ------------------------------------------------------------------------

CyclotronValley::CyclotronValley():
    Component(),
    filename_m(""),
    fieldmap_m(NULL),
    scale_m(1.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    fast_m(false) {
    setElType(isRF);
}


CyclotronValley::CyclotronValley(const CyclotronValley &right):
    Component(right),
    filename_m(right.filename_m),
    fieldmap_m(right.fieldmap_m),
    scale_m(right.scale_m),
    ElementEdge_m(right.ElementEdge_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    fast_m(right.fast_m) {
    setElType(isRF);
}


CyclotronValley::CyclotronValley(const std::string &name):
    Component(name),
    filename_m(""),
    fieldmap_m(NULL),
    scale_m(1.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    fast_m(false) {
    setElType(isRF);
}


CyclotronValley::~CyclotronValley() {
  //    Fieldmap::deleteFieldmap(filename_m);
    /*if(RNormal_m) {
        delete[] RNormal_m;
        delete[] VrNormal_m;
        delete[] DvDr_m;
        }*/
}


void CyclotronValley::accept(BeamlineVisitor &visitor) const {
    visitor.visitCyclotronValley(*this);
}

void CyclotronValley::setFieldMapFN(std::string fn) {
    filename_m = fn;
}

std::string CyclotronValley::getFieldMapFN() const {
    return filename_m;
}

void CyclotronValley::setFast(bool fast) {
    fast_m = fast;
}


bool CyclotronValley::getFast() const {
    return fast_m;
}

bool CyclotronValley::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t tmpR = RefPartBunch_m->R[i];
    return applyToReferenceParticle(tmpR, RefPartBunch_m->P[i], t, E, B);
}

bool CyclotronValley::apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t tmpR = R;
    return applyToReferenceParticle(tmpR, P, t, E, B);
}

bool CyclotronValley::applyToReferenceParticle(const Vector_t &tmpR, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {
    Vector_t  tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);

    if(!fieldmap_m->getFieldstrength(tmpR, tmpE, tmpB)) {
        B +=scale_m * tmpB;
        return false;
    }
    return true;
}

void CyclotronValley::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {// called by ParallelTTracker::visitCyclotronValley --> OpalBeamline::visit
    Inform msg("CyclotronValley ");
    RefPartBunch_m = bunch;

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    if(fieldmap_m != NULL) {
        double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
        if(zEnd > zBegin) {
            msg << getName() << " using file: " << "  ";
            fieldmap_m->getInfo(&msg);

            ElementEdge_m = startField;
            startField_m = startField = ElementEdge_m + zBegin;
            endField_m = endField = ElementEdge_m + zEnd;

        } else {
            endField = startField - 1e-3;
        }
    } else {
        endField = startField - 1e-3;
    }
}


void CyclotronValley::finalise()
{}

bool CyclotronValley::bends() const {
    return false;
}


void CyclotronValley::goOnline(const double &) {
    Fieldmap::readMap(filename_m);
    online_m = true;
}

void CyclotronValley::goOffline() {
    Fieldmap::freeMap(filename_m);
    online_m = false;
}





void CyclotronValley::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementBase::ElementType CyclotronValley::getType() const {
    return CYCLOTRONVALLEY;
}