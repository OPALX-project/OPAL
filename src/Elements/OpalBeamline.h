//
// Class OpalBeamline
//   :FIXME: add class description
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
#ifndef OPAL_BEAMLINE_H
#define OPAL_BEAMLINE_H

#include <set>
#include <string>

#include "Beamlines/Beamline.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Source.h"

#include "Utilities/ClassicField.h"

#include "Algorithms/CoordinateSystemTrafo.h"

template <class T, unsigned Dim>
class PartBunchBase;
class ParticleMatterInteractionHandler;
class BoundaryGeometry;

class OpalBeamline {

public:
    OpalBeamline();
    OpalBeamline(const Vector_t& origin,
                 const Quaternion& rotation);
    ~OpalBeamline();

    void activateElements();
    std::set<std::shared_ptr<Component>> getElements(const Vector_t &x);
    Vector_t transformTo(const Vector_t &r) const;
    Vector_t transformFrom(const Vector_t &r) const;
    Vector_t rotateTo(const Vector_t &r) const;
    Vector_t rotateFrom(const Vector_t &r) const;

    Vector_t transformToLocalCS(const std::shared_ptr<Component> &comp,
                                const Vector_t &r) const;
    Vector_t transformFromLocalCS(const std::shared_ptr<Component> &comp,
                                  const Vector_t &r) const;
    Vector_t rotateToLocalCS(const std::shared_ptr<Component> &comp,
                             const Vector_t &r) const;
    Vector_t rotateFromLocalCS(const std::shared_ptr<Component> &comp,
                               const Vector_t &r) const;
    CoordinateSystemTrafo getCSTrafoLab2Local(const std::shared_ptr<Component> &comp) const;
    CoordinateSystemTrafo getCSTrafoLab2Local() const;
    CoordinateSystemTrafo getMisalignment(const std::shared_ptr<Component> &comp) const;

    double getStart(const Vector_t &) const;
    double getEnd(const Vector_t &) const;

    void switchElements(const double &, const double &, const double &kineticEnergy, const bool &nomonitors = false);
    void switchElementsOff();

    ParticleMatterInteractionHandler *getParticleMatterInteractionHandler(const unsigned int &);

    BoundaryGeometry *getBoundaryGeometry(const unsigned int &);

    unsigned long getFieldAt(const unsigned int &, const Vector_t &, const long &, const double &, Vector_t &, Vector_t &);
    unsigned long getFieldAt(const Vector_t &, const Vector_t &, const double &, Vector_t &, Vector_t &);

    template<class T>
    void visit(const T &, BeamlineVisitor &, PartBunchBase<double, 3> *);

    void prepareSections();
    void positionElementRelative(std::shared_ptr<ElementBase>);
    void compute3DLattice();
    void save3DLattice();
    void save3DInput();
    void print(Inform &) const;

    FieldList getElementByType(ElementType);

    void swap(OpalBeamline & rhs);
    void merge(OpalBeamline &rhs);

    bool containsSource();
private:
    FieldList elements_m;
    bool prepared_m;
    bool containsSource_m;

    CoordinateSystemTrafo coordTransformationTo_m;
};

template<class T> inline
void OpalBeamline::visit(const T &element, BeamlineVisitor &, PartBunchBase<double, 3> *bunch) {
    Inform msg("OPAL ");
    double startField = 0.0;
    double endField = 0.0;
    std::shared_ptr<T> elptr(dynamic_cast<T *>(element.clone()));

    positionElementRelative(elptr);

    if (elptr->isElementPositionSet())
        startField = elptr->getElementPosition();

    elptr->initialise(bunch, startField, endField);
    elements_m.push_back(ClassicField(elptr, startField, endField));
}

template<> inline
void OpalBeamline::visit<Source>(const Source &element, BeamlineVisitor &, PartBunchBase<double, 3> *bunch) {
    containsSource_m = true;
    double startField = 0.0;
    double endField = 0.0;
    std::shared_ptr<Source> elptr(dynamic_cast<Source *>(element.clone()));

    positionElementRelative(elptr);

    if (elptr->isElementPositionSet())
        startField = elptr->getElementPosition();

    elptr->initialise(bunch, startField, endField);
    elements_m.push_back(ClassicField(elptr, startField, endField));
}

template<> inline
void OpalBeamline::visit<Marker>(const Marker &/*element*/, BeamlineVisitor &, PartBunchBase<double, 3> *) {
}

template<> inline
void OpalBeamline::visit<Septum>(const Septum &element, BeamlineVisitor &, PartBunchBase<double, 3> *) {
    WARNMSG(element.getTypeString() << " not implemented yet!" << endl);
}

inline
Vector_t OpalBeamline::transformTo(const Vector_t &r) const {
    return coordTransformationTo_m.transformTo(r);
}

inline
Vector_t OpalBeamline::transformFrom(const Vector_t &r) const {
    return coordTransformationTo_m.transformFrom(r);
}

inline
Vector_t OpalBeamline::rotateTo(const Vector_t &r) const {
    return coordTransformationTo_m.rotateTo(r);
}

inline
Vector_t OpalBeamline::rotateFrom(const Vector_t &r) const {
    return coordTransformationTo_m.rotateFrom(r);
}

inline
Vector_t OpalBeamline::transformToLocalCS(const std::shared_ptr<Component> &comp,
                                          const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().transformTo(r);
}

inline
Vector_t OpalBeamline::transformFromLocalCS(const std::shared_ptr<Component> &comp,
                                            const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().transformFrom(r);
}

inline
Vector_t OpalBeamline::rotateToLocalCS(const std::shared_ptr<Component> &comp,
                                       const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().rotateTo(r);
}

inline
Vector_t OpalBeamline::rotateFromLocalCS(const std::shared_ptr<Component> &comp,
                                         const Vector_t &r) const {
    return comp->getCSTrafoGlobal2Local().rotateFrom(r);
}

inline
CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local(const std::shared_ptr<Component> &comp) const {
    return comp->getCSTrafoGlobal2Local();
}

inline
CoordinateSystemTrafo OpalBeamline::getCSTrafoLab2Local() const {
    return coordTransformationTo_m;
}

inline
CoordinateSystemTrafo OpalBeamline::getMisalignment(const std::shared_ptr<Component> &comp) const {
    return comp->getMisalignment();
}

inline
bool OpalBeamline::containsSource() {
    return containsSource_m;
}
#endif // OPAL_BEAMLINE_H