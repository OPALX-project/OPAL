//
// Class SpecificElementVisitor
//   :FIXME: Add file description
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
#ifndef SPECIFICELEMENTVISITOR_H
#define SPECIFICELEMENTVISITOR_H

#include <list>

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/FlexibleCollimator.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Ring.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/Source.h"
#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/Vacuum.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "AbsBeamline/VerticalFFAMagnet.h"

#ifdef ENABLE_OPAL_FEL
#include "AbsBeamline/Undulator.h"
#endif

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

template <class ELEM1, class ELEM2>
struct CastsTrait {
    typedef std::list<const ELEM1*> ElementList_t;

    static void apply(ElementList_t &, const ELEM2 &)
    { }
};

template <class ELEM>
struct CastsTrait<ELEM,ELEM> {
    typedef std::list<const ELEM*> ElementList_t;

    static void apply(ElementList_t &allElements, const ELEM &element)
    {
        allElements.push_back(dynamic_cast<const ELEM*>(&element));
    }
};

template <class ELEM>
class SpecificElementVisitor: public BeamlineVisitor {
public:
    SpecificElementVisitor(const Beamline &beamline);

    virtual void execute();

    /// Apply the algorithm to a beam line.
    virtual void visitBeamline(const Beamline &);
    
    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator &);

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &);

    /// Apply the algorithm to a closed orbit corrector.
    virtual void visitCorrector(const Corrector &);

    /// Apply the algorithm to an cyclotron
    virtual void visitCyclotron(const Cyclotron &);

    /// Apply the algorithm to a degrader.
    virtual void visitDegrader(const Degrader &);

    /// Apply the algorithm to a drift.
    virtual void visitDrift(const Drift &);

    /// Apply the algorithm to a FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator &);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &);

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor &);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT &);

    /// Apply the algorithm to an arbitrary straight multipole.
    virtual void visitMultipoleTStraight(const MultipoleTStraight &);

    /// Apply the algorithm to an arbitrary curved multipole of constant radius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &);

    /// Apply the algorithm to an arbitrary curved multipole of variable radius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &);

    /// Apply the algorithm to an offset (placement).
    virtual void visitOffset(const Offset &);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend3D(const RBend3D &);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &);

    /// Apply the algorithm to a ring.
    virtual void visitRing(const Ring &);

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &);

    /// Apply the algorithm to a sector bend with 3D field map.
    virtual void visitSBend3D(const SBend3D &);

    /// Apply the algorithm to a scaling FFA magnet.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet &);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &);

    /// Apply the algorithm to a source.
    virtual void visitSource(const Source &);

   /// Apply the algorithm to a particle stripper.
    virtual void visitStripper(const Stripper &);

    /// Apply the algorithm to a traveling wave
    virtual void visitTravelingWave(const TravelingWave &);

#ifdef ENABLE_OPAL_FEL
    /// Apply the algorithm to an undulator.
    virtual void visitUndulator(const Undulator &);
#endif

    /// Apply the algorithm to a vacuum space.
    virtual void visitVacuum(const Vacuum &);

    /// Apply the algorithm to a variable RF cavity.
    virtual void visitVariableRFCavity(const VariableRFCavity &vcav);

    /// Apply the algorithm to a variable RF cavity with Fringe Field..
    virtual void visitVariableRFCavityFringeField
                                      (const VariableRFCavityFringeField &vcav);

    /// Apply the algorithm to a vertical FFA magnet.
    virtual void visitVerticalFFAMagnet(const VerticalFFAMagnet &);


    size_t size() const;

    typedef std::list<const ELEM*> ElementList_t;
    typedef typename ElementList_t::iterator iterator_t;
    typedef typename ElementList_t::const_iterator const_iterator_t;

    typedef typename ElementList_t::reference reference_t;
    typedef typename ElementList_t::const_reference const_reference_t;

    iterator_t begin();
    const_iterator_t begin() const;

    iterator_t end();
    const_iterator_t end() const;

    reference_t front();
    const_reference_t front() const;

private:
    ElementList_t allElementsOfTypeE;
};

template<class ELEM>
SpecificElementVisitor<ELEM>::SpecificElementVisitor(const Beamline &beamline):
    BeamlineVisitor(),
    allElementsOfTypeE()
{
    beamline.iterate(*this, false);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::execute()
{ }

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitBeamline(const Beamline &element) {
    element.iterate(*this, false);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitCCollimator(const CCollimator &element) {
    CastsTrait<ELEM, CCollimator>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitComponent(const Component &element) {
    CastsTrait<ELEM, Component>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitCorrector(const Corrector &element) {
    CastsTrait<ELEM, Corrector>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitCyclotron(const Cyclotron &element) {
    CastsTrait<ELEM, Cyclotron>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitDegrader(const Degrader &element) {
    CastsTrait<ELEM, Degrader>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitDrift(const Drift &element) {
    CastsTrait<ELEM, Drift>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitFlaggedElmPtr(const FlaggedElmPtr &element) {
    const ElementBase* wrappedElement = element.getElement();
    wrappedElement->accept(*this);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitFlexibleCollimator(const FlexibleCollimator &element) {
    CastsTrait<ELEM, FlexibleCollimator>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMarker(const Marker &element) {
    CastsTrait<ELEM, Marker>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMonitor(const Monitor &element) {
    CastsTrait<ELEM, Monitor>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMultipole(const Multipole &element) {
    CastsTrait<ELEM, Multipole>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMultipoleT(const MultipoleT &element) {
    CastsTrait<ELEM, MultipoleT>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMultipoleTStraight(const MultipoleTStraight &element) {
    CastsTrait<ELEM, MultipoleTStraight>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius &element) {
    CastsTrait<ELEM, MultipoleTCurvedConstRadius>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius &element) {
    CastsTrait<ELEM, MultipoleTCurvedVarRadius>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitOffset(const Offset &element) {
    CastsTrait<ELEM, Offset>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitProbe(const Probe &element) {
    CastsTrait<ELEM, Probe>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitRBend(const RBend &element) {
    CastsTrait<ELEM, RBend>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitRBend3D(const RBend3D &element) {
    CastsTrait<ELEM, RBend3D>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitRFCavity(const RFCavity &element) {
    CastsTrait<ELEM, RFCavity>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitRing(const Ring &element) {
    CastsTrait<ELEM, Ring>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitSBend(const SBend &element) {
    CastsTrait<ELEM, SBend>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitSBend3D(const SBend3D &element) {
    CastsTrait<ELEM, SBend3D>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitScalingFFAMagnet(const ScalingFFAMagnet &element) {
    CastsTrait<ELEM, ScalingFFAMagnet>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitSeptum(const Septum &element) {
    CastsTrait<ELEM, Septum>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitSolenoid(const Solenoid &element) {
    CastsTrait<ELEM, Solenoid>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitSource(const Source &element) {
    CastsTrait<ELEM, Source>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitStripper(const Stripper &element) {
    CastsTrait<ELEM, Stripper>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitTravelingWave(const TravelingWave &element) {
    CastsTrait<ELEM, TravelingWave>::apply(allElementsOfTypeE, element);
}

#ifdef ENABLE_OPAL_FEL
template<class ELEM>
void SpecificElementVisitor<ELEM>::visitUndulator(const Undulator &element) {
    CastsTrait<ELEM, Undulator>::apply(allElementsOfTypeE, element);
}
#endif

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitVacuum(const Vacuum &element) {
    CastsTrait<ELEM, Vacuum>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitVariableRFCavity
                                             (const VariableRFCavity &element) {
    CastsTrait<ELEM, VariableRFCavity>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitVariableRFCavityFringeField
                                  (const VariableRFCavityFringeField &element) {
    CastsTrait<ELEM, VariableRFCavityFringeField>::apply(allElementsOfTypeE,
                                                         element);
}

template<class ELEM>
void SpecificElementVisitor<ELEM>::visitVerticalFFAMagnet(const VerticalFFAMagnet &element) {
    CastsTrait<ELEM, VerticalFFAMagnet>::apply(allElementsOfTypeE, element);
}

template<class ELEM>
size_t SpecificElementVisitor<ELEM>::size() const{
    return allElementsOfTypeE.size();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::iterator_t SpecificElementVisitor<ELEM>::begin(){
    return allElementsOfTypeE.begin();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::const_iterator_t SpecificElementVisitor<ELEM>::begin() const{
    return allElementsOfTypeE.begin();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::iterator_t SpecificElementVisitor<ELEM>::end(){
    return allElementsOfTypeE.end();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::const_iterator_t SpecificElementVisitor<ELEM>::end() const{
    return allElementsOfTypeE.end();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::reference_t SpecificElementVisitor<ELEM>::front() {
    return allElementsOfTypeE.front();
}

template<class ELEM>
typename SpecificElementVisitor<ELEM>::const_reference_t SpecificElementVisitor<ELEM>::front() const{
    return allElementsOfTypeE.front();
}

#endif
