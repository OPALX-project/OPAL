// ------------------------------------------------------------------------
// $RCSfile: BeamlineGeometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamlineGeometry
//   A class for implementing the composite geometry of a beamline.
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/BeamlineGeometry.h"
#include "Beamlines/Beamline.h"


// Class BeamlineGeometry
// ------------------------------------------------------------------------

BeamlineGeometry::BeamlineGeometry(const Beamline &line):
    itsLine(line)
{}


BeamlineGeometry::~BeamlineGeometry()
{}


double BeamlineGeometry::getArcLength() const {
    return itsLine.getArcLength();
}


double BeamlineGeometry::getElementLength() const {
    return itsLine.getElementLength();
}


double BeamlineGeometry::getOrigin() const {
    return itsLine.getOrigin();
}


double BeamlineGeometry::getEntrance() const {
    return itsLine.getEntrance();
}


double BeamlineGeometry::getExit() const {
    return itsLine.getExit();
}


Euclid3D BeamlineGeometry::getTransform(double fromS, double toS) const {
    return itsLine.getTransform(fromS, toS);
}


Euclid3D BeamlineGeometry::getTotalTransform() const {
    return itsLine.getTotalTransform();
}


Euclid3D BeamlineGeometry::getTransform(double s) const {
    return itsLine.getTransform(s);
}


Euclid3D BeamlineGeometry::getEntranceFrame() const {
    return itsLine.getEntranceFrame();
}


Euclid3D BeamlineGeometry::getExitFrame() const {
    return itsLine.getExitFrame();
}

