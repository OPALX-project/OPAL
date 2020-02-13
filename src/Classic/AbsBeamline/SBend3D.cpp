/*
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "AbsBeamline/SBend3D.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"

SBend3D::SBend3D(const std::string &name)
        : Component(name), map_m(NULL),
         planarArcGeometry_m(1., 1.), fieldUnits_m(1.), lengthUnits_m(1.),
         polyOrder_m(1), smoothOrder_m(1), dummy() {
    setElType(isDrift);
}

SBend3D::SBend3D(const SBend3D &right)
        : Component(right), map_m(NULL),
          planarArcGeometry_m(right.planarArcGeometry_m),
          fieldUnits_m(right.fieldUnits_m), lengthUnits_m(right.lengthUnits_m),
          polyOrder_m(right.polyOrder_m), smoothOrder_m(right.smoothOrder_m),
          dummy() {
    RefPartBunch_m = right.RefPartBunch_m;
    if (right.map_m != NULL)
        map_m = new SectorMagneticFieldMap(*right.map_m);
    setElType(isDrift);
}

SBend3D::~SBend3D() {
    delete map_m;
}

ElementBase* SBend3D::clone() const {
    return new SBend3D(*this);
}


EMField &SBend3D::getField() {
    return dummy;
}

const EMField &SBend3D::getField() const {
    return dummy;
}

bool SBend3D::apply(const size_t &i, const double &t,
                    Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool SBend3D::apply(const Vector_t &R, const Vector_t &/*P*/,
                    const double &/*t*/, Vector_t &E, Vector_t &B) {
    //std::cerr << "ROGERS SBend3D::apply " << R << " " << B << std::endl;
    return map_m->getFieldstrength(R, E, B);
}

void SBend3D::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void SBend3D::finalise() {
    RefPartBunch_m = NULL;
}

bool SBend3D::bends() const {
    return true;
}

BGeometryBase& SBend3D::getGeometry() {
    return planarArcGeometry_m;
}

const BGeometryBase& SBend3D::getGeometry() const {
    return planarArcGeometry_m;
}

void SBend3D::setFieldMapFileName(std::string name) {
    delete map_m;
    map_m = NULL;
    
    if (name != "") {
        map_m = new SectorMagneticFieldMap(
                        name,
                        "Dipole",
                        lengthUnits_m,
                        fieldUnits_m,
                        polyOrder_m,
                        smoothOrder_m);
        double r_curv = (map_m->getPolarBoundingBoxMax()[0]+
                         map_m->getPolarBoundingBoxMin()[0])/2.;
        double delta_phi = map_m->getDeltaPhi();
        planarArcGeometry_m.setElementLength(r_curv*delta_phi);
        planarArcGeometry_m.setCurvature(1./r_curv);
    }
}

void SBend3D::accept(BeamlineVisitor& visitor) const {
    visitor.visitSBend3D(*this);
}
