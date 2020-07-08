//
// Class Mapper
//   Build transfer map. The visitor class for building a VpsMap for a beamline
//   using a thin-lens approximation for all elements.
//   A visitor class implementing the default behaviour for all
//   visitors capable of tracking a transfer map through a beam line.
//   It implements some default behaviour for such visitors.
//   [P]
//   Phase space coordinates (in this order):
//   [DL]
//   [DT]x:[DD]
//     horizontal displacement (metres).
//   [DT]p_x/p_r:[DD]
//     horizontal canonical momentum (no dimension).
//   [DT]y:[DD]
//     vertical displacement (metres).
//   [DT]p_y/p_r:[DD]
//     vertical canonical momentum (no dimension).
//   [DT]delta_p/p_r:[DD]
//     relative momentum error (no dimension).
//   [DT]v*delta_t:[DD]
//     time difference delta_t w.r.t. the reference frame which moves with
//     uniform velocity
//   [P]
//     v_r = c*beta_r = p_r/m
//   [P]
//     along the design orbit, multiplied by the instantaneous velocity v of
//     the particle (metres).
//   [/DL]
//   Where
//   [DL]
//   [DT]p_r:[DD]
//     is the constant reference momentum defining the reference frame velocity.
//   [DT]m:[DD]
//     is the rest mass of the particles.
//   [/DL]
//   Other units used:
//   [DL]
//   [DT]reference momentum:[DD]
//     electron-volts.
//   [DT]accelerating voltage:[DD]
//     volts.
//   [DT]separator voltage:[DD]
//     volts.
//   [DT]frequencies:[DD]
//     hertz.
//   [DT]phase lags:[DD]
//     multiples of (2*pi).
//   [/DL]
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
#include "Algorithms/Mapper.h"
#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/LinearMap.h"
#include "FixedAlgebra/TransportMap.h"
#include "FixedAlgebra/FTps.h"

typedef FTps<double, 6> Series;
typedef FVps<double, 6> Map;


Mapper::Mapper(const Beamline &beamline, const PartData &reference,
               bool backBeam, bool backTrack):
    AbstractMapper(beamline, reference, backBeam, backTrack),
    itsMap()
{}


Mapper::~Mapper()
{}


void Mapper::getMap(LinearMap<double, 6> &map) const {
    map = LinearMap<double, 6>(itsMap);
}


void Mapper::getMap(TransportMap<double, 6> &map) const {
    map = TransportMap<double, 6>(itsMap);
}


void Mapper::getMap(Map &map) const {
    map = itsMap;
}


void Mapper::setMap(const LinearMap<double, 6> &map) {
    itsMap = Map(map);
}


void Mapper::setMap(const TransportMap<double, 6> &map) {
    itsMap = Map(map);
}


void Mapper::setMap(const Map &map) {
    itsMap = map;
}


void Mapper::visitComponent(const Component &comp) {
    comp.trackMap(itsMap, itsReference, back_beam, back_track);
}


void Mapper::visitPatch(const Patch &patch) {
    Euclid3D transform = patch.getPatch();
    if(back_path) transform = Inverse(transform);
    applyTransform(transform);
}


void Mapper::applyDrift(double length) {
    double kin = itsReference.getM() / itsReference.getP();
    double refTime = length / itsReference.getBeta();

    Series px = itsMap[PX];
    Series py = itsMap[PY];
    Series pt = itsMap[PT] + 1.0;
    Series pz = sqrt(pt * pt - px * px - py * py);
    Series E = sqrt(pt * pt + kin * kin);

    itsMap[X] += length * px / pz;
    itsMap[Y] += length * py / pz;
    itsMap[T] += pt * (refTime / E - length / pz);
}


void Mapper::applyThinMultipole
(const BMultipoleField &field, double scale) {
    int order = field.order();

    if(order > 0) {
        Series x = itsMap[X];
        Series y = itsMap[Y];
        Series kx = + field.normal(order);
        Series ky = - field.skew(order);

        while(--order > 0) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            kx = kxt + field.normal(order);
            ky = kyt - field.skew(order);
        }

        itsMap[PX] -= kx * scale;
        itsMap[PY] += ky * scale;
    }
}


void Mapper::applyThinSBend
(const BMultipoleField &field, double scale, double h) {
    Series As = buildSBendVectorPotential(field, h) * scale;

    // These substitutions work because As depends on x and y only,
    // and not on px or py.
    itsMap[PX] -= As.derivative(X).substitute(itsMap);
    itsMap[PY] -= As.derivative(Y).substitute(itsMap);
}


void Mapper::applyTransform(const Euclid3D &euclid, double refLength) {
    if(! euclid.isIdentity()) {
        double kin = itsReference.getM() / itsReference.getP();
        double refTime = refLength / itsReference.getBeta();
        Series px = itsMap[PX];
        Series py = itsMap[PY];
        Series pt = itsMap[PT] + 1.0;
        Series pz = sqrt(pt * pt - px * px - py * py);

        itsMap[PX] = euclid.M(0, 0) * px + euclid.M(1, 0) * py + euclid.M(2, 0) * pz;
        itsMap[PY] = euclid.M(0, 1) * px + euclid.M(1, 1) * py + euclid.M(2, 1) * pz;
        pz = euclid.M(0, 2) * px + euclid.M(1, 2) * py + euclid.M(2, 2) * pz;

        Series x = itsMap[X] - euclid.getX();
        Series y = itsMap[Y] - euclid.getY();
        Series x2 =
            euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
        Series y2 =
            euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
        Series s2 =
            euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
        Series sByPz = s2 / pz;

        Series E = sqrt(pt * pt + kin * kin);
        itsMap[X] = x2 - sByPz * itsMap[PX];
        itsMap[Y] = y2 - sByPz * itsMap[PY];
        itsMap[T] += pt * (refTime / E  + sByPz);
    }
}
