//
// Class H5Writer
//   Interface for H5 writers.
//
// Copyright (c) 2019-2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "H5Writer.h"

H5Writer::H5Writer(H5PartWrapper* h5wrapper, bool restart)
    : H5PartTimer_m(IpplTimings::getTimer("Write H5-File"))
    , h5wrapper_m(h5wrapper)
    , H5call_m(0)
{
    if ( !restart ) {
        h5wrapper->writeHeader();
    }
    h5wrapper->close();
}


void H5Writer::writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[]) {
    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("B-ref_x", FDext[0](0)),
        std::make_pair("B-ref_z", FDext[0](1)),
        std::make_pair("B-ref_y", FDext[0](2)),
        std::make_pair("E-ref_x", FDext[1](0)),
        std::make_pair("E-ref_z", FDext[1](1)),
        std::make_pair("E-ref_y", FDext[1](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);
}


int H5Writer::writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[], double /*meanEnergy*/,
                              double refPr, double refPt, double refPz,
                              double refR, double refTheta, double refZ,
                              double azimuth, double elevation, bool /*local*/) {

    if (beam->getTotalNum() < 3) return -1; // in single particle mode and tune calculation (2 particles) we do not need h5 data

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", refPr),
        std::make_pair("REFPT", refPt),
        std::make_pair("REFPZ", refPz),
        std::make_pair("REFR", refR),
        std::make_pair("REFTHETA", refTheta),
        std::make_pair("REFZ", refZ),
        std::make_pair("AZIMUTH", azimuth),
        std::make_pair("ELEVATION", elevation),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x",  FDext[2](0)),
        std::make_pair("B-ref_z",  FDext[2](1)),
        std::make_pair("B-ref_y",  FDext[2](2)),
        std::make_pair("E-ref_x",  FDext[3](0)),
        std::make_pair("E-ref_z",  FDext[3](1)),
        std::make_pair("E-ref_y",  FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    ++ H5call_m;
    return H5call_m - 1;
}