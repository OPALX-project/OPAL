//
// Units
//   A namespace defining various units conversions.
//
// Copyright (c) 2021, Carl Jolly, STFC, RAL
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
#ifndef CLASSIC_Units_HH
#define CLASSIC_Units_HH

#include "Physics.h"

namespace Units {

    // metre to millimetre
    constexpr double m2mm = 1e+3;

    // millimetre to metre
    constexpr double mm2m = 1 / m2mm;

    // meter to centimeter
    constexpr double m2cm = 1e+2;

    // centimeter to meter
    constexpr double cm2m = 1 / m2cm;

    //seconds to micro seconds
    constexpr double s2us = 1e+6;

    //micro seconds to seconds
    constexpr double us2s = 1 / s2us;

    //seconds to nano seconds
    constexpr double s2ns = 1e+9;

    //nano seconds to seconds
    constexpr double ns2s = 1 / s2ns;

    //seconds to pico seconds
    constexpr double s2ps = 1e+12;

    //pico seconds to seconds
    constexpr double ps2s = 1 / s2ps;

    //Tesla to kilo Gauss
    constexpr double T2kG = 1e+1;

    //kilo Gauss to Tesla
    constexpr double kG2T = 1 / T2kG;

    //kilo volts to volts
    constexpr double kV2V = 1e+3;

    //volts to kilo volts
    constexpr double V2kV = 1 / kV2V;

    //GeV to eV
    constexpr double GeV2eV = 1e+9;

    //eV to GeV
    constexpr double eV2GeV = 1 / GeV2eV;

    // MeV to eV
    constexpr double MeV2eV = 1e+6;

    // eV to MeV
    constexpr double eV2MeV = 1 / MeV2eV;

    // GeV to MeV
    constexpr double GeV2MeV = GeV2eV * eV2MeV;

    // MeV to GeV
    constexpr double MeV2GeV = 1 / GeV2MeV;

    // keV to eV
    constexpr double keV2eV = 1e+3;

    // eV to keV
    constexpr double eV2keV = 1 / keV2eV;

    // GeV to keV
    constexpr double GeV2keV = GeV2eV * eV2keV;

    // keV t to GeV
    constexpr double keV2GeV = keV2eV * eV2GeV;

    // MeV to keV
    constexpr double MeV2keV = MeV2eV * eV2keV;

    // keV to MeV
    constexpr double keV2MeV = 1 / MeV2keV;

    //GeV/c^2 to V*C*s^2/m^2 (ie GeV/c^2 to kg conversion)
    constexpr double GeV2kg = 1.0e+9 * Physics::q_e / Physics::c / Physics::c;

    //V*C*s^2/m^2 to GeV/c^2 (ie GeV/c^2 to kg conversion)
    constexpr double kg2GeV = 1 / GeV2kg;

    //MHz to Hz
    constexpr double MHz2Hz = 1e+6;

    //Hz to MHz
    constexpr double Hz2MHz = 1 / MHz2Hz;

    //GHz to Hz
    constexpr double GHz2Hz = 1e+9;

    //Hz to GHz
    constexpr double Hz2GHz = 1 / GHz2Hz;

    //V/m to MV/m
    constexpr double Vpm2MVpm = 1e-6;

    //MV/m to V/m
    constexpr double MVpm2Vpm = 1 / Vpm2MVpm;

    //A to mA
    constexpr double A2mA = 1e+3;

    //mA to A
    constexpr double mA2A = 1 / A2mA;

    // rad to mrad
    constexpr double rad2mrad = 1e3;

    // mrad to rad
    constexpr double mrad2rad = 1 / rad2mrad;

    // deg to rad
    constexpr double deg2rad = Physics::pi / 180;

    // rad to deg
    constexpr double rad2deg = 1 / deg2rad;
};

#endif // CLASSIC_Units_HH