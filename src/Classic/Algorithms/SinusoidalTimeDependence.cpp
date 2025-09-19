//
// Class SinusoidalTimeDependence
//   A time dependence class that generates sine waves
//
// Copyright (c) 2025, Jon Thompson, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#include "SinusoidalTimeDependence.h"
#include <Physics/Physics.h>
#include "PolynomialTimeDependence.h"
#include "Utility/Inform.h"


SinusoidalTimeDependence::SinusoidalTimeDependence(
    const std::vector<double>& f, const std::vector<double>& p, const std::vector<double>& a,
    const std::vector<double>& o)
    : f_m(f), p_m(p), a_m(a), o_m(o) {
}

double SinusoidalTimeDependence::getValue(const double time) {
    double result{};
    for (size_t i = 0; i < f_m.size(); i++) {
        const auto f = f_m[i];
        auto p       = 0.0;
        auto a       = 1.0;
        auto o       = 0.0;
        if (i < p_m.size()) {
            p = p_m[i];
        }
        if (i < a_m.size()) {
            a = a_m[i];
        }
        if (i < o_m.size()) {
            o = o_m[i];
        }
        const auto angle = 2.0 * Physics::pi * f * time + p;
        result += a / 2.0 * std::sin(angle) + o;
    }
    return result;
}

SinusoidalTimeDependence* SinusoidalTimeDependence::clone() {
    return new SinusoidalTimeDependence(f_m, p_m, a_m, o_m);
}

Inform& SinusoidalTimeDependence::print(Inform& os) const {
    const Inform::FmtFlags_t ff = os.flags();
    os << std::scientific;
    os << "f=[";
    for (size_t i = 0; i < this->f_m.size(); i++) {
        if (i != 0) {
            os << ", ";
        }
        os << this->f_m[i];
    }
    os << "], p=[";
    for (size_t i = 0; i < this->p_m.size(); i++) {
        if (i != 0) {
            os << ", ";
        }
        os << this->p_m[i];
    }
    os << "], a=[";
    for (size_t i = 0; i < this->a_m.size(); i++) {
        if (i != 0) {
            os << ", ";
        }
        os << this->a_m[i];
    }
    os << "], o=[";
    for (size_t i = 0; i < this->o_m.size(); i++) {
        if (i != 0) {
            os << ", ";
        }
        os << this->o_m[i];
    }
    os << endl;
    os.flags(ff);
    return os;
}
