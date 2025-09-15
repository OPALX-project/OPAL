/*
 *  Copyright (c) 2025, Jon Thompson
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

#include "SinusoidalTimeDependence.h"
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
        const auto angle = 2.0 * M_PI * f * time + p;
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
