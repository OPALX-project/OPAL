//
// Class AmrObject
//   The AMR interface to OPAL. A new AMR library needs
//   to inherit from this class in order to work properly
//   with OPAL. Among other things it specifies the refinement
//   strategies.
//
// Copyright (c) 2016 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#include "Amr/AmrObject.h"

#include "Utilities/OpalException.h"

AmrObject::AmrObject()
    : AmrObject(TaggingCriteria::CHARGE_DENSITY, 0.75, 1.0e-15)
{ }


AmrObject::AmrObject(TaggingCriteria tagging,
                     double scaling,
                     double chargedensity)
    : tagging_m(tagging)
    , scaling_m(scaling)
    , chargedensity_m(chargedensity)
    , maxNumPart_m(1)
    , minNumPart_m(1)
    , refined_m(false)
    , amrSolveTimer_m(IpplTimings::getTimer("AMR solve"))
    , amrRegridTimer_m(IpplTimings::getTimer("AMR regrid"))
{ }


AmrObject::~AmrObject()
{ }


void AmrObject::setTagging(TaggingCriteria tagging) {
    tagging_m = tagging;
}


void AmrObject::setTagging(const std::string& tagging) {
    static const std::map<std::string, TaggingCriteria> stringTaggingCriteria_s = {
        {"CHARGE_DENSITY",     TaggingCriteria::CHARGE_DENSITY},
        {"POTENTIAL",          TaggingCriteria::POTENTIAL},
        {"EFIELD",             TaggingCriteria::EFIELD},
        {"MOMENTA",            TaggingCriteria::MOMENTA},
        {"MIN_NUM_PARTICLES",  TaggingCriteria::MIN_NUM_PARTICLES},
        {"MAX_NUM_PARTICLES",  TaggingCriteria::MAX_NUM_PARTICLES }
    };

    if (stringTaggingCriteria_s.count(tagging)) {
        tagging_m = stringTaggingCriteria_s.at(tagging);
    } else {
        throw OpalException("AmrObject::setTagging",
                            "Not supported refinement criteria.\n"
                            "Check the accepted values: "
                            "[CHARGE_DENSITY | POTENTIAL | EFIELD | "
                            "MOMENTA | MIN_NUM_PARTICLES | MAX_NUM_PARTICLES].");
    }
}


void AmrObject::setScalingFactor(double scaling) {
    scaling_m = scaling;
}


void AmrObject::setChargeDensity(double chargedensity) {
    chargedensity_m = chargedensity;
}


void AmrObject::setMaxNumParticles(size_t maxNumPart) {
    maxNumPart_m = maxNumPart;
}


void AmrObject::setMinNumParticles(size_t minNumPart) {
    minNumPart_m = minNumPart;
}


const bool& AmrObject::isRefined() const {
    return refined_m;
}


std::string AmrObject::getTaggingString(int number) {
    std::string tagging;
    switch ( number ) {
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::CHARGE_DENSITY):
            tagging = "CHARGE_DENSITY";
            break;
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::POTENTIAL):
            tagging = "POTENTIAL";
            break;
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::EFIELD):
            tagging = "EFIELD";
            break;
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::MOMENTA):
            tagging = "MOMENTA";
            break;
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::MIN_NUM_PARTICLES):
            tagging = "MIN_NUM_PARTICLES";
            break;
        case static_cast<std::underlying_type_t<TaggingCriteria>>(TaggingCriteria::MAX_NUM_PARTICLES):
            tagging = "MAX_NUM_PARTICLES";
            break;
        default:
            throw OpalException("AmrObject::getTaggingString",
                                "Only numbers between 0 and 5 allowed.");
    }
    return tagging;
}
