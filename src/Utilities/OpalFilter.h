//
// Class OpalFilter
//   The class for the Filter Object.
//   A FILTER definition is used to define a filter which can be applied
//   to a 1D histogram in order to get rid of noise.
//
// Copyright (c) 2008 - 2022, Christof Metzger-Kraus
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
#ifndef OPAL_FILTER_HH
#define OPAL_FILTER_HH

#include "AbstractObjects/Definition.h"
#include "Filters/Filter.h"

#include <memory>
#include <vector>


class OpalFilter: public Definition {

public:
    /// Exemplar constructor.
    OpalFilter();

    virtual ~OpalFilter();

    /// Test if replacement is allowed.
    //  Can replace only by another OpalFilter
    virtual bool canReplaceBy(Object* object);

    /// Make clone.
    virtual OpalFilter* clone(const std::string& name);

    /// Check the OpalFilter data.
    virtual void execute();

    /// Find named FILTER.
    static OpalFilter* find(const std::string& name);

    /// Update the OpalFilter data.
    virtual void update();

    void print(std::ostream& os) const;

    void initOpalFilter();

    inline void apply(std::vector<double>& histogram);
    inline void calc_derivative(std::vector<double>& histogram, const double& hz);

    Filter* filter_m;

private:
    enum class FilterType: unsigned short {
        SAVITZKYGOLAY,
        FIXEDFFTLOWPASS,
        RELATIVEFFTLOWPASS,
        STENCIL
    };

    // Not implemented.
    OpalFilter(const OpalFilter&);
    void operator=(const OpalFilter&);

    // Clone constructor.
    OpalFilter(const std::string &name, OpalFilter* parent);
};

void OpalFilter::apply(std::vector<double>& histogram) {
    if (filter_m)
        filter_m->apply(histogram);
}

void OpalFilter::calc_derivative(std::vector<double>& histogram, const double& hz) {
    if (filter_m)
        filter_m->calc_derivative(histogram, hz);
}

inline std::ostream& operator<<(std::ostream& os, const OpalFilter& b) {
    b.print(os);
    return os;
}

#endif // OPAL_FILTER_HH
