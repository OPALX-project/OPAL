//
// Class OpalSinusoidalTimeDependence
//   User interface for the time dependence class that generates sine waves
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

#ifndef OPAL_OPALSINUSOIDALTIMEDEPENDENCE_H
#define OPAL_OPALSINUSOIDALTIMEDEPENDENCE_H

#include "Elements/OpalElement.h"

/** The UI wrapper for SinusoidalTimeDependence */
class OpalSinusoidalTimeDependence final : public OpalElement {
public:
    /** Enumeration maps to UI parameters */
    enum {
        FREQUENCIES = COMMON,
        PHASE_OFFSETS,
        AMPLITUDES,
        DC_OFFSETS,
        SIZE // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalSinusoidalTimeDependence();

    /** Default destructor */
    ~OpalSinusoidalTimeDependence() override = default;

    /** Inherited copy constructor */
    OpalSinusoidalTimeDependence *clone(const std::string &name) override;

    /** Receive parameters from the parser and hand them off to the OpalSinusoidalTimeDependence*/
    void update() override;

    /** Calls print on the OpalElement */
    void print(std::ostream &) const override;

private:
    // Not implemented.
    OpalSinusoidalTimeDependence(const OpalSinusoidalTimeDependence&) = delete;
    OpalSinusoidalTimeDependence& operator=(const OpalSinusoidalTimeDependence &) = delete;

    // Clone constructor.
    OpalSinusoidalTimeDependence(const std::string &name, OpalSinusoidalTimeDependence *parent);

    static const std::string doc_string;
};

#endif  // OPAL_OPALSINUSOIDALTIMEDEPENDENCE_H
