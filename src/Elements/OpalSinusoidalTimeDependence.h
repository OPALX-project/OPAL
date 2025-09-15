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
