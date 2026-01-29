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

#ifndef OPAL_OpalPolynomialTimeDependence_HH
#define OPAL_OpalPolynomialTimeDependence_HH

#include "Elements/OpalElement.h"

/** OpalPolynomialTimeDependence provides UI wrapper for the
 *  PolynomialTimeDependence
 */

class OpalPolynomialTimeDependence : public OpalElement {
  public:
    /** Enumeration maps to UI parameters */
    enum {
        P0 = COMMON,
        P1,
        P2,
        P3,
        SIZE // size of the enum
    };

    /** Define mapping from enum variables to string UI parameter names */
    OpalPolynomialTimeDependence();

    /** No memory allocated so does nothing */
    virtual ~OpalPolynomialTimeDependence();

    /** Inherited copy constructor */
    virtual OpalPolynomialTimeDependence *clone(const std::string &name);

    /** Receive parameters from the parser and hand them off to the
     *  PolynomialTimeDependence
     */
    void update();

    /** Calls print on the OpalElement */
    virtual void print(std::ostream &) const;
  private:
    // Not implemented.
    OpalPolynomialTimeDependence(const OpalPolynomialTimeDependence &);
    void operator=(const OpalPolynomialTimeDependence &);

    // Clone constructor.
    OpalPolynomialTimeDependence(const std::string &name,
                                 OpalPolynomialTimeDependence *parent);

    static const std::string doc_string;
};

#endif // OPAL_OpalPolynomialTimeDependence_HH
