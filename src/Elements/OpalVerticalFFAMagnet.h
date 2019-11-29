/* 
 *  Copyright (c) 2017, Chris Rogers
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

#ifndef OPAL_OPALVERTICALFFAMAGNET_H
#define OPAL_OPALVERTICALFFAMAGNET_H

#include "Elements/OpalBend.h"

/** OpalVerticalFFAMagnet provides the user interface for the VERTICALFFA object
 */
class OpalVerticalFFAMagnet : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        B0 = COMMON,
        FIELD_INDEX,
        WIDTH,
        MAX_HORIZONTAL_POWER,
        END_LENGTH,
        CENTRE_LENGTH,
        BB_LENGTH,
        HEIGHT_NEG_EXTENT,
        HEIGHT_POS_EXTENT,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalVerticalFFAMagnet();

    /** Destructor does nothing */
    virtual ~OpalVerticalFFAMagnet();

    /** Inherited copy constructor */
    virtual OpalVerticalFFAMagnet *clone(const std::string &name);

    /** Fill in all registered attributes
     *
     *  Just calls fillRegisteredAttributes on the base class
     */
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /** Update the VerticalFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalVerticalFFAMagnet(const OpalVerticalFFAMagnet &);
    void operator=(const OpalVerticalFFAMagnet &);

    // Clone constructor.
    OpalVerticalFFAMagnet(const std::string &name, OpalVerticalFFAMagnet *parent);

    static std::string docstring_m;
};

#endif // OPAL_OPALVERTICALFFAMAGNET_H

