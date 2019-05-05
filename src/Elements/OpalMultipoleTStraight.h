/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
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


#ifndef OPAL_OPALMULTIPOLET_STRAIGHT_HH
#define OPAL_OPALMULTIPOLET_STRAIGHT_HH

#include "Elements/OpalElement.h"
#include "AbsBeamline/MultipoleTStraight.h"

/** OpalMultipoleTStraight provides user interface information for the MultipoleTStraight object
 *
 *  Defines input parameters
 */
class OpalMultipoleTStraight: public OpalElement {

public:

    // The attributes of class OpalMultipoleTStraight
    enum {
        TP = COMMON,     // Transverse field components
        RFRINGE,         // Length of right fringe field
        LFRINGE,         // Length of left fringe field
        HAPERT,          // Aperture horizontal dimension
        VAPERT,          // Aperture vertical dimension
        MAXFORDER,       // Maximum order in the field expansion
        ROTATION,        // Rotation angle about central axis for skew elements
        EANGLE,          // Entrance angle
        BBLENGTH,        // Distance between centre of magnet and entrance
        SIZE             // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalMultipoleTStraight();

    /** Destructor does nothing */
    virtual ~OpalMultipoleTStraight();

    /** Inherited copy constructor */
    virtual OpalMultipoleTStraight *clone(const std::string &name);

    /** Fill in all registered attributes
     *
     *  Just calls fillRegisteredAttributes on the base class
     */
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /** Update the MultipoleT with new parameters from UI parser */
    virtual void update();

    void print(std::ostream &os) const;

  private:
    // Not implemented.
    OpalMultipoleTStraight(const OpalMultipoleTStraight &);
    void operator=(const OpalMultipoleTStraight &);

    // Clone constructor.
    OpalMultipoleTStraight(const std::string &name, OpalMultipoleTStraight *parent);
};

#endif
