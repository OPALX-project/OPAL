//
// Class OpalVariableRFCavityFringeField
//   The class provides the user interface for
//   the VARIABLE_RF_CAVITY_FRINGE_FIELD object.
//
// Copyright (c) 2018 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#ifndef OPAL_OPALVARIABLERFCAVITYFRINGEFIELD_H
#define OPAL_OPALVARIABLERFCAVITYFRINGEFIELD_H

#include "Elements/OpalElement.h"

class OpalVariableRFCavityFringeField: public OpalElement {

public:
    /** enum maps string to integer value for UI definitions */
    enum {
        PHASE_MODEL = COMMON,
        AMPLITUDE_MODEL,
        FREQUENCY_MODEL,
        WIDTH,
        HEIGHT,
        CENTRE_LENGTH,
        END_LENGTH,
        CAVITY_CENTRE,
        MAX_ORDER,
        SIZE // size of the enum
    };

    /** Copy constructor **/
    OpalVariableRFCavityFringeField(const std::string& name,
                                    OpalVariableRFCavityFringeField* parent);

    /** Default constructor **/
    OpalVariableRFCavityFringeField();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
    **/
    OpalVariableRFCavityFringeField* clone();

    /** Inherited copy constructor
     *
     *  Call on a base class to instantiate an object of derived class's type
     */
    virtual OpalVariableRFCavityFringeField* clone(const std::string& name);

    /** Destructor does nothing */
    virtual ~OpalVariableRFCavityFringeField();

    /** Update the OpalVariableRFCavity with new parameters from UI parser */
    virtual void update();

private:
    // Not implemented.
    OpalVariableRFCavityFringeField(const OpalVariableRFCavityFringeField&);
    void operator=(const OpalVariableRFCavityFringeField&);

    /** Check for conversion to unsigned int */
    inline static size_t convertToUnsigned(double value, const std::string& name);

    static const std::string doc_string;
    static constexpr double unsignedTolerance = 1e-9;
};

#endif // OPAL_OPALVARIABLERFCAVITY_H
