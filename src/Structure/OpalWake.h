//
// Class OpalWake
//   The class for the OPAL WAKE command.
//
// Copyright (c) 2008 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_Wake_HH
#define OPAL_Wake_HH

#include "AbstractObjects/Definition.h"

class ElementBase;
class WakeFunction;

class OpalWake: public Definition {

public:
    /// Exemplar constructor.
    OpalWake();

    virtual ~OpalWake();

    /// Test if replacement is allowed.
    //  Can replace only by another WAKE.
    virtual bool canReplaceBy(Object* object);

    /// Make clone.
    virtual OpalWake* clone(const std::string& name);

    /// Check the WAKE data.
    virtual void execute();

    /// Find named WAKE.
    static OpalWake* find(const std::string& name);

    /// Update the WAKE data.
    virtual void update();

    void print(std::ostream& os) const;

    int getNumberOfBins();

    void initWakefunction(const ElementBase& element);

    WakeFunction* wf_m;

private:
    enum class OpalWakeType: unsigned short {
        CSR,
        CSRIGF,
        LONGSHORTRANGE,
        TRANSVSHORTRANGE
    };

    // Not implemented.
    OpalWake(const OpalWake&);
    void operator=(const OpalWake&);

    // Clone constructor.
    OpalWake(const std::string& name, OpalWake* parent);
};

inline std::ostream& operator<<(std::ostream& os, const OpalWake& b) {
    b.print(os);
    return os;
}

#endif // OPAL_Wake_HH
