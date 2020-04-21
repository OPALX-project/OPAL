//
// Class OpalSample
//   The SAMPLING definition.
//   A SAMPLING definition is used to run the optimizer in sample mode.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_SAMPLE_H
#define OPAL_SAMPLE_H

#include <string>
#include <memory>
#include "AbstractObjects/Definition.h"

#include "Sample/SamplingMethod.h"

class OpalSample: public Definition {

public:

    /// Exemplar constructor.
    OpalSample();

    virtual ~OpalSample() {};

    /// Make clone.
    virtual OpalSample *clone(const std::string &name);

    /// Check the OpalSample data.
    virtual void execute();

    /// Find sampling method
    static OpalSample *find(const std::string &name);

    void initialize(const std::string &dvarName,
                    double lower,
                    double upper,
                    size_t modulo = 1,
                    bool sequence = false);

    std::string getVariable() const;

    unsigned int getSize() const;

    std::shared_ptr<SamplingMethod> sampleMethod_m;

private:

    ///@{ Not implemented.
    OpalSample  (const OpalSample &) = delete;
    void operator=(const OpalSample &) = delete;
    ///@}
    /// Private copy constructor, called by clone
    OpalSample(const std::string &name, OpalSample *parent);

    unsigned int size_m;
};

inline
unsigned int OpalSample::getSize() const{
    return size_m;
}
#endif
