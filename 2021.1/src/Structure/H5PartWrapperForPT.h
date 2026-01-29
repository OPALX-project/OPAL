//
// Class H5PartWrapperForPT
//   A class that manages all calls to H5Part for the Parallel-T tracker.
//
// Copyright (c) 200x-2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_H5PARTWRAPPERFORPT_H
#define OPAL_H5PARTWRAPPERFORPT_H

#include "Structure/H5PartWrapper.h"

#include "H5hut.h"

class H5PartWrapperForPT: public H5PartWrapper {
public:
    H5PartWrapperForPT(const std::string& fileName, h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapperForPT(const std::string& fileName, int restartStep, std::string sourceFile, h5_int32_t flags = H5_O_RDWR);
    virtual ~H5PartWrapperForPT();

    virtual void readHeader();
    virtual void readStep(PartBunchBase<double, 3>*, h5_ssize_t firstParticle, h5_ssize_t lastParticle);

    virtual void writeHeader();
    virtual void writeStep(PartBunchBase<double, 3>*, const std::map<std::string, double>& additionalStepAttributes);

    virtual bool predecessorIsSameFlavour() const;

private:
    void readStepHeader(PartBunchBase<double, 3>*);
    void readStepData(PartBunchBase<double, 3>*, h5_ssize_t, h5_ssize_t);

    void writeStepHeader(PartBunchBase<double, 3>*, const std::map<std::string, double>&);
    void writeStepData(PartBunchBase<double, 3>*);
};

inline
bool H5PartWrapperForPT::predecessorIsSameFlavour() const {
    return (predecessorOPALFlavour_m == "opal-t");
}

#endif //OPAL_H5PARTWRAPPERFORPT_H
