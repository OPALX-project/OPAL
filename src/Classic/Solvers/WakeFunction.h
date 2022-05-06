//
// Class WakeFunction
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef WAKEFUNCTION_HH
#define WAKEFUNCTION_HH

#include <string>
#include <vector>

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;

enum class WakeType: unsigned short {
    CSRWakeFunction,
    CSRIGFWakeFunction,
    GreenWakeFunction
};

class WakeFunction {

public:
    WakeFunction(std::string name, unsigned int n):
        nBins_m(n),
        name_m(name) { };

    virtual ~WakeFunction(){ };

    virtual void initialize(const ElementBase* /*elref*/) { };

    virtual void apply(PartBunchBase<double, 3>* bunch) = 0;

    virtual WakeType getType() const = 0;

    const std::string& getName() const { return name_m; }

protected:
    const unsigned int nBins_m;

private:
    const std::string name_m;
};

class LineDensity: public std::vector<double> {
public:
    LineDensity(int size = 0, double defaultValue = 0.0) : std::vector<double>(size, defaultValue) {}
    void getFirstDerivative(std::vector<double>& firstDerivative, const double& hz);
};

inline void LineDensity::getFirstDerivative(std::vector<double>& firstDerivative, const double& hz) {
    const size_t size = this->size();
    if(firstDerivative.size() != size)
        firstDerivative.resize(size, 0.0);

    firstDerivative[0] = ((*this)[1] - (*this)[0]) / hz;
    for(unsigned int i = 1; i + 1 < size; ++i)
        firstDerivative[i] = ((*this)[i + 1] - (*this)[i - 1]) / hz;
    firstDerivative[size - 1] = ((*this)[size - 1] - (*this)[size - 2]) / hz;
}

#endif // WAKEFUNCTION_HH
