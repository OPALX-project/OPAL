//
// Class Astra1DDynamic_fast
//
// This class provides a reader for Astra style field maps. It pre-computes the field
// on a lattice to increase the performance during simulation.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#ifndef CLASSIC_AstraFIELDMAP1DDYNAMICFAST_HH
#define CLASSIC_AstraFIELDMAP1DDYNAMICFAST_HH

#include "Fields/Astra1D_fast.h"

class Astra1DDynamic_fast: public Astra1D_fast {

public:
    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & F);

private:
    Astra1DDynamic_fast(std::string aFilename);
    virtual ~Astra1DDynamic_fast();

    virtual void readMap();

    bool readFileHeader(std::ifstream &file);
    int stripFileHeader(std::ifstream &file);

    double frequency_m;
    double xlrep_m;

    friend class Fieldmap;
};

#endif