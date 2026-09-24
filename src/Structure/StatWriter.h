//
// Class StatWriter
//   This class writes bunch statistics (*.stat).
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_STAT_WRITER_H
#define OPAL_STAT_WRITER_H

#include "StatBaseWriter.h"

class StatWriter : public StatBaseWriter {

public:
    typedef std::vector<std::pair<std::string, unsigned int> > losses_t;

    StatWriter(const std::string& fname, bool restart);

    using SDDSWriter::write;
    /** \brief Write statistical data.
     *
     * Writes statistical beam data to proper output file. This is information such as RMS beam parameters
     * etc.
     *
     * Also gathers and writes load balancing data to load balance statistics file.
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     * \param timeStep The integration time step [s] written as 'dt'. This is beam->getdT(), except in
     * OPAL-cycl, which integrates with getdT() times the harmonic number.
     */
    void write(const PartBunchBase<double, 3> *beam, Vector_t FDext[],
               const double& timeStep,
               const losses_t &losses = losses_t(), const double& azimuth = -1,
               const size_t npOutside = 0);

private:
    void fillHeader(const PartBunchBase<double, 3> *beam, const losses_t &losses = losses_t());

    /// The header has the columns of FIELDSOLVER, MESHFIT=CORE. Like all columns they are
    /// fixed by the first write: a later TRACK/RUN with MESHFIT=ALL writes zeros into them,
    /// and one with MESHFIT=CORE goes without them if the first had MESHFIT=ALL.
    bool hasMeshFitColumns_m = false;
    /// The warning that the core selection is not written has been given.
    bool warnedMeshFitColumns_m = false;
};

#endif
