//
// Class AbstractAmrWriter
//   Abstract base class for writing AMR data to output files.
//
// Copyright (c) 2016 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef ABSTRACT_AMR_WRITER_H
#define ABSTRACT_AMR_WRITER_H

#include "Amr/AmrObject.h"
#include "Amr/AmrDefs.h"
#include "Algorithms/AmrPartBunch.h"

class AbstractAmrWriter {
    
public:
    /*!
     * @param rho is the charge density on all levels
     * @param phi is the electrostatic potential on all levels
     * @param efield are the electric field components on all levels
     * @param refRatio are the refinement ratios among the levels
     * @param geom are the geometries of all levels
     * @param nLevel available
     * @param time specifies the step.
     * @param scale used for mapping
     */
    virtual void writeFields(const amr::AmrScalarFieldContainer_t& rho,
                             const amr::AmrScalarFieldContainer_t& phi,
                             const amr::AmrVectorFieldContainer_t& efield,
                             const amr::AmrIntArray_t& refRatio,
                             const amr::AmrGeomContainer_t& geom,
                             const int& nLevel,
                             const double& time,
                             const double& scale = 1.0) = 0;
    
    /*!
     * @param bunch_p
     * @param time
     * @param scale used for mapping
     */
    virtual void writeBunch(const AmrPartBunch* bunch_p,
                            const double& time,
                            const double& scale = 1.0) = 0;
    
    virtual ~AbstractAmrWriter() { }
    
};

#endif
