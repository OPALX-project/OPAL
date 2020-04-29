//
// Class AmrYtWriter
//   This concrete class writes output files that are readable by yt
//   (cf. http://yt-project.org/). We have a fork of yt in
//   the repository at https://gitlab.psi.ch/frey_m/yt.
//   The functions of this class are copied from AMReX and modified to fit
//   our needs.
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
#ifndef AMR_YT_WRITER_H
#define AMR_YT_WRITER_H

#include "Amr/AbstractAmrWriter.h"

#include <boost/filesystem.hpp>

#include <vector>

class AmrYtWriter : public AbstractAmrWriter {
    
public:
    
    /*!
     * @param step we write
     * @param bin energy bin we write (multi-bunch simulation)
     */
    AmrYtWriter(int step, int bin = 0);
    
    /*!
     * Write yt files to the simulation subdirectory
     * data/amr/yt.
     * The data can be visualized using the python script
     * pyOPALTools/amrPlots/visualize.py. Use the help
     * to find out how to call the script.
     */
    void writeFields(const amr::AmrScalarFieldContainer_t& rho,
                     const amr::AmrScalarFieldContainer_t& phi,
                     const amr::AmrVectorFieldContainer_t& efield,
                     const amr::AmrIntArray_t& refRatio,
                     const amr::AmrGeomContainer_t& geom,
                     const int& nLevel,
                     const double& time,
                     const double& scale);
    
    
    void writeBunch(const AmrPartBunch* bunch_p,
                    const double& time,
                    const double& gamma);
    
private:
    /* Copied and slightely modified version of
     * AMReX_ParticleContainerI.H
     * 
     * @param level to write
     * @param ofs out stream
     * @param fnum file number
     * @param which file
     * @param count how many particles on this grid
     * @param where file offset
     * @param bunch_p to get data from
     * @param gamma is the Lorentz factor
     */
    void writeParticles_m(int level,
                          std::ofstream& ofs,
                          int fnum,
                          amrex::Vector<int>& which,
                          amrex::Vector<int>& count,
                          amrex::Vector<long>& where,
                          const AmrPartBunch* bunch_p,
                          const double gamma) const;
    
private:
    std::string dir_m;                      ///< directory where to write files
    
    std::vector<std::string> intData_m;     ///< integer bunch data
    std::vector<std::string> realData_m;    ///< real bunch data
};

#endif
