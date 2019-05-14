#include "AmrYtWriter.h"

#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>
#include <AMReX_VectorIO.H>
#include <AMReX_NFiles.H>

#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"


AmrYtWriter::AmrYtWriter(int step, int bin)
{
    intData_m.resize(1);
//     intData_m[0] = "id";
//     intData_m[1] = "cpu";
    intData_m[0] = "energy_bin";
    
    realData_m.resize(//3 +  // coordinates
                      3 +  // momenta
                      3 +  // charge + mass + timestep
//                       1 +  // potential at particle location
                      3 +  // electric field at particle location
                      3);  // magnetic field at particle location
    
    int idx = 0;
//     realData_m[idx++] ="position_x";
//     realData_m[idx++] ="position_y";
//     realData_m[idx++] ="position_z";
    realData_m[idx++] ="momentum_x";
    realData_m[idx++] ="momentum_y";
    realData_m[idx++] ="momentum_z";
    realData_m[idx++] ="charge";
    realData_m[idx++] ="mass";
    realData_m[idx++] ="timestep";
//     realData_m[idx++] ="potential";
    realData_m[idx++] ="electric_field_x";
    realData_m[idx++] ="electric_field_y";
    realData_m[idx++] ="electric_field_z";
    realData_m[idx++] ="magnetic_field_x";
    realData_m[idx++] ="magnetic_field_y";
    realData_m[idx++] ="magnetic_field_z";
    
    namespace fs = boost::filesystem;
    
    fs::path dir = OpalData::getInstance()->getInputBasename();
    boost::filesystem::path path = dir.parent_path() / "data" / "amr" / "yt";
    dir_m = amrex::Concatenate((path / "plt").string(), step, 10);
    dir_m += "-";
    dir_m = amrex::Concatenate(dir_m, bin, 3);
    
    if ( Ippl::myNode() == 0 && !fs::exists(path) ) {
        try {
            fs::create_directories( path );
        } catch(const fs::filesystem_error& ex) {
            throw OpalException("AmrYtWriter::AmrYtWriter()",
                                ex.what());
        }
    }
    
    Ippl::Comm->barrier();
}


void AmrYtWriter::writeFields(const amr::AmrScalarFieldContainer_t& rho,
                              const amr::AmrScalarFieldContainer_t& phi,
                              const amr::AmrVectorFieldContainer_t& efield,
                              const amr::AmrIntArray_t& refRatio,
                              const amr::AmrGeomContainer_t& geom,
                              const int& nLevel,
                              const double& time,
                              const double& scale)
{
    /* We need to scale the geometry and cell sizes according to the
     * particle mapping
     */
    
    //
    // Only let 64 CPUs be writing at any one time.
    //
    amrex::VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if ( Ippl::myNode() == 0 )
        if (!amrex::UtilCreateDirectory(dir_m, 0755))
            amrex::CreateDirectoryFailed(dir_m);
    //
    // Force other processors to wait till directory is built.
    //
    Ippl::Comm->barrier();

    std::string HeaderFileName = dir_m + "/Header";

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int nData = rho[0]->nComp()
              + phi[0]->nComp()
              + efield[0][0]->nComp()
              + efield[0][1]->nComp()
              + efield[0][2]->nComp();
    
    if ( Ippl::myNode() == 0 )
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);
        HeaderFile << "HyperCLaw-V1.1\n";

        HeaderFile << nData << '\n';

        // variable names
        for (int ivar = 1; ivar <= rho[0]->nComp(); ivar++)
            HeaderFile << "rho\n";
        
        for (int ivar = 1; ivar <= phi[0]->nComp(); ivar++)
            HeaderFile << "phi\n";
        
        HeaderFile << "Ex\nEy\nEz\n";
        
        // dimensionality
        HeaderFile << AMREX_SPACEDIM << '\n';
        
        // time
        HeaderFile << time << '\n';
        HeaderFile << nLevel - 1 << '\n'; // maximum level number (0=single level)
        
        // physical domain
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbLo(i) / scale << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            HeaderFile << geom[0].ProbHi(i) / scale << ' ';
        HeaderFile << '\n';
        
        // reference ratio
        for (int i = 0; i < refRatio.size(); ++i)
            HeaderFile << refRatio[i] << ' ';
        HeaderFile << '\n';
        
        // geometry domain for all levels
        for (int i = 0; i < nLevel; ++i)
            HeaderFile << geom[i].Domain() << ' ';
        HeaderFile << '\n';
        
        // number of time steps
        for (int i = 0; i < nLevel; ++i)
            HeaderFile << 0 << ' ';
        HeaderFile << '\n';
        
        // cell sizes for all level
        for (int i = 0; i < nLevel; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                HeaderFile << geom[i].CellSize()[k] / scale << ' ';
            HeaderFile << '\n';
        }
        
        // coordinate system
        HeaderFile << geom[0].Coord() << '\n';
        HeaderFile << "0\n"; // write boundary data
    }
    
    for (int lev = 0; lev < nLevel; ++lev) {
        // Build the directory to hold the MultiFab at this level.
        // The name is relative to the directory containing the Header file.
        //
        static const std::string BaseName = "/Cell";
        char buf[64];
        sprintf(buf, "Level_%d", lev);
        std::string sLevel = buf;
    
        //
        // Now for the full pathname of that directory.
        //
        std::string FullPath = dir_m;
        if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
            FullPath += '/';
        FullPath += sLevel;
        //
        // Only the I/O processor makes the directory if it doesn't already exist.
        //
        if ( Ippl::myNode() == 0 )
            if (!amrex::UtilCreateDirectory(FullPath, 0755))
                amrex::CreateDirectoryFailed(FullPath);
        //
        // Force other processors to wait till directory is built.
        //
        Ippl::Comm->barrier();
        
        if ( Ippl::myNode() == 0 )
        {
            HeaderFile << lev << ' ' << rho[lev]->boxArray().size() << ' ' << 0 /*time*/ << '\n';
            HeaderFile << 0 /* level steps */ << '\n';
    
            for (int i = 0; i < rho[lev]->boxArray().size(); ++i)
            {
                amrex::Real dx[3] = {
                    geom[lev].CellSize(0),
                    geom[lev].CellSize(1),
                    geom[lev].CellSize(2)
                };
                dx[0] /= scale;
                dx[1] /= scale;
                dx[2] /= scale;
                amrex::Real lo[3] = {
                    geom[lev].ProbLo(0),
                    geom[lev].ProbLo(1),
                    geom[lev].ProbLo(2)
                };
                lo[0] /= scale;
                lo[1] /= scale;
                lo[2] /= scale;
                amr::AmrDomain_t loc = amr::AmrDomain_t(rho[lev]->boxArray()[i],
                                                        &dx[0], &lo[0]);
//                                                         geom[lev].CellSize(),
//                                                         geom[lev].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; n++)
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
    
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            HeaderFile << PathNameInHeader << '\n';
        }
        
        //
        // We combine all of the multifabs 
        //
        amr::AmrField_t data(rho[lev]->boxArray(),
                             rho[lev]->DistributionMap(),
                             nData, 0,
                             amrex::MFInfo());
        
        //
        // Cull data -- use no ghost cells.
        //
        // dst, src, srccomp, dstcomp, numcomp, nghost
        /*
        * srccomp: the component to copy
        * dstcmop: the component where to copy
        * numcomp: how many components to copy
        */
        amr::AmrField_t::Copy(data, *rho[lev],       0, 0, 1, 0);
        amr::AmrField_t::Copy(data, *phi[lev],       0, 1, 1, 0);
        amr::AmrField_t::Copy(data, *efield[lev][0], 0, 2, 1, 0); // Ex
        amr::AmrField_t::Copy(data, *efield[lev][1], 0, 3, 1, 0); // Ey
        amr::AmrField_t::Copy(data, *efield[lev][2], 0, 4, 1, 0); // Ez
        
        //
        // Use the Full pathname when naming the MultiFab.
        //
        std::string TheFullPath = FullPath;
        TheFullPath += BaseName;
    
        amrex::VisMF::Write(data, TheFullPath, amrex::VisMF::NFiles, true);
    }
    
    if ( Ippl::myNode() == 0 )
        HeaderFile.close();
}


void AmrYtWriter::writeBunch(const AmrPartBunch* bunch_p,
                             const double& time,
                             const double& gamma)
{
    /* According to
     * 
     * template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
     * void
     * ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
     * ::Checkpoint (const std::string&        dir,
     *               const std::string&        name,
     *               bool                      is_checkpoint,
     *               const Vector<std::string>& real_comp_names,
     *               const Vector<std::string>& int_comp_names) const
     * 
     * in AMReX_ParticleContainerI.H with AMReX_Particles.H and AMReX_ParticleI.H.
     * 
     * ATTENTION: We need to Lorentz transform the particles.
     */
    
    const AmrLayout_t* layout_p = static_cast<const AmrLayout_t*>(&bunch_p->getLayout());
    const AmrPartBunch::pbase_t* amrpbase_p = bunch_p->getAmrParticleBase();
    
    const int  NProcs       = amrex::ParallelDescriptor::NProcs();
    const int  IOProcNumber = amrex::ParallelDescriptor::IOProcessorNumber();
    
    bool doUnlink = true;
    
    //
    // We store the particles in a subdirectory of "dir".
    //
    std::string pdir = dir_m;
    
    if ( ! pdir.empty() && pdir[pdir.size()-1] != '/') {
        pdir += '/';
    }
    
    pdir += "opal";
    //
    // Make the particle directories if they don't already exist.
    //
    if ( Ippl::myNode() == 0 ) {
        if ( ! amrex::UtilCreateDirectory(pdir, 0755)) {
            amrex::CreateDirectoryFailed(pdir);
        }
    }
    
    // Force other processors to wait until directory is built.
    Ippl::Comm->barrier();
    
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the I/O processor writes to the header file.
    //
    std::ofstream HdrFile;
    
    long nParticles = bunch_p->getTotalNum();
    
    // do not modify LocalNumPerLevel in here!!!
    auto LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();
    
    int nLevel = (&amrpbase_p->getAmrLayout())->maxLevel() + 1;
    
    std::unique_ptr<size_t[]> partPerLevel( new size_t[nLevel] );
    std::unique_ptr<size_t[]> globalPartPerLevel( new size_t[nLevel] );
    
    for (size_t i = 0; i < LocalNumPerLevel.size(); ++i)
        partPerLevel[i] = LocalNumPerLevel[i];
    
    allreduce(*partPerLevel.get(),
              *globalPartPerLevel.get(),
              nLevel, std::plus<size_t>());
    
    int finest_level = layout_p->finestLevel();

    if ( Ippl::myNode() == 0 )
    {
        std::string HdrFileName = pdir;
        
        if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/') {
            HdrFileName += '/';
            
        }
        
        HdrFileName += "Header";
        
        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);
        
        if ( ! HdrFile.good()) {
            amrex::FileOpenFailed(HdrFileName);
        }
        //
        // First thing written is our Checkpoint/Restart version string.
        // 
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        HdrFile << "Version_Two_Dot_Zero_double" << '\n';
        //
        // AMREX_SPACEDIM and N for sanity checking.
        //
        HdrFile << AMREX_SPACEDIM << '\n';
        
        // The number of extra real parameters
        HdrFile << realData_m.size() << '\n';
        for (std::size_t j = 0; j < realData_m.size(); ++j)
            HdrFile << realData_m[j] << '\n';
        
        // The number of extra int parameters
        HdrFile << intData_m.size() << '\n';
        for (std::size_t j = 0; j < intData_m.size(); ++j)
            HdrFile << intData_m[j] << '\n';
        
        // is_checkpoint
        HdrFile << true << '\n';
        
        //
        // The total number of particles.
        //
        HdrFile << nParticles << '\n';
        //
        // The value of nextid that we need to restore on restart.
        // --> we don's use this feature
        //
        HdrFile << 0 << '\n';
        //
        // Then the finest level of the AMR hierarchy.
        //
        HdrFile << finest_level << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= finest_level; ++lev) {
            HdrFile << layout_p->ParticleBoxArray(lev).size() << '\n';
        }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles(256);

    nOutFiles = std::max(1, std::min(nOutFiles, NProcs));
    
    for (int lev = 0; lev <= finest_level; ++lev) {
        bool gotsome = (globalPartPerLevel[lev] > 0);
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;
        
        if (gotsome) {
            if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') {
                LevelDir += '/';
            }
            
            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
            
            if ( Ippl::myNode() == 0 ) {
                if ( ! amrex::UtilCreateDirectory(LevelDir, 0755)) {
                    amrex::CreateDirectoryFailed(LevelDir);
                }
            }
            //
            // Force other processors to wait until directory is built.
            //
            Ippl::Comm->barrier();
        }
        
        // Write out the header for each particle
        if ( gotsome && Ippl::myNode() == 0 ) {
            std::string HeaderFileName = LevelDir;
            HeaderFileName += "/Particle_H";
            std::ofstream ParticleHeader(HeaderFileName);
            
            layout_p->ParticleBoxArray(lev).writeOn(ParticleHeader);
            ParticleHeader << '\n';

            ParticleHeader.flush();
            ParticleHeader.close();
        }
        
        amrex::MFInfo info;
        info.SetAlloc(false);
        amr::AmrField_t state(layout_p->ParticleBoxArray(lev),
                              layout_p->ParticleDistributionMap(lev),
                              1, 0, info);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        amrex::Vector<int>  which(state.size(),0);
        amrex::Vector<int > count(state.size(),0);
        amrex::Vector<long> where(state.size(),0);
        
        std::string filePrefix(LevelDir);
        filePrefix += '/';
        filePrefix += "DATA_";
        bool groupSets(false), setBuf(true);

        if (gotsome) {
            for(amrex::NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi) {
                std::ofstream& myStream = (std::ofstream&) nfi.Stream();
                //
                // Write out all the valid particles we own at the specified level.
                // Do it grid block by grid block remembering the seek offset
                // for the start of writing of each block of data.
                //
                writeParticles_m(lev, myStream, nfi.FileNumber(), which, count, where, bunch_p, gamma);
            }

            amrex::ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
            amrex::ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
            amrex::ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
        }

        if ( Ippl::myNode() == 0 ) {
            for (int j = 0; j < state.size(); j++) {
                //
                // We now write the which file, the particle count, and the
                // file offset into which the data for each grid was written,
                // to the header file.
                //
                HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
            }

            if (gotsome && doUnlink) {
                //
                // Unlink any zero-length data files.
                //
                amrex::Vector<long> cnt(nOutFiles,0);

                for (int i = 0, N=count.size(); i < N; i++) {
                    cnt[which[i]] += count[i];
                }

                for (int i = 0, N=cnt.size(); i < N; i++) {
                    if (cnt[i] == 0) {
                        std::string FullFileName = amrex::NFilesIter::FileName(i, filePrefix);
                        amrex::UnlinkFile(FullFileName.c_str());
                    }
                }
            }
        }
    }            // ---- end for(lev...)

    if ( Ippl::myNode() == 0 ) {
        HdrFile.flush();
        HdrFile.close();
        if ( ! HdrFile.good()) {
            throw OpalException("AmrYtWriter:writeBunch()",
                                "Problem writing HdrFile");
        }
    }
}


void AmrYtWriter::writeParticles_m(int level,
                                   std::ofstream& ofs,
                                   int fnum,
                                   amrex::Vector<int>& which,
                                   amrex::Vector<int>& count,
                                   amrex::Vector<long>& where,
                                   const AmrPartBunch* bunch_p,
                                   const double gamma) const
{
    /* Copied and modified of
     * 
     * template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
     * void
     * ParticleContainer<NStructReal,
     *                   NStructInt,
     *                   NArrayReal,
     *                   NArrayInt>
     * ::WriteParticles(int              lev,
     *                  std::ofstream&   ofs,
     *                  int              fnum,
     *                  Vector<int>&     which,
     *                  Vector<int>&     count,
     *                  Vector<long>&    where,
     *                  bool             is_checkpoint) const
     * 
     * in AMReX_ParticleContainerI.H
     */
    
    const AmrLayout_t* layout_p = static_cast<const AmrLayout_t*>(&bunch_p->getLayout());
    const AmrPartBunch::pbase_t* amrpbase_p = bunch_p->getAmrParticleBase();
    
    const auto& LocalNumPerLevel = amrpbase_p->getLocalNumPerLevel();
    size_t lBegin = LocalNumPerLevel.begin(level);
    size_t lEnd   = LocalNumPerLevel.end(level);
    
    for (size_t ip = lBegin; ip < lEnd; ++ip) {
        const int grid = amrpbase_p->Grid[ip];
        
        count[grid] += 1;
    }
    
    amrex::MFInfo info;
    info.SetAlloc(false);
    amr::AmrField_t state(layout_p->ParticleBoxArray(level),
                          layout_p->ParticleDistributionMap(level),
                          1, 0, info);

    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        
        which[grid] = fnum;
        where[grid] = amrex::VisMF::FileOffset(ofs);
      
        if (count[grid] == 0) {
            continue;
        }
        
        
        // First write out the integer data in binary.
        const int iChunkSize = 2 + intData_m.size();
        amrex::Vector<int> istuff(count[grid]*iChunkSize);
        int* iptr = istuff.dataPtr();
        
        // inefficient
        for (size_t ip = lBegin; ip < lEnd; ++ip) {
            const int pGrid = amrpbase_p->Grid[ip];
            
            if ( grid == pGrid ) {
                
                iptr[0] = bunch_p->ID[ip];
                iptr[1] = Ippl::myNode();
                iptr[2] = bunch_p->Bin[ip];
                iptr += iChunkSize;
            }
        }

        amrex::writeIntData(istuff.dataPtr(), istuff.size(), ofs);
      
        // Write the Real data in binary.
        const int rChunkSize = AMREX_SPACEDIM + realData_m.size();
        amrex::Vector<double> rstuff(count[grid]*rChunkSize);
        double* rptr = rstuff.dataPtr();
        
        
        // inefficient
        for (size_t ip = lBegin; ip < lEnd; ++ip) {
            const int pGrid = amrpbase_p->Grid[ip];
            
            if ( grid == pGrid ) {
                
                int idx = 0;
                
                // coordinates
                for (int j = 0; j < AMREX_SPACEDIM; ++j)
                    rptr[idx++] = bunch_p->R[ip](j);
                
                // Lorentz transform
                if ( OpalData::getInstance()->isInOPALCyclMode() )
                    rptr[1] *= gamma; // y is longitudinal direction
                else
                    rptr[2] *= gamma; // z is longitudinal direction
                
                // momenta
                for (int j = 0; j < AMREX_SPACEDIM; ++j)
                    rptr[idx++] = bunch_p->P[ip](j);
                
                // charge
                rptr[idx++] = bunch_p->Q[ip];
                
                // mass
                rptr[idx++] = bunch_p->M[ip];
                
                // timestep
                rptr[idx++] = bunch_p->dt[ip];
                
//                 // potential
//                 rptr[idx++] = bunch_p->Phi[ip];
                
                // electric field
                for (int j = 0; j < AMREX_SPACEDIM; ++j)
                    rptr[idx++] = bunch_p->Ef[ip](j);
                
                // magnetic field
                for (int j = 0; j < AMREX_SPACEDIM; ++j)
                    rptr[idx++] = bunch_p->Bf[ip](j);
                
                rptr += idx; //realData_m.size();
            }
        }
        amrex::writeDoubleData(rstuff.dataPtr(), rstuff.size(), ofs);
    }
}
