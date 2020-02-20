//
//  Copyright & License: See Copyright.readme in src directory
//

/**
   \brief Class: DataSink

   This class acts as an observer during the calculation. It generates diagnostic
   output of the accelerated beam such as statistical beam descriptors of particle
   positions, momenta, beam phase space (emittance) etc. These are written to file
   at periodic time steps during the calculation.

   This class also writes the full beam phase space to an H5 file at periodic time
   steps in the calculation (this period is different from that of the statistical
   numbers).

   Class also writes processor load balancing data to file to track parallel
   calculation efficiency.
*/

#ifndef _OPAL_DATA_SINK_H
#define _OPAL_DATA_SINK_H

#include "Algorithms/PBunchDefs.h"

#include "StatWriter.h"
#include "H5Writer.h"
#include "MultiBunchDump.h"

template <class T, unsigned Dim>
class PartBunchBase;
class BoundaryGeometry;
class H5PartWrapper;

class DataSink {
public:
    typedef StatWriter::losses_t            losses_t;
    typedef std::unique_ptr<StatWriter>     statWriter_t;
    typedef std::unique_ptr<SDDSWriter>     sddsWriter_t;
    typedef std::unique_ptr<H5Writer>       h5Writer_t;
    typedef std::unique_ptr<MultiBunchDump> mbWriter_t;
    
    
    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();
    DataSink(H5PartWrapper *h5wrapper, bool restart, short numBunch);
    DataSink(H5PartWrapper *h5wrapper, short numBunch);

    void dumpH5(PartBunchBase<double, 3> *beam, Vector_t FDext[]) const;
    
    int dumpH5(PartBunchBase<double, 3> *beam, Vector_t FDext[], double meanEnergy,
               double refPr, double refPt, double refPz,
               double refR, double refTheta, double refZ,
               double azimuth, double elevation, bool local) const;
    
    //FIXME https://gitlab.psi.ch/OPAL/src/issues/245
    void dumpH5(EnvelopeBunch &beam, Vector_t FDext[],
                double sposHead, double sposRef,
                double sposTail) const;
    
    
    void dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                  const double& azimuth = -1) const;
    
    void dumpSDDS(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                  const losses_t &losses = losses_t(), const double& azimuth = -1) const;
    
    //FIXME https://gitlab.psi.ch/OPAL/src/issues/245
    void dumpSDDS(EnvelopeBunch &beam, Vector_t FDext[],
                  double sposHead, double sposRef, double sposTail) const;

    /** \brief Write cavity information from  H5 file
     */
    void storeCavityInformation();
    
    void changeH5Wrapper(H5PartWrapper *h5wrapper);
    
    
    /**
     * Write particle loss data to an ASCII fille for histogram
     * @param fn specifies the name of ASCII file
     * @param beam
     */
    void writePartlossZASCII(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg, std::string fn);
    
    /**
     * Write geometry points and surface triangles to vtk file
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeGeomToVtk(BoundaryGeometry &bg, std::string fn);
    //void writeGeoContourToVtk(BoundaryGeometry &bg, std::string fn);
    
    
    /**
     * Write impact number and outgoing secondaries in each time step
     *
     * @param fn specifies the name of vtk file contains the geometry
     *
     */
    void writeImpactStatistics(PartBunchBase<double, 3> *beam,
                               long long int &step,
                               size_t &impact,
                               double &sey_num,
                               size_t numberOfFieldEmittedParticles,
                               bool nEmissionMode,
                               std::string fn);

    /** no statWriter_m dump
     * @param beam
     * @param mbhandler is the multi-bunch handler
     */
    void writeMultiBunchStatistics(PartBunchBase<double, 3> *beam,
                                   MultiBunchHandler* mbhandler);

    /**
     * In restart mode we need to set the correct path length
     * of each bunch
     * @param mbhandler is the multi-bunch handler
     */
    void setMultiBunchInitialPathLengh(MultiBunchHandler* mbhandler_p);

private:
    DataSink(const DataSink& ds) = delete;
    DataSink &operator = (const DataSink &) = delete;

    void rewindLines();

    void init(bool restart = false,
              H5PartWrapper* h5wrapper = nullptr,
              short numBunch = 1);


    h5Writer_t      h5Writer_m;
    statWriter_t    statWriter_m;
    std::vector<sddsWriter_t> sddsWriter_m;
    std::vector<mbWriter_t> mbWriter_m;

    static std::string convertToString(int number, int setw = 5);

    /// needed to create index for vtk file
    unsigned int lossWrCounter_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;

    const bool isMultiBunch_m;

    void initMultiBunchDump(short numBunch);
};


inline
std::string DataSink::convertToString(int number, int setw) {
    std::stringstream ss;
    ss << std::setw(setw) << std::setfill('0') <<  number;
    return ss.str();
}


#endif // DataSink_H_

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
