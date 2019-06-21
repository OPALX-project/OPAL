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

#include <fstream>

#include "Algorithms/PBunchDefs.h"
#include "Utilities/Util.h"
#include "H5hut.h"

#ifdef __linux__
    #include "MemoryProfiler.h"
#endif

template <class T, unsigned Dim>
class PartBunchBase;
class EnvelopeBunch;
class BoundaryGeometry;
class H5PartWrapper;

class DataSink {
public:

    typedef std::vector<std::pair<std::string, unsigned int> > losses_t;

    /** \brief Default constructor.
     *
     * The default constructor is called at the start of a new calculation (as
     * opposed to a calculation restart).
     */
    DataSink();
    DataSink(H5PartWrapper *h5wrapper);

    /** \brief Restart constructor.
     *
     * This constructor is called when a calculation is restarted using data from
     * an existing H5 file.
     */
    DataSink(H5PartWrapper *h5wrapper, int restartStep);

    ~DataSink();

    void reset();

    void rewindLinesLBal(size_t numberOfLines) const;
    unsigned int rewindSDDStoSPos(double maxSpos) const;

    bool doHDF5();

    void changeH5Wrapper(H5PartWrapper *h5wrapper);

    /** \brief Write cavity information from  H5 file
     *
     *
     *
     */
    void storeCavityInformation();

    void setMaxNumBunches(int nBunches);

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
     */
    void doWriteStatData(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                         double E, const losses_t &losses, const double& azimuth = -1);

    /** \brief for OPAL-t

     */
    void writeStatData(PartBunchBase<double, 3> *beam, Vector_t FDext[],
                       const losses_t &losses = losses_t(), const double& azimuth = -1);

    // /** \brief for OPAL-cycl

    //  */
    void writeStatData(PartBunchBase<double, 3> *beam, Vector_t FDext[], double E, const double& azimuth = -1);


    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information to beam statistics file so the SDDS tools can be used
     * for plotting data.
     * \param outputFile Name of file to write to.
     *
     */
    void writeSDDSHeader(std::ofstream &outputFile);

    void writeSDDSHeader(std::ofstream &outputFile, const losses_t &losses);

    /** \brief Dumps Phase Space to H5 file.
     *
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     */
    void writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[]);

    /** \brief Dumps phase space to H5 file in OPAL cyclotron calculation.
     *
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     *  \param E average energy (MeB)
     *  \return Returns the number of the time step just written.
     */
    int writePhaseSpace_cycl(PartBunchBase<double, 3> *beam, Vector_t FDext[], double E,
                             double refPr, double refPt, double refPz,
                             double refR, double refTheta, double refZ,
                             double azimuth, double elevation, bool local);

    /** \brief Dumps Phase Space for Envelope trakcer to H5 file.
     *
     * \param beam The beam.
     * \param FDext The external E and B field for the head, reference and tail particles. The vector array
     * has the following layout:
     *  - FDext[0] = B at head particle location (in x, y and z).
     *  - FDext[1] = E at head particle location (in x, y and z).
     *  - FDext[2] = B at reference particle location (in x, y and z).
     *  - FDext[3] = E at reference particle location (in x, y and z).
     *  - FDext[4] = B at tail particle location (in x, y, and z).
     *  - FDext[5] = E at tail particle location (in x, y, and z).
     *  \param sposHead Longitudinal position of the head particle.
     *  \param sposRef Longitudinal position of the reference particle.
     *  \param sposTail Longitudinal position of the tail particles.
     */
    void writePhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
    void stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
    void dumpStashedPhaseSpaceEnvelope();
    void writeStatData(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);

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
    void writeImpactStatistics(PartBunchBase<double, 3> *beam, long long int &step, size_t &impact, double &sey_num,
                               size_t numberOfFieldEmittedParticles, bool nEmissionMode, std::string fn);

    void writeSurfaceInteraction(PartBunchBase<double, 3> *beam, long long int &step, BoundaryGeometry &bg, std::string fn);

    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information to processor statistics file so the SDDS tools can be used
     * for plotting data.
     * \param outputFile Name of file to write to.
     *
     */
    void writeLBalHeader(PartBunchBase<double, 3> *beam, std::ofstream &outputFile);

    void writeLBalData(PartBunchBase<double, 3> *beam,
                       std::ofstream &outputFile,
                       unsigned int pwi);


    /** \brief Write SDDS header.
     *
     * Writes the appropriate SDDS format header information to processor memory so the SDDS tools can be used
     * for plotting data.
     * \param outputFile Name of file to write to.
     *
     */
    void writeMemoryHeader(std::ofstream &outputFile);

    void writeMemoryData(PartBunchBase<double, 3> *beam,
                         std::ofstream &outputFile,
                         unsigned int pwi);

#ifdef ENABLE_AMR
    /** \brief Write SDDS header. (AMR only)
     *
     * Writes the appropriate SDDS format header information to grid load balancing so the SDDS tools can be used
     * for plotting data.
     * \param outputFile Name of file to write to.
     *
     */
    void writeGridLBalHeader(PartBunchBase<double, 3> *beam,
                             std::ofstream &outputFile);

    void writeGridLBalData(PartBunchBase<double, 3> *beam,
                           std::ofstream &outputFile,
                           unsigned int pwi);

    bool writeAmrStatistics(PartBunchBase<double, 3> *beam);
    
    void noAmrDump(PartBunchBase<double, 3> *beam);
#endif


private:

    DataSink(const DataSink &) { }
    DataSink &operator = (const DataSink &) { return *this; }

    static std::string convertToString(int number);

    void rewindLines(const std::string &fileName, size_t numberOfLines) const;
    void replaceVersionString(const std::string &fileName) const;

    void open_m(std::ofstream& os, const std::string& fileName) const;

    /** \brief First write to the statistics output file.
     *
     * Initially set to std::ios::out so that SDDS format header information is written to file
     * during the first write call to the statistics output file. Variable is then
     * reset to std::ios::app so that header information is only written once.
     */
    std::ios_base::openmode mode_m;

    int nMaxBunches_m;

    /** \brief First write to the H5 surface loss file.
     *
     * If true, file name will be assigned and file will be prepared to write.
     * Variable is then reset to false so that H5 file is only initialized once.
     */
    bool firstWriteH5Surface_m;

    /// Name of output file for beam statistics.
    std::string statFileName_m;

    /// Name of output file for processor load balancing information.
    std::string lBalFileName_m;

    /// Name of output file for processor memory information
    std::string memFileName_m;

#ifdef ENABLE_AMR
    /// Name of output file for grid load balancing information
    std::string gridLBalFileName_m;
#endif

    /// Name of output file for surface loss data.
    std::string surfaceLossFileName_m;

    /// H5 file for surface loss data.
    h5_file_t H5fileS_m;

    /// Current record, or time step, of H5 file.
    int H5call_m;

    /// Timer to track statistics write time.
    IpplTimings::TimerRef StatMarkerTimer_m;

    /// Timer to track particle data/H5 file write time.
    IpplTimings::TimerRef H5PartTimer_m;

    /// needed to create index for vtk file
    unsigned int lossWrCounter_m;

    /// flag to discable all HDF5 output
    bool doHDF5_m;

    H5PartWrapper *h5wrapper_m;

#ifdef __linux__
    std::unique_ptr<MemoryProfiler> memprof_m;
#endif
};

inline
void DataSink::reset() {
    H5call_m = 0;
}

/** \brief
 *   delete the last 'numberOfLines' lines of the load balance file
 */
inline
void DataSink::rewindLinesLBal(size_t numberOfLines) const {
    if (Ippl::myNode() == 0) {
        rewindLines(lBalFileName_m, numberOfLines);
    }
}

/** \brief
 *  delete the last 'numberOfLines' lines of the statistics file
 */
inline
unsigned int DataSink::rewindSDDStoSPos(double maxSPos) const {
    if (Ippl::myNode() == 0) {
        return Util::rewindLinesSDDS(statFileName_m, maxSPos);
    }

    return 0;
}

inline
void DataSink::changeH5Wrapper(H5PartWrapper *h5wrapper) {
    h5wrapper_m = h5wrapper;
}

inline
std::string DataSink::convertToString(int number) {
    std::stringstream ss;
    ss << std::setw(5) << std::setfill('0') <<  number;
    return ss.str();
}

#endif // DataSink_H_

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: