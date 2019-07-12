#ifndef OPAL_H5_WRITER_H
#define OPAL_H5_WRITER_H

#include "Structure/H5PartWrapper.h"
#include "Algorithms/PartBunchBase.h"
#include "H5hut.h"

class H5Writer {
    
public:
    
    H5Writer(H5PartWrapper* h5wrapper, bool restart);

    void close();
    
    void changeH5Wrapper(H5PartWrapper *h5wrapper);
    
    
    void storeCavityInformation();
    
    double getLastPosition();
    
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
    int writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[], double E,
                        double refPr, double refPt, double refPz,
                        double refR, double refTheta, double refZ,
                        double azimuth, double elevation, bool local);

    /**
     * FIXME https://gitlab.psi.ch/OPAL/src/issues/245
     * \brief Dumps Phase Space for Envelope trakcer to H5 file.
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
    void writePhaseSpace(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
//     void stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail);
//     void dumpStashedPhaseSpaceEnvelope();
    
private:
    /// Timer to track particle data/H5 file write time.
    IpplTimings::TimerRef H5PartTimer_m;
    
    H5PartWrapper *h5wrapper_m;
    
    /// Current record, or time step, of H5 file.
    int H5call_m;
};


inline
void H5Writer::close() {
    h5wrapper_m->close();
}


inline
void H5Writer::changeH5Wrapper(H5PartWrapper *h5wrapper) {
    h5wrapper_m = h5wrapper;
}


inline
void H5Writer::storeCavityInformation() {
    h5wrapper_m->storeCavityInformation();
}


inline
double H5Writer::getLastPosition() {
    return h5wrapper_m->getLastPosition();
}

#endif
