#include "H5Writer.h"

H5Writer::H5Writer(H5PartWrapper* h5wrapper, bool restart)
    : H5PartTimer_m(IpplTimings::getTimer("Write H5-File"))
    , h5wrapper_m(h5wrapper)
    , H5call_m(0)
{
    if ( !restart ) {
        h5wrapper->writeHeader();
    }
    h5wrapper->close();
}


void H5Writer::writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[]) {
    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("B-ref_x", FDext[0](0)),
        std::make_pair("B-ref_z", FDext[0](1)),
        std::make_pair("B-ref_y", FDext[0](2)),
        std::make_pair("E-ref_x", FDext[1](0)),
        std::make_pair("E-ref_z", FDext[1](1)),
        std::make_pair("E-ref_y", FDext[1](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);
}


int H5Writer::writePhaseSpace(PartBunchBase<double, 3> *beam, Vector_t FDext[], double meanEnergy,
                              double refPr, double refPt, double refPz,
                              double refR, double refTheta, double refZ,
                              double azimuth, double elevation, bool local) {

    if (beam->getTotalNum() < 3) return -1; // in single particle mode and tune calculation (2 particles) we do not need h5 data

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", refPr),
        std::make_pair("REFPT", refPt),
        std::make_pair("REFPZ", refPz),
        std::make_pair("REFR", refR),
        std::make_pair("REFTHETA", refTheta),
        std::make_pair("REFZ", refZ),
        std::make_pair("AZIMUTH", azimuth),
        std::make_pair("ELEVATION", elevation),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x",  FDext[2](0)),
        std::make_pair("B-ref_z",  FDext[2](1)),
        std::make_pair("B-ref_y",  FDext[2](2)),
        std::make_pair("E-ref_x",  FDext[3](0)),
        std::make_pair("E-ref_z",  FDext[3](1)),
        std::make_pair("E-ref_y",  FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    ++ H5call_m;
    return H5call_m - 1;
}


void H5Writer::writePhaseSpace(EnvelopeBunch &beam, Vector_t FDext[],
                               double sposHead, double sposRef,
                               double sposTail)
{
    // FIXME https://gitlab.psi.ch/OPAL/src/issues/245

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("sposHead", sposHead),
        std::make_pair("sposRef",  sposRef),
        std::make_pair("sposTail", sposTail),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x",  FDext[2](0)),
        std::make_pair("B-ref_z",  FDext[2](1)),
        std::make_pair("B-ref_y",  FDext[2](2)),
        std::make_pair("E-ref_x",  FDext[3](0)),
        std::make_pair("E-ref_z",  FDext[3](1)),
        std::make_pair("E-ref_y",  FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    // does not compile with that line --> should be fixed with https://gitlab.psi.ch/OPAL/src/issues/245
//     h5wrapper_m->writeStep(&beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);
}


//FIXME https://gitlab.psi.ch/OPAL/src/issues/245
// void H5Writer::stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {
// 
//     if (!doHDF5_m) return;
// 
//     /// Start timer.
//     IpplTimings::startTimer(H5PartTimer_m);
// 
//     static_cast<H5PartWrapperForPS*>(h5wrapper_m)->stashPhaseSpaceEnvelope(beam,
//                                                                            FDext,
//                                                                            sposHead,
//                                                                            sposRef,
//                                                                            sposTail);
//     H5call_m++;
// 
//     /// %Stop timer.
//     IpplTimings::stopTimer(H5PartTimer_m);
// }
// 
// void H5Writer::dumpStashedPhaseSpaceEnvelope() {
// 
//     if (!doHDF5_m) return;
// 
//     static_cast<H5PartWrapperForPS*>(h5wrapper_m)->dumpStashedPhaseSpaceEnvelope();
// 
//     /// %Stop timer.
//     IpplTimings::stopTimer(H5PartTimer_m);
// }