#include "opal.h"

extern Ippl *ippl;
extern Inform *gmsg;


#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Utilities/OpalException.h"
#include "Fields/Fieldmap.h"
#include "Structure/IpplInfoWrapper.h"

#include "OPALconfig.h"

#ifdef ENABLE_AMR
    #include <AMReX.H>
#endif


int run_opal(char */*args*/[], std::string inputfile, int restartStep,
             int infoLevel, int warnLevel, MPI_Comm comm)
{
    std::string::size_type startExtension    = inputfile.find_last_of('.');
    std::string outputFileName = inputfile.substr(0,startExtension) + ".out";
    std::ofstream output(outputFileName.c_str());

    MPI_Barrier(comm);

    IpplInfoWrapper *newippl = new IpplInfoWrapper(inputfile, infoLevel, warnLevel, comm);
    gmsg = new Inform("OPAL ", output);
    IpplInfo::Info->setDestination(output);
    IpplInfo::Error->setDestination(output);
    IpplInfo::Warn->setDestination(output);

#ifdef ENABLE_AMR
    amrex::Initialize(comm);
#endif

    OpalData *opal = OpalData::getInstance();

    Configure::configure();
    opal->storeInputFn(inputfile);

    // FileStream is a RCObject
    FileStream *is = 0;
    try {
        is = new FileStream(inputfile);
    } catch(...) {
        is = 0;
        throw new OpalException("run_opal", "Could not open inputfile: " + inputfile);
    }

    // run simulation
    OpalParser *parser = new OpalParser();

    if (restartStep > std::numeric_limits<int>::min()) {
        opal->setRestartRun();
        opal->setRestartStep(restartStep);
        opal->setRestartFileName(inputfile.substr(0,startExtension) + ".h5");
    }

    if(is) parser->run(is);

    Ippl::Comm->barrier();

    IpplInfo::Info->setDestination(std::cout);
    IpplInfo::Error->setDestination(std::cout);
    IpplInfo::Warn->setDestination(std::cout);

    // cleanup
    //OPAL->reset();
    OpalData::deleteInstance();
    Fieldmap::clearDictionary();
    delete parser;
    delete gmsg;

#ifdef ENABLE_AMR
    amrex::Finalize(true);
#endif

    //FIXME: strange side effects
    //ippl = 0;
    //delete aippl;

    //XXX: seems like Ippl is always returning the same instance after the
    //     initial instantiation.
    delete newippl;

    output.close();
    return 0;
}
