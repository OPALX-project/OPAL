// ------------------------------------------------------------------------
// $RCSfile: Main.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.9.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Main program for OPAL
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:10 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "opal.h"

Ippl *ippl;
Inform *gmsg;

#include "H5hut.h"

#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/TerminalStream.h"
#include "Utilities/Timer.h"
#include "Fields/Fieldmap.h"
#include "FixedAlgebra/FTps.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"
#include "Utilities/EarlyLeaveException.h"
#include "Utilities/Util.h"

#include "OPALconfig.h"

#ifdef ENABLE_AMR
#include <AMReX_ParallelDescriptor.H>
#endif
/*
  Includes related to the optimizer
*/
#include "boost/smart_ptr.hpp"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"

#include "Optimizer/EA/FixedPisaNsga2.h"
#include "Optimizer/EA/BlendCrossover.h"
#include "Optimizer/EA/IndependentBitMutation.h"

#include "Optimize/OpalSimulation.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"
#include "Expression/SumErrSq.h"
#include "Expression/SDDSVariable.h"
#include "Expression/RadialPeak.h"
#include "Expression/SumErrSqRadialPeak.h"
#include "Expression/ProbeVariable.h"

#include <gsl/gsl_errno.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <cstring>
#include <set>
#include <algorithm>

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int line,
                         int gsl_errno) {
        throw OpalException(file, reason);
    }
}


bool checkInitAmrFlag(int argc, char* argv[]) {
    std::string noamr = "noInitAMR";
    bool initAMR = true;
    for (int i = 0; i < argc; ++i) {
        std::string sargv = std::string(argv[i]);
        if ( sargv.find(noamr) != std::string::npos ) {
            initAMR = false;
            break;
        }
    }
    return initAMR;
}


int main(int argc, char *argv[]) {
    Ippl *ippl = new Ippl(argc, argv);
    gmsg = new  Inform("OPAL");

    namespace fs = boost::filesystem;

#ifdef ENABLE_AMR
    bool initAMR = checkInitAmrFlag(argc, argv);
    if ( initAMR ) {
        // false: build no parmparse, we use the OPAL parser instead.
        amrex::Initialize(argc, argv, false, Ippl::getComm());
    }
#endif

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    H5SetVerbosityLevel(1); //65535);

    gsl_set_error_handler(&errorHandlerGSL);

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
    IpplTimings::startTimer(mainTimer);

    Inform hmsg("");
    std::string mySpace("            ");

    if(Ippl::myNode() == 0) remove("errormsg.txt");

    hmsg << mySpace <<  "   ____  _____       ___ " << endl;
    hmsg << mySpace <<  "  / __ \\|  __ \\ /\\   | | " << endl;
    hmsg << mySpace <<  " | |  | | |__) /  \\  | |" << endl;
    hmsg << mySpace <<  " | |  | |  ___/ /\\ \\ | |" << endl ;
    hmsg << mySpace <<  " | |__| | |  / ____ \\| |____" << endl;
    hmsg << mySpace <<  "  \\____/|_| /_/    \\_\\______|" << endl;


    std::string gitRevision = "git rev. " + Util::getGitRevision();
    std::string copyRight = "(c) PSI, http://amas.web.psi.ch";
    *gmsg << endl
          << "This is OPAL (Object Oriented Parallel Accelerator Library) Version " << OPAL_PROJECT_VERSION << "\n"
          << std::setw(37 + gitRevision.length() / 2) << std::right << gitRevision << "\n\n" << endl
          << std::setw(37 + copyRight.length() / 2) << std::right << copyRight << "\n\n" << endl
          << "The optimiser (former opt-Pilot) is integrated " << endl
          << endl;

#ifdef OPAL_DKS
    *gmsg << "OPAL compiled with DKS (Dynamic Kernel Scheduler) Version "
          << DKS_VERSION;
    if (IpplInfo::DKSEnabled)
      *gmsg << " GPU present" << endl << endl;
    else
      *gmsg << " GPU not present" << endl << endl;
#endif

    *gmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) \nto the OPAL developers " << PACKAGE_BUGREPORT << "\n" << endl;
    *gmsg << "Time: " << timeStr << " date: " << dateStr << "\n" << endl;


    /*
      Make a directory data for some of the output
    */
    if(Ippl::myNode() == 0) {
        if (!fs::exists("data")) {
            boost::system::error_code error_code;
            if (!fs::create_directory("data", error_code)) {
                std::cerr << error_code.message() << std::endl;
                // use error code to prevent create_directory from throwing an exception
            }
        }
    }
    Ippl::Comm->barrier();
    if (!fs::is_directory("data")) {
        std::cerr << "unable to create directory; aborting" << std::endl;
        abort();
    }

    const OpalParser parser;

    //  DTA
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cerr.precision(16);
    std::cerr.setf(std::ios::scientific, std::ios::floatfield);
    // /DTA

    // Set global truncation orders.
    FTps<double, 2>::setGlobalTruncOrder(20);
    FTps<double, 4>::setGlobalTruncOrder(15);
    FTps<double, 6>::setGlobalTruncOrder(10);

    OpalData *opal = OpalData::getInstance();
    opal->storeArguments(argc, argv);
    try {
        Configure::configure();

        // Read startup file.
        FileStream::setEcho(Options::echo);

        char *startup = getenv("HOME");
        boost::filesystem::path p = strncat(startup, "/init.opal", 20);
        if (startup != NULL && is_regular_file(p)) {

            FileStream::setEcho(false);
            FileStream *is;

            try {
                is = new FileStream(startup);
            } catch(...) {
                is = 0;
                ERRORMSG("Could not open startup file \"" << startup << "\".\n"
                         << "Note: this is not mandatory for an OPAL simulation!\n");
            }

            if(is) {
                *gmsg << "Reading startup file \"" << startup << "\"." << endl;
                parser.run(is);
                *gmsg << "Finished reading startup file." << endl;
            }
            FileStream::setEcho(Options::echo);
        } else {
            *gmsg << "Couldn't find startup file \"" << startup << "\".\n"
                  << "Note: this is not mandatory for an OPAL simulation!\n" << endl;
        }

        if(argc <= 1) {
            // Run commands from standard input
            parser.run(new TerminalStream("OPAL"));
        } else {
            int arg = -1;
            std::string fname;
            std::string restartFileName;
            //         // will write dumping date into a new h5 file
            for(int ii = 1; ii < argc; ++ ii) {
                std::string argStr = std::string(argv[ii]);
                // The sequence of the two arguments is free
                if (argStr == std::string("--input")) {
                    ++ ii;
                    arg = ii;
                    INFOMSG(argv[ii] << endl);
                    continue;
                } else if (argStr == std::string("-restart") ||
                           argStr == std::string("--restart")) {
                    opal->setRestartRun();
                    opal->setRestartStep(atoi(argv[++ ii]));
                    opal->setRestartFileName(argv[1]);
                    continue;
                } else if (argStr == std::string("-restartfn") ||
                           argStr == std::string("--restartfn")) {
                    restartFileName = std::string(argv[++ ii]);
                    continue;
                } else if (argStr == std::string("-version") ||
                           argStr == std::string("--version")) {
                    INFOMSG("OPAL Version " << OPAL_PROJECT_VERSION << ", git rev. " << Util::getGitRevision() << endl);
                    IpplInfo::printVersion(true);
                    std::string options = (IpplInfo::compileOptions() +
                                           std::string(" ") +
                                           std::string(OPAL_COMPILE_OPTIONS) +
                                           std::string(" "));
                    std::set<std::string> uniqOptions;
                    while (options.length() > 0) {
                        size_t n = options.find_first_of(' ');
                        while (n == 0) {
                            options = options.substr(n + 1);
                            n = options.find_first_of(' ');
                        }

                        uniqOptions.insert(options.substr(0, n));
                        options = options.substr(n + 1);
                    }
                    for (auto it: uniqOptions) {
                        options += it + " ";
                    }

                    std::string header("Compile-time options: ");
                    while (options.length() > 58) {
                        std::string line = options.substr(0, 58);
                        size_t n = line.find_last_of(' ');
                        INFOMSG(header << line.substr(0, n) << "\n");

                        header = std::string(22, ' ');
                        options = options.substr(n + 1);
                    }
                    INFOMSG(header << options << endl);
                    exit(0);
                } else if ( argStr.find("noInitAMR") != std::string::npos) {
                    // do nothing here
                } else if (argStr == std::string("-help") ||
                           argStr == std::string("--help")) {
                    IpplInfo::printHelp(argv);
                    INFOMSG("   --version            : Print a brief version summary.\n");
                    INFOMSG("   --input <fname>      : Specifies the input file <fname>.\n");
                    INFOMSG("   --restart <n>        : Performes a restart from step <n>.\n");
                    INFOMSG("   --restartfn <fname>  : Uses the file <fname> to restart from.\n");
                    INFOMSG("   --help               : Display this command-line summary.\n");
                    INFOMSG(endl);
                    exit(0);
                } else {
                    if (arg == -1 &&
                        (ii == 1 || ii + 1 == argc) &&
                        argv[ii][0] != '-') {
                        arg = ii;
                        continue;
                    } else {
                        INFOMSG("Unknown argument \"" << argStr << "\"" << endl);
                        IpplInfo::printHelp(argv);
                        INFOMSG("   --version            : Print a brief version summary.\n");
                        INFOMSG("   --input <fname>      : Specifies the input file <fname>.\n");
                        INFOMSG("   --restart <n>        : Performes a restart from step <n>.\n");
                        INFOMSG("   --restartfn <fname>  : Uses the file <fname> to restart from.\n");
                        INFOMSG("   --help               : Display this command-line summary.\n");
                        INFOMSG(endl);
                        exit(1);
                    }
                }
            }

            if (arg == -1) {
                INFOMSG("No input file provided!" << endl);
                exit(1);
            }

            fname = std::string(argv[arg]);
            if (!fs::exists(fname)) {
                INFOMSG("Input file \"" << fname << "\" doesn't exist!" << endl);
                exit(1);
            }

            opal->storeInputFn(fname);

            if (opal->inRestartRun()) {
                if (restartFileName == "")
                    restartFileName = opal->getInputBasename() + std::string(".h5");
                if (!fs::exists(restartFileName)) {
                    INFOMSG("Restart file \"" << restartFileName << "\" doesn't exist!" << endl);
                    exit(1);
                }
                opal->setRestartFileName(restartFileName);
            }

            FileStream *is;

            try {
                is = new FileStream(fname);
            } catch(...) {
                is = 0;
                *gmsg << "Input file \"" << fname << "\" not found." << endl;
            }

            if(is) {
                *gmsg << "* Reading input stream \"" << fname << "\"." << endl;
                parser.run(is);
                *gmsg << "* End of input stream \"" << fname << "\"." << endl;
            }
        }

        if(Ippl::myNode() == 0) {
            std::ifstream errormsg("errormsg.txt");
            if(errormsg.good()) {
                char buffer[256];
                std::string closure("                                                                                 *\n");
                ERRORMSG("\n"
                         << "* **********************************************************************************\n"
                         << "* ************** W A R N I N G / E R R O R * * M E S S A G E S *********************\n"
                         << "* **********************************************************************************"
                         << endl);
                errormsg.getline(buffer, 256);
                while(errormsg.good()) {
                    ERRORMSG("* ");
                    if(errormsg.gcount() == 1) {
                        ERRORMSG(closure);
                    } else if ((size_t)errormsg.gcount() <= closure.size()) {
                        ERRORMSG(buffer << closure.substr(errormsg.gcount() - 1));
                    } else {
                        ERRORMSG(buffer << endl);
                    }
                    errormsg.getline(buffer, 256);
                }
                ERRORMSG("* " << closure
                         << "* **********************************************************************************\n"
                         << "* **********************************************************************************"
                         << endl);
            }
            errormsg.close();
        }

    } catch(EarlyLeaveException& ex) {
        // do nothing here
    } catch(OpalException &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n*** User error detected by function \""
                 << ex.where() << "\"\n";
        // stat->printWhere(errorMsg, true);
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(ClassicException &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n*** User error detected by function \""
                 << ex.where() << "\"\n";
        // stat->printWhere(errorMsg, true);
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(SDDSParserException &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);

        std::stringstream msg;
        errorMsg << "\n*** Error detected by function \""
                 << ex.where() << "\"\n";
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(IpplException &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);

        std::stringstream msg;
        errorMsg << "\n*** Error detected by function \""
                 << ex.where() << "\"\n";
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(std::bad_alloc &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n*** Error:\n";
        errorMsg << "    Sorry, virtual memory exhausted.\n"
                 << ex.what()
                 << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(assertion &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n*** Runtime-error ******************\n";
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        errorMsg << "\n************************************\n" << endl;
        throw std::runtime_error("in Parser");
    } catch(std::exception &ex) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n"
                 << "*** Error:\n"
                 << "    Internal OPAL error: \n";
        std::string what = ex.what();
        size_t pos = what.find_first_of('\n');
        do {
            errorMsg << "    " << what.substr(0, pos) << endl;
            what = what.substr(pos + 1, std::string::npos);
            pos = what.find_first_of('\n');
        } while (pos != std::string::npos);
        errorMsg << "    " << what << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    } catch(...) {
        Inform errorMsg("Error", std::cerr, INFORM_ALL_NODES);
        errorMsg << "\n*** Error:\n"
                 << "    Unexpected exception caught.\n" << endl;

        MPI_Abort(MPI_COMM_WORLD, -100);
    }


    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    IpplTimings::print(std::string("timing.dat"),
                       OpalData::getInstance()->getProblemCharacteristicValues());

    Ippl::Comm->barrier();
    Fieldmap::clearDictionary();
    OpalData::deleteInstance();
    delete gmsg;

#ifdef ENABLE_AMR
    if ( initAMR ) {
        amrex::Finalize(true);
    }
#endif

    delete ippl;
    delete Ippl::Info;
    delete Ippl::Warn;
    delete Ippl::Error;
    delete Ippl::Debug;

    return 0;
}