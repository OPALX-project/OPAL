#include "MemoryProfiler.h"

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#endif

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "Utilities/OpalException.h"

#include <boost/filesystem.hpp>

#include "Ippl.h"

#include <iomanip>


MemoryProfiler::MemoryProfiler()
    : fname_m(OpalData::getInstance()->getInputBasename() + ".mem")
{
    procinfo_m = {
        {"VmPeak:", VirtualMemory::VMPEAK},
        {"VmSize:", VirtualMemory::VMSIZE},
        {"VmHWM:",  VirtualMemory::VMHWM},
        {"VmRSS:",  VirtualMemory::VMRSS},
        {"VmStk:",  VirtualMemory::VMSTK},
        {"VmData:", VirtualMemory::VMDATA},
        {"VmExe:",  VirtualMemory::VMEXE},
        {"VmLck:",  VirtualMemory::VMLCK},
        {"VmPin:",  VirtualMemory::VMPIN},
        {"VmLib:",  VirtualMemory::VMLIB},
        {"VmPTE:",  VirtualMemory::VMPTE},
        {"VmPMD:",  VirtualMemory::VMPMD},
        {"VmSwap:", VirtualMemory::VMSWAP}
    };
    
    vmem_m.resize(procinfo_m.size());
    unit_m.resize(procinfo_m.size());
    
    if ( !boost::filesystem::exists(fname_m) ) {
        this->header_m();
    }
}


MemoryProfiler::~MemoryProfiler() {
    if ( ofstream_m.is_open() )
        ofstream_m.close();
}


void MemoryProfiler::header_m() {
    
    if ( Ippl::myNode() != 0 )
        return;
    
    ofstream_m.open(fname_m, std::ios::out);
    
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");

    ofstream_m << "SDDS1" << std::endl;
    ofstream_m << "&description\n"
               << indent << "text=\"Memory statistics '"
               << OpalData::getInstance()->getInputFn() << "' "
               << dateStr << "" << timeStr << "\",\n"
               << indent << "contents=\"memory info\"\n"
               << "&end\n";
    ofstream_m << "&parameter\n"
               << indent << "name=processors,\n"
               << indent << "type=long,\n"
               << indent << "description=\"Number of Cores used\"\n"
               << "&end\n";
    ofstream_m << "&parameter\n"
               << indent << "name=revision,\n"
               << indent << "type=string,\n"
               << indent << "description=\"git revision of opal\"\n"
               << "&end\n";
    ofstream_m << "&parameter\n"
               << indent << "name=flavor,\n"
               << indent << "type=string,\n"
               << indent << "description=\"OPAL flavor that wrote file\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=t,\n"
               << indent << "type=double,\n"
               << indent << "units=ns,\n"
               << indent << "description=\"1 Time\"\n"
               << "&end\n";

    // peak virtual memory size
    ofstream_m << "&column\n"
               << indent << "name=VmPeak-Min,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMPEAK] + ",\n"
               << indent << "description=\"2 Minimum peak virtual memory size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmPeak-Max,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMPEAK] + ",\n"
               << indent << "description=\"3 Maximum peak virtual memory size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmPeak-Avg,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMPEAK] + ",\n"
               << indent << "description=\"4 Average peak virtual memory size\"\n"
               << "&end\n";

    // virtual memory size
    ofstream_m << "&column\n"
               << indent << "name=VmSize-Min,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSIZE] + ",\n"
               << indent << "description=\"5 Minimum virtual memory size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmSize-Max,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSIZE] + ",\n"
               << indent << "description=\"6 Maximum virtual memory size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmSize-Avg,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSIZE] + ",\n"
               << indent << "description=\"7 Average virtual memory size\"\n"
               << "&end\n";

    // peak resident set size ("high water mark")
    ofstream_m << "&column\n"
               << indent << "name=VmHWM-Min,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMHWM] + ",\n"
               << indent << "description=\"8 Minimum peak resident set size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmHWM-Max,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMHWM] + ",\n"
               << indent << "description=\"9 Maximum peak resident set size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmHWM-Avg,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMHWM] + ",\n"
               << indent << "description=\"10 Average peak resident set size\"\n"
               << "&end\n";

    // resident set size
    ofstream_m << "&column\n"
               << indent << "name=VmRSS-Min,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMRSS] + ",\n"
               << indent << "description=\"11 Minimum resident set size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmRSS-Max,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMRSS] + ",\n"
               << indent << "description=\"12 Maximum resident set size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmRSS-Avg,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMRSS] + ",\n"
               << indent << "description=\"13 Average resident set size\"\n"
               << "&end\n";

    // stack size
    ofstream_m << "&column\n"
               << indent << "name=VmStk-Min,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSTK] + ",\n"
               << indent << "description=\"14 Minimum stack size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmStk-Max,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSTK] + ",\n"
               << indent << "description=\"15 Maximum stack size\"\n"
               << "&end\n";
    ofstream_m << "&column\n"
               << indent << "name=VmStk-Avg,\n"
               << indent << "type=double,\n"
               << indent << "units=" + unit_m[VirtualMemory::VMSTK] + ",\n"
               << indent << "description=\"16 Average stack size\"\n"
               << "&end\n";
    
    ofstream_m << "&data\n"
               << indent << "mode=ascii,\n"
               << indent << "no_row_counts=1\n"
               << "&end\n";

    ofstream_m << Ippl::getNodes() << std::endl;
    ofstream_m << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl;
    ofstream_m << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                   (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
    ofstream_m.close();
}


void MemoryProfiler::update_m() {
#ifdef __linux__
    static pid_t pid = getpid();
    std::string fname = "/proc/" + std::to_string(pid) + "/status";
    
    if ( !boost::filesystem::exists(fname) ) {
        throw OpalException("MemoryProfiler::update_m()",
                            "File '" + fname + "' doesn't exist.");
    }
    
    std::ifstream ifs(fname.c_str());
    
    if ( !ifs.is_open() ) {
        throw OpalException("MemoryProfiler::update_m()",
                            "Failed to open '" + fname + "'.");
    }
    
    std::string token = "";
    while (ifs >> token) {
        if ( !procinfo_m.count(token) ) {
            continue;
        }
        int idx = procinfo_m[token];
        ifs >> vmem_m[idx]
            >> unit_m[idx];
    }
    
    ifs.close();
#endif
}


void MemoryProfiler::compute_m(vm_t& vmMin,
                               vm_t& vmMax,
                               vm_t& vmAvg)
{
    if ( Ippl::getNodes() == 1 ) {
        for (unsigned int i = 0; i < vmem_m.size(); ++i) {
            vmMin[i] = vmMax[i] = vmAvg[i] = vmem_m[i];
        }
        return;
    }
    
    new_reduce(vmem_m.data(), vmAvg.data(), vmem_m.size(), std::plus<double>());
    
    double inodes = 1.0 / double(Ippl::getNodes());
    for (auto& vm : vmAvg) {
        vm *= inodes;
    }

    new_reduce(vmem_m.data(), vmMin.data(), vmem_m.size(), std::less<double>());
    new_reduce(vmem_m.data(), vmMax.data(), vmem_m.size(), std::greater<double>());
}


void MemoryProfiler::write(double time) {
    
    this->update_m();
    
    vm_t vmMin(vmem_m.size());
    vm_t vmMax(vmem_m.size());
    vm_t vmAvg(vmem_m.size());
    
    this->compute_m(vmMin, vmMax, vmAvg);
    
    if ( Ippl::myNode() != 0 ) {
        return;
    }
    
    int width = 12;
    
    ofstream_m.open(fname_m, std::ios::app);
    
    ofstream_m << time          << std::setw(width) << "\t"
               << vmMin[VMPEAK] << std::setw(width) << "\t"
               << vmMax[VMPEAK] << std::setw(width) << "\t"
               << vmAvg[VMPEAK] << std::setw(width) << "\t"
               
               << vmMin[VMSIZE] << std::setw(width) << "\t"
               << vmMax[VMSIZE] << std::setw(width) << "\t"
               << vmAvg[VMSIZE] << std::setw(width) << "\t"
               
               << vmMin[VMHWM]  << std::setw(width) << "\t"
               << vmMax[VMHWM]  << std::setw(width) << "\t"
               << vmAvg[VMHWM]  << std::setw(width) << "\t"
               
               << vmMin[VMRSS]  << std::setw(width) << "\t"
               << vmMax[VMRSS]  << std::setw(width) << "\t"
               << vmAvg[VMRSS]  << std::setw(width) << "\t"
               
               << vmMin[VMSTK]  << std::setw(width) << "\t"
               << vmMax[VMSTK]  << std::setw(width) << "\t"
               << vmAvg[VMSTK]  << std::setw(width) << "\n";

    ofstream_m.close();
}
