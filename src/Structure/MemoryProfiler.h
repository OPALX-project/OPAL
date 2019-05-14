#ifndef OPAL_MEMORY_PROFILER_H
#define OPAL_MEMORY_PROFILER_H

#include <fstream>
#include <string>
#include <map>
#include <vector>


class MemoryProfiler {
    /* Pay attention with units. /proc/[pid]/status returns values in
     * KiB (Kibibyte) although the units say kB.
     * KiB has base 2 not base 10
     */
    
public:
    typedef std::vector<long double> vm_t;
    typedef std::vector<std::string> units_t;
    
    MemoryProfiler();
    
    ~MemoryProfiler();

    enum VirtualMemory {
        VMPEAK = 0, // VmPeak: Peak virtual memory size.
        VMSIZE,     // VmSize: Virtual memory size.
        VMHWM,      // VmHWM: Peak resident set size ("high water mark").
        VMRSS,      // VmRSS: Resident set size.
        VMDATA,     // VmData: Size of data.
        VMSTK,      // VmStk: Size of stack.
        VMEXE,      // VmExe: Size of text segments.
        VMLCK,      // VmLck: Locked memory size (see mlock(3)).
        VMPIN,      // VmPin: Pinned memory size (since Linux 3.2).  These are pages that can't be moved because something
                    // needs to directly access physical memory.
        VMLIB,      // VmLib: Shared library code size.
        VMPTE,      // VmPTE: Page table entries size (since Linux 2.6.10).
        VMPMD,      // VmPMD: Size of second-level page tables (since Linux 4.0).
        VMSWAP      // VmSwap: Swapped-out virtual memory size by  anonymous  private  pages;  shmem  swap  usage  is  not
                    // included (since Linux 2.6.34).
    };
    
    void write(double time);
    
private:
    void header_m();
    void update_m();
    void compute_m(vm_t& vmMin, vm_t& vmMax, vm_t& vmAvg);
    
private:
    std::string fname_m;
    std::ofstream ofstream_m;
    std::map<std::string, int> procinfo_m;
    vm_t vmem_m;
    units_t unit_m;
};

#endif
