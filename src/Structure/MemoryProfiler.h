//
// Class MemoryProfiler
//   This class writes a SDDS file with virtual memory usage information (Linux only).
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_MEMORY_PROFILER_H
#define OPAL_MEMORY_PROFILER_H

#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "SDDSWriter.h"

class MemoryProfiler : public SDDSWriter {
    /* Pay attention with units. /proc/[pid]/status returns values in
     * KiB (Kibibyte) although the units say kB.
     * KiB has base 2 not base 10
     */
    
public:
    typedef std::vector<long double> vm_t;
    typedef std::vector<std::string> units_t;
    
    MemoryProfiler(const std::string& fname, bool restart);

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
    
    void write(PartBunchBase<double, 3> *beam) override;
    
private:
    void header();
    void update();
    void compute(vm_t& vmMin, vm_t& vmMax, vm_t& vmAvg);
    
private:
    std::string fname_m;
    std::map<std::string, int> procinfo_m;
    vm_t vmem_m;
    units_t unit_m;
};

#endif
