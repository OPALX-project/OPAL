#include "config.h"

#include "Structure/SDDSParser.h"

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>

void printInfo(const std::string &input);
void printUsage(char **argv);

int main(int argc, char **argv) {
    std::string inputFile(""), outputFile("/dev/stdout");

    if (argc == 1) {
        printUsage(argv);
        return 0;
    }
    if (argc == 2 && std::string(argv[1]).substr(0, 2) != "--") {
        inputFile = std::string(argv[1]);

    } else {
        for (int i = 1; i < argc; ++ i) {
            if (std::string(argv[i]) == "--input" && i + 1 < argc) {
                inputFile = std::string(argv[++ i]);
            } else if (std::string(argv[i]) == "--output" && i + 1 < argc) {
                outputFile = std::string(argv[++ i]);
            } else if (std::string(argv[i]) == "--help") {
                printUsage(argv);
                return 0;
            } else {
                inputFile = std::string(argv[i]);
            }
        }
    }

    if (inputFile == "") {
        std::cerr << "Error: no input file provided!\n" << std::endl;
        printUsage(argv);
        return 1;
    }

    SDDS::SDDSParser parser;
    try {
        parser.setInput(inputFile);
        parser.run();
    } catch (OpalException ex) {
        std::cerr << ex.where() << "\n"
                  << ex.what() << std::endl;
        return 1;
    }

    double mins = 0.0, maxs = 0.0;
    try {
        auto sdata = parser.getColumnData("s");
        mins = boost::get<double>(sdata.front());
        maxs = boost::get<double>(sdata.back());
    } catch (OpalException ex) {
        std::cerr << ex.where() << "\n"
                  << ex.what() << std::endl;
        return 1;
    }

    std::cout << "data range: [" << mins << " m, " << maxs << " m]" << std::endl;
    return 0;
}


void printUsage(char **argv) {
    std::string name(argv[0]), indent("  ");
    name = name.substr(name.find_last_of('/') + 1);
    unsigned int width = 25;
    std::cout << name << " version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n"
              << "Usage\n\n"
              << indent << name << " [options] <path to input>\n\n"
              << "Options\n\n"
              << indent << std::setw(width) << std::left << "--input <path to input>" << "= path to input file\n"
              << indent << std::setw(width) << std::left << "--output <file name>"    << "= name of output. If omitted, directed\n"
              << indent << std::setw(width + 2)          << " "                       <<   "to stdout\n"
              << indent << std::setw(width) << std::left << "--step <step number>"    << "= step of the H5Hut file that should be\n"
              << indent << std::setw(width + 2) <<          " "                       <<   "exported. Default: last step\n"
              << indent << std::setw(width) << std::left << "--binary"                << "= whether SDDS output should be in \n"
              << indent << std::setw(width + 2) << std::left << " "                   <<   "binary format. Default: ASCII\n"
              << indent << std::setw(width) << std::left << "--info"                  << "= printing path length corresponding to step\n"
              << indent << std::setw(width) << std::left << "--help"                  << "= this help\n"
              << std::endl;
}