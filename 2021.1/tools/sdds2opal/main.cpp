#include "mdb.h"
#include "SDDS.h"

#include "config.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <unistd.h>

std::vector<std::vector<double> > readSDDSFile(std::string fname, bool temporal);
void printUsage(char **argv);

int main(int argc, char **argv) {
    std::string inputFile(""), outputFile("/dev/stdout");
    bool temporal = false;

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
            } else if (std::string(argv[i]) == "--temporal") {
                temporal = true;
            } else if (std::string(argv[i]) == "--help") {
                printUsage(argv);
                return 0;
            }
        }
    }

    if (inputFile == "") {
        if (!isatty(fileno(stdin))) {
            inputFile = "/dev/stdin";
        } else {
            std::cerr << "Error: no input file provided!\n" << std::endl;
            printUsage(argv);
            return 1;
        }
    }

    auto data = readSDDSFile(inputFile, temporal);
    std::ofstream out(outputFile);
    out.precision(8);
    out << data.size() << "\n";
    for (auto row: data) {
        for (auto column: row)
            out << std::setw(18) << column;
        out << "\n";
    }
    out.close();

    return 0;
}

std::vector<std::vector<double> > readSDDSFile(std::string fname, bool temporal) {

    SDDS_DATASET SDDS_dataset;
    void *columnData;
    long page, rows;
    std::vector<std::vector<double> > fileData;

    char buffer[256];
    strcpy(buffer, fname.c_str());
    if (SDDS_InitializeInput(&SDDS_dataset, buffer ) != 1) {
        std::cerr << "Error: couldn't initialize SDDS file\n" << std::endl;
        std::exit(1);
    }

    char *columnNames[6];
    for (unsigned int i = 0; i < 6u; ++ i) {
        columnNames[i] = new char[3];
    }
    strcpy(columnNames[0], "x");
    strcpy(columnNames[1], "xp");
    strcpy(columnNames[2], "y");
    strcpy(columnNames[3], "yp");
    strcpy(columnNames[4], "t");
    strcpy(columnNames[5], "p");

    long ids[6];
    for (unsigned int i = 0; i < 6; ++ i) {
        ids[i] = SDDS_GetColumnIndex(&SDDS_dataset, columnNames[i]);

        if (ids[i] < 0) {
            std::cerr << "Error: couldn't get column index of column '" << columnNames[i] << "'\n" << std::endl;
            std::exit(1);
        }
    }

    SDDS_SetColumnsOfInterest(&SDDS_dataset, SDDS_NAME_ARRAY, 6, columnNames);
    SDDS_DeleteUnsetColumns(&SDDS_dataset);

    page=SDDS_ReadPage(&SDDS_dataset);
    while (page >= 1) {
        rows = SDDS_RowCount(&SDDS_dataset);
        std::vector<std::vector<double> > pageData(rows);

        for (unsigned int i = 0; i < 6; ++ i) {

            columnData = SDDS_GetColumn(&SDDS_dataset, columnNames[i]);

            for (unsigned int j = 0; j < rows; ++ j) {
                pageData[j].push_back(((double*)columnData)[j]);
            }

            if (columnData) free(columnData);
        }
        fileData.insert(fileData.end(), pageData.begin(), pageData.end());

        page=SDDS_ReadPage(&SDDS_dataset);
    }

    for (std::vector<double> &row: fileData) {
        row[1] *= row[5];
        row[3] *= row[5];
        if (temporal) continue;

        row[4] *= -299792458.0;
    }

    if (SDDS_Terminate(&SDDS_dataset)!=1) {
        std::cerr << "Error: couldn't terminate SDDS properly\n" << std::endl;
        std::exit(1);
    }

    for (unsigned int i = 0; i < 6u; ++ i) {
        delete[] columnNames[i];
    }

    return fileData;
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
              << indent << std::setw(width) << std::left << "--time-of-arrival"       << "= extract and save the time of arrival instead of the longitudinal position\n"
              << indent << std::setw(width + 2)          << " "                       <<   "to stdout\n"
              << indent << std::setw(width) << std::left << "--help"                  << "= this help\n"
              << std::endl;
}