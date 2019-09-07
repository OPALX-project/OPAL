#include "config.h"
#include "mdb.h"
#include "SDDS.h"
#include "H5hut.h"

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>

struct Attribute {
    std::string name;
    std::string unit;
    std::vector<double> value;

    bool operator< (const Attribute& at) const {
        return this->name < at.name;
    }
};

typedef std::map<std::string, std::vector<double> > data_t;
typedef std::set<Attribute> attributes_t;

#define REPORTONERROR(rc) reportOnError(rc, __FILE__, __LINE__)
#define READDATA(type, file, name, value) REPORTONERROR(H5PartReadData##type(file, name, value));
#define READFILEATTRIB(type, file, name, value) REPORTONERROR(H5ReadFileAttrib##type(file, name, value));
#define READSTEPATTRIB(type, file, name, value) REPORTONERROR(H5ReadStepAttrib##type(file, name, value));

enum FORMAT {
    BINARY = SDDS_BINARY,
    ASCII = SDDS_ASCII
};

typedef h5_file_t file_t;

data_t readStepData(file_t file, bool temporal);
attributes_t readStepAttributes(file_t file);
void readH5HutFile(const std::string &fname,
                   size_t step,
                   bool temporal,
                   data_t &data,
                   attributes_t &attr);
void convertToElegantUnits(data_t &data);
void writeSDDSFile(const std::string &fname, const data_t &data, const attributes_t &attr, FORMAT form);
void printInfo(const std::string &input);
void printUsage(char **argv);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    std::string inputFile(""), outputFile("/dev/stdout");
    bool doPrintInfo = false;
    bool temporal = false;
    int step = -1;
    FORMAT form = ASCII;

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
            } else if (std::string(argv[i]) == "--step" && i + 1 < argc) {
                step = atoi(argv[++ i]);
            } else if (std::string(argv[i]) == "--binary") {
                form = BINARY;
            } else if (std::string(argv[i]) == "--help") {
                printUsage(argv);
                return 0;
            } else if (std::string(argv[i]) == "--info") {
                doPrintInfo = true;
            } else if (std::string(argv[i]) == "--temporal") {
                temporal = true;
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

    if (doPrintInfo) {
        printInfo(inputFile);
        return 1;
    }

    data_t data;
    attributes_t attr;
    readH5HutFile(inputFile, step, temporal, data, attr);
    convertToElegantUnits(data);
    writeSDDSFile(outputFile, data, attr, form);

    return 0;
}

void reportOnError(h5_int64_t rc, const char* file, int line) {
    if (rc != H5_SUCCESS)
        std::cerr << "H5 rc= " << rc << " in " << file << " @ line " << line << std::endl;
}

void readH5HutFile(const std::string &fname,
                   size_t step,
                   bool temporal,
                   data_t &data,
                   attributes_t &attr) {
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = MPI_COMM_WORLD;
    H5SetPropFileMPIOCollective (props, &comm);

    file_t file = H5OpenFile(fname.c_str(), H5_O_RDONLY, props);
    H5CloseProp (props);
    h5_ssize_t numStepsInSource = H5GetNumSteps(file);
    h5_ssize_t readStep = (step > (size_t)(numStepsInSource - 1)? numStepsInSource - 1: step);

    REPORTONERROR(H5SetStep(file, readStep));

    data = readStepData(file, temporal);
    attr = readStepAttributes(file);
    REPORTONERROR(H5CloseFile(file));
}

data_t readStepData(file_t file, bool temporal) {
    data_t data;

    h5_ssize_t numParticles = H5PartGetNumParticles(file);

    std::vector<h5_float64_t> buffer(numParticles);
    h5_float64_t *f64buffer = &buffer[0];
    data_t::mapped_type dData(numParticles);

    READDATA(Float64, file, "x", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        dData[n] = f64buffer[n];
    }
    data.insert(std::make_pair(std::string("x"), dData));

    READDATA(Float64, file, "y", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        dData[n] = f64buffer[n];
    }
    data.insert(std::make_pair(std::string("y"), dData));

    if (temporal) {
        READDATA(Float64, file, "time", f64buffer);
        for(long int n = 0; n < numParticles; ++ n) {
            dData[n] = f64buffer[n];
        }
        data.insert(std::make_pair(std::string("time"), dData));
    } else {
        READDATA(Float64, file, "z", f64buffer);
        for(long int n = 0; n < numParticles; ++ n) {
            dData[n] = f64buffer[n];
        }
        data.insert(std::make_pair(std::string("z"), dData));
    }

    READDATA(Float64, file, "px", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        dData[n] = f64buffer[n];
    }
    data.insert(std::make_pair(std::string("px"), dData));

    READDATA(Float64, file, "py", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        dData[n] = f64buffer[n];
    }
    data.insert(std::make_pair(std::string("py"), dData));

    READDATA(Float64, file, "pz", f64buffer);
    for(long int n = 0; n < numParticles; ++ n) {
        dData[n] = f64buffer[n];
    }
    data.insert(std::make_pair(std::string("pz"), dData));

    return data;
}
Attribute readAttribute(file_t file,
                        const std::string &name,
                        const std::string &H5name,
                        unsigned int numComponents) {
    h5_float64_t f64buffer[3];
    char cbuffer[128];
    std::string H5unit = H5name + "Unit";

    Attribute attr;
    attr.name = name;
    READFILEATTRIB(String, file, H5unit.c_str(), cbuffer);
    attr.unit = std::string(cbuffer);
    READSTEPATTRIB(Float64, file, H5name.c_str(), f64buffer);
    for (unsigned int d = 0; d < numComponents; ++ d)
        attr.value.push_back(f64buffer[d]);

    return attr;
}

attributes_t readStepAttributes(file_t file) {
    attributes_t attr;

    attr.insert(readAttribute(file, "energy", "ENERGY", 1));
    attr.insert(readAttribute(file, "charge", "CHARGE", 1));
    attr.insert(readAttribute(file, "path_length", "SPOS", 1));
    attr.insert(readAttribute(file, "time", "TIME", 1));

    return attr;
}

void convertToElegantUnits(data_t &data) {
    size_t size;
    bool temporal = false;

    try {
        size = data["z"].size();
    } catch (...) {
        temporal = true;
        size = data["x"].size();
    }

    for (size_t i = 0; i < size; ++ i) {
        data["px"][i] /= data["pz"][i];
        data["py"][i] /= data["pz"][i];
        if (temporal) continue;

        data["z"][i] /= -299792458.0;
    }
}

void writeSDDSFile(const std::string &fname, const data_t &data, const attributes_t &attr, FORMAT form) {
    std::map<std::string, std::string> nameConversion {{"x", "x"},
                                                       {"y", "y"},
                                                       {"t", "z"},
                                                       {"xp", "px"},
                                                       {"yp", "py"},
                                                       {"p", "pz"}};
    const std::map<std::string, std::string> nameUnits {{"x", "m"},
                                                        {"y", "m"},
                                                        {"t", "s"},
                                                        {"xp", "1"},
                                                        {"yp", "1"},
                                                        {"p", "beta * gamma"}};
    const std::map<std::string, std::string> nameSymbol {{"x", "x"},
                                                        {"y", "y"},
                                                        {"t", "t"},
                                                        {"xp", "x'"},
                                                        {"yp", "y'"},
                                                        {"p", "p"}};
    SDDS_DATASET SDDS_dataset;
    const long rows = data.at("x").size();
    std::vector<std::vector<double> > fileData;

    char buffer0[256];
    char buffer1[64];
    char buffer2[8];

    for (const auto entry: data) {
        if (entry.first == "time") {
            nameConversion["t"] = "time";
            break;
        }
    }

    strcpy(buffer0, fname.c_str());
    strcpy(buffer1, "extracted OPAL data");
    if (SDDS_InitializeOutput(&SDDS_dataset, form, 1, buffer1, NULL, buffer0 ) != 1) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fflush(stdout);
        std::cerr << "Error: couldn't open file '" << fname << "'\n" << std::endl;
        std::exit(1);
    }

    for (auto at: attr) {
        strcpy(buffer0, at.name.c_str());
        strcpy(buffer1, at.unit.c_str());
        if (SDDS_DefineSimpleParameter(&SDDS_dataset, buffer0, buffer1, SDDS_DOUBLE) != 1) {
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            fflush(stdout);
            std::cerr << "Error: couldn't append parameter '" << at.name << "'\n" << std::endl;
            std::exit(1);
        }
    }

    for (auto names: nameConversion) {
        strcpy(buffer0, names.first.c_str());
        strcpy(buffer2, nameSymbol.at(names.first).c_str());
        strcpy(buffer1, nameUnits.at(names.first).c_str());
        if (SDDS_DefineColumn(&SDDS_dataset,
                              buffer0,
                              buffer2,
                              buffer1,
                              NULL,
                              NULL,
                              SDDS_DOUBLE,
                              0) == -1) {
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            fflush(stdout);
            std::cerr << "Error: couldn't append column '" << names.first << "'\n" << std::endl;
            std::exit(1);
        }
    }

    if (SDDS_WriteLayout(&SDDS_dataset) != 1) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fflush(stdout);
        std::cerr << "Error: couldn't write layout\n" << std::endl;
        std::exit(1);
    }

    if (SDDS_StartPage(&SDDS_dataset, rows) != 1) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fflush(stdout);
        std::cerr << "Error: couldn't start page\n" << std::endl;
        std::exit(1);
    }

    for (auto at: attr) {
        if (SDDS_SetParameters(&SDDS_dataset, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE,
                               at.name.c_str(), at.value[0], NULL) != 1) {
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            fflush(stdout);
            std::cerr << "Error: couldn't write parameter '" << at.name << "'\n" << std::endl;
            std::exit(1);
        }
    }

    for (long i = 0; i < rows; ++ i) {
        if (SDDS_SetRowValues(&SDDS_dataset, SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE, i,
                              "x", data.at(nameConversion.at("x"))[i],
                              "y", data.at(nameConversion.at("y"))[i],
                              "t", data.at(nameConversion.at("t"))[i],
                              "xp", data.at(nameConversion.at("xp"))[i],
                              "yp", data.at(nameConversion.at("yp"))[i],
                              "p", data.at(nameConversion.at("p"))[i],
                              NULL) != 1) {
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            fflush(stdout);
            std::cerr << "Error: couldn't add row\n" << std::endl;
            std::exit(1);
        }
    }

    if (SDDS_WritePage(&SDDS_dataset) != 1) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fflush(stdout);
        std::cerr << "Error: couldn't write page\n" << std::endl;
        std::exit(1);
    }

    if (SDDS_Terminate(&SDDS_dataset) != 1) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        fflush(stdout);
        std::cerr << "Error: couldn't terminate properly\n" << std::endl;
        std::exit(1);
    }
}

void printInfo(const std::string &input) {
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = MPI_COMM_WORLD;
    H5SetPropFileMPIOCollective (props, &comm);

    file_t file = H5OpenFile(input.c_str(), H5_O_RDONLY, props);
    H5CloseProp (props);
    h5_ssize_t numStepsInSource = H5GetNumSteps(file);

    std::cout << std::left << std::setw(15) << "Step number" << std::setw(15) << "Position [m]" << std::endl;
    std::cout << std::setfill('-') << std::setw(29) << "-" << std::endl;
    std::cout << std::setfill(' ');
    for (h5_ssize_t i = 0; i < numStepsInSource; ++ i) {
        REPORTONERROR(H5SetStep(file, i));
        double spos;
        READSTEPATTRIB(Float64, file, "SPOS", &spos);
        std::cout << std::setw(15) << i << std::setw(15) << std::setprecision(5) << spos << std::endl;
    }
    REPORTONERROR(H5CloseFile(file));
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
              << indent << std::setw(width) << std::left << "--time-of-arrival"       << "= extract and save the time of arrival instead of the longitudinal position\n"
              << indent << std::setw(width) << std::left << "--help"                  << "= this help\n"
              << std::endl;
}