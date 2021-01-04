#include "Distribution/FromFile.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <strstream>
#include <string>
extern Inform* gmsg;

FromFile::FromFile()
{

}

size_t FromFile::getNumberOfParticlesInFile(std::ifstream &inputFile) {

    size_t numberOfParticlesRead = 0;
    if (Ippl::myNode() == 0) {
        const boost::regex commentExpr("[[:space:]]*#.*");
        const std::string repl("");
        std::string line;
        std::stringstream linestream;
        long tempInt = 0;

        do {
            std::getline(inputFile, line);
            line = boost::regex_replace(line, commentExpr, repl);
        } while (line.length() == 0 && !inputFile.fail());

        linestream.str(line);
        linestream >> tempInt;
        if (tempInt <= 0) {
            throw OpalException("FromFile::getNumberOfParticlesInFile",
                                "The file '" +
                                Attributes::getString(itsAttr[Attrib::FromFile::FNAME]) +
                                "' does not seem to be an ASCII file containing a distribution.");
        }
        numberOfParticlesRead = static_cast<size_t>(tempInt);
    }
    reduce(numberOfParticlesRead, numberOfParticlesRead, OpAddAssign());

    return numberOfParticlesRead;
}

std::vector<DistributionBase::ParticlePhaseSpace> FromFile::create(size_t numberOfParticles,
                                                                   double massInEv,
                                                                   double charge)
{
    
    // Data input file is only read by node 0.
    std::ifstream inputFile;
    std::string fileName = Attributes::getString(itsAttr[Attrib::FromFile::FNAME]);
    if (!boost::filesystem::exists(fileName)) {
        throw OpalException("FromFile::create",
                            "Open file operation failed, please check if \""
                            + fileName +
                            "\" really exists.");
    }
    if (Ippl::myNode() == 0) {
        inputFile.open(fileName.c_str());
    }

    *gmsg << level3 << "\n"
          << "------------------------------------------------------------------------------------\n";
    *gmsg << "READ INITIAL DISTRIBUTION FROM FILE \""
          << Attributes::getString(itsAttr[Attrib::FromFile::FNAME])
          << "\"\n";
    *gmsg << "------------------------------------------------------------------------------------\n" << endl;

    size_t numberOfParticlesRead = getNumberOfParticlesInFile(inputFile);
    /*
     * We read in the particle information using node zero, but save the particle
     * data to each processor in turn.
     */
    unsigned int saveProcessor = 0;

    pmean_m = 0.0;

    size_t numPartsToSend = 0;
    unsigned int distributeFrequency = 1000;
    size_t singleDataSize = (/*sizeof(int) +*/ 6 * sizeof(double));
    unsigned int dataSize = distributeFrequency * singleDataSize;
    std::vector<char> data;

    data.reserve(dataSize);
    std::vector<ParticlePhaseSpace> result;
    const char* buffer;
    if (Ippl::myNode() == 0) {
        char lineBuffer[1024];
        unsigned int numParts = 0;
        while (!inputFile.eof()) {
            inputFile.getline(lineBuffer, 1024);

            Vector_t R(0.0), P(0.0);

            std::istringstream line(lineBuffer);
            line >> R(0);
            if (line.rdstate()) break;
            line >> P(0);
            line >> R(1);
            line >> P(1);
            line >> R(2);
            line >> P(2);

            if (saveProcessor >= (unsigned)Ippl::getNodes())
                saveProcessor = 0;

            if (inputMoUnits_m == InputMomentumUnitsT::EV) {
                P(0) = Util::convertMomentumeVToBetaGamma(P(0), massIneV);
                P(1) = Util::convertMomentumeVToBetaGamma(P(1), massIneV);
                P(2) = Util::convertMomentumeVToBetaGamma(P(2), massIneV);
            }

            pmean_m += P;

            if (saveProcessor > 0u) {
                buffer = reinterpret_cast<const char*>(&R[0]);
                data.insert(data.end(), buffer, buffer + 3 * sizeof(double));
                buffer = reinterpret_cast<const char*>(&P[0]);
                data.insert(data.end(), buffer, buffer + 3 * sizeof(double));
                ++ numPartsToSend;

                if (numPartsToSend % distributeFrequency == 0) {
                    MPI_Bcast(&dataSize, 1, MPI_UNSIGNED, 0, Ippl::getComm());
                    MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());
                    numPartsToSend = 0;

                    std::vector<char>().swap(data);
                    data.reserve(dataSize);
                }
            } else {
                ParticlePhaseSpace particle;
                particle.R_m = R;
                particle.P_m = P;
                result.push_back(particle);
            }

            ++ numParts;
            ++ saveProcessor;
        }

        dataSize = (numberOfParticlesRead == numParts? data.size(): std::numeric_limits<unsigned int>::max());
        MPI_Bcast(&dataSize, 1, MPI_UNSIGNED, 0, Ippl::getComm());
        if (numberOfParticlesRead != numParts) {
            throw OpalException("FromFile::create",
                                "Found " +
                                std::to_string(numParts) +
                                " particles in file '" +
                                fileName +
                                "' instead of " +
                                std::to_string(numberOfParticlesRead));
        }
        MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());

    } else {
        do {
            MPI_Bcast(&dataSize, 1, MPI_UNSIGNED, 0, Ippl::getComm());
            if (dataSize == std::numeric_limits<unsigned int>::max()) {
                throw OpalException("FromFile::create",
                                    "Couldn't find " +
                                    std::to_string(numberOfParticlesRead) +
                                    " particles in file '" +
                                    fileName + "'");
            }
            MPI_Bcast(&data[0], dataSize, MPI_CHAR, 0, Ippl::getComm());

            size_t i = 0;
            while (i < dataSize) {

                if (saveProcessor + 1 == (unsigned) Ippl::myNode()) {
                    const double *tmp = reinterpret_cast<const double*>(&data[i]);
                    ParticlePhaseSpace particle;
                    particle.R_m = Vector_t(tmp[0], tmp[1], tmp[2]);
                    particle.P_m = Vector_t(tmp[3], tmp[4], tmp[5]);
                    result.push_back(particle);
                    i += 6 * sizeof(double);
                } else {
                    i += singleDataSize;
                }

                ++ saveProcessor;
                if (saveProcessor + 1 >= (unsigned) Ippl::getNodes()) {
                    saveProcessor = 0;
                }
            }
        } while (dataSize == distributeFrequency * singleDataSize);
    }

    pmean_m /= numberOfParticlesRead;
    reduce(pmean_m, pmean_m, OpAddAssign());

    if (Ippl::myNode() == 0)
        inputFile.close();

    return result;
}

Inform& FromFile::print(Inform& os) const
{

}

double FromFile::getEmissionTime() const
{

}