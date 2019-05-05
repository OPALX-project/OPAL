#include "Fields/Fieldmap.h"
#include "Fields/FM3DDynamic.h"
#include "Fields/FM3DH5Block.h"
#include "Fields/FM3DH5Block_nonscale.h"
#include "Fields/FM3DMagnetoStaticH5Block.h"
#include "Fields/FM3DMagnetoStatic.h"
#include "Fields/FM3DMagnetoStaticExtended.h"
#include "Fields/FM2DDynamic.h"
#include "Fields/FM2DElectroStatic.h"
#include "Fields/FM2DMagnetoStatic.h"
#include "Fields/FM1DDynamic.h"
#include "Fields/FM1DDynamic_fast.h"
#include "Fields/Astra1DDynamic.h"
#include "Fields/Astra1DDynamic_fast.h"
#include "Fields/FM1DElectroStatic.h"
#include "Fields/FM1DElectroStatic_fast.h"
#include "Fields/Astra1DElectroStatic.h"
#include "Fields/Astra1DElectroStatic_fast.h"
#include "Fields/FM1DMagnetoStatic.h"
#include "Fields/FM1DMagnetoStatic_fast.h"
#include "Fields/Astra1DMagnetoStatic.h"
#include "Fields/Astra1DMagnetoStatic_fast.h"
#include "Fields/FM1DProfile1.h"
#include "Fields/FM1DProfile2.h"
#include "Fields/FMDummy.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Options.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"

#include "H5hut.h"

#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <ios>

#include <assert.h>

namespace fs = boost::filesystem;

#define REGISTER_PARSE_TYPE(X) template <> struct Fieldmap::TypeParseTraits<X> \
    { static const char* name; } ; const char* Fieldmap::TypeParseTraits<X>::name = #X

Fieldmap *Fieldmap::getFieldmap(std::string Filename, bool fast) {
    std::map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if(position != FieldmapDictionary.end()) {
        (*position).second.RefCounter++;
        return (*position).second.Map;
    } else {
        MapType type;
        std::pair<std::map<std::string, FieldmapDescription>::iterator, bool> position;
        type = readHeader(Filename);
        switch(type) {
        case T1DDynamic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DDynamic, new FM1DDynamic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DDynamic, new FM1DDynamic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case TAstraDynamic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraDynamic, new Astra1DDynamic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraDynamic, new Astra1DDynamic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case T1DElectroStatic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DElectroStatic, new FM1DElectroStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DElectroStatic, new FM1DElectroStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case TAstraElectroStatic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraElectroStatic, new Astra1DElectroStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraElectroStatic, new Astra1DElectroStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case T1DMagnetoStatic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DMagnetoStatic, new FM1DMagnetoStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T1DMagnetoStatic, new FM1DMagnetoStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case TAstraMagnetoStatic:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraMagnetoStatic, new Astra1DMagnetoStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            TAstraMagnetoStatic, new Astra1DMagnetoStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;

        case T1DProfile1:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(T1DProfile1, new FM1DProfile1(Filename))));
            return (*position.first).second.Map;
            break;

        case T1DProfile2:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(
                        T1DProfile2, new FM1DProfile2(Filename))));
            return (*position.first).second.Map;
            break;

        case T2DDynamic:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(
                        T2DDynamic, new FM2DDynamic(Filename))));
            return (*position.first).second.Map;
            break;

        case T2DElectroStatic:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(
                        T2DElectroStatic, new FM2DElectroStatic(Filename))));
            return (*position.first).second.Map;
            break;

        case T2DMagnetoStatic:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(
                        T2DMagnetoStatic, new FM2DMagnetoStatic(Filename))));
            return (*position.first).second.Map;
            break;

        case T3DDynamic:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename,
                    FieldmapDescription(
                        T3DDynamic, new FM3DDynamic(Filename))));
            return (*position.first).second.Map;
            break;

        case T3DMagnetoStaticH5Block:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename, FieldmapDescription(
                        T3DMagnetoStaticH5Block, new FM3DMagnetoStaticH5Block(Filename))));
            return (*position.first).second.Map;
            break;

        case T3DMagnetoStatic:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename, FieldmapDescription(
                        T3DMagnetoStatic, new FM3DMagnetoStatic(Filename))));
            return (*position.first).second.Map;
            break;

        case T3DMagnetoStatic_Extended:
            position = FieldmapDictionary.insert(
                std::make_pair(
                    Filename, FieldmapDescription(
                        T3DMagnetoStatic_Extended, new FM3DMagnetoStaticExtended(Filename))));
            return (*position.first).second.Map;
            break;

        case T3DDynamicH5Block:
            if(fast) {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T3DDynamic, new FM3DH5Block_nonscale(Filename))));
            } else {
                position = FieldmapDictionary.insert(
                    std::make_pair(
                        Filename,
                        FieldmapDescription(
                            T3DDynamic, new FM3DH5Block(Filename))));
            }
            return (*position.first).second.Map;
            break;

        default:
            throw GeneralClassicException("Fieldmap::getFieldmap()",
                                          "Couldn't determine type of fieldmap in file \"" + Filename + "\"");
        }
    }
}

std::vector<std::string> Fieldmap::getListFieldmapNames() {
    std::vector<std::string> name_list;
    for(std::map<std::string, FieldmapDescription>::const_iterator it = FieldmapDictionary.begin(); it != FieldmapDictionary.end(); ++ it) {
        name_list.push_back((*it).first);
    }
    return name_list;
}

void Fieldmap::deleteFieldmap(std::string Filename) {
    freeMap(Filename);
}

void Fieldmap::clearDictionary() {
    std::map<std::string, FieldmapDescription>::iterator it = FieldmapDictionary.begin();
    for (;it != FieldmapDictionary.end(); ++ it) {
        delete it->second.Map;
        it->second.Map = NULL;
    }
    FieldmapDictionary.clear();
}

MapType Fieldmap::readHeader(std::string Filename) {
    char magicnumber[5] = "    ";
    std::string buffer;
    int lines_read_m = 0;

    // Check for default map(s).
    if(Filename == "1DPROFILE1-DEFAULT")
        return T1DProfile1;

    if (!fs::exists(Filename))
        throw GeneralClassicException("Fieldmap::readHeader()",
                                      "File \"" + Filename + "\" doesn't exist");

    std::ifstream File(Filename.c_str());
    if(!File.good()) {
        std::cerr << "could not open file " << Filename << std::endl;
        return UNKNOWN;
    }

    getLine(File, lines_read_m, buffer);
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter.read(magicnumber, 4);

    if(strcmp(magicnumber, "3DDy") == 0)
        return T3DDynamic;

    if(strcmp(magicnumber, "3DMa") == 0) {
        char tmpString[21] = "                    ";
        interpreter.read(tmpString, 20);

        if(strcmp(tmpString, "gnetoStatic_Extended") == 0)
            return T3DMagnetoStatic_Extended;
        else
            return T3DMagnetoStatic;
    }

    if(strcmp(magicnumber, "3DEl") == 0)
        return T3DElectroStatic;

    if(strcmp(magicnumber, "2DDy") == 0) {
        // char tmpString[14] = "             ";
        // interpreter.read(tmpString, 13);
        return T2DDynamic;
    }

    if(strcmp(magicnumber, "2DMa") == 0) {
        // char tmpString[20] = "                   ";
        // interpreter.read(tmpString, 19);
        return T2DMagnetoStatic;
    }

    if(strcmp(magicnumber, "2DEl") == 0) {
        // char tmpString[20] = "                   ";
        // interpreter.read(tmpString, 19);
        return T2DElectroStatic;
    }

    if(strcmp(magicnumber, "1DDy") == 0)
        return T1DDynamic;

    if(strcmp(magicnumber, "1DMa") == 0)
        return T1DMagnetoStatic;

    if(strcmp(magicnumber, "1DPr") == 0) {
        // char tmpString[7] = "      ";
        // interpreter.read(tmpString, 6);
        // if(strcmp(tmpString, "ofile1") == 0)
        return T1DProfile1;
        // if(strcmp(tmpString, "ofile2") == 0)
        //     return T1DProfile2;
    }

    if(strcmp(magicnumber, "1DEl") == 0)
        return T1DElectroStatic;

    if(strcmp(magicnumber, "\211HDF") == 0) {
        h5_err_t h5err = 0;
#if defined (NDEBUG)
        // mark variable as unused
        (void)h5err;
#endif
        char name[20];
        h5_size_t len_name = sizeof(name);
        h5_prop_t props = H5CreateFileProp ();
        MPI_Comm comm = Ippl::getComm();
        h5err = H5SetPropFileMPIOCollective (props, &comm);
        assert (h5err != H5_ERR);
        h5_file_t file = H5OpenFile (Filename.c_str(), H5_O_RDONLY, props);
        assert (file != (h5_file_t)H5_ERR);
        H5CloseProp (props);

        h5err = H5SetStep(file, 0);
        assert (h5err != H5_ERR);

        h5_int64_t num_fields = H5BlockGetNumFields(file);
        assert (num_fields != H5_ERR);
        MapType maptype = UNKNOWN;

        for(h5_ssize_t i = 0; i < num_fields; ++ i) {
            h5err = H5BlockGetFieldInfo(
                file, (h5_size_t)i, name, len_name, NULL, NULL, NULL, NULL);
            assert (h5err != H5_ERR);
            // using field name "Bfield" and "Hfield" to distinguish the type
            if(strcmp(name, "Bfield") == 0) {
                maptype = T3DMagnetoStaticH5Block;
                break;
            } else if(strcmp(name, "Hfield") == 0) {
                maptype = T3DDynamicH5Block;
                break;
            }
        }
        h5err = H5CloseFile(file);
        assert (h5err != H5_ERR);
        if(maptype != UNKNOWN)
            return maptype;
    }

    if(strcmp(magicnumber, "Astr") == 0) {
        char tmpString[3] = "  ";
        interpreter.read(tmpString, 2);
        if(strcmp(tmpString, "aE") == 0) {
            return TAstraElectroStatic;
        }
        if(strcmp(tmpString, "aM") == 0) {
            return TAstraMagnetoStatic;
        }
        if(strcmp(tmpString, "aD") == 0) {
            return TAstraDynamic;
        }
    }


    return UNKNOWN;
}

void Fieldmap::readMap(std::string Filename) {
    std::map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if(position != FieldmapDictionary.end())
        if(!(*position).second.read) {
            (*position).second.Map->readMap();
            (*position).second.read = true;
        }
}

void Fieldmap::freeMap(std::string Filename) {
    std::map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    /*
      FIXME: find( ) make problem, crashes
    */
    if(position != FieldmapDictionary.end()) {
        if((*position).second.RefCounter > 0) {
            (*position).second.RefCounter--;
        }

        if ((*position).second.RefCounter == 0) {
            delete (*position).second.Map;
            (*position).second.Map = NULL;
            FieldmapDictionary.erase(position);
        }
    }
}

void Fieldmap::checkMap(unsigned int accuracy,
                        std::pair<double, double> fieldDimensions,
                        double deltaZ,
                        const std::vector<double> &fourierCoefficients,
                        gsl_spline *splineCoefficients,
                        gsl_interp_accel *splineAccelerator) {
    double length = fieldDimensions.second - fieldDimensions.first;
    unsigned int sizeSampling = std::floor(length / deltaZ + 0.5);
    std::vector<double> zSampling(sizeSampling);
    zSampling[0] = fieldDimensions.first;
    for (unsigned int i = 1; i < sizeSampling; ++ i) {
        zSampling[i] = zSampling[i-1] + deltaZ;
    }
    checkMap(accuracy, length, zSampling, fourierCoefficients, splineCoefficients, splineAccelerator);
}

void Fieldmap::checkMap(unsigned int accuracy,
                        double length,
                        const std::vector<double> &zSampling,
                        const std::vector<double> &fourierCoefficients,
                        gsl_spline *splineCoefficients,
                        gsl_interp_accel *splineAccelerator) {
    double error = 0.0;
    double maxDiff = 0.0;
    double ezMax = 0.0;
    double ezSquare = 0.0;
    size_t lastDot = Filename_m.find_last_of(".");
    size_t lastSlash = Filename_m.find_last_of("/");
    lastSlash = (lastSlash == std::string::npos)? 0: lastSlash + 1;

    auto opal = OpalData::getInstance();
    std::ofstream out;
    if (Ippl::myNode() == 0 && !opal->isOptimizerRun()) {
        out.open("data/" + Filename_m.substr(lastSlash, lastDot) + ".check");
        out << "# z  original reproduced\n";
    }
    auto it = zSampling.begin();
    auto end = zSampling.end();
    for (; it != end; ++ it) {
        const double kz = Physics::two_pi * (*it / length + 0.5);
        double onAxisFieldCheck = fourierCoefficients[0];
        unsigned int n = 1;
        for(unsigned int l = 1; l < accuracy; ++l, n += 2) {
            double coskzl = cos(kz * l);
            double sinkzl = sin(kz * l);

            onAxisFieldCheck += (fourierCoefficients[n] * coskzl - fourierCoefficients[n + 1] * sinkzl);
        }
        double ez = gsl_spline_eval(splineCoefficients, *it, splineAccelerator);
        double difference = std::abs(ez - onAxisFieldCheck);
        maxDiff = difference > maxDiff? difference: maxDiff;
        ezMax = std::abs(ez) > ezMax? std::abs(ez): ezMax;
        error += std::pow(difference, 2.0);
        ezSquare += std::pow(ez, 2.0);

        if (Ippl::myNode() == 0 && !opal->isOptimizerRun()) {
            out << std::setw(16) << std::setprecision(8) << *it
                << std::setw(16) << std::setprecision(8) << ez
                << std::setw(16) << std::setprecision(8) << onAxisFieldCheck
                << std::endl;
        }
    }
    out.close();

    if (sqrt(error / ezSquare) > 1e-1 || maxDiff > 1e-1 * ezMax) {
        lowResolutionWarning(sqrt(error / ezSquare), maxDiff / ezMax);

        throw GeneralClassicException("Astra2DDynamic_fast::readMap()",
                                      "Field map can't be reproduced properly with the given number of fourier components");
    }
    if (sqrt(error / ezSquare) > 1e-2 || maxDiff > 1e-2 * ezMax) {
        lowResolutionWarning(sqrt(error / ezSquare), maxDiff / ezMax);
    }
}

void Fieldmap::setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle)
{};

void Fieldmap::setFieldLength(const double &)
{};

void Fieldmap::getLine(std::ifstream &in, int &lines_read, std::string &buffer) {
    size_t firstof = 0;
    size_t lastof;
    size_t comment;

    do {
        ++ lines_read;
        in.getline(buffer_m, READ_BUFFER_LENGTH);

        buffer = std::string(buffer_m);

        comment = buffer.find("#");
        buffer = buffer.substr(0, comment);

        lastof = buffer.find_last_of(alpha_numeric);
        firstof = buffer.find_first_of(alpha_numeric);
    } while(!in.eof() && lastof == std::string::npos);

    if(firstof != std::string::npos) {
        buffer = buffer.substr(firstof, lastof - firstof + 1);
    }
}

bool Fieldmap::interpreteEOF(std::ifstream &in) {
    while(!in.eof()) {
        ++lines_read_m;
        in.getline(buffer_m, READ_BUFFER_LENGTH);
        std::string buffer(buffer_m);
        size_t comment = buffer.find_first_of("#");
        buffer = buffer.substr(0, comment);
        size_t lasto = buffer.find_first_of(alpha_numeric);
        if(lasto != std::string::npos) {
            exceedingValuesWarning();
            return false;
        }
    }
    return true;
}

void Fieldmap::interpreteWarning(const std::string &error_msg,
                                 const std::string &expecting,
                                 const std::string &found) {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "expecting: '" << expecting << "' on line " << lines_read_m << ",\n"
             << "found instead: '" << found << "'.";
    //    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    throw GeneralClassicException("Fieldmap::interpretWarning()",
                                  errormsg.str());
}

void Fieldmap::interpreteWarning(const std::ios_base::iostate &state,
                                 const bool &read_all,
                                 const std::string &expecting,
                                 const std::string &found) {
    std::string error_msg;
    if(!read_all) {
        error_msg = std::string("Didn't find enough values!");
    } else if(state & std::ios_base::eofbit) {
        error_msg = std::string("Found more values than expected!");
    } else if(state & std::ios_base::failbit) {
        error_msg = std::string("Found wrong type of values!");
    }
    interpreteWarning(error_msg, expecting, found);
}

void Fieldmap::missingValuesWarning() {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are only " << lines_read_m - 1 << " lines in the file, expecting more.\n"
             << "Please check the section about field maps in the user manual.";
    //    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    throw GeneralClassicException("Fieldmap::missingValuesWarning()",
                                  errormsg.str());
}

void Fieldmap::exceedingValuesWarning() {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are too many lines in the file, expecting only " << lines_read_m << " lines.\n"
             << "Please check the section about field maps in the user manual.";
    //    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    throw GeneralClassicException("Fieldmap::exceedingValuesWarning()",
                                  errormsg.str());
}

void Fieldmap::disableFieldmapWarning() {
    std::stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" + Filename_m + "' DUE TO PARSING ERRORS." ;
    //    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    throw GeneralClassicException("Fieldmap::disableFieldmapsWarning()",
                                  errormsg.str());
}

void Fieldmap::noFieldmapWarning() {
    std::stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" << Filename_m << "' SINCE FILE COULDN'T BE FOUND!";
    //    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    throw GeneralClassicException("Fieldmap::noFieldmapsWarning()",
                                  errormsg.str());
}

void Fieldmap::lowResolutionWarning(double squareError, double maxError) {
    std::stringstream errormsg;
    errormsg << "IT SEEMS THAT YOU USE TOO FEW FOURIER COMPONENTS TO SUFFICIENTLY WELL\n"
             << "RESOLVE THE FIELD MAP '" << Filename_m << "'.\n"
             << "PLEASE INCREASE THE NUMBER OF FOURIER COMPONENTS!\n"
             << "The ratio (e_i - E_i)^2 / E_i^2 is " << std::to_string(squareError) << " and\n"
             << "the ratio (max_i(|e_i - E_i|) / max_i(|E_i|) is " << std::to_string(maxError) << ".\n"
             << "Here E_i is the field as in the field map and e_i is the reconstructed field.\n"
             << "The lower limit for the two ratios is 1e-2\n"
             << "Have a look into the directory data/ for a reconstruction of the field map.\n";
    std::string errormsg_str = typeset_msg(errormsg.str(), "warning");

    ERRORMSG(errormsg_str << "\n" << endl);

    if(Ippl::myNode() == 0) {
        std::ofstream omsg("errormsg.txt", std::ios_base::app);
        omsg << errormsg.str() << std::endl;
        omsg.close();
    }
}

std::string Fieldmap::typeset_msg(const std::string &msg, const std::string &title) {
    static std::string frame("* ******************************************************************************\n");
    static unsigned int frame_width = frame.length() - 5;
    static std::string closure("                                                                               *\n");

    std::string return_string("\n" + frame);

    int remaining_length = msg.length();
    unsigned int current_position = 0;

    unsigned int ii = 0;
    for(; ii < title.length(); ++ ii) {
        char c = title[ii];
        c = toupper(c);
        return_string.replace(15 + 2 * ii, 1, " ");
        return_string.replace(16 + 2 * ii, 1, &c, 1);
    }
    return_string.replace(15 + 2 * ii, 1, " ");

    while(remaining_length > 0) {
        size_t eol = msg.find("\n", current_position);
        std::string next_to_process;
        if(eol != std::string::npos) {
            next_to_process = msg.substr(current_position, eol - current_position);
        } else {
            next_to_process = msg.substr(current_position);
            eol = msg.length();
        }

        if(eol - current_position < frame_width) {
            return_string += "* " + next_to_process + closure.substr(eol - current_position + 2);
        } else {
            unsigned int last_space = next_to_process.rfind(" ", frame_width);
            if(last_space > 0) {
                if(last_space < frame_width) {
                    return_string += "* " + next_to_process.substr(0, last_space) + closure.substr(last_space + 2);
                } else {
                    return_string += "* " + next_to_process.substr(0, last_space) + " *\n";
                }
                if(next_to_process.length() - last_space + 1 < frame_width) {
                    return_string += "* " + next_to_process.substr(last_space + 1) + closure.substr(next_to_process.length() - last_space + 1);
                } else {
                    return_string += "* " + next_to_process.substr(last_space + 1) + " *\n";
                }
            } else {
                return_string += "* " + next_to_process + " *\n";
            }
        }

        current_position = eol + 1;
        remaining_length = msg.length() - current_position;
    }
    return_string += frame;

    return return_string;
}

void Fieldmap::getOnaxisEz(std::vector<std::pair<double, double> > &onaxis)
{ }

void Fieldmap::get1DProfile1EngeCoeffs(std::vector<double> &engeCoeffsEntry,
                                       std::vector<double> &engeCoeffsExit) {

}

void Fieldmap::get1DProfile1EntranceParam(double &entranceParameter1,
        double &entranceParameter2,
        double &entranceParameter3) {

}

void Fieldmap::get1DProfile1ExitParam(double &exitParameter1,
                                      double &exitParameter2,
                                      double &exitParameter3) {

}

double Fieldmap::getFieldGap() {
    return 0.0;
}

void Fieldmap::setFieldGap(double gap) {

}

void Fieldmap::write3DField(unsigned int nx,
                            unsigned int ny,
                            unsigned int nz,
                            const std::pair<double, double> &xrange,
                            const std::pair<double, double> &yrange,
                            const std::pair<double, double> &zrange,
                            const std::vector<Vector_t> &ef,
                            const std::vector<Vector_t> &bf) {
    const size_t numpoints = nx * ny * nz;
    if (Ippl::myNode() != 0 ||
        (ef.size() != numpoints && bf.size() != numpoints)) return;

    size_t extensionStart = Filename_m.find_last_of('.');
    std::ofstream of;
    of.open (std::string ("data/" + Filename_m.substr(0,extensionStart) + ".vtk").c_str ());
    assert (of.is_open ());
    of.precision (6);

    const double hx = (xrange.second - xrange.first) / (nx - 1);
    const double hy = (yrange.second - yrange.first) / (ny - 1);
    const double hz = (zrange.second - zrange.first) / (nz - 1);

    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated by 3D fieldmaps" << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET RECTILINEAR_GRID" << std::endl;
    of << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;

    of << "X_COORDINATES " << nx << " float" << std::endl;
    of << xrange.first;
    for (unsigned int i = 1; i < nx - 1; ++ i) {
        of << " " << xrange.first + i * hx;
    }
    of << " " << xrange.second << std::endl;

    of << "Y_COORDINATES " << ny << " float" << std::endl;
    of << yrange.first;
    for (unsigned int i = 1; i < ny - 1; ++ i) {
        of << " " << yrange.first + i * hy;
    }
    of << " " << yrange.second << std::endl;

    of << "Z_COORDINATES " << nz << " float" << std::endl;
    of << zrange.first;
    for (unsigned int i = 1; i < nz - 1; ++ i) {
        of << " " << zrange.first + i * hz;
    }
    of << " " << zrange.second << std::endl;

    of << "POINT_DATA " << numpoints << std::endl;

    if (ef.size() == numpoints) {
        of << "VECTORS EField float" << std::endl;
        // of << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < numpoints; ++ i) {
            of << ef[i](0) << " " << ef[i](1) << " " << ef[i](2) << std::endl;
        }
        // of << std::endl;
    }

    if (bf.size() == numpoints) {
        of << "VECTORS BField float" << std::endl;
        // of << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < numpoints; ++ i) {
            of << bf[i](0) << " " << bf[i](1) << " " << bf[i](2) << std::endl;
        }
        // of << std::endl;
    }
}

REGISTER_PARSE_TYPE(int);
REGISTER_PARSE_TYPE(unsigned int);
REGISTER_PARSE_TYPE(double);
REGISTER_PARSE_TYPE(std::string);

std::string Fieldmap::alpha_numeric("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-+\211");
std::map<std::string, Fieldmap::FieldmapDescription> Fieldmap::FieldmapDictionary = std::map<std::string, Fieldmap::FieldmapDescription>();
char Fieldmap::buffer_m[READ_BUFFER_LENGTH];