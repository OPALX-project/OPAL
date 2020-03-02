#ifndef __GENERATE_SIMULATION_H__
#define __GENERATE_SIMULATION_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cstdlib>

#include <vector>
#include <boost/algorithm/string.hpp>
#include "boost/filesystem.hpp"

#include "Util/OptPilotException.h"

/**
 *  \class GenerateOpalSimulation
 *  \brief Generates an OPAL input file from data and template file.
 *
 *  When running the optimizer with the OPAL forward solvers this utility
 *  class helps generating OPAL input files corresponding to the design
 *  variable values of the requested evaluation.
 *  Using the data file we start by generating the dictionary of known
 *  templates.
 *  When the input file is requested, all template variables in the template
 *  input file are replaced with values found in the dictionary.
 *
 *  TODO:
 *    - pass stream directly to opal (detour over file not necessary)
 */
class GenerateOpalSimulation {

public:

    /**
     *  Sets up the scaled variables and fills the dictionary with content of
     *  data file.
     *
     *  @param[in] tmplFile filename of template file
     *  @param[in] varDictionary filename of data file
     *  @param[in] userValues template variables specified by the optimizer
     */
    GenerateOpalSimulation(std::string tmplFile, std::string varDictionary,
                           std::map<std::string,std::string> userValues) {
        varDictionary_ = varDictionary;
        tmplFile_ = tmplFile;

        // some variables need to be scaled
        scaleVars_.insert(std::pair<std::string, double>("GUNSOLB", 1.));

        // add user values and then fill using the data file without
        // overwriting user values..
        dictionary_.insert(userValues.begin(),userValues.end());
        fillDictionary();
    }

    ~GenerateOpalSimulation() {
        dictionary_.clear();
        scaleVars_.clear();
    }

    /**
     *  Replace all template variables in template file with corresponding
     *  values from the dictionary and write file to disk
     *
     *  @param[in] outputFile write resulting input file to this file
     */
    void writeInputFile(std::string outputFile) {

        std::ifstream infile(tmplFile_.c_str());
        std::ostringstream outdata;
        outdata.precision(15);

        while(infile.good()) {
            std::string line;
            std::getline(infile, line, '\n');

            //XXX doing the inverse would be better
            std::map<std::string, std::string>::iterator itr = dictionary_.begin();
            for(;itr != dictionary_.end(); itr++) {
                size_t pos = line.find("_" + itr->first + "_");
                while(pos != std::string::npos) {
                    line.replace(pos, itr->first.length() + 2, itr->second);
                    pos = line.find("_" + itr->first + "_");
                }
            }

            outdata << line << std::endl;
        }
        infile.close();

        // ensure the contents are written to disk
        std::ofstream outfile(outputFile.c_str());
        outfile.precision(15);
        outfile << outdata.str();
        outfile.flush();
        outfile.rdbuf()->pubsync();

        outfile.close();
    }


private:

    /// template filename
    std::string tmplFile_;
    /// data filename
    std::string varDictionary_;
    /// holds parsed template data
    std::map<std::string, std::string> dictionary_;
    /// holds all variables that need scaling
    std::map<std::string, double> scaleVars_;

    /**
     *  Parses the data file and insert all values into the dictionary.
     *
     *  The data file needs to conform to the specification:
     *    - comments start with a hash symbol (end of line comments possible),
     *    - a data row starts with the name of the template variable and the
     *      value is separated by at least one whitespace.
     *
     *  In some cases the simple parser may fail, the characters scanned when
     *  splitting data rows are
     *
     *    ' ', \t
    */
    void fillDictionary() {

        namespace fs = boost::filesystem;

        fs::path pwd = fs::current_path();
        if (!fs::exists(varDictionary_) || fs::is_empty(varDictionary_))
            return;

        std::ifstream infile;
        infile.open(varDictionary_.c_str(), std::ifstream::in);

        char tmp[1024];
        unsigned int line_nr = 0;
        while(infile.good()) {
            std::fill_n(tmp, 1024, '\0');
            infile.getline(tmp, 1024);
            line_nr++;
            if(tmp[0] != '#') {
                std::string stmp(tmp);
                boost::trim(stmp);
                if(stmp.size() == 0)
                    continue;

                std::vector<std::string> all_strings;
                boost::split(all_strings, stmp,
                             boost::is_any_of("\r\n\v\f\t "),
                             boost::token_compress_on);

                if(all_strings.size() < 2) {
                    std::cout << "PROBLEM with the following line "
                              << "(at least name and value required)!"
                              << std::endl;
                    std::cout << stmp << std::endl;
                    std::ostringstream ex;
                    ex << "Invalid data file on line " << line_nr;
                    throw OptPilotException(
                            "GenerateOpalSimulation::fillDictionary()",
                            ex.str());
                }

                std::string varname = all_strings[0];
                std::string value   = all_strings[1];
                scale(varname, &value);
                dictionary_.insert(std::pair<std::string, std::string>(
                            varname, value));
            }
        }

        infile.close();
    }


    /// Helper method to scale variable if necessary.
    void scale(std::string name, std::string *value) {
        if(scaleVars_.count(name) > 0) {
            std::istringstream instr(*value);
            double val;
            instr >> val;
            val *= scaleVars_[name];

            // we keep values as strings, since we have to write them to
            // output streams anyways..
            std::ostringstream of;
            of << val;
            *value = of.str();
        }
    }

};

#endif