#ifndef __SDDSVARIABLE_H__
#define __SDDSVARIABLE_H__

#include <string>

#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/SDDSReader.h"
#include "Util/SDDSParser/SDDSParserException.h"
#include "Expression/Parser/function.hpp"


/**
 *  A simple expression to get SDDS (filename) value near a
 *  specific position (ref_val, default: spos) of a reference
 *  variable (ref_name) for a variable (var_name). Possible
 *  argument orders:
 *      args = [var_name, ref_val, filename]
 *      args = [var_name, ref_name, ref_val, filename]
 */
struct SDDSVariable {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        switch ( args.size() ) {
            case 3: {
                var_name_      = boost::get<std::string>(args[0]);
                ref_name_      = "s";
                ref_val_       = boost::get<double>(args[1]);
                stat_filename_ = boost::get<std::string>(args[2]);
                break;
            }
            case 4: {
                var_name_      = boost::get<std::string>(args[0]);
                ref_name_      = boost::get<std::string>(args[1]);
                ref_val_       = boost::get<double>(args[2]);
                stat_filename_ = boost::get<std::string>(args[3]);
                break;
            }
            default: {
                throw OptPilotException("SDDSVariable::operator()",
                                        "sddsVariableAt expects 3 or 4 arguments, " +
                                        std::to_string(args.size()) + " given");
            }
        }

        bool is_valid = true;

        boost::scoped_ptr<SDDSReader> sim_stats(new SDDSReader(stat_filename_));
        try {
            sim_stats->parseFile();
        } catch (SDDSParserException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sim_value = 0.0;
        try {
            sim_stats->getInterpolatedValue(ref_name_, ref_val_, var_name_, sim_value);
        } catch(SDDSParserException &e) {
            std::cout << "Exception while getting value "
                      << "from SDDS file: " << e.what()
                      << std::endl;
            is_valid = false;
        } catch(...) {
            std::cout << "Exception while getting '" + var_name_ + "' "
                      << "from SDDS file. "
                      << std::endl;
            is_valid = false;
        }

        return boost::make_tuple(sim_value, is_valid);
    }

private:

    std::string var_name_;
    std::string stat_filename_;
    std::string ref_name_;
    double ref_val_;
};

/**
 *  A simple expression to get value from stat file near a
 *  specific position (ref_val, default: spos) of a reference
 *  variable (ref_var) for a variable (var_name). Possible
 *  argument orders:
 *      args = [var_name, ref_val]
 *      args = [var_name, ref_name, ref_val]
 */

struct sameSDDSVariable {
    sameSDDSVariable(const std::string & base_filename) {
        size_t pos = base_filename.find_last_of("/");
        std::string tmplfile = base_filename;
        if(pos != std::string::npos)
            tmplfile = base_filename.substr(pos+1);
        pos = tmplfile.find_last_of(".");
        // std::string simName =
        stat_filename_ = tmplfile.substr(0,pos) + ".stat";
    }

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() < 2 || args.size() > 3) {
            throw OptPilotException("sameSDDSVariable::operator()",
                                    "statVariableAt expects 2 or 3 arguments, " +
                                    std::to_string(args.size()) + " given");
        }

        args.push_back(stat_filename_);

        return var_(args);
    }

private:
    client::function::argument_t stat_filename_;
    SDDSVariable var_;
};

#endif
