//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/file.hpp"
#include "Util/SDDSParser/skipper.hpp"
#include "Util/SDDSParser/array.hpp"
#include "Util/SDDSParser/associate.hpp"
#include "Util/SDDSParser/column.hpp"
#include "Util/SDDSParser/data.hpp"
#include "Util/SDDSParser/description.hpp"
#include "Util/SDDSParser/include.hpp"
#include "Util/SDDSParser/parameter.hpp"

#include "Util/SDDSParser/SDDSParserException.h"


#include "boost/optional/optional_io.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#ifndef __SDDSPARSER_H__
#define __SDDSPARSER_H__

namespace SDDS {

    class SDDSParser {
    private:
        std::string readFile();
        static void fixCaseSensitivity(std::string &for_string);
        static std::string fixCaseSensitivity(const std::string &for_string) {
            std::string retval(for_string);
            fixCaseSensitivity(retval);
            return retval;
        }
        std::string sddsFileName_m;

        /// mapping from parameter name to offset in params_m
        std::map<std::string, int> paramNameToID_m;
        /// mapping from column name to ID in columns_m
        std::map<std::string, int> columnNameToID_m;

        SDDS::file sddsData_m;

    public:
        SDDSParser();
        SDDSParser(const std::string &input);
        void setInput(const std::string &input);
        file run();

        file getData();
        ast::columnData_t getColumnData(const std::string &columnName);

        ast::datatype getColumnType(const std::string &col_name) {
            int index = getColumnIndex(col_name);
            return *sddsData_m.sddsColumns_m[index].type_m;
        }

        /**
         *  Converts the string value of a parameter at timestep t to a value of
         *  type T.
         *
         *  @param t timestep (beginning at 1, -1 means last)
         *  @param column_name parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getValue(int t, std::string column_name, T& nval) {

            fixCaseSensitivity(column_name);

            int col_idx = getColumnIndex(column_name);

            // round timestep to last if not in range
            size_t row_idx = 0;
            size_t num_rows = sddsData_m.sddsColumns_m[col_idx].values_m.size();
            if(t <= 0 || static_cast<size_t>(t) > num_rows)
                row_idx = num_rows - 1;
            else
                row_idx = static_cast<size_t>(t) - 1;

            ast::variant_t val = sddsData_m.sddsColumns_m[col_idx].values_m[row_idx];
            nval = getBoostVariantValue<T>(val, (int)getColumnType(column_name));
        }


        /**
         *  Converts the string value of a parameter to a value
         *  of type T.
         *
         *  @param ref_name reference quantity (e.g. spos)
         *  @param ref_val interpolate value of reference quantity (e.g. spos)
         *  @param col_name parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getInterpolatedValue(std::string ref_name, double ref_val,
                                  std::string col_name, T& nval) {
            T value_before = 0;
            T value_after  = 0;
            double value_before_ref = 0;
            double value_after_ref  = 0;


            size_t col_idx_ref = getColumnIndex(ref_name);
            ast::columnData_t &ref_values = sddsData_m.sddsColumns_m[col_idx_ref].values_m;
            int index = getColumnIndex(col_name);
            ast::columnData_t &col_values = sddsData_m.sddsColumns_m[index].values_m;

            size_t this_row = 0;
            size_t num_rows = ref_values.size();
            int datatype = (int)getColumnType(col_name);
            for(this_row = 0; this_row < num_rows; this_row++) {
                value_after_ref = boost::get<double>(ref_values[this_row]);

                if(ref_val < value_after_ref) {

                    size_t prev_row = 0;
                    if(this_row > 0) prev_row = this_row - 1;

                    value_before = getBoostVariantValue<T>(col_values[prev_row], datatype);
                    value_after  = getBoostVariantValue<T>(col_values[this_row], datatype);

                    value_before_ref = boost::get<double>(ref_values[prev_row]);
                    value_after_ref  = boost::get<double>(ref_values[this_row]);

                    break;
                }
            }

            if(this_row == num_rows)
                throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                          "all values < specified reference value");

            // simple linear interpolation
            if(ref_val - value_before_ref < 1e-8)
                nval = value_before;
            else
                nval = value_before + (ref_val - value_before_ref)
                    * (value_after - value_before)
                    / (value_after_ref - value_before_ref);

            if (!std::isfinite(nval))
                throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                          "Interpolated value either NaN or Inf.");
        }

        /**
         *  Converts the string value of a parameter at a position spos to a value
         *  of type T.
         *
         *  @param spos interpolate value at spos
         *  @param col_name parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getInterpolatedValue(double spos, std::string col_name, T& nval) {
            getInterpolatedValue("s", spos, col_name, nval);
        }

        /**
         *  Converts the string value of a parameter to a value
         *  of type T.
         *
         *  @param parameter_name parameter name
         *  @param nval store result of type T in nval
         */
        template <typename T>
        void getParameterValue(std::string parameter_name, T& nval) {
            fixCaseSensitivity(parameter_name);

            if (paramNameToID_m.count(parameter_name) > 0) {
                size_t id = paramNameToID_m[parameter_name];
                auto value = sddsData_m.sddsParameters_m[id].value_m;
                nval = boost::get<T>(value);
            } else {
                throw SDDSParserException("SDDSParser::getParameterValue",
                                        "unknown parameter name: '" + parameter_name + "'!");
            }
        }

        /// Convert value from boost variant (only numeric types) to a value of type T
        // use integer instead of ast::datatype enum since otherwise boost has ambigious overloads
        // as tested on 8-1-2019, boost 1.68, gcc 7.3
        template <typename T>
            T getBoostVariantValue(const ast::variant_t& val, int datatype) const {
            T value;
            try {
                switch (datatype) {
                case ast::FLOAT:
                    value = boost::get<float>(val);
                    break;
                case ast::DOUBLE:
                    value = boost::get<double>(val);
                    break;
                case ast::SHORT:
                    value = boost::get<short>(val);
                    break;
                case ast::LONG:
                    value = boost::get<long>(val);
                    break;
                default:
                    throw SDDSParserException("SDDSParser::getBoostVariantValue",
                                              "can't convert value to type T");
                }
            }
            catch (...) {
                throw SDDSParserException("SDDSParser::getBoostVariantValue",
                                          "can't convert value");
            }
            return value;
        }

    private:

        int getColumnIndex(std::string col_name) const;
    };

    inline
    file SDDSParser::getData() {
        return sddsData_m;
    }
}

#endif