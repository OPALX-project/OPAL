//
// Class SDDSParser
//   This class writes column entries of SDDS files.
//
// Copyright (c) 2015, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin
// All rights reserved
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
#ifndef __SDDSPARSER_H__
#define __SDDSPARSER_H__

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/file.hpp"
#include "Util/SDDSParser/SDDSParserException.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <variant>

namespace SDDS {

    class SDDSParser {
    private:
        std::string readFile();
        static std::string normalizeKey(std::string_view value);
        std::string sddsFileName_m;

        /// mapping from parameter name to offset in params_m
        std::unordered_map<std::string, std::size_t> paramNameToID_m;
        /// mapping from column name to ID in columns_m
        std::unordered_map<std::string, std::size_t> columnNameToID_m;

        SDDS::file sddsData_m;

    public:
        SDDSParser();
        SDDSParser(const std::string& input);
        void setInput(const std::string& input);
        file run();

        const file& getData() const;
        const ast::columnData_t& getColumnData(const std::string& columnName) const;

        ast::dataType getColumnType(const std::string& col_name) const {
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
        void getValue(int t, const std::string& column_name, T& nval) {

            int col_idx = getColumnIndex(column_name);

            // round timestep to last if not in range
            std::size_t row_idx = 0;
            std::size_t num_rows = sddsData_m.sddsColumns_m[col_idx].values_m.size();
            if (num_rows == 0) {
                throw SDDSParserException("SDDSParser::getValue",
                                          "requested column has no parsed rows");
            }
            if (t <= 0 || static_cast<std::size_t>(t) > num_rows) {
                row_idx = num_rows - 1;
            } else {
                row_idx = static_cast<std::size_t>(t) - 1;
            }

            const ast::variant_t& val = sddsData_m.sddsColumns_m[col_idx].values_m[row_idx];
            nval = getVariantValue<T>(val, getColumnType(column_name));
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
        void getInterpolatedValue(const std::string& ref_name,
                                  double ref_val,
                                  const std::string& col_name,
                                  T& nval) {
            T value_before = 0;
            T value_after  = 0;
            double value_before_ref = 0;
            double value_after_ref  = 0;

            std::size_t col_idx_ref = getColumnIndex(ref_name);
            ast::columnData_t &ref_values = sddsData_m.sddsColumns_m[col_idx_ref].values_m;
            int index = getColumnIndex(col_name);
            ast::columnData_t &col_values = sddsData_m.sddsColumns_m[index].values_m;

            std::size_t this_row = 0;
            std::size_t num_rows = ref_values.size();
            ast::dataType dataType = getColumnType(col_name);
            for (this_row = 0; this_row < num_rows; this_row++) {
                value_after_ref = std::get<double>(ref_values[this_row]);

                if (ref_val < value_after_ref) {

                    std::size_t prev_row = 0;
                    if (this_row > 0) prev_row = this_row - 1;

                    value_before = getVariantValue<T>(col_values[prev_row], dataType);
                    value_after  = getVariantValue<T>(col_values[this_row], dataType);

                    value_before_ref = std::get<double>(ref_values[prev_row]);
                    value_after_ref  = std::get<double>(ref_values[this_row]);

                    break;
                }
            }

            if (this_row == num_rows) {
                throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                          "all values < specified reference value");
            }

            // simple linear interpolation
            if (ref_val - value_before_ref < 1e-8) {
                nval = value_before;
            } else {
                nval = value_before + (ref_val - value_before_ref)
                    * (value_after - value_before)
                    / (value_after_ref - value_before_ref);
            }

            if (!std::isfinite(nval)) {
                throw SDDSParserException("SDDSParser::getInterpolatedValue",
                                          "Interpolated value either NaN or Inf.");
            }

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
        void getInterpolatedValue(double spos, const std::string& col_name, T& nval) {
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
        void getParameterValue(const std::string& parameter_name, T& nval) {
            auto key = normalizeKey(parameter_name);
            auto it = paramNameToID_m.find(key);
            if (it != paramNameToID_m.end()) {
                const auto& value = sddsData_m.sddsParameters_m[it->second].value_m;
                nval = std::get<T>(value);
            } else {
                throw SDDSParserException("SDDSParser::getParameterValue",
                                        "unknown parameter name: '" + key + "'!");
            }
        }

        /// Convert value from variant (only numeric types) to a value of type T
        template <typename T>
        T getVariantValue(const ast::variant_t& val, ast::dataType type) const {
            static_assert(
                std::is_same_v<T, float>  ||
                std::is_same_v<T, double> ||
                std::is_same_v<T, short>  ||
                std::is_same_v<T, long>,
                "Unsupported type T in SDDSParser::getVariantValue"
            );

            try {
                switch (type) {
                    case ast::dataType::FLOAT: {
                        return static_cast<T>(std::get<float>(val));
                    }
                    case ast::dataType::DOUBLE: {
                        return static_cast<T>(std::get<double>(val));
                    }
                    case ast::dataType::SHORT: {
                        return static_cast<T>(std::get<short>(val));
                    }
                    case ast::dataType::LONG: {
                        return static_cast<T>(std::get<long>(val));
                    }
                    default: {
                        throw SDDSParserException("SDDSParser::getVariantValue",
                                                "unsupported ast::dataType"
                        );
                    }
                }
            }
            catch (const std::bad_variant_access&) {
                throw SDDSParserException("SDDSParser::getVariantValue",
                                          "variant and dataType mismatch"
                );
            }
        }

    private:
        int getColumnIndex(const std::string& col_name) const;
    };

    inline
    const file& SDDSParser::getData() const {
        return sddsData_m;
    }
}

#endif