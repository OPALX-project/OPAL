//
// Namespace Matheval
//   Implements a recursive-descent parser for mathematical expressions,
//   building an AST that can be optimized and evaluated at runtime.
//
// Copyright (C) 2011-2013, Rhys Ulerich
// Copyright (C) ???, Martin Bauer
// Copyright (C) 2017, Henri Menke
// Copyright (c) 2026, PSI, Villigen, Switzerland
//
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
#pragma once

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

namespace matheval {

    namespace detail {

        namespace math {

            /** @brief Sign function
             *
             * Missing function in the STL. This calculates the mathematical sign
             * function.
             *
             * @f[
             *   \mathop{\mathrm{sgn}}(x) =
             *   \begin{cases}
             *      1 & x > 0 \\
             *      0 & x = 0 \\
             *     -1 & x < 0 \\
             *   \end{cases}
             * @f]
             *
             * @param[in] x number
             * @returns the sign of x
             */
            template < typename T >
            constexpr T sgn(T x) { return (T{0} < x) - (x < T{0}); }

            /** @brief isnan function with adjusted return type */
            template < typename T >
            T isnan(T x) { return std::isnan(x); }

            /** @brief isinf function with adjusted return type */
            template < typename T >
            T isinf(T x) { return std::isinf(x); }

            /** @brief Convert radians to degrees */
            template < typename T >
            inline constexpr T deg(T x) { return x * Units::rad2deg; }

            /** @brief Convert degrees to radians */
            template < typename T >
            inline constexpr T rad(T x) { return x * Units::deg2rad; }
        }

        // Forward declarations
        template < typename real_t > struct expr_ast;
        template < typename real_t > struct unary_op;
        template < typename real_t > struct binary_op;

        struct nil {};

        // Simple recursive wrapper for std::variant using unique_ptr
        template <typename T>
        struct recursive_wrapper {
            recursive_wrapper() : ptr(std::make_unique<T>()) {}
            recursive_wrapper(const T& t) : ptr(std::make_unique<T>(t)) {}
            recursive_wrapper(T&& t) : ptr(std::make_unique<T>(std::move(t))) {}
            recursive_wrapper(const recursive_wrapper& other) : ptr(std::make_unique<T>(*other.ptr)) {}
            recursive_wrapper(recursive_wrapper&& other) noexcept = default;
            
            recursive_wrapper& operator=(const recursive_wrapper& other) {
                if (this != &other) {
                    ptr = std::make_unique<T>(*other.ptr);
                }
                return *this;
            }

            recursive_wrapper& operator=(recursive_wrapper&& other) noexcept = default;

            T& get() { return *ptr; }
            const T& get() const { return *ptr; }

            explicit operator T&() { return *ptr; }
            explicit operator const T&() const { return *ptr; }

        private:
            std::unique_ptr<T> ptr;
        };

        // AST

        /** @brief Abstract Syntax Tree
         *
         * Stores the abstract syntax tree (AST) of the parsed mathematical
         * expression.
         */
        template < typename real_t >
        struct expr_ast
        {
            using tree_t = std::variant<
                nil // can't happen!
                , real_t
                , std::string
                , recursive_wrapper<expr_ast<real_t>>
                , recursive_wrapper<binary_op<real_t>>
                , recursive_wrapper<unary_op<real_t>>
                >;
        public:
            /** @brief AST storage
             *
             * The syntax tree can hold various types.  Numbers (`real_t`),
             * variables (`std::string`), the recursive tree itself
             * (`expr_ast`), binary operators (`binary_op`), and unary
             * operators (`unary_op`).
             */
            tree_t tree;

            /** @brief Default constructor
             *
             * Initializes the tree to a nil value to indicate inconsistent
             * state.
             */
            expr_ast() : tree(nil{}) {}

            /** @brief Copy constructor
             *
             * Deep copies the syntax tree.
             */
            template <typename Expr>
            expr_ast(Expr other) : tree(std::move(other)) {} // NOLINT

            /** @brief Add a tree */
            expr_ast& operator+=(expr_ast const& rhs);
            /** @brief subtract a tree */
            expr_ast& operator-=(expr_ast const& rhs);
            /** @brief Multiply by a tree */
            expr_ast& operator*=(expr_ast const& rhs);
            /** @brief Divide by a tree */
            expr_ast& operator/=(expr_ast const& rhs);
        };

        /** @brief Store a unary operator and its argument tree */
        template < typename real_t >
        struct unary_op
        {
            /** @brief Signature of a unary operator: op(x) */
            using op_t = std::function<real_t(real_t)>;

            /** @brief Save the operator and the argument tree */
            unary_op(op_t op, expr_ast<real_t> rhs)
                : op(std::move(op)),
                  rhs(std::make_unique<expr_ast<real_t>>(std::move(rhs)))
            {}

            unary_op(const unary_op& other)
                : op(other.op),
                  rhs(std::make_unique<expr_ast<real_t>>(*other.rhs))
            {}

            unary_op(unary_op&&) noexcept = default;

            unary_op& operator=(const unary_op& other) {
                if (this != &other) {
                    op = other.op;
                    rhs = std::make_unique<expr_ast<real_t>>(*other.rhs);
                }
                return *this;
            }

            unary_op& operator=(unary_op&&) noexcept = default;

            /** @brief Stored operator */
            op_t op;

            /** @brief Stored argument tree */
            std::unique_ptr<expr_ast<real_t>> rhs;
        };

        /** @brief Store a binary operator and its argument trees */
        template < typename real_t >
        struct binary_op
        {
            /** @brief Signature of a binary operator: op(x,y) */
            using op_t = std::function<real_t(real_t, real_t)>;

            /** @brief Save the operator and the argument trees */
            binary_op(op_t op, expr_ast<real_t> lhs, expr_ast<real_t> rhs)
                : op(std::move(op)),
                  lhs(std::make_unique<expr_ast<real_t>>(std::move(lhs))),
                  rhs(std::make_unique<expr_ast<real_t>>(std::move(rhs)))
            {}

            binary_op(const binary_op& other)
                : op(other.op), 
                  lhs(std::make_unique<expr_ast<real_t>>(*other.lhs)),
                  rhs(std::make_unique<expr_ast<real_t>>(*other.rhs))
            {}

            binary_op(binary_op&&) noexcept = default;

            binary_op& operator=(const binary_op& other) {
                if (this != &other) {
                    op = other.op;
                    lhs = std::make_unique<expr_ast<real_t>>(*other.lhs);
                    rhs = std::make_unique<expr_ast<real_t>>(*other.rhs);
                }
                return *this;
            }

            binary_op& operator=(binary_op&&) noexcept = default;

            /** @brief Stored operator */
            op_t op;

            /** @brief Stored argument tree of first argument */
            std::unique_ptr<expr_ast<real_t>> lhs;

            /** @brief Stored argument tree of second argument */
            std::unique_ptr<expr_ast<real_t>> rhs;
        };

        template < typename real_t >
        expr_ast<real_t>& expr_ast<real_t>::operator+=(expr_ast<real_t> const& rhs)
        {
            tree = binary_op<real_t>(std::plus<real_t>{}, tree, rhs);
            return *this;
        }
        template < typename real_t >
        expr_ast<real_t>& expr_ast<real_t>::operator-=(expr_ast<real_t> const& rhs)
        {
            tree = binary_op<real_t>(std::minus<real_t>{}, tree, rhs);
            return *this;
        }
        template < typename real_t >
        expr_ast<real_t>& expr_ast<real_t>::operator*=(expr_ast<real_t> const& rhs)
        {
            tree = binary_op<real_t>(std::multiplies<real_t>{}, tree, rhs);
            return *this;
        }
        template < typename real_t >
        expr_ast<real_t>& expr_ast<real_t>::operator/=(expr_ast<real_t> const& rhs)
        {
            tree = binary_op<real_t>(std::divides<real_t>{}, tree, rhs);
            return *this;
        }

        /** @brief Evaluate the Abstract Syntax Tree
         *
         * This visits all the variants of the AST and applies the stored
         * operators.
         */
        template < typename real_t >
        class eval_ast
        {
        public:
            /** @brief Necessary typedef for visitor */
            using result_type = real_t;

            /** @brief Type of the symbol table */
            using symbol_table_t = std::unordered_map<std::string, result_type>;

            /** @brief Constructor
             *
             * Saves the symbol table to apply variables.
             */
            explicit eval_ast(symbol_table_t sym) : st(std::move(sym)) {}

            /** @brief Empty nodes in the tree evaluate to 0 */
            result_type operator()(nil /*unused*/) const { return 0; }

            /** @brief Numbers evaluate to themselves */
            result_type operator()(result_type n)  const { return n; }

            /** @brief Variables evaluate to their value in the symbol table */
            result_type operator()(std::string const &c) const
            {
                auto it = st.find(c);
                if (it == st.end()) {
                    throw std::invalid_argument("Unknown variable " + c); // NOLINT
                }
                return it->second;
            }

            /** @brief Recursively evaluate the AST */
            result_type operator()(expr_ast<real_t> const& ast) const
            {
                return std::visit(*this, ast.tree);
            }

            /** @brief Evaluate a binary operator and optionally recurse its operands */
            result_type operator()(binary_op<real_t> const& tree) const
            {
                return tree.op(std::visit(*this, tree.lhs->tree),
                               std::visit(*this, tree.rhs->tree));
            }

            /** @brief Evaluate a unary operator and optionally recurse its operand */
            result_type operator()(unary_op<real_t> const& tree) const
            {
                return tree.op(std::visit(*this, tree.rhs->tree));
            }

            /** @brief Handle recursive_wrapper for expr_ast */
            result_type operator()(recursive_wrapper<expr_ast<real_t>> const& w) const
            {
                return std::visit(*this, w.get().tree);
            }

            /** @brief Handle recursive_wrapper for binary_op */
            result_type operator()(recursive_wrapper<binary_op<real_t>> const& w) const
            {
                return (*this)(w.get());
            }

            /** @brief Handle recursive_wrapper for unary_op */
            result_type operator()(recursive_wrapper<unary_op<real_t>> const& w) const
            {
                return (*this)(w.get());
            }

        private:
            symbol_table_t st;
        };

        template <typename real_t> struct ConstantFolder {
            /** @brief Necessary typedef for visitor */
            using result_type = typename expr_ast<real_t>::tree_t;

            /** @brief Empty nodes in the tree evaluate to 0 */
            result_type operator()(nil /*unused*/) const { return 0; }

            /** @brief Numbers evaluate to themselves */
            result_type operator()(real_t n) const { return n; }

            /** @brief Variables do not evaluate */
            result_type operator()(std::string const& c) const { return c; }

            /** @brief Recursively evaluate the AST */
            result_type operator()(expr_ast<real_t> const& ast) const {
                return std::visit(*this, ast.tree);
            }

            /** @brief Evaluate a binary operator and optionally recurse its operands */
            result_type operator()(binary_op<real_t> const& tree) const {
                auto lhs = std::visit(*this, tree.lhs->tree);
                auto rhs = std::visit(*this, tree.rhs->tree);

                /* If both operands are known, we can directly evaluate the function,
                 * else we just update the children with the new expressions. */
                if (auto l = std::get_if<real_t>(&lhs)) {
                    if (auto r = std::get_if<real_t>(&rhs)) {
                        return tree.op(*l, *r);
                    }
                }
                return binary_op<real_t>(tree.op, expr_ast<real_t>(lhs), expr_ast<real_t>(rhs));
            }

            /** @brief Evaluate a unary operator and optionally recurse its operand */
            result_type operator()(unary_op<real_t> const& tree) const {
                auto rhs = std::visit(*this, tree.rhs->tree);
                /* If the operand is known, we can directly evaluate the function. */
                if (auto r = std::get_if<real_t>(&rhs)) {
                    return tree.op(*r);
                }
                return unary_op<real_t>(tree.op, expr_ast<real_t>(rhs));
            }

            /** @brief Handle recursive_wrapper for expr_ast */
            result_type operator()(recursive_wrapper<expr_ast<real_t>> const& w) const {
                return std::visit(*this, w.get().tree);
            }

            /** @brief Handle recursive_wrapper for binary_op */
            result_type operator()(recursive_wrapper<binary_op<real_t>> const& w) const {
                return (*this)(w.get());
            }

            /** @brief Handle recursive_wrapper for unary_op */
            result_type operator()(recursive_wrapper<unary_op<real_t>> const& w) const {
                return (*this)(w.get());
            }
        };

        // Simple recursive descent parser
        template <typename real_t>
        class SimpleMathParser {
        public:
            SimpleMathParser(std::string_view input) 
                : input_(input), pos_(0) {}

            expr_ast<real_t> parse() {
                pos_ = 0;
                skipWhitespace();
                auto result = parseExpression();
                skipWhitespace();
                if (pos_ < input_.length()) {
                    throw std::runtime_error("Unexpected characters at end of expression");
                }
                return result;
            }

        private:
            std::string_view input_;
            size_t pos_;

            void skipWhitespace() {
                while (pos_ < input_.length() && std::isspace(static_cast<unsigned char>(input_[pos_]))) {
                    pos_++;
                }
            }

            bool match(char c) {
                skipWhitespace();
                if (pos_ < input_.length() && input_[pos_] == c) {
                    pos_++;
                    return true;
                }
                return false;
            }

            bool match(const std::string& str) {
                skipWhitespace();
                if (input_.compare(pos_, str.length(), str) == 0) {
                    pos_ += str.length();
                    return true;
                }
                return false;
            }

            void expect(char c) {
                if (!match(c)) {
                    throw std::runtime_error(std::string("Expected '") + c + "'");
                }
            }

            expr_ast<real_t> parseExpression() {
                auto result = parseTerm();
                while (true) {
                    if (match('+')) {
                        result += parseTerm();
                    } else if (match('-')) {
                        result -= parseTerm();
                    } else {
                        break;
                    }
                }
                return result;
            }

            expr_ast<real_t> parseTerm() {
                auto result = parseFactor();
                while (true) {
                    if (match('*')) {
                        result *= parseFactor();
                    } else if (match('/')) {
                        result = expr_ast<real_t>(binary_op<real_t>(std::divides<real_t>{}, result, parseFactor()));
                    } else if (match('%')) {
                        auto fmod = static_cast<real_t(*)(real_t,real_t)>(&std::fmod);
                        result = expr_ast<real_t>(binary_op<real_t>(fmod, result, parseFactor()));
                    } else {
                        break;
                    }
                }
                return result;
            }

            expr_ast<real_t> parseFactor() {
                auto lhs = parsePrimary();
                skipWhitespace();
                if (match("**")) {
                    auto rhs = parseFactor();
                    auto pow = static_cast<real_t(*)(real_t,real_t)>(&std::pow);
                    lhs = expr_ast<real_t>(binary_op<real_t>(pow, lhs, rhs));
                }
                return lhs;
            }

            expr_ast<real_t> parsePrimary() {
                skipWhitespace();

                // Number
                if (std::isdigit(static_cast<unsigned char>(input_[pos_])) || 
                    (pos_ < input_.length() && input_[pos_] == '.')) {
                    return parseNumber();
                }

                // Parentheses
                if (match('(')) {
                    auto result = parseExpression();
                    expect(')');
                    return result;
                }

                // Unary operators
                if (match('-')) {
                    return expr_ast<real_t>(unary_op<real_t>(std::negate<real_t>{}, parsePrimary()));
                }
                if (match('+')) {
                    return parsePrimary();
                }

                // Functions
                std::string identifier = parseIdentifier();
                if (identifier.empty()) {
                    throw std::runtime_error("Unexpected character");
                }

                // Check for constants
                for (const auto& c : getConstants()) {
                    if (identifier == c.first) {
                        return expr_ast<real_t>(c.second);
                    }
                }

                // Check for binary functions
                for (const auto& f : getBinaryFunctions()) {
                    if (identifier == f.first) {
                        expect('(');
                        auto arg1 = parseExpression();
                        expect(',');
                        auto arg2 = parseExpression();
                        expect(')');
                        return expr_ast<real_t>(binary_op<real_t>(f.second, arg1, arg2));
                    }
                }

                // Check for unary functions
                for (const auto& f : getUnaryFunctions()) {
                    if (identifier == f.first) {
                        expect('(');
                        auto arg = parseExpression();
                        expect(')');
                        return expr_ast<real_t>(unary_op<real_t>(f.second, arg));
                    }
                }

                // Variable
                return expr_ast<real_t>(identifier);
            }

            expr_ast<real_t> parseNumber() {
                skipWhitespace();
                size_t start = pos_;
                bool has_dot = false;
                bool has_exp = false;
                
                if (pos_ < input_.length() && (input_[pos_] == '+' || input_[pos_] == '-')) {
                    pos_++;
                }

                while (pos_ < input_.length()) {
                    char c = input_[pos_];
                    if (std::isdigit(static_cast<unsigned char>(c))) {
                        pos_++;
                    } else if (c == '.' && !has_dot) {
                        has_dot = true;
                        pos_++;
                    } else if ((c == 'e' || c == 'E') && !has_exp) {
                        has_exp = true;
                        pos_++;
                        if (pos_ >= input_.length() || 
                            (!std::isdigit(static_cast<unsigned char>(input_[pos_])) && 
                            input_[pos_] != '+' && input_[pos_] != '-')) {
                            throw std::runtime_error("Expected exponent digits at position " + std::to_string(pos_));
                        }
                        if (input_[pos_] == '+' || input_[pos_] == '-') pos_++;
                        if (pos_ >= input_.length() || !std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                            throw std::runtime_error("Expected exponent digits at position " + std::to_string(pos_));
                        }
                    } else {
                        break;
                    }
                }

                if (pos_ > start) {
                    try {
                        std::string_view sv(input_.data() + start, pos_ - start);
                        real_t value = static_cast<real_t>(std::stod(std::string(sv)));
                        return expr_ast<real_t>(value);
                    } catch (...) {
                        throw std::runtime_error("Invalid number");
                    }
                }
                throw std::runtime_error("Expected number");
            }

            std::string parseIdentifier() {
                skipWhitespace();
                size_t start = pos_;
                if (pos_ < input_.length() && (std::isalpha(static_cast<unsigned char>(input_[pos_])) || input_[pos_] == '_')) {
                    pos_++;
                    while (pos_ < input_.length() && 
                           (std::isalnum(static_cast<unsigned char>(input_[pos_])) || input_[pos_] == '_')) {
                        pos_++;
                    }
                    return std::string(input_.substr(start, pos_ - start));
                }
                return "";
            }

            static const std::array<std::pair<std::string_view, real_t>, 4>& getConstants() {
                static const std::array<std::pair<std::string_view, real_t>, 4> constants = {{
                    {"e",       Physics::e},
                    {"epsilon", std::numeric_limits<real_t>::epsilon()},
                    {"phi",     (1.0 + std::sqrt(5.0)) / 2.0},
                    {"pi",      Physics::pi}
                }};
                return constants;
            }

            static const std::array<std::pair<std::string_view, typename unary_op<real_t>::op_t>, 31>& getUnaryFunctions() {
                static const std::array<std::pair<std::string_view, typename unary_op<real_t>::op_t>, 31> funcs = {{
                    {"abs",     static_cast<real_t(*)(real_t)>(&std::abs)},
                    {"acos",    static_cast<real_t(*)(real_t)>(&std::acos)},
                    {"acosh",   static_cast<real_t(*)(real_t)>(&std::acosh)},
                    {"asin",    static_cast<real_t(*)(real_t)>(&std::asin)},
                    {"asinh",   static_cast<real_t(*)(real_t)>(&std::asinh)},
                    {"atan",    static_cast<real_t(*)(real_t)>(&std::atan)},
                    {"atanh",   static_cast<real_t(*)(real_t)>(&std::atanh)},
                    {"cbrt",    static_cast<real_t(*)(real_t)>(&std::cbrt)},
                    {"ceil",    static_cast<real_t(*)(real_t)>(&std::ceil)},
                    {"cos",     static_cast<real_t(*)(real_t)>(&std::cos)},
                    {"cosh",    static_cast<real_t(*)(real_t)>(&std::cosh)},
                    {"deg2rad", static_cast<real_t(*)(real_t)>(&math::deg)},
                    {"erf",     static_cast<real_t(*)(real_t)>(&std::erf)},
                    {"erfc",    static_cast<real_t(*)(real_t)>(&std::erfc)},
                    {"exp",     static_cast<real_t(*)(real_t)>(&std::exp)},
                    {"exp2",    static_cast<real_t(*)(real_t)>(&std::exp2)},
                    {"floor",   static_cast<real_t(*)(real_t)>(&std::floor)},
                    {"isinf",   static_cast<real_t(*)(real_t)>(&math::isinf)},
                    {"isnan",   static_cast<real_t(*)(real_t)>(&math::isnan)},
                    {"log",     static_cast<real_t(*)(real_t)>(&std::log)},
                    {"log2",    static_cast<real_t(*)(real_t)>(&std::log2)},
                    {"log10",   static_cast<real_t(*)(real_t)>(&std::log10)},
                    {"rad2deg", static_cast<real_t(*)(real_t)>(&math::rad)},
                    {"round",   static_cast<real_t(*)(real_t)>(&std::round)},
                    {"sgn",     static_cast<real_t(*)(real_t)>(&math::sgn)},
                    {"sin",     static_cast<real_t(*)(real_t)>(&std::sin)},
                    {"sinh",    static_cast<real_t(*)(real_t)>(&std::sinh)},
                    {"sqrt",    static_cast<real_t(*)(real_t)>(&std::sqrt)},
                    {"tan",     static_cast<real_t(*)(real_t)>(&std::tan)},
                    {"tanh",    static_cast<real_t(*)(real_t)>(&std::tanh)},
                    {"tgamma",  static_cast<real_t(*)(real_t)>(&std::tgamma)}
                }};
                return funcs;
            }

            static const std::array<std::pair<std::string_view, typename binary_op<real_t>::op_t>, 4>& getBinaryFunctions() {
                static const std::array<std::pair<std::string_view, typename binary_op<real_t>::op_t>, 4> funcs = {{
                    {"atan2", static_cast<real_t(*)(real_t,real_t)>(&std::atan2)},
                    {"max",   static_cast<real_t(*)(real_t,real_t)>(&std::fmax)},
                    {"min",   static_cast<real_t(*)(real_t,real_t)>(&std::fmin)},
                    {"pow",   static_cast<real_t(*)(real_t,real_t)>(&std::pow)}
                }};
                return funcs;
            }
        };

        /** @brief Parse an expression
         *
         * This function builds the parser and parses the string into
         * an AST.
         *
         * @param[in] first iterator to the start of the input sequence
         * @param[in] last  iterator to the end of the input sequence
         */
        template <typename real_t, typename Iterator>
        expr_ast<real_t> parse(Iterator first, Iterator last) {
            std::string input(first, last);
            SimpleMathParser<real_t> parser(input);
            return parser.parse();
        }

    } // namespace detail


    /** @brief Class interface
     *
     * This class hides the parser, AST, and AST traversal behind some
     * member functions.
     *
     * @tparam real_t datatype of the result
     */
    template < typename real_t >
    class Parser
    {
        detail::expr_ast<real_t> ast;
    public:
        /** @brief Parse an expression
         *
         * This function builds the parser and parses the iterator into
         * an AST.
         *
         * @param[in] first iterator to the start of the input sequence
         * @param[in] last  iterator to the end of the input sequence
         */
        template < typename Iterator >
        void parse(Iterator first, Iterator last)
        {
            ast = detail::parse<real_t>(first, last);
        }

        /** @overload parse(Iterator first, Iterator last) */
        void parse(std::string const &str)
        {
            parse(str.begin(), str.end());
        }

        void optimize() {
            ast.tree = std::visit(detail::ConstantFolder<real_t>{}, ast.tree);
        }

        /** @brief Evaluate the AST with a given symbol table
         *
         * @param[in] st the symbol table for variables
         */
        real_t evaluate(typename detail::eval_ast<real_t>::symbol_table_t const& st) const
        {
            detail::eval_ast<real_t> solver(st);
            return solver(ast);
        }
    };


    /** @brief Convenience function
     *
     * This function builds the parser, parses the iterator to an AST,
     * evaluates it, and returns the result.
     *
     * @param[in] first iterator to the start of the input sequence
     * @param[in] last  iterator to the end of the input sequence
     * @param[in] st    the symbol table for variables
     */
    template < typename real_t, typename Iterator >
    real_t parse(Iterator first, Iterator last,
                 typename detail::eval_ast<real_t>::symbol_table_t const& st)
    {
        Parser<real_t> parser;
        parser.parse(first, last);
        return parser.evaluate(st);
    }

    /** @overload parse(Iterator first, Iterator last, typename detail::eval_ast<real_t>::symbol_table_t const &st) */
    template < typename real_t >
    real_t parse(std::string const& str,
                 typename detail::eval_ast<real_t>::symbol_table_t const& st)
    {
        return parse<real_t>(str.begin(), str.end(), st);
    }

} // namespace matheval