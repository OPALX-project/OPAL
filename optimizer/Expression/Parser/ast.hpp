//
// Namespace ast
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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
#ifndef AST_HPP
#define AST_HPP

#include <list>
#include <memory>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace client { namespace ast
{
    /// Local, C++20 replacement for boost::recursive_wrapper.
    /// Holds a heap-allocated T so that recursive/incomplete types can be
    /// stored as an alternative of a std::variant.
    template <typename T>
    class recursive_wrapper
    {
    public:
        recursive_wrapper() : ptr_(std::make_unique<T>()) {}
        recursive_wrapper(T const& value) : ptr_(std::make_unique<T>(value)) {}
        recursive_wrapper(T&& value) : ptr_(std::make_unique<T>(std::move(value))) {}

        recursive_wrapper(recursive_wrapper const& other) : ptr_(std::make_unique<T>(*other.ptr_)) {}
        recursive_wrapper(recursive_wrapper&&) noexcept = default;

        recursive_wrapper& operator=(recursive_wrapper const& other) {
            if (this != &other) {
                ptr_ = std::make_unique<T>(*other.ptr_);
            }
            return *this;
        }
        recursive_wrapper& operator=(recursive_wrapper&&) noexcept = default;

        recursive_wrapper& operator=(T const& value) { *ptr_ = value; return *this; }
        recursive_wrapper& operator=(T&& value) { *ptr_ = std::move(value); return *this; }

        T& get() { return *ptr_; }
        T const& get() const { return *ptr_; }

        explicit operator T&() { return *ptr_; }
        explicit operator T const&() const { return *ptr_; }

    private:
        std::unique_ptr<T> ptr_;
    };

    template <typename T>
    struct is_recursive_wrapper : std::false_type {};
    template <typename T>
    struct is_recursive_wrapper<recursive_wrapper<T>> : std::true_type {};
    template <typename T>
    inline constexpr bool is_recursive_wrapper_v = is_recursive_wrapper<T>::value;

    /// Local replacement for boost::apply_visitor: dispatches on a
    /// std::variant while transparently unwrapping recursive_wrapper<T>
    /// alternatives to T, matching boost::variant's visitation semantics.
    template <typename Visitor, typename Variant>
    decltype(auto) apply_visitor(Visitor&& vis, Variant&& var)
    {
        return std::visit(
            [&vis](auto&& alternative) -> decltype(auto) {
                using Alternative = std::decay_t<decltype(alternative)>;
                if constexpr (is_recursive_wrapper_v<Alternative>) {
                    return vis(alternative.get());
                } else {
                    return vis(alternative);
                }
            },
            std::forward<Variant>(var));
    }

    struct tagged
    {
        int id; // Used to annotate the AST with the iterator position.
                // This id is used as a key to a map<int, Iterator>
                // (not really part of the AST.)
    };

    struct nil {};
    struct unary;
    struct function_call;
    struct expression;

    struct identifier : tagged
    {
        identifier(std::string const& name = "") : name(name) {}
        std::string name;
    };

    struct quoted_string
    {
        quoted_string(std::string const& value = "") : value(value) {}
        std::string value;
    };

    typedef std::variant<
            nil
          , bool
          , unsigned int
          , double
          , identifier
          , recursive_wrapper<unary>
          , recursive_wrapper<function_call>
          , recursive_wrapper<expression>
        >
    operand;

    typedef std::variant<
            expression
          , quoted_string
        >
    function_call_argument;

    enum optoken
    {
        op_plus,
        op_minus,
        op_times,
        op_divide,
        op_positive,
        op_negative,
        op_not,
        op_equal,
        op_not_equal,
        op_less,
        op_less_equal,
        op_greater,
        op_greater_equal,
        op_and,
        op_or
    };

    struct unary
    {
        optoken operator_;
        operand operand_;
    };

    struct operation
    {
        optoken operator_;
        operand operand_;
    };

    struct function_call
    {
        identifier function_name;
        std::list<function_call_argument> args;
    };

    struct expression
    {
        operand first;
        std::list<operation> rest;
    };
#if 0
    // print functions for debugging
    inline std::ostream& operator<<(std::ostream& out, nil)
    {
        out << "nil"; return out;
    }
#endif
    inline std::ostream& operator<<(std::ostream& out, identifier const& id)
    {
        out << id.name; return out;
    }
}}

#endif