/*=============================================================================
    Adapted from boost spirit mini_c example.

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(REQUIREMENTS_HPP)
#define REQUIREMENTS_HPP

#include "ast.hpp"
#include "error_handler.hpp"
#include <set>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace client { namespace code_gen
{
    struct requirements
    {
        typedef bool result_type;

        template <typename ErrorHandler>
        requirements(ErrorHandler& error_handler_)
        {
            namespace phx = boost::phoenix;
            using boost::phoenix::function;

            error_handler = function<ErrorHandler>(error_handler_)
                (std::string("Error! "),
                 phx::arg_names::_2,
                 phx::cref(error_handler_.iters)[phx::arg_names::_1]);
        }

        bool operator()(ast::nil) { BOOST_ASSERT(0); return false; }
        bool operator()(unsigned int x) { return true; }
        bool operator()(double x) { return true; }
        bool operator()(bool x) { return true; }
        bool operator()(ast::quoted_string const &x) { return true; }

        bool operator()(ast::operation const& x) {
            if (!boost::apply_visitor(*this, x.operand_))
                return false;
            return true;
        }

        bool operator()(ast::unary const& x) {
            if (!boost::apply_visitor(*this, x.operand_))
                return false;
            return true;
        }

        bool operator()(ast::identifier const& x) {
            variables_.insert(x.name);

            return true;
        }

        bool operator()(ast::function_call const& x) {
            functions_.insert(x.function_name.name);

            for(ast::function_call_argument const& arg: x.args) {
                if (!boost::apply_visitor(*this, arg))
                    return false;
                //if (!(*this)(arg))
                    //return false;
            }
            return true;
        }

        bool operator()(ast::expression const& x) {

            if (!boost::apply_visitor(*this, x.first))
                return false;

            for(ast::operation const& oper: x.rest) {
                if (!(*this)(oper))
                    return false;
        }

        return true;
    }

    std::set<std::string> variables() { return variables_; }
    std::set<std::string> functions() { return functions_; }

    private:

        boost::function<
            void(int tag, std::string const& what)>
        error_handler;

        std::set<std::string> variables_;
        std::set<std::string> functions_;
    };
}}

#endif