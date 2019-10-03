/*=============================================================================
    Adapted from boost spirit mini_c example.

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(STACKEVALUATOR_HPP)
#define STACKEVALUATOR_HPP

#include <map>
#include <vector>

#include <boost/function.hpp>
#include <boost/assert.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/variant.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "ast.hpp"
#include "function.hpp"
#include "error_handler.hpp"

namespace client { namespace code_gen
{
    struct StackEvaluator {

        typedef bool result_type;

        template <typename ErrorHandler>
        StackEvaluator(ErrorHandler& error_handler_)
        {
            namespace phx = boost::phoenix;
            using boost::phoenix::function;

            error_handler = function<ErrorHandler>(error_handler_)
                (std::string("Error! "),
                 phx::arg_names::_2,
                 phx::cref(error_handler_.iters)[phx::arg_names::_1]);
        }

        double result() {
            BOOST_ASSERT(evaluation_stack_.size() == 1);
            client::function::argument_t res = evaluation_stack_.back();
            double result = boost::get<double>(res);
            evaluation_stack_.pop_back();
            return result;
        }

        void registerFunction(std::string name,
                              client::function::type callback) {
            functions_.insert(client::function::named_t(name, callback));
        }

        void registerFunctions (
                std::map<std::string, client::function::type> functions) {
            functions_.insert(functions.begin(), functions.end());
        }

        void registerVariables(
                std::map<std::string, double> variableDictionary) {
            variableDictionary_ = variableDictionary;
        }

        // visitor
        bool operator()(ast::nil) { BOOST_ASSERT(0); return false; }
        bool operator()(unsigned int x);
        bool operator()(double x);
        bool operator()(bool x);
        bool operator()(ast::quoted_string const& x);
        bool operator()(ast::identifier const& x);
        bool operator()(ast::operation const& x);
        bool operator()(ast::unary const& x);
        bool operator()(ast::function_call const& x);
        bool operator()(ast::expression const& x);

    private:

        boost::function<void(int tag, std::string const& what)> error_handler;

        std::map<std::string, double> variableDictionary_;
        std::map<std::string, client::function::type> functions_;

        //our stack is conform to function call arguments
        client::function::arguments_t evaluation_stack_;
    };
}}

#endif