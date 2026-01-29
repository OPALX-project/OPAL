/*=============================================================================
    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#include "evaluator.hpp"

#include <boost/variant/apply_visitor.hpp>

namespace client { namespace code_gen {

    bool StackEvaluator::operator()(unsigned int x) {

        evaluation_stack_.push_back(static_cast<double>(x));
        return true;
    }

    bool StackEvaluator::operator()(double x) {

        evaluation_stack_.push_back(x);
        return true;
    }

    bool StackEvaluator::operator()(bool x) {

        evaluation_stack_.push_back(static_cast<double>(x));
        return true;
    }

    bool StackEvaluator::operator()(ast::quoted_string const& x) {

        evaluation_stack_.push_back(x.value);
        return true;
    }

    bool StackEvaluator::operator()(ast::identifier const& x) {

        std::map<std::string, double>::iterator i =
            variableDictionary_.find(x.name);
        if(i == variableDictionary_.end()) {
            std::cout << "Undefined variable " << x.name << std::endl;
            return false;
        }

        evaluation_stack_.push_back(i->second);
        return true;
    }

    bool StackEvaluator::operator()(ast::operation const& x) {

        if (!boost::apply_visitor(*this, x.operand_))
            return false;

        double op2 = boost::get<double>(evaluation_stack_.back());
        evaluation_stack_.pop_back();
        double op1 = boost::get<double>(evaluation_stack_.back());
        evaluation_stack_.pop_back();
        double res = 0.0;

        switch (x.operator_) {
            case ast::op_plus          : res = op1 + op2;  break;
            case ast::op_minus         : res = op1 - op2;  break;
            case ast::op_times         : res = op1 * op2;  break;
            case ast::op_divide        : res = op1 / op2;  break;

            case ast::op_equal         : res = op1 == op2; break;
            case ast::op_not_equal     : res = op1 != op2; break;
            case ast::op_less          : res = op1 < op2;  break;
            case ast::op_less_equal    : res = op1 <= op2; break;
            case ast::op_greater       : res = op1 > op2;  break;
            case ast::op_greater_equal : res = op1 >= op2; break;

            case ast::op_and           : res = op1 && op2; break;
            case ast::op_or            : res = op1 || op2; break;
            default                    : BOOST_ASSERT(0);  return false;
        }

        evaluation_stack_.push_back(res);
        return true;
    }

    bool StackEvaluator::operator()(ast::unary const& x) {

        if (!boost::apply_visitor(*this, x.operand_))
            return false;

        double op = boost::get<double>(evaluation_stack_.back());
        evaluation_stack_.pop_back();

        switch (x.operator_) {
            case ast::op_negative : op = -op; break;
            case ast::op_not      : op = !op; break;
            case ast::op_positive :           break;
            default               : BOOST_ASSERT(0); return false;
        }


        evaluation_stack_.push_back(op);
        return true;
    }

    bool StackEvaluator::operator()(ast::function_call const& x) {

        for(ast::function_call_argument const& arg : x.args) {
            if (!boost::apply_visitor(*this, arg))
                return false;
        }

        std::vector<client::function::argument_t> args(x.args.size());
        for(size_t i=0; i < x.args.size(); i++) {
            client::function::argument_t arg = evaluation_stack_.back();
            args[x.args.size()-1-i] = arg;
            evaluation_stack_.pop_back();
        }

        std::map<std::string, client::function::type>::iterator itr =
            functions_.find(x.function_name.name);
        if(itr == functions_.end()) {
            std::cout << "Undefined function "
                      << x.function_name.name << std::endl;
            return false;
        }

        if(! boost::get<1>(itr->second(args))) return false;

        double function_result = boost::get<0>(itr->second(args));
        evaluation_stack_.push_back(function_result);

        return true;
    }

    bool StackEvaluator::operator()(ast::expression const& x) {

        if (!boost::apply_visitor(*this, x.first))
            return false;

        for(ast::operation const& oper : x.rest) {
            if (!(*this)(oper))
                return false;
        }

        return true;
    }

}}

