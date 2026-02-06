/*=============================================================================
    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#include "evaluator.hpp"

#include <boost/variant/apply_visitor.hpp>

namespace client { namespace code_gen {

    bool StackEvaluator::operator()(unsigned int x) {
        evaluation_stack_.emplace_back(std::in_place_type<double>, static_cast<double>(x));
        return true;
    }

    bool StackEvaluator::operator()(double x) {
        evaluation_stack_.emplace_back(std::in_place_type<double>, x);
        return true;
    }

    bool StackEvaluator::operator()(bool x) {
        evaluation_stack_.emplace_back(std::in_place_type<double>, static_cast<double>(x));
        return true;
    }

    bool StackEvaluator::operator()(ast::quoted_string const& x) {
        evaluation_stack_.emplace_back(std::in_place_type<std::string>, x.value);
        return true;
    }

    bool StackEvaluator::operator()(ast::identifier const& x) {

        auto i = variableDictionary_.find(x.name);
        if (i == variableDictionary_.end()) {
            std::cout << "Undefined variable " << x.name << std::endl;
            return false;
        }

        evaluation_stack_.emplace_back(std::in_place_type<double>, i->second);
        return true;
    }

    bool StackEvaluator::operator()(ast::operation const& x) {

        if (!boost::apply_visitor(*this, x.operand_)) {
            return false;
        }

        double op2 = std::get<double>(evaluation_stack_.back());
        evaluation_stack_.pop_back();
        double op1 = std::get<double>(evaluation_stack_.back());
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

        evaluation_stack_.emplace_back(std::in_place_type<double>, res);
        return true;
    }

    bool StackEvaluator::operator()(ast::unary const& x) {

        if (!boost::apply_visitor(*this, x.operand_))
            return false;

        double op = std::get<double>(evaluation_stack_.back());
        evaluation_stack_.pop_back();

        switch (x.operator_) {
            case ast::op_negative : op = -op; break;
            case ast::op_not      : op = !op; break;
            case ast::op_positive :           break;
            default               : BOOST_ASSERT(0); return false;
        }

        evaluation_stack_.emplace_back(std::in_place_type<double>, op);
        return true;
    }

    bool StackEvaluator::operator()(ast::function_call const& x) {

        for(ast::function_call_argument const& arg : x.args) {
            if (!boost::apply_visitor(*this, arg))
                return false;
        }

        std::vector<client::function::argument_t> args;
        for (size_t i = x.args.size(); i-- > 0; ) {
            args.emplace_back(std::move(evaluation_stack_.back()));
            evaluation_stack_.pop_back();
        }

        std::map<std::string, client::function::type>::iterator itr =
            functions_.find(x.function_name.name);
        if (itr == functions_.end()) {
            std::cout << "Undefined function "
                      << x.function_name.name << std::endl;
            return false;
        }

        if (! std::get<1>(itr->second(args))) {
            return false;
        }

        double function_result = std::get<0>(itr->second(args));
        evaluation_stack_.emplace_back(std::in_place_type<double>, function_result);

        return true;
    }

    bool StackEvaluator::operator()(ast::expression const& x) {

        if (!boost::apply_visitor(*this, x.first))
            return false;

        for (ast::operation const& oper : x.rest) {
            if (!(*this)(oper))
                return false;
        }

        return true;
    }

}}

