/*=============================================================================
    Adapted from boost spirit mini_c example.

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(EXPRESSION_HPP)
#define EXPRESSION_HPP

///////////////////////////////////////////////////////////////////////////////
// Spirit v2.5 allows you to suppress automatic generation
// of predefined terminals to speed up complation. With
// BOOST_SPIRIT_NO_PREDEFINED_TERMINALS defined, you are
// responsible in creating instances of the terminals that
// you need (e.g. see qi::uint_type uint_ below).
#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Uncomment this if you want to enable debugging
// #define BOOST_SPIRIT_QI_DEBUG
///////////////////////////////////////////////////////////////////////////////

#include <boost/spirit/include/qi.hpp>
#include "ast.hpp"
#include "error_handler.hpp"
#include "skipper.hpp"
#include <vector>

namespace client { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    ///////////////////////////////////////////////////////////////////////////////
    //  The expression grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct expression : qi::grammar<Iterator, ast::expression(), qi::locals<char>, skipper<Iterator> >
    {
        expression(error_handler<Iterator>& error_handler);

        qi::rule<Iterator, ast::expression(), qi::locals<char>, skipper<Iterator> >
            expr, equality_expr, relational_expr,
            logical_or_expr, logical_and_expr,
            additive_expr, multiplicative_expr
            ;

        qi::rule<Iterator, ast::operand(), qi::locals<char>, skipper<Iterator> >
            unary_expr, primary_expr, constant_expr
            ;

        qi::rule<Iterator, ast::function_call(), qi::locals<char>, skipper<Iterator> >
            function_call
            ;

        qi::rule<Iterator, std::list<ast::function_call_argument>(), qi::locals<char>, skipper<Iterator> >
            argument_list
            ;

        qi::rule<Iterator, std::string(), qi::locals<char>, skipper<Iterator> >
            identifier
            ;

        qi::rule<Iterator, std::string(), qi::locals<char>, skipper<Iterator> >
            quoted_string
            ;

        qi::symbols<char, ast::optoken>
            logical_or_op, logical_and_op,
            equality_op, relational_op,
            additive_op, multiplicative_op, unary_op
            ;

        qi::symbols<char>
            keywords
            ;
    };
}}

#endif