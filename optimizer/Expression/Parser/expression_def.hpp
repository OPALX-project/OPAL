#include "expression.hpp"
#include "error_handler.hpp"

#include <cctype>
#include <cstdlib>
#include <iterator>
#include <optional>

namespace client { namespace parser
{
    template <typename Iterator>
    expression<Iterator>::expression(error_handler<Iterator>& error_handler)
      : error_handler_(error_handler)
      , annotation_(error_handler.iters)
    {
    }

    template <typename Iterator>
    void expression<Iterator>::skip_ws(Iterator& it, Iterator last) const
    {
        while (it != last) {
            if (std::isspace(static_cast<unsigned char>(*it))) {
                ++it;
            } else if (*it == '/' && std::next(it) != last && *std::next(it) == '*') {
                std::advance(it, 2);
                while (it != last && !(*it == '*' && std::next(it) != last && *std::next(it) == '/')) {
                    ++it;
                }
                if (it != last) {
                    std::advance(it, 2); // consume closing "*/"
                }
            } else {
                break;
            }
        }
    }

    template <typename Iterator>
    bool expression<Iterator>::match(Iterator& it, Iterator last, char c) const
    {
        skip_ws(it, last);
        if (it != last && *it == c) {
            ++it;
            return true;
        }
        return false;
    }

    template <typename Iterator>
    bool expression<Iterator>::match_op_token(Iterator& it, Iterator last, char const* token) const
    {
        skip_ws(it, last);
        Iterator cur = it;
        char const* w = token;
        while (*w != '\0' && cur != last && *cur == *w) {
            ++cur;
            ++w;
        }
        if (*w != '\0') {
            return false;
        }
        it = cur;
        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::match_word(Iterator& it, Iterator last, char const* word) const
    {
        skip_ws(it, last);
        Iterator cur = it;
        char const* w = word;
        while (*w != '\0' && cur != last && *cur == *w) {
            ++cur;
            ++w;
        }
        if (*w != '\0') {
            return false; // did not fully match the keyword
        }
        if (cur != last && (std::isalnum(static_cast<unsigned char>(*cur)) || *cur == '_')) {
            return false; // e.g. "trueish" is not the keyword "true"
        }
        it = cur;
        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::try_identifier_token(Iterator& it, Iterator last, std::string& word) const
    {
        skip_ws(it, last);
        Iterator save = it;
        if (it == last || !(std::isalpha(static_cast<unsigned char>(*it)) || *it == '_')) {
            return false;
        }
        Iterator cur = it;
        ++cur;
        while (cur != last && (std::isalnum(static_cast<unsigned char>(*cur)) || *cur == '_')) {
            ++cur;
        }
        std::string token(save, cur);
        if (token == "true" || token == "false") {
            return false; // keywords cannot be used as identifiers
        }
        it = cur;
        word = std::move(token);
        return true;
    }

    template <typename Iterator>
    ast::operand expression<Iterator>::to_operand(ast::expression&& e)
    {
        // Optimization over the original Spirit grammar: an expression with
        // no trailing operations carries no extra information, so we avoid
        // the heap allocation of wrapping it and reuse its operand directly
        // (evaluating the wrapped form would just forward to `first`).
        if (e.rest.empty()) {
            return std::move(e.first);
        }
        return ast::operand(std::move(e));
    }

    ///////////////////////////////////////////////////////////////////////////
    // Binary, left-associative precedence level:
    //   level = next (op next)*
    ///////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    template <typename NextLevel, typename OpMatcher>
    bool expression<Iterator>::parse_level(
        Iterator& it, Iterator last, ast::expression& out,
        NextLevel next, OpMatcher match_op)
    {
        Iterator save = it;

        ast::expression firstExpr;
        if (!next(it, last, firstExpr)) {
            it = save;
            return false;
        }

        out.first = to_operand(std::move(firstExpr));
        out.rest.clear();

        while (true) {
            Iterator opPos = it;
            std::optional<ast::optoken> op = match_op(it, last);
            if (!op) {
                it = opPos;
                break;
            }

            ast::expression rhsExpr;
            if (!next(it, last, rhsExpr)) {
                error_handler_(std::string("Error! "), std::string("Expecting operand"), it);
                return false;
            }

            ast::operation operation;
            operation.operator_ = *op;
            operation.operand_ = to_operand(std::move(rhsExpr));
            out.rest.push_back(std::move(operation));
        }

        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_logical_or(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) { return parse_logical_and(i, l, o); },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "||")) return ast::op_or;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_logical_and(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) { return parse_equality(i, l, o); },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "&&")) return ast::op_and;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_equality(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) { return parse_relational(i, l, o); },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "==")) return ast::op_equal;
                if (match_op_token(i, l, "!=")) return ast::op_not_equal;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_relational(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) { return parse_additive(i, l, o); },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "<=")) return ast::op_less_equal;
                if (match_op_token(i, l, ">=")) return ast::op_greater_equal;
                if (match_op_token(i, l, "<"))  return ast::op_less;
                if (match_op_token(i, l, ">"))  return ast::op_greater;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_additive(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) { return parse_multiplicative(i, l, o); },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "+")) return ast::op_plus;
                if (match_op_token(i, l, "-")) return ast::op_minus;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_multiplicative(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_level(it, last, out,
            [this](Iterator& i, Iterator l, ast::expression& o) -> bool {
                ast::operand operand;
                if (!parse_unary(i, l, operand)) return false;
                o.first = std::move(operand);
                o.rest.clear();
                return true;
            },
            [this](Iterator& i, Iterator l) -> std::optional<ast::optoken> {
                if (match_op_token(i, l, "*")) return ast::op_times;
                if (match_op_token(i, l, "/")) return ast::op_divide;
                return std::nullopt;
            });
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_unary(Iterator& it, Iterator last, ast::operand& out)
    {
        Iterator save = it;
        ast::optoken op;
        bool isUnary = true;
        if (match(it, last, '+')) {
            op = ast::op_positive;
        } else if (match(it, last, '-')) {
            op = ast::op_negative;
        } else if (match(it, last, '!')) {
            op = ast::op_not;
        } else {
            isUnary = false;
        }

        if (!isUnary) {
            return parse_primary(it, last, out);
        }

        ast::operand sub;
        if (!parse_unary(it, last, sub)) {
            it = save;
            return false;
        }

        ast::unary u;
        u.operator_ = op;
        u.operand_ = std::move(sub);
        out = std::move(u);
        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_constant(Iterator& it, Iterator last, ast::operand& out)
    {
        skip_ws(it, last);
        Iterator save = it;
        Iterator cur = it;

        bool hasDigits = false;
        while (cur != last && std::isdigit(static_cast<unsigned char>(*cur))) {
            ++cur;
            hasDigits = true;
        }
        if (cur != last && *cur == '.') {
            Iterator afterDot = cur;
            ++afterDot;
            if (afterDot != last && std::isdigit(static_cast<unsigned char>(*afterDot))) {
                cur = afterDot;
                while (cur != last && std::isdigit(static_cast<unsigned char>(*cur))) {
                    ++cur;
                    hasDigits = true;
                }
            } else if (hasDigits) {
                cur = afterDot; // trailing dot after digits, e.g. "3."
            }
        }
        if (!hasDigits) {
            return false;
        }
        if (cur != last && (*cur == 'e' || *cur == 'E')) {
            Iterator expIt = cur;
            ++expIt;
            if (expIt != last && (*expIt == '+' || *expIt == '-')) {
                ++expIt;
            }
            if (expIt != last && std::isdigit(static_cast<unsigned char>(*expIt))) {
                while (expIt != last && std::isdigit(static_cast<unsigned char>(*expIt))) {
                    ++expIt;
                }
                cur = expIt;
            }
        }

        std::string token(save, cur);
        char* end = nullptr;
        double value = std::strtod(token.c_str(), &end);
        if (end != token.c_str() + token.size()) {
            return false;
        }

        it = cur;
        out = value;
        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_primary(Iterator& it, Iterator last, ast::operand& out)
    {
        if (parse_constant(it, last, out)) {
            return true;
        }

        Iterator save = it;
        std::string word;
        if (try_identifier_token(it, last, word)) {
            Iterator afterWord = it;
            if (match(it, last, '(')) {
                ast::function_call fc;
                fc.function_name = ast::identifier(word);
                std::list<ast::function_call_argument> args;
                if (!parse_argument_list(it, last, args)) {
                    return false;
                }
                fc.args = std::move(args);
                if (!match(it, last, ')')) {
                    error_handler_(std::string("Error! "), std::string("Expecting \")\""), it);
                    return false;
                }
                out = std::move(fc);
                annotation_(out, save);
                return true;
            }

            it = afterWord;
            ast::identifier id(word);
            out = id;
            annotation_(out, save);
            return true;
        }
        it = save;

        if (match_word(it, last, "true")) {
            out = true;
            annotation_(out, save);
            return true;
        }
        if (match_word(it, last, "false")) {
            out = false;
            annotation_(out, save);
            return true;
        }

        if (match(it, last, '(')) {
            ast::expression e;
            if (!parse_expr(it, last, e)) {
                return false;
            }
            if (!match(it, last, ')')) {
                error_handler_(std::string("Error! "), std::string("Expecting \")\""), it);
                return false;
            }
            out = to_operand(std::move(e));
            annotation_(out, save);
            return true;
        }

        it = save;
        return false;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_call_argument(Iterator& it, Iterator last, ast::function_call_argument& out)
    {
        Iterator save = it;

        ast::expression e;
        if (parse_expr(it, last, e)) {
            out = std::move(e);
            return true;
        }
        it = save;

        std::string s;
        if (parse_quoted_string(it, last, s)) {
            out = ast::quoted_string(s);
            return true;
        }
        it = save;

        return false;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_argument_list(Iterator& it, Iterator last, std::list<ast::function_call_argument>& out)
    {
        out.clear();

        ast::function_call_argument firstArg;
        if (!parse_call_argument(it, last, firstArg)) {
            return true; // the argument list is optional
        }
        out.push_back(std::move(firstArg));

        while (true) {
            Iterator commaPos = it;
            if (!match(it, last, ',')) {
                it = commaPos;
                break;
            }

            ast::function_call_argument arg;
            if (!parse_call_argument(it, last, arg)) {
                error_handler_(std::string("Error! "), std::string("Expecting argument"), it);
                return false;
            }
            out.push_back(std::move(arg));
        }

        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_quoted_string(Iterator& it, Iterator last, std::string& out)
    {
        skip_ws(it, last);
        if (it == last || (*it != '\'' && *it != '"')) {
            return false;
        }

        char quote = *it;
        Iterator cur = it;
        ++cur;
        std::string value;
        while (cur != last && *cur != quote) {
            value += *cur;
            ++cur;
        }
        if (cur == last) {
            return false; // unterminated string
        }
        ++cur; // consume closing quote
        it = cur;
        out = std::move(value);
        return true;
    }

    template <typename Iterator>
    bool expression<Iterator>::parse_expr(Iterator& it, Iterator last, ast::expression& out)
    {
        return parse_logical_or(it, last, out);
    }

    template <typename Iterator>
    bool expression<Iterator>::parse(Iterator& first, Iterator last, ast::expression& attr)
    {
        Iterator it = first;
        if (!parse_expr(it, last, attr)) {
            error_handler_(std::string("Error! "), std::string("Expecting expression"), it);
            return false;
        }
        skip_ws(it, last);
        first = it;
        return true;
    }
}}