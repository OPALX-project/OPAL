// ------------------------------------------------------------------------
// $RCSfile: StringStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: StringStream
//   Implements an input buffer for reading from a string.
//   This string is intended for storing valid MAD-9 expressions,
//   it must not contain comments.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/StringStream.h"
#include "Parser/Token.h"
#include "Utilities/FormatError.h"
#include <cctype>

// Class StringStream
// ------------------------------------------------------------------------

StringStream::StringStream(const std::string &str):
    TokenStream("expression"),
    line_m(str + '\n'),
    currentChar_m(0)
{}


StringStream::~StringStream()
{}


Token StringStream::readToken() {
    if(put_back_flag) {
        put_back_flag = false;
        return put_back;
    }

    while(true) {
        if(currentChar_m >= line_m.length()  ||  line_m[currentChar_m] == '\n') {
            return Token("string", 1, Token::IS_EOF, "EOF");
        } else if(isspace(line_m[currentChar_m])) {
            currentChar_m++;
        } else {
            break;
        }
    }

    // First character.
    char ch = line_m[currentChar_m];

    if(ch == '"'  ||  ch == '\'') {
        // String token.
        return readString();
    } else if(isalpha(ch)) {
        // Word token.
        return readWord();
    } else if(isdigit(ch) ||
              (ch == '.' && isdigit(line_m[currentChar_m+1]))) {
        // Numeric token.
        return readNumber();
    } else {
        // Delimiter tokens.
        if(ch == '<'  &&  line_m[currentChar_m+1] == '=') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "<=");
        } else if(ch == '>'  &&  line_m[currentChar_m+1] == '=') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, ">=");
        } else if(ch == '='  &&  line_m[currentChar_m+1] == '=') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "==");
        } else if(ch == '!'  &&  line_m[currentChar_m+1] == '=') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "!=");
        } else if(ch == '|'  &&  line_m[currentChar_m+1] == '|') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "||");
        } else if(ch == '&'  &&  line_m[currentChar_m+1] == '&') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "&&");
        } else if(ch == ':'  &&  line_m[currentChar_m+1] == '=') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, ":=");
        } else if(ch == '-'  &&  line_m[currentChar_m+1] == '>') {
            currentChar_m += 2;
            return Token("string", 1, Token::IS_DELIMITER, "->");
        } else {
            currentChar_m++;
            return Token("string", 1, Token::IS_DELIMITER, ch);
        }
    }

    return Token(stream_name, curr_line, Token::IS_ERROR, "ERROR");
}


Token StringStream::readNumber() {
    bool digit = false;
    bool eflag = false;
    double value = 0.0;
    int expsig = 1;
    int expval = 0;
    int places = 0;
    int lex_pos = currentChar_m;

    while(isdigit(line_m[currentChar_m])) {
        // Digits preceding decimal point.
        value = 10.0 * value + double(line_m[currentChar_m] - '0');
        digit = true;
        currentChar_m++;
    }

    if(digit && line_m[currentChar_m] != '.' && toupper(line_m[currentChar_m]) != 'E') {
        // Unsigned integer seen.
        std::string lexeme(line_m.data() + lex_pos, currentChar_m - lex_pos);
        return Token("string", 1, lexeme, int(value + 0.5));
    }

    // Decimal point.
    if(line_m[currentChar_m] == '.') {
        currentChar_m++;

        // Digits following decimal point.
        while(isdigit(line_m[currentChar_m])) {
            value = 10.0 * value + double(line_m[currentChar_m++] - '0');
            digit = true;
            places++;
        }
    }

    if(! digit) eflag = true;

    // Exponent ?
    if(toupper(line_m[currentChar_m]) == 'E') {
        currentChar_m++;
        digit = false;

        if(line_m[currentChar_m] == '+') {
            currentChar_m++;
        } else if(line_m[currentChar_m] == '-') {
            currentChar_m++;
            expsig = -1;
        }

        while(isdigit(line_m[currentChar_m])) {
            expval = 10 * expval + (line_m[currentChar_m++] - '0');
            digit = true;
        }

        if(! digit) eflag = true;

        // Skip over any non-punctuation characters.
        char ch = line_m[currentChar_m];

        while(! isspace(ch)  &&  ! ispunct(ch)) {
            eflag = true;
            currentChar_m++;
            ch = line_m[currentChar_m];
        }
    }

    // Put pieces together.
    std::string lexeme(line_m.data() + lex_pos, currentChar_m - lex_pos);

    if(eflag) {
        return Token("string", 1, Token::IS_ERROR,
                     "Invalid numeric token \"" + lexeme + "\".");
    } else {
        int power = expsig * expval - places;

        if(power > 0) {
            for(places = power; places > 0; places--)
                value *= 10.0;
        } else {
            for(places = - power; places > 0; places--)
                value /= 10.0;
        }

        return Token("string", 1, lexeme, value);
    }
}


Token StringStream::readString() {
    std::string lexeme;

    if(line_m[currentChar_m] == '"'  ||  line_m[currentChar_m] == '\'') {
        char quote = line_m[currentChar_m];
        currentChar_m++;

        while(true) {
            if(currentChar_m >= line_m.length()) {
                throw FormatError("StringStream::readString()",
                                  "String not terminated.");
            }

            if(line_m[currentChar_m] == quote) {
                currentChar_m++;
                if(line_m[currentChar_m] != quote) break;
            } else if(line_m[currentChar_m] == '\\') {
                currentChar_m++;
            }

            lexeme += line_m[currentChar_m++];
        }
    }

    return Token("string", 1, Token::IS_STRING, lexeme);
}


Token StringStream::readWord() {
    std::string lexeme;
    char ch = line_m[currentChar_m];

    if(isalpha(line_m[currentChar_m])) {
        lexeme += toupper(ch);
        char ch = line_m[++currentChar_m];

        while(isalnum(ch) || ch == '_' || ch == '.' || ch == '\'') {
            lexeme += toupper(ch);
            ch = line_m[++currentChar_m];
        }
    }

    return Token("string", 1, Token::IS_WORD, lexeme);
}