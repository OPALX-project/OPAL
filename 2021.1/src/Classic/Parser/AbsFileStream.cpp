// ------------------------------------------------------------------------
// $RCSfile: AbsFileStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbsFileStream
//   Implements an abstract input buffer for reading from a file.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:37 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Parser/AbsFileStream.h"
#include "Parser/Token.h"
#include "Utilities/FormatError.h"
#include <cctype>
#include <iostream>
#include <algorithm>

// Class AbsFileStream
// ------------------------------------------------------------------------

AbsFileStream::AbsFileStream(const std::string &name):
    TokenStream(name),
    line("\n"),
    curr_char(0)
{}


AbsFileStream::~AbsFileStream()
{}


Token AbsFileStream::readToken() {
    if(put_back_flag) {
        put_back_flag = false;
        return put_back;
    }

    // Skip white space.
    while(true) {
        char ch = line[curr_char];

        if(ch == '\n') {
            // End of line.
            if(! fillLine()) {
                // End of file.
                return Token(stream_name, curr_line, Token::IS_EOF, "EOF");
            }
        } else if(isspace(ch)) {
            // Skip white space.
            curr_char++;
        } else {
            break;
        }
    }

    // First character.
    char ch = line[curr_char];

    if(ch == '/') {
        // Skip comments.
        char ch1 = line[curr_char+1];
        if(ch1 == '/') {
            // C++ like comment.
            if(fillLine()) {
                return readToken();
            } else {
                return Token(stream_name, curr_line, Token::IS_EOF, "EOF");
            }
        } else if(ch1 == '*') {
            // C like comment.
            if(skipComments()) {
                return readToken();
            } else {
                return Token(stream_name, curr_line, Token::IS_ERROR,
                             "Comment not closed.");
            }
        }
    }

    if(ch == '"'  ||  ch == '\'') {
        // String token.
        return readString();
    } else if(isalpha(ch)) {
        // Word token.
        return readWord();
    } else if(isdigit(ch) ||
              (ch == '.' && isdigit(line[curr_char+1]))) {
        // Numeric token.
        return readNumber();
    } else {
        // Delimiter tokens.
        if(ch == '<'  &&  line[curr_char+1] == '=') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "<=");
        } else if(ch == '>'  &&  line[curr_char+1] == '=') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, ">=");
        } else if(ch == '='  &&  line[curr_char+1] == '=') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "==");
        } else if(ch == '!'  &&  line[curr_char+1] == '=') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "!=");
        } else if(ch == '|'  &&  line[curr_char+1] == '|') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "||");
        } else if(ch == '&'  &&  line[curr_char+1] == '&') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "&&");
        } else if(ch == ':'  &&  line[curr_char+1] == '=') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, ":=");
        } else if(ch == '-'  &&  line[curr_char+1] == '>') {
            curr_char += 2;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, "->");
        } else if(ch == '/'  &&  line[curr_char+1] == '*') {
            curr_char += 2;
            if(skipComments()) {
                return readToken();
            } else {
                return Token(stream_name, curr_line, Token::IS_ERROR,
                             "Comment not closed.");
            }
        } else {
            curr_char++;
            return Token(stream_name, curr_line, Token::IS_DELIMITER, ch);
        }
    }

    return Token(stream_name, curr_line, Token::IS_ERROR, "ERROR");
}


Token AbsFileStream::readNumber() {
    bool digit = false;
    bool eflag = false;
    double value = 0.0;
    int expsig = 1;
    int expval = 0;
    int places = 0;
    int lex_pos = curr_char;

    while(isdigit(line[curr_char])) {
        // Digits preceding decimal point.
        value = 10.0 * value + double(line[curr_char] - '0');
        digit = true;
        curr_char++;
    }

    if(digit  &&  line[curr_char] != '.'  &&
       toupper(line[curr_char]) != 'E') {
        // Unsigned integer seen.
        std::string lexeme(line.data() + lex_pos, curr_char - lex_pos);
        return Token(stream_name, curr_line, lexeme, int(value + 0.5));
    }

    // Decimal point.
    if(line[curr_char] == '.') {
        curr_char++;

        // Digits following decimal point.
        while(isdigit(line[curr_char])) {
            value = 10.0 * value + double(line[curr_char++] - '0');
            digit = true;
            places++;
        }
    }

    if(! digit) eflag = true;

    // Exponent ?
    if(toupper(line[curr_char]) == 'E') {
        curr_char++;
        digit = false;

        if(line[curr_char] == '+') {
            curr_char++;
        } else if(line[curr_char] == '-') {
            curr_char++;
            expsig = -1;
        }

        while(isdigit(line[curr_char])) {
            expval = 10 * expval + (line[curr_char++] - '0');
            digit = true;
        }

        if(! digit) eflag = true;

        // Skip over any non-punctuation characters.
        char ch = line[curr_char];

        while(! isspace(ch)  &&  ! ispunct(ch)) {
            eflag = true;
            curr_char++;
            ch = line[curr_char];
        }
    }

    // Put pieces together.
    std::string lexeme(line.data() + lex_pos, curr_char - lex_pos);

    if(eflag) {
        return Token(stream_name, curr_line, Token::IS_ERROR,
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

        return Token(stream_name, curr_line, lexeme, value);
    }
}


Token AbsFileStream::readString() {
    std::string lexeme;

    if(line[curr_char] == '"'  ||  line[curr_char] == '\'') {
        char quote = line[curr_char];
        curr_char++;

        while(true) {
            if(line[curr_char] == '\n') {
                return Token(stream_name, curr_line, Token::IS_ERROR, "String \"" +
                             lexeme + "\" not terminated at end of line.");
            }

            if(line[curr_char] == quote) {
                curr_char++;
                if(line[curr_char] != quote) break;
            } else if(line[curr_char] == '\\') {
                curr_char++;
            }

            lexeme += line[curr_char++];
        }
    }

    return Token(stream_name, curr_line, Token::IS_STRING, lexeme);
}


Token AbsFileStream::readWord() {
    std::string lexeme;
    char ch = line[curr_char];

    if(isalpha(line[curr_char])) {
        lexeme += toupper(ch);
        char ch = line[++curr_char];

        while(isalnum(ch) || ch == '_' || ch == '.') {
            lexeme += toupper(ch);
            ch = line[++curr_char];
        }
    }

    return Token(stream_name, curr_line, Token::IS_WORD, lexeme);
}


bool AbsFileStream::skipComments() {
    curr_char += 2;

    while(true) {
        while(line[curr_char] != '\n') {
            if(line[curr_char] == '*'  &&  line[curr_char+1] == '/') {
                curr_char += 2;
                return true;
            } else if(line[curr_char] == '/'  &&  line[curr_char+1] == '*') {
                skipComments();
            } else {
                curr_char++;
            }
        }

        // End of line
        if(! fillLine()) return false;
    }
}