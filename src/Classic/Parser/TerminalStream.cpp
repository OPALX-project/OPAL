// ------------------------------------------------------------------------
// $RCSfile: TerminalStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TerminalStream
//   Implements an input buffer for reading tokens from a terminal.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:33:11 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Parser/TerminalStream.h"

#ifdef READLINE
// Special version using GNU readline().
// This version does not work with I/O redirection.
extern "C"
{
#include <stdio.h>
#ifdef __P
    // Corrects a problem in the readline header files
#undef __P
#endif
#include <readline/readline.h>
#include <readline/history.h>
}
#endif // READLINE
#include <iomanip>
#include <iostream>


// Class TerminalStream
// ------------------------------------------------------------------------

// :FIXME: can READLINE support be removed? Why is it not used?
TerminalStream::TerminalStream(
#ifdef READLINE
    const char program[]
#else
    const char[]
#endif
    ): AbsFileStream("standard input") {
#ifdef READLINE
    // Set up the readline() function.
    rl_readline_name = new char[strlen(program) + 1];
    strcpy(const_cast<char *>(rl_readline_name), program);
    rl_initialize();
#endif
}


TerminalStream::~TerminalStream() {
#ifdef READLINE
    delete [] rl_readline_name;
#endif
}


bool TerminalStream::fillLine() {
#ifdef READLINE
    char *p = readline("==>");
    line = std::string(p) + '\n';
    add_history(p);
    line += "\n";
    curr_line++;
    std::cerr.width(5);
    std::cerr << curr_line << " " << line;
    curr_char = 0;
    return true;
#else
    if(std::cin.eof()) {
        // We must test for end of file, even on terminal, as it may be rerouted.
        return false;
    } else {
        std::cerr << "==>";
        std::getline(std::cin, line, '\n');
        line += "\n";
        curr_line++;
        std::cerr.width(5);
        std::cerr << curr_line << " " << line;
        curr_char = 0;
        return true;
    }
#endif
}
