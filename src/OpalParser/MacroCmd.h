//
// Class MacroCmd
//
//   This class parses the MACRO command.
//   Encapsulate the buffer for the ``archetypes'' of all macros.
//   The macro is stored as a MacroStream.  For execution, first the
//   parameters are replaced, then the resulting stream is sent to the parser.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#ifndef CLASSIC_MacroCmd_HH
#define CLASSIC_MacroCmd_HH

#include "OpalParser/Macro.h"
#include "OpalParser/MacroStream.h"
#include "MemoryManagement/Pointer.h"
#include <iosfwd>
#include <string>

// Class MacroCmd
// ------------------------------------------------------------------------

class MacroCmd: public Macro {

public:

    MacroCmd();
    MacroCmd(const std::string &name, MacroCmd *parent);
    virtual ~MacroCmd();

    /// Execute the macro command.
    virtual void execute();

    /// Make a macro instance.
    //  Expects parse pointer in the statement to be set on the first argument.
    //  The parser is used to determine the parse mode
    //  (normal, error, match, edit, track).
    virtual Object *makeInstance
    (const std::string &name, Statement &, const Parser *);

    /// Make a macro template.
    //  Expects parse pointer in the statement to be set on the first argument.
    virtual Object *makeTemplate(const std::string &, TokenStream &, Statement &);

private:

    // Not implemented.
    MacroCmd(const MacroCmd &);
    void operator=(const MacroCmd &);

    // The stream of tokens representing the macro command.
    Pointer<MacroStream> body;

    // Pointer to the parser to be used in execution.
    const Parser *itsParser;
};

#endif // CLASSIC_MacroCmd_HH