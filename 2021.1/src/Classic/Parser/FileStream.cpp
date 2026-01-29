// ------------------------------------------------------------------------
// $RCSfile: FileStream.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FileStream
//   Implements an input buffer for reading from a file.
//
// ------------------------------------------------------------------------
// Class category: Parser
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:33:11 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Parser/FileStream.h"
#include "Utilities/ParseError.h"
#include "Utility/IpplInfo.h"
#include <iomanip>
#include <iostream>

// Class FileStream
// ------------------------------------------------------------------------

bool FileStream::echoFlag = false;


FileStream::FileStream(const std::string &name):
    AbsFileStream(name),
    is(name.c_str()) {
    if(is.fail()) {
        throw ParseError("FileStream::FileStream()",
                         "Cannot open file \"" + name + "\".");
    }
}


FileStream::~FileStream()
{}


void FileStream::setEcho(bool flag) {
    echoFlag = flag;
}


bool FileStream::getEcho() {
    return echoFlag;
}


bool FileStream::fillLine() {
    if(is.eof()) {
        // End of file.
        return false;
    } else if(is.bad()) {
        // Read error.
        throw ParseError("FileStream::fillLine()",
                         "Read error on file \"" + stream_name + "\".");
    } else {
        // Stream OK, read next line.
        line = "";
        std::getline(is, line, '\n');
        line += "\n";
        curr_line++;
        if(echoFlag && Ippl::myNode() == 0) {
            std::cerr.width(5);
            std::cerr << curr_line << " " << line;
        }
        curr_char = 0;
        return true;
    }
}