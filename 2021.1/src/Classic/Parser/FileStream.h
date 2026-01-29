#ifndef CLASSIC_FileStream_HH
#define CLASSIC_FileStream_HH

// ------------------------------------------------------------------------
// $RCSfile: FileStream.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FileStream
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
#include <fstream>


// Class FileStream
// ------------------------------------------------------------------------
/// A stream of input tokens.
//  The source of tokens is a named disk file.

class FileStream: public AbsFileStream {

public:

    /// Constructor.
    //  Open associated file.
    explicit FileStream(const std::string &name);

    /// Destructor.
    //  Close associated file.
    virtual ~FileStream();

    /// Read next input line.
    virtual bool fillLine();

    /// Set echo flag.
    //  If [b]flag[/b] is true, the subsequent input lines are echoed to
    //  the standard error stream.
    static void setEcho(bool flag);

    /// Return echo flag.
    static bool getEcho();

private:

    // Not implemented.
    FileStream();
    FileStream(const FileStream &);
    void operator=(const FileStream &);

    // Current file stream.
    std::ifstream is;

    // Flag for statement echo.
    static bool echoFlag;
};

#endif // CLASSIC_FileStream_HH
