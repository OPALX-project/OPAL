#ifndef SWITCHERERROR_H
#define SWITCHERERROR_H

#include "Utilities/ClassicException.h"

#include <string>

class SwitcherError:public ClassicException
{
public:
    SwitcherError(const std::string &meth, const std::string &msg);

    SwitcherError(const SwitcherError &);
    virtual ~SwitcherError();

private:

    // Not implemented.
    SwitcherError();    
};

#endif
