#ifndef GENERALCLASSICEXCEPTION_H
#define GENERALCLASSICEXCEPTION_H

#include "Utilities/ClassicException.h"

class GeneralClassicException:public ClassicException
{
public:
    GeneralClassicException(const std::string &meth, const std::string &msg);

    GeneralClassicException(const GeneralClassicException &);
    virtual ~GeneralClassicException();

private:

    // Not implemented.
    GeneralClassicException();
};

#endif
