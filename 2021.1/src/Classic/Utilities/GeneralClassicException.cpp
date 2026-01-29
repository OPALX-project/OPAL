#include "Utilities/GeneralClassicException.h"

GeneralClassicException::GeneralClassicException(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}

GeneralClassicException::GeneralClassicException(const GeneralClassicException &rhs):
    ClassicException(rhs)
{}

GeneralClassicException::~GeneralClassicException()
{}
