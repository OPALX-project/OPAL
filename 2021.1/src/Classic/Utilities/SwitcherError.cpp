#include "SwitcherError.h"

SwitcherError::SwitcherError(const std::string &meth, const std::string &msg):
    ClassicException(meth, msg)
{}

SwitcherError::SwitcherError(const SwitcherError &rhs):
    ClassicException(rhs)
{}

SwitcherError::~SwitcherError()
{}
