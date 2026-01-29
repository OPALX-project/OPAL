#ifndef CHANGES_H
#define CHANGES_H

#include <map>
#include <string>
namespace Versions {
    extern std::map<unsigned int, std::string> changes;

    void fillChanges();
}

#endif