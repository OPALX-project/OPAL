#ifndef CLASSIC_FIELDMAP_ICC
#define CLASSIC_FIELDMAP_ICC


#include <string>
#include <iostream>
#include <fstream>
#include <boost/type_index.hpp>
#include "Fields/Fieldmap.h"

template<class T>
bool Fieldmap::interpreteLine(std::ifstream & in, T & value, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(boost::typeindex::type_id<T>().pretty_name());
        interpreteWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

template<class S, class T>
bool Fieldmap::interpreteLine(std::ifstream & in, S & value1, T & value2, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value2;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(boost::typeindex::type_id<S>().pretty_name());
        expecting += std::string(" ") + boost::typeindex::type_id<T>().pretty_name();
        interpreteWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

template<class S, class T, class U>
bool Fieldmap::interpreteLine(std::ifstream & in, S & value1, T & value2, U & value3, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value3;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<T>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<U>().pretty_name());
        interpreteWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

template<class S, class T, class U, class V>
bool Fieldmap::interpreteLine(std::ifstream & in, S & value1, T & value2, U & value3, V & value4, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;

    interpreter >> value3;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value4;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<T>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<U>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<V>().pretty_name());
        interpreteWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

template<class S>
bool Fieldmap::interpreteLine(std::ifstream & in, S & value1, S & value2, S & value3, S & value4, S & value5, S & value6, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;

    interpreter >> value3;

    interpreter >> value4;

    interpreter >> value5;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value6;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<S>().pretty_name() + std::string(" ") +
                              boost::typeindex::type_id<S>().pretty_name());
        interpreteWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}


#endif