#include "SDDSParser.h"

#include <boost/algorithm/string.hpp>

SDDS::SDDSParser::SDDSParser():
    sddsFileName_m("")
{ }

SDDS::SDDSParser::SDDSParser(const std::string &input):
    sddsFileName_m(input)
{ }

void SDDS::SDDSParser::setInput(const std::string &input) {
    sddsFileName_m = input;
}

SDDS::file SDDS::SDDSParser::run() {
    typedef std::string::const_iterator iterator_t;
    typedef SDDS::parser::file_parser<iterator_t> file_parser_t;
    typedef SDDS::parser::skipper<iterator_t> skipper_t;
    typedef SDDS::error_handler<iterator_t> error_handler_t;

    sddsData_m.clear();
    paramNameToID_m.clear();
    columnNameToID_m.clear();

    skipper_t skipper;
    std::string contents = readFile();
    iterator_t contentsIter = contents.begin();
    iterator_t contentsEnd = contents.end();
    error_handler_t error_handler(contentsIter, contentsEnd);
    file_parser_t parser(error_handler);

    bool success = phrase_parse(contentsIter, contentsEnd, parser, skipper, sddsData_m);
    {
        SDDS::parameterList::iterator piter = sddsData_m.sddsParameters_m.begin();
        SDDS::parameterList::iterator pend = sddsData_m.sddsParameters_m.end();
        for (; piter != pend && success; ++ piter) {
            success = piter->parse(contentsIter, contentsEnd, skipper);
        }
        while (success && contentsIter != contentsEnd) {
            SDDS::columnList::iterator citer = sddsData_m.sddsColumns_m.begin();
            SDDS::columnList::iterator cend = sddsData_m.sddsColumns_m.end();
            for (; citer != cend && success; ++ citer) {
                success = citer->parse(contentsIter, contentsEnd, skipper);
            }
        }
    }

    if (!success || contentsIter != contentsEnd)
        {
            throw SDDSParserException("SDDSParser::parseSDDSFile",
                                      "could not parse SDDS file");
        }

    unsigned int param_order = 0;
    for (const SDDS::parameter &param: sddsData_m.sddsParameters_m) {
        std::string name = *param.name_m;
        fixCaseSensitivity(name);
        paramNameToID_m.insert(std::make_pair(name,
                                              param_order));
        ++ param_order;
    }

    unsigned int col_order = 0;
    for (const SDDS::column &col: sddsData_m.sddsColumns_m) {
        std::string name = *col.name_m;
        fixCaseSensitivity(name);
        columnNameToID_m.insert(std::make_pair(name,
                                               col_order));
        ++ col_order;
    }

    return sddsData_m;
}

std::string SDDS::SDDSParser::readFile() {
    std::ifstream in(sddsFileName_m.c_str());

    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);

        in.read(&contents[0], contents.size());

        in.close();

        return contents;
    }

    throw SDDSParserException("SDDSParser::readSDDSFile",
                              "could not open file '" + sddsFileName_m + "'");

    return std::string("");
}

SDDS::ast::columnData_t SDDS::SDDSParser::getColumnData(const std::string &columnName) {
    int idx = getColumnIndex(columnName);

    return sddsData_m.sddsColumns_m[idx].values_m;
}


int SDDS::SDDSParser::getColumnIndex(std::string col_name) const {
    fixCaseSensitivity(col_name);
    auto it = columnNameToID_m.find(col_name);
    if (it != columnNameToID_m.end()) {
        return it->second;
    }

    throw SDDSParserException("SDDSParser::getColumnIndex",
                              "could not find column '" + col_name + "'");

}

//XXX use either all upper, or all lower case chars
void SDDS::SDDSParser::fixCaseSensitivity(std::string &for_string) {

    boost::to_lower(for_string);
}