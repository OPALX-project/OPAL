#ifndef __SAMPLE_INDIVIDUAL_H__
#define __SAMPLE_INDIVIDUAL_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#include "Utilities/OpalException.h"


#include "boost/smart_ptr.hpp"

/**
 * \class SampleIndividual
 *  Structure for an individual in the population holding genes
 *  values.
 *
 *  @see Types.h
 */
class SampleIndividual {

public:

    /// representation of genes
    typedef std::vector<double> genes_t;
    /// gene names
    typedef std::vector<std::string> names_t;
    /// objectives array
    typedef std::vector<double> objectives_t;

    SampleIndividual()
    {}

    SampleIndividual(names_t names)
        : names_m(names)
    {
        genes.resize(names.size(), 0.0);
    }

    /// serialization of structure
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/) {
        ar & genes;
        ar & objectives;
        ar & id;
    }

    /// genes of an individual
    genes_t      genes;
    /// values of objectives of an individual
    objectives_t objectives;
    /// id
    unsigned int id = 0;

    int getIndex(std::string name) {
        auto res = std::find(std::begin(names_m), std::end(names_m), name);

        if (res == std::end(names_m)) {
            throw OpalException("SampleIndividual::getIndex()",
                                "Variable '" + name + "' not contained.");
        }
        return std::distance(std::begin(names_m), res);
    }


    std::string getName(size_t i) {
        return names_m[i];
    }

    void print(std::ostream &out) const {
        out << std::setw(8) << id << std::endl;
        for (unsigned int i = 0; i < genes.size(); ++ i) {
            out << names_m[i] << ": " << genes[i] << std::endl;
        }
    }
private:
    /// gene names
    names_t names_m;
};

inline
std::ostream & operator<<(std::ostream & out, const SampleIndividual &ind) {
    ind.print(out);

    return out;
}
#endif