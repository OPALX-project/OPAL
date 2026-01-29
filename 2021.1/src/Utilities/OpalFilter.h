#ifndef OPAL_FILTER_HH
#define OPAL_FILTER_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalFilter.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalFilter
//
// ------------------------------------------------------------------------
//
// $Date: 2008/10/14 09:33:44 $
// $Author: Christof Kraus $
//
// ------------------------------------------------------------------------

#include <vector>
#include <memory>
#include "AbstractObjects/Definition.h"
#include "Filters/Filter.h"

// Class OpalFilter
// ------------------------------------------------------------------------
/// The FILTER definition.
//  A FILTER definition is used to define a filter which can be applied
//  to a 1D histogram in order to get rid of noise.

class OpalFilter: public Definition {

public:

    /// Exemplar constructor.
    OpalFilter();

    virtual ~OpalFilter();

    /// Test if replacement is allowed.
    //  Can replace only by another OpalFilter
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalFilter *clone(const std::string &name);

    /// Check the OpalFilter data.
    virtual void execute();

    /// Find named FILTER.
    static OpalFilter *find(const std::string &name);

    /// Update the OpalFilter data.
    virtual void update();

    void print(std::ostream &os) const;

    void initOpalFilter();

    inline void apply(std::vector<double> &histogram);
    inline void calc_derivative(std::vector<double> &histogram, const double &hz);

    Filter *filter_m;
private:

    // Not implemented.
    OpalFilter(const OpalFilter &);
    void operator=(const OpalFilter &);

    // Clone constructor.
    OpalFilter(const std::string &name, OpalFilter *parent);

};

void OpalFilter::apply(std::vector<double> &histogram) {
    if(filter_m)
        filter_m->apply(histogram);
}

void OpalFilter::calc_derivative(std::vector<double> &histogram, const double &hz) {
    if(filter_m)
        filter_m->calc_derivative(histogram, hz);
}

inline std::ostream &operator<<(std::ostream &os, const OpalFilter &b) {
    b.print(os);
    return os;
}

#endif // OPAL_FILTER_HH
