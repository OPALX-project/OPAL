#ifndef OPAL_INDEXMAP_H
#define OPAL_INDEXMAP_H

#include <ostream>
#include <map>

#include "AbsBeamline/Component.h"
#include "Utilities/OpalException.h"

#include <set>
#include <utility>


class IndexMap
{
public:
    typedef std::pair<double, double> key_t;
    typedef std::set<std::shared_ptr<Component> > value_t;

    IndexMap();

    void add(key_t::first_type initialStep, key_t::second_type finalStep, const value_t &val);

    value_t query(key_t::first_type s, key_t::second_type ds);

    void tidyUp(double zstop);

    void print(std::ostream&) const;
    void saveSDDS(double startS) const;
    size_t size() const;

    size_t numElements() const;
    std::pair<double, double> getRange(const IndexMap::value_t::value_type &element,
                                       double position) const;
    IndexMap::value_t getTouchingElements(const std::pair<double, double> &range);

    class OutOfBounds: public OpalException {
    public:
        OutOfBounds(const std::string &meth, const std::string &msg):
            OpalException(meth, msg) { }

        OutOfBounds(const OutOfBounds &rhs):
            OpalException(rhs) { }

        virtual ~OutOfBounds() { }

    private:
        OutOfBounds();
    };

private:
    class myCompare {
    public:
        bool operator()(const key_t x , const key_t y) const
        {
            if (x.first < y.first) return true;

            if (x.first == y.first) {
                if (x.second < y.second) return true;
            }

            return false;
        }
    };

    typedef std::map<key_t, value_t, myCompare> map_t;
    typedef std::multimap<value_t::value_type, key_t> invertedMap_t;
    map_t mapRange2Element_m;
    invertedMap_t mapElement2Range_m;

    double totalPathLength_m;

    static bool almostEqual(double, double);
    static const double oneMinusEpsilon_m;
};

inline
size_t IndexMap::size() const {
    return mapRange2Element_m.size();
}

inline
std::ostream& operator<< (std::ostream &out, const IndexMap &im)
{
    im.print(out);
    return out;
}

inline
Inform& operator<< (Inform &out, const IndexMap &im) {
    im.print(out.getStream());
    return out;
}

#endif