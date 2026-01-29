//
// Class IndexMap
//
// This class stores and prints the sequence of elements that the referenc particle passes.
// Each time the reference particle enters or leaves an element an entry is added to the map.
// With help of this map one can determine which element can be found at a given position.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
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
    struct Range
    {
        typedef double first_type;
        typedef double second_type;
        first_type begin;
        second_type end;
    };
    typedef Range key_t;
    typedef std::set<std::shared_ptr<Component> > value_t;

    IndexMap();

    void add(key_t::first_type initialStep, key_t::second_type finalStep, const value_t &val);

    value_t query(key_t::first_type s, key_t::second_type ds);

    void tidyUp(double zstop);

    void print(std::ostream&) const;
    void saveSDDS(double startS) const;
    size_t size() const;

    size_t numElements() const;
    key_t getRange(const IndexMap::value_t::value_type &element, double position) const;
    value_t getTouchingElements(const key_t &range) const;

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
            if (x.begin < y.begin) return true;

            if (x.begin == y.begin) {
                if (x.end < y.end) return true;
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