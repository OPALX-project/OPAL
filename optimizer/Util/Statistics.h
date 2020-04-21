//
// Class Statistics
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <map>
#include <string>
#include <iostream>

template<typename T>
class Statistics {

public:

    Statistics(std::string name) : stat_name_(name) {}
    ~Statistics() {}

    void registerStatistic(std::string name, T initial_value = 0) {
        std::pair<statistics_iterator_t, bool> statistic_position;
        statistic_position = statistics_.insert(std::pair<std::string, T>(name, initial_value));

        if(statistic_position.second == false)
            std::cout << "Statistic " << statistic_position.first->second << " already exists!" << std::endl;
    }

    void changeStatisticBy(std::string name, T change_by_value) {
        statistics_iterator_t name_at;
        name_at = statistics_.find(name);

        if(name_at != statistics_.end())
            statistics_[name] += change_by_value;
        else
            std::cout << "Statistic " << name << " not registered!" << std::endl;
    }

    T getStatisticValue(std::string name) {
        return statistics_[name];
    }

    void dumpStatistics() {
        std::cout << "Statistics: " << stat_name_ << std::endl;

        T sum = 0;
        for(std::pair<std::string, T> stat : statistics_) {
            sum += stat.second;
            std::cout << "\t" << stat.first << " = " << stat.second << std::endl;
        }

        std::cout << "_________________________" << std::endl;
        std::cout << "Total: " << sum << std::endl;
    }

    void dumpStatistics(std::ostringstream &stream) {
        stream << "Statistics: " << stat_name_ << std::endl;

        T sum = 0;
        for (std::pair<std::string, T> stat : statistics_) {
            sum += stat.second;
            stream << "\t" << stat.first << " = " << stat.second << std::endl;
        }

        stream << "_________________________" << std::endl;
        stream << "Total: " << sum << std::endl;
    }


private:

    typedef typename std::map<std::string, T> statistics_t;
    typedef typename std::map<std::string, T>::iterator statistics_iterator_t;

    statistics_t statistics_;
    std::string stat_name_;

};

#endif
