#ifndef __TRACE_H__
#define __TRACE_H__

#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "boost/smart_ptr.hpp"

#include "Util/Trace/TraceComponent.h"

class Trace {

public:

    Trace(std::string name)
        : name_(name)
    {}

    ~Trace()
    {}

    void registerComponent(std::string name,
            boost::shared_ptr<TraceComponent> component) {
        nameToIdx_.insert(
            std::pair<std::string, size_t>(name, pipeline_.size()));
        pipeline_.push_back(component);
    }

    void unregisterComponent(std::string /*name*/) {
        //TODO: set null @ idx
    }

    void log(std::ostringstream &dump) {
        for(boost::shared_ptr<TraceComponent> component : pipeline_) {
            component->execute(dump);
        }
    }

private:

    std::string name_;

    std::vector< boost::shared_ptr<TraceComponent> > pipeline_;
    std::map< std::string, size_t > nameToIdx_;

};

#endif
