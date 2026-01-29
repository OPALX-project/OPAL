#include "Utilities/ClassicField.h"

extern Inform *gmsg;

ClassicField::ClassicField(std::shared_ptr<Component> element, const double &start, const double &end):
    element_m(element),
    start_m(start),
    end_m(end),
    is_on_m(false)
{ }

ClassicField::~ClassicField() {
    element_m = NULL;
}

void ClassicField::setOn(const double &kineticEnergy) {
    if(!is_on_m) {
        element_m->goOnline(kineticEnergy);
        INFOMSG(level3 << element_m->getName() << " gone live" << endl);
        is_on_m = true;
    }
}

void ClassicField::setOff() {
    if(is_on_m) {
        element_m->goOffline();
        INFOMSG(level3 << element_m->getName() << " gone off" << endl);
        is_on_m = false;
    }
}
