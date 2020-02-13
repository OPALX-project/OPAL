#include "Fields/FMDummy.h"
#include "Fields/Fieldmap.hpp"

using namespace std;

FMDummy::FMDummy(std::string aFilename):
    Fieldmap(aFilename),
    zbegin_m(0.0),
    zend_m(-1e-3) {

    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "Could not determine the file type.\n"
             << "Please check the section about field maps in the user manual.\n";
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    WARNMSG(errormsg_str << "\n" << endl);

    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
    disableFieldmapWarning();
}

FMDummy::~FMDummy()
{ }

void FMDummy::readMap()
{ }

void FMDummy::freeMap()
{ }

bool FMDummy::getFieldstrength(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/) const {
    return false;
}

bool FMDummy::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FMDummy::getFieldDimensions(double &zBegin, double &zEnd, double &/*rBegin*/, double &/*rEnd*/) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
}
void FMDummy::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FMDummy::swap()
{ }

void FMDummy::getInfo(Inform */*msg*/)
{ }

double FMDummy::getFrequency() const {
    static double dummy = 0.0;
    return dummy;
}

void FMDummy::setFrequency(double /*freq*/)
{ }
