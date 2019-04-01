#ifndef OPAL_SAMPLE_H
#define OPAL_SAMPLE_H

#include <string>
#include <memory>
#include "AbstractObjects/Definition.h"

#include "Sample/SamplingMethod.h"

// Class OpalSample
// ------------------------------------------------------------------------
/// The SAMPLING definition.
//  A SAMPLING definition is used run the optimizer in sample mode.

class OpalSample: public Definition {

public:

    /// Exemplar constructor.
    OpalSample();

    virtual ~OpalSample() {};

    /// Make clone.
    virtual OpalSample *clone(const std::string &name);

    /// Check the OpalSample data.
    virtual void execute();

    /// Find sampling method
    static OpalSample *find(const std::string &name);

    void initialize(const std::string &dvarName,
                    double lower,
                    double upper,
                    size_t modulo = 1,
                    bool sequence = false);

    std::string getVariable() const;

    unsigned int getSize() const;

    std::shared_ptr<SamplingMethod> sampleMethod_m;

private:

    ///@{ Not implemented.
    OpalSample  (const OpalSample &) = delete;
    void operator=(const OpalSample &) = delete;
    ///@}
    /// Private copy constructor, called by clone
    OpalSample(const std::string &name, OpalSample *parent);

    unsigned int size_m;
};

inline
unsigned int OpalSample::getSize() const{
    return size_m;
}
#endif
