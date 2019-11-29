#include <memory>
#include <vector>
#include <iostream>

class Component;

/** MaxwellTest is a small utility for testing field objects
 * 
 *  Calculated numerical derivatives for calculating maxwell's equations and
 *  checking that field maps are physical.
 */
class MaxwellTest {

  public:
    MaxwellTest(Vector_t dR, double dt, Component* field);

    Component* getField() const {return field_m.get();}
    void setField(Component* field) {field_m.reset(field);}

    Vector_t getDR() const {return dR_m;}
    void setDR(Vector_t dR) {dR_m = dR;}
    
    std::vector< std::vector<double> > partialsDerivB(const Vector_t &R, double t) const;
    std::vector< std::vector<double> > partialsDerivA(const Vector_t &R, double t) const;
    double divB(const Vector_t &R, double t) const;
    Vector_t curlB(const Vector_t &R, double t) const;
    Vector_t curlA(const Vector_t &R, double t) const;
    Vector_t getB(const Vector_t &R, double t) const;

    std::ostream& printHeading(std::ostream& out) const;
    std::ostream& printLine(std::ostream& out, const Vector_t& R, double t) const;
  private:
    std::unique_ptr<Component> field_m;
    Vector_t dR_m;
    double dt_m;
};
