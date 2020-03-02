#include <sstream>

#include "gtest/gtest.h"
#include "Elements/OpalPolynomialTimeDependence.h"

#include "opal_test_utilities/SilenceTest.h"

class BeamlineVisitor;
// some comment
class TestElement : public ElementBase {
  public:

    TestElement() : ElementBase(), base(NULL), type("") {}
    ElementBase::ElementType getType() const {return ElementBase::ANY;}
    BGeometryBase  &getGeometry() {
        return *base;
    }

    const BGeometryBase  &getGeometry() const {
        return *base;
    }

    ElementBase* clone() const {return NULL;}

    void accept(BeamlineVisitor& /*visitor*/) const {}

  private:
    BGeometryBase* base;
    std::string type;
};

TEST(OpalPolynomialTimeDependenceTest, ConstructorTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalPolynomialTimeDependence dep;
    OpalPolynomialTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(OpalPolynomialTimeDependenceTest, PrintTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalPolynomialTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "POLYNOMIAL_TIME_DEPENDENCE;\n");
}

TEST(OpalPolynomialTimeDependenceTest, UpdateTest) {
    OpalTestUtilities::SilenceTest silencer;

    // std::cerr << "HELLO" << std::endl;
    TestElement element;
    // std::cerr << "WORLD" << std::endl;
    element.setAttribute("P0", 99.);
    // std::cerr << "ATTR" << std::endl;
    OpalPolynomialTimeDependence dependence;
    // std::cerr << "FILL" << std::endl;
    // makes a segmentation fault...
    // dependence.fillRegisteredAttributes(element, OpalElement::IDEAL_FLAG);
    // std::cerr << "DONE" << std::endl;

}