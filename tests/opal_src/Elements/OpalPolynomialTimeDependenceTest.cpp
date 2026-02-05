#include <sstream>

#include "Algorithms/AbstractTimeDependence.h"
#include "Attributes/Attributes.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

class BeamlineVisitor;
// some comment
class TestElement final : public ElementBase {
public:
    TestElement() : base(nullptr) {
    }
    ElementType getType() const override {
        return ElementType::ANY;
    }
    BGeometryBase& getGeometry() override {
        return *base;
    }

    const BGeometryBase& getGeometry() const override {
        return *base;
    }

    ElementBase* clone() const override {
        return nullptr;
    }

    void accept(BeamlineVisitor& /*visitor*/) const override {
    }

private:
    BGeometryBase* base;
    std::string type;
};

TEST(OpalPolynomialTimeDependenceTest, ConstructorTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalPolynomialTimeDependence dep;
    const OpalPolynomialTimeDependence* dep_clone = dep.clone("new name");
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
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P0], 1.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P1], 2.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P2], 3.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P3], 4.0);
    dependence.update();
    const auto dep      = AbstractTimeDependence::getTimeDependence("SCALE1");
    const auto value    = dep->getValue(0.1);
    constexpr auto shouldBe = 1.0 + 2.0 * 0.1 + 3.0 * 0.1 * 0.1 + 4.0 * 0.1 * 0.1 * 0.1;
    EXPECT_DOUBLE_EQ(value, shouldBe);
}

TEST(OpalPolynomialTimeDependenceTest, UpdateTest2) {
    OpalTestUtilities::SilenceTest silencer;
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setRealArray(dependence.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS],
        {1.0, 2.0, 3.0, 4.0});
    dependence.update();
    const auto dep      = AbstractTimeDependence::getTimeDependence("SCALE1");
    const auto value    = dep->getValue(0.1);
    constexpr auto shouldBe = 1.0 + 2.0 * 0.1 + 3.0 * 0.1 * 0.1 + 4.0 * 0.1 * 0.1 * 0.1;
    EXPECT_DOUBLE_EQ(value, shouldBe);
}

TEST(OpalPolynomialTimeDependenceTest, UpdateTest3) {
    OpalTestUtilities::SilenceTest silencer;
    OpalPolynomialTimeDependence dependence;
    dependence.setOpalName("SCALE1");
    Attributes::setRealArray(dependence.itsAttr[OpalPolynomialTimeDependence::COEFFICIENTS],
        {1.0, 2.0, 3.0, 4.0});
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P0], 1.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P1], 2.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P2], 3.0);
    Attributes::setReal(dependence.itsAttr[OpalPolynomialTimeDependence::P3], 4.0);
    EXPECT_ANY_THROW(dependence.update());
}
