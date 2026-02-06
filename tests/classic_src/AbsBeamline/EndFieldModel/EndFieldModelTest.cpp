#include "gtest/gtest.h"
#include "Utilities/GeneralClassicException.h"
#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace endfieldmodel {

class EndFieldMockup : public EndFieldModel {
public:
    ~EndFieldMockup() {;}
    double function(double x, int n) const {return x*n;}
    EndFieldMockup* clone() const {return nullptr;}
    std::ostream& print(std::ostream& out) const {return out;}
    void setMaximumDerivative(size_t n) {if(n) {}}
    void rescale(double scaleFactor) { if(scaleFactor) {}}
    double getEndLength() const {return 0;}
    double getCentreLength() const {return 0;}
    double getStartOffset() const {return 0;}
private:

};

TEST(EndFieldModelTest, TestGetSetEndFieldModel) {
    EXPECT_THROW(EndFieldModel::getEndFieldModel("Cheese"), GeneralClassicException);
    EndFieldModel* mockup = dynamic_cast<EndFieldModel*>(new EndFieldMockup());
    std::shared_ptr<endfieldmodel::EndFieldModel> mockupPtr(mockup);
    EndFieldModel::setEndFieldModel("Cheese", mockupPtr);
    std::shared_ptr<endfieldmodel::EndFieldModel> mockupTest =
                            EndFieldModel::getEndFieldModel("Cheese");
    EXPECT_EQ(mockupTest, mockupPtr);
    EXPECT_EQ(EndFieldModel::getName(mockupPtr), "Cheese");
}
}
