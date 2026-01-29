#include "gtest/gtest.h"

#include "Distribution/Distribution.h"
#include "Attributes/Attributes.h"
#include "Attributes/PredefinedString.h"
#include "Parser/SimpleStatement.h"
#include "Utilities/ParseError.h"
#include "ValueDefinitions/StringConstant.h"
#include "AbstractObjects/OpalData.h"

#include "opal_test_utilities/SilenceTest.h"

using namespace Attributes;

TEST(PredefinedStringTest, TestDistributionType) {
    OpalTestUtilities::SilenceTest silencer;

    StringConstant stringConst;
    Distribution dist;

    PredefinedString *typeAttribute = dynamic_cast<PredefinedString*>(&dist.itsAttr[Attrib::Distribution::TYPE].getHandler());
    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "GAUSS");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true));
        EXPECT_EQ(Attributes::getString(dist.itsAttr[Attrib::Distribution::TYPE]), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "GUNGAUSSFLATTOPTH");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true));
        EXPECT_EQ(Attributes::getString(dist.itsAttr[Attrib::Distribution::TYPE]), "GUNGAUSSFLATTOPTH");
    }

    {
        StringConstant *astraFlattop = dynamic_cast<StringConstant *>(OpalData::getInstance()->find("ASTRAFLATTOPTH"));
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, astraFlattop->getString());
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true));
        EXPECT_EQ(Attributes::getString(dist.itsAttr[Attrib::Distribution::TYPE]), astraFlattop->getString());
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "gauss");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true));
        EXPECT_EQ(Attributes::getString(dist.itsAttr[Attrib::Distribution::TYPE]), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "Gauss");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_NO_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true));
        EXPECT_EQ(Attributes::getString(dist.itsAttr[Attrib::Distribution::TYPE]), "GAUSS");
    }

    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "guass");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true), ParseError);
    }
    {
        Token token("PredefinedStringTest.in", 1, Token::IS_STRING, "");
        Statement::TokenList tokenList({token});
        SimpleStatement statement("PredefinedString", tokenList);

        EXPECT_THROW(typeAttribute->parse(dist.itsAttr[Attrib::Distribution::TYPE], statement, true), ParseError);
    }
}