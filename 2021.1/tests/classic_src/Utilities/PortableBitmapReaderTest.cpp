#include "gtest/gtest.h"
#include "OpalSourcePath.h"

#include "Utilities/PortableBitmapReader.h"
#include "opal_test_utilities/SilenceTest.h"

#include <iostream>

TEST(PBMReaderTest, SimpleAsciiTest) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled.pbm";
    PortableBitmapReader reader(pathToBitmapFile);

    bool pixel = reader.isBlack(111, 300);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(111, 319);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(111, 320);
    ASSERT_FALSE(pixel);

    pixel = reader.isBlack(199, 300);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(200, 300);
    ASSERT_FALSE(pixel);
}

TEST(PBMReaderTest, SimpleBinaryTest) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled_binary.pbm";
    PortableBitmapReader reader(pathToBitmapFile);

    bool pixel = reader.isBlack(111, 300);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(111, 319);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(111, 320);
    ASSERT_FALSE(pixel);

    pixel = reader.isBlack(199, 300);
    ASSERT_TRUE(pixel);

    pixel = reader.isBlack(200, 300);
    ASSERT_FALSE(pixel);
}