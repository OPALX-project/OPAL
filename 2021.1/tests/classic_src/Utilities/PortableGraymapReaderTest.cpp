#include "gtest/gtest.h"
#include "OpalSourcePath.h"

#include "Utilities/PortableGraymapReader.h"
#include "opal_test_utilities/SilenceTest.h"

#include <iostream>

TEST(PGMReaderTest, SimpleAscii8Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled8.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    unsigned short pixel = reader.getPixel(199, 299);
    unsigned short expected = 103;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 299);
    expected = 68;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 300);
    expected = 219;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(199, 300);
    expected = 155;
    ASSERT_EQ(expected, pixel);
}

TEST(PGMReaderTest, SimpleAscii16Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled16.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    unsigned short pixel = reader.getPixel(199, 299);
    unsigned short expected = 26610;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 299);
    expected = 17689;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 300);
    expected = 56453;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(199, 300);
    expected = 39879;
    ASSERT_EQ(expected, pixel);
}

TEST(PGMReaderTest, SimpleBinary8Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled_binary8.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    unsigned short pixel = reader.getPixel(199, 299);
    unsigned short expected = 103;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 299);
    expected = 68;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 300);
    expected = 219;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(199, 300);
    expected = 155;
    ASSERT_EQ(expected, pixel);
}

TEST(PGMReaderTest, SimpleBinary16Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled_binary16.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    unsigned short pixel = reader.getPixel(199, 299);
    unsigned short expected = 26610;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 299);
    expected = 17689;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(200, 300);
    expected = 56453;
    ASSERT_EQ(expected, pixel);

    pixel = reader.getPixel(199, 300);
    expected = 39879;
    ASSERT_EQ(expected, pixel);
}