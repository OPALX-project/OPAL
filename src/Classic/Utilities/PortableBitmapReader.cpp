#include "Utilities/PortableBitmapReader.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

PortableBitmapReader::PortableBitmapReader(const std::string &input) {

    std::ifstream in(input);
    readHeader(in);
    pixels_m.resize(width_m * height_m);

    if (type_m == ASCII) {
        readImageAscii(in);
    } else {
        readImageBinary(in);
    }
}

std::string PortableBitmapReader::getNextPart(std::istream &in) {
    char c;
    do {
        c = in.get();
        if (c == '#') {
            do {
                c = in.get();
            } while (c != '\n');
        } else if (!(c == ' ' ||
                     c == '\t' ||
                     c == '\n' ||
                     c == '\r')) {
            in.putback(c);
            break;
        }
    } while (true);

    std::string nextPart;
    in >> nextPart;

    return nextPart;
}

void PortableBitmapReader::readHeader(std::istream &in) {
    std::string magicValue = getNextPart(in);

    if (magicValue == "P1") {
        type_m = ASCII;
    } else if (magicValue == "P4") {
        type_m = BINARY;
    } else {
        throw OpalException("PortableBitmapReader::readHeader",
                            "Unknown magic value: '" + magicValue + "'");
    }

    {
        std::string tmp = getNextPart(in);
        std::istringstream conv;
        conv.str(tmp);
        conv >> width_m;
    }

    {
        std::string tmp = getNextPart(in);
        std::istringstream conv;
        conv.str(tmp);
        conv >> height_m;
    }

    char tmp;
    in.read(&tmp, 1);
}

void PortableBitmapReader::readImageAscii(std::istream &in) {
    unsigned int size = height_m * width_m;
    unsigned int i = 0;
    while (i < size) {
        char c;
        in >> c;

        if (!(c == ' ' ||
              c == '\n' ||
              c == '\t' ||
              c == '\r')) {
            pixels_m[i] = (c == '1');
            ++ i;
        }
    }
}

void PortableBitmapReader::readImageBinary(std::istream &in) {
    static const unsigned int sizeChar = sizeof(char) * 8;

    unsigned int numPixels = 0;
    unsigned char c = 0;
    for (unsigned int row = 0; row < height_m; ++ row) {
        for (unsigned int col = 0; col < width_m; ++ col) {
            if ( col % sizeChar == 0) {
                char c2 = 0;
                in.read(&c2, 1);
                c = (unsigned char) c2;
            }
            unsigned int k = sizeChar - 1 - (col % sizeChar);
            pixels_m[numPixels] = (c >> k & 1);
            ++ numPixels;
        }
    }
}

void PortableBitmapReader::print(std::ostream &/*out*/) const {
    for (unsigned int i = 0; i < height_m; ++ i) {
        for (unsigned int j = 0; j < width_m; ++ j) {
            if (j % 8 == 0) {
                unsigned int byte = 0;
                for (unsigned int k = 0; k < 8 && j + k < width_m; ++ k) {
                    unsigned int idx = getIdx(i, j + k);
                    if (pixels_m[idx])
                        byte = byte | (1 << (7 - k));
                }
                std::cout << "  " << std::hex << std::setw(2) << std::setfill('0') << byte << ": ";
            }
            unsigned int idx = getIdx(i, j);
            std::cout << pixels_m[idx];
        }
        std::cout << std::endl;
    }
}
