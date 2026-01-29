#ifndef PORTABLEBITMAPRREADER_H
#define PORTABLEBITMAPRREADER_H

#include "Utilities/OpalException.h"

#include <string>
#include <vector>
#include <istream>

class PortableBitmapReader {
public:
    PortableBitmapReader(const std::string & input);

    unsigned int getWidth() const;
    unsigned int getHeight() const;

    bool isBlack(unsigned int i, unsigned int j) const;
    std::vector<bool> getPixels() const;
    void print(std::ostream &out) const;

private:
    void readHeader(std::istream &in);
    void readImageAscii(std::istream &in);
    void readImageBinary(std::istream &in);
    std::string getNextPart(std::istream &in);

    unsigned int getIdx(unsigned int h, unsigned int w) const;

    unsigned int width_m;
    unsigned int height_m;

    enum FileType {
        ASCII,
        BINARY
    };

    FileType type_m;

    std::vector<bool> pixels_m;
};

inline
unsigned int PortableBitmapReader::getWidth() const {
    return width_m;
}

inline
unsigned int PortableBitmapReader::getHeight() const {
    return height_m;
}

inline
bool PortableBitmapReader::isBlack(unsigned int i, unsigned int j) const {
    return pixels_m[getIdx(i, j)];
}

inline
std::vector<bool> PortableBitmapReader::getPixels() const {
    return pixels_m;
}

inline
unsigned int PortableBitmapReader::getIdx(unsigned int h, unsigned int w) const {
    if (h >= height_m || w >= width_m) throw OpalException("PortableBitmapReader::getIdx",
                                                           "Pixel number out of bounds");
    return h * width_m + w;
}

#endif