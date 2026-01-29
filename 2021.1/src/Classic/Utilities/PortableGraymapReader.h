#ifndef PORTABLEGRAYMAPRREADER_H
#define PORTABLEGRAYMAPRREADER_H

#include "Utilities/OpalException.h"

#include <string>
#include <vector>
#include <istream>

class PortableGraymapReader {
public:
    PortableGraymapReader(const std::string & input);

    unsigned int getWidth() const;
    unsigned int getHeight() const;

    unsigned short getPixel(unsigned int i, unsigned int j) const;
    std::vector<unsigned short> getPixels() const;
    void print(std::ostream &out) const;

private:
    void readHeader(std::istream &in);
    void readImageAscii(std::istream &in);
    void readImageBinary(std::istream &in);
    std::string getNextPart(std::istream &in);

    unsigned int getIdx(unsigned int h, unsigned int w) const;

    unsigned int width_m;
    unsigned int height_m;
    unsigned short depth_m;

    enum FileType {
        ASCII,
        BINARY
    };

    FileType type_m;

    std::vector<unsigned short> pixels_m;
};

inline
unsigned int PortableGraymapReader::getWidth() const {
    return width_m;
}

inline
unsigned int PortableGraymapReader::getHeight() const {
    return height_m;
}

inline
unsigned short PortableGraymapReader::getPixel(unsigned int i, unsigned int j) const {
    return pixels_m[getIdx(i, j)];
}

inline
std::vector<unsigned short> PortableGraymapReader::getPixels() const {
    return pixels_m;
}

inline
unsigned int PortableGraymapReader::getIdx(unsigned int h, unsigned int w) const {
    if (h >= height_m || w >= width_m) throw OpalException("PortableGraymapReader::getIdx",
                                                           "Pixel number out of bounds");
    return h * width_m + w;
}

#endif