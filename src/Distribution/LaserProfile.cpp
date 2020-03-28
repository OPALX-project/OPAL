// ------------------------------------------------------------------------
// $RCSfile: LaserProfile.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LaserProfile
//
// ------------------------------------------------------------------------

#include "Distribution/LaserProfile.h"
#include "Utility/Inform.h"
#include "Utilities/OpalException.h"
#include "Utilities/PortableGraymapReader.h"

#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <cmath>

//#define TESTLASEREMISSION

extern Inform *gmsg;

LaserProfile::LaserProfile(const std::string &fileName,
                           const std::string &imageName,
                           double intensityCut,
                           short flags):
    sizeX_m(0),
    sizeY_m(0),
    hist2d_m(NULL),
    rng_m(NULL),
    pdf_m(NULL),
    centerMass_m(0.0),
    standardDeviation_m(0.0){

    unsigned short *image = readFile(fileName, imageName);
    // saveData("originalLaserProfile", image);

    if (flags & FLIPX) flipX(image);
    if (flags & FLIPY) flipY(image);
    if (flags & ROTATE90) {
        swapXY(image);
        flipX(image);
    }
    if (flags & ROTATE180) {
        flipX(image);
        flipY(image);
    }
    if (flags & ROTATE270) {
        swapXY(image);
        flipY(image);
    }

    filterSpikes(image);
    normalizeProfileData(intensityCut, image);
    computeProfileStatistics(image);
    // saveData("processedLaserProfile", image);
    fillHistrogram(image);

    delete[] image;

    setupRNG();

#ifdef TESTLASEREMISSION
    saveHistogram();
    sampleDist();
#endif

    printInfo();
}

LaserProfile::~LaserProfile() {
    gsl_histogram2d_pdf_free(pdf_m);
    gsl_histogram2d_free(hist2d_m);
    gsl_rng_free(rng_m);
}

unsigned short * LaserProfile::readFile(const std::string &fileName,
                                        const std::string &imageName) {
    namespace fs = boost::filesystem;
    if (!fs::exists(fileName)) {
        throw OpalException("LaserProfile::readFile",
                            "given file '" + fileName + "' does not exist");
    }

    size_t npos = fileName.find_last_of('.');
    std::string ext = fileName.substr(npos + 1);

    unsigned short *image;
    if (ext == "pgm") {
        image = readPGMFile(fileName);
    } else {
        image = readHDF5File(fileName, imageName);
    }

    return image;
}

unsigned short * LaserProfile::readPGMFile(const std::string &fileName) {
    PortableGraymapReader reader(fileName);

    sizeX_m = reader.getWidth();
    sizeY_m = reader.getHeight();

    unsigned short *image = new unsigned short[sizeX_m * sizeY_m];
    unsigned int idx = 0;
    for (unsigned int i = 0; i < sizeX_m; ++ i) {
        for (unsigned int j = 0; j < sizeY_m; ++ j, ++ idx) {
            image[idx] = reader.getPixel(j, i);
        }
    }

    return image;
}

unsigned short * LaserProfile::readHDF5File(const std::string &fileName,
                                            const std::string &imageName) {

    hid_t h5 = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group = H5Gopen2 (h5, imageName.c_str(), H5P_DEFAULT);

    if (group < 0) {
        throw OpalException("LaserProfile::readFile",
                            "given image name '" + imageName + "' does not exist");
    }

    hsize_t dim[2];
    hid_t   dataSet = H5Dopen2 (group, "CameraData", H5P_DEFAULT);

    if (dataSet < 0) {
        throw OpalException("LaserProfile::readFile",
                            "data set 'CameraData' does not exist");
    }

    hid_t dataSetSpace = H5Dget_space(dataSet);
    H5Sget_simple_extent_dims(dataSetSpace, dim, NULL);
    hid_t filespace = H5Screate_simple(2, dim, NULL);

    sizeX_m        = dim[0];
    sizeY_m        = dim[1];

    hsize_t startHyperslab[]  = {0, 0};
    hsize_t blockCount[]  = {sizeX_m, sizeY_m};

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, startHyperslab, NULL, blockCount, NULL);

    unsigned short  *image = new unsigned short  [sizeX_m * sizeY_m];
    hid_t mem = H5Screate_simple(2, blockCount, NULL);

    H5Dread(dataSet, H5T_NATIVE_USHORT, mem, filespace, H5P_DEFAULT, image);

    H5Sclose(mem);
    H5Sclose(filespace);
    H5Dclose(dataSetSpace);
    H5Dclose(dataSet);
    H5Gclose(group);
    H5Fclose(h5);

    return image;
}

void LaserProfile::flipX(unsigned short *image) {
    unsigned int col2 = sizeX_m - 1;
    for (unsigned int col1 = 0; col1 < col2; ++ col1, -- col2) {

        unsigned int pixel1 = col1 * sizeY_m;
        unsigned int pixel2 = col2 * sizeY_m;
        for (unsigned int row = 0; row < sizeY_m; ++ row) {
            std::swap(image[pixel1 ++], image[pixel2 ++]);
        }
    }
}

void LaserProfile::flipY(unsigned short *image) {
    for (unsigned int col = 0; col < sizeX_m; ++ col) {

        unsigned int pixel1 = col * sizeY_m;
        unsigned int pixel2 = (col + 1) * sizeY_m - 1;

        unsigned int row2 = sizeY_m - 1;
        for (unsigned int row1 = 0; row1 < row2; ++ row1, -- row2) {
            std::swap(image[pixel1 ++], image[pixel2 --]);
        }
    }
}

void LaserProfile::swapXY(unsigned short *image) {
    unsigned short *copy = new unsigned short [sizeX_m * sizeY_m];
    for (unsigned int col = 0; col < sizeX_m; ++ col) {
        for (unsigned int row = 0; row < sizeY_m; ++ row) {
            unsigned int pixel1 = col * sizeY_m + row;
            unsigned int pixel2 = row * sizeX_m + col;

            copy[pixel2] = image[pixel1];
        }
    }

    for (unsigned int pixel = 0; pixel < sizeX_m * sizeY_m; ++ pixel) {
        image[pixel] = copy[pixel];
    }

    delete[] copy;

    std::swap(sizeX_m, sizeY_m);
}

void LaserProfile::filterSpikes(unsigned short *image) {
    hsize_t idx = sizeY_m + 1;
    hsize_t idxN = idx - sizeY_m;
    hsize_t idxNW = idxN - 1, idxNE = idxN + 1;
    hsize_t idxW = idx - 1, idxE = idx + 1;
    hsize_t idxS = idx + sizeY_m;
    hsize_t idxSW = idxS - 1, idxSE = idxS + 1;

    for (hsize_t i = 1; i < sizeX_m - 1; ++ i) {
        for (hsize_t j = 1; j < sizeY_m - 1; ++ j) {

            if (image[idx] > 0) {
                if (image[idxNW] == 0 &&
                    image[idxN] == 0 &&
                    image[idxNE] == 0 &&
                    image[idxW] == 0 &&
                    image[idxE] == 0 &&
                    image[idxSW] == 0 &&
                    image[idxS] == 0 &&
                    image[idxSE] == 0) {

                    image[idx] = 0;
                }
            }

            ++ idxNW; ++ idxN; ++ idxNE;
            ++ idxW;  ++ idx;  ++ idxE;
            ++ idxSW; ++ idxS; ++ idxSE;
        }
    }
}

void LaserProfile::normalizeProfileData(double intensityCut, unsigned short *image) {
    unsigned short int profileMax = getProfileMax(image);

    unsigned int pixel = 0;
    for (unsigned int col = 0; col < sizeX_m; ++ col) {
        for(unsigned int row = 0; row < sizeY_m; ++ row, ++ pixel) {
            double val = (double(image[pixel]) / profileMax - intensityCut) / (1.0 - intensityCut);

            val = std::max(0.0, val);
            image[pixel] = std::round(val * profileMax);
        }
    }
}

void LaserProfile::computeProfileStatistics(unsigned short *image) {
    double totalMass = 0.0;
    centerMass_m = Vector_t(0.0);
    standardDeviation_m = Vector_t(0.0);

    unsigned int pixel = 0;
    for(unsigned int col = 0; col < sizeX_m; ++ col) {
        for(unsigned int row = 0; row < sizeY_m; ++ row, ++ pixel) {
            double val = image[pixel];

            centerMass_m(0) += val * col;
            centerMass_m(1) += val * row;
            standardDeviation_m(0) += val * col*col;
            standardDeviation_m(1) += val * row*row;
            totalMass += val;
        }
    }

    centerMass_m /= totalMass;
    standardDeviation_m(0) = sqrt(standardDeviation_m(0) / totalMass -
                                  centerMass_m(0)*centerMass_m(0));
    standardDeviation_m(1) = sqrt(standardDeviation_m(1) / totalMass -
                                  centerMass_m(1)*centerMass_m(1));
}

void LaserProfile::fillHistrogram(unsigned short *image) {
    unsigned int histSizeX = sizeX_m;
    unsigned int histSizeY = sizeY_m;
    double histRangeX = histSizeX / standardDeviation_m(0) / 2;
    double histRangeY = histSizeY / standardDeviation_m(1) / 2;
    double binSizeX = 1.0 / standardDeviation_m(0);
    double binSizeY = 1.0 / standardDeviation_m(1);

    hist2d_m = gsl_histogram2d_alloc(histSizeX, histSizeY);
    gsl_histogram2d_set_ranges_uniform(hist2d_m,
                                       -histRangeX, histRangeX,
                                       -histRangeY, histRangeY);


    unsigned int pixel = 0;
    for (unsigned int col = 0; col < sizeX_m; ++ col) {
        double x = (col - 0.5 * sizeX_m) / standardDeviation_m(0);
        x = x + 0.5 * binSizeX;
        // std::cout << "x = " << x << std::endl;

        for (unsigned int row = 0; row < sizeY_m; ++ row, ++ pixel) {
            double val = image[pixel];
            double y = -(row - 0.5 * sizeY_m) / standardDeviation_m(1);
            y = y + 0.5 * binSizeY;
            // if (col == 0)
            //     std::cout << "y = " << y << std::endl;

            gsl_histogram2d_accumulate(hist2d_m, x, y, val);
        }
    }
    saveHistogram();
}

void LaserProfile::setupRNG() {
    unsigned int histSizeX = hist2d_m->nx, histSizeY = hist2d_m->ny;

    gsl_rng_env_setup();

    const gsl_rng_type *T = gsl_rng_default;
    rng_m = gsl_rng_alloc(T);

    pdf_m = gsl_histogram2d_pdf_alloc(histSizeX, histSizeY);
    gsl_histogram2d_pdf_init(pdf_m, hist2d_m);
}

void LaserProfile::printInfo() {
    INFOMSG(level3 << "* **********************************************************************************\n");
    INFOMSG("* LASER PROFILE \n");
    INFOMSG("* size = " << sizeX_m << " x " << sizeY_m << " pixels \n");
    INFOMSG("* center of mass: x = " << centerMass_m(0) << ", y = " << centerMass_m(1) << "\n");
    INFOMSG("* standard deviation: x = " << standardDeviation_m(0) << ", y = " << standardDeviation_m(1) << "\n");
    INFOMSG("* **********************************************************************************" << endl);
}

void LaserProfile::saveData(const std::string &fname, unsigned short *image) {
    std::ofstream out("data/" + fname + ".pgm");
    out << "P2" << std::endl;
    out << sizeX_m << " " << sizeY_m << std::endl;
    out << getProfileMax(image) << std::endl;

    for (unsigned int j = 0; j < sizeY_m; ++ j) {
        for (unsigned int i = 0; i < sizeX_m; ++ i) {
            out << image[i * sizeY_m + j] << " ";
        }
        out << std::endl;
    }
}

void LaserProfile::saveHistogram() {
    FILE  *fh = fopen("data/LaserHistogram.dat", "w");
    gsl_histogram2d_fprintf(fh, hist2d_m, "%g", "%g");
    fclose(fh);
}

void LaserProfile::sampleDist() {
    std::ofstream fh("data/LaserEmissionSampled.dat");
    double x, y;

    for(int i = 0; i < 1000000; i++) {
        getXY(x, y);
        fh << x << "\t" << y << "\n";
    }

    fh.close();
}

void LaserProfile::getXY(double &x, double &y) {
    double u = gsl_rng_uniform(rng_m);
    double v = gsl_rng_uniform(rng_m);
    gsl_histogram2d_pdf_sample(pdf_m, u, v, &x, &y);
}

unsigned short LaserProfile::getProfileMax(unsigned short *image) {
    unsigned int numberPixels = sizeX_m * sizeY_m;

    unsigned short maxIntensity = 0;

    for(unsigned int i = 0; i < numberPixels; i++) {
        if(image[i] > maxIntensity)
            maxIntensity = image[i];
    }

    return maxIntensity;
}
