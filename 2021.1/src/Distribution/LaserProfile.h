#ifndef OPAL_LASERPROFILE_HH
#define OPAL_LASERPROFILE_HH
// ------------------------------------------------------------------------
// $RCSfile: LaserProfile.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LaserProfile
//
// ------------------------------------------------------------------------

#include "Algorithms/Vektor.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>
#include <string>
#include "hdf5.h"


class LaserProfile {

public:
    LaserProfile(const std::string &fileName,
                 const std::string &imageName,
                 double intensityCut,
                 short flags);

    ~LaserProfile();

    void getXY(double &x, double &y);

    enum {FLIPX = 1,
          FLIPY = 2,
          ROTATE90 = 4,
          ROTATE180 = 8,
          ROTATE270 = 16};
private:
    unsigned short * readFile(const std::string &fileName,
                              const std::string &imageName);
    unsigned short * readPGMFile(const std::string &fileName);
    unsigned short * readHDF5File(const std::string &fileName,
                                  const std::string &imageName);
    void flipX(unsigned short *image);
    void flipY(unsigned short *image);
    void swapXY(unsigned short *image);
    void filterSpikes(unsigned short *image);
    void normalizeProfileData(double intensityCut, unsigned short *image);
    void computeProfileStatistics(unsigned short *image);
    void fillHistrogram(unsigned short *image);
    void setupRNG();
    void printInfo();

    void saveData(const std::string &fname, unsigned short *image);
    void saveHistogram();
    void sampleDist();

    unsigned short getProfileMax(unsigned short *image);

    hsize_t sizeX_m, sizeY_m;
    gsl_histogram2d *hist2d_m;
    gsl_rng *rng_m;
    gsl_histogram2d_pdf *pdf_m;

    Vector_t centerMass_m;
    Vector_t standardDeviation_m;
};
#endif