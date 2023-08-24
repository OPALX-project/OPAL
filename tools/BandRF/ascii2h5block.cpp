//
// ascii2h5block tool
//
// Copyright (c) 2011 - 2023, Ch. Wang, A. Adelmann, Achim Gsell, Jochem Snuverink,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
//                            Daniel Winklehner, MIT, Cambridge, MA, USA
//                            Chris van Herwaarden and Hui Zhang
//                            Pedro Calvo, CIEMAT, Spain
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "H5hut.h"

#include "Physics/Units.h"

/*
  Purpose: Convert ASCII E & B-Field data into H5hut (H5block)
  format for usage in OPAL.

  Usage: ascii2h5block efield.txt hfield.txt ehfout

  To visualize use Visit: https://wci.llnl.gov/codes/visit/

  The first three rows of a field map that you wish to combine should look like this:

  int1 int2 int3

  the integers are the amount of steps (or different values) in x y and z

*/
int main(int argc,char* argv[]) {

    if (argc != 4) {
        std::cout << "Wrong number of arguments: ascii2h5block efield.txt (or \"\") hfield.txt (or \"\")  ehfout" << std::endl;
        //--commlib mpi" << std::endl;
        std::exit(1);
    }

    // // initialize MPI & H5hut
    // MPI_Init (&argc, &argv);
    // MPI_Comm comm = MPI_COMM_WORLD;
    // int comm_size = 1;
    // MPI_Comm_size (comm, &comm_size);
    // int comm_rank = 0;
    // MPI_Comm_rank (comm, &comm_rank);
    // H5AbortOnError ();
    // H5SetVerbosityLevel (h5_verbosity);

    std::string efin(argv[1]);
    std::string hfin(argv[2]);
    std::string ehfout(argv[3]);
    ehfout += std::string(".h5part");

    h5_float64_t freq = 72.615 * Units::MHz2Hz;

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Combine " << efin << " and " << hfin << " to " << ehfout << std::endl;
    std::cout << "Frequency " << freq << " [Hz]" << std::endl;

    std::ifstream finE, finH;
    finH.open(hfin);
    finE.open(efin);
    bool efield = false, hfield = false;
    if (finE.is_open()) {
        efield = true;
    } else if (efin.empty() == false) {
        std::cout << "E-field \"" << efin << "\" could not be opened" << std::endl;
        std::exit(1);
    }
    if (finH.is_open()) {
        hfield = true;
    } else if (hfin.empty() == false) {
        std::cout << "H-field \"" << hfin << "\" could not be opened" << std::endl;
        std::exit(1);
    }

    if (efield == false && hfield == false) {
        std::cerr << "Neither E-field \"" << efin
                  << "\" nor H-field \""  << hfin
                  << "\" could be opened" << std::endl;
        std::exit(1);
    }


    int linesE = std::count(std::istreambuf_iterator<char>(finE),
                 std::istreambuf_iterator<char>(), '\n');
    linesE = linesE - 1;
    finE.seekg(0, finE.beg); // Reset iterator

    int linesH = std::count(std::istreambuf_iterator<char>(finH),
                 std::istreambuf_iterator<char>(), '\n');
    linesH = linesH - 1;
    finH.seekg(0, finH.beg); // Reset iterator

    int gridPx = 0, gridPy = 0, gridPz = 0;
    int HgridPx = 0, HgridPy = 0, HgridPz = 0;
    char temp[256];
    /* Header and grid info */
    if (efield) {
        finE >> gridPx >> gridPy >> gridPz;
        finE.getline(temp,256);
    }
    if (hfield) {
        finH >> HgridPx >> HgridPy >> HgridPz;
        finH.getline(temp,256);
    }

    int ngridE =  gridPx *  gridPy *  gridPz; /* number of rows in a column */
    int ngridH = HgridPx * HgridPy * HgridPz;

    if (linesE != ngridE) {
        std::cerr << "The number of lines in the E-file (" << linesE
                  <<  ") doesn't match the number of lines specified by the grid ("
                  << ngridE << ")." << std::endl;
        std::exit(1);
    }
    if (linesH != ngridH) {
        std::cerr << "The number of lines in the H-file (" << linesH
                  <<  ") doesn't match the number of lines specified by the grid ("
                  << ngridH << ")." << std::endl;
        std::exit(1);
    }


    h5_file_t file = H5OpenFile(ehfout.c_str(), H5_O_WRONLY, H5_PROP_DEFAULT);
    if (!file) {
        std::cerr << "Could not open output file " << ehfout << std::endl;
        std::exit(1);
    }

    H5SetStep(file, 0);
    H5Block3dSetView(file,
                     0, gridPx - 1,
                     0, gridPy - 1,
                     0, gridPz - 1);

    if (efield) {
        std::cout << "Number of electric field data " << linesE << std::endl;

        h5_float64_t* sEx = new h5_float64_t[linesE]; /* redefines the sE and sH variables as arrays */
        h5_float64_t* sEy = new h5_float64_t[linesE];
        h5_float64_t* sEz = new h5_float64_t[linesE];

        h5_float64_t* FieldstrengthEz = new h5_float64_t[linesE]; /* redefines the fieldstrength variables as arrays */
        h5_float64_t* FieldstrengthEx = new h5_float64_t[linesE];
        h5_float64_t* FieldstrengthEy = new h5_float64_t[linesE];

        double* Ex = new double[linesE]; /* redefines the E and H variables as arrays */
        double* Ey = new double[linesE];
        double* Ez = new double[linesE];

        for (int i = 0; i < linesE; i++) {
            finE >> sEx[i] >> sEy[i] >> sEz[i] >> Ex[i] >> Ey[i] >> Ez[i];
        }
        finE.close();

        h5_float64_t stepEx = (sEx[linesE-1] - sEx[0]) / (gridPx - 1); /* calculates the stepsizes of the x,y,z of the efield and hfield*/
        h5_float64_t stepEy = (sEy[linesE-1] - sEy[0]) / (gridPy - 1);
        h5_float64_t stepEz = (sEz[linesE-1] - sEz[0]) / (gridPz - 1);

        if (stepEx < 0 || stepEy < 0 || stepEz < 0) {
            std::cerr << "Wrong step size of the efield!" << std::endl;
            std::exit(1);
        }
        std::cout << "gridPx = " << gridPx << " --- stepEx = " << stepEx << std::endl;
        std::cout << "gridPy = " << gridPy << " --- stepEy = " << stepEy << std::endl;
        std::cout << "gridPz = " << gridPz << " --- stepEz = " << stepEz << std::endl;
        std::cout << "sEx limits = (" << sEx[0] << ", " << sEx[linesE-1] << ") m" << std::endl;
        std::cout << "sEy limits = (" << sEy[0] << ", " << sEy[linesE-1] << ") m" << std::endl;
        std::cout << "sEz limits = (" << sEz[0] << ", " << sEz[linesE-1] << ") m" << std::endl;

        for (int i = 0; i < gridPz; i++) {
            for (int j = 0; j < gridPy; j++) {
                for (int k = 0; k < gridPx; k++) {
                    FieldstrengthEx[k + j * gridPx + i * gridPx * gridPy] = static_cast<h5_float64_t>(Ex[i + j * gridPz + k * gridPz * gridPy]);
                    FieldstrengthEy[k + j * gridPx + i * gridPx * gridPy] = static_cast<h5_float64_t>(Ey[i + j * gridPz + k * gridPz * gridPy]);
                    FieldstrengthEz[k + j * gridPx + i * gridPx * gridPy] = static_cast<h5_float64_t>(Ez[i + j * gridPz + k * gridPz * gridPy]);
                }
            }
        }

        H5Block3dWriteVector3dFieldFloat64 (
                                            file,            /*!< IN: file handle */
                                            "Efield",        /*!< IN: name of dataset to write */
                                            FieldstrengthEx, /*!< IN: X axis data */
                                            FieldstrengthEy, /*!< IN: Y axis data */
                                            FieldstrengthEz  /*!< IN: Z axis data */
                                            );
        H5Block3dSetFieldSpacing(file, "Efield", stepEx * Units::m2mm, stepEy * Units::m2mm, stepEz * Units::m2mm);
        H5Block3dSetFieldOrigin (file, "Efield", sEx[0] * Units::m2mm, sEy[0] * Units::m2mm, sEz[0] * Units::m2mm);
    }


    if (hfield) {
        std::cout << "Number of magnetic field data " << linesH << std::endl;

        h5_float64_t* sHx = new h5_float64_t[linesH];
        h5_float64_t* sHy = new h5_float64_t[linesH];
        h5_float64_t* sHz = new h5_float64_t[linesH];

        h5_float64_t* FieldstrengthHz = new h5_float64_t[linesH];
        h5_float64_t* FieldstrengthHx = new h5_float64_t[linesH];
        h5_float64_t* FieldstrengthHy = new h5_float64_t[linesH];

        double* Hx = new double[linesH];
        double* Hy = new double[linesH];
        double* Hz = new double[linesH];

        for (int i = 0; i < linesH; i++) {
            finH >> sHx[i] >> sHy[i] >> sHz[i] >> Hx[i] >> Hy[i] >> Hz[i];
        }
        finH.close();

        h5_float64_t stepHx = (sHx[linesH-1] - sHx[0]) / (HgridPx - 1);
        h5_float64_t stepHy = (sHy[linesH-1] - sHy[0]) / (HgridPy - 1);
        h5_float64_t stepHz = (sHz[linesH-1] - sHz[0]) / (HgridPz - 1);

        if (stepHx < 0 || stepHy < 0 || stepHz < 0) {
            std::cerr << "Wrong step size of the efield!" << std::endl;
            std::exit(1);
        }
        std::cout << "HgridPx = " << HgridPx << " --- stepHx = " << stepHx << std::endl;
        std::cout << "HgridPy = " << HgridPy << " --- stepHy = " << stepHy << std::endl;
        std::cout << "HgridPz = " << HgridPz << " --- stepHz = " << stepHz << std::endl;
        std::cout << "sHx limits = (" << sHx[0] << ", " << sHx[linesH-1] << ") m" << std::endl;
        std::cout << "sHy limits = (" << sHy[0] << ", " << sHy[linesH-1] << ") m" << std::endl;
        std::cout << "sHz limits = (" << sHz[0] << ", " << sHz[linesH-1] << ") m" << std::endl;

        for (int i = 0; i < HgridPz; i++) {
            for (int j = 0; j < HgridPy; j++) {
                for (int k = 0; k < HgridPx; k++) {
                    FieldstrengthHx[k + j * HgridPx + i * HgridPx * HgridPy] = static_cast<h5_float64_t>((Hx[i + j * HgridPz + k * HgridPz * HgridPy]));
                    FieldstrengthHy[k + j * HgridPx + i * HgridPx * HgridPy] = static_cast<h5_float64_t>((Hy[i + j * HgridPz + k * HgridPz * HgridPy]));
                    FieldstrengthHz[k + j * HgridPx + i * HgridPx * HgridPy] = static_cast<h5_float64_t>((Hz[i + j * HgridPz + k * HgridPz * HgridPy]));
                }
            }
        }
        H5Block3dWriteVector3dFieldFloat64 (
                                            file,            /*!< IN: file handle */
                                            "Hfield",        /*!< IN: name of dataset to write */
                                            FieldstrengthHx, /*!< IN: X axis data */
                                            FieldstrengthHy, /*!< IN: Y axis data */
                                            FieldstrengthHz  /*!< IN: Z axis data */
                                            );
        H5Block3dSetFieldSpacing(file, "Hfield", stepHx * Units::m2mm, stepHy * Units::m2mm, stepHz * Units::m2mm);
        H5Block3dSetFieldOrigin (file, "Hfield", sHx[0] * Units::m2mm, sHy[0] * Units::m2mm, sHz[0] * Units::m2mm);
    }

    H5WriteFileAttribFloat64 (
                              file,                      /*!< [in] Handle to open file */
                              "Resonance Frequency(Hz)", /*!< [in] Name of attribute */
                              &freq,                     /*!< [in] Array of attribute values */
                              1                          /*!< [in] Number of array elements */
                              );
    H5CloseFile(file);

    std::cout << "Done bye ..." << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}
