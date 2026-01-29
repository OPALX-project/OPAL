/*

  Purpose: Convert ANSYS E & B-Field data into H5hut (H5block)
  format for usage in OPAL.

  Usage: ascii2h5block efield.txt hfield.txt ehfout

  To visualize use Visit: https://wci.llnl.gov/codes/visit/

  Ch. Wang & A. Adelmann, 2011

  ToDo: make it more generic

  Modification by Chris van Herwaarden and Hui Zhang

  The first three rows of a field map that you wish to combine should look like this:

  int1 int2 int3
  int1 int2 int3
  int1 int2 int3

  the integers are the amount of steps (or different values) in x y and z

*/


#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "H5hut.h"

int main(int argc,char *argv[]) {

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

    h5_float64_t freq = 72.615 * 1.0e6;

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Combine " << efin << " and " << hfin << " to " << ehfout << std::endl;

    std::cout << "Frequency " << freq << " [Hz]" << std::endl;

    std::ifstream finE, finH;
    finH.open(hfin);
    finE.open(efin);
    bool efield = false, hfield = false;
    if (finE.is_open())
        efield = true;
    else if (efin.empty() == false) {
        std::cout << "E-field \"" << efin << "\" could not be opened" << std::endl;
        std::exit(1);
    }
    if (finH.is_open())
        hfield = true;
    else if (hfin.empty() == false) {
        std::cout << "H-field \"" << hfin << "\" could not be opened" << std::endl;
        std::exit(1);
    }

    if (efield == false && hfield == false) {
        std::cerr << "Neither E-field \"" << efin
                  << "\" nor H-field \""  << hfin
                  << "\" could be opened" << std::endl;
        std::exit(1);
    }

    int gridPx = 0, gridPy = 0, gridPz = 0;
    char temp[256];
    /* Header and grid info */
    /* Deletes the first two rows in finE and finH to get rid of header info*/
    if (efield) {
        finE.getline(temp,256);
        finE >> gridPx >> gridPy >> gridPz;
        finE.getline(temp,256);
        finE.getline(temp,256);
    }

    int HgridPx = 0, HgridPy = 0, HgridPz = 0;
    if (hfield) {
        finH.getline(temp,256);
        finH >> HgridPx >> HgridPy >> HgridPz;
        finH.getline(temp,256);
        finH.getline(temp,256);
    }

    int n =  gridPx *  gridPy *  gridPz; /* number of rows in a column */
    int m = HgridPx * HgridPy * HgridPz;

    int H5gridPx = std::max(gridPx, HgridPx);
    int H5gridPy = std::max(gridPy, HgridPy);
    int H5gridPz = std::max(gridPz, HgridPz);

    std::cout << "H5block grid" << std::endl;
    std::cout << "H5gridPx " << H5gridPx << std::endl;
    std::cout << "H5gridPy " << H5gridPy << std::endl;
    std::cout << "H5gridPz " << H5gridPz << std::endl;

    h5_file_t file = H5OpenFile(ehfout.c_str(), H5_O_WRONLY, H5_PROP_DEFAULT);
    if (!file) {
        std::cerr << "Could not open output file " << ehfout << std::endl;
        std::exit(1);
    }

    H5SetStep(file, 0);
    H5Block3dSetView(file,
                     0, H5gridPx - 1,
                     0, H5gridPy - 1,
                     0, H5gridPz - 1);

    std::cout << "number Edata " << n << std::endl;
    if (efield) {
        h5_float64_t* sEx = new h5_float64_t[n]; /* redefines the sE and sH variables as arrays */
        h5_float64_t* sEy = new h5_float64_t[n];
        h5_float64_t* sEz = new h5_float64_t[n];

        h5_float64_t* FieldstrengthEz = new h5_float64_t[n]; /* redefines the fieldstrength variables as arrays */
        h5_float64_t* FieldstrengthEx = new h5_float64_t[n];
        h5_float64_t* FieldstrengthEy = new h5_float64_t[n];

        double* Ex = new double[n]; /* redefines the E and H variables as arrays */
        double* Ey = new double[n];
        double* Ez = new double[n];

        for (int i = 0; i < n; i++) {
            finE >> sEx[i] >> sEy[i] >> sEz[i] >> Ex[i] >> Ey[i] >> Ez[i];
        }
        finE.close();

        h5_float64_t stepEx = (sEx[n-1] - sEx[0]) / (gridPx - 1); /* calculates the stepsizes of the x,y,z of the efield and hfield*/
        h5_float64_t stepEy = (sEy[n-1] - sEy[0]) / (gridPy - 1);
        h5_float64_t stepEz = (sEz[n-1] - sEz[0]) / (gridPz - 1);

        std::cout << "gridPx " << gridPx << " stepEx " << stepEx << std::endl;
        std::cout << "gridPy " << gridPy << " stepEy " << stepEy << std::endl;
        std::cout << "gridPz " << gridPz << " stepEz " << stepEz << std::endl;

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
        H5Block3dSetFieldSpacing(file, "Efield", stepEx, stepEy, stepEz);
        H5Block3dSetFieldOrigin (file, "Efield", sEx[0], sEy[0], sEz[0]);
    }

    std::cout << "number Bdata " << m << std::endl;

    if (hfield) {

        h5_float64_t* sHx = new h5_float64_t[m];
        h5_float64_t* sHy = new h5_float64_t[m];
        h5_float64_t* sHz = new h5_float64_t[m];

        h5_float64_t* FieldstrengthHz = new h5_float64_t[m];
        h5_float64_t* FieldstrengthHx = new h5_float64_t[m];
        h5_float64_t* FieldstrengthHy = new h5_float64_t[m];

        double* Hx = new double[m];
        double* Hy = new double[m];
        double* Hz = new double[m];

        for (int i = 0; i < m; i++) {
            finH >> sHx[i] >> sHy[i] >> sHz[i] >> Hx[i] >> Hy[i] >> Hz[i];
        }
        finH.close();

        h5_float64_t stepHx = (sHx[m-1] - sHx[0]) / (HgridPx - 1);
        h5_float64_t stepHy = (sHy[m-1] - sHy[0]) / (HgridPy - 1);
        h5_float64_t stepHz = (sHz[m-1] - sHz[0]) / (HgridPz - 1);

        std::cout << "HgridPx " << HgridPx << " stepHx " << stepHx << std::endl;
        std::cout << "HgridPy " << HgridPy << " stepHy " << stepHy << std::endl;
        std::cout << "HgridPz " << HgridPz << " stepHz " << stepHz << std::endl;

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
        H5Block3dSetFieldSpacing(file, "Hfield", stepHx, stepHy, stepHz);
        H5Block3dSetFieldOrigin (file, "Hfield", sHx[0], sHy[0], sHz[0]);
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