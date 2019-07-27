/*

Purpose: Convert ANSIS E & B-Field data into H5hut (H5block)
         format for usage in OPAL.

Usage: ascii2h5block efield.txt hfield.txt ehfield

To visulize use Visit: https://wci.llnl.gov/codes/visit/

Ch. Wang & A. Adelmann, 2011
D. Winklehner, 2013

ToDo: make it more generic / duh -DW
static2: Changed it to reflect new field files from Daniela on Sep. 15 2014
static2a: 0-Hfield for IsoDAR central region RF
 
*/


#include <fstream>
#include <ios>
#include <iostream>
#include "H5hut.h"
#include <cassert>
#include <string>
#include <algorithm>
#include <cmath>


int main(int argc,char *argv[]) {
    
    if (argc != 3) {
        std::cout << "Wrong number of arguments: ascii2h5block efield.txt outfield" << std::endl;
        exit(1);
    }

    std::string efin(argv[1]);
    std::string ehfout(argv[2]);
    std::string ehfout_c = ehfout + std::string("_CYC.h5part");

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Using " << efin << " to create " << ehfout_c << std::endl;

    // Open file streams
    std::ifstream finE;
    finE.open(efin.c_str());
    assert(finE.is_open());

    // Get number of lines
    int nlinesE = std::count(std::istreambuf_iterator<char>(finE), 
			     std::istreambuf_iterator<char>(), '\n');
    
    std::cout << "Lines in finE: " << nlinesE << std::endl;

    // Header has 5 lines
    int nlines = nlinesE - 5;

    // Reset iterator
    finE.seekg(0, finE.beg);
    
    std::string templine;

    // Skip the 5 header lines
    for (int i = 0; i < 5; i++){
      std::getline(finE, templine);
    }

    /* Set frequency (TODO: ask AA if this is the right thing 
    to do for static fields -DW) */
    h5_float64_t freq = 49.2e6; //49.2 MHz ->Hz

    h5_float64_t *FieldstrengthEz = new h5_float64_t[nlines]; 
    h5_float64_t *FieldstrengthEx = new h5_float64_t[nlines]; 
    h5_float64_t *FieldstrengthEy = new h5_float64_t[nlines]; 
    h5_float64_t *FieldstrengthHz = new h5_float64_t[nlines];  
    h5_float64_t *FieldstrengthHx = new h5_float64_t[nlines];
    h5_float64_t *FieldstrengthHy = new h5_float64_t[nlines];

    // Daniela's latest files don't have the x,y,z coordinates anymore
    // Got the following from the readme file accompanying the fields:
    h5_float64_t xbegin = -90.0;
    h5_float64_t xend   = -10.0;
    h5_float64_t ybegin = -30.0;
    h5_float64_t yend   = 30.0;
    h5_float64_t zbegin = -1.0;
    h5_float64_t zend   = 1.0;

    // Init the arrays for fields
    double *Ex = new double[nlines];
    double *Ey = new double[nlines];
    double *Ez = new double[nlines];
    
    /*N.B.: Daniela's files now have the structure: Bx,By,Bz; no complex numbers
     units are cm, V/cm and Gauss */
    for(int i = 0; i < nlines; i++) {
        finE >> Ex[i] >> Ey[i] >> Ez[i];
    }

    finE.close();

    double Emax = 0.0;
    double E_temp;
    int i_temp = 0;

    for (int i = 0; i < nlines; i++){
        E_temp = std::sqrt(Ex[i] * Ex[i] + Ey[i] * Ey[i] + Ez[i] * Ez[i]);
        if (E_temp > Emax) {
            Emax = E_temp;
            i_temp = i;
        }
    }

    std::cout << "Hardcoded limits: x(" << xbegin << "/" << xend << ") cm" << std::endl;
    std::cout << "Hardcoded limits: y(" << ybegin << "/" << yend << ") cm" << std::endl;
    std::cout << "Hardcoded limits: z(" << zbegin << "/" << zend << ") cm" << std::endl;

    // Set spacing 
    // TODO: Make program find spacing automatically -DW
    double spacing = 0.2;

    std::cout << "Hardcoded spacing: " << spacing << " cm" << std::endl;

    double gridPx_temp = (xend - xbegin) / spacing + 1.0;
    double gridPy_temp = (yend - ybegin) / spacing + 1.0;
    double gridPz_temp = (zend - zbegin) / spacing + 1.0;

    int gridPx = (int) gridPx_temp;
    int gridPy = (int) gridPy_temp;
    int gridPz = (int) gridPz_temp;

    int nlines_hc = gridPx * gridPy * gridPz;

    std::cout << "Hardcoded nlines: " << nlines_hc << std::endl;
    std::cout << "File nlines: " << nlines << std::endl; 

    std::cout << "Grid dimensions: Px = " << gridPx << " , Py = " << gridPy << " , Pz = " << gridPz << std::endl;

    std::cout << "E_max = (" << Ex[i_temp] << ", "<< Ey[i_temp] << ", " << Ez[i_temp] << ") V/cm at index " << i_temp << "." << std::endl;

    std::cout << "Converting from V/cm and cm to kV/mm and mm before saving h5part" << std::endl;


    // Here we also convert from from G to kG
    for (int i = 0; i < gridPz; i++) {
      for (int j = 0; j < gridPy; j++) {
	for (int k = 0; k < gridPx; k++) {
                FieldstrengthEx[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Ex[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthEy[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Ey[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthEz[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Ez[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthHx[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ex[i+j*gridPz+k*gridPz*gridPy])*1e-3;
		     FieldstrengthHy[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ey[i+j*gridPz+k*gridPz*gridPy])*1e-3;
		     FieldstrengthHz[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ez[i+j*gridPz+k*gridPz*gridPy])*1e-3;
            }
        }
    }

/*
    // Here we also convert from V/cm to kV/mm (MV/m)
    for (int i = 0; i < gridPz; i++) {
      for (int j = 0; j < gridPy; j++) {
	for (int k = 0; k < gridPx; k++) {
                FieldstrengthEx[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ex[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthEy[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ey[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthEz[k+j*gridPx+i*gridPx*gridPy]=static_cast<h5_float64_t>(Ez[i+j*gridPz+k*gridPz*gridPy])*1e-4;
                FieldstrengthHx[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Hx[i+j*gridPz+k*gridPz*gridPy])*1e-3;
		     FieldstrengthHy[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Hy[i+j*gridPz+k*gridPz*gridPy])*1e-3;
		     FieldstrengthHz[k+j*gridPx+i*gridPx*gridPy]=0.0; //static_cast<h5_float64_t>(Hz[i+j*gridPz+k*gridPz*gridPy])*1e-3;
            }
        }
    }
*/
    // Change spacing and limits from cm to mm
    spacing *= 10.0;
    xbegin *= 10; ybegin *= 10; zbegin *= 10;
    xend *= 10; yend *= 10; zend *= 10;

    // Write h5part file for OPAL-CYC
    h5_err_t h5err;
    h5_file_t file = H5OpenFile(ehfout_c.c_str(), H5_O_WRONLY, H5_PROP_DEFAULT);
    h5err = H5Block3dSetView(file,
                             0, gridPx - 1,
                             0, gridPy - 1,
                             0, gridPz - 1);
    if(file) {

        H5SetStep(file, 0);

        H5Block3dWriteVector3dFieldFloat64 (
	    file,            /*!< IN: file handle */
	    "Efield",        /*!< IN: name of dataset to write */
	    FieldstrengthEx, /*!< IN: X axis data */
	    FieldstrengthEy, /*!< IN: Y axis data */
	    FieldstrengthEz  /*!< IN: Z axis data */
	);

        h5err = H5Block3dSetFieldSpacing(file, "Efield", spacing, spacing, spacing);
        h5err = H5Block3dSetFieldOrigin(file, "Efield", xbegin, ybegin, zbegin);

        H5Block3dWriteVector3dFieldFloat64 (
	    file,            /*!< IN: file handle */
	    "Hfield",	     /*!< IN: name of dataset to write */
	    FieldstrengthHx, /*!< IN: X axis data */
	    FieldstrengthHy, /*!< IN: Y axis data */
	    FieldstrengthHz  /*!< IN: Z axis data */
	);

        h5err = H5Block3dSetFieldSpacing(file, "Hfield", spacing, spacing, spacing);
        h5err = H5Block3dSetFieldOrigin(file, "Hfield", xbegin, ybegin, zbegin);

	// Frequency here really in Hz ??? -DW
        H5WriteFileAttribFloat64 (
            file,            /*!< [in] Handle to open file */
            "Resonance Frequency(Hz)", /*!< [in] Name of attribute */
            &freq,           /*!< [in] Array of attribute values */ 
            1	             /*!< [in] Number of array elements */
        );

        H5CloseFile(file);
    } 

    /*
    // Write text file for OPAL-T
    std::ofstream foutEH;
    foutEH.open(ehfout_t.c_str());
    assert(foutEH.is_open());

    // Change spacing and limits to cm
    spacing = spacing * 1.0e-1;
    xbegin *= 1.0e-1; ybegin *= 1.0e-1; zbegin *= 1.0e-1;
    xend *= 1.0e-1; yend *= 1.0e-1; zend *= 1.0e-1;

    // Header
    foutEH << "3DDynamic\tXYZ\n";
    foutEH << freq*1e-6 << "\n"; // frequency in MHz
    foutEH << xbegin << "\t" << xend << "\t" << gridPx-1 << "\n"; //cm / cm / #
    foutEH << ybegin << "\t" << yend << "\t" << gridPy-1 << "\n"; //cm / cm / #
    foutEH << zbegin << "\t" << zend << "\t" << gridPz-1 << "\n"; //cm / cm / #
 
    // Fielddata
    // Here we also convert from V/cm to kV/mm (MV/m) and from G to T (NB the difference to the CYC field!)
    foutEH.setf(ios::scientific,ios::floatfield);
    foutEH.precision(5);

    for (int i=0; i<nlines; i++) {
      foutEH << Ex[i]*1e-4 << "\t" << Ey[i]*1e-4 << "\t" << Ez[i]*1e-4 << "\t";
      foutEH << Hx[i]*1e-4 << "\t" << Hy[i]*1e-4 << "\t" << Hz[i]*1e-4 << "\n";
    }

    foutEH.close();
    */

    std::cout << "Done, bye ..." << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}
