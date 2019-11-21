# This is OPAL 2.

The documentation has been moved to the [Wiki](https://gitlab.psi.ch/OPAL/src/wikis/home).

# OPAL-FEL
A new undulator/wiggler element is added to OPAL-t. When a particle bunch approaches it, the solver will change to a full-wave solver that will simulate the run through the undulator.


Based on [MITHRA](https://github.com/aryafallahi/mithra) full-wave Maxwell code for FELs, by Arya Fallahi.

## To compile OPAL-FEL

* Clone repository `git@gitlab.psi.ch:albajacas_a/mithra.git`
* Set environment variable
  `MITHRA_INCLUDE_DIR=/path/to/mithra`
* Build OPAL with `cmake -DENABLE_OPAL_FEL=yes ..`

## To use UNDULATOR element
Example to instantiate an undulator/wiggler
```  // Parameters for full-wave simulation of an undulator/wiggler, based on MITHRA FEL solver by Arya Fallahi.
   // The full-wave solver will start at ELEMEDGE, but the undulator will start at ELEMEDGE + LFRINGE
   // L:  	    	   LFRINGE + undulator length
   // LFRINGE:		   Distance in front of the undulator at which the full-wave solver starts
   // LAMBDA:		   Undulator period
   // K:	   	   Undulator parameter, ~93.37*B0*LAMBDA, with B0 maximum magnetic field
   // MESHLENGTH:  	   Mesh size for full-wave simulation
   // MESHRESOLUTION:	   Mesh discretisation
   // TRUNORDER:	   Truncation order of absorbing boundaries. Can be 1 or 2
   // SPACECHARGE:	   Boolean to include space-charge in undulator simulation
   // FNAME:		   File to indicate desired output from undulator simulation
   // TOTALTIME:	   Optional. Indicate total time of full-wave simulation
   // Undulator of length 85cm + an entrance fringe of 50cm
   UND: UNDULATOR, L = 1.35, LFRINGE = .5, ELEMEDGE = und_edge, K = 11.1, LAMBDA = 8.5e-2,
   	  MESHLENGTH = { 5e-3, 5e-3, 18e-3 }, MESHRESOLUTION = { .1e-3, .1e-3, 2e-6},
     	  TRUNORDER = 2, SPACECHARGE = 1, 
	  FNAME = "undulator_output.job", TOTALTIME = 8.1e-9;
```

Example files can be found in [/examples-OPAL-FEL](https://gitlab.psi.ch/albajacas_a/src/blob/master/examples-OPAL-FEL)