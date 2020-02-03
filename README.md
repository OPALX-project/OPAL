# This is OPAL 2.

The documentation has been moved to the [Wiki](https://gitlab.psi.ch/OPAL/src/wikis/home).

# OPAL-FEL
A new undulator/wiggler element is added to OPAL-t. When a particle bunch approaches it, the solver will change to a full-wave solver that will simulate the run through the undulator.


Based on [MITHRA](https://github.com/aryafallahi/mithra) full-wave Maxwell code for FELs, by Arya Fallahi.

## To compile OPAL-FEL

* Clone repository `git@github.com:arnaualba/mithra.git`
* Set environment variables
```
export MITHRA_INCLUDE_DIR=$HOME/data/githubMithra/mithra
export PATH=$MITHRA_INCLUDE_DIR:$PATH
```

* Build the MITHRA library with 
```
make lib
```

* Build OPAL with 
```
cmake -DENABLE_OPAL_FEL=yes ..
make
```

## To use UNDULATOR element
Example to instantiate an undulator/wiggler
```
	// Parameters for full-wave simulation of an undulator/wiggler, based on MITHRA FEL solver by Arya Fallahi.
   // The full-wave solver will start at ELEMEDGE, but the undulator will start at ELEMEDGE + LFRINGE
   // L:  	    	   LFRINGE + undulator length
   // LFRINGE:		   Distance in front of the undulator at which the full-wave solver starts
   // LAMBDA:		   Undulator period
   // K:	   	   Undulator parameter, ~93.37*B0*LAMBDA, with B0 maximum magnetic field
   // MESHLENGTH:  	   Mesh size for full-wave simulation
   // MESHRESOLUTION:	   Mesh discretisation
   // TRUNORDER:	   Truncation order of absorbing boundaries. Can be 1 or 2
   // SPACECHARGE:	   Boolean to include space-charge in undulator simulation
   // EMITPARTICLES:   Boolean to emit particles as time passes, rather than injecting the whole bunch at once
   // FNAME:		   File to indicate desired output from undulator simulation
   // TOTALTIME:	   Total time of full-wave simulation
   // Undulator of length 85cm + an entrance fringe of 50cm
   UND: UNDULATOR, L = 1.35, LFRINGE = .5, ELEMEDGE = und_edge, K = 11.1, LAMBDA = 8.5e-2,
	MESHLENGTH = { 5e-3, 5e-3, 18e-3 }, MESHRESOLUTION = { .1e-3, .1e-3, 2e-6},
	TRUNORDER = 2, SPACECHARGE = 1, EMITPARTICLES = 1,
	FNAME = "output_parameters.job", TOTALTIME = 8.1e-9;
```

## Output of full wave simulation
The details of the output produced from the full-wave simulation can be found in the manual [Manual](https://github.com/arnaualba/mithra/tree/master/doc), and can be viewd with the Python functions [plotMithraPy](https://github.com/arnaualba/plotMithraPy).
Here we provide an example of the output_parametes.job file:
```
BUNCH
{ 
  bunch-sampling
  {
    sample							= true
    directory						= ./
    base-name 						= bunch-sampling/bunch
    rhythm							= 4e-11
  }
    
  bunch-visualization
  {
    sample							= false
    directory						= ./
    base-name						= bunch-visualization/bunch
    rhythm							= 5.0
  }

  bunch-profile
  {
    sample							= false
    directory						= ./
    base-name						= bunch-profile/bunch
    rhythm							= 1e-9
  }
}


FIELD
{
  field-sampling
  {
    sample							= false
    type							= at-point
    field							= Ex
    field                  			= Ey
    field                 			= Ez
    directory						= ./
    base-name						= field-sampling/field
    rhythm							= 3.2e-12
    position						= (0.0, 0.0, 110.0e-6)
#    line-begin						= (0.0, 0.0, 100.0)
#    line-end						= (0.0, 0.0, 130.0)
#    number-of-points				= 100
  }
  
  field-visualization
  {
    sample							= false
    type							= in-plane
    field							= Ex
    field							= Ey
    field							= Ez
    plane							= yz
    position						= (0.0, 0.0, 110.0)
    directory						= ./
    base-name						= field-visualization/field1
    rhythm							= 5.0
  }
    
  field-profile
  {
    sample							= false
    field                   		= Q
    directory						= ./
    base-name						= field-profile/field
    rhythm							= 100
  }
}


FEL-OUTPUT
{
  radiation-power
  {
    sample							= true
    type							= at-point
    directory						= ./
    base-name						= power-sampling/power-NSFD
    plane-position					= 8e-3
    plane-position					= 4e-3
    normalized-frequency			= 1.00
  }
  screen-profile
  {
    sample							= true
    directory						= ./
    base-name						= screen-profile/bunch
    position						= -.4
    position						= .88
  }
}
```