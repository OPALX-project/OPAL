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
An example to run and test can be found in the folder `src/examples-OPAL-FEL/`
Here is an example on ho to instantiate an undulator/wiggler
```
    // Parameters for full-wave simulation of an undulator/wiggler, based on MITHRA FEL solver by Arya Fallahi.
    // The full-wave solver will start at ELEMEDGE.
    //The full undulator length is LAMBDA * (NUMPERIODS + 4) because fringes measure 2*LAMBDA
    // ELEMEDGE:        Start of undulator fringe, and transfer to FW solver
    // LAMBDA:          Undulator period
    // NUMPERIODS       Number of periods in the undulator
    // K:               Undulator parameter, ~93.37*B0*LAMBDA, with B0 maximum magnetic field
    // MESHLENGTH:      Mesh size for FW solver
    // MESHRESOLUTION:  Mesh discretisation for FW solver
    // TRUNORDER:       Truncation order of absorbing boundaries. Can be 1 or 2
    // FNAME:          File to indicate desired output from undulator simulation
    // TOTALTIME:      Total time of full-wave simulation
    // Undulator of length 2.1m + 70cm fringe field
    UND: UNDULATOR, ELEMEDGE = .01, K = 51.5, LAMBDA = .35, NUMPERIODS = 6,
        MESHLENGTH = { 1.8e-3, 4e-3, .1e-3 }, MESHRESOLUTION = {10e-6, 10e-6, .07e-6},
        TRUNORDER = 2, TOTALTIME = 1.4e-8, FNAME = "output.job";
	
```

NOTE: At the moment the simulation stops right after the FW simulation is done running.
In the future the particles will be transferred back from the FW solver to the static solver, and the simulation will continue.

## Output of full wave simulation
The details of the output produced from the full-wave simulation can be found in the manual [Manual](https://github.com/arnaualba/mithra/tree/master/doc), and can be viewed with the Python functions [plotMithraPy](https://github.com/arnaualba/plotMithraPy).
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
  bunch-profile-lab-frame
  {
    sample							= true
    directory						= ./
    base-name						= screen-profile/bunch
    position						= -.4
    position						= .88
  }
}
```