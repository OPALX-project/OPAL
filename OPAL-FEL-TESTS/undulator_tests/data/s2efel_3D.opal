OPTION, PSDUMPFREQ = 1000;
   OPTION, STATDUMPFREQ = 1000;
  OPTION, BOUNDPDESTROYFQ=10;
 OPTION, AUTOPHASE=4;
        
OPTION, VERSION=10900;

REAL REPARTFREQ=1000000;

Title, string="UNDULATOR";

REAL rf_freq             = 1.0;
     REAL n_particles         = 1224;
      REAL beam_bunch_charge   =  1.6e-19 * 1e6 * 1.846e8;
      
REAL gamma   = 100.41;
 REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;
    
value , {gamma, P0};

UND: UNDULATOR, L = 9.0, , K = 1.417, LAMBDA = .03,      MESHLENGTH = { 3.2e-3, 3.2e-3, .28e-3 }, MESHRESOLUTION = {100e-6, 100e-6, .4e-6},     TRUNORDER = 2, SPACECHARGE = 1,     FNAME = "undulator_output.job",     TOTALTIME = 3.003e-8, TIMESTEPRATIO = 5;

DRIVE: Line = (UND);

Dist: DISTRIBUTION, TYPE = FROMFILE,      	FNAME = "initial-profile-for-OPAL.tsv",	EMITTED = FALSE;

FS_SC: Fieldsolver, FSTYPE = FFT,             MX = 8, MY = 8, MT = 8,             PARFFTX = true,            PARFFTY = true,             PARFFTT = true,              BCFFTX = open,             BCFFTY = open,             BCFFTT = open,            BBOXINCR = 1,             GREENSF = INTEGRATED;

REAL mithra_time = 3e-8;
BEAM1:  BEAM, PARTICLE = ELECTRON, pc = P0, NPART = n_particles,	BFREQ = 1 / mithra_time, BCURRENT = beam_bunch_charge / mithra_time , CHARGE = -1;

TRACK, LINE = DRIVE, BEAM = BEAM1, MAXSTEPS = 1000000, DT = 3.E-10, ZSTOP= 10.0 ;
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,     FIELDSOLVER = FS_SC, DISTRIBUTION = Dist;
ENDTRACK;

Quit;

