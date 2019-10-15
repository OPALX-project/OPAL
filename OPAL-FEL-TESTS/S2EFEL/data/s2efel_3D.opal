OPTION, PSDUMPFREQ = 1000;
   OPTION, STATDUMPFREQ = 1000;
  OPTION, BOUNDPDESTROYFQ=10;
 OPTION, AUTOPHASE=4;
        
OPTION, VERSION=10900;

REAL REPARTFREQ=1000000;

Title, string="AWA Photoinjector + undulator";

REAL n_particles         = 1e6;
      REAL beam_bunch_charge   =  1.602e-19 * 1.846e08 * 1e6 ;
      
REAL gamma = 100.41;
REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;
    
value , {gamma, P0};

REAL gun_inj_phase = 0.0;
  
GUN:    RFCavity, L = 0.2927, VOLT = 60.0 * 6.6, , TYPE = "STANDING",        FMAPFN = "fieldmaps/DriveGun.T7", FREQ = 1300.0, LAG = (gun_inj_phase*Pi)/180.0;
 
REAL KSBF = 0.162544398;
if (OPALVERSION>15000) KSBF = KSBF/1.3528;
 
BF: Solenoid, L = 0.5, , KS = KSBF,    FMAPFN = "fieldmaps/BF_550.T7";

REAL I  = 273;
 REAL SF = (I/440)*1.973966;
if (OPALVERSION>15000) SF = SF/3.2306;

M:  Solenoid, L = 0.5, , KS = SF,    FMAPFN = "fieldmaps/M_440.T7";

value,{KSBF,SF};

U1: UNDULATOR, L = 10.0, LFRINGE = 1.0, Z = 0.5, K = 1.417, LAMBDA = .03,     MESHLENGTH = { 3200e-06, 3200e-06, 280e-06 }, MESHRESOLUTION = { 50e-06, 50e-06, .07e-06 },     TRUNORDER = 2, SPACECHARGE = 1,     FNAME = "undulator_output.job", TIMESTEPRATIO = 3;

GS:  Line = (GUN, BF, M);

DRIVE: Line = (GS, U1);

REAL sigmaz = 50e-05;
Dist: DISTRIBUTION, TYPE = FLATTOP,        SIGMAX = 260e-05,        SIGMAY = 260e-05,        TRISE = 6.0e-12,               TFALL = 6.0e-12,               TPULSEFWHM = sigmaz / ( beta * CLIGHT ),         CUTOFFLONG = 4.0,        NBIN = 5,        EMISSIONSTEPS = 100,        EMISSIONMODEL = ASTRA,        EKIN = 0.55,                   EMITTED = True,                WRITETOFILE = False;
    
FS_SC: Fieldsolver, FSTYPE = FFT,             MX = 16, MY = 16, MT = 16,             PARFFTX = true,             PARFFTY = true,             PARFFTT = true,             BCFFTX = open,             BCFFTY = open,             BCFFTT = open,            BBOXINCR = 1,             GREENSF = INTEGRATED;

REAL rf_freq =  beta * CLIGHT / sigmaz;
BEAM1:  BEAM, PARTICLE = ELECTRON, GAMMA = gamma, NPART = n_particles,        BFREQ = rf_freq, BCURRENT = beam_bunch_charge * rf_freq, CHARGE = -1;

TRACK, LINE = DRIVE, BEAM = BEAM1, MAXSTEPS = 1000000, DT = 1.E-11, ZSTOP = 12.0;
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,     FIELDSOLVER = FS_SC, DISTRIBUTION = Dist;
ENDTRACK;

Quit;

