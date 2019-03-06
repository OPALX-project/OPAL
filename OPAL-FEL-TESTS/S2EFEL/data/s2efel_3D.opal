OPTION, PSDUMPFREQ = 1000;
   OPTION, STATDUMPFREQ = 1000;
  OPTION, BOUNDPDESTROYFQ=10;
 OPTION, AUTOPHASE=4;
        
OPTION, VERSION=10900;

REAL REPARTFREQ=1000000;

Title, string="AWA Photoinjector";

REAL rf_freq             = 1.3e9;
     REAL n_particles         = 1E+5;
      REAL beam_bunch_charge   = 1e-9;
      
REAL Edes    = 1.4e-9;
 REAL gamma   = (Edes+EMASS)/EMASS;
 REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;
    
value , {Edes, P0};

REAL gun_inj_phase = 0.0;
  
GUN:    RFCavity, L = 0.2927, VOLT = 60.0, , TYPE = "STANDING",        FMAPFN = "fieldmaps/DriveGun.T7", FREQ = 1300.0, LAG = (gun_inj_phase*Pi)/180.0;
 
REAL KSBF = 0.162544398;
if (OPALVERSION>15000) KSBF = KSBF/1.3528;
 
BF: Solenoid, L = 0.5, , KS = KSBF,    FMAPFN = "fieldmaps/BF_550.T7";

REAL I  = 273;
 REAL SF = (I/440)*1.973966;
if (OPALVERSION>15000) SF = SF/3.2306;

M:  Solenoid, L = 0.5, , KS = SF,    FMAPFN = "fieldmaps/M_440.T7";

value,{KSBF,SF};

DR1: DRIFT, L = 10.0, Z = 0.5;

GS:  Line = (GUN, BF, M);

DRIVE: Line = (GS, DR1);

Dist: DISTRIBUTION, TYPE = FLATTOP,        SIGMAX = 0.00075,        SIGMAY = 0.00075,        TRISE = 6.0e-12,               TFALL = 6.0e-12,               TPULSEFWHM = 20.0e-12,         CUTOFFLONG = 4.0,        NBIN = 9,        EMISSIONSTEPS = 100,        EMISSIONMODEL = ASTRA,        EKIN = 0.55,                   EMITTED = True,                WRITETOFILE = False;
    
FS_SC: Fieldsolver, FSTYPE = FFT,             MX = 32, MY = 32, MT = 32,             PARFFTX = true,             PARFFTY = true,             PARFFTT = true,              BCFFTX = open,             BCFFTY = open,             BCFFTT = open,            BBOXINCR = 1,             GREENSF = INTEGRATED;

BEAM1:  BEAM, PARTICLE = ELECTRON, pc = P0, NPART = n_particles,        BFREQ = rf_freq,BCURRENT = beam_bunch_charge * rf_freq, CHARGE = -1;

TRACK, LINE = DRIVE, BEAM = BEAM1, MAXSTEPS = 1000000, DT = {1.E-12, 3.E-11}, ZSTOP={0.4, 5.05};
RUN, METHOD = "PARALLEL-T", BEAM = BEAM1,     FIELDSOLVER = FS_SC, DISTRIBUTION = Dist;
ENDTRACK;

Quit;

