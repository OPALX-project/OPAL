;;; opal-mode-el -- Major mode for editing OPAL files

;; Authors: Oscar Roberto Blanco Garcia, Christof Metzger-Kraus
;; email : <oscar.roberto.blanco.garcia@cern.ch>
;; Version: 1.0
;; Created: 17.05.2012
;; Keywords: OPAL major-mode

;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation; either version 2 of
;; the License, or (at your option) any later version.

;; This program is distributed in the hope that it will be
;; useful, but WITHOUT ANY WARRANTY; without even the implied
;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
;; PURPOSE.  See the GNU General Public License for more details.

;; You should have received a copy of the GNU General Public
;; License along with this program; if not, write to the Free
;; Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
;; MA 02111-1307 USA

;;; Commentary:
;;
;; This mode is modified from an example used in a tutorial about Emacs
;; mode creation. The tutorial can be found here:
;; https://www.emacswiki.org/emacs/ModeTutorial

;; Add this to your .emacs file to load and bind it to files with extension
;; .opal

;;; Code:

(defgroup opal nil
 "Major mode to edit OPAL files scripts in emacs"
 :group 'languages
)

(defvar opal-mode-hook nil)

;(defvar opal-mode-map
;  (let ((opal-mode-map (make-keymap)))
;    (define-key opal-mode-map "\C-j" 'newline-and-indent)
;    opal-mode-map)
;  "Keymap for OPAL major mode")

(add-to-list 'auto-mode-alist '("\\.opal\\'" . opal-mode))

; optimiser keywords
;(concat "\\<" (regexp-opt '("DVAR" "OBJECTIVE" "CONSTRAINT" "OPTIMIZE" "SAMPLING") t) "\\>")

(defconst opal-font-lock-keywords-optimise
  (list
  '("\\<\\(CONSTRAINT\\|DVAR\\|O\\(?:\\(?:BJECTIV\\|PTIMIZ\\)E\\)\\|SAMPLING\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (matchingmet).")


;(concat "\\<" (regexp-opt '("LINE" "RUN" "START" "TWISS") t) "\\>")

(defconst opal-font-lock-keywords-simul
  (list
   ; These define the beginning and end of each OPAL entity definition
  '("\\<\\(LINE\\|RUN\\|START\\|TWISS\\)\\>"
 . font-lock-builtin-face)
 )
 "Highlighting expressions for OPAL mode (simul).")

(defconst opal-font-lock-keywords-programflow
  (list
  '("\\<\\(ELSE\\(?:IF\\)?\\|IF\\|MACRO\\|WHILE\\)\\>"
  . font-lock-keyword-face)
  )
  "Highlighting expressions for OPAL mode (programflow).")

;(concat "\\<" (regexp-opt '("CALL" "CONST" "EXIT" "HELP" "OPTION" "PRINT" "QUIT" "REAL" "SELECT" "STOP" "SYSTEM" "TITLE" "VALUE") t) "\\>")

(defconst opal-font-lock-keywords-controlstm
  (list
  '("\\<\\(C\\(?:ALL\\|ONST\\)\\|EXIT\\|HELP\\|OPTION\\|PRINT\\|QUIT\\|REAL\\|S\\(?:ELECT\\|TOP\\|YSTEM\\)\\|\\(?:TITL\\|VALU\\)E\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (controlstm).")

;(concat "\\<" (regexp-opt '("CCOLLIMATOR" "CYCLOTRON" "DEGRADER" "DRIFT" "ECOLLIMATOR" "FLEXIBLECOLLIMATOR" "HKICKER" "KICKER" "MARKER" "MATRIX" "MONITOR" "MULTIPOLE" "OCTUPOLE" "PROBE" "QUADRUPOLE" "RBEND" "RCOLLIMATOR" "RFCAVITY" "RINGDEFINITION" "SBEND" "SBEND3D" "SEPTUM" "SEXTUPOLE" "SOLENOID" "SOURCE" "STRIPPER" "TRAVELINGWAVE" "TRIMCOIL" "UNDULATOR" "VARIABLE_RF_CAVITY" "VKICKER") t) "\\>")

(defconst opal-font-lock-keywords-elements
  (list
  '("\\<\\(C\\(?:COLLIMATOR\\|YCLOTRON\\)\\|D\\(?:EGRADER\\|RIFT\\)\\|ECOLLIMATOR\\|FLEXIBLECOLLIMATOR\\|HKICKER\\|KICKER\\|M\\(?:A\\(?:RKER\\|TRIX\\)\\|ONITOR\\|ULTIPOLE\\)\\|OCTUPOLE\\|PROBE\\|QUADRUPOLE\\|R\\(?:BEND\\|COLLIMATOR\\|FCAVITY\\|INGDEFINITION\\)\\|S\\(?:BEND\\(?:3D\\)?\\|E\\(?:PTUM\\|XTUPOLE\\)\\|O\\(?:LENOID\\|URCE\\)\\|TRIPPER\\)\\|TR\\(?:AVELINGWAVE\\|IMCOIL\\)\\|UNDULATOR\\|V\\(?:ARIABLE_RF_CAVITY\\|KICKER\\)\\)\\>"
  . font-lock-type-face)
  )
  "Highlighting expressions for OPAL mode (elements).")

;(concat "\\<" (regexp-opt '("BEAM" "DISTRIBUTION" "FIELDSOLVER" "POLYNOMIAL_TIME_DEPENDENCE") t) "\\>")

(defconst opal-font-lock-keywords-beamspec
  (list
  '("\\<\\(BEAM\\|DISTRIBUTION\\|FIELDSOLVER\\|POLYNOMIAL_TIME_DEPENDENCE\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (beamspec).")

;(concat "\\<" (regexp-opt '("ENDTRACK" "TRACK") t) "\\>")

(defconst opal-font-lock-keywords-orbit_corr
  (list
  '("\\<\\(ENDTRACK\\|TRACK\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (orbit_corr).")

;; attribute list can be generated in terminal with:
;ack -h '^\s*\("[A-Z0-9_]+"' src/ |sed 's/\s*("\([^"]\+\)".*/"\1"/' |sort|uniq >attributes.txt
;ack -h 'Attributes::make.*\("[A-Z0-9_]+"' |sed 's/.*Attributes::make.*("\([^"]\+\)".*$/"\1"/' |sort|uniq >> attributes.txt
;cat attributes.txt |sort|uniq >attributes2.txt

;(concat "\\<" (regexp-opt '("A" "ALL" "ALPHAX" "ALPHAY" "AMPLITUDE_MODEL" "AMR" "AMR_BFX" "AMR_BFY" "AMR_BFZ" "AMR_DENSITY" "AMR_DOMAIN_RATIO" "AMR_MAXGRIDX" "AMR_MAXGRIDY" "AMR_MAXGRIDZ" "AMR_MAXLEVEL" "AMR_MAX_NUM_PART" "AMR_MG_INTERP" "AMR_MG_NORM" "AMR_MG_NSWEEPS" "AMR_MG_PREC" "AMR_MG_REBALANCE" "AMR_MG_REUSE" "AMR_MG_SMOOTHER" "AMR_MG_TOL" "AMR_MG_VERBOSE" "AMR_MIN_NUM_PART" "AMR_REFX" "AMR_REFY" "AMR_REFZ" "AMR_REGRID_FREQ" "AMR_SCALING" "AMR_TAGGING" "AMR_YT_DUMP_FREQ" "ANGLE" "APERTURE" "APVETO" "ASCIIDUMP" "AUTOPHASE" "AZIMUTHAL_ANGLE" "AZIMUTHAL_EXTENT" "B" "B0" "BB_LENGTH" "BBLENGTH" "BBOXINCR" "BCFFTT" "BCFFTX" "BCFFTY" "BCFFTZ" "BCURRENT" "BEAM" "BEAMHALOBOUNDARY" "BEAM_PHIINIT" "BEAM_PRINIT" "BEAM_RINIT" "BETAX" "BETAY" "BFREQ" "BIRTH_CONTROL" "BMAX" "BOUNDARYGEOMETRY" "BOUNDPDESTROYFQ" "BSCALE" "C" "CATHTEMP" "CAVITY_CENTRE" "CENTRE_LENGTH" "CHARGE" "CLASS" "CLEAR" "CLOTUNEONLY" "CMD" "COEFDENOM" "COEFDENOMPHI" "COEFNUM" "COEFNUMPHI" "COLUMN" "CONDUCT" "CONST_LENGTH" "CONSTRAINTS" "CONV_HVOL_PROG" "COORDINATE_SYSTEM" "CORRT" "CORRX" "CORRY" "CORRZ" "CROSSOVER" "CSRDUMP" "CUTOFF" "CUTOFFLONG" "CUTOFFPX" "CUTOFFPY" "CUTOFFPZ" "CUTOFFR" "CUTOFFX" "CUTOFFY" "CYHARMON" "CZERO" "DDX" "DDY" "DEBIN" "DELETEONTRANSVERSEEXIT" "DELPARTFREQ" "DENERGY" "DESCRIPTION" "DESIGNENERGY" "DISTDIR" "DISTRIBUTION" "DK1" "DK1S" "DK2" "DK2S" "DK3" "DK3S" "DKN" "DKS" "DLAG" "DPHI" "DPSI" "DR" "DT" "DTAU" "DTBUNCH" "DTHETA" "DTSCINIT" "DUMP" "DUMP_DAT" "DUMP_FREQ" "DUMP_OFFSPRING" "DVARS" "DVOLT" "DX" "DY" "DZ" "E1" "E2" "EANGLE" "EBDUMP" "ECHO" "EKIN" "ELASER" "ELEMEDGE" "EMISSIONMODEL" "EMISSIONSTEPS" "EMITTED" "ENABLEHDF5" "ENABLERUTHERFORD" "END_LENGTH" "END_NORMAL_X" "END_NORMAL_Y" "END_POSITION_X" "END_POSITION_Y" "ENERGY" "EPSILON" "ESCALE" "ET" "EX" "EXPECTED_HYPERVOL" "EXPR" "EY" "FAST" "FE" "FGEOM" "FIELD_INDEX" "FIELDMAPDIR" "FIELDSOLVER" "FIELD_UNITS" "FILE" "FILE_NAME" "FILTERS" "FINT" "FLIPX" "FLIPY" "FMAPFN" "FMHIGHE" "FMLOWE" "FNAME" "FREQ" "FREQUENCY_MODEL" "FSTYPE" "FTOSCAMPLITUDE" "FTOSCPERIODS" "FULL" "GAP" "GAPWIDTH" "GAS" "GENE_MUTATION_PROBABILITY" "GEOMETRY" "GREATERTHANPI" "GREENSF" "H1" "H2" "HALOSHIFT" "HAPERT" "HARMONIC_NUMBER" "HEIGHT" "HEIGHT_NEG_EXTENT" "HEIGHT_POS_EXTENT" "HGAP" "HKICK" "HYPERVOLREFERENCE" "ID1" "ID2" "IDEALIZED" "IMAGENAME" "INFO" "INITIAL_OPTIMIZATION" "INITIALPOPULATION" "INPUT" "INPUTMOUNITS" "INSIDEPOINT" "INTENSITYCUT" "INTERPL" "IS_CLOSED" "ITSOLVER" "JSON_DUMP_FREQ" "K" "K0" "K0S" "K1" "K1S" "K2" "K2S" "K3" "K3S" "KEEP" "KICK" "KN" "KS" "L" "LAG" "LAMBDA" "LASERPROFFN" "LAT_PHIINIT" "LAT_RINIT" "LAT_THETAINIT" "LENGTH" "LENGTH_UNITS" "LFRINGE" "LINE" "LOGBENDTRAJECTORY" "LOWERBOUND" "MAGNET_END" "MAGNET_START" "MAP_ORDER" "MASS" "MATERIAL" "MAXFORDER" "MAXGENERATIONS" "MAX_HORIZONTAL_POWER" "MAXITERS" "MAX_ORDER" "MAX_R" "MAXR" "MAXSTEPS" "MAXSTEPSCO" "MAXSTEPSSI" "MAXXORDER" "MAX_Y_POWER" "MAXZ" "MB_BINNING" "MB_ETA" "MBMODE" "MEMORYDUMP" "MESHLENGTH" "MESHRESOLUTION" "MESSAGE" "METHOD" "MINBINEMITTED" "MIN_R" "MINR" "MINSTEPFORREBIN" "MINZ" "MODE" "MT" "MTSSUBSTEPS" "MUTATION" "MUTATION_PROBABILITY" "MX" "MY" "MZ" "N" "NAME" "NBIN" "NFREQ" "NHOLX" "NHOLY" "NLEFT" "NLHS" "NPART" "NPEAKS" "NPOINTS" "NRIGHT" "NSECTORS" "NSLICES" "NSTEPS" "NUMBLOCKS" "NUMCELLS" "NUM_COWORKERS" "NUM_IND_GEN" "NUM_MASTERS" "NUMPERIODS" "OBJECTIVES" "OFFSETP" "OFFSETPX" "OFFSETPY" "OFFSETPZ" "OFFSETT" "OFFSETX" "OFFSETY" "OFFSETZ" "ONE_PILOT_CONVERGE" "OPCHARGE" "OPMASS" "OPYIELD" "ORDER" "ORDERMAPS" "ORIENTATION" "ORIGIN" "OUTDIR" "OUTFN" "OUTPUT" "P1" "P2" "P3" "PARAMB" "PARFFTT" "PARFFTX" "PARFFTY" "PARTICLE" "PARTICLEMATTERINTERACTION" "PATTERN" "PC" "PDIS" "PHASE_MODEL" "PHI" "PHI0" "PHIINIT" "PHIMAX" "PHIMIN" "PHI_START" "PHI_STEPS" "PMAPFN" "POLYORDER" "PRECMODE" "PRESSURE" "PRINIT" "PSCALE" "PSDUMPEACHTURN" "PSDUMPFRAME" "PSDUMPFREQ" "PSI" "PT" "PXMULT" "PYMULT" "PZINIT" "PZMULT" "R" "R0" "R51" "R52" "R61" "R62" "RADIAL_NEG_EXTENT" "RADIAL_POS_EXTENT" "RADIUS" "RANDOM" "RANGE" "RASTER" "RC" "REBINFREQ" "RECOMBINATION_PROBABILITY" "RECYCLEBLOCKS" "REFER" "REFPOS" "REMOTEPARTDEL" "REPARTFREQ" "RESIDUUM" "RESTART_FILE" "RESTART_STEP" "RFFCFN" "RFFREQ" "RFMAPFN" "RFPHI" "RFRINGE" "RFVCFN" "RGUESS" "RHODUMP" "RINIT" "RMAX" "RMIN" "RNGTYPE" "ROTATE180" "ROTATE270" "ROTATE90" "ROTATION" "R_START" "R_STEPS" "S" "SAMPLINGS" "SBIN" "SCALABLE" "SCALE" "SCSOLVEFREQ" "SECTOR" "SEED" "SEPPEAKS" "SIGMA" "SIGMAPT" "SIGMAPX" "SIGMAPY" "SIGMAPZ" "SIGMAR" "SIGMAT" "SIGMAX" "SIGMAY" "SIGMAZ" "SIMBIN_CROSSOVER_NU" "SIMTMPDIR" "SLICES" "SLPTC" "SOL_SYNCH" "SPIRAL" "SPTDUMPFREQ" "STARTPOPULATION" "STATDUMPFREQ" "STEP" "STEPSIZE" "STEPSPERTURN" "STOP" "STOREOBJECTIVES" "STRING" "SUPERPOSE" "SURFDUMPFREQ" "SYMMETRY" "T" "T0" "TABLE" "TAN_DELTA" "TANGENTIAL_OFFSET" "TAU" "TELL" "TEMPERATURE" "TEMPLATEDIR" "TFALL" "THETA" "THETA_IN" "THETA_OUT" "THRESHOLD" "TIMEINTEGRATOR" "TIMES" "TMULT" "TOL" "TOPO" "TOTALTIME" "TP" "TPULSEFWHM" "TRACE" "TRACKBACK" "TRANSPARENT" "TRIMCOIL" "TRIMCOILTHRESHOLD" "TRISE" "TRUNORDER" "T_START" "T_STEPS" "TURNS" "TYPE" "UPPERBOUND" "VALUE" "VALUES" "VAPERT" "VARIABLE" "VARRADIUS" "VERSION" "VKICK" "VOLT" "W" "WAKEF" "WARN" "WEIGHT" "WIDTH" "WRITETOFILE" "X" "XEND" "XMULT" "XSCALE" "XSIZE" "X_START" "XSTART" "X_STEPS" "XYZSCALE" "Y" "YEND" "YMULT" "YSCALE" "YSIZE" "Y_START" "YSTART" "Y_STEPS" "Z" "Z0" "ZEND" "ZINIT" "ZSCALE" "ZSHIFT" "Z_START" "ZSTART" "Z_STEPS" "ZSTOP") t) "\\>")

(defconst opal-font-lock-keywords-parameters
  (list
 '("\\<\\(A\\(?:L\\(?:L\\|PHA[XY]\\)\\|M\\(?:PLITUDE_MODEL\\|R\\(?:_\\(?:BF[XYZ]\\|D\\(?:ENSITY\\|OMAIN_RATIO\\)\\|M\\(?:AX\\(?:GRID[XYZ]\\|LEVEL\\|_NUM_PART\\)\\|G_\\(?:INTERP\\|N\\(?:ORM\\|SWEEPS\\)\\|PREC\\|RE\\(?:\\(?:BALANC\\|US\\)E\\)\\|SMOOTHER\\|TOL\\|VERBOSE\\)\\|IN_NUM_PART\\)\\|RE\\(?:F[XYZ]\\|GRID_FREQ\\)\\|SCALING\\|TAGGING\\|YT_DUMP_FREQ\\)\\)?\\)\\|NGLE\\|P\\(?:ERTURE\\|VETO\\)\\|SCIIDUMP\\|UTOPHASE\\|ZIMUTHAL_\\(?:ANGLE\\|EXTENT\\)\\)\\|B\\(?:0\\|B\\(?:LENGTH\\|OXINCR\\|_LENGTH\\)\\|C\\(?:FFT[TXYZ]\\|URRENT\\)\\|E\\(?:AM\\(?:HALOBOUNDARY\\|_\\(?:\\(?:P\\(?:HI\\|R\\)\\|R\\)INIT\\)\\)?\\|TA[XY]\\)\\|FREQ\\|IRTH_CONTROL\\|MAX\\|OUND\\(?:ARYGEOMETRY\\|PDESTROYFQ\\)\\|SCALE\\)\\|C\\(?:A\\(?:THTEMP\\|VITY_CENTRE\\)\\|ENTRE_LENGTH\\|HARGE\\|L\\(?:ASS\\|EAR\\|OTUNEONLY\\)\\|MD\\|O\\(?:EF\\(?:DENOM\\(?:PHI\\)?\\|NUM\\(?:PHI\\)?\\)\\|LUMN\\|N\\(?:DUCT\\|ST\\(?:RAINTS\\|_LENGTH\\)\\|V_HVOL_PROG\\)\\|ORDINATE_SYSTEM\\|RR[TXYZ]\\)\\|ROSSOVER\\|SRDUMP\\|UTOFF\\(?:LONG\\|P[XYZ]\\|[RXY]\\)?\\|YHARMON\\|ZERO\\)\\|D\\(?:D[XY]\\|E\\(?:BIN\\|L\\(?:ETEONTRANSVERSEEXIT\\|PARTFREQ\\)\\|NERGY\\|S\\(?:CRIPTION\\|IGNENERGY\\)\\)\\|IST\\(?:DIR\\|RIBUTION\\)\\|K\\(?:[123]S\\|[123NS]\\)\\|LAG\\|P\\(?:[HS]I\\)\\|T\\(?:AU\\|BUNCH\\|HETA\\|SCINIT\\)\\|UMP\\(?:_\\(?:DAT\\|FREQ\\|OFFSPRING\\)\\)?\\|V\\(?:ARS\\|OLT\\)\\|[RTXYZ]\\)\\|E\\(?:ANGLE\\|BDUMP\\|CHO\\|KIN\\|L\\(?:ASER\\|EMEDGE\\)\\|MI\\(?:SSION\\(?:MODEL\\|STEPS\\)\\|TTED\\)\\|N\\(?:ABLE\\(?:HDF5\\|RUTHERFORD\\)\\|D_\\(?:LENGTH\\|NORMAL_[XY]\\|POSITION_[XY]\\)\\|ERGY\\)\\|PSILON\\|SCALE\\|XP\\(?:ECTED_HYPERVOL\\|R\\)\\|[12TXY]\\)\\|F\\(?:AST\\|E\\|GEOM\\|I\\(?:ELD\\(?:MAPDIR\\|SOLVER\\|_\\(?:INDEX\\|UNITS\\)\\)\\|L\\(?:E\\(?:_NAME\\)?\\|TERS\\)\\|NT\\)\\|LIP[XY]\\|M\\(?:APFN\\|\\(?:HIGH\\|LOW\\)E\\)\\|NAME\\|REQ\\(?:UENCY_MODEL\\)?\\|STYPE\\|TOSC\\(?:AMPLITUDE\\|PERIODS\\)\\|ULL\\)\\|G\\(?:A\\(?:PWIDTH\\|[PS]\\)\\|E\\(?:\\(?:NE_MUTATION_PROBABILIT\\|OMETR\\)Y\\)\\|RE\\(?:ATERTHANPI\\|ENSF\\)\\)\\|H\\(?:A\\(?:LOSHIFT\\|PERT\\|RMONIC_NUMBER\\)\\|EIGHT\\(?:_\\(?:\\(?:NEG\\|POS\\)_EXTENT\\)\\)?\\|GAP\\|KICK\\|YPERVOLREFERENCE\\|[12]\\)\\|I\\(?:D\\(?:EALIZED\\|[12]\\)\\|MAGENAME\\|N\\(?:FO\\|ITIAL\\(?:\\(?:POPUL\\|_OPTIMIZ\\)ATION\\)\\|PUT\\(?:MOUNITS\\)?\\|SIDEPOINT\\|TE\\(?:NSITYCUT\\|RPL\\)\\)\\|S_CLOSED\\|TSOLVER\\)\\|JSON_DUMP_FREQ\\|K\\(?:0S\\|1S\\|2S\\|3S\\|EEP\\|ICK\\|[0-3NS]\\)\\|L\\(?:A\\(?:G\\|MBDA\\|SERPROFFN\\|T_\\(?:\\(?:PHI\\|R\\|THETA\\)INIT\\)\\)\\|ENGTH\\(?:_UNITS\\)?\\|FRINGE\\|INE\\|O\\(?:GBENDTRAJECTORY\\|WERBOUND\\)\\)\\|M\\(?:A\\(?:GNET_\\(?:END\\|START\\)\\|P_ORDER\\|SS\\|TERIAL\\|X\\(?:FORDER\\|GENERATIONS\\|ITERS\\|STEPS\\(?:CO\\|SI\\)?\\|\\(?:XORDE\\|_\\(?:\\(?:HORIZONTAL_POW\\|ORD\\|Y_POW\\)E\\)?\\)R\\|[RZ]\\)\\)\\|B\\(?:MODE\\|_\\(?:BINNING\\|ETA\\)\\)\\|E\\(?:MORYDUMP\\|S\\(?:H\\(?:LENGTH\\|RESOLUTION\\)\\|SAGE\\)\\|THOD\\)\\|IN\\(?:BINEMITTED\\|STEPFORREBIN\\|_R\\|[RZ]\\)\\|ODE\\|TSSUBSTEPS\\|UTATION\\(?:_PROBABILITY\\)?\\|[TXYZ]\\)\\|N\\(?:AME\\|BIN\\|FREQ\\|HOL[XY]\\|L\\(?:EFT\\|HS\\)\\|P\\(?:ART\\|\\(?:EAK\\|OINT\\)S\\)\\|RIGHT\\|S\\(?:\\(?:ECTOR\\|LICE\\|TEP\\)S\\)\\|UM\\(?:BLOCKS\\|CELLS\\|PERIODS\\|_\\(?:COWORKERS\\|IND_GEN\\|MASTERS\\)\\)\\)\\|O\\(?:BJECTIVES\\|FFSET\\(?:P[XYZ]\\|[PTXYZ]\\)\\|NE_PILOT_CONVERGE\\|P\\(?:CHARGE\\|MASS\\|YIELD\\)\\|R\\(?:DER\\(?:MAPS\\)?\\|I\\(?:\\(?:ENTATIO\\|GI\\)N\\)\\)\\|UT\\(?:DIR\\|FN\\|PUT\\)\\)\\|P\\(?:A\\(?:R\\(?:AMB\\|FFT[TXY]\\|TICLE\\(?:MATTERINTERACTION\\)?\\)\\|TTERN\\)\\|DIS\\|H\\(?:ASE_MODEL\\|I\\(?:0\\|INIT\\|M\\(?:AX\\|IN\\)\\|_ST\\(?:ART\\|EPS\\)\\)?\\)\\|MAPFN\\|OLYORDER\\|R\\(?:E\\(?:\\(?:CMOD\\|SSUR\\)E\\)\\|INIT\\)\\|S\\(?:CALE\\|DUMP\\(?:EACHTURN\\|FR\\(?:AME\\|EQ\\)\\)\\|I\\)\\|\\(?:XMUL\\|YMUL\\|Z\\(?:INI\\|MUL\\)\\)T\\|[123CT]\\)\\|R\\(?:5[12]\\|6[12]\\|A\\(?:DI\\(?:AL_\\(?:\\(?:NEG\\|POS\\)_EXTENT\\)\\|US\\)\\|N\\(?:DOM\\|GE\\)\\|STER\\)\\|E\\(?:BINFREQ\\|C\\(?:OMBINATION_PROBABILITY\\|YCLEBLOCKS\\)\\|F\\(?:ER\\|POS\\)\\|MOTEPARTDEL\\|PARTFREQ\\|S\\(?:IDUUM\\|TART_\\(?:FILE\\|STEP\\)\\)\\)\\|F\\(?:F\\(?:CFN\\|REQ\\)\\|MAPFN\\|PHI\\|RINGE\\|VCFN\\)\\|GUESS\\|HODUMP\\|INIT\\|M\\(?:AX\\|IN\\)\\|NGTYPE\\|OTAT\\(?:E\\(?:\\(?:18\\|27\\|9\\)0\\)\\|ION\\)\\|_ST\\(?:ART\\|EPS\\)\\|[0C]\\)\\|S\\(?:AMPLINGS\\|BIN\\|C\\(?:AL\\(?:\\(?:ABL\\)?E\\)\\|SOLVEFREQ\\)\\|E\\(?:CTOR\\|ED\\|PPEAKS\\)\\|I\\(?:GMA\\(?:P[TXYZ]\\|[RTXYZ]\\)?\\|M\\(?:BIN_CROSSOVER_NU\\|TMPDIR\\)\\)\\|L\\(?:ICES\\|PTC\\)\\|OL_SYNCH\\|P\\(?:IRAL\\|TDUMPFREQ\\)\\|T\\(?:A\\(?:RTPOPULATION\\|TDUMPFREQ\\)\\|EP\\(?:S\\(?:IZE\\|PERTURN\\)\\)?\\|O\\(?:P\\|REOBJECTIVES\\)\\|RING\\)\\|U\\(?:PERPOSE\\|RFDUMPFREQ\\)\\|YMMETRY\\)\\|T\\(?:A\\(?:BLE\\|N\\(?:GENTIAL_OFFSET\\|_DELTA\\)\\|U\\)\\|E\\(?:LL\\|MP\\(?:ERATURE\\|LATEDIR\\)\\)\\|FALL\\|H\\(?:ETA\\(?:_\\(?:IN\\|OUT\\)\\)?\\|RESHOLD\\)\\|IME\\(?:INTEGRATOR\\|S\\)\\|MULT\\|O\\(?:L\\|PO\\|TALTIME\\)\\|PULSEFWHM\\|R\\(?:A\\(?:C\\(?:E\\|KBACK\\)\\|NSPARENT\\)\\|I\\(?:MCOIL\\(?:THRESHOLD\\)?\\|SE\\)\\|UNORDER\\)\\|URNS\\|YPE\\|_ST\\(?:ART\\|EPS\\)\\|[0P]\\)\\|UPPERBOUND\\|V\\(?:A\\(?:LUES?\\|PERT\\|R\\(?:IABLE\\|RADIUS\\)\\)\\|ERSION\\|KICK\\|OLT\\)\\|W\\(?:A\\(?:KEF\\|RN\\)\\|EIGHT\\|IDTH\\|RITETOFILE\\)\\|X\\(?:END\\|MULT\\|S\\(?:CALE\\|IZE\\|TART\\)\\|YZSCALE\\|_ST\\(?:ART\\|EPS\\)\\)\\|Y\\(?:END\\|MULT\\|S\\(?:CALE\\|IZE\\|TART\\)\\|_ST\\(?:ART\\|EPS\\)\\)\\|Z\\(?:0\\|END\\|INIT\\|S\\(?:CALE\\|HIFT\\|T\\(?:ART\\|OP\\)\\)\\|_ST\\(?:ART\\|EPS\\)\\)\\|[ABCKLNRSTW-Z]\\)\\>"
  . font-lock-variable-name-face)
  )
  "Highlighting expressions for OPAL mode (parameters).")

;(concat "\\<" (regexp-opt '("EALIGN" "EFCOMP" "ERROR") t) "\\>")

(defconst opal-font-lock-keywords-errordef
  (list
  '("\\<\\(E\\(?:ALIGN\\|FCOMP\\|RROR\\)\\)\\>"
  . font-lock-warning-face)
  )
  "Highlighting expressions for OPAL mode (errordef).")

;(concat "\\<" (regexp-opt '("AC" "AIR" "ALPHA" "ALUMINAAL2O3" "ALUMINUM" "ANTIPROTON" "ASTRA" "ASTRAFLATTOPTH" "AUTO" "AVFEQ" "BANDRF" "BEAMSTRIPPING" "BERYLLIUM" "BICGSTAB" "BINOMIAL" "BLEND" "BLOCK_GS" "BLOCK_JACOBI" "BORONCARBIDE" "BOXCORNER" "BUNCH_BINNING" "BUNCH_MEAN" "CARBON" "CARBONCYCL" "CARTESIAN" "CENTRE" "CG" "CHEBYSHEV" "CLIGHT" "CMASS" "COLLIM" "CONSTANT" "COPPER" "CYCIAE" "CYLINDRICAL" "DC" "DEGRAD" "DEUTERON" "DIRICHLET" "DMASS" "E" "ELECTRON" "ELLIPTIC" "EMASS" "ENTRY" "EVOVERC" "EXIT" "FALSE" "FFA" "FFT" "FFTPERIODIC" "FIXEDFFTLOWPASS" "FLATTOP" "FORCE" "FROMFILE" "GAMMA_BINNING" "GAUSS" "GAUSSIAN" "GAUSSMATCHED" "GITREVISION" "GLOBAL" "GMRES" "GOLD" "GRAPHITE" "GRAPHITER6710" "GS" "GUNGAUSSFLATTOPTH" "H2" "H2P" "HALTON" "HIERARCHY" "HMINUS" "HMMASS" "ILUT" "INDEPENDENTBIT" "INTEGRATED" "JACOBI" "KAPTON" "L1_NORM" "L2_NORM" "LAGRANGE" "LATIN_HYPERCUBE" "LF2" "LINEAR" "LINF_NORM" "MMASS" "MOLYBDENUM" "MTS" "MULTIGAUSS" "MUON" "MYLAR" "NAIVEONEPOINT" "NAIVEUNIFORM" "NIEDERREITER" "NONE" "NONEQUIL" "ONEBIT" "OPEN" "P0" "P3M" "PC" "PERIODIC" "PI" "PMASS" "POSITRON" "PROTON" "QUADRATIC" "RADDEG" "RANDOM_SEQUENCE_UNIFORM" "RANDOM_SEQUENCE_UNIFORM_INT" "RAP" "RECTANGULAR" "REFERENCE" "RELATIVEFFTLOWPASS" "REUSE" "RILUK" "RING" "RK4" "RP" "SA" "SAAMG" "SCATTERING" "SGS" "SIMULATEDBINARY" "SINGLEGAP" "SOBOL" "SPATIAL" "STANDARD" "STANDING" "STD" "STENCIL" "SYNCHROCYCLOTRON" "TEMPORAL" "THICK" "TITANIUM" "TRILINEAR" "TRUE" "TWOPI" "UMASS" "UNIFORM" "UNIFORM_INT" "URANIUM" "WATER" "XEMASS" "XENON") t) "\\>")

(defconst opal-font-lock-keywords-constants
  (list
  '("\\<\\(A\\(?:C\\|IR\\|L\\(?:PHA\\|UMIN\\(?:AAL2O3\\|UM\\)\\)\\|NTIPROTON\\|STRA\\(?:FLATTOPTH\\)?\\|UTO\\|VFEQ\\)\\|B\\(?:ANDRF\\|E\\(?:AMSTRIPPING\\|RYLLIUM\\)\\|I\\(?:CGSTAB\\|NOMIAL\\)\\|L\\(?:END\\|OCK_\\(?:GS\\|JACOBI\\)\\)\\|O\\(?:RONCARBIDE\\|XCORNER\\)\\|UNCH_\\(?:BINNING\\|MEAN\\)\\)\\|C\\(?:AR\\(?:BON\\(?:CYCL\\)?\\|TESIAN\\)\\|ENTRE\\|G\\|HEBYSHEV\\|LIGHT\\|MASS\\|O\\(?:LLIM\\|NSTANT\\|PPER\\)\\|Y\\(?:CIAE\\|LINDRICAL\\)\\)\\|D\\(?:C\\|E\\(?:GRAD\\|UTERON\\)\\|IRICHLET\\|MASS\\)\\|E\\(?:L\\(?:ECTRON\\|LIPTIC\\)\\|MASS\\|NTRY\\|VOVERC\\|XIT\\)?\\|F\\(?:ALSE\\|F\\(?:TPERIODIC\\|[AT]\\)\\|IXEDFFTLOWPASS\\|LATTOP\\|\\(?:ORC\\|ROMFIL\\)E\\)\\|G\\(?:A\\(?:MMA_BINNING\\|USS\\(?:IAN\\|MATCHED\\)?\\)\\|ITREVISION\\|LOBAL\\|MRES\\|OLD\\|RAPHITE\\(?:R6710\\)?\\|S\\|UNGAUSSFLATTOPTH\\)\\|H\\(?:2P?\\|ALTON\\|IERARCHY\\|M\\(?:\\(?:INU\\|MAS\\)S\\)\\)\\|I\\(?:LUT\\|N\\(?:DEPENDENTBIT\\|TEGRATED\\)\\)\\|JACOBI\\|KAPTON\\|L\\(?:1_NORM\\|2_NORM\\|A\\(?:\\(?:GRANG\\|TIN_HYPERCUB\\)E\\)\\|F_2\\|IN\\(?:EAR\\|F_NORM\\)\\)\\|M\\(?:MASS\\|OLYBDENUM\\|TS\\|U\\(?:LTIGAUSS\\|ON\\)\\|YLAR\\)\\|N\\(?:AIVE\\(?:ONEPOINT\\|UNIFORM\\)\\|IEDERREITER\\|ONE\\(?:QUIL\\)?\\)\\|O\\(?:NEBIT\\|PEN\\)\\|P\\(?:3M\\|ERIODIC\\|MASS\\|\\(?:OSITR\\|ROT\\)ON\\|[0CI]\\)\\|QUADRATIC\\|R\\(?:A\\(?:DDEG\\|P\\)\\|E\\(?:CTANGULAR\\|FERENCE\\|LATIVEFFTLOWPASS\\|USE\\)\\|I\\(?:LUK\\|NG\\)\\|K_4\\|P\\)\\|S\\(?:A\\(?:AMG\\)?\\|CATTERING\\|GS\\|I\\(?:MULATEDBINARY\\|NGLEGAP\\)\\|OBOL\\|PATIAL\\|T\\(?:AND\\(?:ARD\\|ING\\)\\|D\\|ENCIL\\)\\|YNCHROCYCLOTRON\\)\\|T\\(?:EMPORAL\\|HICK\\|ITANIUM\\|R\\(?:ILINEAR\\|UE\\)\\|WOPI\\)\\|U\\(?:MASS\\|NIFORM\\(?:_INT\\)?\\|RANIUM\\)\\|WATER\\|XE\\(?:MASS\\|NON\\)\\)\\>"
  . font-lock-constant-face)
  )
  "Highlighting expressions for OPAL mode (constants).")

;(regexp-opt '("TITLE") t)

(defconst opal-font-lock-keywords-stringatt
  (list
  '("\\<\\(TITLE\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (stringatt).")

;(concat "\\<" (regexp-opt '("ABS" "ACOS" "ASIN" "ATAN" "ATAN2" "COS" "COSH" "EXP" "GAUSS" "LOG" "LOG10" "RANF" "SIN" "SINH" "SQRT" "TAN" "TANH" "TGAUSS") t) "\\>")

(defconst opal-font-lock-keywords-functions
  (list
  '("\\<\\(A\\(?:BS\\|COS\\|SIN\\|TAN2?\\)\\|COSH?\\|EXP\\|GAUSS\\|LOG\\(?:10\\)?\\|RANF\\|S\\(?:INH?\\|QRT\\)\\|T\\(?:ANH?\\|GAUSS\\)\\)\\>"
  . font-lock-function-name-face)
  )
  "Highlighting expressions for OPAL mode (functions).")

(defconst opal-font-lock-special_operators
  (list
   '("\\(->\\|:=\\)"
  . font-lock-warning-face)
  )
  "Highlighting expressions for OPAL mode (variables_opal).")

(defconst opal-font-lock-special_constants
  (list
   '("\\(#[es]\\)"
  . font-lock-constant-face)
  )
  "Highlighting expressions for OPAL mode (variables_opal).")


(defconst opal-font-lock-keywords-3
  (append
     opal-font-lock-keywords-optimise
     opal-font-lock-keywords-programflow
     opal-font-lock-keywords-simul
     opal-font-lock-keywords-controlstm
     opal-font-lock-keywords-elements
     opal-font-lock-keywords-beamspec
     opal-font-lock-keywords-orbit_corr
     opal-font-lock-keywords-parameters
     opal-font-lock-keywords-errordef
     opal-font-lock-keywords-constants
     opal-font-lock-keywords-stringatt
     opal-font-lock-keywords-functions
     opal-font-lock-special_operators
     opal-font-lock-special_constants
  )
 "Balls-out highlighting in OPAL mode.")

(defvar opal-font-lock-keywords opal-font-lock-keywords-3
  "Default highlighting expressions for OPAL mode.")

(defvar opal-mode-syntax-table
  (let ((opal-mode-syntax-table (make-syntax-table c-mode-syntax-table)))

    ; This is added so entity names with underscores can be more easily parsed
        (modify-syntax-entry ?_ "w" opal-mode-syntax-table)
        (modify-syntax-entry ?. "w" opal-mode-syntax-table)

        ;  Comment styles are same as C++
        (modify-syntax-entry ?/ ". 124 b" opal-mode-syntax-table)
        (modify-syntax-entry ?* ". 23" opal-mode-syntax-table)
        (modify-syntax-entry ?\n "> b" opal-mode-syntax-table)
        (modify-syntax-entry ?! "< b" opal-mode-syntax-table)
        (modify-syntax-entry ?' "|" opal-mode-syntax-table)
        opal-mode-syntax-table)
  "Syntax table for opal-mode")

;;; ### autoload
(defun opal-mode ()
  "Major mode for editing OPAL script files"
  (interactive)
  (kill-all-local-variables)
  (setq mode-name "OPAL")
  (setq major-mode 'opal-mode)
  (setq comment-start "//")
;  (use-local-map opal-mode-map)
  (set-syntax-table opal-mode-syntax-table)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(opal-font-lock-keywords nil t))
;; Set up search
  (add-hook 'opal-mode-hook
     (lambda ()  (setq case-fold-search t)))
  (run-hooks 'opal-mode-hook)
)
(provide 'opal-mode)

;;; opal-mode.el ends here
