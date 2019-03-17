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


;(concat "\\<" (regexp-opt '("CONSTRAINT" "ENDMATCH" "LINE" "MATCH" "RUN" "START" "TWISS" "VARY") t) "\\>")

(defconst opal-font-lock-keywords-simul
  (list
   ; These define the beginning and end of each OPAL entity definition
  '("\\<\\(ENDMATCH\\|LINE\\|MATCH\\|RUN\\|START\\|TWISS\\|VARY\\)\\>"
 . font-lock-builtin-face)
 )
 "Highlighting expressions for OPAL mode (simul).")

(defconst opal-font-lock-keywords-programflow
  (list
  '("\\<\\(ELSE\\(?:IF\\)?\\|IF\\|MACRO\\|WHILE\\)\\>"
  . font-lock-keyword-face)
  )
  "Highlighting expressions for OPAL mode (programflow).")

;(concat "\\<" (regexp-opt '("CALL" "CONST" "ESAVE" "EXIT" "HELP" "OPTION" "PRINT" "QUIT" "REAL" "SAVE" "SELECT" "SHOW" "STOP" "SYSTEM" "TITLE" "VALUE") t) "\\>")

(defconst opal-font-lock-keywords-controlstm
  (list
  '("\\<\\(C\\(?:ALL\\|ONST\\)\\|E\\(?:SAVE\\|XIT\\)\\|HELP\\|OPTION\\|PRINT\\|QUIT\\|REAL\\|S\\(?:AVE\\|ELECT\\|HOW\\|TOP\\|YSTEM\\)\\|\\(?:TITL\\|VALU\\)E\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (controlstm).")

;(concat "\\<" (regexp-opt '("CCOLLIMATOR" "CYCLOTRON" "DEGRADER" "DRIFT" "ECOLLIMATOR" "FLEXIBLECOLLIMATOR" "HKICKER" "KICKER" "MARKER" "MATRIX" "MONITOR" "MULTIPOLE" "OCTUPOLE" "PROBE" "QUADRUPOLE" "RBEND" "RCOLLIMATOR" "RFCAVITY" "RINGDEFINITION" "SBEND" "SBEND3D" "SEPTUM" "SEXTUPOLE" "SOLENOID" "SOURCE" "STRIPPER" "TRAVELINGWAVE" "TRIMCOIL" "VARIABLE_RF_CAVITY" "VKICKER") t) "\\>")

(defconst opal-font-lock-keywords-elements
  (list
  '("\\<\\(C\\(?:COLLIMATOR\\|YCLOTRON\\)\\|D\\(?:EGRADER\\|RIFT\\)\\|ECOLLIMATOR\\|FLEXIBLECOLLIMATOR\\|HKICKER\\|KICKER\\|M\\(?:A\\(?:RKER\\|TRIX\\)\\|ONITOR\\|ULTIPOLE\\)\\|OCTUPOLE\\|PROBE\\|QUADRUPOLE\\|R\\(?:BEND\\|COLLIMATOR\\|FCAVITY\\|INGDEFINITION\\)\\|S\\(?:BEND\\(?:3D\\)?\\|E\\(?:PTUM\\|XTUPOLE\\)\\|O\\(?:LENOID\\|URCE\\)\\|TRIPPER\\)\\|TR\\(?:AVELINGWAVE\\|IMCOIL\\)\\|V\\(?:ARIABLE_RF_CAVITY\\|KICKER\\)\\)\\>"
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

;(concat "\\<" (regexp-opt '("LMDIF" "MIGRAD" "SIMPLEX") t) "\\>")

(defconst opal-font-lock-keywords-matchingmet
  (list
  '("\\<\\(LMDIF\\|MIGRAD\\|SIMPLEX\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (matchingmet).")

;(concat "\\<" (regexp-opt '("ENDTRACK" "SURVEY" "TRACK") t) "\\>")

(defconst opal-font-lock-keywords-orbit_corr
  (list
  '("\\<\\(ENDTRACK\\|SURVEY\\|TRACK\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (orbit_corr).")

;(concat "\\<" (regexp-opt '("EPRINT") t) "\\>")

(defconst opal-font-lock-keywords-plot
  (list
  '("\\<\\(EPRINT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (plot).")

;(concat "\\<" (regexp-opt '("CYCLE" "ENDEDIT" "FLATTEN" "INSTALL" "MOVE" "REFLECT" "REMOVE" "SEQEDIT") t) "\\>")

(defconst opal-font-lock-keywords-seqediting
  (list
  '("\\<\\(CYCLE\\|ENDEDIT\\|FLATTEN\\|INSTALL\\|MOVE\\|RE\\(?:FLECT\\|MOVE\\)\\|SEQEDIT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (seqediting).")

;(concat "\\<" (regexp-opt '("A" "ADD" "AMPLITUDE_MODEL" "ANGLE" "APERTURE" "APVETO" "AT" "AUTOPHASE" "B" "BBOXINCR" "BCFFTZ" "BCFFTT" "BCFFTX" "BCFFTY" "BCURRENT" "BEAM_PHIINIT" "BEAM_PRINIT" "BEAM_RINIT" "BFREQ" "BMAX" "BOUNDPDESTROYFQ" "BSCALE" "BY" "CALLS" "CATHTEMP" "CENTRE" "CHARGE" "CLASS" "CLEAR" "CMD" "COEFDENOM" "COEFDENOMPHI" "COEFNUM" "COEFNUMPHI" "COLUMN" "CONDUCT" "CONSTRAINTS" "CONST_LENGTH" "CONV_HVOL_PROG" "CORRX" "CORRY" "CORRZ" "CSRDUMP" "CROSSOVER" "CUTOFFLONG" "CUTOFFPX" "CUTOFFPY" "CUTOFFPZ" "CUTOFFR" "CUTOFFX" "CUTOFFY" "CYHARMON" "CZERO" "DESIGNENERGY" "DK1" "DK1S" "DK2" "DK2S" "DKN" "DKNR" "DKS" "DKSR" "DLAG" "DPHI" "DPSI" "DS" "DT" "DTHETA" "DUMP" "DUMP_DAT" "DUMP_FREQ" "DUMP_OFFSPRING" "DVARS" "DVOLT" "DX" "DY" "DZ" "E1" "E2" "EBDUMP" "ECHO" "EKIN" "ELASER" "ELEMEDGE" "ELEMENT" "EMISSIONMODEL" "EMISSIONSTEPS" "EMITTED" "ENABLEHDF5" "ENERGY" "ENDSEQUENCE" "END_NORMAL_X" "END_NORMAL_Y" "END_POSITION_X" "END_POSITION_Y" "EPSILON" "ESCALE" "ET" "EVERYSTEP" "EX" "EXPECTED_HYPERVOL" "EXPR" "EY" "FE" "FGEOM" "FIELDMAPDIR" "FIELDMAPDIR" "FILE" "FINT" "FMAPFN" "FNAME" "FORM" "FREQ" "FREQUENCY_MODEL" "FROM" "FSTYPE" "FTOSCAMPLITUDE" "FTOSCPERIODS" "FULL" "GAMMA" "GAP" "GAPWIDTH" "GENE_MUTATION_PROBABILITY" "GEOMETRY" "GREENSF" "H1" "H2" "HAPERT" "HARMON" "HARMONIC_NUMBER" "HGAP" "HKICK" "HYPERVOLREFERENCE" "IDEALIZED" "INITIAL_OPTIMIZATION" "IMAGENAME" "INFO" "INITIALPOPULATION" "INPUT" "INPUTMOUNITS" "INTENSITYCUT" "INTERPL" "IS_CLOSED" "ITSOLVER" "K0" "K0S" "K1" "K1S" "K2" "K2S" "K3" "K3S" "KEYWORD" "KICK" "KN" "KS" "L" "LAG" "LASERPROFFN" "LATTICE_PHIINIT" "LATTICE_RINIT" "LATTICE_THETAINIT" "LENGTH" "LEVEL" "LOGBENDTRAJECTORY" "LOWER" "LOWERBOUND" "MASS" "MAXGENERATIONS" "MAXITERS" "MAXR" "MAXSTEPS" "MAXZ" "METHOD" "MINR" "MINZ" "MODE" "MREX" "MREY" "MSCALX" "MSCALY" "MT" "MUTATION_PROBABILITY" "MX" "MY" "NBIN" "NFREQ" "NLEFT" "NO" "NPART" "NPOINTS" "NRIGHT" "NUMCELLS" "NUM_COWORKERS" "NUM_IND_GEN" "NUM_MASTERS" "OBJECTIVES" "OFFSETPX" "OFFSETPY" "OFFSETPZ" "OFFSETT" "OFFSETX" "OFFSETY" "OFFSETZ" "ONE_PILOT_CONVERGE" "OPCHARGE" "OPMASS" "OPYIELD" "ORDER" "ORIENTATION" "ORIGIN" "OUTDIR" "OUTFN" "OUTPUT" "P1" "P2" "P3" "P4" "PARFFTT" "PARFFTX" "PARFFTY" "PARTICLE" "PARTICLEMATTERINTERACTION" "PATTERN" "PC" "PDIS" "PHI" "PHI0" "PHIINIT" "PHIMAX" "PHIMIN" "PLANE" "POLYORDER" "PRECMODE" "PRINIT" "PSDUMPEACHTURN" "PSDUMPFRAME" "PSDUMPFREQ" "PSI" "PTC" "PYMULT" "PZINIT" "PZMULT" "R51" "R52" "R61" "R62" "RADIUS" "RANDOM" "RANGE" "RASTER" "RECOMBINATION_PROBABILITY" "REFER" "REFPOS" "REPARTFREQ" "RESET" "RFFREQ" "RFMAPFN" "RFPHI" "RINIT" "RMAX" "RMIN" "ROTATION" "ROW" "S" "SAMPLINGS" "SCALABLE" "SEED" "SELECTED" "SEQUENCE" "SIGMA" "SIGMAPX" "SIGMAPY" "SIGMAPZ" "SIGMAR" "SIGMAT" "SIGMAX" "SIGMAY" "SIGMAZ" "SIGX" "SIGY" "SIMBIN_CROSSOVER_NU" "SIMTMPDIR" "SLPTC" "SOL_SYNCH" "SPLIT" "SPTDUMPFREQ" "STATDUMPFREQ" "STEP" "STEPSPERTURN" "STOP" "STRING" "SUPERPOSE" "SYMMETRY" "T0" "TABLE" "TAU" "TELL" "TEMPLATEDIR" "TFALL" "THETA" "THIN" "THRESHOLD" "TIME" "TIMEINTEGRATOR" "TMULT" "TO" "TOL" "TOLERANCE" "TPULSEFWHM" "TRACE" "TRIMCOILTHRESHOLD" "TRISE" "TURNS" "TYPE" "UPPER" "UPPERBOUND" "VARIABLE" "VERIFY" "VERSION" "VKICK" "VMAX" "VMIN" "VOLT" "W" "WAKEF" "WARN" "WARP" "WEIGHT" "WIDTH" "WRITETOFILE" "X" "XEND" "XMA" "XMULT" "XSIZE" "XSTART" "Y" "YEND" "YMA" "YMULT" "YSIZE" "YSTART" "Z" "Z0" "ZEND" "ZINIT" "ZSTART" "ZSTOP") t) "\\>")

(defconst opal-font-lock-keywords-parameters
  (list
 '("\\<\\(A\\(?:DD\\|MPLITUDE_MODEL\\|NGLE\\|P\\(?:ERTURE\\|VETO\\)\\|T\\|UTOPHASE\\)\\|B\\(?:BOXINCR\\|C\\(?:FFT[TXYZ]\\|URRENT\\)\\|EAM_\\(?:\\(?:P\\(?:HI\\|R\\)\\|R\\)INIT\\)\\|FREQ\\|MAX\\|OUNDPDESTROYFQ\\|SCALE\\|Y\\)\\|C\\(?:A\\(?:LLS\\|THTEMP\\)\\|ENTRE\\|HARGE\\|L\\(?:ASS\\|EAR\\)\\|MD\\|O\\(?:EF\\(?:DENOM\\(?:PHI\\)?\\|NUM\\(?:PHI\\)?\\)\\|LUMN\\|N\\(?:DUCT\\|ST\\(?:RAINTS\\|_LENGTH\\)\\|V_HVOL_PROG\\)\\|RR[XYZ]\\)\\|ROSSOVER\\|SRDUMP\\|UTOFF\\(?:LONG\\|P[XYZ]\\|[RXY]\\)\\|YHARMON\\|ZERO\\)\\|D\\(?:ESIGNENERGY\\|K\\(?:1S\\|2S\\|[NS]R\\|[12NS]\\)\\|LAG\\|P\\(?:[HS]I\\)\\|THETA\\|UMP\\(?:_\\(?:DAT\\|FREQ\\|OFFSPRING\\)\\)?\\|V\\(?:ARS\\|OLT\\)\\|[STXYZ]\\)\\|E\\(?:BDUMP\\|CHO\\|KIN\\|L\\(?:ASER\\|EME\\(?:DGE\\|NT\\)\\)\\|MI\\(?:SSION\\(?:MODEL\\|STEPS\\)\\|TTED\\)\\|N\\(?:ABLEHDF5\\|D\\(?:SEQUENCE\\|_\\(?:NORMAL_[XY]\\|POSITION_[XY]\\)\\)\\|ERGY\\)\\|PSILON\\|SCALE\\|VERYSTEP\\|XP\\(?:ECTED_HYPERVOL\\|R\\)\\|[12TXY]\\)\\|F\\(?:E\\|GEOM\\|I\\(?:ELDMAPDIR\\|LE\\|NT\\)\\|MAPFN\\|NAME\\|ORM\\|R\\(?:EQ\\(?:UENCY_MODEL\\)?\\|OM\\)\\|STYPE\\|TOSC\\(?:AMPLITUDE\\|PERIODS\\)\\|ULL\\)\\|G\\(?:A\\(?:MMA\\|P\\(?:WIDTH\\)?\\)\\|E\\(?:\\(?:NE_MUTATION_PROBABILIT\\|OMETR\\)Y\\)\\|REENSF\\)\\|H\\(?:A\\(?:PERT\\|RMON\\(?:IC_NUMBER\\)?\\)\\|GAP\\|KICK\\|YPERVOLREFERENCE\\|[12]\\)\\|I\\(?:DEALIZED\\|MAGENAME\\|N\\(?:FO\\|ITIAL\\(?:\\(?:POPUL\\|_OPTIMIZ\\)ATION\\)\\|PUT\\(?:MOUNITS\\)?\\|TE\\(?:NSITYCUT\\|RPL\\)\\)\\|S_CLOSED\\|TSOLVER\\)\\|K\\(?:0S\\|1S\\|2S\\|3S\\|EYWORD\\|ICK\\|[0-3NS]\\)\\|L\\(?:A\\(?:G\\|SERPROFFN\\|TTICE_\\(?:\\(?:PHI\\|R\\|THETA\\)INIT\\)\\)\\|E\\(?:NGTH\\|VEL\\)\\|O\\(?:GBENDTRAJECTORY\\|WER\\(?:BOUND\\)?\\)\\)\\|M\\(?:A\\(?:SS\\|X\\(?:GENERATIONS\\|ITERS\\|STEPS\\|[RZ]\\)\\)\\|ETHOD\\|IN[RZ]\\|ODE\\|RE[XY]\\|SCAL[XY]\\|UTATION_PROBABILITY\\|[TXY]\\)\\|N\\(?:BIN\\|FREQ\\|LEFT\\|O\\|P\\(?:ART\\|OINTS\\)\\|RIGHT\\|UM\\(?:CELLS\\|_\\(?:COWORKERS\\|IND_GEN\\|MASTERS\\)\\)\\)\\|O\\(?:BJECTIVES\\|FFSET\\(?:P[XYZ]\\|[TXYZ]\\)\\|NE_PILOT_CONVERGE\\|P\\(?:CHARGE\\|MASS\\|YIELD\\)\\|R\\(?:DER\\|I\\(?:\\(?:ENTATIO\\|GI\\)N\\)\\)\\|UT\\(?:DIR\\|FN\\|PUT\\)\\)\\|P\\(?:A\\(?:R\\(?:FFT[TXY]\\|TICLE\\(?:MATTERINTERACTION\\)?\\)\\|TTERN\\)\\|DIS\\|HI\\(?:0\\|INIT\\|M\\(?:AX\\|IN\\)\\)?\\|LANE\\|OLYORDER\\|R\\(?:ECMODE\\|INIT\\)\\|S\\(?:DUMP\\(?:EACHTURN\\|FR\\(?:AME\\|EQ\\)\\)\\|I\\)\\|TC\\|\\(?:YMUL\\|Z\\(?:INI\\|MUL\\)\\)T\\|[1-4C]\\)\\|R\\(?:5[12]\\|6[12]\\|A\\(?:DIUS\\|N\\(?:DOM\\|GE\\)\\|STER\\)\\|E\\(?:COMBINATION_PROBABILITY\\|F\\(?:ER\\|POS\\)\\|PARTFREQ\\|SET\\)\\|F\\(?:FREQ\\|MAPFN\\|PHI\\)\\|INIT\\|M\\(?:AX\\|IN\\)\\|O\\(?:TATION\\|W\\)\\)\\|S\\(?:AMPLINGS\\|CALABLE\\|E\\(?:ED\\|LECTED\\|QUENCE\\)\\|I\\(?:G\\(?:MA\\(?:P[XYZ]\\|[RTXYZ]\\)?\\|[XY]\\)\\|M\\(?:BIN_CROSSOVER_NU\\|TMPDIR\\)\\)\\|LPTC\\|OL_SYNCH\\|P\\(?:LIT\\|TDUMPFREQ\\)\\|T\\(?:ATDUMPFREQ\\|EP\\(?:SPERTURN\\)?\\|OP\\|RING\\)\\|UPERPOSE\\|YMMETRY\\)\\|T\\(?:A\\(?:BLE\\|U\\)\\|E\\(?:LL\\|MPLATEDIR\\)\\|FALL\\|H\\(?:ETA\\|IN\\|RESHOLD\\)\\|IME\\(?:INTEGRATOR\\)?\\|MULT\\|OL\\(?:ERANCE\\)?\\|PULSEFWHM\\|R\\(?:ACE\\|I\\(?:MCOILTHRESHOLD\\|SE\\)\\)\\|URNS\\|YPE\\|[0O]\\)\\|UPPER\\(?:BOUND\\)?\\|V\\(?:ARIABLE\\|ER\\(?:IFY\\|SION\\)\\|KICK\\|M\\(?:AX\\|IN\\)\\|OLT\\)\\|W\\(?:A\\(?:KEF\\|R[NP]\\)\\|EIGHT\\|IDTH\\|RITETOFILE\\)\\|X\\(?:END\\|M\\(?:A\\|ULT\\)\\|S\\(?:IZE\\|TART\\)\\)\\|Y\\(?:END\\|M\\(?:A\\|ULT\\)\\|S\\(?:IZE\\|TART\\)\\)\\|Z\\(?:0\\|END\\|INIT\\|ST\\(?:ART\\|OP\\)\\)\\|[ABLSW-Z]\\)\\>"
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

;(concat "\\<" (regexp-opt '("AMR" "ANTIPROTON" "ASTRA" "CARBON" "CENTRE" "CLIGHT" "CMASS" "COLLIM" "DEGRAD" "DEUTERON" "DMASS" "E" "ELECTRON" "EMASS" "ENTRY" "EXIT" "FALSE" "FFT" "FFTPERIODIC" "HMINUS" "HMMASS" "INTEGRATED" "LONG-SHORT-RANGE" "MMASS" "MUON" "NONE" "NONEQUIL" "OPEN" "P0" "P3M" "PERIODIC" "PI" "PMASS" "POSITRON" "PROTON" "RADDEG" "SAAMG" "STANDARD" "TRANSV-SHORT-RANGE" "TRUE" "TWOPI" "UMASS" "URANIUM" "XEMASS" "XENON" "1D-CSR") t) "\\>")

(defconst opal-font-lock-keywords-constants
  (list
  '("\\<\\(1D-CSR\\|A\\(?:MR\\|NTIPROTON\\|STRA\\)\\|C\\(?:ARBON\\|ENTRE\\|LIGHT\\|MASS\\|OLLIM\\)\\|D\\(?:E\\(?:GRAD\\|UTERON\\)\\|MASS\\)\\|E\\(?:LECTRON\\|MASS\\|NTRY\\|XIT\\)?\\|F\\(?:ALSE\\|FT\\(?:PERIODIC\\)?\\)\\|HM\\(?:\\(?:INU\\|MAS\\)S\\)\\|INTEGRATED\\|LONG-SHORT-RANGE\\|M\\(?:MASS\\|UON\\)\\|NONE\\(?:QUIL\\)?\\|OPEN\\|P\\(?:3M\\|ERIODIC\\|MASS\\|\\(?:OSITR\\|ROT\\)ON\\|[0I]\\)\\|RADDEG\\|S\\(?:AAMG\\|TANDARD\\)\\|T\\(?:R\\(?:\\(?:ANSV-SHORT-RANG\\|U\\)E\\)\\|WOPI\\)\\|U\\(?:MASS\\|RANIUM\\)\\|XE\\(?:MASS\\|NON\\)\\)\\>"
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

;(concat "\\<" (regexp-opt '("ABS" "ACOS" "ASIN" "ATAN" "COS" "COSH" "EXP" "GAUSS" "LOG" "LOG10" "RANF" "SIN" "SINH" "SQRT" "TAN" "TANH" "TGAUSS") t) "\\>")

(defconst opal-font-lock-keywords-functions
  (list
  '("\\<\\(A\\(?:BS\\|COS\\|\\(?:SI\\|TA\\)N\\)\\|COSH?\\|EXP\\|GAUSS\\|LOG\\(?:10\\)?\\|RANF\\|S\\(?:INH?\\|QRT\\)\\|T\\(?:ANH?\\|GAUSS\\)\\)\\>"
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
     opal-font-lock-special_operators
     opal-font-lock-special_constants
     opal-font-lock-keywords-programflow
     opal-font-lock-keywords-simul
     opal-font-lock-keywords-controlstm
     opal-font-lock-keywords-elements
     opal-font-lock-keywords-beamspec
     opal-font-lock-keywords-matchingmet
     opal-font-lock-keywords-orbit_corr
     opal-font-lock-keywords-plot
     opal-font-lock-keywords-seqediting
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

    ; This is added so entity names with unde rscores can be more easily parsed
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