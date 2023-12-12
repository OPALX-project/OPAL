# Copyright (c) 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
#
# This file is part of OPAL.
#
# OPAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with OPAL.  If not, see <https://www.gnu.org/licenses/>.

import unittest
import tempfile
import subprocess
import pyopal.objects.parser
import pyopal.objects.encapsulated_test_case


class ParserTest(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    def make_temp(self, a_string):
        """Dump string to a temporary file for use in testing"""
        my_temp = tempfile.NamedTemporaryFile(mode='w+')
        my_temp.write(a_string)
        my_temp.flush()
        return my_temp

    def encapsulated_test_parser(self):
        """OPAL can kill python execution so we hide the test in a subprocess"""
        temp_file = self.make_temp(self.good_lattice)
        pyopal.objects.parser.initialise_from_opal_file(temp_file.name)

    command = """
import pyopal.objects.parser
pyopal.objects.parser.initialise_from_opal_file(
"""
    good_lattice = """
///////////////////////////////////////
Title,string="Dummy lattice for testing";
///

OPTION, VERSION=20210100;
OPTION, ECHO=False;
OPTION, INFO=False;
OPTION, WARN=False;
OPTION, STATDUMPFREQ = 10; 
OPTION, AUTOPHASE = 4; 

Title, string="OPAL-T Parser Test";

//Initial Momentum Calculation 
REAL Edes    = 1.4e-9; //initial energy in GeV 
REAL gamma   = (Edes+EMASS)/EMASS; 
REAL beta    = sqrt(1-(1/gamma^2)); 
REAL P0      = gamma*beta*EMASS;    //inital z momentum

REAL rf_freq             = 1.3e3;     //RF frequency. (Hz)
REAL n_particles         = 1E2;       //Number of particles in simulation.
REAL beam_bunch_charge   = 1e-9;      //Charge of bunch. (C)
REAL beam_current = beam_bunch_charge*rf_freq*1e6;

dr1: DRIFT, L=1.0, ELEMEDGE=0.5; 

l1: Line = (dr1);

Dist1: DISTRIBUTION, TYPE=FLATTOP, 
	SIGMAX = 0.0005, 
	SIGMAY = 0.0005, 
	EMITTED= True;

Fs1:FIELDSOLVER, FSTYPE=None, MX=5, MY=5, MT=5,
                 PARFFTX=true, PARFFTY=false, PARFFTT=false,
                 BCFFTX=open, BCFFTY=open, BCFFTT=open,BBOXINCR=1;

beam1: BEAM, PARTICLE=ELECTRON, pc=P0, NPART=n_particles, BCURRENT=beam_current, CHARGE=-1.0, BFREQ=rf_freq;

TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=1, DT=1e-10, ZSTOP=0.5;
RUN, METHOD="PARALLEL-T", BEAM=beam1, FIELDSOLVER=Fs1, DISTRIBUTION=Dist1;
ENDTRACK;
STOP;
"""

if __name__ == "__main__":
    unittest.main()
