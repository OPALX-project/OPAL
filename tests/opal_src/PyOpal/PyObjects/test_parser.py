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

null: LOCAL_CARTESIAN_OFFSET,
                end_position_x=0., end_position_y=0.,
                end_normal_x=1.0, end_normal_y=0.;

ringdef: RINGDEFINITION, HARMONIC_NUMBER=1, LAT_RINIT=1, LAT_PHIINIT=0,
         LAT_THETAINIT=0.0, BEAM_PHIINIT=0, BEAM_PRINIT=0,
         BEAM_RINIT=1, SYMMETRY=1, RFFREQ=1, IS_CLOSED=false,
         MIN_R=0.1, MAX_R=2;

l1: Line = (ringdef, null);

Dist1: DISTRIBUTION, TYPE=gauss;

Fs1:FIELDSOLVER, FSTYPE=None, MX=5, MY=5, MT=5,
                 PARFFTX=true, PARFFTY=false, PARFFTT=false,
                 BCFFTX=open, BCFFTY=open, BCFFTT=open,BBOXINCR=2;

beam1: BEAM, PARTICLE=PROTON, pc=1, NPART=1, BCURRENT=1.6e-19, CHARGE=1.0, BFREQ=1;

TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=1, STEPSPERTURN=1;
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=Fs1, DISTRIBUTION=Dist1;
ENDTRACK;
STOP;
"""

if __name__ == "__main__":
    unittest.main()
