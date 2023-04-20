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

"""Test the field module"""
import unittest
import tempfile
import sys
import os
import pyopal.objects.parser
import pyopal.objects.field
import pyopal.objects.encapsulated_test_case

class FieldTest(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Test that we can get fields out"""
    def setUp(self):
        """Dump string to a temporary file for use in testing"""
        self.write_temp(self.good_lattice)

    def write_temp(self, a_string):
        self.temp = tempfile.NamedTemporaryFile(mode='w+')
        self.temp.write(a_string)
        self.temp.flush()
        self.file_name = self.temp.name

    def encapsulated_test_field(self):
        """Test that we can get out a field value"""
        print("TEST FIELD")
        pyopal.objects.parser.initialise_from_opal_file(self.file_name)
        float_tolerance = 1e-3
        test_fail = False
        print("Test result")

        value = pyopal.objects.field.get_field_value(1, 0, 0, 0)
        reference = (False, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0)
        print("Field", value, "reference", reference)
        for i, ref_value in enumerate(reference):
            test_fail = abs(value[i]-ref_value) > float_tolerance

        value = pyopal.objects.field.get_number_of_elements()
        ref_value = 2
        test_fail = test_fail or value != ref_value
        print("Number of elements", value, "reference", ref_value, "Test", test_fail)

        value = pyopal.objects.field.get_element_name(1)
        ref_value = "OFFSET_TEST"
        test_fail = test_fail or value != ref_value
        print("Element name", value, "reference", ref_value, "Test", test_fail)


        value = pyopal.objects.field.get_element_start_position(1)
        ref_value = [1.0, 10.0, 0.0]
        for i in range(3):
            test_fail = test_fail or abs(value[i] - ref_value[i]) > float_tolerance
        print("Offset start", value, "reference", ref_value, "Test", test_fail)

        value = pyopal.objects.field.get_element_start_normal(1)
        ref_value = [0.0, 1.0, 0.0]
        for i in range(3):
            test_fail = test_fail or abs(value[i] - ref_value[i]) > float_tolerance
        print("Offset start normal", value, "reference", ref_value, "Test", test_fail)

        value = pyopal.objects.field.get_element_end_position(1)
        ref_value = [0.0, 11.0, 0.0]
        for i in range(3):
            test_fail = test_fail or abs(value[i] - ref_value[i]) > float_tolerance
        print("Offset end position", value, "reference", ref_value, "Test", test_fail)

        value = pyopal.objects.field.get_element_end_normal(1)
        ref_value = [2**-0.5, 2**-0.5, 0.0]
        for i in range(3):
            test_fail = test_fail or abs(value[i] - ref_value[i]) > float_tolerance
        print("Offset end normal", value, "reference", ref_value, "Test", test_fail)

        return test_fail

    debug = False
    command = """
import pyopal.objects.parser
pyopal.objects.parser.initialise_from_opal_file(
"""
    good_lattice = """
///////////////////////////////////////
Title,string="Dummy lattice for testing";
///

OPTION, VERSION=20210100;

field: VERTICALFFAMAGNET,
        B0=2,
        FIELD_INDEX=2,
        MAX_HORIZONTAL_POWER=4,
        END_LENGTH=0.01,
        CENTRE_LENGTH=1000,
        HEIGHT_NEG_EXTENT=4,
        HEIGHT_POS_EXTENT=4,
        WIDTH=2000,
        BB_LENGTH=10;

offset_test: LOCAL_CARTESIAN_OFFSET, END_POSITION_X=-1.0, END_POSITION_Y=1.0,
                                END_NORMAL_X=2,      END_NORMAL_Y=2;

ringdef: RINGDEFINITION, HARMONIC_NUMBER=1, LAT_RINIT=1, LAT_PHIINIT=0,
         LAT_THETAINIT=0.0, BEAM_PHIINIT=0, BEAM_PRINIT=0,
         BEAM_RINIT=1, SYMMETRY=1, RFFREQ=1, IS_CLOSED=false,
         MIN_R=0.1, MAX_R=2;

l1: Line = (ringdef, field, offset_test);

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
    FieldTest.debug = False
    unittest.main()
