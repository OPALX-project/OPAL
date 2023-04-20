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

"""Test that distribution parses okay"""
import unittest
import pyopal.objects.distribution

class TestDistribution(pyopal.objects.encapsulated_test_case.EncapsulatedTestCase):
    """Test that distribution parses okay"""
    def encapsulated_test_init(self):
        """Check we can initialise variables without a problem"""
        my_distribution = pyopal.objects.distribution.Distribution()
        my_distribution.type = "cheese"
        self.assertEqual(my_distribution.type, "CHEESE")
        my_distribution.fname = "disttest.dat"
        self.assertEqual(my_distribution.fname, "disttest.dat")
        my_distribution.momentum_units = "asdasd"
        self.assertEqual(my_distribution.momentum_units, "ASDASD")

    def encapsulated_test_register(self):
        """Check we can call register function"""
        my_distribution = pyopal.objects.distribution.Distribution()
        my_distribution.register()

if __name__ == "__main__":
    unittest.main()
