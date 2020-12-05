import os
import unittest

import pyopal
import pyopal.parser
import pyopal.field

class TestField(unittest.TestCase):
    def setUp(self):
        self.start = os.getcwd()
        test_dir = os.path.split(__file__)[0]
        self.lattice_dir = os.path.join(test_dir, "dummy_lattice")
        os.chdir(self.lattice_dir)
        self.lattice_fname = "dummy_lattice.opal"

    def tearDown(self):
        os.chdir(self.start)

    def test_a_throws_on_no_lattice(self):
        # assumes that this test is run first; default TestLoader behaviour is
        # sort by method name, so that is okay here.
        try:
            field = pyopal.field.get_field_value(0.5, 5, 0, 0)
            self.assertTrue(false, "Should have thrown an exception")
        except RuntimeError:
            pass # Lattice not initialised - correctly threw error

    def test_initialise_from_opal_file(self):
        pyopal.parser.initialise_from_opal_file(self.lattice_fname)

        field = pyopal.field.get_field_value(0.5, 5, 0, 0)
        self.assertFalse(field[0]) # inBounds
        self.assertAlmostEqual(field[1], 0.0, 1e-9)
        self.assertAlmostEqual(field[2], 0.0, 1e-9)
        self.assertAlmostEqual(field[3], 0.46211715726000974, 1e-9)
        self.assertAlmostEqual(field[4], 0.0, 1e-9)
        self.assertAlmostEqual(field[5], 0.0, 1e-9)
        self.assertAlmostEqual(field[6], 0.0, 1e-9)

        # Really I am testing the boost framework here; normally this is not
        # necessary, I am just doing it as a sort of experiment.
        try:
            field = pyopal.field.get_field_value("fish", 5, 0, 0)
            self.assertTrue(false, "Should have thrown an exception")
        except TypeError:
            pass # bad value - correctly threw an error

if __name__ == "__main__":
	unittest.main() # for running directly from the command line
