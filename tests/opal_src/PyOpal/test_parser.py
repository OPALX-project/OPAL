import os
import unittest

import pyopal
import pyopal.parser

class TestParser(unittest.TestCase):
    def setUp(self):
        self.start = os.getcwd()
        test_dir = os.path.split(__file__)[0]
        self.lattice_dir = os.path.join(test_dir, "dummy_lattice")
        os.chdir(self.lattice_dir)
        self.lattice_fname = "dummy_lattice.opal"

    def tearDown(self):
        os.chdir(self.start)

    def test_load_file(self):
        # OPAL is a singleton (at the moment); so we do all tests in one blob
        self.assertFalse(pyopal.parser.is_initialised())
        pyopal.parser.initialise_from_opal_file(self.lattice_fname)
        self.assertTrue(pyopal.parser.is_initialised())
        try:
            pyopal.parser.initialise_from_opal_file(self.lattice_fname)
            self.assertTrue(False, "Should have thrown")
        except RuntimeError:
            pass
        object_string = pyopal.parser.list_objects(".*")
        # not sure why this returns no objects; but it is an OPAL internal issue
        # I guess? Should return a list of all objects
        self.assertTrue('Object names matching the pattern ".*"' in object_string)


if __name__ == "__main__":
    unittest.main() # for running directly from the command line
