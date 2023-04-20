# Routine to encapsulate a test case, so that all memory initialisations are done
# in an forked process
#
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
#

import unittest
import sys
import os

class EncapsulatedTestCase(unittest.TestCase):
    """
    Because OPAL is sloppy with memory, we run many tests in a forked process.
    This is done by using EncapsulatedTestCase instead of the regular unittest
    TestCase.

    TestCases should inherit from this class. By default make install will put
    the class in pyopal.objects. Individual test methods should have a name like
    "encapsulated_test_<name>" to indicate that the method needs to be run in a
    forked process. EncapsulatedTestCase will search for these tests at run time
    and execute them.
    """
    def test_run_encapsulated(self):
        """
        Search through the methods and look for functions with name
        "encapsulated_test_<name>". If one is found, it will be executed in a
        forked process.

        This method replaces all of the encapsulated_test_<name> methods, so if
        multiple iterations fail it still gets counted as a single failure.
        """
        for name, test_method in vars(self.__class__).items():
            if name[:18] != "encapsulated_test_":
                continue
            if not callable(test_method):
                continue
            a_pid = os.fork()
            if a_pid == 0: # the child process
                try:
                    retvalue = test_method(self)
                except:
                    print("Encapsulated test failed with:")
                    sys.excepthook(*sys.exc_info())
                    os._exit(1)
                os._exit(0)
            else:
                return_value = os.waitpid(a_pid, 0)[1]
                classname = self.__class__.__name__
                self.assertEqual(return_value, 0, 
                                 msg=f"Failed test {classname}.{name}")

