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

"""Module to drive pylint tests"""

import tempfile
import os
import unittest
import subprocess
import inspect
import argparse
import pyopal

class PyLintTest(unittest.TestCase):
    """
    Because there is no actual python code in PyOpal, we use a pylint run
    against the tests as a proxy for checking that code is reasonably pythonic.
    The aim is to ensure that devs understand the basic principles of python
    style and have done a reasonable effort of implementing them - but it is not
    required to be religious about it. Hence require 8.0/10.0 as a score.
    """
    def setUp(self):
        """set up the test"""
        self.log = tempfile.NamedTemporaryFile("w+")
        self.opal_test_path = \
                        os.path.dirname(inspect.getfile(inspect.currentframe()))
        self.opal_install_path = os.path.abspath(os.path.dirname(pyopal.__file__))
        self.pass_score = 8.0
        self.disables = ["missing-module-docstring", "missing-class-docstring"]
        self.pylint_cmd = ["pylint"]

    def run_pylint(self, dirpath, fname):
        """Run pylint and fill the log tempfile"""
        dis_list = []
        for disable_item in self.disables:
            dis_list += ["--disable", disable_item]
        full_path = os.path.join(dirpath, fname)
        subprocess.run(self.pylint_cmd+dis_list+[full_path],
            stdout=self.log, stderr=subprocess.STDOUT)

    def get_score_from_line(self, line):
        """Extract the pylint score by scraping the text output"""
        if "Your code has been rated at" not in line:
            return -1
        score = line.split("rated at ")[1:]
        score = score[0].split("/10")[0]
        score = float(score)
        return score

    def check_logfile(self):
        """Read the log and look for errors"""
        self.log.flush()
        self.log.seek(0)
        buffer = ""
        for line in self.log.readlines():
            buffer += line
            score = self.get_score_from_line(line)
            if score < 0:
                continue
            msg=f"""

Failed static code analysis with output

{buffer}

Score must be greater than {self.pass_score} to pass. Try running pylint to check manually.
            """
            self.assertGreater(score, self.pass_score, msg)
            buffer = ""

    def will_not_run_arg(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--do_not_run_pylint", action='store_true', help="disable pylint execution")
        args = parser.parse_args()
        return args.do_not_run_pylint

    def test_check_pylint_tests(self):
        """Recurse up the directory structure and run pylint on each .py file in
        the tests"""
        if self.will_not_run_arg():
            print("Pylint execution disabled - skipping")
            return
        for path_root in [self.opal_install_path, self.opal_test_path]:
            print(f"Checking pylint starting at {path_root}")
            for dirpath, dirnames, filenames in os.walk(top=path_root):
                print("  Entering", dirpath)
                for fname in filenames:
                    if fname[-3:] != ".py":
                        continue
                    print("    Testing", os.path.join(dirpath, fname))
                    self.log.truncate(0)
                    self.run_pylint(dirpath, fname)
                    self.check_logfile()

    verbose = False

if __name__ == "__main__":
    unittest.main()
