import os
import unittest
import pyopal.objects.minimal_runner


class MinimalRunnerOverload(pyopal.objects.minimal_runner.MinimalRunner):
    def __init__(self):
        """Initialise"""
        super().__init__()
        self.preproc = os.path.join(self.tmp_dir, "preproc")
        self.postproc = os.path.join(self.tmp_dir, "postproc")

    def make_element_iterable(self):
        """Iterable is added to the line"""
        return [self.null_drift(), self.null_drift()]

    def preprocess(self):
        """Preprocess is called just before tracking"""
        # check that the element iterable was added to the line
        if len(self.line) != 4: # RingDef, minimal runner null drift, two more null drifts
            raise RuntimeError("Expected line of length 3, was actually "+str(len(self.line)))

        with open(self.preproc, "a") as fout:
            fout.write("a")

    def postprocess(self):
        """Preprocess should be called just after tracking"""
        with open(self.postproc, "a") as fout:
            fout.write("b")



class MinimalRunnerTest(unittest.TestCase):
    def test_minimal_runner(self):
        """Test that minimal runner calls overloaded functions okay"""
        runner = MinimalRunnerOverload()
        for i in range(5):
            runner.execute_fork()
        # check that the runner preproc and postproc were called
        with open(runner.preproc) as fin:
            self.assertEqual(len(fin.read()), 5)
        with open(runner.postproc) as fin:
            self.assertEqual(len(fin.read()), 5)
 
if __name__ == "__main__":
    unittest.main()


