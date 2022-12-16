"""Test Line implementation"""

# we are testing the dunder operator, hence disable pylint warning
# pylint: disable=unnecessary-dunder-call

import unittest
import pyopal.elements.vertical_ffa_magnet
import pyopal.objects.line

class TestLine(unittest.TestCase):
    """Test Line implementation"""
    def setUp(self):
        """Define a few default variables"""
        self.magnet1 = pyopal.elements.vertical_ffa_magnet.VerticalFFAMagnet()
        self.magnet2 = pyopal.elements.vertical_ffa_magnet.VerticalFFAMagnet()
        self.line = pyopal.objects.line.Line()

    def test_initialisation(self):
        """Test that we can set member data okay"""
        self.line.origin = "bob"
        self.line.orientation = "fred"
        self.line.length = 1.0
        self.line.x = 2.0
        self.line.y = 3.0
        self.line.z = 4.0

        self.line.theta = 5.0
        self.line.phi = 6.0
        self.line.psi = 7.0
        self.assertEqual(self.line.origin, "bob")
        self.assertEqual(self.line.orientation, "fred")
        self.assertEqual(self.line.length, 1.0)
        self.assertEqual(self.line.x, 2.0)
        self.assertEqual(self.line.y, 3.0)
        self.assertEqual(self.line.z, 4.0)

        self.assertEqual(self.line.theta, 5.0)
        self.assertEqual(self.line.phi, 6.0)
        self.assertEqual(self.line.psi, 7.0)

        # Note, if debugging, this is the name of the Opal/Lines/Line.h line
        # that is defined in PyLine::register()
        self.line.set_opal_name("A OPAL NAME")
        self.assertEqual(self.line.get_opal_name(), "A OPAL NAME")


    def test_append(self):
        """Check that we can append elements"""
        self.line.append(self.magnet1)
        self.assertEqual(self.line[0], self.magnet1)
        self.line.append(self.magnet2)
        self.assertEqual(self.line[1], self.magnet2)
        self.line.append(self.magnet1)
        self.assertEqual(self.line[2], self.magnet1)
        with self.assertRaises(AttributeError):
            self.line.append("not an element")
 
    def test_set(self):
        """Check that we can set elements"""
        self.line.append(self.magnet1)
        self.assertEqual(self.line[0], self.magnet1)
        self.line.append(self.magnet2)
        self.assertEqual(self.line[1], self.magnet2)
        self.line[1]  = self.magnet1
        self.assertEqual(self.line[1], self.magnet1)
        self.line.__setitem__(1, self.magnet2)
        self.assertEqual(self.line[1], self.magnet2)
        self.line.__setitem__(-1, self.magnet1)
        self.assertEqual(self.line[1], self.magnet1)
        self.line.__setitem__(-2, self.magnet2)
        self.assertEqual(self.line[0], self.magnet2)
        with self.assertRaises(AttributeError):
            self.line[1]  = "not an element"
        with self.assertRaises(RuntimeError):
            self.line.__setitem__(2, self.magnet2)

    def test_get(self):
        """Check that we can get elements"""
        with self.assertRaises(RuntimeError):
            self.line[0]
        self.line.append(self.magnet1)
        self.line.append(self.magnet2)
        self.assertEqual(self.line.__getitem__(0), self.magnet1)
        self.assertEqual(self.line.__getitem__(1), self.magnet2)
        self.assertEqual(self.line.__getitem__(-1), self.magnet2)
        self.assertEqual(self.line.__getitem__(-2), self.magnet1)
        with self.assertRaises(RuntimeError):
            self.line.__getitem__(3)

    def test_len(self):
        """Check len function"""
        self.assertEqual(len(self.line), 0)
        self.line.append(self.magnet1)
        self.line.append(self.magnet2)
        self.assertEqual(len(self.line), 2)

if __name__ == "__main__":
    unittest.main()
