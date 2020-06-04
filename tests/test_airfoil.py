
import unittest

from pybemt.airfoil import Airfoil, load_airfoil
from math import radians,pi,degrees
import numpy as np
import matplotlib.pyplot as pl

class TestAirfoil(unittest.TestCase):

    def setUp(self):
        self.a = load_airfoil('NACA_4412')
        self.alphas = np.array([radians(a) for a in [-180,-10,1,10,90]])

    def test_load(self):
        self.assertEqual('NACA_4412', self.a.name)

    def test_normalize(self):
        """
        Test that angles are kept within -pi,pi range
        """
        degs = [-270,360,450,540,720]
        rads_correct = [pi/2,0,pi/2,pi,0]
        for alpha,alpha2 in zip(degs, rads_correct):
            self.assertAlmostEqual(alpha2, self.a._normalize_angle(radians(alpha)), places=5)

    def test_drag(self):
        """
        Test that correct drag coefficients are given.
        """
        Cds = [0.0060, 0.0431, 0.0180, 0.0368, 1.7759] 
        for alpha,Cd in zip(self.alphas, Cds):
            self.assertAlmostEqual(Cd, self.a.Cd(alpha), places=5)

    def test_lift(self):
        """
        Test that correct lift coefficients are given.
        """
        Cls = [-0.0922,-0.4011,0.4939,1.1991,0.2213] 
        for alpha,Cl in zip(self.alphas, Cls):
            self.assertAlmostEqual(Cl, self.a.Cl(alpha), places=5)

    def test_interpolation(self):
        """
        Test interpolated drag and lift values at zero angle of attack.
        """
        self.assertAlmostEqual(0.01887, self.a.Cd(0), places=5)
        self.assertAlmostEqual(0.38715, self.a.Cl(0), places=5)

    def test_plot(self):
        """
        Test that airfoil data is plotted.
        """
        self.a.plot()
        ax = pl.gca()
        x,y = self.a.alpha_, self.a.Cl_
        x_plot, y_plot = ax.lines[0].get_xydata().T
        self.assertIsNone(np.testing.assert_array_equal(y_plot, y))
        y = self.a.Cd_
        x_plot, y_plot = ax.lines[1].get_xydata().T
        self.assertIsNone(np.testing.assert_array_equal(y_plot, y))

if __name__ == '__main__':
    unittest.main()
