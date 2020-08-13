
import unittest
from configparser import SafeConfigParser
from pybemt.rotor import Rotor, Section
from pybemt.airfoil import load_airfoil
from pybemt.fluid import Fluid
from math import pi,radians
from scipy import optimize

class TestRotor(unittest.TestCase):
    def setUp(self):
        cfg = SafeConfigParser()
        cfg.read('test_config.ini')

        self.rotor = Rotor(cfg, 'rotor', 'rotor')

    def test_rotor_parameters(self):
        self.assertEqual(2, self.rotor.n_blades)
        self.assertEqual(0.8, self.rotor.diameter)
        self.assertEqual(8, len(self.rotor.sections))

    def test_precalc(self):
        self.rotor.blade_radius = 0.0
        self.rotor.area = 0.0
        twist = 1.0
        self.rotor.precalc(twist=twist)
        self.assertEqual(0.4, self.rotor.blade_radius)
        self.assertAlmostEqual(0.50265482, self.rotor.area, places=8)
        self.assertAlmostEqual(radians(19.6+twist), self.rotor.sections[0].pitch, places=8)
        self.assertAlmostEqual(radians(twist), self.rotor.sections[-1].pitch, places=8)

    def test_dataframe(self):
        df = self.rotor.sections_dataframe()
        columns = ['radius','chord','pitch','Cl','Cd','dT','dQ','F','a','ap','Re']
        for c in columns:
            self.assertTrue(c in df.columns)

        chords = [0.056, 0.07, 0.07, 0.065, 0.058, 0.05, 0.043, 0.01]
        for i,c in enumerate(chords):
            self.assertEqual(c, df['chord'][i])


class TestSection(unittest.TestCase):
    def setUp(self):
        cfg = SafeConfigParser()
        cfg.read('test_config.ini')
        # Correct phi and induction factors for hand-calculated case
        self.v = 1.0
        self.omega = 2000*2*pi/60.
        self.a = 1.017166762893497
        self.ap = 0.002340435194091869
        self.phi = 0.03216841958789685
        self.F = 0.9999799142439774

        self.rotor = Rotor(cfg, 'rotor', 'rotor')
        self.sec = self.rotor.sections[-1]

        self.fluid = Fluid(cfg)

    def test_section_parameters(self):
        self.assertEqual('GOE_408', self.sec.airfoil.name)
        self.assertEqual(1, self.sec.C)

    def test_precalc(self):
        self.sec.sigma = 0
        self.sec.precalc()
        self.assertAlmostEqual(1.0/(30*pi), self.sec.sigma, places=6)

    def test_tiploss(self):
        self.assertEqual(1, self.sec.tip_loss(0))
        self.assertEqual(1, self.sec.tip_loss(0.001))
        self.assertEqual(1, self.sec.tip_loss(pi))
        self.assertEqual(0.0, self.sec.tip_loss(-pi/10))
        self.assertAlmostEqual(0.7521512, self.sec.tip_loss(pi/10), places=7)
        self.assertAlmostEqual(0.9999799142439774, self.sec.tip_loss(self.phi), places=7)

    def test_airfoil_forces(self):
        # Test against hand-calculated values for phi=0.
        Cl = 0.4002
        Cd = 0.0220
        CT,CQ = self.sec.airfoil_forces(0)
        self.assertAlmostEqual(Cl, CT)
        self.assertAlmostEqual(Cd, CQ)
        # Test turbine mode
        self.sec.C = -1
        CT,CQ = self.sec.airfoil_forces(0)
        self.assertAlmostEqual(Cl, CT)
        self.assertAlmostEqual(-Cd, CQ)

    def test_induction_factors(self):
        # Test against hand-calculated values for phi=0.
        CT = 0.4002
        CQ = 0.0220
        a, ap = self.sec.induction_factors(0)
        self.assertAlmostEqual(-1, a)
        self.assertAlmostEqual(1, ap)

        a, ap = self.sec.induction_factors(self.phi)
        self.assertAlmostEqual(self.a, a)
        self.assertAlmostEqual(self.ap, ap)
        
    def test_func(self):
        # Should give 0 in residual
        self.assertAlmostEqual(0.0, self.sec.func(self.phi, self.v, self.omega))

    def test_forces(self):
        # Calculating forces using momentum theory should give same results as the implemented blade element theory
        dT, dQ = self.sec.forces(self.phi, self.v, self.omega, self.fluid)

        dTmom = 4*pi*self.fluid.rho*self.sec.radius*self.v**2*(1 + self.sec.C*self.a)*self.a*self.F*self.sec.width
        dQmom = 4*pi*self.fluid.rho*self.sec.radius**3*self.v*(1 + self.sec.C*self.a)*self.ap*self.omega*self.F*self.sec.width

        self.assertAlmostEqual(dTmom, dT)
        self.assertAlmostEqual(dQmom, dQ)


