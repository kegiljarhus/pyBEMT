
import unittest
from configparser import SafeConfigParser
from pybemt.fluid import Fluid

class TestFluid(unittest.TestCase):
    def setUp(self):
        cfg = SafeConfigParser()
        cfg.read('test_config.ini')

        self.fluid = Fluid(cfg)

    def test_fluid_parameters(self):
        self.assertEqual(1.225, self.fluid.rho)
        self.assertEqual(1.81e-5, self.fluid.mu)
