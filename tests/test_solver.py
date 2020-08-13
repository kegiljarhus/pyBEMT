
import unittest
from pybemt.solver import Solver

class TestSolver(unittest.TestCase):
    def setUp(self):

        self.solver = Solver('test_config.ini')
        self.T,self.Q,self.P,self.df = self.solver.run()

    def test_solvers(self):
        # Brute and bisect solver should give roughly same answer
        self.solver.solver = 'brute'
        Tb,Qb,Pb,df = self.solver.run()
        self.assertAlmostEqual(self.T/self.P, Tb/Pb, places=3)

    def test_coeffs(self):
        J,CT,CQ,CP,eta = self.solver.rotor_coeffs(self.T, self.Q, self.P)
        self.assertAlmostEqual(0.0375, J)
        self.assertAlmostEqual(0.0089192, CP)
 
