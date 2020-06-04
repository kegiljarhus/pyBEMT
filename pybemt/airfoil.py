# -*- coding: utf-8 -*-

"""
Module for holding airfoil data, and providing drag and lift coefficients to the solver.

Airfoil data is stored in the folder pybemt/airfoils. Currently, the Aerodyn format
is supported, with only a single airfoil table. Data must be available from -180 to 180 degrees,
and a quadratic function is built to interpolate between data in the airfoil table.

This module can also be executed from the command-line to plot drag and lift coefficients 
for a single airfoil, e.g. 

.. code-block:: console

    python airfoil.py NACA_4412

"""

import os
import sys
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from math import degrees, radians, atan2, sin, cos


class Airfoil:
    """
    Class for storing airfoil drag and lift coefficients. Should be initialized using 
    the load_airfoil() function.
    """

    def __init__(self):
        self.name = None
        self.alpha_ = None
        self.Cl_ = None
        self.Cd_ = None
        self.Cl_func = None
        self.Cd_func = None
        self.zero_lift = 0.0
        
    def _normalize_angle(self, alpha):
        """
        Ensure that the angle fulfils :math:`\pi < \alpha < \pi`

        :param float alpha: Angle in radians
        :return: Normalized angle
        :rtype: float
        """

        return atan2(sin(alpha), cos(alpha))
    
    def Cd(self, alpha): 
        """
        Provide drag coefficent for a given angle of attack.

        :param float alpha: Angle in radians
        :return: Drag coefficient
        :rtype: float
        """
        # NB! The stored data uses degrees for the angle
        return self.Cd_func(-self.zero_lift + degrees(self._normalize_angle(alpha)))

    def Cl(self, alpha): # NB! The stored data uses degrees for the angle
        """
        Provide drag coefficent for a given angle of attack.

        :param float alpha: Angle in radians
        :return: Drag coefficient
        :rtype: float
        """
        return self.Cl_func(-self.zero_lift + degrees(self._normalize_angle(alpha)))

    def plot(self, color='k'):
        """
        Plot lift and drag coefficients.

        :param string color: Matplotlib color key. Default value is k, i.e. black.
        """

        pl.plot(self.alpha_, self.Cl_, color + '-')
        pl.plot(self.alpha_, self.Cd_, color + '--')
        pl.title('Airfoil characteristics for ' + self.name)
        pl.xlabel('Angle of attack')
        pl.ylabel('Drag and lift coefficients')
        pl.legend(('$C_l$','$C_d$'))
        
        
def load_airfoil(name):
    """
    Load airfoil data from data file into an Airfoil object. 

    :param string name: name of the airfoil to load, e.g. 'NACA_4412'.

    :return: Airfoil object
    :rtype: pybemt.Airfoil
    """
    a = Airfoil()
    
    this_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(this_dir, 'airfoils', name + '.dat')
    
    a.name = name
    
    a.alpha_, a.Cl_, a.Cd_ = np.loadtxt(path, skiprows=14, unpack=True)
    
    a.Cl_func = interp1d(a.alpha_, a.Cl_, kind='quadratic')
    a.Cd_func = interp1d(a.alpha_, a.Cd_, kind='quadratic')
            
    return a


if __name__ == '__main__':
    # Load and plot airfoil from command line
    name = sys.argv[-1]
    a = load_airfoil(name)
    a.plot()

    pl.show()
