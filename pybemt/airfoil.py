# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from math import degrees, radians, atan2, sin, cos



class Airfoil:
    def __init__(self):
        self.name = None
        self.alpha_ = None
        self.Cl_ = None
        self.Cd_ = None
        self.Cl_func = None
        self.Cd_func = None
        self.zero_lift = 0.0
        
    def normalize_angle(self, alpha):
        return atan2(sin(alpha), cos(alpha))
    
    def Cd(self, alpha): # NB! The stored data uses degrees for the angle
        return self.Cd_func(-self.zero_lift + degrees(self.normalize_angle(alpha)))

    def Cl(self, alpha): # NB! The stored data uses degrees for the angle
        return self.Cl_func(-self.zero_lift + degrees(self.normalize_angle(alpha)))

    def plot(self,col='k'):
        pl.plot(self.alpha_, self.Cl_, col+'-')
        pl.plot(self.alpha_, self.Cd_, col+'--')
        pl.title('Airfoil characteristics for ' + self.name)
        pl.xlabel('Angle of attack')
        pl.ylabel('Drag and lift coefficients')
        pl.legend(('$C_l$','$C_d$'))
        #alphas = np.linspace(self.alpha_[0], self.alpha_[-1], 2*len(self.alpha_))
        #pl.plot(alphas, [self.Cl(radians(a)) for a in alphas], 'C1--')
        #pl.plot(alphas, [self.Cd(radians(a)) for a in alphas], 'C2--')
        
        
def load_airfoil(name):
    
    a = Airfoil()
    
    this_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(this_dir, 'airfoils', name + '.dat')
    
    a.name = name
    
    a.alpha_, a.Cl_, a.Cd_ = np.loadtxt(path, skiprows=14, unpack=True)
    
    a.Cl_func = interp1d(a.alpha_, a.Cl_, kind='quadratic')
    a.Cd_func = interp1d(a.alpha_, a.Cd_, kind='quadratic')

    for i in range(5,len(a.alpha_)):
        C0 = a.Cl_[i-1]
        C1 = a.Cl_[i]
        if C0 < 0 and C1 >= 0:
            a0 = a.alpha_[i-1]
            a1 = a.alpha_[i]
            a.zero_lift =0.0 # a0 + (0.0 - C0)*(a1 - a0)/(C1 - C0)
            break
            
    # print('zero',a.zero_lift)
    
    return a


if __name__ == '__main__':
    #name = sys.argv[-1]
    for i,name in enumerate(sys.argv[1:]):
        print('name',name)
        a = load_airfoil(name)
        a.plot('C%i'%i)

    pl.show()
