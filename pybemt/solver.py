# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
from scipy import optimize
from configparser import SafeConfigParser
from math import radians, degrees, sqrt, pi
from .fluid import Fluid
from .rotor import Rotor


class Solver: 
    def __init__(self, config_file):

        # Read configuration file
        cfg = SafeConfigParser()
        cfg.read(config_file)
        
        # Case
        self.v_inf = cfg.getfloat('case', 'v_inf')
        self.rpm = cfg.getfloat('case', 'rpm')
        if cfg.has_option('case', 'coaxial'):
            self.coaxial = cfg.getboolean('case', 'coaxial')
        else:
            self.coaxial = False
        

        # Rotor
        if cfg.has_section('turbine'):
            self.mode = 'turbine'
            self.rotor = Rotor(cfg, 'turbine', self.mode)
        else:
            self.mode = 'rotor'
            self.rotor = Rotor(cfg, 'rotor', self.mode)

        # Fluid
        self.fluid = Fluid(cfg)
        
        # Output
        self.T = 0 # Thrust
        self.Q = 0 # Torque
        self.P = 0 # Power
        
        # Coaxial
        if self.coaxial:
            self.rpm2 = cfg.getfloat('case','rpm2')
            self.rotor2 = Rotor(cfg, 'rotor2', self.mode)
            self.zD = cfg.getfloat('case','dz')/self.rotor.diameter
            self.T2 = 0
            self.Q2 = 0
            self.P2 = 0

        # Solver
        self.solver = 'bisect'
        self.Cs = 0.625
        if cfg.has_section('solver'):
            self.solver = cfg.get('solver','solver')
            if cfg.has_option('solver', 'Cs'):
                self.Cs = cfg.getfloat('solver','Cs')
       
    def rotor_coeffs(self, T, Q, P):
        """ Calculate non-dimensional coefficients for rotor. 
        """
        D = self.rotor.diameter
        R = 0.5*D
        rho = self.fluid.rho
        n = self.rpm/60.0
        J = self.v_inf/(n*D)
        omega = self.rpm*2*pi/60.0
 
        CT = T/(rho*n**2*D**4)
        CQ = Q/(rho*n**2*D**5)
        CP = 2*pi*CQ

        if J==0.0:
            eta = (CT/CP)
        else:
            eta = (CT/CP)*J

        return J, CT, CQ, CP, eta

    def turbine_coeffs(self, T, Q, P):
        rho = self.fluid.rho
        V = self.v_inf
        omega = self.rpm*2*pi/60.0
        TSR = omega*self.rotor.blade_radius/V
        CT = T/(0.5*rho*self.rotor.area*V**2)
        CP = P/(0.5*rho*self.rotor.area*V**3)

        return TSR, CP, CT

        
    def run_sweep(self, parameter, n, low, high):
        """Utility function to run a sweep of a single parameter."""

        if self.mode == 'turbine':
            df = pd.DataFrame(columns = [parameter, 'T', 'Q', 'P', 'TSR', 'CT', 'CP'], index=range(n))
        else:
            if self.coaxial:
                cols = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta', 
                        'CT2', 'CQ2', 'CP2', 'eta2']
            else:
                cols = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta']

            df = pd.DataFrame(columns = cols, index=range(n))

        sections = []
        for i,p in enumerate((np.linspace(low, high, n))):
            setattr(self, parameter, p)
            
            if self.mode == 'turbine':
                T,Q,P,sec_df = self.run()
                TSR,CP,CT = self.turbine_coeffs(T, Q, P)
                df.iloc[i] = [p, T, Q, P, TSR, CT, CP]
            else:
                if self.coaxial:
                    T,Q,P,sec_df,T2,Q2,P2,sec_df2 = self.run()
                    J,CT,CQ,CP,eta = self.rotor_coeffs(T, Q, P)
                    J,CT2,CQ2,CP2,eta = self.rotor_coeffs(T2, Q2, P2)
                    df.iloc[i] = [p, T, Q, P, T2, Q2, P2, J, CT, CQ, CP, eta, CT2, CP2, eta2]
                else:
                    T,Q,P,sec_df = self.run()
                    J,CT,CQ,CP,eta = self.rotor_coeffs(T, Q, P)
                    df.iloc[i] = [p, T, Q, P, J, CT, CQ, CP, eta]

            sections.append(sec_df)
        
        return df, sections
    
    def solve(self, rotor, rpm, v_inflow, r_inflow):
        rotor.precalc()

        omega = rpm*2*pi/60.0
        # Axial momentum (thrust)
        T = 0.0
        # Angular momentum
        Q = 0.0

        for sec in rotor.sections:
            if sec.radius < r_inflow:
                v = v_inflow
            else:
                v = 0.0
                
            if self.solver == 'brute':
                phi = self.brute_solve(sec, v, omega)
            else:
                try:
                    phi = optimize.bisect(sec.func, 0.01*pi, 0.9*pi, args=(v, omega))
                except ValueError as e:
                    print(e)
                    print('Bisect failed, switching to brute solver')
                    phi = self.brute_solve(sec, v, omega)

            
            dT, dQ = sec.forces(phi, v, omega, self.fluid)

            # Integrate
            T += dT
            Q += dQ
        
        # Power
        P = Q*omega  

        return T, Q, P

    def slipstream(self):
        # Slipstream properties

        r_s = self.rotor.blade_radius/sqrt(2.0)
        v_s = self.Cs*sqrt(2*self.T/(self.fluid.rho*self.rotor.area))

        return r_s, v_s

 
    def run(self):
        self.T, self.Q, self.P = self.solve(self.rotor, self.rpm, self.v_inf, self.rotor.diameter)
       
        print('--- Results ---')
        print('Trust (N):\t',self.T)
        print('Torque (Nm):\t',self.Q)
        print('Power (W):\t',self.P)

        # Coaxial calculaction
        if self.coaxial:
            self.r_s, self.v_s = self.slipstream()
           
            self.T2, self.Q2, self.P2 = self.solve(self.rotor2, self.rpm2, self.v_s, self.r_s)

            print('Trust 2 (N):\t',self.T2)
            print('Torque 2 (Nm):\t',self.Q2)
            print('Power 2 (W):\t',self.P2)

            return self.T, self.Q, self.P, self.rotor.sections_dataframe(), self.T2, self.Q2, self.P2, self.rotor2.sections_dataframe()
        
        else:
            return self.T, self.Q, self.P, self.rotor.sections_dataframe()


    def brute_solve(self, sec, v, omega, n=3600):
        """ 
        Solve by a simple brute force procedure, iterating through all
        possible angles and selecting the one with lowest residual.
        """
        resid = np.zeros(n)
        phis = np.linspace(-0.9*np.pi,0.9*np.pi,n)
        for i,phi in enumerate(phis):
            res = sec.func(phi, v, omega)
            if not np.isnan(res):
                resid[i] = res
            else:
                resid[i] = 1e30
        i = np.argmin(abs(resid))
        return phis[i]

    def optimize_pitch(self):
        """
        Optimize rotor pitch for either maximum thrust (propeller) or maximum power (turbine)
        using a genetic evolution algorithm.

        This is intended as an example of how optimization can be done using the scipy.optimize
        package. The overall procedure can be readily modified to optimize for other parameters,
        e.g. a parametrized function for the pitch, or for a parameter sweep instead of 
        a single parameter set.

        return: Array of optimized pitches
        """

        def run_bemt(x):
            print('Current iteration:',x)
            for sec,pitch in zip(self.rotor.sections, x):
                sec.pitch = np.radians(pitch)

            T,Q,P,df = self.run()
            if self.mode == 'turbine':
                return -P
            else:
                return -T
 
        x = [sec.pitch for sec in self.rotor.sections]
        bounds = [(0,30)]*len(x)

        result = optimize.differential_evolution(run_bemt, bounds, tol=1e-1)

        return result 
        
        
