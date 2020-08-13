# -*- coding: utf-8 -*-
"""
Module for the solver class.
"""

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
    """
    The Solver object loads the config file and contains functions for running a single simulation,
    parameter sweeps and optimization.

    :param string config_path: Path to config file
    """
    def __init__(self, config_path):

        # Read configuration file
        cfg = SafeConfigParser()
        cfg.read(config_path)
        
        # Case
        self.v_inf = cfg.getfloat('case', 'v_inf')
        self.rpm = cfg.getfloat('case', 'rpm')
        if cfg.has_option('case', 'twist'):
            self.twist = cfg.getfloat('case', 'twist')
        else:
            self.twist = 0.0
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
            if cfg.has_option('case', 'twist2'):
                self.twist2 = cfg.getfloat('case', 'twist2')
            else:
                self.twist2 = 0.0
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
        """ 
        Dimensionless coefficients for a rotor. 

        .. math::
            \\text{J} = \\frac{V_\\infty}{nD} \\\\
            C_T = \\frac{T}{\\rho n^2 D^4} \\\\
            C_Q = \\frac{Q}{\\rho n^2 D^5} \\\\
            C_P = 2\\pi C_Q \\\\
            \\eta = \\frac{C_T}{C_P}J \\\\

        :param float T: Thrust
        :param float Q: Torque
        :param float P: Power
        :return: Advance ratio, thrust coefficient, torque coefficient, power coefficient and efficiency
        :rtype: tuple
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
        """
        Dimensionless coefficients for a turbine.

        .. math::
            \\text{TSR} = \\frac{\\Omega R}{V_\\infty} \\\\
            C_T = \\frac{2T}{\\rho A V_\\infty^2} \\\\
            C_P = \\frac{2P}{\\rho A V_\\infty^3} \\\\

        :param float T: Thrust
        :param float Q: Torque
        :param float P: Power
        :return: Tip-speed ratio, power coefficient and thrust coefficient
        :rtype: tuple
        """

        rho = self.fluid.rho
        V = self.v_inf
        omega = self.rpm*2*pi/60.0
        TSR = omega*self.rotor.blade_radius/V
        CT = T/(0.5*rho*self.rotor.area*V**2)
        CP = P/(0.5*rho*self.rotor.area*V**3)

        return TSR, CP, CT

        
    def run_sweep(self, parameter, n, low, high):
        """
        Utility function to run a sweep of a single parameter.

        :param string parameter: Parameter to sweep, must be a member of the Solver class.
        :param int n: Number of runs
        :param float low: Minimum parameter value
        :param float high: Maximum parameter value

        :return: DataFrame of results and list of sections for each run
        :rtype: tuple
        """

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
    
    def solve(self, rotor, twist, rpm, v_inflow, r_inflow):
        """
        Find inflow angle and calculate forces for a single rotor given rotational speed, inflow velocity and radius.

        :param Rotor rotor: Rotor to solve for
        :param float twist: Angle to adjust rotor pitch
        :param float rpm: Rotations per minute
        :param float v_inflow: Inflow velocity
        :param float r_inflow: Inflow radius (equal to blade radius for single rotors)
        :return: Calculated thrust, torque and power for the rotor
        :rtype: tuple
        """

        rotor.precalc(twist)

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
        """
        For coaxial calculations. Calculates slipstream radius and velocity for the upper rotor according to
        momentum theory. Currently only the static case is included.

        .. math::
            r_s = \\frac{R}{\\sqrt{2}} \\\\
            v_s = C_s\\sqrt{\\frac{2 T}{\\rho A}} \\\\

        :return: Radius and velocity of the slipstream
        :rtype: tuple
        """

        r_s = self.rotor.blade_radius/sqrt(2.0)
        v_s = self.Cs*sqrt(2*self.T/(self.fluid.rho*self.rotor.area))

        return r_s, v_s

 
    def run(self):
        """
        Runs the solver, i.e. finds the forces for each rotor.

        :return: Calculated thrust, torque, power and DataFrame with properties for all sections.
        :rtype: tuple
        """
        self.T, self.Q, self.P = self.solve(self.rotor, self.twist, self.rpm, self.v_inf, self.rotor.diameter)
       
        print('--- Results ---')
        print('Trust (N):\t',self.T)
        print('Torque (Nm):\t',self.Q)
        print('Power (W):\t',self.P)

        # Coaxial calculaction
        if self.coaxial:
            self.r_s, self.v_s = self.slipstream()
           
            self.T2, self.Q2, self.P2 = self.solve(self.rotor2, self.twist2, self.rpm2, self.v_s, self.r_s)

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

        :param Section sec: Section to solve for
        :param float v: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :param int n: Number of angles to test for, optional
        :return: Inflow angle with lowest residual
        :rtype: float
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
        
        
