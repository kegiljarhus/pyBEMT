# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
from configparser import SafeConfigParser, NoOptionError
from math import radians, degrees, sqrt, cos, sin, atan2, atan, pi, acos, exp
from scipy.optimize import bisect
from .airfoil import load_airfoil


class Section: # TODO - is this stupid? Just have arrays instead?
    def __init__(self, airfoil, radius, width, pitch, chord, rotor, mode):
        self.airfoil = airfoil
        self.radius = radius
        self.width = width
        self.pitch = pitch
        self.chord = chord
        self.rotor = rotor

        if mode == 'turbine':
            self.C = -1
        else:
            self.C = 1
        
        self.v = 0.0
        self.v_theta = 0.0
        self.v_rel = 0.0
        self.a=0.0
        self.ap=0.0
        self.Re = 0.0
        self.alpha = 0.0
        self.dT = 0.0
        self.dQ = 0.0
        self.F = 0.0
        self.Cl = 0.0
        self.Cd = 0.0

        self.precalc()
        
    def precalc(self):
        self.sigma = self.rotor.n_blades*self.chord/(2*pi*self.radius)

    def tip_loss(self, phi):
        def prandtl(dr, r, phi):
            f = self.rotor.n_blades*dr/(2*r*(sin(phi)))
            if (-f > 500): # exp can overflow for very large numbers
                F = 1.0
            else:
                F = 2*acos(min(1.0, exp(-f)))/pi
                
            return F
        
        if phi == 0:
            F = 1.0
        else:    
            r = self.radius
            Ftip = prandtl(self.rotor.blade_radius - r, r, phi)
            Fhub = prandtl(r - self.rotor.radius_hub, r, phi)
            F = Ftip*Fhub
            
        self.F = F
        return F
 
                    
    def airfoil_forces(self, phi):
        C = self.C

        alpha = C*(self.pitch - phi)
                
        Cl = self.airfoil.Cl(alpha)
        Cd = self.airfoil.Cd(alpha)
                
        CT = Cl*cos(phi) - C*Cd*sin(phi)
        CQ = Cl*sin(phi) + C*Cd*cos(phi)
        
        return CT, CQ
    
    def induction_factors(self, phi):
        C = self.C
        
        F = self.tip_loss(phi)
        
        CT, CQ = self.airfoil_forces(phi)
        
        kappa = 4*F*sin(phi)**2/(self.sigma*CT)
        kappap = 4*F*sin(phi)*cos(phi)/(self.sigma*CQ)

        a = 1.0/(kappa - C)
        ap = 1.0/(kappap + C)
        
        return a, ap
        
    def func(self, phi, v_inf, omega):
        # Function to solve for a single blade element
        C = self.C

        a, ap = self.induction_factors(phi)
        
        resid = sin(phi)/(1 + C*a) - v_inf*cos(phi)/(omega*self.radius*(1 - C*ap))
        
        self.a = a
        self.ap = ap
        
        return resid
    
    def forces(self, phi, v_inf, omega, fluid):
        C = self.C
        r = self.radius
        rho = fluid.rho
        
        a, ap = self.induction_factors(phi)
        CT, CQ = self.airfoil_forces(phi)
        
        v = (1 + C*a)*v_inf
        vp = (1 - C*ap)*omega*r
        U = sqrt(v**2 + vp**2)   
        
        self.Re = rho*U*self.chord/fluid.mu
            
        # From blade element theory
        self.dT = self.sigma*pi*rho*U**2*CT*r*self.width
        self.dQ = self.sigma*pi*rho*U**2*CQ*r**2*self.width

        # From momentum theory
        # dT = 4*pi*rho*r*self.v_inf**2*(1 + a)*a*F
        # dQ = 4*pi*rho*r**3*self.v_inf*(1 + a)*a*self.omega*F
                
        return self.dT, self.dQ
        
    
class Blade: # struct of arrays instead of array of structs
    pass

class Rotor: # struct of arrays instead of array of structs
    def __init__(self, cfg, name, mode):
        self.n_blades = cfg.getint(name, 'nblades')
        self.diameter = cfg.getfloat(name, 'diameter')

        s = cfg.get(name, 'section').split()
        c = cfg.get(name, 'chord').split()
        r = cfg.get(name, 'radius').split()
        self.n_sections = len(s)
        try:
            dr = cfg.get(name,'dr').split()
        except NoOptionError:
            dr = self.n_sections*[float(r[1]) - float(r[0])]
        
        p   = cfg.get(name, 'pitch').split()
        self.sections = []
        for i in range(self.n_sections): 
            sec = Section(load_airfoil(s[i]), float(r[i]), float(dr[i]), radians(float(p[i])), float(c[i]), self, mode)
            self.sections.append(sec)
        
        self.radius_hub = cfg.getfloat(name,'radius_hub')

        self.precalc()

    def precalc(self):
        self.blade_radius = 0.5*self.diameter

    def sections_dataframe(self):
        columns = ['radius','chord','pitch','Cl','Cd','dT','dQ','F','a','ap','Re']
        data = {}
        for param in columns:
            array = [getattr(sec, param) for sec in self.sections]
            data[param] = array
        
        return pd.DataFrame(data)
 
class Fluid:
    def __init__(self, cfg):
        self.rho = cfg.getfloat('fluid','rho')
        self.mu = cfg.getfloat('fluid','mu')


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
            self.rotor2 = Rotor(cfg, 'rotor2')
            self.zD = cfg.getfloat('case','dz')/self.rotor.diameter
            self.T2 = 0
            self.Q2 = 0
            self.P2 = 0

        # Solver
        self.solver = 'bisect'
       
    def rotor_coeffs(self, T, Q, P):
        """ Calculate non-dimensional coefficients. For propellers,
        we typically use the advance ratio,
        to 
        """
        D = self.rotor.diameter
        rho = self.fluid.rho
        n = self.rpm/60.0
        J = self.v_inf/(n*D)
 
        CT = T/(rho*n**2*D**4)
        CQ = Q/(rho*n**2*D**5)
        CP = 2*pi*CQ

        if J==0.0:
            eta = (CT/CP)
        else:
            eta = (CT/CP)*J

        # TODO - coaxial???
        
        return J, CT, CQ, CP, eta

    def turbine_coeffs(self, T, Q, P):
        R = 0.5*self.rotor.diameter
        A = pi*R**2
        rho = self.fluid.rho
        V = self.v_inf
        omega = self.rpm*2*pi/60.0
        TSR = omega*R/V
        CT = T/(0.5*rho*A*V**2)
        CP = P/(0.5*rho*A*V**3)

        return TSR, CP, CT
        
    def run_sweep(self, parameter, n, low, high):
        """Utility function to run a sweep of a single parameter."""

        if self.mode == 'turbine':
            df = pd.DataFrame(columns = [parameter, 'T', 'Q', 'P', 'TSR', 'CT', 'CP'], index=range(n))
        else:
            df = pd.DataFrame(columns = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta'], index=range(n))
            if self.coaxial:
                df2 = pd.DataFrame(columns = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta'], index=range(n))

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
                    df.iloc[i] = [p, T, Q, P, T2, Q2, P2, J, CT, CQ, CP, eta]
                else:
                    T,Q,P,sec_df = self.run()
                    J,CT,CQ,CP,eta = self.rotor_coeffs(T, Q, P)
                    df.iloc[i] = [p, T, Q, P, J, CT, CQ, CP, eta]

            sections.append(sec_df)
        
        return df, sections
    
    def solve(self, rotor, rpm, v_inflow, r_inflow):

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
                phi = bisect(sec.func, 0.01*pi, 0.9*pi, args=(v, omega))
            
            dT, dQ = sec.forces(phi, v, omega, self.fluid)

            # Integrate
            T += dT
            Q += dQ
        
        # Power
        P = Q*omega  

        return T, Q, P
 
    def run(self):

        R = 0.5*self.rotor.diameter
        self.T, self.Q, self.P = self.solve(self.rotor, self.rpm, self.v_inf, R)
       
        print('--- Results ---')
        print('Trust (N):\t',self.T)
        print('Torque (Nm):\t',self.Q)
        print('Power (W):\t',self.P)

        # Coaxial calculaction
        if self.coaxial:
            # Slipstream properties
            Cs = 1.2
            Cs = 1.0
            # zD 32: 0.14
            # zD 28: 0.16
            self.r_s = Cs*R/sqrt(2.0)
            #A = pi*(self.r_s**2 - self.rotor.radius_hub**2)
            A = pi*R**2
            #Cs = 5.0 # 0.75 bra for 2828, 1.1 for 3232
            #Cs = 3.0
            Cs= 0.5 # awesome for 28/32
            Cs= 0.75 # awesome for 28/28
            Cs=0.625
            # velocity from momentum theory
            self.v_s = Cs*sqrt(2*self.T/(self.rho*A))
            print('Coaxial rad/vel',self.r_s, self.v_s)
            
            self.T2, self.Q2, self.P2 = self.solve(self.rotor2, self.rpm2, self.v_s, self.r_s)

            print('Trust 2 (N):\t',self.T2)
            print('Torque 2 (Nm):\t',self.Q2)
            print('Power 2 (W):\t',self.P2)

        
        if self.coaxial:
            return self.T, self.Q, self.P, self.rotor.sections_dataframe(), self.T2, self.Q2, self.P2, self.rotor2.sections_dataframe()
        else:
            return self.T, self.Q, self.P, self.rotor.sections_dataframe()


    def brute_solve(self, sec, v, omega, n=3600):
        """ Solve by a simple brute force procedure, iterating through all
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

 
       
        
        
        
