# -*- coding: utf-8 -*-

"""
Module for storing rotor properties and calculation of induction factors and forces for the airfoil sections.
"""
import pandas as pd
from configparser import NoOptionError
from math import radians, degrees, sqrt, cos, sin, atan2, atan, pi, acos, exp
from .airfoil import load_airfoil

class Rotor: 
    """
    Holds rotor properties and a list of all airfoil sections.

    :param configparser.SafeConfigParser cfg: Configuration object
    :param string name: Name of rotor
    :param string mode: Solver mode
    """
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
        
        self.alpha = [float(p) for p in cfg.get(name, 'pitch').split()]
        self.sections = []
        for i in range(self.n_sections): 
            sec = Section(load_airfoil(s[i]), float(r[i]), float(dr[i]), radians(self.alpha[i]), float(c[i]), self, mode)
            self.sections.append(sec)
        
        self.radius_hub = cfg.getfloat(name,'radius_hub')

        self.precalc(twist=0.0)

    def precalc(self, twist):
        """
        Calculation of properties before each solver run, to ensure all parameters are correct for parameter sweeps.

        :return: None
        """
        self.blade_radius = 0.5*self.diameter
        self.area = pi*self.blade_radius**2

        # Apply twist
        for i,sec in enumerate(self.sections):
            sec.pitch = radians(self.alpha[i] + twist)


    def sections_dataframe(self):
        """
        Creates a pandas DataFrame with all calculated section properties.

        :return: DataFrame with section properties
        :rtype: pd.DataFrame
        """

        columns = ['radius','chord','pitch','Cl','Cd','dT','dQ','F','a','ap','Re','AoA']
        data = {}
        for param in columns:
            array = [getattr(sec, param) for sec in self.sections]
            data[param] = array
        
        return pd.DataFrame(data)
 

class Section: 
    """
    Class for calculating induction factors and forces according to the BEM theory for a single airfoil section.

    :param Airfoil airfoil: Airfoil of the section
    :param float radius: Distance from center to mid of section
    :param float width: Width of section
    :param float pitch: Pitch angle in radians
    :param float chord: Chord length of section
    :param Rotor rotor: Rotor that section belongs to
    :param string mode: Solver mode
    """
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
        self.AoA = 0.0
        self.dT = 0.0
        self.dQ = 0.0
        self.F = 0.0
        self.Cl = 0.0
        self.Cd = 0.0

        self.precalc()
        
    def precalc(self):
        """
        Calculation of properties before each solver run, to ensure all parameters are correct for parameter sweeps.

        :return: None
        """
        self.sigma = self.rotor.n_blades*self.chord/(2*pi*self.radius)

    def tip_loss(self, phi):
        """
        Prandtl tip loss factor, defined as

        .. math::
            F = \\frac{2}{\\pi}\\cos^{-1}e^{-f} \\\\
            f = \\frac{B}{2}\\frac{R-r}{r\\sin\\phi}

        A hub loss is also caluclated in the same manner.

        :param float phi: Inflow angle
        :return: Combined tip and hub loss factor
        :rtype: float
        """
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
    
    def stall_delay_model(self, phi, alpha, Cl, Cd):
        """
        The 3D correction model based on Chaviaropoulos and Hansen ref:
        
        .. bib::
            @article{chaviaropoulos2000investigating,
                     title={Investigating three-dimensional and rotational effects on wind turbine blades by means of a quasi-3D Navier-Stokes solver},
                     author={Chaviaropoulos, PK and Hansen, Martin OL},
                     journal={J. Fluids Eng.},
                     volume={122},
                     number={2},
                     pages={330--336},
                     year={2000}
                     }

        .. math::
            Cl_3D = Cl_2D + a (c / r)^h \cos^n{twist} (Cl_inv - Cl_2d) \\
            Cl_inv = \sqrt{Cl_2d^2 + Cd^2}
        where:
        a = 2.2, h = 1.3 and n = 4

        :param float phi: Inflow angle
        :return: Lift coefficient with 3D correction
        :rtype: float
        """
        Cl_inv = sqrt(Cl**2+Cd**2)
        twist = alpha - self.pitch
        r = self.radius - self.rotor.radius_hub
        c = self.chord
        
        a = 2.2
        h = 1.3
        n = 4
        return Cl + a * (c / r) ** h * (cos(twist)) ** n * (Cl_inv - Cl)
                
    def airfoil_forces(self, phi):
        """
        Force coefficients on an airfoil, decomposed in axial and tangential directions:

        .. math::
            C_T = C_l\\cos{\\phi} - CC_d\\sin{\\phi} \\\\
            C_Q = C_l\\sin{\\phi} + CC_d\\cos{\\phi} \\\\

        where drag and lift coefficients come from
        airfoil tables.

        :param float phi: Inflow angle
        :return: Axial and tangential force coefficients
        :rtype: tuple
        """

        C = self.C

        alpha = C*(self.pitch - phi)
                
        Cl = self.airfoil.Cl(alpha)
        Cd = self.airfoil.Cd(alpha)
                
        CT = Cl*cos(phi) - C*Cd*sin(phi)
        CQ = Cl*sin(phi) + C*Cd*cos(phi)

        self.AoA = degrees(alpha)
        self.Cl = float(Cl)
        self.Cd = float(Cd)
        
        return CT, CQ
    
    def induction_factors(self, phi):
        """
        Calculation of axial and tangential induction factors,

        .. math::
            a = \\frac{1}{\\kappa - C} \\\\
            a\' = \\frac{1}{\\kappa\' + C} \\\\
            \\kappa = \\frac{4F\\sin^2{\\phi}}{\\sigma C_T} \\\\
            \\kappa\' = \\frac{4F\\sin{\\phi}\\cos{\\phi}}{\\sigma C_Q} \\\\
            
        :param float phi: Inflow angle
        :return: Axial and tangential induction factors
        :rtype: tuple
        """

        C = self.C
        
        F = self.tip_loss(phi)
        
        CT, CQ = self.airfoil_forces(phi)
        
        kappa = 4*F*sin(phi)**2/(self.sigma*CT)
        kappap = 4*F*sin(phi)*cos(phi)/(self.sigma*CQ)

        a = 1.0/(kappa - C)
        ap = 1.0/(kappap + C)
        
        return a, ap
        
    def func(self, phi, v_inf, omega):
        """
        Residual function used in root-finding functions to find the inflow angle for the current section.

        .. math::
            \\frac{\\sin\\phi}{1+Ca} - \\frac{V_\\infty\\cos\\phi}{\\Omega R (1 - Ca\')} = 0\\\\

        :param float phi: Estimated inflow angle
        :param float v_inf: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :return: Residual
        :rtype: float
        """
        # Function to solve for a single blade element
        C = self.C

        a, ap = self.induction_factors(phi)
        
        resid = sin(phi)/(1 + C*a) - v_inf*cos(phi)/(omega*self.radius*(1 - C*ap))
        
        self.a = a
        self.ap = ap
        
        return resid
    
    def forces(self, phi, v_inf, omega, fluid):
        """
        Calculation of axial and tangential forces (thrust and torque) on airfoil section. 

        The definition of blade element theory is used,

        .. math::
            \\Delta T = \\sigma\\pi\\rho U^2C_T r\\Delta r \\\\
            \\Delta Q = \\sigma\\pi\\rho U^2C_Q r^2\\Delta r \\\\
            U = \\sqrt{v^2+v\'^2} \\\\
            v = (1 + Ca)V_\\infty \\\\
            v\' = (1 - Ca\')\\Omega R \\\\

        Note that this is equivalent to the momentum theory definition,

        .. math::
            \\Delta T = 4\\pi\\rho r V_\\infty^2(1 + Ca)aF\\Delta r \\\\
            \\Delta Q = 4\\pi\\rho r^3 V_\\infty\\Omega(1 + Ca)a\'F\\Delta r \\\\


        :param float phi: Inflow angle
        :param float v_inf: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :param Fluid fluid: Fluid 
        :return: Axial and tangential forces
        :rtype: tuple
        """

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
        # dQ = 4*pi*rho*r**3*self.v_inf*(1 + a)*ap*self.omega*F
                
        return self.dT, self.dQ
        
    

