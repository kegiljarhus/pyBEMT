import pandas as pd
from configparser import NoOptionError
from math import radians, degrees, sqrt, cos, sin, atan2, atan, pi, acos, exp
from .airfoil import load_airfoil

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
        
        self.alpha = [float(p) for p in cfg.get(name, 'pitch').split()]
        self.sections = []
        for i in range(self.n_sections): 
            sec = Section(load_airfoil(s[i]), float(r[i]), float(dr[i]), radians(self.alpha[i]), float(c[i]), self, mode)
            self.sections.append(sec)
        
        self.radius_hub = cfg.getfloat(name,'radius_hub')

        self.precalc()

    def precalc(self):
        self.blade_radius = 0.5*self.diameter
        self.area = pi*self.blade_radius**2

    def sections_dataframe(self):
        columns = ['radius','chord','pitch','Cl','Cd','dT','dQ','F','a','ap','Re']
        data = {}
        for param in columns:
            array = [getattr(sec, param) for sec in self.sections]
            data[param] = array
        
        return pd.DataFrame(data)
 

class Section: 
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
        
    

