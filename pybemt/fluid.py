# -*- coding: utf-8 -*-

"""
Module for holding and calculating fluid properties. Currently, only viscosity and density are included.
"""

class Fluid:
    """
    Class for loading fluid properties from configuration file and providing them to the solver.

    :param configparser.SafeConfigParser cfg: Configuration object
    """

    def __init__(self, cfg):
        self.rho = cfg.getfloat('fluid','rho')
        self.mu = cfg.getfloat('fluid','mu')


