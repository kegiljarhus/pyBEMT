---
title: 'pyBEMT: An implementation of the Blade Element Momentum Theory in Python'
tags:
  - Python
  - BEMT
  - Wind energy
  - Tidal energy
  - Aerodynamics
  - Aerospace
  - Aeronautics
authors:
  - name: Knut Erik T. Giljarhus
    orcid: 0000-0002-4144-0454 
    affiliation: 1
affiliations:
 - name: Department of Mechanical and Structural Engineering and Materials Science, University of Stavanger, Stavanger, Norway
   index: 1
date: 04 June 2020
bibliography: paper.bib
---

# Summary

The use of rotating blades to generate thrust in the form of propellers or
torque in the form of turbines is of great significance for transportation and
energy generation. This has led to extensive research and development of
mathematical models to better understand, predict and optimize the performance
of these machines. The blade element momentum theory (BEMT) is one such method
with a long history dating back to @glauert1935airplane:1935.  

Despite the development of more sophisticated methods, such as vortex
methods[@gohard1978free:1978], the blade element momentum theory is still
considered relevant for the study of rotor design. Its simple formulation lends
itself to use in education and to quickly analyse new ideas. Several
open-source packages include an implementation of BEMT, for instance the
`AeroDyn`[@moriarty2005aerodyn:2005] solver used in the whole-turbine
simulation software `OpenFAST`[@OpenFAST:2020], and the
`QBlade`[@marten2013qblade:2013] software for wind turbine blade design.
`QBlade` has later been forked to develop `JBLADE`[@silvestre2013jblade:2013],
a software focusing on propeller design.

`pyBEMT` is a unified implementation of the blade element moomentum theory,
supporting both propellers and turbines. The software is designed as a
stand-alone Python implementation with emphasis on readability and
extensibility.  Its modular design and permissive license also makes it
suitable for integration into other simulation tools. Other notable features of
the package are a model for coaxial rotors and optimization of rotor parameters
using the differential evolution algorithm in `SciPy`. `pyBEMT` is currently
applied in research projects on rotor design for unmanned aerial vehicles and
turbine design for tidal stream turbines, as well as used in education within
fluid dynamics and computational engineering.  

# Acknowledgements

J\o rgen Apeland, Vetle B. Ingebretsen and Stian R. Hidle from the University
of Stavanger are acknowledged for their contribution of experimental validation
data. 

# References
