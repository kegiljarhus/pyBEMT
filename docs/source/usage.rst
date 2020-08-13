.. highlight:: python

.. _usage:

Usage
=====


Configuration file (.ini)
-------------------------

The parameters to the solver is given in a configuration file
using the `INI <http://en.wikipedia.org/wiki/INI_file>`_ format. 

Each section and available parameters are described below.

[case]
------

- ``v_inf`` --- The inflow velocity in m/s.
- ``rpm`` --- Rotational speed of the rotor in rotations per minute.
- ``twist`` --- Optional global pitch change in degrees.
- ``coaxial`` --- Whether to add a second rotor. Set to False by default.
- ``rpm2`` --- Rotational speed of the second rotor in rotations per minute.
- ``twist2`` --- Optional global pitch change of the second rotor in degrees.
- ``dz`` --- Distance to second rotor.

[rotor]
-------

- ``nblades`` --- Number of blades.
- ``diameter`` --- Diameter in m.
- ``radius_hub`` --- Radius of hub in m.
- ``section`` --- List of airfoil sections given by airfoil name.
- ``radius`` --- List of distance to each section in m.
- ``chord`` --- List of chord length of each section in m.
- ``pitch`` --- List of pitch of each section in degrees.

[rotor2]
--------

Description of second rotor for coaxial simulations, same parameters as 
first rotor.

[fluid]
-------

- ``rho`` --- Fluid density in kg/m^3.
- ``mu`` --- Fluid dynamic viscosity in Pa s.


[solver]
--------

Optional solver settings.

- ``solver`` --- The default solver is a bisection solver. This can be replaced by a more stable brute force solver using 'brute' here.
- ``Cs`` --- Slipstream coefficient for coaxial solver.


Running a single simulation
---------------------------

To run a single simulation, we create a solver object with the 
configuration file as parameter and then execute the :meth:`solver.Solver.run` method:

.. code-block:: python

        from pybemt.solver import Solver

        s = Solver('rotor.ini')
        T, Q, P, section_df = s.run()

The solver returns the thrust [N], the torque [Nm] and the power [W]. Additionally, a pandas DataFrame with the result for each rotor section is provided.


Running a parameter sweep
-------------------------

A typical use case is to run a sweep of a parameter, for instance the 
rotational speed of the rotor. A utility function is provided for this, called
:meth:`solver.Solver.run_sweep`.

This is used in the following way:

.. code-block:: python

        from pybemt.solver import Solver

        s = Solver('rotor.ini')
        df, section_df = s.run_sweep('rpm', 20, 150, 350)

The `run_sweep` function takes as arguments the parameter to sweep, the number 
of points and the minimum and maximum values. It returns a data frame with 
values for each parameter value, as well as a list of dataframes for each
section.


Running an optimization
-----------------------

Optimizaton of parameters can easily be done using the scipy.optimize
library. Currently, only optimization of pitch is supported directly by
the library:

.. code-block:: python

        from pybemt.solver import Solver

        s = Solver('rotor.ini')
 
        pitches = s.optimize_pitch()

The differential evolution algorithm is used in the current implementation,
as it has been found to give the best results.  Each section is considered as
a separate parameter for the optimization. Using a parameterized function instead can
lead to significant speedups, but this is not currently directly supported
by the package.


Adding a new airfoil
--------------------

pyBEMT uses the same file format as the AeroDyn software. However, note 
that currently only the drag and lift tables are used and only a single 
Reynolds number is supported.

To add a new airfoil, either add the file to the airfoils directory and 
re-install the software or add the airfoil file directly to the installed
airfoils directory. The installation location can be found by running
the following code snippet: ::

        import os
        import pybemt
        print(os.path.dirname(pybemt.__file__))

