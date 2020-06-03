 # -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import pandas as pd
from math import pi

from pybemt.solver import Solver

s = Solver('tmotor28.ini')

df, sections = s.run_sweep('rpm', 20, 1000.0, 3200.0)
ax = df.plot(x='rpm', y='T') 

p = s.optimize_pitch()

df, sections = s.run_sweep('rpm', 20, 1000.0, 3200.0)
df.plot(x='rpm', y='T', ax=ax) 

pl.xlabel('RPM')
pl.ylabel('Thrust (N)')
pl.legend(('Baseline','Optimized'))


pl.show()

