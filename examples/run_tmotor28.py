 # -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
import pandas as pd
from math import pi

from pybemt.solver import Solver

s = Solver('tmotor28.ini')

df, sections = s.run_sweep('rpm', 20, 1000.0, 3200.0)

ax = df.plot(x='rpm', y='T') 
ax2 = df.plot(x='rpm', y='P') 

df_exp = pd.read_csv("tmotor28_data.csv", delimiter=';')

df_exp.plot(x='RPM',y='T(N)',style='o',ax=ax)
pl.figure(1)
pl.ylabel('Thrust (N)')
pl.legend(('BEMT','Experiment'))
df_exp.plot(x='RPM',y='P(W)',style='o',ax=ax2)
pl.figure(2)
pl.ylabel('Power (W)')
pl.legend(('BEMT','Experiment'))


pl.show()

