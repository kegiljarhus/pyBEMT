 # -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import sys

from pybemt.solver import Solver

s = Solver('tmotor28_coaxial.ini')

dfe = pd.read_csv("tmotor28_coaxial_data.csv", delimiter=',')

Ts=[]
T2s=[]
Ps=[]
P2s=[]
for rpm_lower, rpm_upper in zip(dfe['RPM'], dfe['RPM_B']):
    s.rpm = rpm_upper
    s.rpm2 = rpm_lower
    T,Q,P,dfU,T2,Q2,P2,dfL = s.run()

    Ts.append(T)
    T2s.append(T2)
    Ps.append(P)
    P2s.append(P2)

pl.figure()
ax=pl.gca()
rpm = dfe['RPM_B']
pl.plot(rpm, Ts) 
pl.plot(rpm, T2s) 
pl.figure()
ax2=pl.gca()
pl.plot(rpm, Ps) 
pl.plot(rpm, P2s) 

df28 = pd.read_csv("tmotor28_data.csv", delimiter=';')

dfe.plot(x='RPM_B',y='T_B(N)',style='C0o',ax=ax)
dfe.plot(x='RPM_B',y='T_A(N)',style='C1o',ax=ax)
df28.plot(x='RPM',y='T(N)',style='--',ax=ax,color='lightgrey')
pl.figure(1)
pl.ylabel('Thrust (N)')
pl.legend(('Upper BEMT','Lower BEMT', 'Upper Experiment', 'Lower Experiment', '28 isolated'))

dfe.plot(x='RPM_B',y='P_B(W)',style='C0o',ax=ax2)
dfe.plot(x='RPM_B',y='P_A(W)',style='C1o',ax=ax2)
df28.plot(x='RPM',y='P(W)',style='--',ax=ax2,color='lightgrey')
pl.figure(2)
pl.ylabel('Power (W)')
pl.legend(('Upper BEMT','Lower BEMT', 'Upper Experiment', 'Lower Experiment'))
pl.legend(('Upper BEMT','Lower BEMT', 'Upper Experiment', 'Lower Experiment', '28 isolated'))

pl.show()


