"""
Run analysis of power output and drag of a three-bladed tidal turbine.

The results are compared against cavitation tunnel tests presented in
    Bahaj, A. S., et al. "Power and thrust measurements of marine current 
    turbines under various hydrodynamic flow conditions in a cavitation 
    tunnel and a towing tank." Renewable energy 32.3 (2007): 407-426.
"""

import matplotlib.pyplot as pl
from pybemt.solver import Solver

# Run sweep of tip speed ratio with the BEMT method
s = Solver('tidal.ini')
df, section_df = s.run_sweep('rpm', 20, 150.0, 330.0)

# Experimental data 
cp_tsr = [4.170616, 4.423381,4.660348, 4.897314,5.134281,5.371248, 
        5.371248,5.592417, 5.845182, 6.082148,6.303318,6.540284, 6.777251,
        7.014218,7.219589, 7.440758, 7.693523]
cp = [0.413793, 0.430885, 0.437181, 0.446177,0.445277, 0.454273,
        0.457871,0.452474,0.454273,0.452474,0.449775,0.447976,0.441679,
        0.435382, 0.429085, 0.412894, 0.409295]
ct_tsr = [4.184953, 4.435737,4.670846, 4.905956, 5.141066,5.39185, 
        5.376176,5.611285, 5.862069,5.877743,6.097179,6.097179,6.316614,
        6.551724,6.786834,7.021944,7.257053, 7.460815,7.711599]
ct = [0.64451, 0.672997, 0.701484,0.726409, 0.740653, 0.753116,
        0.77270, 0.778042,0.795846,0.802967, 0.811869,0.817211,0.820772,
        0.836795, 0.851039,0.859941, 0.874184, 0.884866,0.890208]

# Plot results
pl.plot(ct_tsr, ct, 'C0o')
ax = pl.gca()
df.plot(x='TSR', y='CT', style='C0-', linewidth=2, ax=ax) 

df.plot(x='TSR', y='CP', style='C1-', linewidth=2, ax=ax) 
pl.plot(cp_tsr, cp,'C1o')
pl.legend(('Exp CT', 'BEMT CT', 'Exp CP', 'BEMT CP'))
pl.xlabel('Blade tip speed ratio')
pl.ylabel('Power and thrust coefficients')

pl.show()



