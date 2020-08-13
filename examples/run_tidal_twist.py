"""
Example on how to use the twist variable to modify the pitch angle globally.
"""

import matplotlib.pyplot as pl
from pybemt.solver import Solver

# Run sweep of tip speed ratio with the BEMT method
s = Solver('tidal.ini')
df, section_df = s.run_sweep('twist', 9, -2, 2)

# Plot results
df.plot(x='twist', y='CT', style='C0-o', linewidth=2) 
pl.legend(loc='upper left')
ax = pl.gca()

ax2 = ax.twinx()
df.plot(x='twist', y='CP', style='C1-o', linewidth=2, ax = ax2) 
ax2.axis((-2.5,2.5,0.46,0.48))
pl.legend(loc='upper right')
ax.set_xlabel('Twist angle')
ax.set_ylabel('Thrust coefficient')
ax2.set_ylabel('Power coefficient')
pl.tight_layout()
pl.show()



