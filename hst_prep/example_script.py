import transits
import batman
import numpy as np

# Last day of Cycle 23 is 30 Sep 2016:
jd_cycle_end = 2457662.500000

# WASP-121 planet parameters:
batpar = batman.TransitParams()
batpar.per = 1.2749255
batpar.ecc = 0.0
batpar.w = 90
batpar.rp = 0.12454
batpar.a = 3.754
batpar.t0 = 2456635.70832
b = 0.160
batpar.inc = np.rad2deg( np.arccos( b/batpar.a ) )

# Ephemeris uncertainties:
P_unc = 2.5e-7
T0_unc = 1.1e-4

# Phase constraints:
phase_range = [ 0.85, 0.86 ]

transits.phase_constraints( jd_cycle_end, batpar, P_unc=P_unc, T0_unc=T0_unc, \
                            phase_range=phase_range, n_hstorb=5, nsigma=3 )
