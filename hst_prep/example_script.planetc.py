import transits
import numpy as np

# Last day of Cycle 23 is 30 Sep 2016:
jd_cycle_end = 2457662.500000

# WASP-121 planet parameters:
syspars = {}
syspars['ld'] = None
syspars['P'] = 1.2749255
syspars['ecc'] = 0.0
syspars['RpRs'] = 0.12454
syspars['aRs'] = 3.754
syspars['T0'] = 2456635.70832
syspars['b'] = 0.160
syspars['incl'] = np.rad2deg( np.arccos( syspars['b']/syspars['aRs'] ) )

# Ephemeris uncertainties:
P_unc = 2.5e-7
T0_unc = 1.1e-4

# Phase constraints:
phase_range = [ 0.85, 0.86 ]

transits.phase_constraints( jd_cycle_end, syspars, P_unc=P_unc, T0_unc=T0_unc, \
                            phase_range=phase_range, n_hstorb=5, nsigma=3 )
