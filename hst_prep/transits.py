import numpy as np
import matplotlib.pyplot as plt
import pdb, sys, os
#from planetc import transit
import batman

HST_ORB_PERIOD_DAYS = 96.36/60./24.

def phase_constraints( jd_cycle_end, batpar, P_unc=0.0, T0_unc=0.0, phase_range=[], phase_exposure=['first',1], n_hstorb=5, tvis_minutes=45, nsigma=3 ):
    """
    Given a particular transit system and a set of orbital phase
    constraints, plot the resulting HST phase coverage over the
    specified number of orbits.

    NOTE: This is currently hard-wired to primary transits only. 
    Wouldn't be too hard to generalise to primary transits and secondary eclipses.

    """

    # Find the last transit of the HST cycle:
    Tmid = batpar.t0
    while Tmid<jd_cycle_end:
        Tmid += batpar.per
    while Tmid>=jd_cycle_end:
        Tmid -= batpar.per

    # Get the number of planet orbits between the 
    # reference epoch and the HST cycle:
    norb = int( np.round( np.abs( Tmid - batpar.t0 )/batpar.per ) )
    
    # Determine the uncertainty on the mid-time by
    # the last transit of the HST cycle:
    unc = T0_unc + norb*P_unc

    # Determine the JD corresponding to the earliest possible
    # observation start time, i.e. the lower phase constraint
    # and assuming the lower limit of the possible Tmid value:
    phase_low = phase_range[0] - nsigma*unc/batpar.per
    phase_upp = phase_range[1] + nsigma*unc/batpar.per
    jd_phase_low = Tmid - batpar.per*( 1-phase_low )
    phase_orb = phase_exposure[1]
    dt = ( phase_orb-1 )*HST_ORB_PERIOD_DAYS
    if phase_exposure[0]=='first':
        jd_start = jd_phase_low-dt
    elif phase_exposure[0]=='last':
        jd_start = ( jd_phase_low-tvis_minutes/60./24. )-dt
        
    # Determine the offset between the earliest possible start
    # time and the latest possible start time:
    delt = batpar.per*( phase_range[1]-phase_range[0] ) + 2*nsigma*unc

    # Amount of time spent observing
    # per HST orbit:
    #tvis = 50./60./24. # i.e. 50 minutes
    tvis = tvis_minutes/60./24. # i.e. 50 minutes

    # Make sure the limb darkening is turned off:
    batpar.limb_dark = 'quadratic'
    batpar.u = [ 0, 0 ]

    # Set an arbitrary number of points to be evaluated
    # per HST orbit, at sufficiently high time-sampling to
    # trace out the shape of the transit:
    npoints = 100
    jd = []
    for i in range( n_hstorb ):
        jd_i = np.r_[ jd_start:jd_start+tvis:1j*npoints ]
        jd += [ jd_i ]
        jd_start += HST_ORB_PERIOD_DAYS
    jd1 = np.concatenate( jd )
    t1 = 1+( jd1-Tmid )/batpar.per
    #psignal1 = transit.ma02_aRs( jd1, **syspars )
    jd2 = jd1+delt
    t2 = 1+( jd2-Tmid )/batpar.per
    #psignal2 = transit.ma02_aRs( jd2, **syspars )
    transittype = 'primary' # currently hard-wired for primary transits only...
    pmodel1 = batman.TransitModel( batpar, jd1, transittype=transittype )
    pmodel2 = batman.TransitModel( batpar, jd2, transittype=transittype )
    jdfull = np.r_[ jd1.min():jd2.max():1j*1000 ]
    tfull = 1+( jdfull-Tmid )/batpar.per
    #psignalfull = transit.ma02_aRs( jdfull, **syspars )
    pmodelfull = batman.TransitModel( batpar, jdfull, transittype=transittype )
    psignal1 = pmodel1.light_curve( batpar )
    psignal2 = pmodel2.light_curve( batpar )
    psignalfull = pmodelfull.light_curve( batpar )
    
    ixs = ( psignalfull<psignalfull.max()-0.01*( psignalfull.max()-psignalfull.min() ) )
    tdur = jdfull[ixs].max()-jdfull[ixs].min()
    print( 'Transit duration = {0:.2f} days = {1:.0f} minutes'.format( tdur, tdur*24*60 ) )

    # Make plot:
    fig = plt.figure( figsize=[14,6] )
    c1 = 'c'
    c2 = 'r'
    plt.plot( tfull, psignalfull, '-k', lw=2, zorder=0 )
    ms = 8
    plt.plot( t1, psignal1, 'o', mfc=c1, mec='none', ms=ms, zorder=1 )
    plt.plot( t2, psignal2, 'o', mfc=c2, mec='none', ms=ms*0.5, zorder=2 )
    dy = psignal1.max()-psignal1.min()
    ymax = psignal1.max()+0.1*dy
    ymin = psignal1.min()-0.1*dy
    plt.ylim( [ ymin, ymax ] )
    xmin = min( [ t1.min(), t2.min() ] )
    xmax = max( [ t1.max(), t2.max() ] )
    dx = xmax-xmin
    plt.xlim( [ xmin-0.1*dx, xmax+0.1*dx ] )
    plt.ylabel( 'Relative flux' )
    plt.xlabel( 'Planet orbital phase' )
    text_fs = 14
    text_str = 'Assuming:\n  HST period = {0:.2f} minutes\n  Visibility = {1:.2f} minutes'\
               .format( HST_ORB_PERIOD_DAYS*24*60, tvis*24*60 )
    #text_str = '{0}\n  P = {1:.8f} days\n  T0 = {2:.8f}'.format( text_str, syspars['P'], syspars['T0'] )
    text_str = '{0}\n  P = {1:.8f} days\n  T0 = {2:.8f}'.format( text_str, batpar.per, batpar.t0 )
    expstr = 'to {0} exposure of orbit {1}:'.format( phase_exposure[0], phase_exposure[1] )
    if phase_range[0]<0: phase_range[0] += 1
    if phase_range[1]<0: phase_range[1] += 1
    text_str = '{0}\n\nSpecified phase range applied\n{1}\n  Lower = {2:.5f}\n  Upper = {3:.5f}'\
               .format( text_str, expstr, phase_range[0], phase_range[1] )
    nmin_window = int( np.ceil( 24*60*batpar.per*( phase_range[1]-phase_range[0] ) ) )
    text_str = ' {0}\n({1} minute window)'.format( text_str, nmin_window )
    text_str = '{0}\n\n{1}-sigma plausible range:\n  Phase lower = {2:.5f}\n  Phase upper = {3:.5f}'\
               .format( text_str, nsigma, phase_low, phase_upp )
    ax = plt.gca()
    ax.text( 0.03, 0.8, text_str, fontsize=text_fs, transform=ax.transAxes, \
             horizontalalignment='left', verticalalignment='top' )

    return fig
