"""
ATCA 7mm continuum reduction pipeline
======================================
Calibrates a multi-target ATCA 7mm (43 GHz) track, extracts point-source
fluxes for every science target via uvmodelfit, writes them to a txt file,
and images only the N brightest sources for a visual check.

Run inside CASA:
    casa -c atca_7mm_pipeline.py
or paste section wise into an interactive CASA session.

"""

import os

vis        = '2026-06-25_1542_C3826/raw.ms' 
bpcal      = '1921-293'      # bandpass calibrator 
fluxcal    = '1934-638'      # flux calibrator
phasecal   = '2355-534'      # phase calibrator
#refant     = 'ca04'          # Using CA04 as reference coz fairly central, 01 is weird w amplitudes, 2 has weird attenuators, 3 out in the field

#flux_results_file = 'target_fluxes_2026-06-25.txt'
#n_brightest_to_image = 2

#List all targets
#obs_summary = listobs(vis=vis)
#targets = sorted(set(
    #f['name'] for k, f in obs_summary.items()
    #if k.startswith('field_') and f['name'] not in (bpcal, fluxcal, phasecal)
#))
#print(f"[INFO] {len(targets)} target fields identified: {targets[:5]} ...")

# Just testing one target
tgt = 'J2301082-593932'
safe = tgt.replace('+', 'p').replace('-', 'm')

# PLOTS - TAKE FOREVER SO SAVE ONCE 
#u-v plot, amplitude vs time plot, *& implies only cross-corr data, amp vs time takes insane memory 
#plotms(vis=vis, xaxis='U', yaxis='V', antenna='*&', coloraxis='field')
#plotms(vis=vis, xaxis='time', yaxis='amp', antenna='*&', coloraxis='baseline')
#plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', coloraxis='baseline')

#Flagging
"""
# save initial flag state
flagmanager(vis=vis, mode='save', versionname='flag_v1')
#Need to flag pointing scans
flagdata(vis=vis, mode='manual', intent='*POINT*', flagbackup=False)
#clip zeros and NaNs
flagdata(vis=vis, mode='clip', clipzeros=True, flagbackup = False) 
#flag antenna shadowing
flagdata(vis=vis, mode='shadow', tolerance = 0.0, flagbackup = False)
# remove first 10 seconds of observations
flagdata(vis=vis, mode='quack', quackinterval=10.0, quackmode='beg', flagbackup=False)
#flags all autocorr
flagdata(vis=vis, autocorr=True)
#Actual RFI flagging - looks good, default options are best for now!
flagdata(vis=vis, mode='tfcrop', extendflags=False, action='apply', display='report', writeflags=True)

#Check in plotted data if everything is well flagged 
#plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', flaggedsymbolshape='circle', customflaggedsymbol=True)
#extend flagged data - if more than 80% of time range is flagged, then flag the entire chunk
flagdata(vis=vis, mode='extend', growtime=80.0, growfreq=80.0, action='apply', display='report', writeflags=True)

#save flags version 
flagmanager(vis=vis, mode='save', versionname='flag_v2')

# calibration

setjy(vis=vis, field=fluxcal, scalebychan=True, standard='Stevens-Reynolds 2016', usescratch=True) #Perley-Butler 2010
# do primary gain calibration - do we need bpcal or fluxcal here - did bpcal in first analysis. now fluxcal
gaincal(vis = vis, caltable = 'cal.G0', field = fluxcal, refant = refant, gaintype = 'G', calmode = 'p', parang = True, solint = '60s') 
#do bandpass response
bandpass(vis=vis, caltable='cal.B0', field=bpcal, spw='', refant=refant, solnorm=True, solint='inf', bandtype='B', gaintable=['cal.G0'], parang=True) 
#secondary gain calibraton??
gaincal(vis=vis, caltable='cal.G1', field=','.join([bpcal,fluxcal,phasecal]), refant=refant, spw='*', gaintype='G', calmode='ap', parang=True, solint='45s', gaintable=['cal.B0'])
#polarization response - do we do phasecal or fluxcal - did phasecal before, now fluxcal
polcal(vis=vis, caltable='cal.D0', field=fluxcal, refant=refant, gaintable=['cal.B0', 'cal.G1'], poltype='Df+QU', solint='inf')

#Need to do better solutions??? so run things again???
bandpass(vis = vis, caltable = 'cal.B1', field = bpcal, spw='', refant = refant, solnorm = True, solint = 'inf', bandtype = 'B', gaintable = ['cal.G1','cal.D0'], parang = True)
gaincal(vis = vis, caltable = 'cal.G2', field=','.join([bpcal,fluxcal,phasecal]), refant = refant, spw = '*', gaintype = 'G', calmode = 'ap', parang = True, solint = '45s', gaintable = ['cal.B1','cal.D0'])
polcal(vis=vis, caltable='cal.D1', field=fluxcal, refant=refant, gaintable=['cal.B1', 'cal.G2'], poltype='Df+QU', solint='inf')

#MAKE BANDPASS PLOTS
#plotms(vis = 'cal.B1', xaxis = 'freq', yaxis = 'amp', coloraxis = 'spw') #plot the bandpass solutions
#plotms(vis='cal.G2', xaxis='time/freq', yaxis='amp/phase', coloraxis='spw', field=bpcal/fluxcal/phasecal, iteraxis='antenna') 

#set the flux
fluxscale(vis=vis, caltable='cal.G2', fluxtable='cal.F0', reference=fluxcal) 

#apply calibration on main cals - default is linear extrapolation - can change using interp. Since D1 uses fluxcal, need fluxcal as second col, and third col is the field for self consistent solns?
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, fluxcal], field=fluxcal, parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, bpcal],   field=bpcal,   parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, phasecal], field=phasecal, parang=True, flagbackup=False)

# apply calibration to all targets

for tgt in targets:
    applycal(vis=vis, field=tgt, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, phasecal], parang=True, flagbackup=False)

#the final set of flagging routines to make sure any outliers are properly flagged - change any options???
flagmanager(vis = vis, mode = 'save', versionname = 'before_rflag')
flagdata(vis=vis, mode='rflag', field=bpcal, spw='', datacolumn='corrected', action='apply', display='report', correlation='ABS_ALL', timedevscale=3.0, freqdevscale=3.0, winsize=3, combinescans=True, ntime='9999999min', extendflags=False, flagbackup=False)
flagdata(vis=vis, mode='rflag', field=fluxcal, spw='', datacolumn='corrected', action='apply', display='report', correlation='ABS_ALL', timedevscale=3.0, freqdevscale=3.0, winsize=3, combinescans=True, ntime='9999999min', extendflags=False, flagbackup=False)
flagdata(vis=vis, mode='rflag', field=phasecal, spw='', datacolumn='corrected', action='apply', display='report', correlation='ABS_ALL', timedevscale=3.0, freqdevscale=3.0, winsize=3, combinescans=True, ntime='9999999min', extendflags=False, flagbackup=False)

for tgt in targets: 
    flagdata(vis=vis, mode='rflag', field=tgt, spw='', datacolumn='corrected', action='apply', display='report', correlation='ABS_ALL', timedevscale=3.0, freqdevscale=3.0, winsize=3, combinescans=True, ntime='9999999min', extendflags=False, flagbackup=False)

#Extend all new flagging
flagdata(vis=vis, mode='extend', field=','.join([bpcal, fluxcal, phasecal]), spw='', action='apply', display='report', flagbackup=False, extendpols=True, correlation='', growtime=95.0, growfreq=95.0, growaround=True, flagneartime=False, flagnearfreq=False, combinescans=True, ntime='9999999min')

for tgt in targets:
    flagdata(vis=vis, mode='extend', field=tgt, spw='', action='apply', display='report', flagbackup=False, extendpols=True, correlation='', growtime=95.0, growfreq=95.0, growaround=True, flagneartime=False, flagnearfreq=False, combinescans=True, ntime='9999999min')

flagmanager(vis=vis, mode='save', versionname='after_rflag')

#Check calibrated flux bp and phase plots - 
#plotms(vis = vis, field = fluxcal, xaxis = 'freq', yaxis = 'amp', correlation ='xx,yy' , ydatacolumn = 'corrected', coloraxis = 'spw')
#plotms(vis=vis, field=bpcal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw')
#plotms(vis=vis, field=phasecal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw')

"""
#Flux extraction via uvmodelfit - ALL targets, ONE combined MS
import os
import sys

# Compute weights with statwt
statwt(vis=vis, datacolumn='corrected')
 
#targets_corrected_ms = 'targets_corrected.ms'
#if os.path.exists(targets_corrected_ms):
#    os.system(f'rm -rf {targets_corrected_ms}')
 
#need to transform for some reason
#mstransform(vis=vis, outputvis=targets_corrected_ms,
#            field=','.join(targets), datacolumn='corrected',
#            keepflags=False)

# generate single MS so it picks up new weights
test_ms = 'target_test_corrected.ms'
if os.path.exists(test_ms):
    os.system(f'rm -rf {test_ms}')

mstransform(vis=vis, outputvis=test_ms, field=tgt, datacolumn='corrected', keepflags=False)

#function to capture the terminal output because what it outputs to a normal txt file was bs
 
def capture_c_level_stdout(logfile, func, **kwargs):
    stdout_fd = sys.stdout.fileno()
    saved_fd = os.dup(stdout_fd)
    with open(logfile, 'w') as f:
        sys.stdout.flush()
        os.dup2(f.fileno(), stdout_fd)
        try:
            result = func(**kwargs)
        finally:
            sys.stdout.flush()
            os.dup2(saved_fd, stdout_fd)
            os.close(saved_fd)
    return result
 
#results = []  # (field, flux_Jy, flux_err_Jy, status)
 
#for tgt in targets[:5]: #Only on first 5 targets for now
#    safe = tgt.replace('+', 'p').replace('-', 'm') #some weird cl business
#    clfile = f'{safe}.cl'
#    logfile = f'{safe}_uvmodelfit.log'

# Save single file now instead of all
clfile = f'{safe}_statwt_test.cl'
logfile = f'{safe}_statwt_test_uvmodelfit.log'
 
if os.path.exists(clfile):
   os.system(f'rm -rf {clfile}')
 
 # Modified to use single ms for this test instead
try:
    capture_c_level_stdout(
        logfile, uvmodelfit, #use uvmodelfit to fit for a ptsrc flux
        vis=test_ms,
        field=tgt,                # select just this target within the combined MS
        comptype='P',              # point source
        sourcepar=[0.005, 0.0, 0.0],  # [flux Jy, RA offset, Dec offset] initial guess
        niter=10,
        spw='',                    # fit across all combined SPWs
        outfile=clfile
    )
        
    print(f"[DONE] Generated new log file for testing.")
   
   
except Exception as e:
    print("Error in uvmodelfit")

