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
refant     = 'ca04'          # Using CA04 as reference coz fairly central, 01 is weird w amplitudes, 2 has weird attenuators, 3 out in the field

flux_results_file = 'target_fluxes_2026-06-25.txt'
n_brightest_to_image = 2

#List all targets
obs_summary = listobs(vis=vis)
targets = sorted(set(
    f['name'] for k, f in obs_summary.items()
    if k.startswith('field_') and f['name'] not in (bpcal, fluxcal, phasecal)
))
print(f"[INFO] {len(targets)} target fields identified: {targets[:5]} ...")

# PLOTS - TAKE FOREVER SO SAVE ONCE 
#u-v plot, amplitude vs time plot, *& implies only cross-corr data, amp vs time takes insane memory 
#plotms(vis=vis, xaxis='U', yaxis='V', antenna='*&', coloraxis='field')
#plotms(vis=vis, xaxis='time', yaxis='amp', antenna='*&', coloraxis='baseline')
#plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', coloraxis='baseline')

#Flagging

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
#Actual RFI flagging - unsure about what extra options to use here
flagdata(vis=vis, mode='tfcrop', extendflags=False, action='apply', display='report', writeflags=True)

#Check in plotted data if everything is well flagged 
plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', flaggedsymbolshape='circle', customflaggedsymbol=True)
#extend flagged data
flagdata(vis=vis, mode='extend', growtime=80.0, growfreq=80.0, action='apply', display='report', writeflags=True)

#save flags version 
flagmanager(vis=vis, mode='save', versionname='flag_v2')

# calibration

setjy(vis=vis, field=fluxcal, scalebychan=True, standard='Perley-Butler 2010', usescratch=True)
# do primary gain calibration - do we need bpcal or fluxcal here???
gaincal(vis = vis, caltable = 'cal.G0', field = bpcal, refant = refant, gaintype = 'G', calmode = 'p', parang = True, solint = '60s') 
#do bandpass response
bandpass(vis=vis, caltable='cal.B0', field=bpcal, spw='', refant=refant, solnorm=True, solint='inf', bandtype='B', gaintable=['cal.G0'], parang=True) 
#secondary gain calibraton??
gaincal(vis=vis, caltable='cal.G1', field=','.join([bpcal,fluxcal,phasecal]), refant=refant, spw='*', gaintype='G', calmode='ap', parang=True, solint='45s', gaintable=['cal.B0'])
#polarization response - do we do phasecal or fluxcal?
polcal(vis=vis, caltable='cal.D0', field=phasecal, refant=refant, gaintable=['cal.B0', 'cal.G1'], poltype='Df+QU', solint='inf')

#Need to do better solutions??? so run things again???
bandpass(vis = vis, caltable = 'cal.B1', field = bpcal, spw='', refant = refant, solnorm = True, solint = 'inf', bandtype = 'B', gaintable = ['cal.G1','cal.D0'], parang = True)
gaincal(vis = vis, caltable = 'cal.G2', field=','.join([bpcal,fluxcal,phasecal]), refant = refant, spw = '*', gaintype = 'G', calmode = 'ap', parang = True, solint = '45s', gaintable = ['cal.B1','cal.D0'])
polcal(vis=vis, caltable='cal.D1', field=phasecal, refant=refant, gaintable=['cal.B1', 'cal.G2'], poltype='Df+QU', solint='inf')

#MAKE BANDPASS PLOTS
#plotms(vis = 'cal.B1', xaxis = 'freq', yaxis = 'amp', coloraxis = 'spw') #plot the bandpass solutions
#plotms(vis='cal.G2', xaxis='time', yaxis='amp/phase', coloraxis='spw', field=bpcal, iteraxis='antenna', plotrange=[0,0,-180.0,180.0]) 

#set the flux
fluxscale(vis=vis, caltable='cal.G2', fluxtable='cal.F0', reference=fluxcal) 

#apply calibration on main cals
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, fluxcal], field=fluxcal, parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, bpcal], field=bpcal, parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, bpcal], field=phasecal, parang=True, flagbackup=False)

# apply calibration to all targets

for tgt in targets:
    applycal(vis=vis, field=tgt, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, phasecal], parang=True, flagbackup=False)

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


#Flux extraction via uvmodelfit - ALL targets, ONE combined MS
import os
import sys
 
targets_corrected_ms = 'targets_corrected.ms'
if os.path.exists(targets_corrected_ms):
    os.system(f'rm -rf {targets_corrected_ms}')
 
#need to transform for some reason
mstransform(vis=vis, outputvis=targets_corrected_ms,
            field=','.join(targets), datacolumn='corrected',
            keepflags=False)

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
 
results = []  # (field, flux_Jy, flux_err_Jy, status)
 
for tgt in targets[:5]: #Only on first 5 targets for now
    safe = tgt.replace('+', 'p').replace('-', 'm') #some weird cl business
    clfile = f'{safe}.cl'
    logfile = f'{safe}_uvmodelfit.log'
 
    if os.path.exists(clfile):
        os.system(f'rm -rf {clfile}')
 
    try:
        capture_c_level_stdout(
            logfile, uvmodelfit, #use uvmodelfit to fit for a ptsrc flux
            vis=targets_corrected_ms,
            field=tgt,                # select just this target within the combined MS
            comptype='P',              # point source
            sourcepar=[0.005, 0.0, 0.0],  # [flux Jy, RA offset, Dec offset] initial guess
            niter=10,
            spw='',                    # fit across all combined SPWs
            outfile=clfile
        )
 
        # Primary flux value
        cl.open(clfile)
        flux_val = cl.getcomponent(0)['flux']['value'][0]
        cl.close()
 
        # Formal uncertainty: parse from the captured fit summary.
        # Real line format confirmed from log: "I = 0.000981418 +/- 1.04456e-07"
        flux_err = float('nan')
        with open(logfile) as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith('I =') and '+/-' in stripped:
                    try:
                        flux_err = float(stripped.split('+/-')[1].strip().split()[0])
                    except (IndexError, ValueError):
                        pass
                    break
 
        status = 'OK' if flux_err == flux_err else 'OK_NO_ERR'  # NaN check without numpy
        results.append((tgt, flux_val, flux_err, status))
        print(f"[{status}]   {tgt}: {flux_val:.6f} +/- {flux_err:.2e} Jy")
 
    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: uvmodelfit failed - {e}")
 

#flux results to a txt file

with open(flux_results_file, 'w') as f:
    f.write(f"{'field':<25}{'flux_Jy':<25}{'flux_err_Jy':<15}{'status':<20}\n")
    for tgt, flux_val, flux_err, status in results:
        f.write(f"{tgt:<25}{flux_val:<15.8f}{flux_err:<15.8f}{status:<20}\n")
 

# image only the N brightest targets using tclean 

good_results = [r for r in results if r[3] == 'OK']
good_results.sort(key=lambda r: r[1], reverse=True)
brightest = good_results[:n_brightest_to_image]
 
print(f"[INFO] Imaging the {n_brightest_to_image} brightest sources:")
for tgt, flux_val, flux_err, _ in brightest:
    print(f"       {tgt}: {flux_val:.4f} Jy")
 
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe outside a GUI session - incase I need to run on a cluster
import matplotlib.pyplot as plt
 
for tgt, flux_val, flux_err, _ in brightest:
    safe = tgt.replace('+', 'p').replace('-', 'm')
    imagename = f'{safe}_cont'
 
    for ext in ['.image', '.model', '.residual', '.psf', '.pb', '.sumwt', '.mask']:
        if os.path.exists(imagename + ext):
            os.system(f'rm -rf {imagename}{ext}')
 
    tclean(vis=vis,
           field=tgt,
           datacolumn='corrected',   
           imagename=imagename,
           specmode='mfs',
           deconvolver='hogbom',
           imsize=512,
           cell='0.3arcsec',        #not sure what goes here
           niter=1000,
           threshold='1mJy',         #is this 50 uJy?
           weighting='briggs',
           robust=0.5,
           interactive=False)
 
    #FITS coz we love it
    fitsname = f'{imagename}.fits'
    if os.path.exists(fitsname):
        os.system(f'rm -f {fitsname}')
    exportfits(imagename=imagename + '.image', fitsimage=fitsname, overwrite=True)
 
    #imview/casaviewer is no longer available on macOS CASA for some fucking reason
    ia.open(imagename + '.image')
    pix = ia.getchunk()[:, :, 0, 0]  # axes are [RA, Dec, Stokes, Freq] -> take a 2D slice
    csys = ia.coordsys()
    increment = csys.increment()['numeric']  # radians per pixel
    ia.close()
 
    cell_arcsec = abs(increment[0]) * 206265.0  # radians -> arcsec
    npix = pix.shape[0]
    half_extent = npix / 2.0 * cell_arcsec
 
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(pix.T, origin='lower', cmap='inferno',
                    extent=[-half_extent, half_extent, -half_extent, half_extent])
    ax.invert_xaxis()  #astronomy
    ax.set_xlabel('RA offset (arcsec)')
    ax.set_ylabel('Dec offset (arcsec)')
    ax.set_title(f'{tgt}  (uvmodelfit flux = {flux_val * 1000:.3f} mJy)')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Jy/beam')
    fig.tight_layout()
    fig.savefig(f'{imagename}.png', dpi=150)
    plt.close(fig)
 
print("[DONE] Flux extraction + imaging complete.")
 
 
 