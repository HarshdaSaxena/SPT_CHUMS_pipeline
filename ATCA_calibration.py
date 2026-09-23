"""
ATCA 7mm continuum reduction pipeline for calibration
======================================
Calibrates a multi-target ATCA 7mm (43 GHz) track.

Commands to run on Etude:

export DISPLAY=:99
export QT_QPA_PLATFORM=offscreen
casa --nogui -c ATCA_calibration.py
"""

import os

vis        = '2026-06-25_1542_C3826/raw.ms' 
bpcal      = '1921-293'      # bandpass calibrator 
fluxcal    = '1934-638'      # flux calibrator
phasecal   = '2355-534'      # phase calibrator
refant     = 'ca04'          # Using CA04 as reference coz fairly central, 01 is weird w amplitudes, 2 has weird attenuators, 3 out in the field

#List all targets
obs_summary = listobs(vis=vis)
targets = sorted(set(
    f['name'] for k, f in obs_summary.items()
    if k.startswith('field_') and f['name'] not in (bpcal, fluxcal, phasecal)
))
print(f"[INFO] {len(targets)} target fields identified: {targets[:5]} ...")

if not os.path.exists('plots'):
    os.makedirs('plots')

# PLOTS - TAKE FOREVER SO SAVE ONCE 
#u-v plot, amplitude vs time plot, *& implies only cross-corr data, amp vs time takes insane memory 
plotms(vis=vis, xaxis='U', yaxis='V', antenna='*&', coloraxis='field', plotfile='plots/uv_coverage_inital.png', showgui=False, overwrite=True, highres=True)
plotms(vis=vis, xaxis='time', yaxis='amp', antenna='*&', coloraxis='baseline', plotfile='plots/ampvstime_initial.png', showgui=False, overwrite=True, highres=True)
plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', coloraxis='baseline', plotfile='plots/ampvsfreq_initial.png', showgui=False, overwrite=True, highres=True)

#Flagging

# save initial flag state
flagmanager(vis=vis, mode='save', versionname='flag_initial')
#Need to flag pointing scans
flagdata(vis=vis, mode='manual', intent='*POINT*', flagbackup=False)
#clip zeros and NaNs
flagdata(vis=vis, mode='clip', clipzeros=True, flagbackup = False) 
#flag antenna shadowing
flagdata(vis=vis, mode='shadow', tolerance = 0.0, flagbackup = False)
# remove first 10 seconds of observations
flagdata(vis=vis, mode='quack', quackinterval=10.0, quackmode='beg', flagbackup=False)
#flags all autocorr
flagdata(vis=vis, autocorr=True, flagbackup=False)
#Actual RFI flagging - looks good, default options are best for now!
flagdata(vis=vis, mode='tfcrop', extendflags=False, action='apply', display='report', writeflags=True, flagbackup=False)

#Check in plotted data if everything is well flagged 
plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', flaggedsymbolshape='circle', customflaggedsymbol=True, plotfile='plots/ampvsfreq_postflag.png', showgui=False, overwrite=True, highres=True)
#extend flagged data - if more than 80% of time range is flagged, then flag the entire chunk
flagdata(vis=vis, mode='extend', growtime=80.0, growfreq=80.0, action='apply', display='report', writeflags=True, flagbackup=False)

#save flags version 
flagmanager(vis=vis, mode='save', versionname='flag_postflag')

# calibration

setjy(vis=vis, field=fluxcal, scalebychan=True, standard='Stevens-Reynolds 2016', usescratch=True) #From Aarons comment, use this table for now
gaincal(vis=vis, caltable='bpcal.K0', field=bpcal, gaintype='K', solint='inf', refant=refant) #delay K type calibration
# do primary gain calibration fluxcal is fine here
gaincal(vis = vis, caltable = 'cal.G0', field = fluxcal, refant = refant, gaintype = 'G', calmode = 'p', parang = True, solint = '60s') 
#do bandpass response - bandtype B is sufficient 
bandpass(vis=vis, caltable='cal.B0', field=bpcal, spw='', refant=refant, solnorm=True, solint='inf', bandtype='B', gaintable=['cal.G0'], parang=True) 
#secondary gain calibraton
gaincal(vis=vis, caltable='cal.G1', field=','.join([bpcal,fluxcal,phasecal]), refant=refant, spw='*', gaintype='G', calmode='ap', parang=True, solint='45s', gaintable=['cal.B0'])
#polarization response - phasecal would be suggested since there is larger parallactic coverage, use Df+QU for phasecal and Df if doing on fluxcal
polcal(vis=vis, caltable='cal.D0', field=phasecal, refant=refant, gaintable=['cal.B0', 'cal.G1'], poltype='Df+QU', solint='inf')

#Need to do better solutions, so run things again apparently
bandpass(vis = vis, caltable = 'cal.B1', field = bpcal, spw='', refant = refant, solnorm = True, solint = 'inf', bandtype = 'B', gaintable = ['cal.G1','cal.D0'], parang = True)
gaincal(vis = vis, caltable = 'cal.G2', field=','.join([bpcal,fluxcal,phasecal]), refant = refant, spw = '*', gaintype = 'G', calmode = 'ap', parang = True, solint = '45s', gaintable = ['cal.B1','cal.D0'])
polcal(vis=vis, caltable='cal.D1', field=phasecal, refant=refant, gaintable=['cal.B1', 'cal.G2'], poltype='Df+QU', solint='inf')

#MAKE BANDPASS PLOTS
plotms(vis='cal.B1', xaxis = 'freq', yaxis = 'amp', coloraxis = 'spw', plotfile='plots/B1_ampvsfreq.png', showgui=False, overwrite=True, highres=True) #plot the bandpass solutions
plotms(vis='cal.G2', xaxis='freq', yaxis='amp', coloraxis='spw', field=bpcal, iteraxis='antenna', plotfile='plots/G2bpcal_ampvsfreq.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='freq', yaxis='amp', coloraxis='spw', field=fluxcal, iteraxis='antenna', plotfile='plots/G2fluxcal_ampvsfreq.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='freq', yaxis='amp', coloraxis='spw', field=phasecal, iteraxis='antenna', plotfile='plots/G2phasecal_ampvsfreq.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='time', yaxis='phase', coloraxis='spw', field=bpcal, iteraxis='antenna', plotfile='plots/G2bpcal_phasevstime.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='time', yaxis='phase', coloraxis='spw', field=fluxcal, iteraxis='antenna', plotfile='plots/G2fluxcal_phasevstime.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='time', yaxis='phase', coloraxis='spw', field=phasecal, iteraxis='antenna', plotfile='plots/G2phasecal_phasevstime.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='freq', yaxis='phase', coloraxis='spw', field=bpcal, iteraxis='antenna', plotfile='plots/G2bpcal_phasevsfreq.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='freq', yaxis='phase', coloraxis='spw', field=fluxcal, iteraxis='antenna', plotfile='plots/G2fluxcal_phasevsfreq.png', showgui=False, overwrite=True, highres=True) 
plotms(vis='cal.G2', xaxis='freq', yaxis='phase', coloraxis='spw', field=phasecal, iteraxis='antenna', plotfile='plots/G2phasecal_phasevsfreq.png', showgui=False, overwrite=True, highres=True) 

#set the flux
fluxscale(vis=vis, caltable='cal.G2', fluxtable='cal.F0', reference=fluxcal) 

#apply calibration on main cals - default is linear extrapolation - can change using interp. Since D1 uses fluxcal, need fluxcal as second col, and third col is the field for self consistent solns i think
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, fluxcal], field=fluxcal, interp=['linear', 'linear', 'linear'], parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, bpcal],   field=bpcal, interp=['linear', 'linear', 'linear'], parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, phasecal], field=phasecal, interp=['linear', 'linear', 'linear'], parang=True, flagbackup=False)

# apply calibration to all targets

for tgt in targets:
    applycal(vis=vis, field=tgt, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, phasecal, phasecal], interp=['linear', 'linear', 'linear'], parang=True, flagbackup=False)

#the final set of flagging routines to make sure any outliers are properly flagged - using options from BIGCAT guide now. 
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
plotms(vis = vis, field = fluxcal, xaxis = 'freq', yaxis = 'amp', correlation ='xx,yy' , ydatacolumn = 'corrected', coloraxis = 'spw', plotfile='plots/fluxcal_postcalampvsfreq.png', showgui=False, overwrite=True, highres=True)
plotms(vis=vis, field=bpcal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw', plotfile='plots/bpcal_postcalampvsfreq.png', showgui=False, overwrite=True, highres=True)
plotms(vis=vis, field=phasecal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw', plotfile='plots/phasecal_postcalampvsfreq.png', showgui=False, overwrite=True, highres=True)


# image bpcal, phasecal and fluxcal
 
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe outside a GUI session - incase I need to run on a cluster
import matplotlib.pyplot as plt

#Image calibrators, dfferent thresholds for these guys compared to targets coz much brighter
calibrators = [
    (bpcal,    'bandpass calibrator', 2000, '1mJy'),
    (fluxcal,  'flux calibrator',     2000, '1mJy'),
    (phasecal, 'phase calibrator',    1500, '1mJy'),
]

for field_name, label, niter, threshold in calibrators:
    safe = field_name.replace('+', 'p').replace('-', 'm')
    imagename = f'{safe}_cont'

    for ext in ['.image', '.model', '.residual', '.psf', '.pb', '.sumwt', '.mask']:
        if os.path.exists(imagename + ext):
            os.system(f'rm -rf {imagename}{ext}')

    tclean(vis=vis,
           field=field_name,
           datacolumn='corrected',
           imagename=imagename,
           specmode='mfs',
           deconvolver='hogbom',
           imsize=256,
           cell='0.1arcsec',       
           niter=niter,
           threshold=threshold,    
           weighting='briggs',
           robust=0.5,
           interactive=False)

    # FITS coz we love it
    fitsname = f'{imagename}.fits'
    if os.path.exists(fitsname):
        os.system(f'rm -f {fitsname}')
    exportfits(imagename=imagename + '.image', fitsimage=fitsname, overwrite=True)

    # imview/casaviewer is no longer available on macOS CASA for some reason
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
    ax.invert_xaxis()  # astronomy
    ax.set_xlabel('RA offset (arcsec)')
    ax.set_ylabel('Dec offset (arcsec)')
    ax.set_title(f'{field_name}  ({label})')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Jy/beam')
    fig.tight_layout()
    fig.savefig(f'plots/{imagename}.png', dpi=150)
    plt.close(fig)

    print(f"[INFO] Saved {field_name} -> {imagename}.png")

print("[DONE] Calibrator imaging complete.")