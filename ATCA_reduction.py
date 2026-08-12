"""
ATCA 7mm continuum reduction pipeline
======================================
Calibrates a multi-target ATCA 7mm (43 GHz) track, extracts point-source
fluxes for every science target via uvmodelfit, writes them to a txt file,
and images only the N brightest sources for a visual check.

Run inside CASA:
    casa -c atca_7mm_pipeline.py
"""

import os
import math

vis        = '2026-06-25_1542_C3826/raw.ms' 
bpcal      = '1921-293'      # bandpass calibrator 
fluxcal    = '1934-638'      # flux calibrator
phasecal   = '2355-534'      # phase calibrator
refant     = 'ca04'          # Using CA04 as reference coz fairly central, 01 is weird w amplitudes, 2 has weird attenuators, 3 out in the field

flux_results_file = 'target_fluxes_2026-06-25.txt'

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
#Actual RFI flagging - looks good, default options are best for now!
flagdata(vis=vis, mode='tfcrop', extendflags=False, action='apply', display='report', writeflags=True)

#Check in plotted data if everything is well flagged 
plotms(vis=vis, xaxis='freq', yaxis='amp', antenna='*&', flaggedsymbolshape='circle', customflaggedsymbol=True)
#extend flagged data - if more than 80% of time range is flagged, then flag the entire chunk
flagdata(vis=vis, mode='extend', growtime=80.0, growfreq=80.0, action='apply', display='report', writeflags=True)

#save flags version 
flagmanager(vis=vis, mode='save', versionname='flag_v2')

# calibration

setjy(vis=vis, field=fluxcal, scalebychan=True, standard='Stevens-Reynolds 2016', usescratch=True) #From Aarons comment, use this table for now
# do primary gain calibration - do we need bpcal or fluxcal here - did bpcal in first analysis. now fluxcal
gaincal(vis = vis, caltable = 'cal.G0', field = fluxcal, refant = refant, gaintype = 'G', calmode = 'p', parang = True, solint = '60s') 
#do bandpass response
bandpass(vis=vis, caltable='cal.B0', field=bpcal, spw='', refant=refant, solnorm=True, solint='inf', bandtype='B', gaintable=['cal.G0'], parang=True) 
#secondary gain calibraton
gaincal(vis=vis, caltable='cal.G1', field=','.join([bpcal,fluxcal,phasecal]), refant=refant, spw='*', gaintype='G', calmode='ap', parang=True, solint='45s', gaintable=['cal.B0'])
#polarization response - do we do phasecal or fluxcal - did phasecal before, now fluxcal
polcal(vis=vis, caltable='cal.D0', field=fluxcal, refant=refant, gaintable=['cal.B0', 'cal.G1'], poltype='Df', solint='inf')

#Need to do better solutions, so run things again apparently
bandpass(vis = vis, caltable = 'cal.B1', field = bpcal, spw='', refant = refant, solnorm = True, solint = 'inf', bandtype = 'B', gaintable = ['cal.G1','cal.D0'], parang = True)
gaincal(vis = vis, caltable = 'cal.G2', field=','.join([bpcal,fluxcal,phasecal]), refant = refant, spw = '*', gaintype = 'G', calmode = 'ap', parang = True, solint = '45s', gaintable = ['cal.B1','cal.D0'])
polcal(vis=vis, caltable='cal.D1', field=fluxcal, refant=refant, gaintable=['cal.B1', 'cal.G2'], poltype='Df', solint='inf')

#MAKE BANDPASS PLOTS
#plotms(vis = 'cal.B1', xaxis = 'freq', yaxis = 'amp', coloraxis = 'spw') #plot the bandpass solutions
#plotms(vis='cal.G2', xaxis='time/freq', yaxis='amp/phase', coloraxis='spw', field=bpcal/fluxcal/phasecal, iteraxis='antenna') 

#set the flux
fluxscale(vis=vis, caltable='cal.G2', fluxtable='cal.F0', reference=fluxcal) 

#apply calibration on main cals - default is linear extrapolation - can change using interp. Since D1 uses fluxcal, need fluxcal as second col, and third col is the field for self consistent solns i think
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, fluxcal], field=fluxcal, parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, bpcal],   field=bpcal,   parang=True, flagbackup=False)
applycal(vis=vis, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, phasecal], field=phasecal, parang=True, flagbackup=False)

# apply calibration to all targets

for tgt in targets:
    applycal(vis=vis, field=tgt, gaintable=['cal.B1', 'cal.D1', 'cal.F0'], gainfield=[bpcal, fluxcal, phasecal], parang=True, flagbackup=False)

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
#plotms(vis = vis, field = fluxcal, xaxis = 'freq', yaxis = 'amp', correlation ='xx,yy' , ydatacolumn = 'corrected', coloraxis = 'spw')
#plotms(vis=vis, field=bpcal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw')
#plotms(vis=vis, field=phasecal, xaxis='freq', yaxis='amp', correlation='xx,yy', ydatacolumn='corrected', coloraxis='spw')


#Flux extraction via uvmodelfit - ALL targets
import os
import sys

statwt(vis=vis, datacolumn='corrected') #do weighting correctly - thanks Zac + Adam
 
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

MAX_OFFSET_ARCSEC = 10.0  #Can change for now this is limit for max peak offset from the centers we have

def run_uvmodelfit_and_parse(vis, field, clfile, logfile, sourcepar, varypar=None):
    """Runs uvmodelfit once, returns (flux, flux_err, chi2, x, x_err, y, y_err) or Nones on failure."""
    kwargs = dict(vis=vis, field=field, comptype='P', sourcepar=sourcepar,
                   niter=10, spw='', outfile=clfile)
    if varypar is not None:
        kwargs['varypar'] = varypar

    capture_c_level_stdout(logfile, uvmodelfit, **kwargs)

    cl.open(clfile)
    flux_val = cl.getcomponent(0)['flux']['value'][0]
    cl.close()

    flux_err = float('nan')
    chi2 = float('nan')
    x_off, x_err = float('nan'), float('nan')
    y_off, y_err = float('nan'), float('nan')

    with open(logfile) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('iter='):
                try:
                    chi2 = float(stripped.split('reduced chi2=')[1].split(':')[0].strip())
                except (IndexError, ValueError):
                    pass
            elif stripped.startswith('I =') and '+/-' in stripped:
                try:
                    flux_err = float(stripped.split('+/-')[1].strip().split()[0])
                except (IndexError, ValueError):
                    pass
            elif stripped.startswith('x =') and '+/-' in stripped:
                try:
                    rhs = stripped.split('=')[1].split('+/-')
                    x_off = float(rhs[0].strip())
                    x_err = float(rhs[1].strip().split()[0])
                except (IndexError, ValueError):
                    pass
            elif stripped.startswith('y =') and '+/-' in stripped:
                try:
                    rhs = stripped.split('=')[1].split('+/-')
                    y_off = float(rhs[0].strip())
                    y_err = float(rhs[1].strip().split()[0])
                except (IndexError, ValueError):
                    pass

    return flux_val, flux_err, chi2, x_off, x_err, y_off, y_err


results = []  # (field, flux_Jy, flux_err_Jy, chi2, x_arcsec, x_err, y_arcsec, y_err, status)

for tgt in targets:
    safe = tgt.replace('+', 'p').replace('-', 'm')
    clfile = f'{safe}.cl'
    logfile = f'{safe}_uvmodelfit.log'

    if os.path.exists(clfile):
        os.system(f'rm -rf {clfile}')

    try:
        #free fit - position allowed to vary
        flux_val, flux_err, chi2, x_off, x_err, y_off, y_err = run_uvmodelfit_and_parse(
            targets_corrected_ms, tgt, clfile, logfile, sourcepar=[0.005, 0.0, 0.0])

        offset_radius = math.sqrt(x_off**2 + y_off**2) if x_off == x_off else float('nan')  # NaN-safe

        if offset_radius == offset_radius and offset_radius > MAX_OFFSET_ARCSEC:
            # Free fit wandered too far - Refit with position FIXED at phase center
            print(f"[REFIT] {tgt}: free fit offset {offset_radius:.1f}\" exceeds "
                  f"{MAX_OFFSET_ARCSEC}\" limit - refitting with position fixed at (0,0)")
            flux_val, flux_err, chi2, x_off, x_err, y_off, y_err = run_uvmodelfit_and_parse(
                targets_corrected_ms, tgt, clfile, logfile,
                sourcepar=[flux_val if flux_val == flux_val else 0.005, 0.0, 0.0],
                varypar=[True, False, False])
            status = 'OK_POS_CONSTRAINED' if flux_err == flux_err else 'OK_POS_CONSTRAINED_NO_ERR'
        else:
            status = 'OK' if flux_err == flux_err else 'OK_NO_ERR'

        results.append((tgt, flux_val, flux_err, chi2, x_off, x_err, y_off, y_err, status))
        print(f"[{status}]   {tgt}: {flux_val:.6f} +/- {flux_err:.2e} Jy, "
              f"chi2={chi2:.3g}, offset=({x_off:.2f},{y_off:.2f}) arcsec")

    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), float('nan'),
                         float('nan'), float('nan'), float('nan'), float('nan'),
                         f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: uvmodelfit failed - {e}")

# ----------------------------------------------------------------------
# Write extended flux results to a txt file
# ----------------------------------------------------------------------
with open(flux_results_file, 'w') as f:
    header = (f"{'field':<25}{'flux_Jy':<15}{'flux_err_Jy':<15}{'chi2':<12}"
              f"{'x_arcsec':<12}{'x_err':<12}{'y_arcsec':<12}{'y_err':<12}{'status':<22}\n")
    f.write(header)
    for tgt, flux_val, flux_err, chi2, x_off, x_err, y_off, y_err, status in results:
        f.write(f"{tgt:<25}{flux_val:<15.6e}{flux_err:<15.6e}{chi2:<12.4g}"
                f"{x_off:<12.4f}{x_err:<12.4f}{y_off:<12.4f}{y_err:<12.4f}{status:<22}\n")

print(f"[INFO] Extended flux results for {len(results)} targets written to {flux_results_file}")


# image bpcal, phasecal and fluxcal
 
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, safe outside a GUI session - incase I need to run on a cluster
import matplotlib.pyplot as plt

#Image calibrators, dfferent thresholds for these guys compared to targets coz much brighter
calibrators = [
    (bpcal,    'bandpass calibrator', 2000, '10mJy'),
    (fluxcal,  'flux calibrator',     2000, '5mJy'),
    (phasecal, 'phase calibrator',    1500, '5mJy'),
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
           cell='1arcsec',        #maybe for now?
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
    fig.savefig(f'{imagename}.png', dpi=150)
    plt.close(fig)

    print(f"[INFO] Saved {field_name} -> {imagename}.png")

print("[DONE] Calibrator imaging complete.")

# Image brightest sources

def plot_casa_image(image_path, title, cbar_label, png_out):
    """Read a CASA image with the ia tool and save a matplotlib PNG.
    Reused for the clean, residual, and dirty beam (.psf) products."""
    ia.open(image_path)
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
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    fig.savefig(png_out, dpi=150)
    plt.close(fig)

n_brightest_to_image = 5
good_results = [r for r in results if r[-1] == 'OK']
good_results.sort(key=lambda r: r[1], reverse=True)
brightest = good_results[:n_brightest_to_image]

print(f"[INFO] Imaging the {n_brightest_to_image} brightest sources:")
for tgt, flux_val, flux_err, _, _, _, _, _, _ in brightest:
    print(f"{tgt}: {flux_val:.4f} Jy")

for tgt, flux_val, flux_err, _, _, _, _, _, _ in brightest:
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
           cell='1arcsec',
           imsize=256,          #letting this and cell size vary makes the image v shitty - check what is the best?
           niter=1000,
           threshold='0.25mJy',
           weighting='briggs',
           robust=0.5,
           interactive=False)

    # FITS coz we love it
    fitsname = f'{imagename}.fits'
    if os.path.exists(fitsname):
        os.system(f'rm -f {fitsname}')
    exportfits(imagename=imagename + '.image', fitsimage=fitsname, overwrite=True)

    #Clean image
    plot_casa_image(imagename + '.image',
                     f'{tgt}  (uvmodelfit flux = {flux_val * 1000:.3f} mJy)',
                     'Jy/beam',
                     f'{imagename}.png')

    #Residual map 
    plot_casa_image(imagename + '.residual',
                     f'{tgt}  (residual)',
                     'Jy/beam',
                     f'{imagename}_residual.png')

    #Dirty beam / PSF
    plot_casa_image(imagename + '.psf',
                     f'{tgt}  (dirty beam / PSF)',
                     'Normalized response',
                     f'{imagename}_psf.png')

    print(f"[INFO] Saved clean/residual/psf images for {tgt}")

print("[DONE] Flux extraction + imaging complete.")
 
 
