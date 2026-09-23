"""
ATCA 7mm flux estimation pipeline
======================================
Estimates the flux via a uvmodel fit and then uses phase closure to estimate the fluxes of 
significant sources (>3sigma). 

Commands to run on Etude:

export DISPLAY=:99
export QT_QPA_PLATFORM=offscreen
casa --nogui -c ATCA_flux_estimation.py
"""
'''
import os
import math
import itertools
import numpy as np
import pandas as pd

vis        = '2026-06-25_1542_C3826/raw.ms' 
bpcal      = '1921-293'      # bandpass calibrator 
fluxcal    = '1934-638'      # flux calibrator
phasecal   = '2355-534'      # phase calibrator
refant     = 'ca04'          # Using CA04 as reference coz fairly central, 01 is weird w amplitudes, 2 has weird attenuators, 3 out in the field

targets_corrected_ms = 'targets_corrected.ms'

obs_summary = listobs(vis=vis)
targets = sorted(set(
    f['name'] for k, f in obs_summary.items()
    if k.startswith('field_') and f['name'] not in (bpcal, fluxcal, phasecal)
))
print(f"[INFO] {len(targets)} target fields identified: {targets[:5]} ...")

# Only redo the (slow) statwt + mstransform step if the output doesn't already exist
if not os.path.exists(targets_corrected_ms):
    statwt(vis=vis, datacolumn='corrected')
    mstransform(vis=vis, outputvis=targets_corrected_ms,
                field=','.join(targets), datacolumn='corrected',
                keepflags=False)
else:
    print(f"[INFO] {targets_corrected_ms} already exists, skipping statwt/mstransform")


#FLUX EXTRACTION VIA UVMODELFIT
#function to capture the terminal output because what it outputs to a normal txt file was bs
flux_results_file = 'target_fluxes_2026-06-25_uvmodelfit.txt'

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

#THIS IS FITTING ALL BASELINES SIMULATENOUSLY 
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

##READ SIGNIFICANT SOURCES FROM UVMODELFIT AND THEN RE-RUN PHASE CLOSURE ON THEM
#This is when we first average on time and then compute triple amplitude  - 
# ----------------------------------------------------------------------
# Filter targets using the uvmodelfit results: only keep sources with
# positive flux and >3sigma significance (flux/flux_err > 3), since
# phase closure on a non-detection is not meaningful to compute.
# ----------------------------------------------------------------------
uvmodelfit_results_file = 'target_fluxes_2026-06-25_uvmodelfit.txt'
uvmodelfit_widths = [25, 15, 15, 12, 12, 12, 12, 12, 22]

df_uvfit = pd.read_fwf(uvmodelfit_results_file, widths=uvmodelfit_widths)
df_uvfit.columns = [c.strip() for c in df_uvfit.columns]
df_uvfit['field'] = df_uvfit['field'].str.strip()

# guard against div-by-zero / NaN errors when computing significance
with np.errstate(divide='ignore', invalid='ignore'):
    sigma = df_uvfit['flux_Jy'] / df_uvfit['flux_err_Jy']

detected_mask = (df_uvfit['flux_Jy'] > 0) & (sigma > 3) & np.isfinite(sigma)
detected_targets = sorted(df_uvfit.loc[detected_mask, 'field'].tolist())

print(f"[INFO] {len(detected_targets)} / {len(df_uvfit)} targets passed "
      f">3sigma, positive-flux cut from uvmodelfit and will be processed "
      f"via phase closure: {detected_targets[:5]} ...")

# ----------------------------------------------------------------------
# Flux extraction via closure triple-product (bispectrum), with
# theoretical (radiometer/weight-propagated) flux errors alongside
# the empirical (across-triangle scatter) errors
# ----------------------------------------------------------------------
DATACOL_MAP = {'data': 'data', 'corrected': 'corrected_data', 'model': 'model_data'}

def get_field_baseline_vis(vis, field, datacolumn='data', spw=''):
    """
    Returns {(ant1, ant2): (weighted_mean_vis, sigma_theoretical, n)}
    averaged over all unflagged channels/times, using the mean of the
    parallel-hand polarizations (assumed to be pol index 0 and -1),
    weighted by the WEIGHT column (1/sigma^2) populated by statwt.
    """
    col = DATACOL_MAP[datacolumn]
    ms.open(vis)
    ms.selectinit(reset=True)
    sel = {'field': field}
    if spw:
        sel['spw'] = spw
    ms.msselect(sel)

    rec = ms.getdata([col, 'flag', 'antenna1', 'antenna2', 'weight'])

    ms.close()

    data = rec[col]          # shape (npol, nchan, nrow)
    flag = rec['flag']
    ant1 = rec['antenna1']
    ant2 = rec['antenna2']
    weight = rec['weight']   # shape (npol, nrow) -- one value per row per corr

    npol = data.shape[0]
    if npol >= 2:
        data_i = 0.5 * (data[0] + data[-1])
        flag_i = flag[0] | flag[-1]
        w_XX = weight[0, :]
        w_YY = weight[-1, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            var_i_row = 0.25 * (np.where(w_XX > 0, 1.0 / w_XX, np.inf) +
                                 np.where(w_YY > 0, 1.0 / w_YY, np.inf))
        w_i_row = np.where(np.isfinite(var_i_row) & (var_i_row > 0), 1.0 / var_i_row, 0.0)
    else:
        data_i = data[0]
        flag_i = flag[0]
        w_i_row = weight[0, :]

    sum_wv, sum_w, counts = {}, {}, {}
    nrow = data_i.shape[-1]
    for row in range(nrow):
        a1, a2 = int(ant1[row]), int(ant2[row])
        if a1 == a2:
            continue
        vals = data_i[:, row]
        good = ~flag_i[:, row]
        w_row = w_i_row[row]
        if not np.any(good) or w_row <= 0:
            continue
        key = (a1, a2)
        n_good = int(np.sum(good))
        v_sum = np.sum(vals[good])
        if key in sum_wv:
            sum_wv[key] += w_row * v_sum
            sum_w[key] += w_row * n_good
            counts[key] += n_good
        else:
            sum_wv[key] = w_row * v_sum
            sum_w[key] = w_row * n_good
            counts[key] = n_good

    out = {}
    for key in sum_wv:
        w_tot = sum_w[key]
        v_mean = sum_wv[key] / w_tot
        sigma_theoretical = 1.0 / np.sqrt(w_tot)   # theoretical noise on the time-averaged vis
        out[key] = (v_mean, sigma_theoretical, counts[key])
    return out


def closure_triples(baseline_vis):
    """
    Given {(a1,a2): (vis, sigma, n)}, returns a list of
    (triangle, flux_estimate_Jy, flux_err_theoretical_Jy,
     closure_phase_deg, n_min) for every closed antenna triangle.
    """
    def get_vis(a, b):
        if (a, b) in baseline_vis:
            v, s, n = baseline_vis[(a, b)]
            return v, s, n
        elif (b, a) in baseline_vis:
            v, s, n = baseline_vis[(b, a)]
            return np.conj(v), s, n
        return None, None, None

    ants = sorted(set(a for pair in baseline_vis for a in pair))
    out = []
    for i, j, k in itertools.combinations(ants, 3):
        v_ij, s_ij, n_ij = get_vis(i, j)
        v_jk, s_jk, n_jk = get_vis(j, k)
        v_ik, s_ik, n_ik = get_vis(i, k)
        if v_ij is None or v_jk is None or v_ik is None:
            continue

        a_ij, a_jk, a_ik = np.abs(v_ij), np.abs(v_jk), np.abs(v_ik)
        if a_ij == 0 or a_jk == 0 or a_ik == 0:
            continue

        bispectrum = v_ij * v_jk * np.conj(v_ik)
        flux_est = np.abs(bispectrum) ** (1.0 / 3.0)
        closure_phase = np.degrees(np.angle(bispectrum))

        frac_err_sq = (s_ij / a_ij) ** 2 + (s_jk / a_jk) ** 2 + (s_ik / a_ik) ** 2
        flux_err_theoretical = (flux_est / 3.0) * np.sqrt(frac_err_sq)

        out.append(((i, j, k), flux_est, flux_err_theoretical, closure_phase, min(n_ij, n_jk, n_ik)))
    return out


def combine_flux_theoretical(triples):
    """
    Combines per-triangle flux estimates two ways:
    - inverse-variance-weighted mean using the theoretical (radiometer) errors
    - median + scatter-based error across triangles (empirical)
    """
    fluxes = np.array([t[1] for t in triples])
    errs_theo = np.array([t[2] for t in triples])
    phases = np.array([t[3] for t in triples])

    w = 1.0 / errs_theo**2
    flux_theoretical = np.sum(w * fluxes) / np.sum(w)
    flux_err_theoretical = 1.0 / np.sqrt(np.sum(w))

    flux_empirical = np.median(fluxes)
    flux_err_empirical = np.std(fluxes) / np.sqrt(len(fluxes))

    return {
        'flux_theoretical': flux_theoretical,
        'flux_err_theoretical': flux_err_theoretical,
        'flux_empirical': flux_empirical,
        'flux_err_empirical': flux_err_empirical,
        'mean_cphase': np.mean(phases),
        'cphase_scatter': np.std(phases),
        'n_triangles': len(triples),
    }


flux_results_file = 'target_fluxes_2026-06-25_phaseclosure.txt'
results = []  # (field, flux_theo_Jy, flux_err_theo_Jy, flux_emp_Jy, flux_err_emp_Jy,
              #  err_ratio, mean_cphase_deg, cphase_scatter_deg, n_triangles, status)

for tgt in detected_targets:
    try:
        baseline_vis = get_field_baseline_vis(targets_corrected_ms, tgt, datacolumn='data')
        triples = closure_triples(baseline_vis)

        if len(triples) == 0:
            results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                             float('nan'), float('nan'), float('nan'), 0, 'NO_TRIANGLES'))
            print(f"[FAIL] {tgt}: no closed antenna triangles found")
            continue

        combined = combine_flux_theoretical(triples)
        err_ratio = combined['flux_err_empirical'] / combined['flux_err_theoretical']

        status = 'OK' if len(triples) == 20 else 'OK_FEW_TRIANGLES'
        results.append((tgt, combined['flux_theoretical'], combined['flux_err_theoretical'],
                         combined['flux_empirical'], combined['flux_err_empirical'], err_ratio,
                         combined['mean_cphase'], combined['cphase_scatter'],
                         combined['n_triangles'], status))
        print(f"[{status}]   {tgt}: theo={combined['flux_theoretical']:.6f} +/- "
              f"{combined['flux_err_theoretical']:.2e} Jy, "
              f"emp={combined['flux_empirical']:.6f} +/- {combined['flux_err_empirical']:.2e} Jy "
              f"(emp/theo err ratio={err_ratio:.2f}), "
              f"closure phase = {combined['mean_cphase']:.2f} +/- {combined['cphase_scatter']:.2f} deg "
              f"({combined['n_triangles']} triangles)")

    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                         float('nan'), float('nan'), float('nan'), 0, f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: closure-phase flux extraction failed - {e}")

# ----------------------------------------------------------------------
# Write closure-phase flux results to a txt file
# ----------------------------------------------------------------------
with open(flux_results_file, 'w') as f:
    header = (f"{'field':<25}{'flux_theo_Jy':<15}{'flux_err_theo_Jy':<18}"
              f"{'flux_emp_Jy':<15}{'flux_err_emp_Jy':<17}{'err_ratio':<12}"
              f"{'mean_cphase_deg':<18}{'cphase_scatter_deg':<20}"
              f"{'n_triangles':<14}{'status':<22}\n")
    f.write(header)
    for (tgt, flux_theo, flux_err_theo, flux_emp, flux_err_emp, err_ratio,
         mean_cphase, cphase_scatter, n_tri, status) in results:
        f.write(f"{tgt:<25}{flux_theo:<15.6e}{flux_err_theo:<18.6e}"
                f"{flux_emp:<15.6e}{flux_err_emp:<17.6e}{err_ratio:<12.3f}"
                f"{mean_cphase:<18.4f}{cphase_scatter:<20.4f}"
                f"{n_tri:<14d}{status:<22}\n")

print(f"[INFO] Closure-phase flux results for {len(results)} targets written to {flux_results_file}")
'''
'''
#THIS IS FITTING SMALL BASELINES FIRST AND THEN INFORMING FIT 
MAX_OFFSET_ARCSEC = 10.0  # limit for max peak offset from center

# Placeholder 
SHORT_BASELINE_UVRANGE = '0~250m'

flux_results_file_shortbaselineprior = 'target_fluxes_2026-06-25_shortb_prior.txt'

def run_uvmodelfit_and_parse(vis, field, clfile, logfile, sourcepar, varypar=None, uvrange=''):
    """Runs uvmodelfit once, returns (flux, flux_err, chi2, x, x_err, y, y_err) or Nones on failure."""
    kwargs = dict(vis=vis, field=field, comptype='P', sourcepar=sourcepar,
                   niter=10, spw='', outfile=clfile)
    if uvrange:
        kwargs['uvrange'] = uvrange
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
        # ------------------------------------------------------------
        # Pass 0: short baselines only - rough position estimate.
        # ------------------------------------------------------------
        short_clfile = f'{safe}_short.cl'
        short_logfile = f'{safe}_short_uvmodelfit.log'
        if os.path.exists(short_clfile):
            os.system(f'rm -rf {short_clfile}')

        _, _, _, x_guess, _, y_guess, _ = run_uvmodelfit_and_parse(
            targets_corrected_ms, tgt, short_clfile, short_logfile,
            sourcepar=[0.005, 0.0, 0.0], uvrange=SHORT_BASELINE_UVRANGE)

        if x_guess != x_guess or y_guess != y_guess:  # NaN check - short-baseline fit itself failed
            print(f"[WARN] {tgt}: short-baseline position fit failed, "
                  f"falling back to (0,0) starting guess for the full fit")
            x_guess, y_guess = 0.0, 0.0

        # ------------------------------------------------------------
        # Pass 1: free fit across all baselines, seeded with the
        # short-baseline position instead of a blind guess
        # ------------------------------------------------------------
        flux_val, flux_err, chi2, x_off, x_err, y_off, y_err = run_uvmodelfit_and_parse(
            targets_corrected_ms, tgt, clfile, logfile,
            sourcepar=[0.005, x_guess, y_guess])

        offset_radius = math.sqrt(x_off**2 + y_off**2) if x_off == x_off else float('nan')  # NaN-safe

        if offset_radius == offset_radius and offset_radius > MAX_OFFSET_ARCSEC:
            # Full-baseline fit still wandered too far - refit with
            # position FIXED at the short-baseline estimate rather than
            # blindly at phase center, since that estimate is likely
            # closer to the truth than (0,0) even if imperfect.
            print(f"[REFIT] {tgt}: full-baseline fit offset {offset_radius:.1f}\" exceeds "
                  f"{MAX_OFFSET_ARCSEC}\" limit - refitting with position fixed at "
                  f"short-baseline estimate ({x_guess:.2f},{y_guess:.2f})")
            flux_val, flux_err, chi2, x_off, x_err, y_off, y_err = run_uvmodelfit_and_parse(
                targets_corrected_ms, tgt, clfile, logfile,
                sourcepar=[flux_val if flux_val == flux_val else 0.005, x_guess, y_guess],
                varypar=[True, False, False])
            status = 'OK_POS_CONSTRAINED' if flux_err == flux_err else 'OK_POS_CONSTRAINED_NO_ERR'
        else:
            status = 'OK' if flux_err == flux_err else 'OK_NO_ERR'

        results.append((tgt, flux_val, flux_err, chi2, x_off, x_err, y_off, y_err, status))
        print(f"[{status}]   {tgt}: {flux_val:.6f} +/- {flux_err:.2e} Jy, "
              f"chi2={chi2:.3g}, offset=({x_off:.2f},{y_off:.2f}) arcsec "
              f"(short-baseline guess was ({x_guess:.2f},{y_guess:.2f}))")

    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), float('nan'),
                         float('nan'), float('nan'), float('nan'), float('nan'),
                         f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: uvmodelfit failed - {e}")

# ----------------------------------------------------------------------
# Write extended flux results to a txt file
# ----------------------------------------------------------------------
with open(flux_results_file_shortbaselineprior, 'w') as f:
    header = (f"{'field':<25}{'flux_Jy':<15}{'flux_err_Jy':<15}{'chi2':<12}"
              f"{'x_arcsec':<12}{'x_err':<12}{'y_arcsec':<12}{'y_err':<12}{'status':<22}\n")
    f.write(header)
    for tgt, flux_val, flux_err, chi2, x_off, x_err, y_off, y_err, status in results:
        f.write(f"{tgt:<25}{flux_val:<15.6e}{flux_err:<15.6e}{chi2:<12.4g}"
                f"{x_off:<12.4f}{x_err:<12.4f}{y_off:<12.4f}{y_err:<12.4f}{status:<22}\n")

print(f"[INFO] Extended flux results for {len(results)} targets written to {flux_results_file_shortbaselineprior}")
#
'''
'''
#THIS IS FOR FITTING FROM THE CLEAN MAPS
#Flux extraction via tclean + imfit (fitting the CLEAN maps) - ALL targets
#Meant to be compared against the uvmodelfit fluxes in target_fluxes_2026-06-25.txt
# No masks/regions restrict the CLEANing or the fitting anywhere below -
# tclean cleans the whole image, imfit searches the whole image for the peak.
# The 10" figure is used ONLY to label the status column afterwards.
# Only ONE image per target is saved: the full, uncropped CLEAN map.
#Flux extraction via tclean + imfit (fitting the CLEAN maps) - ALL targets
#Meant to be compared against the uvmodelfit fluxes in target_fluxes_2026-06-25.txt
#
# No masks/regions restrict the CLEANing or the fitting anywhere below -
# tclean cleans the whole image, imfit searches the whole image for the peak.
# The 10" figure is used ONLY to label the status column afterwards.
# Only ONE image per target is saved: the full, uncropped CLEAN map.
import os
import shutil
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless backend - no display needed, works on macOS
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Config - reuses vis/targets/refant etc from the calibration script,
# and targets_corrected_ms from the uvmodelfit script.
# ----------------------------------------------------------------------
date_str                = '2026-06-25'
clean_flux_results_file = f'target_fluxes_{date_str}_clean.txt'
clean_maps_dir          = f'clean_maps_{date_str}'   # one full, uncropped CLEAN map per target, fits+png

os.makedirs(clean_maps_dir, exist_ok=True)

# imaging params - adjust to match your array config / target field of view
#imsize            = 512
cell_arcsec       = 0.1
cell              = [f'{cell_arcsec}arcsec']
OFFSET_FLAG_ARCSEC = 10.0   # purely a label threshold, not used to restrict anything


def export_png(imagename, png_path, title=None):
    """Raster a CASA image to a png (imview is EOL on macOS), plotting every
    pixel of imagename with no region/crop applied."""
    try:
        ia.open(imagename)
        try:
            shape = ia.shape()
            data  = ia.getchunk(blc=[0] * len(shape), trc=[s - 1 for s in shape])
            unit  = ia.brightnessunit()
        finally:
            ia.close()

        img2d = np.squeeze(data)
        while img2d.ndim > 2:
            img2d = img2d[..., 0]        # take stokes=0 / freq=0 if present
        img2d = img2d.T                  # ia.getchunk is [x,y]; imshow wants [row=y, col=x]
        img2d = np.where(np.isfinite(img2d), img2d, np.nan)

        fig, ax = plt.subplots(figsize=(6, 6))
        im = ax.imshow(img2d, origin='lower', cmap='inferno')
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(unit if unit else 'Jy/beam')
        ax.set_xlabel('pixel x')
        ax.set_ylabel('pixel y')
        ax.set_title(title if title else os.path.basename(imagename))
        fig.tight_layout()
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[WARN] could not export png for {imagename}: {e}")


def fit_peak(imagename):
    """imfit over the WHOLE image, no region restriction. Returns the
    component dict, or None if nothing was found / it failed outright."""
    try:
        fit_result = imfit(imagename=imagename)
    except Exception as e:
        print(f"[INFO] imfit failed on {imagename}: {e}")
        return None
    if not fit_result or 'results' not in fit_result or 'component0' not in fit_result['results']:
        return None
    return fit_result['results']['component0']


def phase_centre(imagename):
    """The image's own reference direction (radians) - i.e. its phase centre."""
    ia.open(imagename)
    try:
        cs = ia.coordsys()
        world = cs.referencevalue()['numeric']
        return world[0], world[1]
    finally:
        ia.close()


def sky_offset_arcsec(ra, dec, ra_ref, dec_ref):
    """Proper tangent-plane offset in arcsec (cos(dec) factor on RA), not a
    naive coordinate subtraction. Wraps the RA difference into (-pi, pi]
    first - imfit's fitted RA and the image's reference RA can end up on
    opposite sides of the 0/2pi branch cut even when the true sky positions
    are almost identical, which otherwise produces a bogus ~360 degree
    "offset" instead of the real (usually tiny) one."""
    RAD2ARCSEC = 206264.806
    dra_raw = ra - ra_ref
    dra_wrapped = (dra_raw + math.pi) % (2 * math.pi) - math.pi
    dra  = dra_wrapped * math.cos(dec_ref) * RAD2ARCSEC
    ddec = (dec - dec_ref) * RAD2ARCSEC
    return dra, ddec


results = []  # (field, flux_Jy, flux_err_Jy, ra_off_arcsec, dec_off_arcsec, status)

for tgt in targets:
    safe       = tgt.replace('+', 'p').replace('-', 'm')
    imagename  = f'{safe}_clean'
    image_file = f'{imagename}.image'
    pbcor_file = f'{imagename}.image.pbcor'

    for ext in ['.image', '.image.pbcor', '.model', '.residual', '.psf', '.pb', '.mask', '.sumwt']:
        if os.path.exists(imagename + ext):
            shutil.rmtree(imagename + ext, ignore_errors=True)

    flux_val, flux_err = float('nan'), float('nan')
    ra_off_arcsec, dec_off_arcsec = float('nan'), float('nan')
    status = 'FIT_FAIL'

    try:
        # CLEAN the whole image - no mask, no restriction
        tclean(vis=targets_corrected_ms,
               imagename=imagename,
               field=tgt,
               spw='',
               specmode='mfs',
               deconvolver='hogbom',
               gridder='standard',
               #imsize=imsize,
               cell=cell,
               weighting='briggs',
               robust=0.5,
               niter=1000,
               threshold='0.1mJy',
               interactive=False)

        use_image = pbcor_file if os.path.exists(pbcor_file) else image_file
        # plot the plain (non-pbcor) image - .pbcor gets NaN-blanked outside
        # the primary beam, which looks like a crop even though it isn't one
        plot_image = image_file if os.path.exists(image_file) else use_image

        clean_fits = os.path.join(clean_maps_dir, f'{safe}_clean.fits')
        exportfits(imagename=use_image, fitsimage=clean_fits, overwrite=True)
        export_png(plot_image, os.path.join(clean_maps_dir, f'{safe}_clean.png'), title=tgt)

        # fit the peak over the WHOLE image - no region restriction
        comp = fit_peak(use_image)

        if comp is not None:
            flux_val = comp['flux']['value'][0]
            flux_err = comp['flux']['error'][0]
            ra  = comp['shape']['direction']['m0']['value']
            dec = comp['shape']['direction']['m1']['value']

            ra_ref, dec_ref = phase_centre(use_image)
            ra_off_arcsec, dec_off_arcsec = sky_offset_arcsec(ra, dec, ra_ref, dec_ref)
            total_offset = math.hypot(ra_off_arcsec, dec_off_arcsec)

            status = 'OK' if total_offset <= OFFSET_FLAG_ARCSEC else f'OK_OFFSET_GT_{OFFSET_FLAG_ARCSEC:.0f}ARCSEC'

            print(f"[{status}]   {tgt}: {flux_val:.6f} +/- {flux_err:.2e} Jy, "
                  f"offset=({ra_off_arcsec:.2f},{dec_off_arcsec:.2f})\" total={total_offset:.2f}\"")
        else:
            status = 'FIT_FAIL: no component found'
            print(f"[FAIL] {tgt}: {status}")

        results.append((tgt, flux_val, flux_err, ra_off_arcsec, dec_off_arcsec, status))

    except Exception as e:
        status = f'FIT_FAIL: {e}'
        results.append((tgt, flux_val, flux_err, ra_off_arcsec, dec_off_arcsec, status))
        print(f"[FAIL] {tgt}: CLEAN/imfit pipeline failed - {e}")

# ----------------------------------------------------------------------
# Write CLEAN-map flux results to a txt file
# ----------------------------------------------------------------------
with open(clean_flux_results_file, 'w') as f:
    header = (f"{'field':<25}{'flux_Jy':<15}{'flux_err_Jy':<15}"
              f"{'ra_off_arcsec':<16}{'dec_off_arcsec':<16}{'status':<30}\n")
    f.write(header)
    for tgt, flux_val, flux_err, ra_off, dec_off, status in results:
        f.write(f"{tgt:<25}{flux_val:<15.6e}{flux_err:<15.6e}"
                f"{ra_off:<16.4f}{dec_off:<16.4f}{status:<30}\n")

print(f"[INFO] CLEAN-map flux results for {len(results)} targets written to {clean_flux_results_file}")
print(f"[INFO] One full, uncropped CLEAN map per target (fits+png) saved to ./{clean_maps_dir}/")
'''
'''
#THIS IS GETTING FLUXES FROM PHASE CLOSURE
#This is we first average on time and then compute triple amplitude  - 
import itertools
import numpy as np

# ----------------------------------------------------------------------
# Flux extraction via closure triple-product (bispectrum), with
# theoretical (radiometer/weight-propagated) flux errors alongside
# the empirical (across-triangle scatter) errors
# ----------------------------------------------------------------------
DATACOL_MAP = {'data': 'data', 'corrected': 'corrected_data', 'model': 'model_data'}

def get_field_baseline_vis(vis, field, datacolumn='data', spw=''):
    """
    Returns {(ant1, ant2): (weighted_mean_vis, sigma_theoretical, n)}
    averaged over all unflagged channels/times, using the mean of the
    parallel-hand polarizations (assumed to be pol index 0 and -1),
    weighted by the WEIGHT column (1/sigma^2) populated by statwt.
    """
    col = DATACOL_MAP[datacolumn]
    ms.open(vis)
    ms.selectinit(reset=True)
    sel = {'field': field}
    if spw:
        sel['spw'] = spw
    ms.msselect(sel)

    rec = ms.getdata([col, 'flag', 'antenna1', 'antenna2', 'weight'])

    ms.close()

    data = rec[col]          # shape (npol, nchan, nrow)
    flag = rec['flag']
    ant1 = rec['antenna1']
    ant2 = rec['antenna2']
    weight = rec['weight']   # shape (npol, nrow) -- one value per row per corr

    npol = data.shape[0]
    if npol >= 2:
        # assumes parallel hands are first/last pol index (XX,YY) -
        # double check against your correlator's pol order if unsure
        # (see rec = ms.getdata([col,'axis_info'], ifraxis=True) -> 'corr_axis')
        data_i = 0.5 * (data[0] + data[-1])
        flag_i = flag[0] | flag[-1]
        w_XX = weight[0, :]
        w_YY = weight[-1, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            var_i_row = 0.25 * (np.where(w_XX > 0, 1.0 / w_XX, np.inf) +
                                 np.where(w_YY > 0, 1.0 / w_YY, np.inf))
        w_i_row = np.where(np.isfinite(var_i_row) & (var_i_row > 0), 1.0 / var_i_row, 0.0)
    else:
        data_i = data[0]
        flag_i = flag[0]
        w_i_row = weight[0, :]

    sum_wv, sum_w, counts = {}, {}, {}
    nrow = data_i.shape[-1]
    for row in range(nrow):
        a1, a2 = int(ant1[row]), int(ant2[row])
        if a1 == a2:
            continue
        vals = data_i[:, row]
        good = ~flag_i[:, row]
        w_row = w_i_row[row]
        if not np.any(good) or w_row <= 0:
            continue
        key = (a1, a2)
        n_good = int(np.sum(good))
        v_sum = np.sum(vals[good])
        if key in sum_wv:
            sum_wv[key] += w_row * v_sum
            sum_w[key] += w_row * n_good
            counts[key] += n_good
        else:
            sum_wv[key] = w_row * v_sum
            sum_w[key] = w_row * n_good
            counts[key] = n_good

    out = {}
    for key in sum_wv:
        w_tot = sum_w[key]
        v_mean = sum_wv[key] / w_tot
        sigma_theoretical = 1.0 / np.sqrt(w_tot)   # theoretical noise on the time-averaged vis
        out[key] = (v_mean, sigma_theoretical, counts[key])
    return out


def closure_triples(baseline_vis):
    """
    Given {(a1,a2): (vis, sigma, n)}, returns a list of
    (triangle, flux_estimate_Jy, flux_err_theoretical_Jy,
     closure_phase_deg, n_min) for every closed antenna triangle.
    """
    def get_vis(a, b):
        if (a, b) in baseline_vis:
            v, s, n = baseline_vis[(a, b)]
            return v, s, n
        elif (b, a) in baseline_vis:
            v, s, n = baseline_vis[(b, a)]
            return np.conj(v), s, n
        return None, None, None

    ants = sorted(set(a for pair in baseline_vis for a in pair))
    out = []
    for i, j, k in itertools.combinations(ants, 3):
        v_ij, s_ij, n_ij = get_vis(i, j)
        v_jk, s_jk, n_jk = get_vis(j, k)
        v_ik, s_ik, n_ik = get_vis(i, k)
        if v_ij is None or v_jk is None or v_ik is None:
            continue

        a_ij, a_jk, a_ik = np.abs(v_ij), np.abs(v_jk), np.abs(v_ik)
        if a_ij == 0 or a_jk == 0 or a_ik == 0:
            continue

        bispectrum = v_ij * v_jk * np.conj(v_ik)
        flux_est = np.abs(bispectrum) ** (1.0 / 3.0)
        closure_phase = np.degrees(np.angle(bispectrum))

        # propagate fractional baseline errors through the triple product:
        # sigma_flux/flux = (1/3) * sqrt(sum of fractional variances)
        frac_err_sq = (s_ij / a_ij) ** 2 + (s_jk / a_jk) ** 2 + (s_ik / a_ik) ** 2
        flux_err_theoretical = (flux_est / 3.0) * np.sqrt(frac_err_sq)

        out.append(((i, j, k), flux_est, flux_err_theoretical, closure_phase, min(n_ij, n_jk, n_ik)))
    return out


def combine_flux_theoretical(triples):
    """
    Combines per-triangle flux estimates two ways:
    - inverse-variance-weighted mean using the theoretical (radiometer) errors
    - median + scatter-based error across triangles (empirical)
    Returns both so the ratio between them can flag noise-dominated
    vs. structure/offset-dominated targets.
    """
    fluxes = np.array([t[1] for t in triples])
    errs_theo = np.array([t[2] for t in triples])
    phases = np.array([t[3] for t in triples])

    w = 1.0 / errs_theo**2
    flux_theoretical = np.sum(w * fluxes) / np.sum(w)
    flux_err_theoretical = 1.0 / np.sqrt(np.sum(w))

    flux_empirical = np.median(fluxes)
    flux_err_empirical = np.std(fluxes) / np.sqrt(len(fluxes))

    return {
        'flux_theoretical': flux_theoretical,
        'flux_err_theoretical': flux_err_theoretical,
        'flux_empirical': flux_empirical,
        'flux_err_empirical': flux_err_empirical,
        'mean_cphase': np.mean(phases),
        'cphase_scatter': np.std(phases),
        'n_triangles': len(triples),
    }


results = []  # (field, flux_theo_Jy, flux_err_theo_Jy, flux_emp_Jy, flux_err_emp_Jy,
              #  err_ratio, mean_cphase_deg, cphase_scatter_deg, n_triangles, status)

for tgt in targets:
    try:
        baseline_vis = get_field_baseline_vis(targets_corrected_ms, tgt, datacolumn='data')
        triples = closure_triples(baseline_vis)

        if len(triples) == 0:
            results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                             float('nan'), float('nan'), float('nan'), 0, 'NO_TRIANGLES'))
            print(f"[FAIL] {tgt}: no closed antenna triangles found")
            continue

        combined = combine_flux_theoretical(triples)
        err_ratio = combined['flux_err_empirical'] / combined['flux_err_theoretical']

        status = 'OK' if len(triples) == 20 else 'OK_FEW_TRIANGLES'
        results.append((tgt, combined['flux_theoretical'], combined['flux_err_theoretical'],
                         combined['flux_empirical'], combined['flux_err_empirical'], err_ratio,
                         combined['mean_cphase'], combined['cphase_scatter'],
                         combined['n_triangles'], status))
        print(f"[{status}]   {tgt}: theo={combined['flux_theoretical']:.6f} +/- "
              f"{combined['flux_err_theoretical']:.2e} Jy, "
              f"emp={combined['flux_empirical']:.6f} +/- {combined['flux_err_empirical']:.2e} Jy "
              f"(emp/theo err ratio={err_ratio:.2f}), "
              f"closure phase = {combined['mean_cphase']:.2f} +/- {combined['cphase_scatter']:.2f} deg "
              f"({combined['n_triangles']} triangles)")

    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                         float('nan'), float('nan'), float('nan'), 0, f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: closure-phase flux extraction failed - {e}")

# ----------------------------------------------------------------------
# Write closure-phase flux results to a txt file
# ----------------------------------------------------------------------
with open(flux_results_file, 'w') as f:
    header = (f"{'field':<25}{'flux_theo_Jy':<15}{'flux_err_theo_Jy':<18}"
              f"{'flux_emp_Jy':<15}{'flux_err_emp_Jy':<17}{'err_ratio':<12}"
              f"{'mean_cphase_deg':<18}{'cphase_scatter_deg':<20}"
              f"{'n_triangles':<14}{'status':<22}\n")
    f.write(header)
    for (tgt, flux_theo, flux_err_theo, flux_emp, flux_err_emp, err_ratio,
         mean_cphase, cphase_scatter, n_tri, status) in results:
        f.write(f"{tgt:<25}{flux_theo:<15.6e}{flux_err_theo:<18.6e}"
                f"{flux_emp:<15.6e}{flux_err_emp:<17.6e}{err_ratio:<12.3f}"
                f"{mean_cphase:<18.4f}{cphase_scatter:<20.4f}"
                f"{n_tri:<14d}{status:<22}\n")

print(f"[INFO] Closure-phase flux results for {len(results)} targets written to {flux_results_file}")


#NOW we compute triple amplitude and then average over time
import itertools
import numpy as np

# ----------------------------------------------------------------------
# Flux extraction via closure triple-product (bispectrum), time-binned
# to avoid decorrelation bias from averaging visibilities across a
# window where antenna phases may have drifted. Bispectra are formed
# per time bin (where phase cancellation is exact), then the complex
# bispectra themselves are vector-averaged across bins, and only then
# is the flux (|.|^(1/3)) and error extracted.
# ----------------------------------------------------------------------
DATACOL_MAP = {'data': 'data', 'corrected': 'corrected_data', 'model': 'model_data'}


def get_field_baseline_vis_timebinned(vis, field, datacolumn='data', spw='', time_bin_sec=30.0):
    """
    Returns {tbin_idx: {(a1,a2): (weighted_mean_vis, sigma_theoretical, n)}}
    i.e. the same per-baseline weighted-mean/sigma as before, but computed
    separately within each time bin instead of across the whole track.
    time_bin_sec should be short relative to the residual atmospheric/
    antenna phase coherence time (a good starting point is your gaincal
    solint, e.g. 45-60s, or shorter if you want to test sensitivity).
    """
    col = DATACOL_MAP[datacolumn]
    ms.open(vis)
    ms.selectinit(reset=True)
    sel = {'field': field}
    if spw:
        sel['spw'] = spw
    ms.msselect(sel)

    rec = ms.getdata([col, 'flag', 'antenna1', 'antenna2', 'weight', 'time'])
    ms.close()

    data = rec[col]          # shape (npol, nchan, nrow)
    flag = rec['flag']
    ant1 = rec['antenna1']
    ant2 = rec['antenna2']
    weight = rec['weight']   # shape (npol, nrow)
    time = rec['time']

    npol = data.shape[0]
    if npol >= 2:
        # assumes parallel hands are first/last pol index (XX,YY) -
        # double check against your correlator's pol order if unsure
        data_i = 0.5 * (data[0] + data[-1])
        flag_i = flag[0] | flag[-1]
        w_XX = weight[0, :]
        w_YY = weight[-1, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            var_i_row = 0.25 * (np.where(w_XX > 0, 1.0 / w_XX, np.inf) +
                                 np.where(w_YY > 0, 1.0 / w_YY, np.inf))
        w_i_row = np.where(np.isfinite(var_i_row) & (var_i_row > 0), 1.0 / var_i_row, 0.0)
    else:
        data_i = data[0]
        flag_i = flag[0]
        w_i_row = weight[0, :]

    t0 = time.min()
    tbin_idx = np.floor((time - t0) / time_bin_sec).astype(int)

    binned = {}  # tbin -> {(a1,a2): [sum_wv, sum_w, count]}
    nrow = data_i.shape[-1]
    for row in range(nrow):
        a1, a2 = int(ant1[row]), int(ant2[row])
        if a1 == a2:
            continue
        vals = data_i[:, row]
        good = ~flag_i[:, row]
        w_row = w_i_row[row]
        if not np.any(good) or w_row <= 0:
            continue

        tb = int(tbin_idx[row])
        key = (a1, a2)
        n_good = int(np.sum(good))
        v_sum = np.sum(vals[good])

        bin_dict = binned.setdefault(tb, {})
        if key in bin_dict:
            bin_dict[key][0] += w_row * v_sum
            bin_dict[key][1] += w_row * n_good
            bin_dict[key][2] += n_good
        else:
            bin_dict[key] = [w_row * v_sum, w_row * n_good, n_good]

    out = {}
    for tb, bl in binned.items():
        out[tb] = {}
        for key, (sum_wv, sum_w, n) in bl.items():
            v_mean = sum_wv / sum_w
            sigma_theoretical = 1.0 / np.sqrt(sum_w)
            out[tb][key] = (v_mean, sigma_theoretical, n)
    return out


def bispectra_per_bin(baseline_vis_bin):
    """
    Given one time bin's {(a1,a2): (vis, sigma, n)}, returns
    {triangle: (B_complex, sigma_B, n_min)} for every closed triangle
    present in this bin. This is where the phase cancellation happens
    exactly, before any time-averaging occurs.
    """
    def get_vis(a, b):
        if (a, b) in baseline_vis_bin:
            v, s, n = baseline_vis_bin[(a, b)]
            return v, s, n
        elif (b, a) in baseline_vis_bin:
            v, s, n = baseline_vis_bin[(b, a)]
            return np.conj(v), s, n
        return None, None, None

    ants = sorted(set(a for pair in baseline_vis_bin for a in pair))
    out = {}
    for i, j, k in itertools.combinations(ants, 3):
        v_ij, s_ij, n_ij = get_vis(i, j)
        v_jk, s_jk, n_jk = get_vis(j, k)
        v_ik, s_ik, n_ik = get_vis(i, k)
        if v_ij is None or v_jk is None or v_ik is None:
            continue

        a_ij, a_jk, a_ik = np.abs(v_ij), np.abs(v_jk), np.abs(v_ik)
        if a_ij == 0 or a_jk == 0 or a_ik == 0:
            continue

        B = v_ij * v_jk * np.conj(v_ik)

        frac_err_sq = (s_ij / a_ij) ** 2 + (s_jk / a_jk) ** 2 + (s_ik / a_ik) ** 2
        sigma_B = np.abs(B) * np.sqrt(frac_err_sq)   # approx magnitude of complex bispectrum error

        out[(i, j, k)] = (B, sigma_B, min(n_ij, n_jk, n_ik))
    return out


def closure_triples_timeavg(binned_baseline_vis):
    """
    Vector-averages the complex bispectrum for each triangle across all
    time bins where it's present (inverse-variance weighted), THEN takes
    |.|^(1/3) for flux and angle() for closure phase. This preserves the
    exact per-instant phase cancellation instead of averaging visibilities
    (and therefore antenna phases) across the whole track first.
    Returns list of (triangle, flux_Jy, flux_err_theoretical_Jy, closure_phase_deg, n_bins).
    """
    per_triangle = {}  # triangle -> [sum_wB, sum_w, n_bins_used]

    for tb, baseline_vis_bin in binned_baseline_vis.items():
        bispectra = bispectra_per_bin(baseline_vis_bin)
        for tri, (B, sigma_B, n_min) in bispectra.items():
            if sigma_B <= 0 or not np.isfinite(sigma_B):
                continue
            w = 1.0 / sigma_B**2
            if tri in per_triangle:
                per_triangle[tri][0] += w * B
                per_triangle[tri][1] += w
                per_triangle[tri][2] += 1
            else:
                per_triangle[tri] = [w * B, w, 1]

    out = []
    for tri, (sum_wB, sum_w, n_bins) in per_triangle.items():
        B_avg = sum_wB / sum_w
        sigma_B_avg = 1.0 / np.sqrt(sum_w)

        flux_est = np.abs(B_avg) ** (1.0 / 3.0)
        closure_phase = np.degrees(np.angle(B_avg))

        frac_err = sigma_B_avg / np.abs(B_avg)
        flux_err_theoretical = (flux_est / 3.0) * frac_err

        out.append((tri, flux_est, flux_err_theoretical, closure_phase, n_bins))
    return out


def combine_flux_theoretical(triples):
    """
    Combines per-triangle flux estimates two ways:
    - inverse-variance-weighted mean using the theoretical (radiometer) errors
    - median + scatter-based error across triangles (empirical)
    """
    fluxes = np.array([t[1] for t in triples])
    errs_theo = np.array([t[2] for t in triples])
    phases = np.array([t[3] for t in triples])

    w = 1.0 / errs_theo**2
    flux_theoretical = np.sum(w * fluxes) / np.sum(w)
    flux_err_theoretical = 1.0 / np.sqrt(np.sum(w))

    flux_empirical = np.median(fluxes)
    flux_err_empirical = np.std(fluxes) / np.sqrt(len(fluxes))

    return {
        'flux_theoretical': flux_theoretical,
        'flux_err_theoretical': flux_err_theoretical,
        'flux_empirical': flux_empirical,
        'flux_err_empirical': flux_err_empirical,
        'mean_cphase': np.mean(phases),
        'cphase_scatter': np.std(phases),
        'n_triangles': len(triples),
    }


TIME_BIN_SEC = 30.0  # tune relative to residual phase coherence time

results = []  # (field, flux_theo_Jy, flux_err_theo_Jy, flux_emp_Jy, flux_err_emp_Jy,
              #  err_ratio, mean_cphase_deg, cphase_scatter_deg, n_triangles, status)

for tgt in targets:
    try:
        binned_baseline_vis = get_field_baseline_vis_timebinned(
            targets_corrected_ms, tgt, datacolumn='data', time_bin_sec=TIME_BIN_SEC)
        triples = closure_triples_timeavg(binned_baseline_vis)

        if len(triples) == 0:
            results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                             float('nan'), float('nan'), float('nan'), 0, 'NO_TRIANGLES'))
            print(f"[FAIL] {tgt}: no closed antenna triangles found")
            continue

        combined = combine_flux_theoretical(triples)
        err_ratio = combined['flux_err_empirical'] / combined['flux_err_theoretical']

        status = 'OK' if len(triples) == 20 else 'OK_FEW_TRIANGLES'
        results.append((tgt, combined['flux_theoretical'], combined['flux_err_theoretical'],
                         combined['flux_empirical'], combined['flux_err_empirical'], err_ratio,
                         combined['mean_cphase'], combined['cphase_scatter'],
                         combined['n_triangles'], status))
        print(f"[{status}]   {tgt}: theo={combined['flux_theoretical']:.6f} +/- "
              f"{combined['flux_err_theoretical']:.2e} Jy, "
              f"emp={combined['flux_empirical']:.6f} +/- {combined['flux_err_empirical']:.2e} Jy "
              f"(emp/theo err ratio={err_ratio:.2f}), "
              f"closure phase = {combined['mean_cphase']:.2f} +/- {combined['cphase_scatter']:.2f} deg "
              f"({combined['n_triangles']} triangles)")

    except Exception as e:
        results.append((tgt, float('nan'), float('nan'), float('nan'), float('nan'),
                         float('nan'), float('nan'), float('nan'), 0, f'FIT_FAIL: {e}'))
        print(f"[FAIL] {tgt}: closure-phase flux extraction failed - {e}")

# ----------------------------------------------------------------------
# Write closure-phase flux results to a txt file
# ----------------------------------------------------------------------
with open(flux_results_file, 'w') as f:
    header = (f"{'field':<25}{'flux_theo_Jy':<15}{'flux_err_theo_Jy':<18}"
              f"{'flux_emp_Jy':<15}{'flux_err_emp_Jy':<17}{'err_ratio':<12}"
              f"{'mean_cphase_deg':<18}{'cphase_scatter_deg':<20}"
              f"{'n_triangles':<14}{'status':<22}\n")
    f.write(header)
    for (tgt, flux_theo, flux_err_theo, flux_emp, flux_err_emp, err_ratio,
         mean_cphase, cphase_scatter, n_tri, status) in results:
        f.write(f"{tgt:<25}{flux_theo:<15.6e}{flux_err_theo:<18.6e}"
                f"{flux_emp:<15.6e}{flux_err_emp:<17.6e}{err_ratio:<12.3f}"
                f"{mean_cphase:<18.4f}{cphase_scatter:<20.4f}"
                f"{n_tri:<14d}{status:<22}\n")

print(f"[INFO] Closure-phase flux results for {len(results)} targets written to {flux_results_file}")
'''

#Code to compare fluxes 
import pandas as pd
import matplotlib.pyplot as plt

widths_uvmodelfit  = [25, 15, 15, 12, 12, 12, 12, 12, 22]
widths_phaseclosure = [25, 15, 18, 15, 17, 12, 18, 20, 14, 22]

df_phaseclosure = pd.read_fwf('target_fluxes_2026-06-25_phaseclosure.txt', widths=widths_phaseclosure)
df_phaseclosure.columns = [c.strip() for c in df_phaseclosure.columns]
df_phaseclosure['field'] = df_phaseclosure['field'].str.strip()

df_uvmodelfit = pd.read_fwf('target_fluxes_2026-06-25_uvmodelfit.txt', widths=widths_uvmodelfit)
df_uvmodelfit.columns = [c.strip() for c in df_uvmodelfit.columns]
df_uvmodelfit['field'] = df_uvmodelfit['field'].str.strip()

# merge on field name so rows actually correspond to the same target,
# and drop any target that failed/has no valid flux in either method
merged = pd.merge(df_uvmodelfit, df_phaseclosure, on='field', suffixes=('_uvmodelfit', '_phaseclosure'))
merged = merged.dropna(subset=['flux_Jy', 'flux_theo_Jy'])

flux_uvmodelfit = merged['flux_Jy']
flux_phaseclosure  = merged['flux_theo_Jy']

plt.figure(figsize=(8, 6))
plt.scatter(flux_uvmodelfit, flux_phaseclosure, s=20, label='Flux fits')

lims = [min(flux_uvmodelfit.min(), flux_phaseclosure.min()), max(flux_uvmodelfit.max(), flux_phaseclosure.max())]
plt.plot(lims, lims, 'k--', label='y=x')

plt.xlabel('Flux from uvmodelfit (Jy)')
plt.ylabel('Flux from phase closure(Jy)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.title('Flux comparison for >3sigma detected sources')
plt.savefig('flux_uvmodelfit_vs_phaseclosure.png')
