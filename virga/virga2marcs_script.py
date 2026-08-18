import numpy as np
import pandas as pd
import astropy.units as u
from molmass import Formula

import virga.justdoit as jdi

import h5py

k_B = 1.380649e-16               #ergs/K


def _read_header_value(lines, tag):
    for l in lines:
        if tag in l:
            return float(l.strip().split(tag)[-1].strip())
    raise ValueError(f"could not find '{tag}' in marcs2virga.dat header")


def lognormal(N, rg, sig, r_array):
    # Ackerman and Marley 2001 lognormal distribution
    # Equation 9, can't find a later reference in the literature
    # units are dn/dr so cm^-3cm^-1
    # N and rg are shape (nz, ncond), sig is a scalar, and r_array is shape (nbins,)
    # want final output to be shape (nz, cond, nbins) so that we can integrate over
    # the size distribution for each layer and then the condensate species separately
    prefactor = N[:,np.newaxis] / (r_array[np.newaxis, np.newaxis, :] * np.log(sig) * np.sqrt(2 * np.pi))
    exponent = - (np.log(r_array[np.newaxis, np.newaxis, :] / rg[:,np.newaxis]) / (np.sqrt(2) * np.log(sig)))**2
    distribution = prefactor * np.exp(exponent)
    return distribution


def lognormal_abs_sca_sum_matrix(r_array, distribution, Qabs, Qsca):
    # Integrate over radius bins first (per condensate), then sum condensates.
    # distribution: (nz, ncond, nbins)
    # Qabs, Qsca expected: (ncond, nbins, nwavelengths)
    # r_array: (nbins,)
    nz, ncond, nbins = distribution.shape
    nwavelengths = Qabs.shape[-1]

    # Build per-bin widths dr (same length as r_array)
    dr = np.empty_like(r_array)
    dr[1:-1] = 0.5 * (r_array[2:] - r_array[:-2])
    dr[0] = r_array[1] - r_array[0]
    dr[-1] = r_array[-1] - r_array[-2]

    # Geometric factor (nz, ncond, nbins)
    factor = distribution * np.pi * r_array**2 * dr

    factor_2D = factor.reshape(nz, ncond*nbins)
    Qabs_2D = Qabs.reshape(ncond*nbins, nwavelengths)
    Qsca_2D = Qsca.reshape(ncond*nbins, nwavelengths)

    kappa_abs = factor_2D @ Qabs_2D  # shape (nz, nwavelengths)
    kappa_sca = factor_2D @ Qsca_2D  # shape (nz, nwavelengths)

    return kappa_abs, kappa_sca


def read_in_mie_h5(file_path):
    """Load precomputed Mie data from HDF5 file."""
    with h5py.File(file_path, "r") as f:
        mie_radii = f["sizes_um"][:]
        wavelengths_um = f["wavelengths_um"][:]
        qabs = f["qabs"][:]
        qsca = f["qsca"][:]

    # Ensure orientation is (nbins, nwavelength)
    if qabs.shape[0] != mie_radii.size and qabs.shape[1] == mie_radii.size:
        qabs = qabs.T
    if qsca.shape[0] != mie_radii.size and qsca.shape[1] == mie_radii.size:
        qsca = qsca.T

    return mie_radii, wavelengths_um, qabs, qsca


def run_virga(input_dir='.', mieff_dir=None, mie_data_dir=None):
    """
    Read marcs2virga.dat + marcs_wnos.dat (both written by marcs.f's
    OSTABLOOK), run virga, and write virga2marcs.dat in the exact format
    and on the exact wavenumber grid marcs.f reads back: one header line,
    then "wavenumber(1/cm) pressure(dyn/cm^2) kappa_abs(cm^2/g)
    kappa_sca(cm^2/g)" rows, nwtot wavelengths x ntau depths, ordered
    wavelength-major (all depths for wavelength 1, then all depths for
    wavelength 2, ...).

    Paths are relative to `input_dir` so this can be run from the SLURM
    job's scratch directory where marcs.f writes its files.
    """
    input_dir = input_dir.rstrip('/') or '.'
    if mieff_dir is None:
        mieff_dir = f'{input_dir}/virga/'
    if mie_data_dir is None:
        mie_data_dir = f'{input_dir}/Mie_data'

    input_data_marcs = f'{input_dir}/marcs2virga.dat'
    wnos_path = f'{input_dir}/marcs_wnos.dat'

    condensates_start = ' # Begin Condensates'
    condensates_end = ' # End Condensates'

    with open(input_data_marcs, 'r') as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if condensates_start in l)
    end = next(i for i, l in enumerate(lines) if condensates_end in l)
    block = lines[start+1:end]

    condensates = {}
    for line in block:
        parts = line.strip().lstrip('#').strip().split(':')
        molecule_name = parts[0].strip()
        condensates[molecule_name] = float(parts[1].strip())

    molecules = list(condensates.keys())
    # marcs.f writes each condensate's partial pressure (dyn/cm^2) at the
    # bottom/reference layer (ppallmol(ntau, ...)), not a mass mixing
    # ratio. Converted to the dimensionless ext_mmr virga expects below,
    # once P (bottom-layer total pressure) and mean_molecular_weight are
    # available.
    partial_pressures = np.array(list(condensates.values()))

    mean_molecular_weight = _read_header_value(lines[:start], 'mean molecular weight')
    metallicity = _read_header_value(lines[:start], 'metallicity')
    # marcs.f writes GRAV from COMMON /CG/ in cm/s^2 (cgs), not m/s^2.
    gravity_cgs = _read_header_value(lines[:start], 'surface gravity')

    T = np.genfromtxt(input_data_marcs, skip_header=end+2)[:,0]             #K
    P = np.genfromtxt(input_data_marcs, skip_header=end+2)[:,1]*1e-6        #bar
    kz = np.genfromtxt(input_data_marcs, skip_header=end+2)[:,2]            #cm^2/s
    ntau = len(P)

    # mole fraction -> mass mixing ratio, at the same bottom layer (ntau)
    # marcs.f sampled the partial pressure from:
    #   mmr = (pp / Ptot_bottom) * (gas molar mass / mean molecular weight)
    Ptot_bottom = P[-1] * 1e6   # bar -> dyn/cm^2
    ext_mmr = np.array([
        (pp / Ptot_bottom) * (Formula(name).mass / mean_molecular_weight)
        for name, pp in zip(molecules, partial_pressures)
    ])

    ############################# SETTING UP VIRGA #################################
    profile = pd.DataFrame({'pressure': P, 'temperature': T, 'kz': kz})

    a = jdi.Atmosphere(molecules,
                    fsed=0.1, mh=metallicity,
                    mmw=mean_molecular_weight)

    a.gravity(gravity=gravity_cgs, gravity_unit=u.Unit('cm/(s**2)'))

    a.ptk(df=profile)

    all_out = jdi.compute(a, as_dict=True,
                        directory=mieff_dir, ext_mmr=ext_mmr)

    # virga treats the ntau (T,P) points given as layer *boundaries* and
    # returns one fewer per-layer quantity (confirmed empirically: a
    # dry-run with ntau=53 input rows got back 52-row output arrays).
    # Derive n_layers from virga's own output rather than assuming
    # ntau-1 always holds.
    n_layers = all_out['mean_particle_r'].shape[0]

    pp = np.zeros((n_layers, len(molecules)))
    n_c = np.zeros((n_layers, len(molecules)))
    rg = np.zeros((n_layers, len(molecules)))
    for i, igas in enumerate(molecules):
        mass_igas = Formula(igas).mass                      #amu
        cgs_mass_H = Formula("H").mass*u.u.to(u.g)          #g

        mmr_c = all_out['condensate_mmr'][:,i]
        mmr_tot = all_out['cond_plus_gas_mmr'][:,i]
        rho_p = all_out['condensate_density'][i]            # g/cm^3
        rg[:,i] = all_out['mean_particle_r'][:,i]*1e-4       #microns -> cm

        for j in range(n_layers):
            pp[j][i] = (mean_molecular_weight/mass_igas) * P[j]*1e6 * (mmr_tot[j] - mmr_c[j])         #dyn/cm^2
            if rg[j][i] > 0:
                n_c[j][i] = (3*mmr_c[j]*mean_molecular_weight*cgs_mass_H*(P[j]*1e6))/(4*np.pi*rg[j][i]**3*rho_p*k_B*T[j])   #number density

    # Diagnostic side-product only -- marcs.f does not read this file.
    particle_props_path = f'{input_dir}/virga_particle_props.dat'
    with open(particle_props_path, 'w') as file:
        file.write("P (bar) r (cm) pp (dyn/cm^2) n_c (1/cm^3)\n")
        for j in range(len(molecules)):
            file.write(f"\nGas: {molecules[j]}\n")
            for i in range(len(P)):
                if i < n_layers:
                    file.write(f"{P[i]:.6e} {rg[i][j]:.6e} {pp[i][j]:.6e} {n_c[i][j]:.6e}\n")
                else:
                    file.write(f"{P[i]:.6e} {0} {0} {0}\n")

    Qabs_all = Qsca_all = mie_radii_cm = wavelength = None
    for i, igas in enumerate(molecules):
        h5_path = f"{mie_data_dir}/Mie_{igas}.h5"
        mie_radii, wavelength, Qabs, Qsca = read_in_mie_h5(h5_path) #TODO need to actually loop over the materials
        mie_radii_cm = mie_radii * 1e-4  # microns -> cm, for the distribution function

        if i == 0:
            Qabs_all = np.zeros((len(molecules), Qabs.shape[0], Qabs.shape[1]))
            Qsca_all = np.zeros((len(molecules), Qsca.shape[0], Qsca.shape[1]))

        Qabs_all[i, :, :] = Qabs
        Qsca_all[i, :, :] = Qsca

    sig = list(all_out['scalar_inputs'].values())[2]  # Geometric standard deviation, is a scalar

    distribution = lognormal(n_c, rg, sig, mie_radii_cm)
    k_abs, k_sca = lognormal_abs_sca_sum_matrix(mie_radii_cm, distribution, Qabs_all, Qsca_all)

    # Convert from cm^-1 (extinction per length) to cm^2/g (mass opacity)
    # by dividing by the gas density: rho_gas = P/(k_B*T) * mmw * m_H,
    # with P in dyn/cm^2, T in K, mmw dimensionless, m_H in g.
    rho_gas = (P[:n_layers]*1e6) / (k_B * T[:n_layers]) * mean_molecular_weight * Formula("H").mass*u.u.to(u.g)  # g/cm^3
    kappa_abs_cm2_per_g = k_abs / rho_gas[:, np.newaxis]
    kappa_sca_cm2_per_g = k_sca / rho_gas[:, np.newaxis]

    # virga's own Mie-table wavenumber grid (cm^-1) -- NOT the same grid as
    # MARCS's own opacity-sampling wavenumbers (WNOS), which is why this
    # gets interpolated onto marcs_wnos.dat below rather than written out
    # directly.
    mie_wavenumber = 1e4 / wavelength   # microns -> cm^-1

    with open(wnos_path) as f:
        nwtot = int(f.readline())
    marcs_wnos = np.genfromtxt(wnos_path, skip_header=1)
    if marcs_wnos.size != nwtot:
        raise ValueError(f"{wnos_path} header says {nwtot} wavelengths "
                          f"but contains {marcs_wnos.size}")

    sort_idx = np.argsort(mie_wavenumber)
    log_mie_wn = np.log(mie_wavenumber[sort_idx])
    log_marcs_wn = np.log(marcs_wnos)

    kappa_abs_regridded = np.empty((n_layers, nwtot))
    kappa_sca_regridded = np.empty((n_layers, nwtot))
    for ip in range(n_layers):
        kappa_abs_regridded[ip, :] = np.exp(np.interp(
            log_marcs_wn, log_mie_wn,
            np.log(np.clip(kappa_abs_cm2_per_g[ip, sort_idx], 1e-99, None))))
        kappa_sca_regridded[ip, :] = np.exp(np.interp(
            log_marcs_wn, log_mie_wn,
            np.log(np.clip(kappa_sca_cm2_per_g[ip, sort_idx], 1e-99, None))))

    # marcs.f's OSTABLOOK expects exactly ntau depths per wavelength
    # (it wrote T(k)/Ptot(k)/K_zz(k) for k=1,ntau into marcs2virga.dat).
    # virga gave back n_layers=ntau-1 layer-averaged values, so pad the
    # deepest layer by repeating the last one -- an approximation, not a
    # physically distinct extra layer.
    if n_layers == ntau - 1:
        kappa_abs_out = np.vstack([kappa_abs_regridded, kappa_abs_regridded[-1:]])
        kappa_sca_out = np.vstack([kappa_sca_regridded, kappa_sca_regridded[-1:]])
    elif n_layers == ntau:
        kappa_abs_out = kappa_abs_regridded
        kappa_sca_out = kappa_sca_regridded
    else:
        raise ValueError(f"virga returned {n_layers} layers, expected "
                          f"{ntau} or {ntau-1} for ntau={ntau} input rows")

    output_path = f'{input_dir}/virga2marcs.dat'
    with open(output_path, 'w') as file:
        file.write("Wavenumber (1/cm) Pressure (dyn/cm^2) kappa_abs (cm^2/g) kappa_sca (cm^2/g)\n")
        for iw in range(nwtot):
            waveno = marcs_wnos[iw]
            for ip in range(ntau):
                file.write(f"{waveno:.6e} {P[ip]*1e6:.6e} "
                           f"{kappa_abs_out[ip, iw]:.6e} "
                           f"{kappa_sca_out[ip, iw]:.6e}\n")

    return output_path


if __name__ == '__main__':
    run_virga()
