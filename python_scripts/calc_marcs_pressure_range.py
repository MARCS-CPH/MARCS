"""
Given a finished MARCS run, predict the range of log(tau_ross) needed so
that the model's gas-pressure structure spans a desired pressure range.

Reads the "M O D E L  A T M O S P H E R E" table out of marcs_output.dat
(columns: K, TAUROSS, TAU(5000), GEOM.DEPTH, T, PE, PG, PRAD, PTURB,
KAPPAROSS -- PG is the total gas pressure in dyn/cm^2 = barye, KAPPAROSS
is the Rosseland mean opacity in cm^2/g), converts PG to bar, and finds
log(tau_ross) at the requested min/max pressure using one of two modes:

  table         Interpolate within the computed grid, or linearly
                extrapolate in log(Pg)-log(tau) space using the nearest
                edge points if the target falls outside the grid. Purely
                data-driven, no physics assumed beyond a local power law.

  hydrostatic   Integrate the hydrostatic-equilibrium relation
                dPg/dtau_ross = g/kappa_ross (g from the model's log(g),
                kappa_ross from the tabulated Rosseland mean opacity)
                starting from the top of the computed grid. Within the
                grid this gives an independent, physics-based cross-check
                of the tabulated Pg; beyond the grid it extrapolates by
                assuming kappa_ross stays at its edge value, which is a
                different (and often better motivated) assumption than
                the table mode's power-law fit.

  both          Run and report both modes side by side.

For now the target range is fixed to 1 microbar -- 1 bar, matching the
convention used for the K2-18b models, but it can be overridden with
--min-pressure / --max-pressure (in bar).

Usage:
    python calc_marcs_pressure_range.py [model_path] [--min-pressure BAR] [--max-pressure BAR] [--mode {table,hydrostatic,both}]

model_path may be either a marcs_output.dat file or the run directory
containing one. Defaults to the z_200 K2-18b run.
"""

import argparse
import os
import re
import sys

import numpy as np

DEFAULT_MODEL = (
    "/groups/astro/tbalduin/MARCS_CPH/MARCS/models/"
    "K2_18b_highZ_CNOS/chain/z_200/marcs_output.dat"
)

DEFAULT_MODEL = ("/groups/astro/tbalduin/MARCS_CPH/MARCS/models/K2_18b_highZ_CNOS/chain/z_200_v6_DMS_noneq_off_virga/marcs_output.dat")

DYN_CM2_PER_BAR = 1.0e6  # 1 bar = 1e6 dyn/cm^2 (Ba); 1 microbar = 1 dyn/cm^2

# column indices in the "M O D E L  A T M O S P H E R E" table
COL_TAUROSS = 1
COL_PG = 6
COL_KAPPAROSS = 9


def resolve_model_path(path):
    if os.path.isdir(path):
        path = os.path.join(path, "marcs_output.dat")
    if not os.path.isfile(path):
        sys.exit(f"error: could not find a MARCS output file at '{path}'")
    return path


def parse_logg(path):
    """Return log(g) [cgs] from the model header, e.g. '... log(g) = 3.0'."""
    with open(path) as f:
        for _, line in zip(range(20), f):
            m = re.search(r"log\(g\)\s*=\s*([-+0-9.eE]+)", line)
            if m:
                return float(m.group(1))
    sys.exit(f"error: could not find 'log(g) = ...' in the header of {path}")


def parse_model_atmosphere(path):
    """Return (logtau, pg_dyncm2, kappa_ross) from the model-atmosphere table."""
    with open(path) as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if "M O D E L" in line and "A T M O S P H E R E" in line:
            header_idx = i
            break
    if header_idx is None:
        sys.exit(f"error: no 'M O D E L  A T M O S P H E R E' table found in {path}")

    # the column-name row (contains "TAUROSS" and "PG") follows shortly after
    col_row = None
    for i in range(header_idx, min(header_idx + 5, len(lines))):
        if "TAUROSS" in lines[i] and "PG" in lines[i]:
            col_row = i
            break
    if col_row is None:
        sys.exit(f"error: could not locate the column-header row in {path}")

    logtau, pg, kappa = [], [], []
    for line in lines[col_row + 1:]:
        if not line.strip():
            break
        parts = line.split()
        try:
            tauross = float(parts[COL_TAUROSS])
            pg_val = float(parts[COL_PG])
            kappa_val = float(parts[COL_KAPPAROSS])
        except (ValueError, IndexError):
            break
        logtau.append(np.log10(tauross))
        pg.append(pg_val)
        kappa.append(kappa_val)

    if len(logtau) < 2:
        sys.exit(f"error: fewer than 2 data rows parsed from {path}")

    return np.array(logtau), np.array(pg), np.array(kappa)


def logtau_at_pressure_table(logtau, log_pg, target_log_pg, n_extrap=3):
    """
    Find log(tau_ross) at a given log10(PG) by interpolating within the
    computed grid, or linearly extrapolating (in log-log space) using the
    nearest n_extrap points if the target falls outside the grid.

    Returns (logtau_value, extrapolated: bool, dex_beyond: float)
    """
    lo, hi = log_pg[0], log_pg[-1]
    within = min(lo, hi) <= target_log_pg <= max(lo, hi)

    if within:
        # np.interp needs increasing xp; log_pg should already increase
        # with depth (== with logtau), but sort defensively just in case.
        order = np.argsort(log_pg)
        value = np.interp(target_log_pg, log_pg[order], logtau[order])
        return value, False, 0.0

    n = min(n_extrap, len(logtau))
    if target_log_pg < lo:
        lt_fit, lp_fit = logtau[:n], log_pg[:n]
        dex_beyond = lo - target_log_pg
    else:
        lt_fit, lp_fit = logtau[-n:], log_pg[-n:]
        dex_beyond = target_log_pg - hi

    # fit logtau as a linear function of log_pg over the edge points, then
    # evaluate it at the target pressure
    slope, intercept = np.polyfit(lp_fit, lt_fit, 1)
    value = slope * target_log_pg + intercept
    return value, True, dex_beyond


def integrate_hydrostatic(logtau, pg_dyncm2, kappa_ross, logg):
    """
    Integrate dPg/dtau_ross = g/kappa_ross from the top of the grid,
    anchored to the tabulated Pg there. Returns P_hydro in dyn/cm^2 on
    the same tau grid, for comparison against the tabulated Pg.
    """
    g = 10.0 ** logg
    tau = 10.0 ** logtau
    integrand = g / kappa_ross
    increments = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(tau)
    return pg_dyncm2[0] + np.concatenate(([0.0], np.cumsum(increments)))


def logtau_at_pressure_hydrostatic(logtau, kappa_ross, p_hydro_dyncm2, target_bar,
                                    logg, n_edge=3):
    """
    Find log(tau_ross) at a given pressure (in bar) using the hydrostatic
    P_hydro(tau) curve for interpolation, or by extending dPg/dtau = g/kappa_ross
    past the grid edge (assuming kappa_ross ~ its mean edge value) for
    extrapolation.

    Returns (logtau_value, extrapolated: bool, dex_beyond: float)
    """
    g = 10.0 ** logg
    tau = 10.0 ** logtau
    target_dyncm2 = target_bar * DYN_CM2_PER_BAR
    lo, hi = p_hydro_dyncm2[0], p_hydro_dyncm2[-1]
    within = min(lo, hi) <= target_dyncm2 <= max(lo, hi)

    if within:
        order = np.argsort(p_hydro_dyncm2)
        tau_value = np.interp(target_dyncm2, p_hydro_dyncm2[order], tau[order])
        return np.log10(tau_value), False, 0.0

    n = min(n_edge, len(kappa_ross))
    if target_dyncm2 < lo:
        kappa_edge = np.mean(kappa_ross[:n])
        tau_target = tau[0] + (target_dyncm2 - lo) * kappa_edge / g
        dex_beyond = np.log10(lo) - np.log10(target_dyncm2)
    else:
        kappa_edge = np.mean(kappa_ross[-n:])
        tau_target = tau[-1] + (target_dyncm2 - hi) * kappa_edge / g
        dex_beyond = np.log10(target_dyncm2) - np.log10(hi)

    if tau_target <= 0:
        sys.exit("error: hydrostatic extrapolation to the requested pressure gives a "
                  "non-positive tau_ross -- the target pressure is too far outside the "
                  "computed range for a local constant-opacity extrapolation to be meaningful.")

    return np.log10(tau_target), True, dex_beyond


def report(label, target_bar, logtau_val, extrap, dex):
    tag = f"EXTRAPOLATED ({dex:.2f} dex beyond computed range)" if extrap else "interpolated"
    print(f"  {label} (Pg = {target_bar:.3e} bar): logtau = {logtau_val:7.3f}   [{tag}]")


def format_logtau_cards(logtau_min, logtau_max, differ=0.1):
    """
    Format the two boundary LOGTAU cards exactly as MARCS's TAUSCA reader
    expects (marcs.f, FORMAT 50 = 5(7X,F8.0)): 'LOGTAU=' occupies the first
    7 columns, the value occupies the next 8, and -- only on the first line
    -- 'DIFFER=' must start at column 16 for the step to be read correctly.
    The value's justification within its 8-column field doesn't matter to
    the Fortran reader (it parses the field for a number, not by position),
    but the field width and the DIFFER= column position do.
    """
    val_min = f"{logtau_min:.2f}"
    if len(val_min) > 8:
        raise ValueError(f"logtau_min={logtau_min} does not fit an F8.0 field")
    val_max = f"{logtau_max:.2f}"
    if len(val_max) > 8:
        raise ValueError(f"logtau_max={logtau_max} does not fit an F8.0 field")
    line1 = f"LOGTAU={val_min:<8}DIFFER={differ:g}"
    line2 = f"LOGTAU={val_max}"
    return line1, line2


def patch_input_file(template_path, out_path, logtau_min, logtau_max, differ=0.1):
    """
    Copy template_path to out_path with the first two 'LOGTAU=' cards (the
    pair TAUSCA actually reads at startup) replaced by the recommended
    range. These .input decks conventionally carry one more inert
    'LOGTAU=...' line further down (dead leftover text, never reached by
    TAUSCA's sequential read) -- that one is intentionally left untouched.
    Refuses to run if the leading pair can't be found, so it never silently
    edits the wrong lines.
    """
    line1, line2 = format_logtau_cards(logtau_min, logtau_max, differ)
    with open(template_path) as f:
        lines = f.readlines()
    replaced = 0
    for i, line in enumerate(lines):
        if line.startswith("LOGTAU="):
            lines[i] = (line1 if replaced == 0 else line2) + "\n"
            replaced += 1
            if replaced == 2:
                break
    if replaced != 2:
        sys.exit(f"error: expected 2 leading LOGTAU= cards in {template_path}, "
                  f"found {replaced} -- refusing to guess which to replace")
    if os.path.abspath(out_path) == os.path.abspath(template_path):
        print(f"warning: overwriting {out_path} in place", file=sys.stderr)
    with open(out_path, "w") as f:
        f.writelines(lines)
    print(f"Wrote {out_path} with:\n  {line1}\n  {line2}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model_path", nargs="?", default=DEFAULT_MODEL,
                         help="marcs_output.dat file, or the run directory containing one "
                              f"(default: {DEFAULT_MODEL})")
    parser.add_argument("--min-pressure", type=float, default=1.0e-6,
                         help="target minimum gas pressure in bar (default: 1e-6 = 1 microbar)")
    parser.add_argument("--max-pressure", type=float, default=1.0,
                         help="target maximum gas pressure in bar (default: 1 bar)")
    parser.add_argument("--mode", choices=["table", "hydrostatic", "both"], default="table",
                         help="'table': interpolate/extrapolate the tabulated Pg curve (default). "
                              "'hydrostatic': integrate dPg/dtau = g/kappa_ross instead. "
                              "'both': report both.")
    parser.add_argument("--emit-cards", action="store_true",
                         help="Also print the LOGTAU/DIFFER cards ready to paste into a .input file")
    parser.add_argument("--write-input", metavar="TEMPLATE:OUTPUT",
                         help="Copy TEMPLATE .input file to OUTPUT with its two leading LOGTAU "
                              "cards replaced by the recommended range (e.g. "
                              "old.input:new.input). Never overwrites TEMPLATE unless OUTPUT "
                              "is the same path.")
    parser.add_argument("--emit-mode", choices=["table", "hydrostatic"], default="hydrostatic",
                         help="Which recommendation to use for --emit-cards/--write-input "
                              "(default: hydrostatic, since that's the same g/kappa_ross "
                              "relation TRYCK/SOLVE actually enforce as the pressure boundary "
                              "condition -- see marcs.f:10198)")
    parser.add_argument("--differ", type=float, default=0.1,
                         help="Grid step in dex to record on the emitted LOGTAU card "
                              "(default: 0.1, matching the existing K2-18b templates)")
    args = parser.parse_args()

    path = resolve_model_path(args.model_path)
    logtau, pg_dyncm2, kappa_ross = parse_model_atmosphere(path)
    pg_bar = pg_dyncm2 / DYN_CM2_PER_BAR
    log_pg_bar = np.log10(pg_bar)

    print(f"Model file:            {path}")
    print(f"Computed grid:          logtau = {logtau[0]:.2f} to {logtau[-1]:.2f} "
          f"({len(logtau)} layers)")
    print(f"Computed Pg range:      {pg_bar[0]:.3e} to {pg_bar[-1]:.3e} bar")
    print()
    print(f"Target pressure range:  {args.min_pressure:.3e} bar  to  {args.max_pressure:.3e} bar")
    print()

    results = {}
    emitting = args.emit_cards or args.write_input
    need_table = args.mode in ("table", "both") or (emitting and args.emit_mode == "table")
    need_hydro = args.mode in ("hydrostatic", "both") or (emitting and args.emit_mode == "hydrostatic")

    report_table = args.mode in ("table", "both")
    report_hydro = args.mode in ("hydrostatic", "both")

    if need_table:
        if report_table:
            print("-- table mode (interpolation / log-log extrapolation of tabulated Pg) --")
        logtau_min, extrap_min, dex_min = logtau_at_pressure_table(
            logtau, log_pg_bar, np.log10(args.min_pressure))
        logtau_max, extrap_max, dex_max = logtau_at_pressure_table(
            logtau, log_pg_bar, np.log10(args.max_pressure))
        if report_table:
            report("minimum logtau needed", args.min_pressure, logtau_min, extrap_min, dex_min)
            report("maximum logtau needed", args.max_pressure, logtau_max, extrap_max, dex_max)
            print(f"  -> recommended logtau range: {logtau_min:.2f}  to  {logtau_max:.2f}")
            if extrap_min or extrap_max:
                print("  Note: extrapolated values are a local power-law (log-log) projection")
                print("  from the edge of the computed model -- treat as a starting guess.")
            print()
        results["table"] = (logtau_min, logtau_max)

    if need_hydro:
        logg = parse_logg(path)
        p_hydro_dyncm2 = integrate_hydrostatic(logtau, pg_dyncm2, kappa_ross, logg)
        rel_dev = np.abs(p_hydro_dyncm2 / pg_dyncm2 - 1.0)

        if report_hydro:
            print(f"-- hydrostatic mode (integrating dPg/dtau = g/kappa_ross, log(g) = {logg:.2f}) --")
            print(f"  Consistency with tabulated Pg: mean deviation {np.mean(rel_dev) * 100:.1f}%, "
                  f"max {np.max(rel_dev) * 100:.1f}% (measures how well grey hydrostatic "
                  f"equilibrium alone reproduces this run's actual Pg)")
        logtau_min, extrap_min, dex_min = logtau_at_pressure_hydrostatic(
            logtau, kappa_ross, p_hydro_dyncm2, args.min_pressure, logg)
        logtau_max, extrap_max, dex_max = logtau_at_pressure_hydrostatic(
            logtau, kappa_ross, p_hydro_dyncm2, args.max_pressure, logg)
        if report_hydro:
            report("minimum logtau needed", args.min_pressure, logtau_min, extrap_min, dex_min)
            report("maximum logtau needed", args.max_pressure, logtau_max, extrap_max, dex_max)
            print(f"  -> recommended logtau range: {logtau_min:.2f}  to  {logtau_max:.2f}")
            if extrap_min or extrap_max:
                print("  Note: extrapolated values assume kappa_ross stays at its mean edge value")
                print("  beyond the computed grid -- treat as a starting guess.")
            print()
        results["hydrostatic"] = (logtau_min, logtau_max)

    if args.mode == "both":
        (t_min, t_max) = results["table"]
        (h_min, h_max) = results["hydrostatic"]
        print(f"Comparison: table gives logtau [{t_min:.2f}, {t_max:.2f}], "
              f"hydrostatic gives [{h_min:.2f}, {h_max:.2f}]")
        print()

    if emitting:
        logtau_min, logtau_max = results[args.emit_mode]
        line1, line2 = format_logtau_cards(logtau_min, logtau_max, args.differ)
        if args.emit_cards:
            print(f"-- LOGTAU cards ({args.emit_mode} mode) --")
            print(f"  {line1}")
            print(f"  {line2}")
            print()
        if args.write_input:
            if ":" not in args.write_input:
                sys.exit("error: --write-input expects TEMPLATE:OUTPUT (e.g. old.input:new.input)")
            template_path, out_path = args.write_input.split(":", 1)
            patch_input_file(template_path, out_path, logtau_min, logtau_max, args.differ)


if __name__ == "__main__":
    main()
