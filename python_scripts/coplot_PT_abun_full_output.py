import numpy as np
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

######## SETTINGS ########

plt.style.use("/groups/astro/tbalduin/python_scripts/MyStyle.mplstyle")

output_path = "/groups/astro/tbalduin/plots/"

altitude_axis = False  # If True, add an altitude axis on the right side of the plot

input_data = '/groups/astro/tbalduin/MARCS_CPH/MARCS/models/K2_18/K2_18b/Marcus/new/non_eq/full_DMS_network/rerun/full_output.dat'
input_data = '/groups/astro/tbalduin/MARCS_CPH/MARCS/models/K2_18b_highZ_CNOS/chain/z_200_v6_DMS/full_output.dat'
#input_data = '/groups/astro/tbalduin/MARCS_CPH/MARCS/models/mixing_test/dtinc_5/full_output.dat'
#input_data = '/groups/astro/tbalduin/MARCS_CPH/MARCS/models/mixing_test/dtinc_5/no_safeguard/full_output.dat'
input_data = '/lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/models/archean/earth_modern/test29_modern_cont25/full_output.dat'
input_data = '/lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/models/archean/earth_modern/test30_grid1bar/full_output.dat'
input_data = '/lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/models/archean/earth_modern/test38_lowteff/full_output.dat'
input_data = '/lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/models/archean/earth_modern/test44_convfix_cont50/full_output.dat'
input_data = '/lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/models/archean/earth_archean/test77_faintsun/full_output.dat'

# Mode: 'manual'    -> plot only species listed in selected_species
#        'threshold' -> plot all species whose max log10-abundance exceeds abundance_threshold
mode = 'manual'

# Manual mode: list of species names exactly as they appear in the file header
selected_species = ['H2', 'H2O', 'CO2', 'CO', 'CH4', 'N2', 'NH3', 'H2S', 'CH3SCH3']
selected_species = ['H2', 'H2O', 'CO2', 'CO', 'CH4', 'N2', 'NH3', 'H2S', 'CH3SCH3']
selected_species=['H2', 'H2O', 'CO2', 'CO', 'CH4', 'N2', 'NH3', 'H2S', 'CH3SCH3','H','S','N2','CH','S2','SH', 'C5H10O2S', 'SO2', 'H2SO4']
selected_species=['H2', 'H2O', 'CO2', 'CO', 'CH4', 'N2', 'NH3', 'H2S','H','S','N2','CH','S2','SH','SO2', 'H2SO4']
selected_species = ['N2', 'O2', 'O3', 'CO2', 'CH4']
# Threshold mode: include species with max log10(abundance) > abundance_threshold
# Values of -99 (GGchem zero-floor) are excluded when computing the maximum.
abundance_threshold = -3.0

# Lower x-limit for the mixing-ratio axis (linear, e.g. 1e-12)
xmin_vmr = 1e-12

input_label = ['TP profile']

save_fig_name = "TP_abund_earth_archean_test77_faintsun"

dyn_to_bar = 1e-6


def find_nearest(array, search_for):
    diff = abs(np.array(array) - search_for)
    return int(np.argmin(diff))


######## READ DATA ########

with open(input_data) as f:
    all_lines = f.readlines()

# --- Locate section 1: "M O D E L  A T M O S P H E R E" ---
sec1_start = None
for i, line in enumerate(all_lines):
    if "M O D E L  A T M O S P H E R E" in line:
        sec1_start = i + 2  # skip model header + column-name line
        break

if sec1_start is None:
    raise ValueError("Could not find 'M O D E L  A T M O S P H E R E' in input file")

# Collect section 1 rows (lines starting with an integer layer index)
# Stop when a non-integer-starting line is encountered -> that is the species header.
sec1_rows = []
sec2_header_idx = None
i = sec1_start
while i < len(all_lines):
    tokens = all_lines[i].strip().split()
    if not tokens:
        i += 1
        continue
    try:
        int(tokens[0])
        sec1_rows.append(tokens)
    except ValueError:
        sec2_header_idx = i
        break
    i += 1

if sec2_header_idx is None:
    raise ValueError("Could not find section-2 species header line")

nr_layers = len(sec1_rows)
print(f"Detected {nr_layers} atmospheric layers")

# Parse section 1: same column layout as marcs_output.dat
# [0]=K  [1]=TAUROSS  [2]=TAU(5000)  [3]=GEOM.DEPTH  [4]=T  [5]=PE  [6]=PG  ...
pressure_data   = np.zeros(nr_layers)   # PG [Ba = dyn/cm^2]
temperature_data = np.zeros(nr_layers)  # T [K]
geom_depth_data = np.zeros(nr_layers)   # Geometric depth from top [cm]

for ii, row in enumerate(sec1_rows):
    geom_depth_data[ii]  = float(row[3])
    temperature_data[ii] = float(row[4])
    pressure_data[ii]    = float(row[6])

# --- Section 2: species abundances ---
# Header line: species names separated by whitespace
species_names = all_lines[sec2_header_idx].strip().split()
n_species = len(species_names)
print(f"Found {n_species} species in section 2")

# Data rows: [layer_index]  [log10(p_i/Ba)]  [log10(p_i/Ba)]  ...
# -99 is the GGchem floor value (effectively log10(0)); stored as NaN.
abun_matrix = np.full((nr_layers, n_species), np.nan)
for ii in range(nr_layers):
    row = all_lines[sec2_header_idx + 1 + ii].strip().split()
    for jj in range(n_species):
        val = float(row[jj + 1])   # +1 to skip layer index
        if val > -98.0:
            abun_matrix[ii, jj] = val

# Volume mixing ratio: VMR_i = p_i / p_tot = 10^(log10 p_i) / p_tot
# pressure_data is in Ba (dyn/cm^2), same unit as the GGchem log10 values.
vmr_matrix = 10**abun_matrix / pressure_data[:, np.newaxis]

######## SELECT SPECIES ########

if mode == 'manual':
    selected_indices = []
    selected_labels  = []
    for name in selected_species:
        if name in species_names:
            selected_indices.append(species_names.index(name))
            selected_labels.append(name)
        else:
            print(f"Warning: species '{name}' not found in file header")

elif mode == 'threshold':
    selected_indices = []
    selected_labels  = []
    for jj, name in enumerate(species_names):
        col_vmr = vmr_matrix[:, jj]
        col_vmr = col_vmr[~np.isnan(col_vmr)]
        if len(col_vmr) > 0 and np.max(col_vmr) > 10**abundance_threshold:
            selected_indices.append(jj)
            selected_labels.append(name)
    print(f"Threshold mode (VMR > 1e{abundance_threshold:.0f}): {len(selected_labels)} species selected")
    print(f"With these species selected:",selected_labels)

else:
    raise ValueError(f"Unknown mode '{mode}'. Use 'manual' or 'threshold'.")

######## PLOT ########

fig, ax1 = plt.subplots()

ax1.set_ylabel(r"Pressure [bar]")
ax1.set_xlabel(r"Volume Mixing Ratio [$p_i / p_{tot}$]")

pmax = np.min(pressure_data) * dyn_to_bar   # top of atmosphere (lowest pressure)
pmin = np.max(pressure_data) * dyn_to_bar   # bottom of atmosphere (highest pressure)
#pmin = 1  # grid already stops at ~1 bar for test30_grid1bar, no floor needed
ax1.set_ylim(1.1 * pmin, 0.9 * pmax)
ax1.set_yscale('log')

# Determine x range from selected VMR data
all_vmr_sel = vmr_matrix[:, selected_indices]
valid_vmr = all_vmr_sel[~np.isnan(all_vmr_sel)]
xmax_vmr = np.max(valid_vmr) * 2 if len(valid_vmr) > 0 else 1.0
ax1.set_xlim(xmin_vmr, xmax_vmr)
ax1.set_xscale('log')

for idx in selected_indices:
    ax1.plot(vmr_matrix[:, idx], pressure_data * dyn_to_bar)

lines_style = [Line2D([0], [0], linestyle='dashed', color='black')]
legend1 = plt.legend(lines_style, input_label, fontsize="12", loc='lower left', bbox_to_anchor=(0., 0.2))
legend2 = plt.legend(selected_labels, fontsize="8", loc='lower left', bbox_to_anchor=(0., 0.025),ncols=3)
ax1.add_artist(legend1)
ax1.add_artist(legend2)

ax2 = ax1.twiny()
ax2.plot(temperature_data, pressure_data * dyn_to_bar, label="Temp", color="black", linestyle='dashed')
ax2.xaxis.label.set_color('grey')
ax2.tick_params(axis='x', colors='grey', labelsize=20)
ax2.set_xlabel(r"Temperature [K]", labelpad=10)
ax2.set_yscale('log')

if altitude_axis:
    ax1.tick_params(right=False, labelright=False, which='both')
    ax3 = ax1.secondary_yaxis('right')
    ind_1bar = find_nearest(pressure_data, 1 / dyn_to_bar)
    altitude_data = (geom_depth_data[ind_1bar] - np.array(geom_depth_data)) / 1e5  # cm -> km
    altitude_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80]
    if max(altitude_ticks) > max(altitude_data):
        print('Chosen altitude ticks exceed data range')
    altitude_ticks = [n for n in altitude_ticks if n <= max(altitude_data)]
    pressure_ticks = [pressure_data[find_nearest(altitude_data, t)] * dyn_to_bar for t in altitude_ticks]
    ax3.set_ylabel(r"Altitude [km]")
    ax3.tick_params(right=False, labelright=False, which='minor')
    ax3.set_yticks(pressure_ticks)
    ax3.set_yticklabels(altitude_ticks, fontsize=20)

plt.savefig(output_path + save_fig_name + ".png")
print(f"Saved to {output_path + save_fig_name}.png")
