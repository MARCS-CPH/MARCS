# MARCS
(A full list of all input parameters can be found in parameters_list.txt)
How to run a marcs model.  
In order to run a marcs model from this repository for the first time, a few steps have to be taken.  
### Prepare input files
One should check that all required files are existing and in the right directory for marcs to work.
Here is a little checklist:  
  - elabund.dat in the **data** directory (make sure it is the right elemental abundace, for example solar or earth like)
  - jonabs.dat in **data**
  - mol_names.dat in **data** (make sure this contains all the molecules and atoms which you want the code to calculate the opacities for)
  - all dispol and DustChem files in **data** (should always be the case)
  - parameter.inc in the main directory
  - Makefile in the main directory
  - the ggchem16 executable in **GGchem/src16** (build it by running "make" inside GGchem/src16/;
    this binary is gitignored, so a fresh checkout won't have it, and marcs.f calls it through the
    GGchem/ggchem symlink, which points at GGchem/src16/ggchem16 and will otherwise be dangling)
  - marcs.input in the **input_files** directory (make sure to align all the input parameters correctly, for example have equal signs be aligned etc.)
        - in addition if you do not want non equilibrium chemistry to be calculated set NONEQ to zero (maybe a good idea for the first test run)
  - runmarcs file in the **runmarcs_files** directory (make sure all the output you want is comment out/in)
  - make sure that you have a compiler installed and available to your shell. MARCS by default uses the ifort compilers, but gfortran is possible if needed.   
    (on the Copenhagen HPC this is done by adding the line "module load intel" to your .bashrc in your home folder. Make sure to restart your terminal/shell after adding the line)

### Run MARCS
Afterwards you can compile marcs either with one of the commands listed in the file compile.txt, or by executing the "make" command if you use KROME.  
Then you can run marcs by either executing runmarcs or by adding it to your HPC queue.  
(Note that the runmarcs file in this repo assumes you'll do the later so you might need to adjust it if you run marcs locally)

# KROME
(For a full explanation of krome specific compiling flags please refer to https://bitbucket.org/tgrassi/krome/wiki/optionsALL.)
Running MARCS with KROME

In order to run KROME with MARCS a few steps have to be added to the above mentioned checklist.  

### Prepare network
Before you can run a model with KROME you have to prepare the network you want to run.  
You can find the relevant networks in the folder krome/networks.  
There you can find some already prepared in the subfolder "noneq". 
(for the example case in this repository it is advised to choose the "DMS_v2" network,
krome/networks/non_eq/full_network_DMS_v2.dat -- the current standard/most up to date network)  
If you know which network you want to run you should go to the file "krome/compile_KROME.sh"
and set the `network=` variable near the top to your network's path (e.g.
`network=networks/non_eq/full_network_DMS_v2.dat`), then give an appropriate
name to the `project=` variable right below it.  

### Compile and run
Then you can execute "krome/compile_KROME.sh" and your krome build should get compiled.     
When KROME is compiled, you can switch back to the main directory and compile marcs with krome.
This is done by simply excuting the "make" command.  
The relevant compiling options can be found in the makefile.  
(Note that at this point the debugging flags are very problematic as they also show problems with krome itself.
It is highly advised to always use the optimised flags for compiling)  
Once you compiled marcs with krome you can run marcs by executing the runmarcs file as usual.  
Make sure to comment in/out the krome output that you wanted to see in the runmarcs file.

# Additional physics

### CIA (collision-induced absorption)
MARCS can include collision-induced absorption via the `LINCIA` input
parameter (0 = off, the default in every existing .input file; 1 = on).

When enabled, MARCS reads its list of active pairs from a fixed path,
`./data/cia_list.dat` -- this file is **not** created automatically, so you
need to put one there yourself before setting LINCIA=1. Format: first line
is the number of pairs, then one line per pair as
`SPECIES1  SPECIES2  /path/to/converted/pair/data/`.
`data/cia_pairs_all.dat` lists every CIA pair currently converted and
available (21 pairs) in exactly this format and can be used as a source to
copy the pairs you want into `data/cia_list.dat` -- it is a reference list
only and is not read directly by MARCS.

Pairs are matched generically by species name against GGchem's species
tables, so any pair present in the converted data can be used, not just the
originally-targeted N2/O2/CO2 set for Earth-type atmospheres.

### N2/O2 Rayleigh scattering
N2 and O2 Rayleigh scattering is included automatically whenever those
species are present in the atmosphere. There is no input switch for this --
it runs independently of LINCIA/CIA.

### conv_crit (chemistry integration length)
Part of the noneq input block (see parameters_list.txt). A value of exactly
`0` switches per-layer KROME chemistry integration to always run for the
full `tMAX`, instead of the default `min(vert_mix_time, tMAX)`. Any nonzero
value keeps the default behavior. This gives direct input-file control over
integration length independent of the vertical mixing timescale.
