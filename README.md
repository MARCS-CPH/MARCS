# MARCS
(A full list of all input parameters can be found in parameters_list.txt.
For explanations of individual physics modules, see PHYSICS_OVERVIEW.md.)
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
  - krome/data/database/leiden_xsecs and krome/data/database/swri_xsecs -- these are symlinks,
    not tracked directories (the LEIDEN/SWRI photodissociation cross-section files are too large/
    numerous to track in git). They point at /groups/astro/tbalduin/krome_leiden_xsecs/ and
    /groups/astro/tbalduin/krome_swri_xsecs/ respectively, so a fresh checkout on a machine
    without access to those paths will need the symlink targets updated to wherever that data
    has been copied to.
    (For now this assumes running on the Copenhagen HPC, where every current user has access to
    those paths -- this repository is not yet set up for use outside that environment. Making
    this data publicly downloadable is a planned future improvement.)
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
