# MARCS
How to run a marcs model.  
In order to run a marcs model from this repository for the first time, a few steps have to be taken.  
### Prepare input files
Firstly all files with "_basic" should be copied/moved to their respecitve names without "_basic".
(This is just in place to avoid git tracking of the actual input files)  
Afterwards one should check that all required files are existing and in the right directory for marcs to work.
Here is a little checklist:  
  - elabund.dat in the **data** directory (make sure it is the right elemental abundace, for example solar or earth like)
  - jonabs.dat in **data**
  - mol_names.dat in **data** (make sure this contains all the molecules and atoms which you want the code to calculate the opacities for)
  - all dispol and DustChem files in **data** (should always be the case)
  - parameter.inc in the main directory
  - marcs.input in the main directory (make sure to align all the input parameters correctly, for example have equal signs be aligned etc.)
        - in addition if you do not want non equilibrium chemistry to be calculated set NONEQ to zero (maybe a good idea for the first test run)
  - runmarcs file in the main directory (make sure all the output you want is comment out/in)

### Run MARCS
Afterwards you can compile marcs either with one of the commands from compile.txt, or by executing the "make" command if you use KROME.  
Then you can run marcs by either executing runmarcs or by adding it to your HPC queue.  
(Note that the runmarcs file in this repo assumes you'll do the later so you might need to adjust it if you run marcs locally)

# KROME
Running MARCS with KROME

In order to run KROME with MARCS a few steps have to be added to the above mentioned checklist.  

### Prepare network
Before you can run a model with KROME you have to prepare the network you want to run.  
You can find the relevant networks in the folder krome/networks.  
There you can find some already prepared in the subfolder "noneq".  
If you know which network you want to run you should go to the file "compile_krome.sh"
and make sure to add your network path to this line "./krome -n networks/ADD_YOUR_NETWORK_PATH_HERE".  
When this is in place make sure the first line "project" also gets an approriate name.  

### Compile and run
Then you can execute "compile_krome.sh" and your krome build should get compiled.
(for the example case in this repository its is advised to choose the "react_Chapman_incl_photo" network)  
When KROME is compiled, you can switch back to the main directory and compile marcs with krome.
This is done by simply excuting the "make" command.  
The relevant compiling options can be found in the makefile.  
(Note that at this point the debugging flags are very problematic as they also show problems with krome itself.
It is highly advised to always use the optimised flags for compiling)  
Once you compiled marcs with krome you can run marcs by executing the runmarcs file as usual.  
Make sure to comment in/out the krome output that you wanted to see in the runmarcs file.

### MARCS Noneq input and KROME compilation flags

This part just serves as a quick explanation of the new noneq input parameters found in "marcs.input" and a short summary of important krome flags.

The noneq input parameters are:
 - NONEQ (Basic on/off switch for non equilibrium chemistry. Expects a numeric value of 0 (off) or 1 (on))
 - PHOTO (Basic on/off switch for including photo rates. Expects a numeric value of 0 (off) or 1 (on).
          currently only able to turn on/off the photolysis module for ozone with more potential features later)
 - DTMIN (Starting timstep for KROME in seconds. Expects a Value in the format of X.XE+/-XX)
 - DTMAX (Maximal timestep size for KROME in seconds. Expects a Value in the format of X.XE+/-XX)
 - tMAX  (Final time the network will be solved for in seconds. Expects a Value in the format of X.XE+/-XX)
 - DTINC (By how much the starting timestep should increase every iteration until DTMAX is reached. Expect a value in the format of X.XX.
         Usually it is advised to have some increase or choose starting timestep and final time wisely to not have your calculation run for too long.)
 - KROMEO (KROME output parameter. Expectes either 1,2 or 3.  
          1 is just a full output at the end. Good for debugging.  
          2 is just the final output. Should be default for normal operations  
          3 is both outputs.)
 - KROMER (Switch that determines whether the krome calculations should be retunred to MARCS itself. 0 is off 1 is on.
           Currently still a bit work in progress.)
 - OUTINT (How often should the output be written out in case you choose to the full output. Expects a value in the format X.XE+/-XX.
           probably only needs to be changed when debugging.)

For a full explanation of krome specific compiling flags please refer to https://bitbucket.org/tgrassi/krome/wiki/optionsALL.
A few set of quick adjustments that can be done if the solver is an unstable are:

-ATOL (Absolute tolerance of the solver. Default 1d-20. Can be used to make the solver more or less accurate to ensure stability if needed or performance if needed.
      also has an option to define custom ATOL for each species.)
-RTOL (Relative tolerance of the solver. Default 1d-4. Can also be used to make the solver more or less accurate to ensure stability if needed or performance if needed.) 



# SW 
This file explains how to run MARCS with Static Weather,
in order to produce MARCS models with clouds.

The driftmarcs.f90 file has the conjugation of the calls to both 
codes as well as important parameters between the codes comunication.

At the beginning of driftmarcs.f90 you will find the following variables 
declared and defined:
. f_opac: defines the fraction of the cloud opacity to be considered in marcs.
in the first instance one should set this to the cloud opacity fraction one wants to start with.
0 means no clouds, 1.0 means the full cloud calculated by driftmarcs.

. delta_t: the relative error between two models after a MARCS+SW run. Initially set to a high value (100.0 fex) 
so that the code always has to perform at least 2 runs of MARCS+SW before convergence.

. delta_opac: the steps in cloud opacity to be taken between each run.

. idriftok: is 1 if SW has ran successfully.

Files to care about:
. opac_frac.in: store the current cloud opacity fraction to be taken.
It is read my MARCS.

.marcs2drift.dat: read by MARCS with SW's cloud information.

.out3_dust.dat: contains main output from SW.

.sw.out: contains log from SW run, where one can check for convergence.
If the statement "regular end of integration" is found, then SW
has converged and successfully ran.

In driftmarcs.f90 we encounter two WHILE loops.

The first while loop:
    Here we run MARCS+SW while slowly adding the cloud to MARCS.
    For example, if f_opac is set to 0.20 at the beginning and delta_opac 
    to 0.20, then we will call MARCS+SW 5 times, each time with 20% more opacity.
    In the final model we will have a model with the full cloud, but this
    does not mean it is fully converged yet. This is where the second loop 
    comes in...

The second while loop:
    Loops through MARCS+SW at 100% cloud opacity until the MARCS model previously
    computed has a relative error in temperature of less than 10% - this value can be changed.
    If convergence is not obtained within 20 loops, the code stops. - this value can also be modified.
