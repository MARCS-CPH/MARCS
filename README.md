# MARCS

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

# KROME
Brief description for using KROME with MARCS

Enter the krome subdirectory
Modify the file run_MARCS_KROME.sh to reflect the network wanted
Return to "main" directory by cd ../
Compile the combined MARCS and KROME by running:
make
Now you can run MARCS via runmarcs as usual

Changes made to the default Makefile for the linking to MARCS to work as intended

Assumption:
marcs.f is in the “main” directory and all files related to krome are in a sub-directory called “krome”.
Inside the krome directory, the files created by running ./krome (or better run_MARCS_KROME.sh) are in MARCS_build 
Exception: the file “reactions_verbatim.dat” must be in the “main” directory
Desired result: marcs executable is created in the "main" directory and everything else is kept in "krome/MARCS_build"

Added at the top of the Makefile to tell the compiler to look in this directory
#Paths
VPATH = ./krome/MARCS_build

Changed the name of the intended executable
#executable name
exec = marcs

Added path for krome_subs.f90 file
GREP = $(shell grep -i 'dgesv' krome/MARCS_build/krome_subs.f90)

Specified location of module files to be “krome/MARCS_build” directory using “-module krome/MARCS_build”
switchOPT = -O3 -ipo -ip -unroll -xHost -g -module krome/MARCS_build -fp-model precise

All the objects have the path “krome/MARCS_build” added in front of their filenames

Modified the default target
#default target
all:  $(objs) marcs.o
      $(fc) $(objs) marcs.o -o $(exec) $(switch) $(lib)

Included specifically the -save compiler option for the marcs compilation to save variables in static memory to aviod segmentation violation
#Special rule for marcs
marcs.o:marcs.f
      $(fc) -save $(switch) $(nowarn) -c $^ -o $@






