#F90          = ifx
#F90LINKER    = ifx
F90          = ifort
F90LINKER    = ifort

# for IFORT-Compiler (matches the toolchain marcs.f itself is built with)
#===================
#DEFS     = -D IFORT -diag-disable 8290,8291
#FFLAGS   = -r8 -i4 -g -traceback -fp-stack-check -check all,noarg_temp_created -gen-interfaces -warn interfaces -fpe0 -mcmodel medium -shared-intel -qopenmp
#FFLAGS   = -r8 -i4 -O3 -msse4.2 -xHOST -prec-div -fp-model source -mcmodel medium -shared-intel -qopenmp
DEFS      =
FFLAGS    = -r8 -g -traceback -O3
INCLUDES  =
LFLAGS    = $(FFLAGS)

# for gfortran-Compiler
#======================
#F90          = gfortran
#F90LINKER    = gfortran
#DEFS      =
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -g -fbounds-check -fbacktrace -std=legacy
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -g -O5 -std=legacy
#INCLUDES  =
#LFLAGS    = $(FFLAGS)

EXECS     = ggchem16

OBJECTS  = datamod.o \
           database.o \
	   init_opac.o \
           read_parameter.o \
           demo_sweep.o \
           demo_structure.o \
	   demo_phasediagram.o \
	   adapt_condensates.o \
	   equil_cond.o \
	   is_nan.o \
	   init.o \
           nasa_polynomial.o \
           init_chemistry.o \
	   init_dustchem.o \
	   ggchem.o \
           gauss16.o \
           gauss8.o \
	   gauss_nm.o \
	   massfrac.o \
	   linpack_q.o \
	   slatec_routines.o \
           smchem8.o \
           smchem16.o \
           supersat.o \
	   stindex.o \
           upper.o \
           nucleation.o \
	   auto_structure.o \
	   MIEX.o \
	   q_mie.o \
	   calc_opac.o \
	   effMedium.o \
	   fastmie.o \
	   disk_problem.o \
	   kappaRoss.o \
           main_marcs.o

.SUFFIXES : .o .f .F .f90

default: ggchem

ggchem:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

# MARCS integration: same OBJECTS, but main_marcs.o (which additionally
# writes pp.dat/ndensity.dat/GGchem_ppel for marcs.f to read back) replaces
# main.o. Built separately since it needs its own object files (main.o vs
# main_marcs.o can't share a link line), same FFLAGS/compiler as above.
OBJECTS_MARCS = $(OBJECTS:main.o=main_marcs.o)
EXECS_MARCS   = ggchem16_marcs

ggchem_marcs: $(OBJECTS_MARCS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS_MARCS) -o $(EXECS_MARCS)

clean:
	rm -f *.o 

realclean:
	rm -f *.o *~ *.mod *_genmod.f90 $(EXECS)

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

