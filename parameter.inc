*  INCLUDE FILE FOR MARCS   version 2005 DJ

* Number of depth points in the model.
      INTEGER   NDP
      PARAMETER ( NDP=300)

* Maximal number of rays for the solution of the transfer equation.
* This should be at least NDP/2 +3.
      INTEGER   NRAYS
      PARAMETER ( NRAYS=120) 

* Equation of state parameters   (for tabs and absko)
* mkomp  = max # of absorption coefficient contributors
* mkompr = max # of temp. indep. absorption coefficient contributors
      integer mkomp,mkompr,ifadim,kfadim
      parameter (mkomp=27,mkompr=16,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp)

* Number of wavelength points in the  OS calculation
      INTEGER NWL
      PARAMETER (NWL=80000)

* Max number of points in any single "WAVELENGTH SETS" block read by
* VAAGL from an .input file. Sizes the XLA/XLA3 arrays in COMMON /CA5/
* (marcs.f, ABSKO/DETABS) and the equivalently-sized local XLA/XLA3
* buffer INABS fills per-set to write/read the ABKOF opacity-table
* scratch file. THESE ARE NOT AUTO-SIZED TO THE INPUT FILE -- a set
* with more than maxwlset points silently reads/writes past the end of
* those fixed arrays rather than raising any error (traced 2026-09-14:
* produces a garbage/occasionally-Infinity Rosseland weight for the
* affected indices, corrupting the Rosseland mean for every layer).
* INITAB checks NL(1)/NL(2) against this and STOPs with a clear message
* if exceeded, rather than letting it corrupt silently -- if you hit
* that STOP, raise this number to fit your largest wavelength set and
* recompile (make) before rerunning.
* In practice only the Rosseland set (set 1, NL(1)) is at real risk of
* exceeding this: the general/detailed set (set 2 onward) is chopped by
* SETDIS into sub-sets of at most LMAX=15 points each (marcs.f, call
* site in INITAB uses LMAX=15), so it self-limits regardless of the
* total general-set size in the .input file. The Rosseland set gets no
* such chunking. Default (100) covers every Rosseland set in
* input_files/ as of 2026-09-14 (stellar: 20, earth: 28, warm_planet:
* 33) with headroom for future ones.
* Also sizes XLR/ABSKA/SPRIDA in LISTMO and ABSKTR/SPRTR in ARCHIV
* (marcs.f) -- same disease, found 2026-09-14 right after the ABSKO/
* DETABS/INABS fix above: ARCHIV writes, and LISTMO reads back, NLP=
* NL(1)+1 elements of the Rosseland wavelength/opacity arrays when
* saving/listing a converged model, and those were ALSO hardcoded to
* 20 regardless of the actual Rosseland set size, corrupting memory
* the same way (confirmed: two independent converged runs using a
* 28/31-point Rosseland set crashed here, one with a Fortran I/O
* error, one with a segfault -- same root cause, different memory-
* layout manifestation) right after otherwise-successful convergence.
      INTEGER MAXWLSET
      PARAMETER (MAXWLSET=100)

* number of wavelengths after MARCS opacity sampling
      INTEGER NWREAL
      parameter (NWREAL=7949)

* For call to INTERP from OSTABLOOK
      integer mtemp,mpe
      parameter (mtemp=12,mpe=1)

* NUMBER OF ELEMENTS TO DELIVER IONS in JON (and joinabs.dat)
      integer n16,n17
      parameter (n16=16,n17=17)

* Include file for testing Tsuji's molecular equilibrium
* maxmol is maximum number of molecules to be stored in 'takemolec'
* maxmet is maximum number of atoms
* npoint is the number of copmuting temperatures
      integer maxmol,maxmet
      parameter (maxmol=350, maxmet=100)

* maxosmol is the maximum number of molecules for which opacity is included
* (including one number for the atomic lines). Bumped 100->160 2026-08-14:
* data/mol_names_cool_WIP_126.dat (126 molecules, updated + new linelists)
* hit "Increase dimension for MAXOSMOL" at 100 -- headroom kept above 126
* for further WIP-list growth.
      integer maxosmol
      parameter (maxosmol=160)

* NOPI and NOTI set the maximum possible number of points in pressure
* and temperature in the opacity tables. Those are regulated and checked by
* OPAMOD
* Note: make sure to change in .github/build.
* NOSPEC bumped 75->160 2026-08-14 alongside maxosmol: was sized exactly
* to the old 75-molecule production mol_names list, so it hit "increase
* dimension NOSPEC to molno" as soon as molno=126 (mol_names_cool_WIP_126.dat).
* Temporarily dropped to 80 to isolate the 75-species checkpoint sweep
* (all passed clean, confirmed by rerunning the identical 75-species file
* against a NOSPEC=160 rebuild -- THAT crashed, isolating the regression
* to the NOSPEC=160 rebuild itself, unrelated to species count/content).
* Testing NOSPEC=126 (exact molno of the full WIP list) to binary-search
* the real working threshold between 80 (known good) and 160 (known bad).
      integer NOPI, NOTI, NOSPEC
      parameter (NOPI = 30, NOTI = 80, NOSPEC = 126)

* natms = number of atoms considered in Gibbs energy minimalization routines
* nspec = total number of species (atoms, ions, molecules) considered in Gibbs 
*         energy minimalization routines (= 0.5 * number of lines in 
*         /home/ast/uffegj/falkesg/gfits.data
      integer natms,nspec, nions
      parameter(natms=99,nions=127,nspec=896)       ! # elements and species


* Virga
* nr_cond number of condesnates considered in VIRGA, condensates species shall be put into virga_condensates.dat
      integer nr_cond
      parameter (nr_cond = 1)

* CIA (collision-induced absorption), repurposing LIN_CIA.
* maxcia   = max number of active pairs any single run's
*            data/cia_pairs_<target>.dat can list (not the number of
*            pairs converted by the python pipeline -- see
*            output/cia/, which has all 24 downloaded HITRAN pairs).
*            Earth uses 21 (all 24 HITRAN pairs minus N2-air/O2-air,
*            excluded since HITRAN's -air coefficients are defined
*            against total air density and would double-count the
*            explicit N2/O2 pairs already listed -- widened 2026-08-10
*            from the original 4 N2/O2/CO2-only set); headroom kept for
*            other targets' pair lists.
* cia_maxt/cia_maxwn = per-pair native-grid ceilings. Checked across all
*            24 converted pairs (2026-08-10, widened from the original
*            4-pair check on 2026-08-08): ktemp peaks at 334 (H2-He),
*            nwnos peaks at 21600 (CO2-CO2). Generous headroom kept.
      integer maxcia, cia_maxt, cia_maxwn
      parameter (maxcia=24, cia_maxt=350, cia_maxwn=25000)