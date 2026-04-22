!6-----------------------------------------------------------------------
C BEM�RK:         RS=0.00465047
C      RAU=1.00
C      STEFF=5770

!         version 2024      by KHM and TB: implementation of NonEq. Chemistry via KROME
!         version 2022      by BCE: update coupling to GGchem and to SW (DRIFT)
!                                   and to add irradiation
!         version 2018      by ERC: coupling to GGCHEM to model down to 100 K
!         version 2015      by DJ: coupling to DRIFT to model clouds
!         version 00-03-22  by UGJ includes JF's equilibrium and C2H opac.
!         version 97-03-14  by Uffe Graae Jorgensen: double precision
!         version 94-12-15  by Christiane Helling: Tsuji's molecular eq.
!  VMS    version 93-01-01  by Uffe Graae Jorgensen (Opacity Sampling, Sphericity)
!         version 90-03-07  by Aake Nordlund (stripped all line code, ODF)
!  Cyber  version 1975      by Gustafsson et al. (ODF)
!
! Iteration control parameters:
! NEWMOD can take the following values:
!     000=> start from a stored model
!     001=> start a new model
!     002=> no iterations; transfer in old model
!     003=> continue with new teff,logg,a/h
!     004=> continue with new tau-scale
!     005=> combine old model with temp from other file 
!     006=> start from a binary stored model
!     009=> one iteration to compute the limb functions, no correction on
!           the model, be careful to use the same parameters in the input file.
! If NEWMOD=1, one reads the last 3 entries of the input file, namely: 
!     taucnv=> approx. tau where convection starts
!     dtblnk=> approx. backwarming and surface cooling
!     taubln=> turnover between backw. and surf. cool.
!
! JONTYP can be:
!     001=> constant partition functions
!     002=> full pf, Bashek et al.
!     003=> full pf, Fischel et al., all possible ionisation states
!  
! JUMP can be:
!     1-3 => subroutine MOL doesnt work => no presmo
!     1 => Tsuji-routines, partryck for molecules, xmettryck for atoms
!     2 => JFF routines after same principle as Tsuji jump=1 routines
!     3 => JFF Gibbs minimising routines => GEM-package
!     4 => ERC added GGCHEM code: Equilibrium chemistry down to 100 K
!-----------------------------------------------------------------------
      program scmarcs
      
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      integer molh, jump
      character molname*4,osfil*60,sampling*3,pp_sph*3
      logical pf,pfe,pfd,fixros,itstop,quit,onemor
      integer krome_on,krome_photo_on
      integer krome_output,krome_debug,krome_return
      real*8 krome_tmax,krome_photo_scale
      common/carciv/ larciv    !=1 if called from arciv, otherwise = 0
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *  kos_step,nwtot,nosmol,newosatom,newosatomlist
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     *  vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     * ntau,iter
      common /carc3/ f1p,f3p,f4p,f5p,hnic,presmo(33)
      common /cit/it,itmax
      common /ci4/dum,idum(96),molh,jump 
      common /cpf/pf,pfe,pfd,fixros,itstop
      common /carc1/istral,idrab1,idrab2,idrab3,idrab4,idrab5,
     *  idrab6,iarch
      common /newmo/newmod
      common /cmolrat/ fold(ndp,8),molold,kl
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      common /cisph/isph
      common /clist/nlte
      common/ci5/abmarcs(18,ndp),anjon(18,5),h(5),part(18,5),
     *  dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      namelist /outlist/ masslinf
      common /coutlist/pplist(843)
      common /clevprint/ prj2(ndp),masslinf
      common/cabinit/abinit(natms),kelem(natms),nelem
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      character atnames*2, molnames*8, molnames2*4
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
      common /ggchembool/ iggcall
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      common /noneq_time/ dt_start,dt_max,dt_inc,krome_tmax
      common /noneq_output/ krome_output,krome_debug,krome_return
      common /photochem/ FLUX_RAD(ndp,nwreal) !second dimension should be nwtot, in most cases 7949
C    

! Initiation
      larciv = 0
      call gettime(0)
      
      open(unit=5,file='mxms7.input',status='old',readonly)
      open(unit=7,file='mxmodel.dat',status='unknown')
      open(unit=9,file='data/jonabs.dat',status='old',readonly)
      open(unit=16,file='arcivaaa.dat',status='old',readonly)

      read(5,'(6(7x,i3,5x))') itmax,nprint,newmod,noarch,jontyp,idust
      if (idust == 1) idustopac = 1
      read(5,'(7x,i3,12x,a3,2(12x,i3),12x,a3)') jump,sampling,molold,
     *  nlte,pp_sph

      if(pp_sph.eq.'sph' .or. pp_sph.eq.'SPH') isph = 1
      
      if(jump.ne.4) stop 
     * 'Error: This MARCS version only works with GGchem eq. chemistry'
      iggcall = 0
      call oldsta
      call mainb

      io=-1
      call initjn(io)
      call modjon(jontyp,io)
      call initab(io)

      if(idust .eq. 1) then
        if(newmod .ne. 0) stop 'Error: Dust only works for NEWMOD=0'
        if(jontyp .ne. 3) stop 'Error: Dust only works for JONTYP=3'
        if(jump .ne. 2) then
            if(jump.ne.4) then 
              stop 'Error: Dust only works for JUMP=2 or 4'
            end if
        end if
        if(isph .eq. 1) stop 'Error: Dust only works for ISPH=0'
        
        print *, 'Reading DRIFT file...'
        open(unit=976, file='f_cloud.in')
        read(976,*) f_opac
        close(976)
        call drift2marcs
        print *, 'Done.'
        print *
        call dust_opac_eps_interp
        if (icloud_conv == 1) stop
      end if

      molh = 0                           ! molh=1 => only h,h2,h2+ in molec. eq.      
      metals = 0
     
      read(5,outlist)
      pfd=itmax.lt.0
      if(pfd) itmax=-itmax
      pfe=itmax.le.nprint      
      
      if(newmod.eq.1) call startm
      if(newmod.eq.2) print*,' newmo=2; no iteration; trans old mod'
      if(newmod.eq.3) call scale(22)
      if(newmod.eq.4) then
        call resume(22,1)
        if (isph.eq.1) then 
          call tryck_sph
        else 
          call tryck
        end if
      end if
      if(newmod.eq.5) call coscal(22)
      if(newmod.eq.6) call oldarc(22)
      if(newmod.eq.7) read(22)
      if(newmod.eq.8) then
        read(22)
        call scale(22)
      end if
      
      if(newmod.ne.1) then
        if(isph.eq.1) then 
          call presnt_sph
        else 
          call presnt
        end if
      end if

      print *, 'Initiation finished.'
      
! Iterate model
      lun=22
      itstop=.false.
      quit=.false.
      onemor=.false.

      do it=1,itmax 
        print *
        write(*,'(a23,i3,a8)')' Iteration #           ',it,' started'
        call gettime(1)
        
        pf=it.gt.itmax-nprint
        if(quit) pf=.true.

        if (isph.eq.1) then 
          call solve_sph(1)
        else 
          call solve(1)
        end if
        
        if(newmod.eq.2) go to 102
        if(newmod.eq.9) goto 101
        
        if(isph.eq.1) then 
          call matrix_sph
        else 
          call matrix
        end if
        
      
        call newsta
        
        if(itstop.ne..True.) then !catches onemor iterations which are not converged
         onemor=.False.
        endif
        if(itstop .and. onemor) quit=.true.
        if(itstop) onemor=.true.
        if(quit) exit
      end do
      
      
      if(noarch.ge.2) go to 101
      call modjon(3,1)
      if(noarch.eq.1) goto 101
      call newsta
102   continue
      call archiv(22)
              call listmo(1,22,isph)


! End program
101   continue
      write(7,*) 'Model written out successfully.'
      call gettime(1)
      print *
      print *, 'Successful termination of execution!'
      print *
      if(it.gt.itmax) then
        print *, 'WARNING: Maximum number of iterations reached'
      end if

      stop
      end
      
C
      SUBROUTINE ABSKO(NEWT,NT,TSKAL,PESKAL,ISETA,J,ABSK,SPRID,nlayer)
      implicit real*8 (a-h,o-z)
C
C        THE ROUTINE ADMINISTERS THE COMPUTATION OF ABSORPTION
C        COEFFICIENTS. IT CALLS THE ROUTINES, GIVING THE PROPER THERMO-
C        DYNAMIC INFORMATION ( J O N ) , THE DETAILS OF THE ABSORPTION
C        MECHANISMS ( D E T A B S ) AND THE FACTORS FOR THE INTERPOLATION
C        IN T  ( T A B S ) . IT CHOOSES (IF NECESSARY READS) THE RIGHT SET
C        OF ABSORPTION-COEFFICIENT DATA (ISETA), STATEMENT NO. 5 AND MAKES
C        THE INTERPOLATION IN T, STATEMENTS NOS. 10-18, AND THE SUMMATION
C        OF A ROSSELAND MEAN, IF INDICATED BY J = 0, STATEMENTS NOS. 25-28
C
C        NEWT SHOULD BE GT 1 THE FIRST TIME THIS ROUTINE IS USED,
C                       EQ 1 WHEN A NEW SET OF T-PE IS USED,
C                       EQ 0 OTHERWISE.
C
C        NT IS THE NUMBER OF T-PE POINTS. THE TEMPERATURES T SHOULD BE EX-
C        PRESSED IN KELVIN, THE ELECTRON PRESSURES PE IN DYNES PER CM2
C        ISETA IS THE WAVELENGTH-SET NUMBER, J THE WAVELENGTH NUMBER IN THAT
C        SET. J EQUAL TO ZERO INDICATES THAT A ROSSELAND MEAN IS WANTED.
C        THIS MEAN IS COMPUTED USING THE WAVELENGTH POINTS OF THE ACTUAL
C        SET (ISETA) AND THE QUADRATURE WEIGHTS GIVEN IN ROSW.
C        IN ABSK AND SPRID THE ABSORPTION AND SCATTERING COEFFICIENTS PER GRAM
C        OF STELLAR MATTER ARE STORED.
C
C        DIMENSIONS NECESSARY
C        AB(NKOMP),ABSK(1),FAKT(NKOMP),FAKTP(IFADIM),NTPO(NTO),PE(NT),PESKAL(1),
C        ROSW(MAX(NL)),SPRID(1),SUMW(NT),T(NT),TSKAL(1),XLA(MAX(NL)),
C        XLA3(MAX(NL))
C        THE DIMENSIONS ARE LOWER LIMITS.
C        DIMENSIONS OF ARRAYS IN COMMONS /CA2/,/CA3/ AND /CFIL/ ARE COMMENTED
C        ON IN SUBROUTINE INABS, THOSE OF ARRAYS IN COMMON /CA4/ IN SUBROUTINE
C        TABS.
C        NKOMP IS THE NUMBER OF 'COMPONENTS'
C        NL(I) IS THE NUMBER OF WAVELENGTHS IN WAVELENGTH SET I
C        NT IS THE NUMBER OF T-PE POINTS SUPPLIED IN TSKAL AND PESKAL
C        NTO IS THE NUMBER OF POINTS IN THOSE SCALES FOR WHICH A DETAILED
C              PRINT-OUT IS WANTED.
C
      include 'parameter.inc'
C
C      PARAMETER (KFADIM=4000,IFADIM=1000)
      DIMENSION TSKAL(NDP),PESKAL(NDP),ABSK(NDP),SPRID(NDP)
      DIMENSION FAKTP(ifadim)
      DIMENSION SUMW(NDP)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA3/ILOGTA(30),NULL
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
      COMMON/CA5/AB(30),FAKT(30),PE(NDP),T(NDP),XLA(20),XLA3(20),RO,
     &           SUMABS,SUMSCA,VIKTR,ISET,NLB
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/COUTR/NTO,NTPO(10)
      COMMON/CROS/ROSW(20)
      COMMON /CARC3/ F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /CARC4/ PROV(30),NPROVA,NPROVS,NPROV
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
     &             ,FT,FTE
      COMMON /DENSTY/ ROTEST(NDP),PRH2O(NDP)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CMETPE/ PPEL(NDP), METPE
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      common /ggchemdetabs / f1_dt(ndp), f5_dt(ndp), 
     >                       rC(ndp), rMg(ndp), rAl(ndp), 
     >                       rSi(ndp), rHe(ndp), ro_dt(ndp)
      INTEGER MOLH, JUMP
C
C
      ISET=ISETA
      IF(NEWT.GT.1)ISETP=-1
      IF(NEWT.EQ.0)GO TO 5
C
C        FACTORS ONLY DEPENDENT ON T-PE
C
      CALL TABS(NT,TSKAL)

      IFAK=1
      KFAK=1
      JP=0
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE FIRST NTP-LOOP')
      !print*, "NT is ", nt
      DO4 NTP=1,NT
      T(NTP)=TSKAL(NTP)
      PE(NTP)=PESKAL(NTP)
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 3
      IF(NTP.EQ.NTPO(KP))GO TO 1
      GO TO 3
    1 IOUTR=1
      KP=KP+1
    3 CONTINUE
C
      if(j.le.0) then
        molhs=molh
        molh =0
      endif


      if(pe(ntp).le.1d-300) then
        print *, "Electron pressure too small in ABSKO"
        write(*,*) ntp,pe(ntp)
        write(*,*) ntp,PESKAL(ntp)
      end if
      
      CALL JON(T(NTP),PE(NTP),1,PG,RO,DUM,IOUTR)
      
      
      if(j.le.0) then
        !molh=molhs
        molh = 0

        IF (JUMP.GE.1) THEN
         prh2o(ntp)=partryck(ntp,4)
         ELSE
         prh2o(ntp)=presmo(4)
        ENDIF
        
      endif
      if (nlayer == -1 ) then
        CALL DETABS(J,0,NTP,IOUTR,NTP)
      else 
        CALL DETABS(J,0,NTP,IOUTR,nlayer)
      end if     
        
C
C        WE STORE THE FAKT ARRAY, MADE IN JON-DETABS IN LONGER ARRAYS NAMELY
C                  IN AFAK FOR TEMPERATURE-INDEPENDENT COMPONENTS
C                  IN FAKTP FOR TEMPERATURE-DEPENDENT ONES

      DO2 KOMP=1,KOMPR
      AFAK(KFAK)=FAKT(KOMP)
    2 KFAK=KFAK+1

      DO4 KOMP=KOMPS,NKOMP
      FAKTP(IFAK)=FAKT(KOMP)
      KFAK=KFAK+NOFAK(IFAK)
    4 IFAK=IFAK+1
C        END OF 'THE FIRST NTP-LOOP'
C
C        READING  OF A NEW WAVELENGTH SET IF INDICATED BY ISET
    5 IF(ISET.EQ.ISETP)GO TO 6
      IREADP=IRESET(ISET)
   51 READ(IREADP,END=52)ISETP,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
      GO TO 5
   52 REWIND IREADP
      GO TO 51
C        ROSSELAND MEAN OR NOT
    6 IF(J.GT.0)GO TO 9
    7 J1=1
      J2=NLB
      DO8 NTP=1,NT
      SUMW(NTP)=0.
    8 ABSK(NTP)=0.
      GO TO 10
    9 J1=J
      J2=J
C
C        INTERPOLATION IN T
C        LOOP OVER ALL THE WAVELENGTHS IN CASE OF ROSSELAND MEAN. THIS
C        LOOP ENDS IN STATEMENT NO. 26
   10 CONTINUE
      DO26 JP=J1,J2
      KFAK=1
      IFAK=1
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE SECOND NTP-LOOP')
      DO26 NTP=1,NT
C
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 93
      IF(NTP.EQ.NTPO(KP))GO TO 92
      GO TO 93
   92 IOUTR=1
      KP=KP+1
      IF(KP.EQ.2)IOUTR=2
   93 CONTINUE
      IU=JP
C
C        COMPONENTS WITH ABSORPTION COEFFICIENTS INDEPENDENT OF THE
C        TEMPERATURE
C
      DO KOMP=1,KOMPR
      
      IF(KOMPLA(IU).LE.0)GO TO 12
C        THE VECTOR KOMPLA IS DETERMINED IN SUBROUTINE INABS.
C        KOMPLA GREATER THAN ZERO GIVES THE INDEX IN ABKOF, WHERE THE TABLE FOR
C        THIS COMPONENT AND WAVELENGTH BEGINS.
C        KOMPLA LESS THAN OR EQUAL TO ZERO INDICATES THAT THE ACTUAL ABSORPTION
C        COEFFICIENT FOR THIS COMPONENT AND WAVELENGTH IS ZERO, AS FOUND IN SUB-
C        ROUTINE INABS.
C
   11 INDEX=KOMPLA(IU)
      AB(KOMP)=AFAK(KFAK)*ABKOF(INDEX)
      
      GO TO 13
   12 AB(KOMP)=0.
   13 KFAK=KFAK+1
      IU=IU+NLB
      
      end do
C
C        COMPONENTS WITH T-DEPENDENT ABSORPTION COEFFICIENTS
      DO19 KOMP=KOMPS,NKOMP

      NOP=NOFAK(IFAK)

      IF(NOP.EQ.0)GO TO 17
      IF(KOMPLA(IU).LE.0)GO TO 17
   15 INDEX=NPLATS(IFAK)-1+KOMPLA(IU)

C        THE VECTOR NPLATS IS DETERMINED BY SUBROUTINE TABS. IT GIVES THE ARRAY
C        INDEX OF THE TEMPERATURE AT WHICH THE INTERPOLATION IN ABKOF
C        BEGINS. NOFAK, GIVING INFORMATION ON THE T-INTERPOLATION AND
C        POSSIBLY INDICATING THAT AB=0 (NOFAK=0) IS ALSO DETERMINED BY TABS
C
C        INTERPOLATION
      DELSUM=0.
      DO16 NP=1,NOP

      DELSUM=DELSUM+AFAK(KFAK)*ABKOF(INDEX)
      KFAK=KFAK+1
   16 INDEX=INDEX+1
      

C
C        HAS THE INTERPOLATION BEEN MADE ON THE LOGARITHM
      IF(ILOGTA(KOMP).GT.0)DELSUM=EXP(DELSUM)
C        MULTIPLICATION BY FACTOR FROM JON-DETABS
      DELSUM=DELSUM*FAKTP(IFAK)
      IF(DELSUM.GE.0)GO TO 162
C
C        A NEGATIVE INTERPOLATION RESULT
  161 IF(NULL.GT.0)WRITE(IWRIT,200)KOMP,DELSUM,JP,ISET,T(NTP)
  200 FORMAT(4H AB(,I4,11H) NEGATIVE=,E12.4,5X,17HFOR WAVELENGTH NO,I5,
     *5X,6HSET NO,I5,5X,2HT=,F10.4,'  AND THEREFORE PUT =0 ***ABSKO***')
      AB(KOMP)=0.
      GO TO 18
  162 AB(KOMP)=DELSUM
      GO TO 18
   17 AB(KOMP)=0.
      KFAK=KFAK+NOP
   18 IU=IU+NLB
   19 IFAK=IFAK+1
C
C        WE MULTIPLY BY WAVELENGTH-DEPENDENT  FACTORS AND ADD UP. THIS IS
C        DONE IN DETABS.
      
      if (nlayer == -1 ) then
            CALL DETABS(J,JP,NTP,IOUTR,NTP)
      else 
            CALL DETABS(J,JP,NTP,IOUTR,nlayer)
      end if 
      
C
      IF(J.LE.0)GO TO 25
   24 ABSK(NTP)=SUMABS
      SPRID(NTP)=SUMSCA
      GO TO 26
C
C        SUMMATION TO GET A ROSSELAND MEAN
   25 CONTINUE
C use only central wavelengthinterval for Rosseland to avoid too
C large "abitrariness" on the value due to very small kap_nu in the
C wings of the planck function.
      if (xla(jp).le.5000. .or. xla(jp).ge.1.e5) go to 26
      IF(J.EQ.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/(SUMABS+SUMSCA)
      IF(J.LT.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/SUMABS
      SUMW(NTP)=SUMW(NTP)+ROSW(JP)*VIKTR
   26 CONTINUE
C
C        END OF 'THE SECOND NTP-LOOP'
C
      IF(J.GT.0)GO TO 29
      DO28 NTP=1,NT
      SPRID(NTP)=0.
   28 ABSK(NTP)=SUMW(NTP)/ABSK(NTP)
C
   29 CONTINUE
      RETURN
      END

************************************************************************
      SUBROUTINE AINV3(A,M)
      implicit real*8 (a-h,o-z)
C
C***** THIS SUBROUTINE EVALUATES THE INVERSE OF A
C***** SQUARE M*M MATRIX A
C***** THE DIMENSIONS OF THE ARRAYS MAY BE INTEGER VARIABLES WHEN USED**
C***** ON THE IBM 7094,BUT THEY MUST BE INTEGER CONSTANTS ON THE IBM 113
C
C
C      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(8,8),C(8),IND(8)
C
CC
C
  100 AMAX=0.0
      DO 2 I=1,M
      IND(I)=I
      IF(DABS(A(I,1))-AMAX)2,2,3
    3 AMAX=DABS(A(I,1))
      I4=I
    2 CONTINUE
      MM=M-1
      DO 111 J=1,MM
      IF(I4-J)6,6,4
    4 ISTO=IND(J)
      IND(J)=IND(I4)
      IND(I4)=ISTO
      DO 5 K=1,M
      STO=A(I4,K)
      A(I4,K)=A(J,K)
      A(J,K)=STO
    5 CONTINUE
    6 AMAX=0.0
      J1=J+1
      DO 11 I=J1,M
      A(I,J)=A(I,J)/A(J,J)
      DO 10 K=J1,M
      A(I,K)=A(I,K)-A(I,J)*A(J,K)
      IF (K-J1)14,14,10
   14 IF(DABS(A(I,K))-AMAX)10,10,17
   17 AMAX=DABS(A(I,K))
      I4=I
   10 CONTINUE
   11 CONTINUE
  111 CONTINUE
      DO 140 I1=1,MM
      I=M+1-I1
      I2=I-1
      DO 41 J1=1,I2
      J=I2+1-J1
      J2=J+1
      W1=-A(I,J)
      IF(I2-J2)141,43,43
   43 DO 42 K=J2,I2
      W1=W1-A(K,J)*C(K)
   42 CONTINUE
  141 C(J)=W1
   41 CONTINUE
      DO 40 K=1,I2
      A(I,K)=C(K)
   40 CONTINUE
  140 CONTINUE
      DO 150 I1=1,M
      I=M+1-I1
      I2=I+1
      W=A(I,I)
      DO 56 J=1,M
      IF (I-J)52,53,54
   52 W1=0.0
      GO TO 55
   53 W1=1.0
      GO TO 55
   54 W1=A(I,J)
   55 IF(I1-1)156,156,57
   57 DO 58 K=I2,M
      W1=W1-A(I,K)*A(K,J)
   58 CONTINUE
  156 C(J)=W1
   56 CONTINUE
      DO 50 J=1,M
      A(I,J)=C(J)/W
   50 CONTINUE
  150 CONTINUE
C        DENNA RUTIN LOESER TAANSPORTEKVATIONEN MED FEAUTRIERS METOD,
      DO 60 I=1,M
   63 IF(IND(I)-I)61,60,61
   61 J=IND(I)
      DO 62 K=1,M
      STO=A(K,I)
      A(K,I)=A(K,J)
      A(K,J)=STO
   62 CONTINUE
      ISTO=IND(J)
      IND(J)=J
      IND(I)=ISTO
      GO TO 63
   60 CONTINUE
      RETURN
      END
**********************************************************************
      subroutine algebn(n)
      implicit real*8 (a-h,o-z)
*
*  algebn does the algebraic manipulations when solving the radiation
*  transoport+flux equation system.  The routine is called once per
*  wavelength point.  The results are subtracted from the matrices O
*  and Q, which describe the dependence of the flux and radiation
*  pressure equations on the mean intensities.  The matrices are
*  arranged as follows:
*
*  a2(1)  a3(1)                               b2(1)  b3(1)
*  a1(2)  a2(2)  a3(2)                        b1(2)  b2(2)  b3(2)
*         a1(3)  a2(3)  a3(3)                        b1(3)  b2(3)  b3(3)
*                  .................                         ............
*                                a1(n)  a2(n)
*
*
*  d2(1)
*  d1(2)  d2(2)
*         d1(3)  d2(3)
*                ..................
*                                d1(n)  d2(n)
*
*  c(1)
*         c(2)
*                c(3)
*                 .................
*                                c(n)
*  i.e.,
*
*                   A*dJ + B *dT + E *dPe = 0
*  sum over wavel:  C*dJ + O1*dT + O2*dPe = 0
*  sum over wavel:  D*dJ + Q1*dT + Q2*dPe = 0
*
*  where A and B are tridiagonal, D is bi-diagonal, C is diagonal,
*  and the O and Q matrices are full.
*
*  n is the rank of the matrices.  Operation count 17 n**2 flops.
*
*  Modifierad nov-72/aake.
*  Subtracts d*ainverse*b from q and c*ainverse*b from o, where c is a
*  diagonal matrix and o is a full matrix.
*
*  Modified oct-78/aake.
*  b is now tridiagonal.
*
*  Modified june-79/aake.
*  Following stein (priv. comm. -79), a2 contains the row sum of
*  elements instead of the diagonal eelement.  This eliminates the
*  need for double precision.
*
*  Modified may-89/aake.
*  Bug in the inversion of A fixed.  Syntax cleaned up to show the
*  function more clearly (fortran 77 syntax).
*
      include 'parameter.inc'
c
      dimension a1(ndp),a2(ndp),a3(ndp),d1(ndp),d2(ndp),c(ndp),h(ndp),
     * b1(ndp),b2(ndp),b3(ndp),g1(ndp),q1(ndp,ndp),o1(ndp,ndp),
     * e1(ndp),e2(ndp),e3(ndp),g2(ndp),q2(ndp,ndp),o2(ndp,ndp)
      common/space1/a1,a2,a3,d1,d2,b1,b2,b3,c,q1,o1,
     * e1,e2,e3,q2,o2
*
*  a1(k) stores the fraction of row k-1 to subtract from row k
*  a2(k) stores the normalization factor for row k
*  a3(k) stores the fraction of (normalized) row k+1 to subtract
*        from (normalized) row k
*
      do 1 k=1,n-1
*
*  scaling factor to apply to row k
*
        f=1./(a2(k)-a3(k))
*
*  fraction of row k to add to row k+1
*
        a1(k+1)=a1(k+1)*f
*
*  the new rowsum in row k+1
*
        a2(k+1)=a2(k+1)-a1(k+1)*a2(k)
*
*  scale super-diagonal, save scale factor
*
        a3(k)=a3(k)*f
        a2(k)=f
    1 continue
*
*  last scaling factor
*
      a2(n)=1./a2(n)
*
      do 7 l=1,n
*
*  We first initialize the l'th column vector in Ainv*B.  At this point
*  we also normalize by multiplying each row with a2(k).
*
      do 3 k=1,l-2
        g1(k)=0.
        g2(k)=0.
    3 continue        
*
      if (l.eq.1) then
        g1(l  )  =b2(l  )
        g2(l  )  =e2(l  )
      else
        g1(l-1)=b3(l-1)
        g2(l-1)=e3(l-1)
        g1(l  )  =b2(l  )-a1(l  )*g1(l-1)
        g2(l  )  =e2(l  )-a1(l  )*g2(l-1)
      endif
      if (l.lt.n) then
        g1(l+1)=b1(l+1)-a1(l+1)*g1(l  )
        g2(l+1)=e1(l+1)-a1(l+1)*g2(l  )
      endif
*
*  Eliminate the sub-diagonal in A by adding a fraction of row k-1 to
*  row k.  [2 flops average].
*
      do 4 k=l+2,n
        g1(k)=-a1(k)*g1(k-1)
        g2(k)=-a1(k)*g2(k-1)
    4 continue
*
*  Eliminate the super-diagonal by adding a fraction of (normailzed)
*  row k+1 to (normalized) row k.  [3 flops average].
*
      g1(n)=a2(n)*g1(n)
      g2(n)=a2(n)*g2(n)
      do 6 k=n-1,1,-1
        g1(k)=a2(k)*g1(k)-a3(k)*g1(k+1)
        g2(k)=a2(k)*g2(k)-a3(k)*g2(k+1)
    6 continue
*
*  Subtract the results from the Q and O matrices.
*  [12 flops]
*
      o1(1,l)=o1(1,l)-c(1)*g1(1)
      o2(1,l)=o2(1,l)-c(1)*g2(1)
      q1(1,l)=q1(1,l)-d2(1)*g1(1)
      q2(1,l)=q2(1,l)-d2(1)*g2(1)
      do 7 k=2,n
        o1(k,l)=o1(k,l)-c(k)*g1(k)
        o2(k,l)=o2(k,l)-c(k)*g2(k)
        q1(k,l)=q1(k,l)-d1(k)*g1(k-1)-d2(k)*g1(k)
        q2(k,l)=q2(k,l)-d1(k)*g2(k-1)-d2(k)*g2(k)
    7 continue
*
      return
      end
**********************************************************************
      SUBROUTINE ARCHIV(LUN)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE STORES ALL INTERESTING INFORMATION ON A MODEL ON
C        FORTRAN FILE IARCH. MOREOVER, IT PUNCES CARDS FOR BELL'S USE.
C
      include 'parameter.inc'
      
C
      CHARACTER*1  DAY(10)
      CHARACTER*10 DAG,KLOCK
      DIMENSION ABSKA(NDP),SPRIDA(NDP),ABSKTR(20),SPRTR(20)
      DIMENSION PRESMP(33)     !345??? was 33 ??? (UGJ/Sep.10.98)
      EQUIVALENCE (DAG,DAY(3))
C
C        COMMONS SHARED BY MAIN PROGRAM
C        THESE SHOULD BE MODIFIED BY NORDLUND FOR HIS PURPOSES.
      common/carciv/ larciv    !=1 if called from arciv, otherwise = 0
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON/STATEC/PRAD(NDP),PTURB(NDP),P(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &              VV(NDP),FLUXC(NDP),PE(NDP),T(NDP),
     &              TAUDUM(NDP),RO(NDP),NTAU,ITMAX
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      COMMON /CTEFF/TEFF,FLUX
      COMMON/CG/G,KONSG
      COMMON /MIXC/PALFA,PBETA,PNY,PY
      COMMON /CSTYR/MIHAL,NOCONV
      COMMON/CLINE4/ILINE
      COMMON/CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &             IARCH
      COMMON/CFIL/IRESET(10),IDUM3,IDUM4
      COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/CXLSET/XL(20,10),IDUM6,NL(10)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE  
      COMMON /CTAUM/TAUM
C
C        COMMON SHARED BY SOLVE
      COMMON /CARC2/TKORRM(NDP),FCCORR(NDP),FLUXME(NWL),TAU5(NDP),INORD
C
C        COMMONS SHARED BY JON
      COMMON/CI5/abmarcs(18,ndp),ANJON(18,5),H(5),PART(18,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,NJ(16),IEL(16),
     *SUMM(ndp),NEL
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CMOL2/PPKDUM(33),NMOL
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C        COMMON SHARED BY DETABS
      COMMON/CARC4/PROV(30),NPROVA,NPROVS,NPROV
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      common /cmtest/pgm1(ndp),pgm2(ndp),pem1(ndp),pem2(ndp)
     *     ,tm1(ndp),tm2(ndp),pgos(ndp)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      COMMON /CMETPE/ PPEL(NDP), METPE
      common /ggchemmu/ggmu(NDP),ggrho(NDP),ppsum(ndp),ppappsum(ndp),
     &   ppnonappsum(ndp),tg(ndp),pges(ndp)
     &  ,ppat1sum(ndp),ppat2sum(ndp),ppmolsum(ndp),ppgs(ndp)
c        COMMON SHARED BY GGCHEM
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      character atnames*2, molnames*8
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)

      common /ggchemdetabs / f1_dt(ndp), f5_dt(ndp), 
     >                       rC(ndp), rMg(ndp), rAl(ndp), 
     >                       rSi(ndp), rHe(ndp), ro_dt(ndp)
      
      CHARACTER*8 SOURCE,ABNAME
      INTEGER MOLH, JUMP
               
C
c      print *, "archiv called"
      larciv = 1
      IARCH=LUN
      ISTAN2=1
      JSTAN2=NL(1)+1
      
      DO 60 K=1,JTAU
60    FCCORR(K)=3.14159*FCCORR(K)
C
C        WHICH PROGRAM WAS USED AND WHEN
C
C NOT SUPPORTED ON APOLLO      CALL ISODAT(DAY,KLOCK)
      WRITE(IARCH)INORD,DAG,KLOCK
C
C        STORE MODEL PARAMETERS
      WRITE(IARCH)TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(abmarcs(I,1),I=1,NEL)
      WRITE(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
C
C        STORE LAST TEMPERATURE-CORRECTION ARRAY
      WRITE(IARCH)JTAU,(TKORRM(K),K=1,JTAU),(FCCORR(K),K=1,JTAU)
C
C        COMPUTE AND STORE THERMODYNAMIC QUANTITIES AND DEPTHSCALES.
C        PUNCH BELL'S CARDS.
C
      

C
C
      DO 20 K=1,JTAU
      KL=K
C
C        JUST PRELIMINARY WE CHANGE TO FULL MOLECULAR EQUILIBRIUM. THIS
C        COMPUTED FOR DEMONSTRATION - IT IS NOT CONSISTENT WITH THE PRES
C        DENSITIES ETC GIVEN IN THE LISTING (WHERE ONLY H2 IS CONSIDERED
C        PORTANT). THE TAU(STANDARD) VALUES ARE NOT QUITE CONSISTENT, HO
C        THIS IS CERTAINLY QUITE UNIMPORTANT FOR MOST MODELS.
      MOLHO=MOLH
      MOLH=0
      
      if(metpe.eq.1) then
      CALL ABSKO(1,1,T(K),PE(K),ISTAN2,JSTAN2,ABSKA(1),SPRIDA(1),k)
      else if(metpe.eq.2) then
      if (T(k).gt.2000.) then
      CALL ABSKO(1,1,T(K),PE(K),ISTAN2,JSTAN2,ABSKA(1),SPRIDA(1),k)
      else
      CALL ABSKO(1,1,T(K),PPEL(K),ISTAN2,JSTAN2,ABSKA(1),SPRIDA(1),k)
      endif
      end if

      if(k.eq.1 .or. k.eq.jtau) write(6,2001) k,t(k),pe(k),ppel(k)
2001  format(' in archiv; k,t,pe,ppel = ',i3,f8.3,2e12.3)
****************18.12.04 Ch.H
*
* if JUMP >= 1 routine MOL dosn't work => no presmo
*    JUMP=1: Tsuji-routins are working instaed of MOL => partryck for molecules
*                                                xmettryck fot atoms
*    JUMP=2: JFF routines after same principle as Tsuji jump=1 routines
*    JUMP=3: JFF Gibs minimerings routines => call gem_package
*    JUMP=4: Added GGchem code by ERC
* ==>> more explanations in routine JON !!!
*
******

C This following piece doesn't really make sense if JUMP .ne.0 !!!
C -- but on the other hand probably it doesn't make harm, but be careful with the 
C  meaning of presmp  -- study later!!! (UGJ 25/1/01):

      IF (JUMP.GE.1) THEN 
 100  DO 101 I=1,33          !345
C  also in Tsujis routine is the number of molecules called NMOL
C  for no confusing I didn't use this name here
       PRESMP(I)=MAX(PARTRYCK(KL,I),1.D-99)
 101  CONTINUE
      ELSE 
C  NMOL=33 is number of molecules considered in the old MARCS chem.equilibrium
      DO 10 I=1,NMOL
       PRESMP(I)=MAX(PRESMO(I),1.D-99)
10    CONTINUE
      ENDIF
     
      HNIP=HNIC
C        CHANGE BACK AGAIN
C UGJ981018      MOLH=MOLHO
      CALL TERMO(k,T(K),PE(K),PRAD(K),PTOT,RRO,CP,CV,AGRAD,Q,U2)
      FORE=(ABSKA(1)+SPRIDA(1))/XKAPR(K)
      FURE=1./(XKAPR(K)*RRO)
      IF(K.GT.1)GO TO 11
      TAUS=FORE*TAU(1)
      Z=0.
      GO TO 12
   11 TAUS=TAUS+(TAU(K)-TAU(K-1))*(FORE+FOREM)*0.5
      Z=Z+(TAU(K)-TAU(K-1))*(FURE+FUREM)*0.5
   12 FOREM=FORE
      FUREM=FURE
C
      K1=MIN0(K+1,JTAU)
      PG=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
      EMU=(1.38*RRO*T(K))/(1.67E-8*PG)
      if (U2.lt.0.) then
         !write(7,*) ' in Archiv for K=',K,' U2= ',U2
         U2 = -U2
      end if
      U=SQRT(U2)
C
C INTERPOLATE THE CONVECTIVE FLUX AND VELOCITY TO THE DEPTH POINTS
C LOGARITHMICALLY
      IF(K.GT.1) GO TO 13
C K=1
      FCONV=0.
      V=0.
      GO TO 15
13    IF(K.EQ.JTAU) GO TO 14
C K"1,<JTAU
      YA=(TAU(K)-TAU(K-1))/(TAU(K+1)-TAU(K-1))
      YB=1.-YA
      FCONV=YA*FLUXC(K+1)+YB*FLUXC(K)
      IF(FLUXC(K).GT.0..AND.FLUXC(K+1).GT.0.) FCONV=
     &   EXP(YA*log(FLUXC(K+1))+YB*log(FLUXC(K)))
      V=YA*VV(K+1)+YB*VV(K)
      IF(VV(K).GT.0..AND.VV(K+1).GT.0.) V=
     &   EXP(YA*log(VV(K+1))+YB*log(VV(K)))
      GO TO 15
14    CONTINUE
C K=JTAU
      YA=(2.*TAU(K)-TAU(K-1)-TAU(K-2))/(TAU(K)-TAU(K-2))
      YB=1.-YA
      FCONV=YA*FLUXC(K)+YB*FLUXC(K-1)
      IF(FCONV.GT.FLUX) FCONV=FLUX
      V=YA*VV(K)+YB*VV(K-1)
15    CONTINUE
      ANCONV=FCONV/FLUX
C
C partial pressure of He is put into presmo(17)

        XNHE = abmarcs(2,k) / (XMH*XMY(k)) 
        PRESMO(17) = 1.38053d-16 * T(K) * XNHE * RRO
        PRESMP(17)=MAX(PRESMO(17),1.D-99)
C      PG=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
C      EMU=(1.38*RRO*T(K))/(1.67E-8*PG)

      WRITE(IARCH)K,TAU(K),TAUS,Z,T(K),PE(K),PG,PRAD(K),PTURB(K),
     &            XKAPR(K),RRO,EMU,CP,CV,AGRAD,Q,U,V,ANCONV,HNIP,
     &            NMOL,(PRESMP(I),I=1,NMOL)
   20 CONTINUE
C
      do  k=1,ndp
            K1=MIN0(K+1,ndp)
            PGm2(k)=P(K)-PRAD(K)-0.5*(PTURB(K)+PTURB(K1))
            pem2(k)=pe(k)
            tm2(k)=t(k)
      end do
C
C        STORE TYPICAL IONIZATION EQUILIBRIA AND ABSORPTION COEFFICIENTS
C
      NLP=NL(1)+1
      WRITE(IARCH)(NJ(I),I=1,NEL),NLP,(XL(J,1),J=1,NLP),
     &   NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
      IF (NLP.GE.21) WRITE(7,*) ' *** ERROR: NLP>20 ****'
C
CUGJ STORE ABSORBTIONKOEFFICIENT, SCATTERING AND LAMBDA FOR SPECTRUM 
CUGJ CALCULATIONS
C
      DO 25 K=1,JTAU
      DO 27 J=1,NLP
      if(metpe.eq.1) then
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1),k)
      else if(metpe.eq.2) then
      if (T(k).gt.2000.) then
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1),k)
      else
      CALL ABSKO(1,1,T(K),PPEL(K),1,J,ABSKA(1),SPRIDA(1),k)
      endif 
      end if
      
      ABSKTR(J)=ABSKA(1)
27    SPRTR(J)=SPRIDA(1)
      WRITE(IARCH) (ABSKTR(J),J=1,NLP)
      WRITE(IARCH) (SPRTR(J),J=1,NLP)
25    CONTINUE
C
      DO40 K=1,JTAU
      KL=K
      TAUK=log10(TAU(K))+10.01
      KTAU=TAUK
      DO32 J=1,NLP
      if(metpe.eq.1) then
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1),k)
      else if(metpe.eq.2) then
      if (T(k).gt.2000.) then
      CALL ABSKO(1,1,T(K),PE(K),1,J,ABSKA(1),SPRIDA(1),k)
      else
      CALL ABSKO(1,1,T(K),PPEL(K),1,J,ABSKA(1),SPRIDA(1),k)
      endif
      end if
      IF(J.GT.1)GO TO 31
      DO30 I=1,NEL
      NJP=NJ(I)
      WRITE(IARCH)K,TAU(K),T(K),PE(K),IEL(I),abmarcs(I,1),
     &    (ANJON(I,JJ),JJ=1,NJP),(PART(I,JJ),JJ=1,NJP)
   30 CONTINUE
   31 WRITE(IARCH)K,TAU(K),(PROV(KP),KP=1,NPROV),ABSKA(1),SPRIDA(1)
   32 CONTINUE
   40 CONTINUE
C
C        STORE FLUXES
C
      WRITE(IARCH)NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C
      OPEN(UNIT=33,FILE='FLUX.DAT',STATUS='unknown')
      DO 324 J=1,NLB
      FLUXME(J)=3.14159*FLUXME(J)
      WRITE(33,333) XLB(J),FLUXME(J)
  324 CONTINUE
      CLOSE(33)
333   FORMAT(1P2E15.5)

      RETURN
      END
C
      FUNCTION BPL(T,X)
      implicit real*8 (a-h,o-z)
C
C T*X must be greater than 1.7e6 to give BPL > 1.e-37
C BPL therefore limited to be > 1.e-30. UGJ 900510
C (generally changed to > 1.e-20 due to conv. problems. UGJ 961230)
C
      COMMON /BPLC/EX,X5
      DATA CP/1.191E27/,C2/1.438E8/  !CP is 2hc2 and C2 hc/k

      X5=((X**2)**2)*(X/CP)
      EX=EXP(-C2/(T*X))
      BPL=EX/((1.-EX)*X5)
      BPL = MAX(1.0D-99,BPL)
      RETURN
C
      ENTRY DIVBP(T,X)
      X6=X5*X
      TEX=T*(1.-EX)
      BPL=C2*(EX/TEX)/(TEX*X6)
      BPL = MAX(1.0D-99,BPL)
      RETURN
      END
C     MARK 4 RELEASE NAG COPYRIGHT 1974
C     MARK 4.5 REVISED
C     THIS ROUTINE ATTEMPTS TO SOLVE A REAL POLYNOMIAL EQUATION
C     HAVING N COEFFICIENTS (DEGREE  EQUALS  N-1) USING THE SEARCH
C     ALGORITHM PROPOSED IN GRANT AND HITCHINS (1971) TO
C     LIMITING MACHINE PRECISION.  ON ENTRY THE COEFFICIENTS
C     OF THE POLYNOMIAL ARE HELD IN THE ARRAY A(N), WITH A(0)
C     HOLDING THE COEFFICIENT OF THE HIGHEST POWER.  ON NORMAL
C     ENTRY THE PARAMETER IFAIL HAS VALUE 0 (HARD FAIL) OR 1
C     (SOFT FAIL) AND WILL BE ZERO ON SUCCESFUL EXIT WITH
C     THE CALCULATED ESTIMATES OF THE ROOTS HELD AS
C     REZ(K)+IIMZ(K), K EQUALS 1(1)N-1, IN APPROXIMATE DECREASING
C     ORDER OF MODULUS.  THE VALUE OF TOL IS OBTAINED BY
C     CALLING THE NAG ROUTINE X02AAF.
C     ABNORMAL EXITS WILL BE INDICATED BY IFAIL HAVING
C     VALUE 1 OR 2.  THE FORMER IMPLIES THAT EITHER A(1) EQUALS 0
C     OR N.LT.2 OR N.GT.100.  FOR IFAIL  EQUALS  2, A POSSIBLE
C     SADDLE
C     POINT HAS BEEN DETECTED.  THE NUMBER OF COEFFICIENTS
C     OF THE REDUCED POLYNOMIAL IS STORED IN N AND ITS
C     COEFFICIENTS ARE STORED IN A(1) TO A(N), THE ROOTS
C     THUS FAR BEING STORED IN THE ARRAYS REZ AND IMZ
C     STARTING WITH REZ(N)+IIMZ(N).  AN IMMEDIATE RE-ENTRY
C     IS POSSIBLE WITH IFAIL UNCHANGED AND WITH A NEW
C     STARTING POINT FOR THE SEARCH HELD IN REZ(1)+IIMZ(1).
C     REF - J.I.M.A., VOL.8., PP122-129 (1971).
CCCC
C
C
      SUBROUTINE C02AEF(A, N, REZ, IMZ, TOL, IFAIL)
      implicit real*8 (a-h,o-z)
C
C      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER IFAIL, IND, N, I, K, II, I2, JTEMP
      REAL*8  IMZ,J,JX,NFUN
      DIMENSION A(N),B(100),C(100),REZ(N),IMZ(N)
      LOGICAL SAT,FLAG
      COMMON /AC02AE/ X, Y, R, RX, J, JX, SAT
C     ETEXT/DTEXT
C     DATA ONE/1.0/,A1P5/1.5/,ZERO/0.0/,P4Z1/1.0E-5/
      DATA ONE /1.0/, A1P5 /1.5/, ZERO /0.0/, P4Z1 /1.0E-5/
C     DATA TWO/2.0/,P5/0.5/,P2Z1/1.0E-3/,P1/0.1/
      DATA TWO /2.0/, P5 /0.5/, P2Z1 /1.0E-3/, P1 /0.1/
C     DATA P3Z2/2.0E-4/,FOUR/4.0/
C     DATA P3Z2/2.0E-4/,FOUR/4.0/
      DATA P3Z2 /2.0E-4/, FOUR /4.0/
      XXX = X02AAF(XXX)
C     THE ABOVE TEST WAS ADDED AT 4.5 TO PREVENT TOL BEING TOO
C     SMALL

      IF (TOL.LT.XXX) TOL = XXX
      FAC = ONE
      FLAG = IFAIL.EQ.2
      IF (FLAG) IFAIL = 1
      IND = 0
      TOL2 = TOL**A1P5
      IF (A(1).NE.ZERO .AND. N.GE.2 .AND. N.LE.100) GO TO 20
      IND = 0
      GO TO 720
   20 IF (A(N).NE.0.0) GO TO 40
      REZ(N-1) = ZERO
      IMZ(N-1) = ZERO
      N = N - 1
      GO TO 20
   40 SCALE = ZERO
C     FUNCTION/DFUNCTION
      DO 60 I=1,N
       IF (DABS(A(I)).GE.P4Z1) SCALE = SCALE + DLOG(DABS(A(I)))
C     FUNCTION/DFUNCTION
   60 CONTINUE
      K = IDINT(SCALE/(DBLE(N)*DLOG(TWO))+P5)
      SCALE = TWO**(-K)
      DO 80 I=1,N
       A(I) = A(I)*SCALE
       B(I) = A(I)
C     TEST FOR LOW ORDER POLYNOMIAL FOR EXPLICIT SOLUTION
   80 CONTINUE
      IF (N.GT.3) GO TO 100
      GO TO (720, 560, 580), N
  100 DO 160 I=2,N
       II = N - I + 2
       DO 120 K=2,II
       I2 = II - K + 1
       C(K-1) = B(II)*B(K) - B(1)*B(I2)
  120  CONTINUE
       IF (C(II-1).LT.-TOL) GO TO 200
       T = ONE
       IF (C(II-1).GE.ONE) T = ONE/C(II-1)
       JTEMP = II - 1
       DO 140 K=1,JTEMP
       B(K) = C(K)*T
  140  CONTINUE
  160 CONTINUE
      FAC = FAC*TWO
      SCALE = ONE
      I = N
  180 I = I - 1
      IF (I.LT.1) GO TO 100
      SCALE = SCALE*TWO
      A(I) = A(I)*SCALE
      B(I) = A(I)
      GO TO 180
  200 IF (.NOT.FLAG) GO TO 220
      X = REZ(1)
      Y = IMZ(1) + TOL
      FLAG = .FALSE.
      GO TO 240
  220 X = P2Z1
      Y = P1
  240 CALL C02AEZ(A, N, TOL)
      FUN = R*R + J*J
  260 G = RX*RX + JX*JX
      IF (G.GE.FUN*TOL2) GO TO 300
      IND = 0
      SCALE = ONE
      I = N
  280 I = I - 1
      IF (I.LT.1) GO TO 720
      SCALE = SCALE*FAC
      A(I) = A(I)/SCALE
      GO TO 280
  300 S1 = -(R*RX+J*JX)/G
      S2 = (R*JX-J*RX)/G
C     FUNCTION/DFUNCTION
      SIG = P3Z2
      S = DSQRT(S1*S1+S2*S2)
      IF (S.LE.ONE) GO TO 320
      S1 = S1/S
      S2 = S2/S
C     VALID DIRECTION OF SEARCH HAS BEEN DETERMINED, NOW
C     PROCEED TO DETERMINE SUITABLE STEP
      SIG = SIG/S
  320 X = X + S1
      Y = Y + S2
  340 CALL C02AEZ(A, N, TOL)
      IF (SAT) GO TO 380
      NFUN = R*R + J*J
      IF (FUN-NFUN.GE.SIG*FUN) GO TO 360
      S1 = P5*S1
      S2 = P5*S2
      S = P5*S
      SIG = P5*SIG
      X = X - S1
      Y = Y - S2
      GO TO 340
  360 FUN = NFUN
      GO TO 260
  380 FUN = ONE/TOL2
      K = 0
C     FUNCTION/DFUNCTION
      IMZ(N-1) = Y*FAC
C     CHECK POSSIBILITY OF REAL ROOT
      IF (DABS(Y).GT.P1) GO TO 420
      S1 = Y
      Y = ZERO
      CALL C02AEZ(A, N, TOL)
      Y = S1
C     REAL ROOT ACCEPTED AND BOTH BACKWARD AND FORWARD DEFLATIONS
C     ARE PERFORMED WITH LINEAR FACTOR
      IF (.NOT.SAT) GO TO 420
      REZ(N-1) = X*FAC
      IMZ(N-1) = ZERO
      N = N - 1
      B(1) = A(1)
      C(N) = -A(N+1)/X
      DO 400 I=2,N
       B(I) = A(I) + X*B(I-1)
       II = N - I + 1
       C(II) = (C(II+1)-A(II+1))/X
  400 CONTINUE
C     COMPLEX ROOT ACCEPTED AND BOTH BACKWARD AND FORWARD
C     DEFLATIONS ARE PERFORMED WITH QUADRATIC FACTOR
      GO TO 460
  420 REZ(N-1) = X*FAC
      REZ(N-2) = X*FAC
      IMZ(N-2) = -IMZ(N-1)
      N = N - 2
      R = TWO*X
      J = -(X*X+Y*Y)
      B(1) = A(1)
      B(2) = A(2) + R*B(1)
      C(N) = -A(N+2)/J
      C(N-1) = -(A(N+1)+R*C(N))/J
      IF (N.EQ.2) GO TO 460
      DO 440 I=3,N
       B(I) = A(I) + R*B(I-1) + J*B(I-2)
       II = N - I + 1
       C(II) = -(A(II+2)-C(II+2)+R*C(II+1))/J
C     MATCHING POINT FOR COMPOSITE DEFLATION
  440 CONTINUE
C     FUNCTION/DFUNCTION
  460 DO 480 I=1,N
       NFUN = DABS(B(I)) + DABS(C(I))
C     FUNCTION/DFUNCTION
       IF (NFUN.LE.TOL) GO TO 480
       NFUN = DABS(B(I)-C(I))/NFUN
       IF (NFUN.GE.FUN) GO TO 480
       FUN = NFUN
       K = I
  480 CONTINUE
      IF (K.EQ.1) GO TO 520
      JTEMP = K - 1
      DO 500 I=1,JTEMP
       A(I) = B(I)
  500 CONTINUE
  520 A(K) = P5*(B(K)+C(K))
      IF (K.EQ.N) GO TO 40
      JTEMP = K + 1
      DO 540 I=JTEMP,N
       A(I) = C(I)
  540 CONTINUE
      GO TO 40
  560 REZ(1) = -A(2)/A(1)*FAC
      IMZ(1) = ZERO
      GO TO 700
  580 R = A(2)*A(2) - FOUR*A(1)*A(3)
      IF (R.GT.ZERO) GO TO 600
      REZ(2) = -P5*A(2)/A(1)*FAC
C     FUNCTION/DFUNCTION
      REZ(1) = REZ(2)
      IMZ(2) = P5*DSQRT(-R)/A(1)*FAC
      IMZ(1) = -IMZ(2)
      GO TO 700
  600 IMZ(1) = ZERO
      IMZ(2) = ZERO
C     FUNCTION/DFUNCTION
      IF (A(2)) 620, 640, 660
  620 REZ(1) = P5*(-A(2)+DSQRT(R))/A(1)*FAC
      GO TO 680
  640 REZ(1) = -P5*A(2)/A(1)*FAC
C     FUNCTION/DFUNCTION
      GO TO 680
  660 REZ(1) = P5*(-A(2)-DSQRT(R))/A(1)*FAC
  680 REZ(2) = A(3)/(REZ(1)*A(1))*FAC*FAC
  700 N = 1
  720 IFAIL = IND
      RETURN
      END
C
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     EVALUATES R,RX,J,JX AT THE POINT X+IY AND APPLIES THE ADAMS
C     TEST.
C     THE BOOLEAN VARIABLE SAT IS GIVEN THE VALUE TRUE IF THE TEST
C     IS
C     SATISFIED.
CCCC
C
C
      SUBROUTINE C02AEZ(A, N, TOL)
      implicit real*8 (a-h,o-z)
C
C      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER N, K
      REAL*8  J,JX
      DIMENSION A(N)
      LOGICAL SAT
C
C     ETEXT/DTEXT
      COMMON /AC02AE/ X, Y, R, RX, J, JX, SAT
      DATA TWO/2.0/,ZERO/0.0/,P8/0.8/,TEN /10.0/,A8/8.0/
C

      P = -TWO*X
C     FUNCTION/DFUNCTION
      Q = X*X + Y*Y
      T = SQRT(Q)
      A2 = ZERO
      B2 = ZERO
      B1 = A(1)
C     FUNCTION/DFUNCTION
      A1 = A(1)
      C = ABS(A1)*P8
      N = N - 2
      DO 20 K=2,N
       A3 = A2
       A2 = A1
C     FUNCTION/DFUNCTION
       A1 = A(K) - P*A2 - Q*A3
       C = T*C + ABS(A1)
       B3 = B2
       B2 = B1
       B1 = A1 - P*B2 - Q*B3
   20 CONTINUE
      N = N + 2
      A3 = A2
      A2 = A1
      A1 = A(N-1) - P*A2 - Q*A3
      R = A(N) + X*A1 - Q*A2
      J = A1*Y
      RX = A1 - TWO*B2*Y*Y
C     FUNCTION/DFUNCTION
      JX = TWO*Y*(B1-X*B2)
C     FUNCTION/DFUNCTION
      C = T*(T*C+ABS(A1)) + ABS(R)
      SAT = (SQRT(R*R+J*J)).LT.((TEN*C-A8*(ABS(R)+ABS(A1)*T)+TWO* ABS(X*
     *A1))*TOL)
      RETURN
      END
C
      SUBROUTINE CLOCK
      implicit real*8 (a-h,o-z)
C
C  TIME SINCE LAST CALL,ACCUMULATED EXECUTION TIME
C  AND TIME REMAINING OF REQUESTED CPU TIME (for Cyber machines ?)
C
C
C     DATA IS/0/
CC
CC
C     IF (IS.NE.1) THEN
C       CALL MSLEFT(MS)
C       IS=1
C       MPP=MS
C     END IF
CC
C     CALL MSLEFT(MP)
C     SEC=(MPP-MP)/1000.
C     ACC=(MS-MP)/1000.0
C     REM=MP/1000.0
C     MPP=MP
CCC
      SEC=0.
      ACC=0.
        REM=0.
CCC

      RETURN  
      END
C
C
C         
      subroutine cclock
      implicit real*8 (a-h,o-z)
      call clock
      END
C
      SUBROUTINE COSCAL(LUN)
      implicit real*8 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /STATEC/DUM1(9*NDP),TT(NDP),DUM2(NDP),RO(NDP),NTAU,ITER
      COMMON /CISPH/ISPH
      DIMENSION TP(NDP)
C
C READ OLD MODEL

      CALL OLDSTA
C
C READ SCALE MODEL
      LUN1=LUN+1
      READ(LUN1) (DUMS,K=1,360),(TP(K),K=1,NTAU)
C
C PRINT
      WRITE(7,50) LUN,LUN1
50    FORMAT('0OLD MODEL FROM LUN',I3,' COMBINED WITH TEMPERATURE '
     &,'FROM LUN',I3,' FOR TAU < 1.0')
C
C ASSEMBLE
      K1=24
      DO 100 K=1,K1
100   TT(K)=TP(K)
      K2=K1+1
      TTK2=TT(K2)
      DO 101 K=K2,NTAU
101   TT(K)=TT(K)+TP(K2)-TTK2
C
C INTEGRATE PRESSURE EQUATION
        if (isph.eq.1) then 
          CALL TRYCK_sph
        else 
          CALL TRYCK
        end if
      ITER=0
      RETURN
      END
C
      SUBROUTINE DETABS(J,JP,NTP,IOUTR,nlayer)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE GIVES THE DETAILS OF THE ABSORPTION MECHANISMS.
C        CHANGES IN THE ABSORPTION-COEFFICIENT PROGRAM ARE EXPECTED TO
C        BE CONFINED TO THE TABLES AND TO THIS ROUTINE.
C        DETABS HAS TWO PURPOSES
C        1. JP=0   DETERMINATION OF WAVELENGTH-INDEPENDENT FACTORS (DEP. ON
C                  T, PE AND THE COMPONENT) STORED IN FAKT.
C        2. JP= THE ACTUAL WAVELENGTH NUMBER.
C                  MULTIPLICATION OF AB, COMPUTED IN SUBROUTINE ABSKO,
C                  BY WAVELENGTH-DEPENDENT FACTORS. SUMMATION OF THE TOTAL
C                  ABSORPTION AND SCATTERING COEFFICIENTS ( SUMABS AND
C                  SUMSCA ).
C
C        N O T E .  BEFORE A CALL ON DETABS FOR PURPOSE 1, SUBROUTINE
C        JON MUST HAVE BEEN CALLED.
C
C        IF J IS LESS THAN OR EQUAL TO ZERO, THE WEIGHT FOR A ROSSELAND MEAN
C        WILL BE COMPUTED AND STORED IN VIKTR (THE WEIGHT BEING 1/VIKTR).
C        NTP IS THE ARRAY INDEX OF THE T-PE POINT.
C        IF IOUTR IS GREATER THAN ZERO AT A CALL WITH JP GREATER THAN ZERO
C        (PART TWO OF THE ROUTINE), DETAILS OF THE ABSORPTION COEFFICIENTS
C        ARE PRINTED. IF IOUTR IS GREATER THAN ONE, A TABLE HEADING IS ALSO
C        PRINTED.
C
C
C        CONTENTS OF COMMON/CI5/, COMMUNICATING PHYSICAL INFORMATION FROM
C        SUBROUTINE JON.
C             ABUND  ABUNDANCES
C             ANJON  FRACTIONS OF IONIZATION
C             H      QUANTUM NUMBER OF THE HIGHEST EXISTING HYDROGENIC LEVEL
C             PART   PARTITION FUNCTIONS
C             DXI    DECREASE OF IONIZATION ENERGY OF HYDROGEN IN ELECTRON-VOLTS
C             F1     N(HI)/N(H)
C             F2     N(HII)/N(H)
C             F3     N(H-)/N(H)
C             F4     N(H2+)/N(H)
C             F5     N(H2)/N(H)
C             XKHM   'DISSOCIATION CONSTANT' OF H-
C             XMH    MASS OF THE HYDROGEN ATOM IN GRAMS
C             XMY    GRAMS OF STELLAR MATTER/GRAMS OF HYDROGEN
C
C        DIMENSIONS NECESSARY
C        ABUND(NEL),ANJON(NEL,MAX(NJ)),ELS(NT),H(5),HREST(NT),PART(NEL,MAX(NJ)),
C        PROV(NKOMP)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS IN COMMON /CA5/ ARE
C        COMMENTED ON IN SUBROUTINE ABSKO.
C        NEL IS THE NUMBER OF CHEMICAL ELEMENTS INITIATED IN SUBROUTINE INJON
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I
C        NKOMP IS THE NUMBER OF COMPONENTS, NOT INCLUDING THOSE ADDED BY
C             ANALYTICAL EXPRESSIONS AFTER STATEMENT NO. 13.
C        NT   IS THE NUMBER OF TEMPERATURES-ELECTRON PRESSURES GIVEN AT THE
C             CALL OF SUBROUTINE ABSKO.
C
C
      include 'parameter.inc'
C
      DIMENSION ELS(NDP),HREST(NDP)
      DIMENSION FAKRAY(NDP)
      DIMENSION PHTVA(NDP),PHEL(NDP),H2RAY(NDP)
      COMMON /CLIN/lin_cia
      COMMON/CI5/abmarcs(18,ndp),ANJON(18,5),H(5),PART(18,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CA2/RCA2DUM(4000),ICA2DUM(602),NKOMP
      COMMON/CA5/AB(30),FAKT(30),PE(NDP),T(NDP),XLA(20),XLA3(20),RO,
     *SUMABS,SUMSCA,VIKTR,ISET,NLB
      COMMON/UTPUT/IREAD,IWRIT
      COMMON /CARC4/ PROV(30),NPROVA,NPROVS,NPROV
      CHARACTER*8 SOURCE,ABNAME
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /ggchembool/ iggcall
      common /ggchemdetabs / f1_dt(ndp), f5_dt(ndp), 
     >                       rC(ndp), rMg(ndp), rAl(ndp), 
     >                       rSi(ndp), rHe(ndp), ro_dt(ndp)
      LOGICAL FIRST
      CHARACTER*8 NHMIN,NH2PR,NHEPR,NELS,NHRAY,NH2RAY
      DATA FIRST/.TRUE./
      DATA NH2PR/'H2PR'/,NHEPR/'HEPR'/,NELS/'ELSC'/,NHRAY/'H-SC'/
     & ,NH2RAY/'H2SC'/,NHMIN/'H-'/

C
C        SAVE ABSORPTION COMPONENT NAMES THE FIRST TIME DETABS IS CALLED


      if (first) then
            do komp=19, nkomp
                  abname(1) = nhmin 
                  abname(komp-16) = abname(komp)
            end do
            nprova=nkomp-16+2
            abname(nprova-1)=nh2pr
            abname(nprova)=nhepr
            nprovs=3
            nprov=nprova+nprovs
            abname(nprov-2)=nels
            abname(nprov-1)=nhray
            abname(nprov)=nh2ray
            ncall=0
            first = .false.

      end if
      ncall=ncall+1
      teta = 5040./t(ntp)
      
      if (jp==0) then
      !iggcall =0
      if (iggcall == 0 ) then
      
C     1. COMPUTATION OF WAVELENGTH-INDEPENDENT QUANTITIES
      hn=1./(xmh*xmy(ntp))
      hnh=f1*hn
      
C        H-
      FAKT(1)=PE(NTP)*HNH*1.E-17/XKHM !bound free H-
      FAKT(18)=PE(NTP)*HNH*2.E-26/PART(1,1) !free-free H-
      
C        HI
      TETA31=31.30364*TETA
      XFAKH=2.0898E-26/PART(1,1)*HNH
      NNIV=15
      XNIV=15.
      IF(H(1).LT.XNIV)NNIV=INT(H(1))
      DO M=1,NNIV
      XM2=M*M
      XM3=XM2*DFLOAT(M)
      FAKT(M+1)=XFAKH*EXP(-TETA31*(1.-1./XM2))/XM3
      end do
      FAKT(NNIV+1)=FAKT(NNIV+1)*MIN(H(1)-NNIV,1.0D+0)
      IF(NNIV.LT.15) then
      N1=NNIV+1
      DO M=N1,15
       FAKT(M+1)=0.
      end do
C
C        FREE-FREE HI ABSORPTION
      end if
      UMC=2.3026*DXI*TETA
      EXPJ=XFAKH*EXP(-TETA31+UMC)/(2.*TETA31)
      ADDF=EXP(TETA31/((DFLOAT(NNIV)+0.5)**2)-UMC)-1.
      IF(H(1).LT.XNIV+0.5)ADDF=0.
      FAKT(17)=EXPJ
      HREST(NTP)=EXPJ*ADDF
            
C        H+H
      FAKT(19)=(HNH*1.E-25)*(HNH*1.E-25)*RO
           
                  
            
C        H2+
      FAKT(20)=(HNH*1.E-20)**2*RO*ANJON(1,2)/ANJON(1,1)
      
            
C        H2-
      FAKT(21)=PE(NTP)*F5*HN
            
C        C I
      FAKT(22)=ANJON(3,1)*abmarcs(3,ntp)*HN*9./PART(3,1)
      
                  
                  
C        MG I
      FAKT(23)=ANJON(8,1)*abmarcs(8,ntp)*HN/PART(8,1)
      
            
C        AL I
      FAKT(24)=ANJON(9,1)*abmarcs(9,ntp)*HN*6./PART(9,1)
      
            
C        SI I
      FAKT(25)=ANJON(10,1)*abmarcs(10,ntp)*HN*9./PART(10,1)
      
            
C        HE I
      FAKT(26)=ANJON(2,1)*abmarcs(2,ntp)*HN/PART(2,1)
       
C        HE-
      FAKT(27)=PE(NTP)*ANJON(2,1)*abmarcs(2,ntp)*HN
      
                  
                  
C        ELECTRON SCATTERING
      ELS(NTP)=4.8206E-9*PE(NTP)/(T(NTP)*RO)
            
      
C        RAYLEIGH SCATTERING
      FAKRAY(NTP)=HNH*2./PART(1,1)
      H2RAY(NTP)=F5*HN
                  
            
      RETURN
      else
      
C     1. COMPUTATION OF WAVELENGTH-INDEPENDENT QUANTITIES
        hn=1./(xmh*xmy(ntp))
        !hnh=f1*hn
        
        hnh=f1_dt(nlayer)*hn
C        H-
      FAKT(1)=PE(NTP)*HNH*1.E-17/XKHM !bound free H-
      FAKT(18)=PE(NTP)*HNH*2.E-26/PART(1,1) !free-free H-
      
      
C        HI
      TETA31=31.30364*TETA
      XFAKH=2.0898E-26/PART(1,1)*HNH
      NNIV=15
      XNIV=15.
      IF(H(1).LT.XNIV)NNIV=INT(H(1))
      DO M=1,NNIV
      XM2=M*M
      XM3=XM2*DFLOAT(M)
       FAKT(M+1)=XFAKH*EXP(-TETA31*(1.-1./XM2))/XM3
      end do
      FAKT(NNIV+1)=FAKT(NNIV+1)*MIN(H(1)-NNIV,1.0D+0)
      IF(NNIV.LT.15) then
      N1=NNIV+1
      DO M=N1,15
       FAKT(M+1)=0.
      end do
      end if
C
C        FREE-FREE HI ABSORPTION
      UMC=2.3026*DXI*TETA
      EXPJ=XFAKH*EXP(-TETA31+UMC)/(2.*TETA31)
      ADDF=EXP(TETA31/((DFLOAT(NNIV)+0.5)**2)-UMC)-1.
      IF(H(1).LT.XNIV+0.5)ADDF=0.
      FAKT(17)=EXPJ
      HREST(NTP)=EXPJ*ADDF

C        H+H
      
      FAKT(19)=(HNH*1.E-25)*(HNH*1.E-25)*ro_dt(nlayer)


C        H2+
      FAKT(20)=(HNH*1.E-20)**2*ro_dt(nlayer)*ANJON(1,2)/ANJON(1,1)
      
C        H2-
       FAKT(21)=PE(NTP)*f5_dt(nlayer)*HN

C        C I
      FAKT(22)=rC(nlayer)*abmarcs(3,ntp)*HN

C        MG I
      FAKT(23)=rMg(nlayer)*abmarcs(8,ntp)*HN
    
      
C        AL I
      !FAKT(24)=ANJON(9,1)*abmarcs(9,ntp)*HN*6./PART(9,1)
      FAKT(24)=rAl(nlayer)*abmarcs(9,ntp)*HN
      !print*, "Si ", fakt(25)
      
      
C        HE I
      !FAKT(26)=ANJON(2,1)*abmarcs(2,ntp)*HN/PART(2,1)
      FAKT(26)=rHe(nlayer)*abmarcs(2,ntp)*HN
      !print*, "He ", fakt(26)
     
C        HE-
      !FAKT(27)=PE(NTP)*ANJON(2,1)*abmarcs(2,ntp)*HN
      FAKT(27)=PE(NTP)*abmarcs(2,ntp)*HN
      !print*, "He- ", fakt(27)

           
C        ELECTRON SCATTERING
      !ELS(NTP)=4.8206E-9*PE(NTP)/(T(NTP)*RO)
      ELS(NTP)=4.8206E-9*PE(NTP)/(T(NTP)*ro_dt(nlayer))

C        RAYLEIGH SCATTERING
      FAKRAY(NTP)=HNH*2./PART(1,1)
      !H2RAY(NTP)=F5*HN
      h2ray(ntp) = f5_dt(nlayer)*hn

      RETURN
      end if
C        N O T E . APART FROM VECTORS HREST AND ELS, NONE OF THE
C        TEMPERATURE- OR PRESSURE-DEPENDENT VARIABLES DEFINED ABOVE CAN
C        GENERALLY BE USED AT THE NEXT VISIT BELOW.
C        ANY SET OF FACTORS WHICH IS WANTED SHOULD BE STORED IN AN ARRAY WITH
C        DIMENSION = NT, LIKE HREST AND ELS, OR IN FAKT, WHERE THE DATA FOR
C        FURTHER USE ARE STORED IN SUBR. ABSKO.
C
C        2. WAVELENGTH-DEPENDENT FACTORS. SUMMATION.
C        CORRECTION FOR STIMULATED EMISSION
      else 
      EXPA=EXP(-28556.*TETA/XLA(JP))

   11 STIM=1.-EXPA
C
C        ABSORPTION
      SUMABS=0.
C        H I
      DO12 KOMP=2,17
      SUMABS=SUMABS+AB(KOMP)

   12 CONTINUE
     
      
      SUMABS=(SUMABS+HREST(NTP))*XLA3(JP)

      PROV(2)=SUMABS
C        H-
      HMIN=AB(1)+AB(18)/STIM
      SUMABS=SUMABS+HMIN
      PROV(1)=HMIN
C        H+H, H2+, HE I, C I, MG I, AL I, SI I
      DO13 KOMP=19,NKOMP
      SUMABS=SUMABS+AB(KOMP)
      PROV(KOMP-16)=AB(KOMP)
   13 CONTINUE

      SUMABS=SUMABS*STIM
    
C        SCATTERING
      XRAY=MAX(XLA(JP),1026.0D+0)
      XRAY2=1./(XRAY*XRAY)
      RAYH=XRAY2*XRAY2*(5.799E-13+XRAY2*(1.422E-6+XRAY2*2.784))*
     *FAKRAY(NTP)
      RAYH2=XRAY2*XRAY2*(8.14E-13+XRAY2*(1.28E-6+XRAY2*1.61))*H2RAY(NTP)
      SUMSCA=ELS(NTP)+RAYH+RAYH2

      

      PROV(NPROV-2)=ELS(NTP)
      PROV(NPROV-1)=RAYH
      PROV(NPROV)=RAYH2
C
      IF(J.GT.0)GO TO 15
C
C        WEIGHT FOR A ROSSELAND MEAN
   14 VIKTR=EXPA/(STIM*STIM*(XLA3(JP)*1E-3)**2)
   15 CONTINUE
C
      IF(IOUTR-1)23,21,20
C
C        **** PRINT-OUT ****
   20 WRITE(IWRIT,200)XLA(JP),(ABNAME(KP),KP=1,NPROV)
  200 FORMAT(' WAVEL.=',F7.0,'    ABS       SCAT  ',6A6,/10A6)
   21 DO22 KP=1,NPROVA
   22 PROV(KP)=PROV(KP)/SUMABS*STIM
      DO24 KP=1,NPROVS
   24 PROV(NPROVA+KP)=PROV(NPROVA+KP)/SUMSCA
      if(lin_cia.ne.1) 
     *   WRITE(IWRIT,*)
     * ' (Linskys H2-H2 and H2-He CIA are not included in SUMABS)'
      WRITE(IWRIT,201)T(NTP),SUMABS,SUMSCA,(PROV(KP),KP=1,NPROV)
  201 FORMAT(' T=',F7.1,1X,1p2E10.3,0p6F7.4,//10F7.4)
   23 CONTINUE
      RETURN
      end if
      END
C
      SUBROUTINE DUBINT(NXSKAL,XSKAL,NYSKAL,YSKAL,IYBEG,IYEND,N,X,Y,
     &                                            FACTOR,IX,IY1,IY2)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES INTERPOLATION FACTORS FOR INTERPOLATION I
C        NOT NECESSARILY EQUIDISTANT, TWO DIMENSIONAL TABLE. SCALES
C           XSKAL (NXSKAL POINTS)
C           YSKAL(NYSKAL POINTS)
C        THE TABLE IS ONLY DEFINED FOR Y-VALUES YSKAL(L), WHERE L IS
C        WITHIN THE INTERVAL IYBEG(NX) TO IYEND(NX) FOR A GIVEN XSKAL(NX
C        ARGUMENTS ARE GIVEN IN X AND Y (N POINTS).
C        RESULTING FACTOR FOR POINT K IS PUT IN FACTOR(K,1-2,1-2)
C        STARTING POINTS AT INTERPOLATION IN IX(K), REFERRING TO THE XSC
C        IN IY1(K) AND IY2(K), REFERRING TO THE RESTRICTED Y SCALES DEFI
C        XSCALE(IX(K)) AND XSCALE(IX(K)+1), RESPECTIVELY.
C        INTERPOLATIONS AND EXTRAPOLATIONS ARE  L I N E A R .
C
      include 'parameter.inc'
C
      DIMENSION XSKAL(30),YSKAL(30),X(NDP),Y(NDP),FACTOR(NDP,2,2),
     &          IYBEG(30),IYEND(30),IY1(NDP),IY2(NDP),IX(NDP)
C
      DO5 K=1,N
      DO1 J=2,NXSKAL
      JMEM1=J
      IF(X(K).LT.XSKAL(J))GO TO 2
    1 CONTINUE
    2 IX(K)=JMEM1-1
      JM1=JMEM1-1
      DO3 J=2,NYSKAL
      JMEM2=J
      IF(Y(K).LT.YSKAL(J))GO TO 4
    3 CONTINUE
    4 IY=JMEM2-1
      IY1(K)=IY+1-IYBEG(JM1)
      IY1(K)=MIN0(IY1(K),IYEND(JM1)-IYBEG(JM1))
      IY1(K)=MAX0(IY1(K),1)
      IY2(K)=IY+1-IYBEG(JMEM1)
      IY2(K)=MIN0(IY2(K),IYEND(JMEM1)-IYBEG(JMEM1))
      IY2(K)=MAX0(IY2(K),1)
      I1=IY1(K)+IYBEG(JM1)-1
      I2=IY2(K)+IYBEG(JMEM1)-1
      DX=(X(K)-XSKAL(JMEM1-1))/(XSKAL(JMEM1)-XSKAL(JMEM1-1))
      DY1=(Y(K)-YSKAL(I1))/(YSKAL(I1+1)-YSKAL(I1))
      DY2=(Y(K)-YSKAL(I2))/(YSKAL(I2+1)-YSKAL(I2))
      FACTOR(K,1,1)=(1.-DX-DY1+DX*DY1)
      FACTOR(K,2,1)=(1.-DY2)*DX
      FACTOR(K,1,2)=(1.-DX)*DY1
    5 FACTOR(K,2,2)=DX*DY2
      RETURN
      END
C
      SUBROUTINE DUMIN
      implicit real*8 (a-h,o-z)
C
C 'DUMIN' READS THE FIRST SETS OF CARDS IN THE JONABS-DATA. THESE DATA
C ARE USED IN ATMOS BUT NOT IN MARCS.  *NORD*
C
      COMMON /UTPUT/IREAD,IWRIT
C
      READ(IREAD,50)I,J
      N=4+I/8+J/16+2*(I/16)
      DO 100 I=1,N
100   READ(IREAD,51)A
      DO 101 J=1,3
      READ(IREAD,50)N
      DO 101 I=1,N
101   READ(IREAD,51)A
      RETURN
50    FORMAT(2I5)
51    FORMAT(A4)
      END
C
      FUNCTION FOUR(Y,X,K,N)
      implicit real*8 (a-h,o-z)
C
      include 'parameter.inc'
C
      DIMENSION X(NDP),Y(NDP)
C
C FOURPOINT LAGRANGE INTERPOLATION TO FIND Y(K-.5). Y AND X OF LENGTH N.
C
C START ADRESS AND NONCENTERING
      IF(K.EQ.2) GOTO 3
      KK=MIN0(MAX0(K-3,0),N-4)
      II=0
      IF(K.LE.2) II=-1
      IF(K.EQ.N) II=1
      XX=.5*( X(KK+II+2)+X(KK+II+3))
C
C I LOOP
      FOUR=0.
      DO 1 I=1,4
      PROD=Y(KK+I)
C
C J LOOP
      DO 2 J=1,4
      IF(J.EQ.I) GO TO 2
      PROD=PROD*(XX-X(KK+J))/(X(KK+I)-X(KK+J))
2     CONTINUE
C
1     FOUR=FOUR+PROD
      RETURN
C
C LINEAR INTERPOLATION AT K=2. 780605.
3     FOUR=.5*(Y(1)+Y(2))
      RETURN
      END
C
      SUBROUTINE GAUSI(K,A,B,AI,XMYI)
      implicit real*8 (a-h,o-z)
C
C        RUTINEN GER VIKTER OCH INTEGRATIONSPUNKTER FOER GAUSSINTEGRATIO
C        MELLAN A OCH B - B AER OEVRE GRAENS , A NEDRE. KAELLA FOER DATA
C        LOWAN, DAVIDS, LEVENSON,  BULL AMER MATH SOC  48 SID 739  (1942
C        AI=VIKTER, XMYI=INTEGRATIONSPUNKTER.
C        INTEGRATIONSORDNING K.  K VAELJES MELLAN 2 OCH 10.
C
      DIMENSION AI(K),XMYI(K),AP(29),XMYP(29),INDOV(9)
      DOUBLE PRECISION AP,XMYP
C               10 DATAKORT FOER AP, 9 FOER XMYP OCH 1 FOER INDOV
      DATA AP/1.0,0.55555555555555,.88888888888888,.347854845137
     *,0.65214515486254,0.23692688505618,0.47862867049936,
     * 0.56888888888888,0.17132449237917,0.36076157304813,
     * 0.46791393457269,0.12948496616887,0.27970539148927,
     * 0.38183005050511,0.41795918367346,0.10122853629037,
     * 0.22238103445337,0.31370664587788,0.36268378337836,
     * 0.08127438836157,0.18064816069485,0.26061069640293,
     * 0.31234707704000,0.33023935500126,0.06667134430868,
     * 0.14945134915058,0.21908636251598,0.26926671930999,
     * 0.29552422471475/
      DATA XMYP/
     *0.57735026918962,.77459666924148,.0,0.86113631159405,
     *0.33998104358485,.90617984593866,.53846931010568,.0,
     *0.93246951420315,.66120938646626,.23861918608319,
     *0.94910791234275,.74153118559939,.40584515137739,.0,
     *0.96028985649753,.79666647741362,.52553240991632,
     *0.18343464249565,.96816023950762,.83603110732663,
     *0.61337143270059,.32425342340380,.0,0.97390652851717,
     *0.86506336668898,.67940956829902,.43339539412924,
     *0.14887433898163/
      DATA INDOV/1,3,5,8,11,15,19,24,29/
      IF(K.EQ.1)GO TO 7
      KUD=0
      FLK=DFLOAT(K)/2.
      K2=K/2
      FK=DFLOAT(K2)
      IF(ABS(FLK-FK)-1.E-7)2,1,1
    1 K2=K2+1
      KUD=1
    2 IOEV=INDOV(K-1)
      INED=IOEV-K2
      DO3 I=1,K2
      IP=INED+I
      XMYI(I)=-XMYP(IP)*(B-A)*0.5+(B+A)*0.5
    3 AI(I)=(B-A)*0.5*AP(IP)
      K2=K2+1
      DO4 I=K2,K
      IP=IOEV+K2-I
      IF(KUD)6,6,5
    5 IP=IP-1
    6 CONTINUE
      XMYI(I)= XMYP(IP)*(B-A)*0.5+(B+A)*0.5
    4 AI(I)=(B-A)*0.5*AP(IP)
      RETURN
    7 XMYI(1)=(B+A)*0.5
      AI(1)=B-A
      RETURN
      END
C


      SUBROUTINE INABS(IOUTS)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE  READS ABSORPTION COEFFICIENT TABLES AND INTER/EXTRA-
C        POLATES THEM TO OUR WAVELENGTHS GIVEN IN XL. THE INTERPOLATION IS
C        PERFORMED SEPARATELY FOR EACH WAVELENGTH SET.
C
C        NKOMP IS THE NUMBER OF COMPONENTS IN THE FULL TABLE.
C        NEXTL SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED ON EXTRA-
C              POLATION IN WAVELENGTH,
C        NUTZL IF PRINT-OUT IS WANTED WHEN WE PUT THE COEFFICIENT =0 OUTSIDE THE
C              WAVELENGTH REGION OF THE TABLES.
C        NEXTT AND NUTZL ARE THE CORRESPONDING QUANTITIES ON INTERPOLATION IN
C              T, MADE IN SUBROUTINE TABS.
C        NULL  SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED (FROM SUB-
C              ROUTINE ABSKO) WHEN A COEFFICIENT IS FOUND TO BE LESS THAN ZERO
C              ON INTERPOLATION IN T AND THEREFORE PUT EQUAL TO ZERO.
C
C        FOR EACH COMPONENT THE FOLLOWING PARAMETERS MUST BE SPECIFIED
C        ABNAME IS THE NAME OF, OR A SYMBOL FOR, THE ABSORPTION MECHANISM.
C        SOURCE INDICATES THE SOURCE OR REFERENCE OF THE DATA
C
C        1. PARAMETERS FOR THE WAVELENGTH INTERPOLATION.
C          ILOGL SHOULD BE GREATER THAN ZERO IF INTERPOLATION IN WAVELENGTH IS
C                TO BE PERFORMED ON THE LOGARITHMIC ABSORPTION COEFFICIENTS
C                (WITH SUBSEQUENT EXPONENTIATION OF THE RESULTS - HERE IF ILOGT
C                IS EQUAL TO ZERO OR IN SUBROUTINE ABSKO IF ILOGT IS GREATER
C                THAN ZERO). OTHERWISE INTERPOLATION IN WAVELENGTH IS MADE
C                DIRECTLY ON THE ABSORPTION COEFFICIENTS THEMSELVES.
C          KVADL SHOULD BE GREATER THAN ZERO IF QUADRATIC INTERPOLATION IN
C                WAVELENGTH IS WANTED. OTHERWISE INTERPOLATION WILL BE LINEAR
C          MINEX SHOULD BE GT 0 IF LINEAR EXTRAPOLATION (INSTEAD OF PUTTING THE
C                COEFFICIENT = 0) IS WANTED TOWARDS SHORTER WAVELENGTHS.
C          MAXEX, CORRESPONDING TOWARDS LONGER WAVELENGTHS.
C          NLATB IS THE NUMBER OF WAVELENGTH POINTS OF THE ABSORPTION COEFFI-
C                CIENT TABLE TO BE READ.
C          XLATB ARE THOSE WAVELENGTHS. THEY SHOULD BE GIVEN IN INCREASING ORDER
C
C        2. PARAMETERS FOR THE TEMPERATURE INTERPOLATION.
C          ILOGT, KVADT, MINET, MAXET AND NTETB ARE THE T-INTERPOLATION
C                ANlogUES TO ILOGL-NLATB.
C          ITETA IS PUT GREATER THAN ZERO WHEN TETA VALUES (TETA=5040./T) ARE
C                GIVEN IN XTET INSTEAD OF TEMPERATURES.
C          XTET ARE THE TEMPERATURE (TETA) VALUES OF THE ABSORPTION
C                COEFFICIENT TABLE TO BE READ. THE XTET VALUES SHOULD BE GIVEN
C                IN INCREASING ORDER AND EQUIDISTANTLY, HOWEVER (IELMAX-1)
C                CHANGES OF THE INTERVAL ARE ALLOWED. THE PROGRAM CHECKS  THAT
C                THIS NUMBER IS NOT EXCEEDED.
C        XKAP IS THE ABSORPTION COEFFICIENT TABLE FOR THE ACTUAL COMPONENT. THE
C                WAVELENGTHS INCREASES MORE RAPIDLY THAN T (TETA).
C
C        THE TABLES FOR  T E M P E R A T U R E - I N D E P E N D E N T
C        C O M P O N E N T S  S H O U L D  B E  P U T  F I R S T .
C           THE RESULTING TABLE IS PUT IN ABKOF. HERE T (TETA) INCREASES MORE
C        RAPIDLY THAN XLA, WHICH INCREASES MORE RAPIDLY THAN KOMP. IF THE RESULT
C        OF THE INTERPOLATION IS ZERO FOR A CERTAIN XLA(J) AND KOMP, THIS IS NOT
C        PUT IN ABKOF. INSTEAD A NOTE IS MADE IN KOMPLA (KOMPLA(NLB*(KOMP-
C        1)+J) IS PUT EQUAL TO ZERO). OTHERWISE THE KOMPLA VALUE TELLS WHERE IN
C        ABKOF THE TABLE FOR THE COMPONENT KOMP AND THE WAVELENGTH J BEGINS.
C
C        A DETAILED PRINT-OUT IS GIVEN IF IOUTS IS GREATER THAN ZERO.
C
C
C        DIMENSIONS NECESSARY
C        ABKOF(NABDIM),ABNAME(NKOMP),DELT(NKOMP,IELMAX),IDEL(NKOMP),
C        IDISKV(MAX(NLATB)),ILOGTA(NKOMP),IRESET(NSET),ISVIT(NKOMP),ITETA(NKOMP)
C        KOMPLA(MAX(NL)*NKOMP),KVADT(NKOMP),MAXET(NKOMP),MINET(NKOMP),
C        NL(NSET),NTAET(NKOMP),NTM(NKOMP,IELMAX),SOURCE(NKOMP),
C        TBOLT(NKOMP,IELMAX),XKAP(MAX(NLATB),MAX(NTETB)),XL(MAX(NL),NSET)
C        XLA(MAX(NL)),XLA3(MAX(NL)),XLATB(MAX(NLATB)),XTET(MAX(NTETB)),
C        XTETP(MAX(NTETB))
C
C        THE DIMENSIONS ARE LOWER LIMITS
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS (GIVEN BELOW) IN
C               ANY ABSORPTION COEFFICIENT TABLE.
C        NABDIM IS THE DIMENSION OF THE ABKOF ARRAY (GIVEN BELOW).
C        NKOMP IS THE NUMBER OF 'COMPONENTS', I.E. EQUAL TO THE NUMBER OF
C               DIFFERENT ABSORPTION COEFFICIENT TABLES TO BE READ.
C        NL(I)  IS THE NUMBER OF WAVELENGTHS IN THE WAVELENGTH SET I.
C        NLATB(KOMP) IS THE NUMBER OF WAVELENGTH POINTS IN THE TABLE TO BE READ
C               FOR THE COMPONENT KOMP.
C        NSET   IS THE NUMBER OF WAVELENGTH SETS.
C        NTETB  IS THE NUMBER OF TEMPERATURE POINTS IN THE TABLE FOR THE COM-
C               PONENT BEING CONSIDERED.
C
C
      DIMENSION IDISKV(40),XLATB(40),XTET(30),NTAET(30),XKAP(40,30),
     *XLA3(20),XLA(20),XTETP(30)
      COMMON/CARC4/PROV(30),NDUM(3)
      COMMON /CHAR/ ABNAME(30),SOURCE(30)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(30,2),TBOT(30,2),IDEL(30),ISVIT(30),ITETA(30),
     *KVADT(30),MAXET(30),MINET(30),NTM(30,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA3/ILOGTA(30),NULL
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      CHARACTER*8 ABNAME,SOURCE
C
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS IN THE XKAP-
C        TABLE. THE DIMENSIONS OF TBOT, DELT AND NTM ARE AFFECTED BY THIS NUMBER
      IELMAX=2
C        THE DIMENSION OF THE ABKOF ARRAY
      NABDIM=4000
      DO705 L=1,30
  705 XTETP(L)=0.
      
C
      READ(IREAT,101)NKOMP,NEXTL,NUTZL,NEXTT,NUTZT,NULL
C
C
      KOMPR=0
      REWIND ISLASK
C
C        LOOP OVER COMPONENTS STARTS (THE 'FIRST KOMP-LOOP')

      DO720 KOMP=1,NKOMP
      READ(IREAT,105)ABNAME(KOMP),SOURCE(KOMP)
      READ(IREAT,102)ILOGL,KVADL,MINEX,MAXEX,NLATB
      READ(IREAT,103)(XLATB(J),J=1,NLATB)

C
C        WE FIND THE DISCONTINUITIES IN WAVELENGTH
C        A DISCONTINUITY IN A TABLE IS DEFINED BY TWO WAVELENGTH POINTS
C        WITHIN LESS THAN TWO ANGSTROEMS.
      IDISK=0
      IDISKV(1)=0
      DO700 J=2,NLATB
      IDISKV(J)=0
      IF((XLATB(J)-XLATB(J-1)).GE.2.)GO TO 700

      IDISKV(J-1)=1
      IDISKV(J)=1
      IDISK=1
  700 CONTINUE
C
C        CONTINUE READING
      READ(IREAT,102)ILOGT,KVADT(KOMP),MINET(KOMP),MAXET(KOMP),NTETB,
     * ITETA(KOMP)

      ILOGTA(KOMP)=ILOGT
      IF(NTETB.GT.1)GO TO 702
  701 KOMPR=KOMPR+1
      GO TO 703
  702 READ(IREAT,103)(XTET(L),L=1,NTETB)
C
C        FINALLY THE ABSORPTION COEFFICIENT TABLE IS READ
  703 DO 704 K=1,NTETB
  704 READ(IREAT,104)(XKAP(JJ,K),JJ=1,NLATB)
C
C        WE TAKE THE LOGARITHMS BEFORE THE WAVELENGTH INTERPOLATION
C        IF ILOGL IS GREATER THAN ZERO.
      IF(ILOGL.LT.1)GO TO 712
  710 DO 711 K=1,NTETB
      DO 711 JJ=1,NLATB
      IF(XKAP(JJ,K).GT.0.)GO TO 711
C
C        A COEFFICIENT FOR WHICH THE LOGARITHM SHOULD BE TAKEN IS ZERO
      
      XKAP(JJ,K)=1.E-20
  711 XKAP(JJ,K)=log(XKAP(JJ,K))
  712 CONTINUE
C
C        PREPARATION OF THE T-INTERPOLATION IN SUBROUTINE TABS
C
C        WE FIND OUT WHETHER ISVIT(KOMP) CAN BE CHOSEN GREATER THAN ZERO. THIS
C        IS THE CASE IF THE T SCALE AND MAXET, MINET AND KVADT ARE IDENTICAL
C        WITH THOSE OF THE PREVIOUS COMPONENT. IF ISVIT IS GREATER THAN ZERO
C        THE TIME SPENT IN SUBR. TABS WILL BE DECREASED.
      ISVIT(KOMP)=0
      IF(NTETB.LE.1)GO TO 719
      DO 721 L=1,NTETB
      IF(XTET(L).NE.XTETP(L))GO TO 722
  721 CONTINUE
      IF(NTETB.NE.NTETBP)GO TO 722
      IF(MAXET(KOMP).NE.MAXETP) GO TO 722
      IF(MINET(KOMP).NE.MINETP) GO TO 722
      IF(KVADT(KOMP).NE.KVADTP)GO TO 722
      ISVIT(KOMP)=1
  722 CONTINUE
C
C        WE REMEMBER TEMPERATURES ETC. FOR NEXT COMPONENT
      DO723 L=1,NTETB
  723 XTETP(L)=XTET(L)
      NTETBP=NTETB
      MAXETP=MAXET(KOMP)
      MINETP=MINET(KOMP)
      KVADTP=KVADT(KOMP)
C
C        WE FIND THE INTERVALS IN THE T (TETA) SCALE
      TBOT(KOMP,1)=XTET(1)
      DELT(KOMP,1)=XTET(2)-XTET(1)
      NTM(KOMP,1)=1
      IDEL(KOMP)=1
      IF(NTETB.EQ.2)GO TO 719
C
      J=1
      LF=1
      DO714 L=3,NTETB
      DIFF=XTET(L)-XTET(L-1)
      IF(ABS(1.-DIFF/DELT(KOMP,J)).LT.1.E-4)GO TO 714
      J=J+1
      IF(J.GT.IELMAX)GO TO 715
      TBOT(KOMP,J)=XTET(L-1)
      DELT(KOMP,J)=DIFF
      NTM(KOMP,J-1)=LF
      LF=0
  714 LF=LF+1
      NTM(KOMP,J)=LF

      IDEL(KOMP)=J
      GO TO 719
C        TOO MANY DIFFERENT INTERVALS IN THE T-TABLE FOR THIS COMPONENT
  715 STOP 'INABS 1'
C
  719 NTAET(KOMP)=NTETB
C        ALL DATA NECESSARY BELOW FOR THIS COMPONENT ARE STORED ON UNIT
C        ISLASK
      WRITE(ISLASK)KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),J=1,
     *NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,(XTET(L),L=1,NTETB)
     *,((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
C
  720 CONTINUE
C        END OF 'THE FIRST KOMP-LOOP'
C
      KOMPS=KOMPR+1
C
C
C        WE BUILD THE ABKOF ARRAY. INTERPOLATION IN WAVELENGTH.
C
C        LOOP OVER WAVELENGTH SETS ('THE ISET-LOOP')
      DO 70 ISET=1,NSET
      REWIND ISLASK
      NLB=NL(ISET)
      DO1 J=1,NLB
      XLA(J)=XL(J,ISET)
    1 XLA3(J)=XLA(J)**3
      INDEX=1
C
C        LOOP OVER COMPONENTS STARTS ('THE SECOND KOMP-LOOP')
      DO60 KOMP=1,NKOMP
      READ(ISLASK)KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),J=1,
     *NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,(XTET(L),L=1,NTETB)
     *,((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
      JI=1
      LAMBI=1
C
C        LOOP OVER WAVELENGTHS ('THE J-LOOP') STARTS
      DO60 J=1,NLB
C        SEARCHING IN WAVELENGTH
      IU=NLB*(KOMP-1)+J
      KOMPLA(IU)=INDEX
      DO24 JJ=1,NLATB
      IHELP=JJ
      IF(XLA(J)-XLATB(JJ))25,24,24
   24 LAMBI=JJ
   25 CONTINUE
      IF(IHELP-1)45,45,26
   26 IF(KVADL)33,33,27
   33 IF(NLATB-LAMBI-1)41,31,31
   27 IF(NLATB-LAMBI-1)41,28,29
C
C        QUADRATIC INTERPOLATION
   28 LAMBI=LAMBI-1
   29 CONTINUE
C        ARE DISCONTINUITIES PRESENT
      IF(IDISK.LE.0)GO TO 299
      IF(IDISKV(LAMBI+1).LE.0)GO TO 299
      IF(XLA(J).GT.XLATB(LAMBI+1))GO TO 292
  291 IF(IDISKV(LAMBI).GT.0)GO TO 31
      IF(LAMBI.EQ.1)GO TO 31
      LAMBI=LAMBI-1
      GO TO 299
  292 LAMBI=LAMBI+1
      IF(IDISKV(LAMBI+1).GT.0)GO TO 31
      IF(LAMBI+1.EQ.NLATB)GO TO 31
  299 CONTINUE

      DXX1=XLA(J)-XLATB(LAMBI)
      DXX2=XLA(J)-XLATB(LAMBI+1)
      DXX3=XLA(J)-XLATB(LAMBI+2)
      DX21=XLATB(LAMBI+1)-XLATB(LAMBI)
      DX32=XLATB(LAMBI+2)-XLATB(LAMBI+1)
      DX31=XLATB(LAMBI+2)-XLATB(LAMBI)
      A1=DXX2*DXX3/(DX21*DX31)
      A2=DXX1*DXX3/(DX21*DX32)
      A3=DXX1*DXX2/(DX31*DX32)
     
C
      DO30 K=1,NTETB
      
      ABKOF(INDEX)=A1*XKAP(LAMBI,K)-A2*XKAP(LAMBI+1,K)+A3*
     &XKAP(LAMBI+2,K)
      

   30 INDEX=INDEX+1
      GO TO 59
C
C        LINEAR INTER- AND EXTRAPOLATION
   31 A2=(XLA(J)-XLATB(LAMBI))/(XLATB(LAMBI+1)-XLATB(LAMBI))
      A1=1.-A2
      DO32 K=1,NTETB
   
      ABKOF(INDEX)=A1*XKAP(LAMBI,K)+A2*XKAP(LAMBI+1,K)
   32 INDEX=INDEX+1
      GO TO 59
C
C        TOO GREAT A WAVELENGTH - OUTSIDE THE TABLE
   41 IF(MAXEX)50,50,42
   42 LAMBI=LAMBI-1
      GO TO 31
C
C        TOO SMALL A WAVELENGTH - OUTSIDE THE TABLE
   45 IF(MINEX)50,50,46
   46 GO TO 31
C
C        ABS. COEFF. IS PUT = ZERO
   50 KOMPLA(IU)=0
      
      GO TO 60
C
   59 IF(ILOGL.LT.1)GO TO 592
      IF(ILOGT.GT.0)GO TO 60
C
C        LOGARITHMIC INTERPOLATION ONLY IN WAVELENGTH
      LIP=INDEX-NTETB
      LAP=INDEX-1
      DO 591 LL=LIP,LAP
  591 ABKOF(LL)=EXP(ABKOF(LL))
C
  592 CONTINUE
C
      IF(ILOGT.LE.0)GO TO 60
C        WE TAKE THE LOGARITHM BEFORE THE T INTERPOLATION IF ILOGT GT 0
      LIP=INDEX-NTETB
      LAP=INDEX-1
      DO593 LL=LIP,LAP
      IF(ABKOF(LL).GT.0.)GO TO 593
C
C        IMPOSSIBLE TO TAKE THE LOGARITHM OF A NEGATIVE COEFFICIENT
      LUS=LL-LIP+1
     
      ABKOF(LL)=1.E-20
  593 ABKOF(LL)=log(ABKOF(LL))
   60 CONTINUE
C        END OF 'THE J-LOOP'
C        END OF 'THE SECOND KOMP-LOOP'
C
C        WRITE THE DATA OF THE SET ISET ON UNIT IRESET(ISET)
      NABKOF=INDEX-1
      NKOMPL=IU
      IREADP=IRESET(ISET)
      WRITE(IREADP)ISET,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
C
      END FILE IREADP
      BACKSPACE IREADP
C
C        CHECK DIMENSION OF ABKOF
     
      IF(NABKOF.LE.NABDIM)GO TO 70
C        TOO SMALL DIMENSION FOR ABKOF
      
      STOP 'INABS 2'
   70 CONTINUE
C

      DO 71 ISET=1,NSET
      IREADP=IRESET(ISET)
      REWIND IREADP
71    CONTINUE
      IF(IOUTS.LE.0) GOTO 74
C

      DO73 M=1,NSET
      NP=NL(M)
   73 CONTINUE
   74 CONTINUE
  101 FORMAT(8X,I2,5(9X,I1))
  102 FORMAT(4(9X,I1),8X,I2,9X,I1)
  103 FORMAT(6F10.0)
  104 FORMAT(6E10.3)
  105 FORMAT(2A8)
      RETURN
      END
C  end of       SUBROUTINE INABS(IOUTS)
C

      ! begin ADS
C     Aaron's routines to read two dimensional OS files:

      MODULE OPAMOD
         implicit real*8 (a-h,o-z)
         private
         public :: get_opac, read_opac  
         include 'parameter.inc'  

         ! WARNING: -mcmodel=medium is needed to compile - array is too large! 
         double precision crossec_data(NOSPEC,NWL,NOPI,NOTI) 
         double precision crossec_pgrid(NOSPEC,NOPI)
         double precision crossec_tgrid(NOSPEC,NOTI)
         integer crossec_ktemp(NOSPEC)
         integer crossec_kpres(NOSPEC)        
         character(len=4) crossec_molid(NOSPEC)
         
         contains
         subroutine get_opac(press, temp, spec_i, opa_out, struc_len)
            ! io
            integer :: spec_i, struc_len
            double precision, intent(in) :: press(struc_len)
            double precision, intent(in) :: temp(struc_len) 
            double precision, intent(inout) :: opa_out(NWL, struc_len)
   
            ! local
            integer :: ktemp, kpres
            double precision, allocatable :: pgrid(:), tgrid(:)
            double precision, allocatable :: opa_in(:, :, :)
   
            ktemp = crossec_ktemp(spec_i)
            kpres = crossec_kpres(spec_i)

            if (allocated(pgrid)) deallocate(pgrid)
            if (allocated(tgrid)) deallocate(tgrid)
            if (allocated(opa_in)) deallocate(opa_in)
            allocate(pgrid(kpres))
            allocate(tgrid(ktemp))
            allocate(opa_in(NWL,kpres,ktemp))
            pgrid(1:kpres) = crossec_pgrid(spec_i, 1:kpres)
            tgrid(1:ktemp) = crossec_tgrid(spec_i, 1:ktemp)
            opa_in(1:NWL, 1:kpres, 1:ktemp) = 
     &         crossec_data(spec_i, 1:NWL, 1:kpres, 1:ktemp)

            ! NOTE: pressure in opacity array is currently in bar !
            call interpol_kappa(press, temp, pgrid, tgrid, 
     &         opa_in, opa_out, 
     &         struc_len, NWL, kpres, ktemp)
            
         end subroutine get_opac

         subroutine read_opac(file_b, spec_i, marcs_wn_grid)   
           ! io
           double precision, intent(in) :: marcs_wn_grid(NWL)
           integer, intent(in) :: spec_i
           character(len=*), intent(in) :: file_b
   
           ! local
           character(len=:), allocatable :: file_nml
           character(len=:), allocatable :: file_data
           character(len=:), allocatable :: file_wnos
   
           character(len=4)          :: molid
   
           integer                   :: kpres, ktemp, nwnos, lwriteos 
           integer                   :: p_i, t_i, freq_i
           integer                   :: temporary_index
           
           double precision  :: pmol_read(100*NOPI)
           double precision  :: tmol_read(100*NOTI)
           double precision :: wn
   
           double precision, allocatable :: cread(:,:,:)
           double precision, allocatable :: wnmol(:)

           character(len=100) :: filename_crossec_out
C       -----------------------------------
C       Set the filenames of the input data
C       -----------------------------------
           file_nml = trim(file_b) // "input.nml"
           file_data = trim(file_b) // "crossec.dat"
           file_wnos = trim(file_b) // "wn.dat"

C       --------------------
C       Read in the namelist
C       --------------------
           !write(*,*) file_nml,file_data,file_wnos
           call read_opac_namelist(file_nml, NOPI, NOTI, kpres, ktemp, 
     &        molid, nwnos, pmol_read, tmol_read)
           crossec_kpres(spec_i) = kpres
           crossec_ktemp(spec_i) = ktemp
           crossec_molid(spec_i) = molid  
           crossec_pgrid(spec_i,1:kpres) = pmol_read(1:kpres)
           crossec_tgrid(spec_i,1:ktemp) = tmol_read(1:ktemp)
C       ----------------------------------------------
C       Allocate the crosssec and the wavelength array
C         and readin the crossections and the wavelength
C       ----------------------------------------------        
           if (ALLOCATED(cread)) THEN
             deallocate(cread)
           endif
           allocate(cread(nwnos, kpres, ktemp))
           if (ALLOCATED(wnmol)) THEN
             deallocate(wnmol)
           endif
           allocate(wnmol(nwnos))
   
           call read_opac_data(file_data, file_wnos, 
     &        nwnos, kpres, ktemp,
     &        cread, wnmol)
   
C       ----------------------------------------------
C       Load opac into array - match closest marcs wn
C       ----------------------------------------------
           do freq_i =1, NWL
             wn = marcs_wn_grid(freq_i)
             ! search for closest marcs wn
             call search_intp_ind(wnmol, nwnos, wn, temporary_index)
             ! assign data to common block:
             crossec_data(spec_i, freq_i, 1:kpres, 1:ktemp) = 
     &          cread(temporary_index,1:kpres,1:ktemp)
     
             ! set opacity to 0 outside range of input data
             if (wn<wnmol(1) .or. wn>wnmol(size(wnmol))) then
                crossec_data(spec_i, freq_i, 1:kpres,1:ktemp) = 0
             endif
            end do 

             lwriteos=0
             if (lwriteos .eq.0) go to 7778
             write(6,*) 'lwriteos = ',lwriteos    
            !START OF CROSSSEC OUT MODULE
            !module to write out crossec data in a plotable and visible format
            !spec_i has to be the number of the desired molecule in mol_names
                        !START OF CROSSSEC OUT MODULE
            !module to write out crossec data in a plotable and visible format
            !spec_i has to be the number of the desired molecule in mol_names
      !       if (spec_i.eq.71) then
 7777 format(A19,I2.2,A4) 
             write(filename_crossec_out,7777) "crosssec_write_out_"
     >        ,spec_i,'.dat'
             write(*,*) filename_crossec_out
             open(unit=777,file=filename_crossec_out)
            ! Write the fixed part of the block
             write(777,'(a)') ' &INPUTOSMOL'
             write(777,'(a)') ' MOLID   = ' // trim(molid) // '  ,'
             write(777,'(a,i12,a)') ' KTEMP   =', ktemp, ','
             write(777,'(a,4(1x,f18.12,a))')
     >         ' TMOL    =', 
     >         tmol_read(1), '     ,', 
     >         tmol_read(2), '     ,', 
     >         tmol_read(3), '     ,', 
     >         tmol_read(4), '     ,'
           write(777,'(5(1x,f18.11,a))')
     >         tmol_read(5), '     ,', 
     >         tmol_read(6), '     ,', 
     >         tmol_read(7), '     ,', 
     >         tmol_read(8), '     ,', 
     >         tmol_read(9), '     ,'
           write(777,'(3(1x,f18.11,a))')
     >         tmol_read(10), '     ,', 
     >         tmol_read(11), '     ,', 
     >         tmol_read(12), '     ,'
           write(777,'(a,i12,a)') ' NWNOS   =', nwnos, ','
           write(777,'(a,f20.14,a)') ' VKMS    =   ', 3.000d0, '     ,'
           write(777,'(a,i11,a)') ' KISO    =', 1, ','
          write(777,'(a,f20.14,a,E22.15,a)')
     >         ' RELISO  =   ', 1.000d0, '     , 14*', 0.00d0,'  ,'
           write(777,'(a,i11)') ' L_PER_STELLAR   =           0,'
           write(777,'(a,i11)') ' LCHROM  =           0'
           write(777,'(a)') ' /'
           do i=1,NWL
             !do j=1,kpres !comment this and the following if in to get an overview of the prid of the data
              !do k=1,ktemp !comment this and the following if in to get an overview of the tgrid of the data      
              !write(777,'(E10.3,A1,E14.6,A1,E10.3,A1,E10.3)') 
      !>          crossec_data(spec_i,i,1,k)," ",marcs_wn_grid(i)," ",
      !>          crossec_pgrid(spec_i,j), " ", crossec_tgrid(spec_i,k) 
              !enddo
             !enddo
             
             write(777, '(F10.3, 48E11.3)') !comment out when going through the whole pgrid and/or tgrid
     >         marcs_wn_grid(i), crossec_data(spec_i,i,1,1:ktemp)
           enddo
           !enddo
            close(777)
      !      stop !this module is designed to just run at the beginning of the simulation to just write out the crossec data
                  !you can also leave this stop commented out if you want to run afterwards
            !endif            
            !END OF CROSSSEC OUT 
7778     continue
         end subroutine read_opac
         
         subroutine read_opac_data(file_data, file_wnos, 
     &        nwnos, kpres, ktemp,
     &        crossec, wnmol)
   
           implicit none
   
           !! Reads Opacity data from given file.
           character(len=*), intent(in) :: file_data, file_wnos
           integer,  intent(in)    :: kpres, ktemp, nwnos
           double precision, intent(out) :: crossec(nwnos,kpres,ktemp)
           double precision, intent(out) :: wnmol(nwnos)
   
   
           ! local
           integer :: fu, rc, readin_size, freq_i, p_i, t_i
           integer :: readin_index
           double precision, allocatable :: crossec_readin(:)
   
C       ----------------
C       OPACITY DATA IO 
C       ----------------        
           readin_size = kpres*ktemp*nwnos
           IF(ALLOCATED( crossec_readin ) ) THEN
            deallocate(crossec_readin)
           ENDIF
           allocate(crossec_readin(readin_size))

C-      Readin kappa to flat array        
           open (action='read', file=file_data, newunit=fu, 
     &        FORM='unformatted')
           read (unit=fu) crossec_readin
           close (fu)
   
C-      Unfold to output array
   
           readin_index=1
           DO freq_i=1,nwnos
             DO p_i=1,kpres
               DO t_i=1,ktemp
                 readin_index = freq_i + nwnos * (p_i-1) +
     &                             nwnos * kpres * (t_i-1)
                 crossec(freq_i,p_i,t_i) =
     &                crossec_readin(readin_index)
               END DO
             END DO
           END DO
           
C       ---------------------
C       OPACITY WAVENUMBER IO 
C       ---------------------
C-      Readin temperature to flat array        
           open (action='read', file=file_wnos, newunit=fu, 
     &        FORM='unformatted')
           read (unit=fu) wnmol
           close (fu)        

           do freq_i =2, nwnos
             if ((wnmol(freq_i)-wnmol(freq_i-1)) .lt. 0) then
             write(*,*) 'ERROR: wavenumbers need to '
             write(*,*) '  be in in increasing sorted order'
             STOP
             end if
           end do  
         end subroutine read_opac_data
   
   
         subroutine read_opac_namelist(file_path,
     &        NOPI, NOTI,          
     &        kpres, ktemp, molid, nwnos, pmol, 
     &        tmol)
           implicit none
           
           !! Reads Namelist from given file.
           character(len=*),  intent(in)    :: file_path
           integer,  intent(in)    :: NOPI, NOTI
   
           integer,  intent(out)                  :: kpres, ktemp, nwnos
           character(len=4),  intent(out)         :: molid
           double precision, intent(out) :: pmol(100*NOPI)
           double precision, intent(out) :: tmol(100*NOTI)
   
           ! local
           integer                          :: fu, rc
           integer                          :: max_ti, max_pi
   
   
           ! Namelist definition.
           namelist / inputosmol / kpres, ktemp, molid, 
     &           nwnos, pmol, tmol

           pmol(:) = 0.0 
           tmol(:) = 0.0 
           
           ! Open and read Namelist file.
           open (action='read', file=file_path, newunit=fu)
           read (nml=inputosmol, unit=fu)
           close (fu)
   
           max_ti = minloc(tmol, 1)-1
           max_pi = minloc(pmol, 1)-1
           
           if (max_pi .gt. NOPI) then
             write(*,*) 'ERROR: increase NOPI'
             write(*,*) 'in parameter.inc'
             STOP
           endif
           if (max_ti .gt. NOTI) then
             write(*,*) 'ERROR: increase NOTI'
             write(*,*) 'in parameter.inc'
             STOP
           endif
           if (max_pi .ne. kpres) then
             write(*,*) 'ERROR: input pressure is wrong'
             STOP
           endif
           if (max_ti .ne. ktemp) then
             write(*,*) 'ERROR: input temperature is wrong'
             STOP
           endif
   
   
         end subroutine read_opac_namelist 
   
         subroutine search_intp_ind(binbord,binbordlen,val,intpint)
           ! ADAPTED FROM PETITRADTRANS
   
           implicit none
         
           INTEGER            :: binbordlen, intpint
           DOUBLE PRECISION   :: binbord(binbordlen),val
           INTEGER            :: i_arr
           INTEGER            :: pivot, k0, km
         
           ! carry out a binary search for the interpolation bin borders
           
         
            if (val >= binbord(binbordlen)) then
               intpint = binbordlen - 1
            else if (val <= binbord(1)) then
               intpint = 1
            else
   
               k0 = 1
               km = binbordlen
               pivot = (km+k0)/2
   
               do while(km-k0>1)
   
                  if (val >= binbord(pivot)) then
                     k0 = pivot
                     pivot = (km+k0)/2
                  else
                     km = pivot
                     pivot = (km+k0)/2
                  end if
   
               end do
   
               intpint = k0
   
            end if
         
         end subroutine search_intp_ind
   
   
         subroutine interpol_kappa(press, temp, pgrid, tgrid, 
     &                  opa_in, opa_out,
     &                  struc_len, freq_len, kpres, ktemp)
            
            ! ADAPTED FROM PETITRADTRANS (interpol_opa_ck)
            implicit none
            ! I/O
            integer, intent(in) :: struc_len
            integer, intent(in) :: freq_len
            integer, intent(in) :: kpres
            integer, intent(in) :: ktemp
   
            double precision, intent(in) :: press(struc_len)
            double precision, intent(in) :: temp(struc_len)
            double precision, intent(in) :: pgrid(kpres)
            double precision, intent(in) :: tgrid(ktemp)
   
            double precision,intent(in)::opa_in(freq_len, kpres, ktemp)
            double precision,intent(out)::opa_out(freq_len, struc_len)
   
            ! local
            double precision :: temp_min, temp_max, PorT
            integer          :: buffer_scalar
            double precision :: slopes(freq_len)
            double precision :: buffer1(freq_len), buffer2(freq_len)
            double precision :: buffer_Ts(freq_len), buffer_Tl(freq_len)
            integer :: i_str, ts, tl, ps, pl
            temp_min = MINVAL(tgrid)
            temp_max = MAXVAL(tgrid)

            do i_str = 1, struc_len
               call search_intp_ind(tgrid,ktemp,temp(i_str),
     &               buffer_scalar)
               ts = buffer_scalar 
               tl = buffer_scalar + 1
   
               call search_intp_ind(pgrid,kpres,press(i_str),
     &               buffer_scalar)     
               ps = buffer_scalar
               pl = buffer_scalar + 1

               ! Interpolate...
         
               !**********************************************************
               ! Interpolation to correct pressure at smaller temperatures
               !**********************************************************
               ! kappas
               if (kpres > 1) then         
                  ! kappa at smaller T and smaller P
                  buffer1 = opa_in(:, ps, ts)
                  ! kappa at smaller T and larger P
                  buffer2 = opa_in(:, pl, ts)
            
                  PorT = log(press(i_str))-log(pgrid(ps))
            
                  slopes = (buffer2-buffer1)/
     &                  (log(pgrid(pl))-log(pgrid(ps)))
            
                  if (press(i_str) >= pgrid(pl)) then
                        buffer_Ts = buffer2
                  else if (press(i_str) <= pgrid(ps)) then
                        buffer_Ts = buffer1
                  else
                        buffer_Ts = buffer1 + slopes*PorT
                  end if
               else
                  ! There is only one pressure value -> 1D interpolation in T instead
                  buffer_Ts = opa_in(:, 1, ts)
               endif
         
               !*********************************************************
               ! Interpolation to correct pressure at larger temperatures
               !*********************************************************
               if (kpres > 1) then
                  ! kappa at larger T and smaller P
                  buffer1 = opa_in(:, ps, tl)
                  ! kappa at larger T and larger P
                  buffer2 = opa_in(:, pl, tl)
            
                  PorT = log(press(i_str))-log(pgrid(ps))
            
                  ! slopes to correct to correct pressure are larger T
                  slopes = (buffer2-buffer1)/
     &                  (log(pgrid(pl))-log(pgrid(ps)))
            
                  ! kappa at larger temperature and correct pressure
                  if (press(i_str) >= pgrid(pl)) then
                        buffer_Tl = buffer2
                  else if (press(i_str) <= pgrid(ps)) then
                        buffer_Tl = buffer1
                  else
                        buffer_Tl = buffer1 + slopes*PorT
                  end if
               else
                  ! There is only one pressure value -> 1D interpolation in T instead
                  buffer_Tl = opa_in(:, 1, tl)
               endif
         
               !***********************************************************
               ! Interpolation to correct pressure and correct temperatures
               !***********************************************************
         
               PorT = temp(i_str)-tgrid(ts)
         
               slopes = (buffer_Tl-buffer_Ts)/(tgrid(tl)-tgrid(ts))
         
               if (temp(i_str) >= temp_max) then
                  opa_out(:,i_str) = buffer_Tl
               else if (temp(i_str) <= temp_min) then
                  opa_out(:,i_str) = buffer_Ts
               else
                  opa_out(:,i_str) = buffer_Ts + slopes*PorT        
               end if
         
            end do
       
         end subroutine interpol_kappa   
      END MODULE OPAMOD

      subroutine opac_wrapper_read(file_b, spec_i, wn_grid)
          use OPAMOD
          implicit real*8 (a-h,o-z)
          include 'parameter.inc'  
 
          integer, intent(in) :: spec_i
          double precision, intent(in) :: wn_grid(NWL)
          character(len=*), intent(in) :: file_b
      
          call read_opac(file_b, spec_i, wn_grid)
         
      end subroutine opac_wrapper_read
      subroutine opac_wrapper_interp(press_cgs, temp,
     &         spec_i, opa, struc_len)
              use OPAMOD
              implicit real*8 (a-h,o-z)
              include 'parameter.inc'  
     
              integer, intent(in) :: spec_i, struc_len
              double precision, intent(in) :: press_cgs(struc_len)
              double precision, intent(in) :: temp(struc_len)

              double precision, intent(out) :: opa(NWL,struc_len)
              ! local:
              double precision :: press_bar(struc_len)
            
              ! NOTE: pressure in opacity array is currently in bar !
              press_bar = press_cgs*1e-6
              call get_opac(press_bar, temp, spec_i, opa, struc_len)
     
      end subroutine opac_wrapper_interp    
      ! end ADS
C     End of Aaron's routines (module OPAMOD and a few more subroutines) to read 
C     two dimensional OS files of molecular cross sections (absorption coefficients)
C     and interpolate linearly to the temperatures and gas pressures
C     of a full MSG model structure. The interpolation in done for the
C     nearest input wavenumbers from the OS files to those in the input
C     wavenumber set for solution of the radiative transfer (because OS
C     is a statistical pick). 
      

      SUBROUTINE INITAB (IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3,INWNFIL*60

      NAMELIST /INPUTOS/ WNB,WNSTEP,WNEND
     *         ,INTVOS,nchrom,masabs
     *         ,losresl,osresl,wnos_first,wnos_last,kos_step
     *         ,LISTWN,INWNFIL,NEWC3
     *         ,newosatom,newosatomlist
      
      character(len=200) :: filebdir(maxosmol)     
      character atnames*2, molnames*8, mol_file*20      
      logical ggchem_mol(maxosmol), ggchem_index_read
      integer ggchem_index(maxosmol), molno
      character(len=5) molnames_new(maxosmol)
      common /molupdate/ molnames_new, 
     * ggchem_index,
     * molno, ggchem_mol, ggchem_index_read     
      COMMON/COSWR/osresl,losresl,listwn
      COMMON/CNEWC3 /NEWC3
      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON/CA1/DUMA(120),IDUM(240),NEXTT,NUTZT
      COMMON/CFIL/IRESET(10),ISLASK,IREAT
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON /CROS/WROS(20)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/COUTR/NTO,NTPO(10)
      COMMON/CVAAGL/XLA(500),W(500),NLB
      COMMON/CLINE4/ILINE
      common /cmasabs/ masabs(3)
      common /cinos/resl,wnos_first,wnos_last,kstep
C
      resl = OSRESL
      kstep = kos_step
C
C LOGICAL UNITS
      ISLASK=11
      IREAT=9
      ISAVE=IWRIT
      DO 1 I=1,10
 1    IRESET(I)=10
C
C ZEROSET WEIGHTS
      DO 4 I=1,500
 4    W(I)=0.
      DO 5 I=1,20
 5    WROS(I)=0.
C
C WAVELENGTH SETS
C ROSSELAND
      CALL VAAGL(NLBRO,XL,WROS)
      NL(1)=NLBRO
      CALL VAAGL(NLB,XLA,W)
      LMAX=15
      IFIRST=2
      ILAST=10
      CALL SETDIS(NLB,XLA,LMAX,IFIRST,ILAST)
C 8000. $NGSTR@M STANDARD
      NL(1)=NL(1)+1
      NL1=NL(1)
      XL(NL1,1)=8000.
      IF(IOUTS.LT.0) IWRIT=4
      CALL INABS(IOUTS)
      IF(IOUTS.GE.-1)IWRIT=ISAVE
      NL(1)=NL(1)-1
C
C CONTROL PARAMETTERS
      NTO=0
      NEXTT=0
      NUTZT=0
C
C  Calculate the wave numbers, WN, for OS opacity and the total number, 
C  NWTOT, of OS frequency points used. This wavenumber scale should be 
C  identical to the one used when the OS-tables were created.
C
      READ(5,INPUTOS)
C     ADS: Note: mol_names.dat can be molecules or atoms
      mol_file = "./data/mol_names.dat"


      open(unit=7397, file=mol_file)
         read(7397, '(i2)') molno
         do n=1,molno
            read(7397, '(A,A)') molnames_new(n), filebdir(n)
            molname(n) = trim(molnames_new(n))
         end do
      close(7397)
      print*, "molecules/atoms included are: "
      print*, molname(1:molno)
      nosmol = molno
      write(7,396) nosmol
396   format("The model used",i4," molecules/atoms in the opacity"
     & ," computation. They were:")
C & /" (from mol_names.dat read in initab) really are: ")
      write(7,397) molname(1:molno)
397   format(18(x,a4))
C     Reminder: We need to match the ggchem_indices to our read in mols later!
      ggchem_index_read = .FALSE.

      IF(molno.GT.MAXOSMOL) STOP ' Increase dimension for MAXOSMOL'
C

C  Calculate (or read, if listwn=1) the wave numbers, WN, for OS opacity
C  If losresl = 1, the os wavenumbers are computed based on a specified
C  and the total number, NWTOT, of OS frequency points for the OS-tables.
C  spectral resolution, osresl. Else it is computed in fixed steps inside a
C  number of prespecified intervals.
C

      
      !IF (LISTWN .eq. 1) go to 225   ! read an existing OS - wn list


C  compute an OS - wn list:

      !if (losresl.eq.0) go to 228     ! compute list with fixed steps

C  use OS list with fixed resolution, osresl, through spectrum:
      wnos(1) = wnos_first
      step = 1.d0 + 1.d0/osresl
      do 240 k = 1,nwl-1
      wnkj = wnos(k)
      do 242 kj = 1,kos_step
242   wnkj = wnkj * step
      nwtot = k + 1
      wnos(nwtot) = wnkj
      if (wnos(nwtot) .gt. wnos_last) then
            go to 241
      end if
240   continue
      write(7,247)nwtot,wnos(1),wnos(wntot)
247   format('nwtot,wnos(1),wnos(wntot):'i6,f6.1,f8.1)
C we come here only if dimension for the OS is too small:
      wnos1 = wnos_first
      step = 1.d0 + 1.d0/osresl
      do 740 k = 1,100000
      wnkj = wnos1
      do 742 kj = 1,kos_step
742   wnkj = wnkj * step
      nwtot = k + 1
      wnos1 = wnkj
      if (wnos1 .gt. wnos_last) go to 741
740   continue
741   continue
      print*,' given,needed dimensions -nwl,nwtot='
     &         ,nwl,nwtot
      stop ' error: increase dimension nwl in parameter.inc for wnos '
241   continue

      if (molno .GT. NOSPEC) then
      print*,' total molno=', molno
      stop ' error: increase dimension NOSPEC to molno in parameter.inc'
      endif

      do nm=1,molno
         call opac_wrapper_read(trim(filebdir(nm)), nm, wnos) 
      end do
      !stop !COMMENT THIS STOP OUT WHEN LWRITEOS IS NOT 1 ANYMORE 
      !(WHICH IS THE CASE WHEN WE WANT TO RUN THE FULL SIMULATION WITH THE NEW OPACITY MODULE, OTHERWISE THIS MODULE IS DESIGNED TO JUST RUN ONCE TO READ IN THE OPACITY DATA AND WRITE IT OUT IN A PLOTABLE FORMAT)
      !write(7,*) "OS done with ADS's routine."     
      oskres = osresl/dfloat(kos_step)
      write(7,245) nwtot,wnos_first,wnos_last,1.e4/wnos_last,
     & 1.e4/wnos_first,kos_step,oskres
245   format('The OS cross sections (cm2/molecule) were read into',i6,
     &' wavenumbers between'/,f6.1,' and',f8.1,' cm-1 (=',f6.3,' to',
     & f5.0,' mu). We used only each',i3,'th in the radiative transfer',
     &/' of this model, resulting in a (statistical spectral) ',
     & 'OS-resolution of',f7.1)

      DO 103 I=1,NWTOT
      L=NWTOT-I+1
103   WLOS(L)=1.D8/WNOS(I)
      DO 104 L=2,NWTOT-1
      LP=L+1
      LM=L-1
104   WLSTEP(L) = ( WLOS(LP)-WLOS(LM) ) / 2.
      WLSTEP(1) = WLOS(2)-WLOS(1)
      WLSTEP(NWTOT) = WLOS(NWTOT)-WLOS(NWTOT-1)

      IF (NWTOT.GT.NWL)  STOP ' DIMENSION NWL TOO SMALL'
      
C
C The wavelengths for calculation of the line-opacity must be inside
C the vavelengths where the continuum-opacity is calculated.
C The first and last line-opacity wavelength are WLOS(1) and 
C WLOS(NWTOT), respectively.
C The first and last continuum-opacity wavelength are XL(1,2) and 
C XL(NL(NLB),NLB), respectively.
C
      IF (WLOS(1) .LE. XL(1,2)) THEN
           PRINT*, ' Error-message from INITAB.FOR:'
           PRINT*, ' WLOS(1) = ',WLOS(1)
           PRINT*, ' XL(1,2) = ',XL(1,2)
           STOP ' WLOS(1) < first cont. point'
      END IF
      IF (WLOS(NWTOT) .GT. XL(NL(NSET),NSET)) THEN
           PRINT*, ' Error-message from INITAB.FOR:'
           PRINT*, 'NLB = ',NLB
           PRINT*, 'NSET = ',NSET
           PRINT*, 'NL(NSET) = ',NL(NSET)
           PRINT*, 'XL(NL(NSET),NSET) = ',XL(NL(NSET),NSET)
           STOP ' WLOS(last) > last cont. point'
      END IF
C
      RETURN
      END

      SUBROUTINE init_thermo()
        implicit real*8 (a-h,o-z)
        include 'parameter.inc'
        character(len=200) line
        character(len=56) comments
        character(len=6) id
        character(len=10) spname
        character(len=2) el(5)
        integer type,ggi,ggindex,N_matched,i,j,nT
        logical species_found, atactive, molactive
        logical is_atom
        dimension st(5),ai(7),ci(7)
        
        character atnames*2, molnames*8
        common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
! Output
        common / thermodata / xmoltha(543,3,7),xmolthTlim(543,4),
     >               xattha(22,3,7),xatthTlim(22,4), thexp(7),
     >               xtemplimlow,xtemplimup,molactive(543),
     >               atactive(22)


C       From nasa9.dat, temperature limits, sets default values
        xtemplimlow = 200.0d0
        xtemplimup = 20000.0d0
        xmolthTlim(:,1) = xtemplimlow
        xmolthTlim(:,2) = 1000.0d0
        xmolthTlim(:,3) = 6000.0d0
        xmolthTlim(:,4) = xtemplimup
        xatthTlim(:,1) = xtemplimlow
        xatthTlim(:,2) = 1000.0d0
        xatthTlim(:,3) = 6000.0d0
        xatthTlim(:,4) = xtemplimup

C      Set temperature exponents (for nasa-9 format):
        thexp(1) = -2.0d0
        thexp(2) = -1.0d0
        thexp(3) = 0.0d0
        thexp(4) = 1.0d0
        thexp(5) = 2.0d0
        thexp(6) = 3.0d0
        thexp(7) = 4.0d0

        xmoltha(:,:,:) = 0.0d0
        xattha(:,:,:) = 0.0d0
        molactive(:)=.FALSE.
        atactive(:)=.FALSE.

C       read comments
        open(8897, file='data/nasa9.dat',status='old')
        do i=1,58
            read(8897,*) line
        end do

c       Iteration over the content blocks
        N_matched = 0
        rloop: do 
          read(8897, '(A10,A56)',end=3923) spname, comments
          read(8897, '(I2,A8,5(A2,F6.2),I2,F13.5,F13.5)') nT, id, 
     &                               el(1), st(1),el(2),st(2),el(3),
     &                               st(3),el(4),st(4), el(5),st(5), 
     &                               type, xm, heatf
C         GGchem matching here 
          species_found=.FALSE.  
          atloop: do ggi=1,22
             if (upper(trim(atnames(ggi))).EQ.upper(trim(spname)))
     &          then
                ggindex = ggi
                species_found=.TRUE.
                is_atom=.TRUE.
                atactive(ggi)=.TRUE.
                exit atloop
             end if
          end do atloop
          if (.NOT. species_found) then
            molloop: do ggi=1,543
               if (upper(trim(molnames(ggi))).EQ.upper(trim(spname))) 
     &            then
                  ggindex = ggi
                  species_found=.TRUE.
                  is_atom=.FALSE.
                  molactive(ggi)=.TRUE.
                  exit molloop
               end if
            end do molloop   
          endif    
C         
C          ADS: debugging purposes          
C          if (species_found) then
C             N_matched = N_matched+1
C          else 
C             print*, trim(spname), 'not found'
C          end if
             

          if (type.NE.0) goto 3923 ! Done with the readin   

          if(nT.GT.3) then 
            print*, 'error readin thermo data: nT>3'
            STOP
          end if

          do i=1,nT
            read(8897, '(2(F11.3),I1,8(F5.1),F15.3)') tlow,tup,ncoeff,
     &                            ci(1),ci(2),ci(3),ci(4),ci(5),
     &                            ci(6),ci(7),dump, dH0
            read(8897, '(5(D16.8))') ai(1),ai(2),ai(3),ai(4),ai(5)
            read(8897, '(5(D16.8))') ai(6),ai(7),dump,b1,b2

            do j=1,7
               if(abs(ci(j)-thexp(j)).GT. 0.1d0) then 
                  print*, 'error readin thermo data, exp mismatch'
                  STOP
               endif
            end do

            if (species_found) then
                if (is_atom) then
                   xattha(ggindex,i,:)=ai(:)
                   if(i.eq.1) xatthTlim(ggindex,1)=tlow
                   xatthTlim(ggindex,i+1)=tup 
                else
                   xmoltha(ggindex,i,:)=ai(:)
                   if(i.eq.1) xmolthTlim(ggindex,1)=tlow
                   xmolthTlim(ggindex,i+1)=tup 
                end if
            end if
          end do
        end do rloop
3923    continue

        close(8897)
C        print*, N_matched

        contains 
        function upper(strIn) result(strOut)
        ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
        ! Original author: Clive Page

        implicit none
  
        character(len=*), intent(in) :: strIn
        character(len=len(strIn)):: strOut
        integer :: i,j
  
        do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j>= iachar("a") .and. j<=iachar("z") ) then
                    strOut(i:i) = achar(iachar(strIn(i:i))-32)
              else
                    strOut(i:i) = strIn(i:i)
              end if
        end do
      end function upper

      end subroutine init_thermo

C
      SUBROUTINE INITJN(IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
      COMMON /UTPUT/IREAD,IWRIT
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac, icloud_conv
C
C INITJN INITIATES THE JON BLOCK FROM LUN 9
C
      ISAVE=IREAD
      IREAD=9              !unit=9 is jonabs.dat
*      CALL DUMIN
      CALL INJON(IOUTS)
      IREAD=ISAVE
      RETURN
      END
C
      SUBROUTINE INJON(IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
C
C        THIS ROUTINE READS DATA NECESSARY FOR THE COMPUTATION OF IONIZATION
C        EQUILIBRIA ETC. (IN SUBR. JON).
C        1. NEL= THE NUMBER OF CHEMICAL ELEMENTS CONSIDERED.
C           A=   THE RATIO OF THE NUMBER OF HYDROGEN NUCLEI TO THE NUMBER OF
C                NUCLEI OF METALLIC ELEMENTS.
C           NMET=THE NUMBER OF METALLIC ELEMENTS IN THE LIST OF CHEMICAL
C                ELEMENTS CONSIDERED. THE LAST NMET ELEMENTS OF THE LIST ARE
C                CONSIDERED TO BE METALLIC, FOR THE CALCULATION OF THE
C                QUANTITY A (DEFINED ABOVE).
C        2. IEL  IS THE ARRAY WHICH WILL CONTAIN THE SYMBOLS FOR THE CHEMICAL
C                ELEMENTS CONSIDERED.
C           ABUND IS THE ARRAY WHICH WILL CONTAIN THE PREVAILING ABUNDANCES
C                THE CHEMICAL ELEMENTS CONSIDERED AT INPUT. THESE ABUNDANCES
C                ARE EXPRESSED ON A LOGARITHMIC SCALE (BASE 10) AND NEED NOT BE
C                NORMALIZED. THE ABUNDANCES ARE MODIFIED IN THIS SUBROUTINE
C                SO THAT THE RIGHT VALUE OF A (DEFINED ABOVE) IS OBTAINED.
C        3. AI   IS THE ARRAY WHICH WILL CONTAIN THE ATOMIC WEIGHTS OF THE
C                ELEMENTS CONSIDERED.
C        4. DATA FOR THE COMPUTATION OF THE PARTITION FUNCTIONS IS READ NEXT
C           NJ(I)= THE NUMBER OF STAGES OF IONIZATION CONSIDERED FOR ELEMENT I.
C         FOR EACH STAGE OF IONIZATION JA THE FOLLOWING QUANTITIES ARE READ
C           G0(JA)=THE STATISTICAL WEIGHT OF THE GROUND LEVEL,
C           NK(JA)=THE NUMBER OF ELECTRON CONFIGURATIONS CONSIDERED.
C         FOR EACH ELECTRON CONFIGURATION JB THE FOLLOWING QUANTITIES ARE READ
C           XION(JB)=THE IONIZATION ENERGY IN ELECTRON VOLTS,
C           G2(JB)=THE STATISTICAL WEIGHT (2L+1)*(2J+1)
C           XL(JB)=THE LOWEST QUANTUM NUMBER OF THE ASYMPTOTIC (HYDROGENIC) PART
C                OF THE PARTITION FUNCTION,
C           NL(JB)=THE NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSION FOR THE
C                'MIDDLE PART' OF THE PARTITION FUNCTION ('QPRIME').
C           ALFA IS AN ARRAY WHICH WILL CONTAIN THE 'STATISTICAL WEIGHTS' OF
C                THE (APPROXIMATE) EXPRESSIONS FOR THE 'MIDDLE PARTS' OF THE
C                PARTITION FUNCTIONS.
C           GAMMA IS AN ARRAY CONTAINING THE CORRESPONDING 'EXCITATION
C                POTENTIALS' (EXPRESSED IN ELECTRON VOLTS).
C         FOR THE METHOD USED SEE TRAVING ET AL., ABH. HAMB. VIII, I (1966).
C        5. ELEMENTS AND STAGES OF IONIZATION THAT SHOULD BE DISREGARDED ARE
C           INDICATED BY IELEM(I)=0 FOR ELEMENT I AND BY ION(I,J)=0 FOR
C           IONIZATION STAGE J. OTHERWISE INDICATORS SHOULD BE = 1.
C        6. NQFIX IS THE NUMBER OF PARTITION FUNCTIONS THAT SHOULD BE CONSTANT.
C                THE VALUES ARE READ INTO THE VECTOR PARCO AND AN INDICATION IS
C                MADE IN IQFIX.  IQFIX(I,J)=0 MEANS THAT THE PARTITION FUNCTION
C                FOR ELEMENT I, STAGE OF IONIZATION J, IS CONSIDERED TO BE
C                CONSTANT.
C           NQTEMP IS THE NUMBER OF PARTITION FUNCTIONS  THAT SHOULD BE
C                PRESSURE-INDEPENDENT AND INTERPOLATED IN T. VALUES OF FOUR
C                TEMPERATURES (TPARF, THE SAME FOR ALL ELEMENTS) AND
C                CORRESPONDING PARTITION FUNCTIONS (PARF) ARE READ. IQFIX(I,J)=1
C                MEANS THAT A PRESSURE-INDEPENDENT PARTITION FUNCTION FOR INTER-
C                POLATION IN T IS GIVEN.
C        7. IFISH IS A PARAMETER FOR THE CHOICE OF THE ASYMPTOTIC PARTITION 
C                FUNCTION. IFISH=0 MEANS THAT THE ASYMPTOTIC PART WILL BE EVALU-
C                ATED FOLLOWING BASCHEK ET AL., ABH. HAMB. VIII,26 (1966). IFISH
C                =1 MEANS THAT IT WILL BE EVALUATED FOLLOWING FISCHEL AN SPARKS
C                ASTROPHYS. J. 164, 356 (1971).
C        8. TMOLIM IS THE HIGHER TEMPERATURE LIMIT BEYOND WHICH MOLECULES WILL
C                NOT BE CONSIDERED
C
C        MOREOVER SOME INITIATING WORK IS DONE FOR SUBR. JON. UNLOGARITHMIC
C        ABUNDANCES ARE NORMALIZED ON HYDROGEN, XMY AND SUMH (DEFINED BELOW)
C        ARE COMPUTED AND SOME FURTHER QUANTITIES ARE EVALUATED AT THE END.
C        A DETAILED PRINTOUT IS GIVEN IF IOUTS IS EQUAL TO ONE. AFTER INJON
C        HAS BEEN CALLED ONCE, A NEW DETAILED PRINTOUT IS OBTAINED IF
C        INJON IS CALLED WITH IOUTS GREATER THAN ONE.
C
C        DIMENSIONS NECESSARY
C        ABUND(NEL),AI(NEL),ALFA(LMAX),ANJON(NEL,MAX(NJ)),FL2(5),F1Q(3),
C        GAMMA(LMAX),G0(JMAX),G2(KMAX),H(5),IEL(NEL),IELEM(NEL),
C        ION(NEL,MAX(NJ)),IQFIX(NEL,MAX(NJ)),JAMEM(NEL,MAX(NJ)),JBBEG(JMAX)
C        JCBEG(JMAX),NJ(NEL),NK(JMAX),NL(KMAX),PARCO(JMAX),PARF(4*JMAX),
C        PARPP(4),PARPT(4),PARQ(4*JMAX),PART(NEL,MAX(NJ)),SHXIJ(5),TPARF
C        XION(KMAX),XIONG(NEL,MAX(NJ)),XL(KMAX)
C        THE DIMENSIONS ARE LOWER LIMITS
C        JMAX IS THE TOTAL NUMBER OF STAGES OF IONIZATION, INCLUDING NEU
C             ATOMS.
C        KMAX IS THE TOTAL NUMBER OF ELECTRON CONFIGURATIONS.
C        LMAX IS THE TOTAL NUMBER OF TERMS IN THE (APPROXIMATE) EXPRESSI
C             FOR THE MIDDLE PART OF THE PARTITION FUNCTIONS ('QPRIME'),
C             ACCORDING TO TRAVING ET AL., CITED ABOVE.
C        NEL  IS THE NUMBER OF CHEMICAL ELEMENTS.
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUT
C             STAGE, FOR ELEMENT I.
C
      DIMENSION AI(16),F1Q(3),F2Q(2),PARF(180),PARPP(4),PARPT(4)
      DIMENSION JAMEM(16,5)
      character*4 abcname
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,NJ(16),IEL(16),
     *SUMM(ndp),NEL
      COMMON/CI9/AI
      COMMON/CI3/ALFA(300),GAMMA(300),G0(45),G2(80),XION(80),XL(80),
     *JBBEG(45),JCBEG(45),NK(45),NL(80),IFISH
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CI5/abmarcs(18,ndp),ANJON(18,5),H(5),PART(18,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI6/TP,IQFIX(16,5),NQTEMP
      COMMON/UTPUT/IREAD,IWRIT
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500)
      INTEGER MOLH, JUMP
      character sunz*1,head_elabund*100
      common/cabinit/abinit(natms),kelem(natms),nelem
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     *  vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     * ntau,iter
      dimension mx_elm(18),abundatms_inp(natms),sum(ndp),fakt(ndp)
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      namelist /abundances/sunz,zscale,abundatms_inp
      data mx_elm /1, 2,6,7,8,10,11,12,13,14,16,19,20,23,25,27,21,17/
C                  H He C N O Ne Na Mg Al Si  S K Ca Cr Fe Ni Ti Cl
      character*3 aifix
      common /pefix/xionfix,aifix
      common/cabnames/abcname(natms)
      common /cisph/isph
      COMMON /CSTYR/MIHAL,NOCONV 
      COMMON /CG/GRAV,KONSG 
      common /CTEFF/TEFF,FLUX
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      character*4 osnames(54)
C The 54 molecules we have absorption coefficients for (March 2026) are:
      data osnames /
     & 'C2  ','CAH ','CH  ','CN  ','CO  ','CO2 ','FEH ','H2O ','HCN ',
     & 'MGH ','NH  ','OH  ','SIH ','SIO ','TIH ','TIO ','CRH ','NO  ',
     & 'LIH ','VO  ','ZRO ','H2  ','C2H2','C3  ','ALCL','ALF ','ALH ',
     & 'ALO ','BEH ','CAF ','CH3F','CH4 ','CP  ','CS  ','H2CO','HCL ',
     & 'HNO3','KCL ','KF  ','LICL','LIF ','MGF ','NACL','NAF ','NAH ',
     & 'NH3 ','SN  ','PH3 ','PN  ','PO  ','PS  ','HS  ','SIS ','SO2 '/

C
C
        IF(IOUTS.GT.1) GO TO 25
C
C        READING OF THE ABUNDANCES AND THEI ASSOCIATED QUANTITIES
C        **** 1 ****

       NEL = 16
       A = 0.
       NMET = 0
C
C        **** 2 ****
C

      READ(IREAD,110)(IEL(I),I=1,NEL)    !the 16 element names H - Ni

C The abundances are now read from elabund.dat i gem_init, and can be
C adjusted collectively or individually in input ($abundances). The 
C 17 abundances usually read here from input file are brought into injon
C and other routines from gem_init by common CI5.
      open(unit=2,file='data/elabund.dat',status='old',readonly)
      read(2,299) head_elabund
299   format(A100)
      natmsact = 0
      do i=1,200
        read(2,*,iostat=io) kelem(i), abinit(i), abcname(i)
        if(io .ne. 0) exit
        natmsact = natmsact + 1
      end do
      nelem = i-1

      WRITE(7,283)
      if (isph==0) write(7,181)teff,log10(grav)
      if (isph==1) write(7,182)teff,log10(grav)
181   format('   Plane-parallel Cool-MARCS (=MSG) model for Teff =',
     &   f6.0,' K, log(g) =',f4.1)
182   format('   Spherical Cool-MARCS (=MSG) model for Teff =',
     &   f6.0,' log(g) =',f4.1)
      WRITE(7,284)
  283 FORMAT(1X,81('*'))
  284 FORMAT(1X,81('*')/)

      write(7,270) jtau,log10(tau(1)),log10(tau(jtau))
270   format('The model is computed in',i3,
     &  ' layers from log_10(tau_ross)=',f6.2,' to',f6.2)

      read(5,abundances)
      !BCE (10.02.23 - introducing metallicity changes)
      !if Z= solar, then sunz=y and nothing changes
      !else add to the abundances of all elements but H and He log(zscale).
      print*, "Metallicity is ", zscale, " time(s) solar."
      !REWRITE HEAD OF ELABUND AND SCREEN OUTPUT WITH GENERALIZED INPUT FORMAT AND OUTPUT FORMAT

      write(7,*) 'This model of the general type ', trim(head_elabund)
      do i=1, nelem
      if (i>2) then
            abinit(i) = abinit(i) + log10(zscale)
      end if
      end do
      if (zscale.ne.1.0) write(7,298) nelem-2,zscale
298   format("but thereafter scaled the",I3,
     &" elements that were not H or He with a factor",1Pe9.2)
      write(7,297)
297   format("meaning that at a scale where epsilon(H) = log(H)=12,",
     &     " epsilon(elm) included:") 
      write(7,296)abinit(1),abinit(2),(abinit(l),l=6,8),
     &  (abinit(l),l=11,14),(abinit(l),l=19,22),abinit(25),abinit(35)
296   format(2f6.2,13f5.2)
      write(7,295) abcname(1),abcname(2),(abcname(l),l=6,8),
     &(abcname(l),l=11,14),(abcname(l),l=19,22),abcname(25),abcname(35)
295   format(2(2x,a4),x,13(x,a4))
      do i=1,18
          abmarcs(i,1:ntau) = abinit(mx_elm(i))
      end do
      
C
C        **** 3 ****
C
      READ(IREAD,102)(AI(I),I=1,NEL)      !the atomic weight 1.0,4.0,...,58.7

      NU=NEL
      SUM(1:ntau)=0.
      SUMM(1:ntau)=0.
      FAKT(1:ntau)=1.
C
C        THE ABUNDANCES ARE CONVERTED FROM A LOGARITHMIC SCALE TO A DIRECT
C        SCALE, AND ARE THEN NORMALIZED ON HYDROGEN. XMY=GRAMS OF STELLAR MATTER
C        /GRAMS OF HYDROGEN. SUMH=NUMBER OF OTHER NUCLEI/NUMBER OF HYDROGEN
C        NUCLEI.
C        SUMM=NUMBER OF NUCLEI OTHER THAN H, C, N, O / NUMBER OF HYDROGEN
C
      if(idust .eq. 0) then
      if(nmet.le.0) go to 22      !->which seems to be always the case
      nu=nel-nmet+1
      do i=nu,nel
            abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau) 
            sum(1:ntau)=abmarcs(i,1:ntau)+sum(1:ntau)
      end do
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)
      abmarcs(18,1:ntau)=10.**abmarcs(18,1:ntau)
      fakt(1:ntau)=sum(1:ntau)*a/10.**abmarcs(1,1:ntau)
      nu=nu-1
   22 do i=1,nu
         abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau)*fakt(1:ntau)  !->fakt==1., so this is just sum abund
         sum(1:ntau)=sum(1:ntau)+abmarcs(i,1:ntau)
      end do
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)
      abmarcs(18,1:ntau)=10.**abmarcs(18,1:ntau)
      xmy(1:ntau)=0.
      aha=abmarcs(1,1)
      do3 i=1,nel
      abmarcs(i,1:ntau)=abmarcs(i,1:ntau)/aha
      summ(1:ntau)=summ(1:ntau)+abmarcs(i,1:ntau)
    3 xmy(1:ntau)=xmy(1:ntau)+abmarcs(i,1:ntau)*ai(i)   !AI(I)=atomic weight, so XMY=#AU/H_nuclei
      abmarcs(17,1:ntau)=abmarcs(17,1:ntau)/aha
      abmarcs(18,1:ntau)=abmarcs(18,1:ntau)/aha

      xmy(1:ntau)=xmy(1:ntau)/ai(1)            !AI(1) = 1.008, so on AU scale.
      sumh(1:ntau)=sum(1:ntau)/aha-1.
      summ(1:ntau)=summ(1:ntau)-abmarcs(1,1:ntau)-abmarcs(3,1:ntau)-
     &  abmarcs(4,1:ntau)-abmarcs(5,1:ntau)
      end if
C
C        **** 4 ****
C
C        READING OF DATA FOR THE PARTITION FUNCTIONS.
C        FOR THE SYMBOLS, SEE ABOVE.
C        28/10/21 as of this date this is only relevant to obtain
c        data of relevant constants for absorption coefficients 
c        computation. The partition functions are no longer used.

      READ(IREAD,103)(NJ(I),I=1,NEL)
      JA=1
      JB=1
      JC1=1
      DO11 I=1,NEL

      NJP=NJ(I)
      DO11 J=1,NJP
      JAMEM(I,J)=JA
      JBBEG(JA)=JB
      JCBEG(JA)=JC1
C        JBBEG AND JCBEG ARE INDICATORS USED BY FUNCTION QTRAV
C
      READ(IREAD,104)G0(JA),NK(JA)
      NKP=NK(JA)

      IQFIX(I,J)=2
C        IQFIX(I,J)=2 MEANS THAT A 'FULL' PARTITION FUNCTION SHOULD BE
C        COMPUTED. THIS MAY BE CHANGED UNDER **** 7 ****.
C
      JA=JA+1
      DO11 K=1,NKP
      READ(IREAD,105)XION(JB),G2(JB),XL(JB),NL(JB)

      IF(K.GT.1)GO TO 9
      XIONG(I,J)=XION(JB)
C        XIONG IS THE IONIZATION ENERGY IN ELECTRON VOLTS FOR THE GROUND STATE,
C        USED IN THE COMPUTATION OF IONIZATION EQUILIBRIA IN SUBROUTINE JON.
      if(i.eq.12 .and. j.eq.1) then
            iafix = iel(12)
C            write(7,1105) iel(12), xion(jb)
C            write(7,1106) iel(12), xiong(i,j)
            xionfix = xiong(i,j)
      end if
1105  format('ionization of element 12:',a3,' was read to',f8.3,' eV')
1106  format('ionization of element 12:',a3,' jon uses ',f8.3,' eV')
C
    9 CONTINUE
      JC2=NL(JB)+JC1-1
      JBM=JB
      JB=JB+1
      IF(NL(JBM).LE.0)GO TO 10
      READ(IREAD,106)(GAMMA(L),ALFA(L),L=JC1,JC2)
   10 JC1=JC2+1
   11 CONTINUE
C
C        **** 5 ****
C
C        READING OF THE INDICATORS OF THE ELEMENTS AND THE STAGES OF IONIZATION
C        TO BE DISREGARDED.
      DO12 I=1,NEL
      NJP=NJ(I)
      READ(IREAD,107)IELEM(I),(ION(I,J),J=1,NJP)

   12 CONTINUE
C
C        **** 6 ****
C
C        SPECIFICATION OF THOSE PARTITION FUNCTIONS GIVEN AS CONSTANTS.
C        INDICATION IN IQFIX.
      READ(IREAD,103)NQFIX
      IF(NQFIX.LE.0)GO TO 15
   13 DO14 I=1,NQFIX
      READ(IREAD,109)I1,J1,PARCOP
      JA=JAMEM(I1,J1)
      PARCO(JA)=PARCOP
   14 IQFIX(I1,J1)=0
   15 CONTINUE

C
C        SPECIFICATION OF THOSE PARTITION FUNCTIONS TO BE INTERPOLATED IN T.
C        INDICATION IN IQFIX.
      READ(IREAD,103)NQTEMP
      IF(NQTEMP.EQ.0)GO TO 20
      READ(IREAD,101)TPARF
      DO17 I=1,NQTEMP
      READ(IREAD,109)I1,J1,(PARPP(K),K=1,4)
      IQFIX(I1,J1)=1
C
C        PREPARATION FOR INTERPOLATION OF PARTITION FUNCTIONS IN T (CONCLUDED
C        IN SUBROUTINE JON).
      DO16 K=1,3
   16 F1Q(K)=(PARPP(K+1)-PARPP(K))/(TPARF(K+1)-TPARF(K))
      DO161 K=1,2
  161 F2Q(K)=(F1Q(K+1)-F1Q(K))/(TPARF(K+2)-TPARF(K))
      F3Q=(F2Q(2)-F2Q(1))/(TPARF(4)-TPARF(1))
      PARPT(1)=PARPP(1)
      PARPT(2)=F1Q(1)
      PARPT(3)=F2Q(1)
      PARPT(4)=F3Q
      JA=JAMEM(I1,J1)
      DO17 K=1,4
      JK=(JA-1)*4+K
      PARQ(JK)=PARPT(K)
   17 PARF(JK)=PARPP(K)
C        PARQ IS IN COMMON/CI1/ AND IS USED IN SUBROUTINE JON. PARF IS JUST
C        USED BELOW.
C
   20 CONTINUE
C
C        **** 7, 8 ****
C
C        THE PARAMETERS IFISH AND TMOLIM. INITIATING WORK FOR SUBROUTINE JON.
C        WHEN MOLH IS GREATER THAN ZERO THE MOLECULAR FORMATION WILL BE
C        IN SUBR. MOLEQ (ONLY H2 AND H2+), ELSE MORE COMPLETE MOLECULAR
C        FORMATION WILL BE EVALUATED IN SUBR. MOL.
C        --> not H too, because it is an atom !!
C
      READ(IREAD,100)IFISH

C TMOLIM, MOLH are given in JONABS.DAT
C here is set : IREAD=9 --> = JONABS.DAT
C we seemingly use (2023) tmolim=15000.   and molh=0

      READ(IREAD,4528) TMOLIM,MOLH
 4528 FORMAT(F10.0,I5)

      DO21 J=1,5
      FLJ=J
      FL2(J)=31.321*FLJ*FLJ
   21 SHXIJ(J)=SQRT(13.595*FLJ)
C
C        EEV=THE ELECTRON VOLT (EXPRESSED IN TERMS OF ERGS)
C        XMH=THE MASS OF THE HYDROGEN ATOM (EXPRESSED IN GRAMS)
C        XKBOL=BOLTZMANN'S CONSTANT (EXPRESSED IN ERGS PER KELVIN
      EEV=1.602095E-12
      XMH=1.67339E-24
      XKBOL=1.38053E-16
      if(idust.eq.0) enamn(1:ntau)=eev/(xmh*xmy(1:ntau))
      TP=0.
C        TP IS THE TEMPERATURE AT THE 'PRECEDING' CALL OF JON.
C
C        
      DO33 I=1,NEL
      NJP=NJ(I)
   33 CONTINUE
      JA=1
      DO32 I=1,NEL
      NJP=NJ(I)
      DO32 J=1,NJP
      JP=J-1
      JK1=(JA-1)*4+1
      JK2=(JA-1)*4+4
   32 JA=JA+1

      IF(IOUTS.LE.0)GO TO 40
   25 CONTINUE

   40 CONTINUE

C
      IWRIT=ISAVE
C
      RETURN
C
  100 FORMAT(I10,F10.4,I10)
  101 FORMAT(6F10.4)
  102 FORMAT(6F10.4)
  103 FORMAT(12I5)
  104 FORMAT(F5.0,I5)
  105 FORMAT(F6.3,F4.0,F5.1,I5)
  106 FORMAT(4(F10.3,F10.4))
  107 FORMAT(I10,5I5)
  108 FORMAT(2F10.4)
  109 FORMAT(2I5,4F10.4)
  110 FORMAT(16A3)
      END
C
      SUBROUTINE INP3(X,Y,XINT,YINT)
      implicit real*8 (a-h,o-z)
C
C        NEWTONINTERPOLATION, TREPUNKTS
C        OBS *** INGEN SPAERR MOT EXTRAPOLATION *****
C
      DIMENSION X(3),Y(3),F1(2)
      DO1 K=1,2
    1 F1(K)=(Y(K+1)-Y(K))/(X(K+1)-X(K))
      F2=(F1(2)-F1(1))/(X(3)-X(1))
      YINT=Y(1)+(XINT-X(1))*F1(1)+(XINT-X(1))*(XINT-X(2))*F2
      RETURN
      END
C
      SUBROUTINE JON(T,PE,IEPRO,PG,RO,E,IOUTR)
      implicit real*8 (a-h,o-z)
C
C
C        THIS ROUTINE COMPUTES IONIZATION EQUILIBRIA FOR A GIVEN TEMPERATURE 
C        (T, EXPRESSED IN KELVIN) AND A GIVEN ELECTRON PRESSURE (PE, IN
C        DYNES PER CM2). THE FRACTIONS OF IONIZATION ARE PUT IN THE ANJON VECTOR
C        AND THE PARTITION FUNCTIONS ARE PUT IN PART. IF IEPRO IS GREATER THAN
C        ZERO, THE GAS PRESSURE (PG,IN DYNES PER CM2), DENSITY (RO, IN GRAMS
C        PER CM3) AND INNER ENERGY (E, IN ERGS PER GRAM) ARE ALSO EVALUATED.
C        N O T E . RADIATION PRESSURE IS NOT INCLUDED IN E.
C
C        THE ENERGIES OF IONIZATION ARE REDUCED BY DXI, FOLLOWING BASCHEK ET 
C        AL., ABH. HAMB. VIII, 26 EQ. (10). THESE REDUCTIONS ARE ALSO MADE IN
C        THE COMPUTATION OF E.
C        THE ENERGY OF DISSOCIATION FOR H- HAS BEEN REDUCED BY 2*DXI, FOLLOWING
C        TARAFDAR AND VARDYA, THIRD HARV. SMITHS. CONF., PAGE 143. THE FORMATION
C        OF MOLECULES IS CONSIDERED FOR T LESS THAN TMOLIM.
C
C        IF IOUTR IS GREATER THAN ZERO, A DETAILED PRINT-OUT WILL BE GIVEN.
C
C
C        THE FUNCTION  QTRAV AND SUBROUTINE MOLEQ ARE CALLED.
C        THEY CALL QAS AND MOLFYS RESPECTIVELY.
C
C        DIMENSIONS NECESSARY
C        A(5),DQ(4),F(MAX(NJ)),PFAK(MAX(NJ)),RFAK(JMAX)
C        DIMENSIONS OF ARRAYS IN COMMONS /CI1/,/CI4/,/CI5/ AND /CI6/ ARE
C        COMMENTED ON IN SUBROUTINE INJON.
C        JMAX IS THE TOTAL NUMBER OF STAGES OF IONIZATION, INCLUDING NEUTRAL
C             ATOMS.
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I.
C
C
      include 'parameter.inc'
      character atnames*2, molnames*8
C
      DIMENSION DQ(4),F(5),PFAK(5),RFAK(45)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN(ndp),SUMH(ndp),XKBOL,NJ(16),IEL(16),
     *SUMM(ndp),NEL
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON/CI5/abmarcs(18,ndp),ANJON(18,5),H(5),PART(18,5),
     *DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON/CI6/TP,IQFIX(16,5),NQTEMP
      COMMON/CI7/A(5),PFISH,ITP
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/RABELL/XXRHO(NDP),XYRHO
      COMMON/CI8/YYPG,YYRHO,YYE
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
     *             ,FT,FTE
      COMMON/CMOL2/PK(33),NMOL
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      common /fullequilibrium/ partryck(ndp,maxmol),
     & xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partpx(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
C      common /fullequilibrium/ partryck(ndp,maxmol),
C     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet)
      COMMON/CPHYDRO/PHYDRO
      common /cg/grav,konsg /cteff/teff,flux
      COMMON /CMETPE/ PPEL(NDP), METPE
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
      
        INTEGER MOLH, JUMP      

C      REAL PHYDRO
C        STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)

CC     NOTE FOR PLANETS:     ! ADS, 31.07.2023
C      IF METPE=2 JON will use a simplified saha equation, which
C      results in a physically wrong electron pressure-pressure relation.
C      The physical electron pressure is the one from ggchem. The electron pressure
C      in solve is only useful for convergence of the model and SHALL NEVER BE USED ELSEWHERE
C
      ITP=1

C
C        IS T=THE TEMPERATURE OF THE PRECEDING CALL
      IF(ABS((T-TP)/T).LT.1.E-8)GO TO 53
   51 ITP=0
C
C        SOME QUANTITIES, ONLY DEPENDENT ON T
      TETA=5040./T
      TETA25=1.202D9/(TETA*TETA*SQRT(TETA))
      DO52 J=1,5
   52 A(J)=FL2(J)*TETA
C        A=ALFA(BASCHEK ET AL., CITED ABOVE)
C
      
      IF(NQTEMP.EQ.0)GO TO 53
C
C        PREPARATION FOR INTERPOLATION OF PARTITION FUNCTIONS IN T
      DQ(1)=1.
      DQ(2)=T-TPARF(1)
      DQ(3)=DQ(2)*(T-TPARF(2))
      DQ(4)=DQ(3)*(T-TPARF(3))
C
C        SOME QUANTITIES ALSO DEPENDENT ON PE
C        THE PFAK FACTORS ARE USED IN THE SAHA EQUATION. H(J) IS THE
C        QUANTUM NUMBER OF THE CUT OF THE PARTITION FUNCTIONS (ACCORDING
C        TO BASCHEK ET AL., CITED ABOVE) FOR J-1 TIMES IONIZED ATOMS. H IS
C        USED IN QAS.
C
C        XNEL= THE ELECTRON (NUMBER) DENSITY (PER CM3)
C        PFISH= P(FISCHEL AND SPARKS, ASTROPHYS. J. 164, 359 (1971)) IS USED IN
C        FUNCTION QAS.
C
   53 DXI=4.98E-4*TETA*SQRT(PE)
      DUM=TETA25/PE
      DIM=EXP10(DXI*TETA)
      PFAK(1)=DIM*DUM
      SQDXI=1./SQRT(DXI)
      H(1)=SHXIJ(1)*SQDXI
      DO54 J=2,5
      PFAK(J)=PFAK(J-1)*DIM
   54 H(J)=SHXIJ(J)*SQDXI
      XNEL=PE/(XKBOL*T)
      PFISH=4.2E3/XNEL**0.166666667
C
C        PARTITION FUNCTIONS AND IONIZATION EQUILIBRIA
C
      XNECNO=0.
      XNENH=0.
      EJON=0.
      JA=1

C     ADS: initialize calculation of pg-pe relation in simplified saha equation
      xpgpesumup = 0.0d0
      xpgpesumdown = 0.0d0
C
C        BEGINNING OF LOOP OVER ELEMENTS ('THE I-LOOP').
      DO24 I=1,NEL
      NJP=NJ(I)


C
C        SHOULD ELEMENT NO. I BE CONSIDERED
      IF(IELEM(I).GT.0)GO TO 9
      DO 55 J=1,NJP
      ANJON(I,J)=0.
      PART(I,J)=0.
   55 CONTINUE
      GO TO 23
C
C        BEGINNING OF LOOP OVER STAGES OF IONIZATION ('THE J-LOOP')
    9 DO19 J=1,NJP
      JM1=J-1

C
C        SHOULD STAGE OF IONIZATION NO. J BE CONSIDERED
      IF(ION(I,J).GT.0)GO TO 10
      ANJON(I,J)=0.
      PART(I,J)=0.
      GO TO 18
C
C        WHICH KIND OF PARTITION FUNCTION SHOULD BE COMPUTED

   10 IF(IQFIX(I,J)-1)14,11,13
   11 IF(T.LT.TPARF(1).OR.T.GT.TPARF(4))GO TO 13
      PARTP=PART(I,J)
      IF(ITP.GT.0)GO TO 15
C
C        PARTITION FUNCTIONS TO BE INTERPOLATED IN T
c      print *, "in interpolation of pf"
      JPARF=(JA-1)*4+1
      PARTP=0.
      DO12 IP=1,4
      PARTP=PARTP+PARQ(JPARF)*DQ(IP)
   12 JPARF=JPARF+1
      GO TO 15
C
C        PARTITION FUNCTIONS FOLLOWING TRAVING ET AL., ABH. HAMB. VIII,1 (1966)
  
   13 PARTP=QTRAV(TETA,H(J),J,JA)
      GO TO 15
C
C        THE PARTITION FUNCTION IS CONSTANT
   14  PARTP=PARCO(JA)
c      print *, "constant partition function"
   15 PART(I,J)=PARTP


C
C        IONIZATION EQUILIBRIA AND TOTAL NUMBER OF ELECTRONS
C
      IF(J.LE.1)GO TO 19
      IF(ITP.GT.0)GO TO 17
      RFAK(JA)=EXP10(-XIONG(I,JM1)*TETA)
   17 F(JM1)=PFAK(JM1)*RFAK(JA)*PARTP/PART(I,J-1)
      GO TO 19
   18 IF(J.GT.1)F(JM1)=0.
   19 JA=JA+1
C        END OF 'THE J-LOOP'
C

C     ADS: start here with code for planets, warning we are in loop of elements here
C     FOLLOWING GRAYS book, p.165  here

      IF(METPE.EQ.2) THEN 
        xphi = TETA25*PART(I,2)/PART(I,1)*EXP10(-TETA*XIONG(I,1)) 
        xmultfactor = (xphi/PE)/(1.0d0+xphi/PE)
C        xmultfactor = 0.01d0
        xpgpesumdown = xpgpesumdown + abmarcs(I,kl)*xmultfactor
        xpgpesumup = xpgpesumup + abmarcs(I,kl)*(1.0d0+xmultfactor)
      ENDIF

      FIL=1.
      DO20 J=2,NJP
      LL=NJP-J+1
   20 FIL=1.+F(LL)*FIL
      
      ANJON(I,1)=1./FIL
      XNEN=0.
      DO21 J=2,NJP
      JM1=J-1
      ANJON(I,J)=ANJON(I,JM1)*F(JM1)
      IF(I.LE.1)GO TO 24
      FLJM1=JM1
   21 XNEN=ANJON(I,J)*FLJM1+XNEN
      IF(I.GT.2.AND.I.LT.6) XNECNO=XNECNO+XNEN*abmarcs(I,kl)
      XNENH=XNEN*abmarcs(I,kl)+XNENH
C        XNENH=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN HYDROGEN (Q IN
C        MIHALAS, METH. COMP. PHYS. 7, 1 (1967), EQ. (35))
C        XNECNO=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN H, C, N, O

C
C
C        COMPUTATION OF THE ENERGY OF IONIZATION (EJON). HYDROGEN IS NOT
C        INCLUDED.
C
      XERG=0.
C        XERG= THE ENERGY OF IONIZATION PER ATOM (IN ELECTRON VOLTS)
C
      DO22 J=2,NJP
      JM1=J-1
      FLJM1=JM1
   22 XERG=ANJON(I,J)*(XIONG(I,JM1)-DXI*FLJM1)+XERG
      EJON=XERG*abmarcs(I,kl)+EJON
      GO TO 24
   23 JA=JA+NJP
   24 CONTINUE
C        END OF 'THE I-LOOP'
C
C
      XNECNO=XNENH-XNECNO
C now xnecno is all that is not from H, C, N, O.
      TP=T
      IF(IEPRO.LE.0)GO TO 71
C
C        COMP. OF PRESSURE, DENSITY AND INNER ENERGY
C
      XIH=XIONG(1,1)-DXI
      XIHM=0.747-2.*DXI
C        XIH AND XIHM ARE THE ENERGIES OF IONIZATION FOR H AND H- RESPECTIVELY
C        (IN ELECTRON VOLTS).
C
      XKHM=TETA25*2.*EXP10(-TETA*XIHM)
C        XKHM = THE 'DISSOCIATION CONSTANT' FOR H-.
C
      HJONH=ANJON(1,2)/ANJON(1,1)
C

*********
C 15.12.94 Ch.Helling
C
C JUMP is set in JONABS.DAT JUMP = 1 => MOLH = 1 :
C                                       call only MOLEQ
C                                       = only H2, H2+ as molecules cosidered
C                                                                            
C( gives the possibility to use only TSUJI-code instaed of MOL for the molecular
C  equilibrium )
C
C                           Jump = 0 : call MOL
C                                      = all molecules cosidered in old MARCS-routine
C                                        MOL
C                                        ( there are problems to calculate lower 
C                                          temperatures )
C
C 25/1/01, UGJ:
C JUMP = 2: JFF routine corresponding to Tsuji's method
C JUMP = 3: JFF Gibs minimalisation method (892 molecules, ions and atoms).
C JUMP = 4: Added GGchem code by ERC
*********

C     ADS: We use a simplified saha equation for planets, when we use metpe=2, better for convergence
      IF(METPE.EQ.2 .and. t.le.2000.0) then
2421  format('teff,t,pe,ppel,pg = ',2f8.1,3e12.3)   
      GO TO 421
      end if
C     JUMP to no molecules, when temperature is highre than tmolim      
      IF(T.GT.TMOLIM)GO TO 42
C      IF(JUMP.ge.1) MOLH=1
* the former step is necessary because often is MOLH=0 on other places !
      IF(MOLH.LE.0) GO TO 45


*********
C
C => if MOLH < or = 0 than is MOLEQ allways " overjumped " ( !! )
C                     and a more complete moleculare formation will
C                     be evaluated in MOL
C => if MOLH > 0  then only MOLEQ is supposed to work and only the 
C           (= 1) moleculare formation of H2 and H2+ is considered 
C                 ( => total hydrogen pressure [H, H2, H2+] 
C                   => possibility to use only Tsujis eqilibrium routines EQMOL,
C                      called in TEST_TSUJI )
*********
C
C        FORMATION OF MOLECULES. ONLY H2 AND H2+
      
   41 PRINT *, "moleq is called"
       CALL MOLEQ(T,PE,HJONH,XIH,XKHM,XIHM,XNENH,F1,F2,F3,F4,F5,FE,FSUM,
     *   EH)

*
* F1 = N(HI)/N(H)   ==> F1P = P(HI)
* F2 = N(HII)/N(H)  ==> F2P = P(HII)
* F3 = N(H-)/N(H)   ==> F3P = P(H-)
* F4 = N(H2+)/N(H)  ==> F4P = P(H2+)
* F5 = N(H2)/N(H)   ==> F5P = P(H2)
*

      FEPE=PE/FE
      F1P=F1*FEPE
      F3P=F3*FEPE
      F4P=F4*FEPE
      F5P=F5*FEPE
      PHYDRO=FSUM*PE/FE
      GO TO 43
C        FORMATION OF MOLECULES COMPOSED OF H,C,N,O
   45 IF(ANJON(3,1).LE.0..OR.ANJON(4,1).LE.0..OR.ANJON(5,1).LE.0.)
     * GOTO 41
*      print*,' now I am in MOL ,  jumped over MOLEQ'
      HJONC=ANJON(3,2)/ANJON(3,1)
      HJONN=ANJON(4,2)/ANJON(4,1)
      HJONO=ANJON(5,2)/ANJON(5,1)
      ABUC=abmarcs(3,kl)/abmarcs(1,kl)
      ABUN=abmarcs(4,kl)/abmarcs(1,kl)
      ABUO=abmarcs(5,kl)/abmarcs(1,kl)

      CALL MOL(T,PE,HJONH,HJONC,HJONN,HJONO,ABUC,ABUO,ABUN,XIH,XKHM,XIHM
     *,XNECNO,F1,F2,F3,F4,F5)

      SUMPMO=0.
      PRESMO(1)=FHE*PK(1)
      PRESMO(2)=FHE*FHE*PK(2)
      PRESMO(3)=FHE*FHE*HJONH*PK(3)
      PRESMO(4)=FHE*FHE*FOE*PK(4)
      PRESMO(5)=FHE*FOE*PK(5)
      PRESMO(6)=FHE*FCE*PK(6)
      PRESMO(7)=FCE*FOE*PK(7)
      PRESMO(8)=FCE*FNE*PK(8)
      PRESMO(9)=FCE*FCE*PK(9)
      PRESMO(10)=FNE*FNE*PK(10)
      PRESMO(11)=FOE*FOE*PK(11)
      PRESMO(12)=FNE*FOE*PK(12)
      PRESMO(13)=FNE*FHE*PK(13)
      PRESMO(14)=FCE*FCE*FHE*FHE*PK(14)
      PRESMO(15)=FHE*FCE*FNE*PK(15)
      PRESMO(16)=FCE*FCE*FHE*PK(16)
      PRESMO(17)=0.0
      PRESMO(18)=FHE*FSE*PK(18)
      PRESMO(19)=FKE*FHE*PK(19)
      PRESMO(20)=FCE*FCE*FCE*FHE*PK(20)
      PRESMO(21)=FCE*FCE*FCE*PK(21)
      PRESMO(22)=FCE*FSE*PK(22)
      PRESMO(23)=FKE*FCE*PK(23)
      PRESMO(24)=FKE*FCE*FCE*PK(24)
      PRESMO(25)=FNE*FSE*PK(25)
      PRESMO(26)=FKE*FNE*PK(26)
      PRESMO(27)=FKE*FOE*PK(27)
      PRESMO(28)=FSE*FOE*PK(28)
      PRESMO(29)=FSE*FSE*PK(29)
      PRESMO(30)=FKE*FSE*PK(30)
      PRESMO(31)=FTE*FOE*PK(31)
      PRESMO(32)=FTE*FOE*FOE*PK(32)
      PRESMO(33)=FTE*FCE*FCE*PK(33)
      DO 30 I=1,NMOL
      PRESMO(I)=PRESMO(I)*PE
   30 SUMPMO=SUMPMO+PRESMO(I)
      
C SUMPMO = sum partial pressures [dyn/cm2] of all molecules
C SUMPA = sum partial pressures of neutral atoms H,C,N,O not in molecules
C SUMPI = sum partial pressure from ionized H,C,N,O
C SUMM=SUMM-ABUND(1)-ABUND(3)-ABUND(4)-ABUND(5) in INJON is number
C of nuclei other than H,C,N,O per H atom (abund are normalized to abund(1) ).
C PE/FE = ro*k_boltz*T/(XMY*XMH) = conversion from #atoms/#H to P[dyn/cm3]
C PE*SUMM/FE = sum partial pressures from all nuclei not HCNO.

      SUMPA=PE*(FHE+FCE+FNE+FOE)
      SUMPI=PE*(FHE*HJONH+FCE*HJONC+FNE*HJONN+FOE*HJONO)
      HNIC=PE*FHE
      HPNIC=HNIC*HJONH
      PG=PE+SUMPMO+SUMPA+SUMPI+PE*SUMM(kl)/FE
      P_MOL(KL) = SUMPMO
      P_NEU_HCNO(KL) = SUMPA
      P_ION_HCNO(KL) = SUMPI
      P_NEU_HE(KL) = 
     &   PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,1)
      P_ION_HE(KL) = 
     &   PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,2)
      PHE = P_ION_HE(KL) + P_NEU_HE(KL)
      P_NON_HHECNO(KL) = PE*SUMM(kl)/FE-PHE
      PG_JON(KL) = PG
      HN_JON(KL) = 1./(XMH*XMY(kl))
      P6_JON(KL) = HNIC+0.42*PHE+0.85*PRESMO(2)
      GOTO 46
C
C     NO MOLECULES (high temperature)
   42 CONTINUE
      F2=ANJON(1,2)
      FE=XNENH+F2
      F1=ANJON(1,1)
      F3=0.
      F4=0.
      F5=0.
      FSUM=1.
      EH=-XIH*F1

      GOTO 43

  421 CONTINUE
C     planet (metpe=2, fake electron pressure)
C     Warning: at low temperature, we use a superficial/simplified saha relation for the pe-pg relation
      F1 = ANJON(1,1)
      F2 = 0.0d0
      F3 = 0.0d0
      F4 = 0.0d0
      F5 = 0.0d0
C This relation will be calculated:      
C PG = PE + PE*FSUM/FE + PE*SUMH/FE   ( = line 43 )      
C but we have pg = pgpesumup/pgpesumdown * pe
C setup FE so that equations match:
      FSUM=0.0d0
      FE = SUMH(kl)/(xpgpesumup/xpgpesumdown+1.0d0)
      EH=-XIH*F1

      PG = xpgpesumup/xpgpesumdown * PE
      !write(*,*) kl,PG, PE, PPEL(kl)
      !write(*,*) xpgpesumup, xpgpesumdown ,xpgpesumup/xpgpesumdown 
      !write(*,*) xpgpesumup/xpgpesumdown * PPEL(kl)   

**********
C
C ==>
C
C PG = PE + PE*FSUM/FE + PE*SUMH/FE   ( = line 43 )
C    = PE + PHYDRO + PE*SUMH/FE
C
C SUMH = number of other nuclei / number of hydrogen nuclei
C   FE = PE/PH   with PH = NH*kT ( NH number of hydrogen nuclei per cm3 )
C   => FE = number of other nuclei * kT
C         = ficticious pressure of the other nuclei, except hydrogen and e-
C
C   => PE*SUMH/FE = gaspressure contributed by all the atomic species
C
C   PHYGRO = P(HI) + P(HII) + P(H-) + P(H2+) + P(H2) 
C          = P(H)  + P(H+)  + P(H-) + P(H2+) + P(H2)
C 
**********

 43   CONTINUE 

46    continue
C
 

      RO=PE*XMY(kl)*(XMH/XKBOL)/(FE*T)

      IF (JUMP.GE.1 .and. MOLH.eq.1) THEN

       PG=PE*(1.+(FSUM+SUMH(kl))/FE)
           XNHE = abmarcs(2,kl) / (XMH*XMY(kl))
           PHE = ro * 1.38053e-16 * T * XNHE
       P6_JON(KL) = xmettryck(kl,1)+0.42*PHE+0.85*partryck(kl,2)
      
      P_NEU_HE(KL) = PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,1)/
     &  (ANJON(2,1)+ANJON(2,2))
      P_ION_HE(KL) = PE*abmarcs(1,kl)*abmarcs(2,kl)/FE*ANJON(2,2)/
     &  (ANJON(2,1)+ANJON(2,2))
      PHE = P_ION_HE(KL) + P_NEU_HE(KL)
      PG_JON(KL) = PG
      HN_JON(KL) = 1./(XMH*XMY(kl))

 
      ENDIF


      XYRHO=RO
      E=1.5*PG/RO+(EH+EJON)*ENAMN(kl)
      YYPG=PG
      YYRHO=RO
      YYE=E
      RO_JON(KL) = RO
C partial pressure of He is put into presmo(17)
        XNHE = abmarcs(2,kl) / (XMH*XMY(kl)) 
        PRESMO(17) = 1.38053d-16 * T * XNHE * RO



      IF(IOUTR.LE.0)GO TO 71
C
C        **** PRINT-OUT ****
C
* iwrit = 6 = mxms7.out
      WRITE(IWRIT,204)T,PE,PG,RO,E
      WRITE(IWRIT,201)
      WRITE(IWRIT,202)
      DO93 I=1,NEL
      NJP=NJ(I)
      WRITE(IWRIT,203)IEL(I),abmarcs(I,1),(ANJON(I,J),J=1,NJP)
      WRITE(IWRIT,207)(PART(I,J),J=1,NJP)
   93 CONTINUE
      IF(T.GT.TMOLIM)GO TO 44
      IF(MOLH.LE.0)GOTO 47
      WRITE(IWRIT,205)F1P,F3P,F5P,F4P
      GO TO 71
   47 CONTINUE

***************18.12.94 Ch.H

      IF (JUMP.EQ.0) THEN
       WRITE(IWRIT,208)HNIC,(PRESMO(I),I=1,13)
       ELSE
       WRITE(IWRIT,209) xmettryck(kl,1),(partryck(kl,I),I=1,2),
     &  (PARTRYCK(kl,I),I=4,13)
      ENDIF

************************

      GOTO 71
   44 WRITE(IWRIT,206)
C
   71 CONTINUE
c
c      print*,'end of jon. fo,fe,foe ',fo,fe,foe
c
      RETURN
  201 FORMAT(1H0,'ELEMENT  ABUNDANCE  IONIZATION FRACTIONS',17X,
     *'PARTITION FUNCTIONS')
  202 FORMAT(1H ,23X,1HI,7X,2HII,6X,3HIII,5X,2HIV,12X,1HI,9X,2HII,8X,
     *3HIII,7X,2HIV)
  203 FORMAT(6H      ,A2,E12.4,4F8.4)
  204 FORMAT(3H0T=,F7.1,5X,3HPE=,E12.4,5X,3HPG=,E12.4,5X,3HRO=,E12.4,
     *5X,2HE=,E12.4)
  205 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+'/1X,
     *4(1PE9.2))
  206 FORMAT(1H0,'NO MOLECULES CONSIDERED IN MARCS-ROUTINES -- ', 
     *' Tsuji routine for molecules used')
  207 FORMAT(1H+,56X,4E10.3) 
  208 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+',6X,
     *'H2O',6X,'OH',7X,'CH',7X,'CO',7X,'CN',7X,'C2',7X,'N2',7X,'O2',7X,
     *'NO',7X,'NH'/1X,14(1PE9.2))
 209  FORMAT(1HO,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',6X,
     *'H2O',6X,'OH',7X,'CH',7X,'CO',7X,'CN',7X,'C2',7X,'N2',7X,'O2',7X,
     *'NO',7X,'NH'/1X,13(1PE9.2))
      END
C
      SUBROUTINE KAP5(T,PE,ABSK,nlayer)
      implicit real*8 (a-h,o-z)
C
      COMMON /CXLSET/XL(20,10),NSET,NL(10)
C
C COMPUTE KAPPA(5000.). 73.10.17 *NORD*.
      !print *, "absko call in kap5 "
      CALL ABSKO(1,1,T,PE,1,NL(1)+1,ABSK,SPRID,nlayer)
      
      RETURN
      END
C
       FUNCTION LENSTR(STRING)
      implicit real*8 (a-h,o-z)
*
* Returns the length of a string not counting trailing blanks
*
       CHARACTER*(*)  STRING
*
       DO 10 I = LEN(STRING), 1, -1
         IF(STRING(I:I) .NE. ' ') THEN
           LENSTR = I
           RETURN
         ENDIF
   10  CONTINUE
       LENSTR = 0
       RETURN
       END

      SUBROUTINE LISTMO(MO,IARCH,ISPH)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE READS AND PRINTS A MODEL IN A READABLE FORMAT
C        AFTER IT HAS BEEN STORED ON FORTRAN UNIT IARCH.
C
      include 'parameter.inc'
      DIMENSION TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(33,NDP),FCONV(NDP),RR(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(20),IEL(16),PROV(20,20),PALL(NDP)
     *,ABSKA(20),SPRIDA(20),XLB(500),PEP(16)
     *,ABNAME(30),SOURCE(30),ABSKTR(NDP),SPRTR(NDP)
C     *,ABNAME(30),SOURCE(30),ABSKTR(NDP),SPRTR(NDP),DUMMY(11*NDP+2)
      DIMENSION W(500),UW(12),BW(21),VW(25),SUMPMOL(NDP)
      dimension xdp(8),kdp(8)
C      REAL*8 ROSSO,PTAUO
      COMMON /CVAAGL/XLB,W,NLAM
      COMMON /CARC2/T(NDP),FC(NDP),FLUXME(NWL),TAU5(NDP),INORD
      CHARACTER*10 DAG,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10)
      COMMON /COPINF/ SUMOP(maxosmol,NDP),SUMKAP(maxosmol,NDP)
      COMMON /UTPUT/IREAD,IWRIT
      COMMON/CI5/abmarcs(18,ndp),ANJON(18,5),H(5),PART(18,5),
     * DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY(ndp)
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLEVPRINT/ PRJ2(NDP),masslinf
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON/COSWR/osresl,losresl,listwn
C      COMMON/COPPR/oppr(15,3,120,3),jvxmax,itxmax  !15mol,10dpt,100wn
C      COMMON/COPPRR/xconop(120,10),xlineop(120,10)    !100wn,10dpt
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CLIST/NLTE
      COMMON /MASSE/RELM
      COMMON /CLIN/lin_cia
      COMMON /CNEWC3 /NEWC3
      COMMON /CG/GRAV,KONSG
      COMMON /CSTYR/MIHAL,NOCONV
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom         !irrin=1~comp.irrad,steff=rad*
      COMMON /CXMAX/XMAX /CTAUM/TAUM
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX                          
      COMMON /CPOLY/FACPLY,MOLTSUJI
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500)
      COMMON/COPsum/ SSUM(NDP),XSUM(NDP),CONSUM(NDP)
      COMMON/CI4/ TMOLIM, IELEM(16), ION(16,5), MOLH, JUMP
      COMMON /ROSSC/CXKAPR(NDP),CROSS(NDP)
      COMMON /COSEXP/ LOPS,NOPS
      DATA A,B/.34785485,.65214515/
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)
      COMMON /CKMOL/KMOL(MAXOSMOL)   !connects OS-molecule with presmo-index
      data xdp / -5., -4., -3., -2., -1., 0., 1., 2./
      character*5 name_mol, name_listmo
      COMMON /CMOLNAME/NAME_MOL(maxmol),NAME_LISTMO(maxmol)
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO_ST(NDP),
     * NTAU,ITER
      common /cgem/pres_gem(ndp,nspec)
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character name_gem*8
      common /cabink/abink(ndp,nspec)
      dimension phe(ndp), trpe(ndp), trphe(ndp)
      dimension printpp(ndp,nspec), printname(nspec)
      character printname*8
      dimension listpp(nspec)
      dimension ptot(ndp),pp_sum(ndp)
      dimension z(ndp)
      dimension pe_gem(ndp)
      COMMON /CMETPE/ PPEL(NDP), METPE
      common /dpeset/ dpein,dtin
      COMMON /CORRECT/TDIFF,TCONV,KORT
      common /ctcorlast/tcorlast   !the maximum temperature correction in last iteration 
      common /cu2warning/nu2warning
      common /ggchemmu/ggmu(NDP),ggrho(NDP),ppsum(ndp),ppappsum(ndp),
     &   ppnonappsum(ndp),tg(ndp),pges(ndp)
     &  ,ppat1sum(ndp),ppat2sum(ndp),ppmolsum(ndp),ppgs(ndp)
      
      character atnames*2, molnames*8
      
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
      dimension ppm(ndp),ppp(ndp),ppnmol(ndp),ppamol(ndp),ppnat(ndp)
     &     ,pp24(ndp),pp54(ndp),pp24m(ndp),pp54m(ndp)
      common /consistlist/pgg(ndp)
      common /consistgginit/ttgg(ndp),ppgg(ndp),kk
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      character*3 aifix
      common /pefix/xionfix,aifix
      dimension idplus(543),idminus(543),idmol(543),idos(543)
      character molnam4*4,molnam8*8
      character*4 osnames(54) !UNHARD CODE THIS
      real*8,dimension(ndp,543):: dummy_ppallmol
      real*8,dimension(ndp,22):: dummy_ppallat

C The 54 molecules we have absorption coefficients for (March 2026) are:
      data osnames /
     & 'C2  ','CAH ','CH  ','CN  ','CO  ','CO2 ','FEH ','H2O ','HCN ',
     & 'MGH ','NH  ','OH  ','SIH ','SIO ','TIH ','TIO ','CRH ','NO  ',
     & 'LIH ','VO  ','ZRO ','H2  ','C2H2','C3  ','ALCL','ALF ','ALH ',
     & 'ALO ','BEH ','CAF ','CH3F','CH4 ','CP  ','CS  ','H2CO','HCL ',
     & 'HNO3','KCL ','KF  ','LICL','LIF ','MGF ','NACL','NAF ','NAH ',
     & 'NH3 ','SN  ','PH3 ','PN  ','PO  ','PS  ','HS  ','SIS ','SO2 '/
C A full list of all molecules and atoms that can be included in the output can be found in at_and_mol_list.txt
C the 15 of the 22 neutral atoms we list in the output (...dat) are:
      dimension kpratoms(15)
      data kpratoms /1, 2,3,4,5,6, 7, 8,10,11,12,13,14,15,17/
C                    H He C N O Na Mg Si Fe Al Ca Cr Ti S  K
      dimension kions(15)
      data kions/15,46,389,414,453,456,462,22,55,447,454,463,466,2*0/
C               H2+ OH+ AL+ CA+ K+ Mg+ Na+ H2- OH- H- K- NA- O- sum+ pe 
C instead of the 3X15 ions we listed until "the-web-output-form" of
C marcs, we now list all the ones from GGchem in two large arrays,
C kionsm and kionsp, the negative (...m) and the positive (...p) ions, respectively.
C We list the sum of the positive (sum+) and negative (sum- and Pe),
C sum+, sum-, pe, separately in orde to make the electron budget clear.
C As a test, sum+ should be equal to Pe + sum-
C  There are 82 positive ions computed in GGchem:
      dimension kionsp(82)
      data kionsp /15,16,17,18,19,20,21,42,43,44,45,46,47,48,49,
     & 50,51,52,53,65,66,78,79,97,98,99,389,391,392,393,
     & 395,396,400,402,405,406,407,409,410,411,412,414,415,416,417,
     & 418,420,421,423,425,426,428,431,433,435,437,439,441,444,446,
     & 448,449,450,451,452,453,456,457,461,462,465,468,470,472,474,
     & 476,479,480,486,517,520,530/
C There are  52 negative ions:
      dimension kionsm(52)
      data kionsm /22,23,54,55,56,57,67,68,100,390,394,397,398,399,401,
     & 403,404,408,413,419,422,424,427,429,430,432,434,436,438,440,
     & 442,443,445,447,454,455,458,459,460,463,464,466,467,469,471,
     & 473,475,477,510,521,525,531/      
      dimension k1mol(15)
      data k1mol/1,3,4,5,26,58,59,192,350,185,354,173,205,214,92/
C   H2 C2 N2 O2 CH CN CO CO2 H2O CH4 NH3 HCN C2H2 C3 TiO
      dimension k2mol(15)
      data k2mol/28,31,32,35,36,122,37,38,222,39,40,41,86,61,212/
C   OH MgH AlH HS HCl NaCl KH CaH CaOH TiH CrH FeH SiO SiC SiC2
      dimension k3mol(15)
      data k3mol/63,93,94,95,96,115,141,337,336,353,351,355,374,376,532/
C   CS VO CrO FeO ZrO SiS AlOH KOH FeS H2S H2SO4 PH3 TiO2 O3 TiC

*
      FLUMAG(I)=-2.5*log10(FLUXME(I))-STMAGN
*
      ispec =1
      if (ispec.eq.1) then
       open (unit=29,file='spectrum.dat',status='unknown')
       do 290 i=1,nwtot
       write(29,295) (spec(i,k),k=1,3)
290    continue
       close(29)
295    format(f10.2,1p2e12.4)
      end if
C     Iniation of natural constants for Kzz calculation
      HPLNCK=6.62554E-27
      BOLTZK=1.38054E-16
      CLIGHT=2.997925E10
      ECHARG=4.80298E-10
      RYDBRG=1.097373E5
      STEFAN=5.675E-5
      CLIGHT=2.99793E10
      PI=3.14159265
      PI4C=PI*4./CLIGHT

C Calling ROSSOS here makes it possible to list X and S as function of
C wavelength and depth in the model, but it comes with the cost that
C OPAC is called for each (continuums)wavelength and depth. OPAC calls
C OSTABLOOK which runs the whole setup of GGchem.
C Probably that's why it takes up all the computing time from the last
C iteration to the printing of the model. Maybe it could be avoided?
      

      
      CALL ROSSOS

      IREAD=5
      IWRIT=7
C
      REWIND IARCH
      DO 1 IMO=1,MO
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(abmarcs(I,1),I=1,NEL)

      GLOG=log10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX

      WRITE(7,1011) palfa,pbeta,pny,py,NOCONV 
1011  format(/'Convection was computed based on the',
     & ' parameter choice:',
     & /'palfa, pbeta, pny, py =',3f4.1,f6.3,
     & /'(but convection was apriori considered irelevant for the',
     & ' uppermost',i3,' layers)')

      relmj = 1048.0*relm !mass of object in units of Jupiter-masses,  1.989e30/1.898e27
      relme = 318.0*relmj !mass of object in units of Earth-masses 
      IF (isph.eq.1) THEN 
      write(7,*) ' '
      write(7,*)'The model is computed in spherical geometry,',
     & ' based on the following parameters:'
      if(relm.ge.0.08) then       !mass greater than 0.08 M_sun = 80 M_jup
         write(7,1531)relm
      else if(relm.ge.0.0003) then !=0.3Mjup = M_saturn; 0.001 Msun = 1 Mjup
         write(7,1532)relmj        !i.e. M_saturn<relm<80M_jup
      else 
         write(7,1533)relme        !i.e. relm<M_saturn
      end if
      write(7,412) MMY,NCORE,KDIFF

412   format('The ray tracings parameters mmy,ncore,kdiff:',3i2)
1531  format('mass (in units of solar mass):', f6.3) 
1532  format('mass (in units of jupiter mass):', f6.3) 
1533  format('mass (in units of earth mass):', f6.3) 
      ENDIF
      IF (irrin == 1) THEN

      write(7,516)
516   FORMAT(/'The model takes irradiation into account,',
     & ' based on the parameters:')
      write(7,401) steff
401   format("Stellar effective temperature ", f8.0, " Kelvin ")
      write(7,402)rstar
402   format("Stellar radius of  ", f7.2, " solar radii")
      write(7,403)semimajor
403   format("Planet at ",f6.3, " AU from star") 
      write(7,406) f_irrad
406   format('irradiation factor:',f6.2)
      write(7,407) teffp
407   format('final planetary Teff:',f8.1)

      if (irrinp == 1) then
      write(7,404) tbottom, reflect
404   format("The model is a rocky planet with surface temperature ", 
     * f5.0, " K and surface albedo ",f4.2)
      else
      write(7,*) 'the irradiated model is a gas-planet, a brown dwarf',
     * ' or star and therefore has no "reflecting surface" '
      end if

      END IF   

      write(7,1540) tconv,abs(tcorlast)
1540  format(/'The convergence criterium was delta(T)=',f4.1,' K and ',
     & 'temp. corr. in the last iteration was',f5.2,' K'//)

      write(7,*) 'Technical settings:'
      write(7,517) dpein,dtin,metpe,xionfix,kort,kpp,mihal,lops,nops
517   FORMAT('dpein,dtin,method-pe,inz.K=12,kort,kpi,mihal,lops,nops:'
     &    2F6.3,i2,f6.3,'eV',i2,2i3,i2,i3)

      if (lops.ne.0  .or.  nops.ne.1) write(7,2016)
2016  format (' the os-absorption was shifted',i3,' os-steps,',
     *    ' and only each',i4,' os-value was used (rest==0)')

      if(metpe.eq.2) write(7,*)
     &  '(listed Pe is from GGchem; Pe for converg. in colm Pturb)'
      write(7,2020)
2020  FORMAT(/1X,81('*')//)
      if (nu2warning .eq. 1) then
         write(7,*) 'WARNING: there were problems with convergence of',
     &   'some thermodynamical values (in TERMO), so do not use, or '
         write(7,*) 'at least use with care, Cp, Cv etc',
     &   '; see output block thermodynamical values below'
      endif
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
      READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &            PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &            AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &            (PRESMO(J,K),J=1,NMOL)
      TAUK=log10(TAU(K))+10.01
      KTAU=TAUK
      IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
      IF(KTAU.EQ.10) K0=K
      NTPO=NTPO+1
   31 CONTINUE
    3 CONTINUE
      WRITE(7,204)
      WRITE(7,205)
      DO 4 I=1,JTAU
      FCONV(I)=ANCONV(I)*FLUX
      WRITE(7,206) I,log10(tau(i)),
     *      TAU(I),T(I),TKORRM(I),FCONV(I),FCORR(I),I
    4 CONTINUE
C
      IF (NLTE.EQ.0) GO TO 4000
C*
C* 90-05-13 START OF MODIFICATIONS (MATS CARLSSON)
C* PRINT MULTI ATMOSPHERIC FILE: ATMOS.MULTI
C*
      OPEN(33,FILE='ATMOS.MULTI',STATUS='NEW',CARRIAGE CONTROL='LIST')
      WRITE(33,400) TEFF
  400 FORMAT(' MARCS MODEL ATMOSPHERE, TEFF=',F10.2/' TAU(5000) SCALE')
      WRITE(33,410) G
  410 FORMAT('*'/'* LG G'/F6.2)
      WRITE(33,420) JTAU
  420 FORMAT('*'/'* NDEP'/I3)
      WRITE(33,430)
  430 FORMAT('*'/'*LG TAU(5000)    TEMPERATURE        NE         V',
     * '              VTURB')
      DO 450 I=1,JTAU
        IF(TAUS(I).GT.0.0) THEN
          TAULG=LOG10(TAUS(I))
        ELSE
          TAULG=2.*LOG10(TAUS(I+1))-LOG10(TAUS(I+2))
        ENDIF
        WRITE(33,440) TAULG,T(I),PE(I)/T(I)/1.380662E-16,0.,2.
  440   FORMAT(1P,5E14.6)
  450 CONTINUE
      CLOSE(33)
C*
C* 90-05-13 END OF MODIFICATIONS
C*
4000  CONTINUE


      WRITE(7,208)
      Z0=Z(1)
      DO 5 I=1,JTAU
      Z(I)=Z(I)-Z0
      if(metpe.eq.1) then
      WRITE(7,209) I,TAU(I),TAUS(I),Z(I),T(I),PE(I),PG(I),PRAD(I),
     &             PTURB(I),XKAPR(I),I
      else if(metpe.eq.2) then
      WRITE(7,209) I,TAU(I),TAUS(I),Z(I),T(I),PPEL(I),PG(I),PRAD(I),
     &             PE(I),XKAPR(I),I
C     &             PTURB(I),XKAPR(I),I
       end if
        !pgx=PP(i)-PPR(i)-PPT(i)
        ppallsum=ppappsum(i)+ppnonappsum(i)+ppat1sum(i)+ppel(i)
     
c       IF (T(I).GT.TEFF) then
c          iint = i
         
c          if ((t(i)-teff .gt. teff-t(i-1)) .and. 
c      &    (iint /= 1)) iint = i-1
         
c       END IF
    5 CONTINUE

2085  FORMAT(' K',1X,'lg(tau)',2X,'T',3X,'Pe-solve',1X,
     & 'Pe-GGchm',1x,'Pg-solve',1x,'PPapsum',2x,'PPnon-ap',1x,
     & 'PPmolsum',1x,'PPatsum',2x,'PPallsum',3x,
     & 'roGG',3x,'muGG')
2095  FORMAT(I2,F6.2,F7.1,1P9E9.2,0pf5.2)
     

C      masslinf = 1                    !now (sept.2006) in namelist outlist 
      if (masslinf.eq.0) go to 4001
      WRITE(7,2071)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,2081)
      DO 51 I=1,JTAU
      IF (I.LE.JTAU-1) THEN
      ROMIDT = (RO(I)+RO(I+1) )/2.
      DPGDZ = ( PG(I+1)-PG(I) ) / ( Z(I)-Z(I+1) ) / ROMIDT
      DPRDZ = ( PRAD(I+1)-PRAD(I) ) / ( Z(I)-Z(I+1) ) / ROMIDT
      END IF
      RO(I) = MAX (RO(I),1.D-99)
      RI = SQRT(RELM/G) * 1.152E13 + Z(JTAU) - Z(I)
      RILOG = LOG10(RI)
      VLOG = log10(AMLOSS) + 24.700 - 2.*RILOG - LOG10(RO(I))
      VLOSS = 10.**VLOG *1.e-5    !velocity in km/s
      ZINT = Z(IINT) - Z(I)
      WRITE(7,2091) I,TAU(I),LOG10(TAU(I)),RO(I),ZINT,T(I),PPRG(I),
     *         DPRDZ,VLOSS,ROSSO(I),PTAUO(I),cross(I)
   51 CONTINUE
4001  CONTINUE

      WRITE(7,210)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,211)
      DO 6 I=1,JTAU
      WRITE(7,212) I,TAU(I),RO(I),EMU(I),CP(I),CV(I),AGRAD(I),Q(I),U(I),
     &             V(I),ANCONV(I),I
    6 CONTINUE
C

      DO 772 I=1,JTAU
        HNIC(I) = max(1.d-99,HNIC(I))
        HNIC(I)=LOG10(HNIC(I))
           SUMPMOL(I) = 0.
        DO 771 J=1,33
           PRESMO(J,I)=max(1.d-99,PRESMO(J,I))
           if(j.eq.1.or.j.eq.17.or.j.eq.32.or.j.eq.33) then
                 PRESMO(J,I)=LOG10(PRESMO(J,I))
                 go to 771   !(P_He is in presmo(17), 1 is H-, 32,33 extr)
           end if
           SUMPMOL(I) = SUMPMOL(I) + PRESMO(J,I)
           PRESMO(J,I)=LOG10(PRESMO(J,I))
771      CONTINUE
772   CONTINUE
C Molecular partial pressures from the Marcs equilibrium
C UGJ 14/11/98: Test with partial pressures from the two
C routines show them almost identical for giants (for very cool
C giants and dwarfs they are probably a bit different from 
C one another). You can write both sets by deleting the first
C "if(jump.ne.1)" loop. but be aware that the spectrum program
C may have difficult finding P(ZrO), P6 and other part.pressures then.
C     
      IF (JUMP.GT.0) GO TO 2101


      WRITE(7,213)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,214)

      DO 7 I=1,JTAU
        WRITE(7,215) I,HNIC(I),(PRESMO(JJ,I),JJ=1,13),PRESMO(31,I),I
7     CONTINUE


      WRITE(7,213)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,216)
      DO 8 I=1,JTAU
      WRITE(7,217) I,(PRESMO(J,I),J=14,16),(PRESMO(J,I),J=18,30),I
8     CONTINUE


      GO TO 2109

2101  CONTINUE     ! go here if jump>0 (i.e. not old marcs chem. equilibrium.

      IF (JUMP.EQ.3) GO TO 2102

C here: JF's/Tsuji's molecular partial pressures from JANAF polynomial fits:

C for Tsuji's eq names are not read in (MOL(J) as real*8...)

      IF (jump.eq.4) THEN

C Here comes output blocks for selected of the 22 atoms and 543 molecules directly from GGchem; 
C Select which by chancing the index numbers in the data statements kions, kpratoms, k1mol, k2mol, k3mol
C in the top of the LISTMO subroutine ; do not put more than 15 species per output block.
C With time these ought to be the only partial pressure listing in the atmosphere data files.
C Unit=90 is an output file with more infor on the partial pressures than stored in the model atmosphere files.


C        open(unit=90,file='allpp.dat',status='unknown')


C The OS-molecules:
        
        nos = 0
        do 4201 k=1,54
        do 4200 j=1,543
        ipos = 0
        molnam8 = molnames(j)
        molnam4 = molnam8(1:4)
        ipos = index(molnames(j),osnames(k))
        if (ipos.ne.1) go to 4200
        nos = nos + 1
        idos(nos) = j
        go to 4201
4200    continue
4201    continue
C
4218    format('there were no GGchem partial pressure for OS molecule:')
        write(90,*) 'noggos,k,ipos, osnames(k)'
4219    format(3i4,2x,a4)
        noggos = 0
        write(90,4218)
        do 4226 k=1,54
        do 4225 j=1,nos
        ipos = 0
        ipos = index(molnames(idos(j)),osnames(k))
        if (ipos.ne.1) go to 4225
        go to 4226
4225    continue
        noggos = noggos + 1
        write(90,4219) noggos,k,ipos, osnames(k)
4226    continue
      WRITE(90,4232) NOSMOL
4232  format(' Following ',i3,
     &       ' molecules are included in opacity for this calculation:')
      WRITE(90,4233) (MOLNAME(I),I=1,NOSMOL)
4233  FORMAT(18(2X,A4))
C
        write(90,*)
        write(90,*)
     &  'Compute partial pressure sums of molecules with known OS:'

        write(90,4202) nos
4202    format('There are',i4,
     &    ' OS molecules with known GGchem partial pressure:')
        WRITE(90,4203) (molnames(idos(j)),j=1,nos)
        WRITE(90,4208) (idos(j),j=1,nos)
4203    format(9A8)
4208    format(9i8)
C
C
        write(90,4205) 
4205    format('There are 54 OS molecules in total that we have '
     &      ,'line lists for:')
        WRITE(90,4206) (osnames(j),j=1,54)
4206    format( 9(a4,4x) )

C
C k should be the depth layer if we should be able to make the test of
C summation in each layer, which would require a listing of the
C pp_j(k) for all j=1,54 OS-molecules; for now we set k=ndp here
C to avoid it being 0 (i.e. out of dimension). ppallmol(k,..) has
C already at this state been filled in by k=1,ntau calls to ggchem and
C corresponding saving in the 2D array.

      open(unit=707,file='pp.dat')
      if (krome_on.eq.1) then !avoid overriding non-eq results with final ggchem call
       read(707,*) (dummy_ppallat(ntau,m),m=1,22)
       read(707,*) (dummy_ppallmol(ntau,m),m=1,543)
      else
       read(707,*) (ppallat(ntau,m),m=1,22)
       read(707,*) (ppallmol(ntau,m),m=1,543)
      endif
      read(707,708) ((km,atnames(m)),m=1,22)   
      read(707,*)
      read(707,709) (molnames(m),m=1,543)

708         format(i4,18x,a2)
709           format(10a8)
      close(707)

C The positive and negative ions:
        natplus = 0
        do 3200 j=1,543
        ipos = 0
        ipos = index(molnames(j),'+')
        if (ipos.eq.0) go to 3200
        natplus = natplus + 1
        idplus(natplus) = j
3200    continue

        write(90,3201) natplus
3201    format('There are',i4,' positive ions:')
        WRITE(90,3253) (molnames(idplus(j)),j=1,natplus)
        WRITE(90,4209) (idplus(j),j=1,natplus)
4209    format(15i8)

        do 3220 i=1,jtau
        ppp(i) = 0.
        do 3221 j=1,natplus
        ppp(i) = ppp(i)+ppallmol(i,idplus(j))
3221    continue
3220    continue

        natminus = 0
        do 3203 j=1,543
        ipos = 0
        ipos = index(molnames(j),'-')
        if (ipos.eq.0) go to 3203
        natminus = natminus + 1
        idminus(natminus) = j
3203    continue

        write(90,3204) natminus
3204    format(/'There are',i4,' negative ions:')
        WRITE(90,3253) (molnames(idminus(j)),j=1,natminus)
        WRITE(90,4209) (idminus(j),j=1,natminus)

        do 3226 i=1,jtau
        ppm(i) = 0.
        do 3223 j=1,natminus
        ppm(i) = ppm(i)+ppallmol(i,idminus(j))
3223    continue
3226    continue


C The neutral molecules:
        neutralmol = 0
        do 3205 j=1,543
        ipos = 0
        ipos = index(molnames(j),'-')
        if (ipos.ne.0) go to 3205
        ipos = index(molnames(j),'+')
        if (ipos.ne.0) go to 3205
        neutralmol = neutralmol + 1
        idmol(neutralmol) = j
3205    continue

        write(90,3206) neutralmol
3206    format(/'There are',i4,' neutral molecules:')
        WRITE(90,3207) (molnames(idmol(j)),j=1,neutralmol)
        WRITE(90,3210) (idmol(j),j=1,neutralmol)
3207    FORMAT(2x,15a8)
3210    FORMAT(2x,15i8)

        do 3227 i=1,jtau
        ppnmol(i) = 0.
        do 3228 j=1,neutralmol
        ppnmol(i) = ppnmol(i)+ppallmol(i,idmol(j))
3228    continue
3227    continue

C The neutral atoms:
        WRITE(90,3202) (atnames(j),j=1,22)
3202    format(/'There are 22 neutral atoms:', 22(1x,a2,1x))

        do 3145 jm=1,2
        jmin=(jm-1)*15 + 1
        jmax = jmin+14
        if(jmin.gt.22) go to 3146
        if(jmax.ge.22) jmax=22
        WRITE(90,2233)
2233    FORMAT(//' P A R T I A L  P R E S S U R E S ',
     &    ' of all the 22 neutral atoms from GGchem')
        WRITE(90,3148) (atnames(j),j=jmin,jmax)
        WRITE(90,3211) (j,j=jmin,jmax)
        DO I=1,JTAU
                WRITE(90,'(I2,15F8.3)')I,
     *          (log10(max(1.e-99,(ppallat(i,j)))),J=jmin,jmax)
        end do

3145    continue
3146    continue

        WRITE(7,2234)
2234    FORMAT(//' P A R T I A L  P R E S S U R E S   ',
     &  ' of 15 selected neutral atoms of the 22 calculated in GGchem')
      
        WRITE(7,3148) (atnames(kpratoms(j)),j=1,15)
        DO I=1,JTAU
                 WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,ppallat(i,kpratoms(j)))),J=1,15)
        END DO
3148    FORMAT(2x,15(3x,a2,3x))
3211    FORMAT(2x,15(1x,i6,1x))

        WRITE(7,2235)
2235    FORMAT(//' P A R T I A L  P R E S S U R E S   of ',
     &  '13 selected ions of the 82(+)+52(-)=134 calculated in GGchem')
        WRITE(7,2236) (molnames(kions(j)),j=1,13)
        DO I=1,JTAU
              patp=   log10(max(1.e-99,ppp(i)))
              patm =  log10(max(1.e-99,ppm(i)))
              ppelk = log10(max(1.e-99,ppel(i)))
              pelbal = log10(max(1.e-99,(ppel(i)+ppm(i))))
           WRITE(7,'(I3,15F8.3)')I,
     *     (log10(max(1.e-99,ppallmol(i,kions(j)))),J=1,13),patp,pelbal
        END DO
2236    FORMAT(5x,12(a8),a6,' ions+  pe+ions-')


        WRITE(7,2239)
2239    FORMAT(//' P A R T I A L  P R E S S U R E S   of 15 iselected ',
     &  'neutral molecules (k1mol) of the 409 calculated in GGchem')
        WRITE(7,2240) (molnames(k1mol(j)),j=1,15)
        DO I=1,JTAU
                 WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,ppallmol(i,k1mol(j)))),J=1,15)
        END DO
2240    FORMAT(5x,14a8,a6)

        WRITE(7,2241)
2241    FORMAT(//' P A R T I A L  P R E S S U R E S   of 15 selected ',
     &  'neutral molecules (k2mol) of the 409 calculated in GGchem')
        WRITE(7,2240) (molnames(k2mol(j)),j=1,15)
        DO I=1,JTAU
                 WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,ppallmol(i,k2mol(j)))),J=1,15)
        END DO

        WRITE(7,2242)
2242    FORMAT(//' P A R T I A L  P R E S S U R E S   of 15 iselected ',
     &  'neutral molecules (k3mol) of the 409 calculated in GGchem')
        WRITE(7,2240) (molnames(k3mol(j)),j=1,15)
        DO I=1,JTAU
                 WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,ppallmol(i,k3mol(j)))),J=1,15)
        END DO

           do 3229 i=1,jtau
           ppnat(i) = 0.
           do 3229 j=1,22
           ppnat(i) = ppnat(i)+ppallat(i,j)
3229       continue

        WRITE(7,2243)
2243    FORMAT(//' P A R T I A L  P R E S S U R E S  -- overview of all'
     &  ,' the 565 neutral and ionized species from GGchem'/
     &  ,10x,'22 neutral atoms (atnames), 409 neutral molecules, '
     &  ,'82 positive ions, 52 negative ions (molnames):')
 
        do 2245 jm=1,2
        jmin=(jm-1)*15 + 1
        jmax = jmin+14
        if(jmin.gt.22) go to 3150
        if(jmax.ge.22) jmax=22
        WRITE(7,3148) (atnames(j),j=jmin,jmax)
                WRITE(7,'(2x,15(2x,i3,3x))')((j),J=jmin,jmax)
        DO I=1,JTAU,jtau-1
                WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,(ppallat(i,j)))),J=jmin,jmax)
        end do
2245    continue
3150    continue

        do 3143 jm=1,37
        jmin=(jm-1)*15 + 1
        jmax = jmin+14
        if(jmin.gt.543) go to 3144
        if(jmax.ge.543) jmax=543
        WRITE(90,2133)
        WRITE(90,3147) TEFF,GLOG,jm,jmin,jmax
        WRITE(90,3153) (molnames(j),j=jmin,jmax)
                WRITE(90,'(2x,15(2x,i3,3x))')((j),J=jmin,jmax)
        WRITE(7,3153) (molnames(j),j=jmin,jmax)
                WRITE(7,'(2x,15(2x,i3,3x))')((j),J=jmin,jmax)
        DO I=1,JTAU
                WRITE(90,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,(ppallmol(i,j)))),J=jmin,jmax)
                if(i.eq.1 .or. i.eq.jtau)
     *          WRITE(7,'(I3,15F8.3)')I,
     *          (log10(max(1.e-99,(ppallmol(i,j)))),J=jmin,jmax)
        end do

3143    continue
3144    continue
2153    FORMAT(' K ',15a8)
3153    FORMAT(5x,15a8)
3253    FORMAT(15a8)
3147    FORMAT(/' TEFF=',F6.0,' LOG G=',F5.1,6X,'jm,jmin,jmax=',3i4/)
        WRITE(7,2247)
2247    FORMAT(//' P A R T I A L  P R E S S U R E S  -- '
     &      ,'illustrative sums')
         write(7,3238)
3238     format(/'positive ions   negative ions  neutral atoms '
     &  ,' neutral molecules  sum{el,at,mol}  electron pres '
     &  ,'  static el.pres' 
     &  ,/'  ppH2  ppCH4  ',
     &   'ppNH3  ppH2O  ppH  ppHe ')

         WRITE(7,2248)
2248     format(/' I   pat+    pat-    pnat    pnmol    psum    pel  ',
     &    'pel-stat   ppH2   ppCH4   ppNH3   ppH2O',
     &    '    ppH    ppHe')

           do 3231 i=1,jtau
           ppamol(i) = 0.
           do 3230 j=1,543
           ppamol(i) = ppamol(i)+ppallmol(i,j)
3230       continue
          patp=   log10(max(1.e-99,ppp(i)))
          patm =  log10(max(1.e-99,ppm(i)))
          pnmolk = log10(max(1.e-99,ppnmol(i)))
          pamolk = log10(max(1.e-99,ppamol(i)))
          pek = log10(max(1.e-99,pe(i)))
          ppelk = log10(max(1.e-99,ppel(i)))
          pelbal = log10(max(1.e-99,(ppel(i)+ppm(i))))
          ppnatk = log10(max(1.e-99,ppnat(i)))
          pgk = log10(max(1.e-99,pg(i)))
          psumk = log10(max(1.e-99,(ppel(i)+ppamol(i)+ppnat(i))))
          pp24k = log10(max(1.e-99,pp24(i)))
          pp54k = log10(max(1.e-99,pp54(i)))
          pp24mk = log10(max(1.e-99,pp24m(i)))
          pp54mk = log10(max(1.e-99,pp54m(i)))
          pph2k = log10(max(1.e-99,ppallmol(i,1)))
          ppch4k = log10(max(1.e-99,ppallmol(i,185)))
          ppnh3k = log10(max(1.e-99,ppallmol(i,354)))
          pph2ok = log10(max(1.e-99,ppallmol(i,350)))
          pphk = log10(max(1.e-99,ppallat(i,1)))
          pphek = log10(max(1.e-99,ppallat(i,2)))
          WRITE(7,'(I2,15F8.3)')I,patp,patm,ppnatk,pnmolk,
     &       psumk,ppelk,pek,pph2k
     &      ,ppch4k,ppnh3k
     &      ,pph2ok,pphk,pphek
3231       continue


      END IF

      write(7,*)
      write(7,*)
      close(90)
      GO TO 2109

!------------ END OF  WRITING PP FOR JUMP = 4 -----------



2102  CONTINUE     ! go here if jump>2 (i.e. not old marcs chem.eq. and not JANAF fits).
      IF (JUMP.GT.4) GO TO 2103   !here one can put other chem.eq. possibilities
!       ERC: was GT3 before, but changed to GT4 to have Pgas printed in
!       SUB routine init_ggchem
C output from GEM:

C

      GO TO 2109
2103  CONTINUE
2109  CONTINUE     ! go here when finished the partial pressure writing



1905  FORMAT(I3,18(2X,A4))
      write(7,*)' '
      write(7,*)' '
      WRITE(7,1904) NOSMOL
1904  FORMAT
     *(I3,' MOLECULES(/+ATOMS) HAVE BEEN CONSIDERED IN THE OPACITY')


      WRITE(7,1900)
1900  format(' A B S O R P T I O N   C O E F F I C I E N T S  [CM/MOL]')
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,1906) (MOLNAME(I),I=1,NOSMOL)
1906  FORMAT(' K',52(2X,A4))
      DO 2074 k=1,JTAU
      do 2076 i=1,nosmol
2076  SUMOP(I,K)= max(1.d-99, SUMOP(I,K) )
      write(7,2078) k,( log10(SUMOP(I,K)),I=1,NOSMOL )
2074  continue
2078  format(i3,52f6.2)

      write(7,*)' '
      write(7,*)' '
      WRITE(7,1901)
1901  format(' I N T E G R A T E D   O P A C I T Y   [CM/G*]')
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6

      WRITE(7,1906) (MOLNAME(I),I=1,NOSMOL)
      DO 2072 k=1,JTAU
      do 2077 i=1,nosmol
2077  SUMKAP(I,K)= max(1.d-99, SUMKAP(I,K) )
      write(7,2073) k,( log10(SUMKAP(I,K)),I=1,NOSMOL )
2072  continue
2073  format(i3,52f6.1)
C End of molecular partial pressures


      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
C
C
CUGJ STORE ABSORBTIONCOEFFICIENT, SCATTERING AND LAMBDA FOR SPECTRUM 
CUGJ CALCULATIONS
C
      WRITE(7,*) ' '
      WRITE(7,*) ' '
      WRITE(7,*) 'LAMBDA FOR ABS-. AND SCAT.COEF. = '
      WRITE(7,*) NLP
      WRITE(7,530) (XLR(KLAM),KLAM=1,NLP)
530   FORMAT(10F12.3)
531   FORMAT(I4,E12.3)
532   FORMAT(1P10E12.5)
      WRITE(7,*) ' '
      WRITE(7,*) 'ABSORBTION COEF. AND SCATTERING COEF. I EACH TAU:'
      WRITE(7,*) ' '
      DO 81 KT=1,JTAU
      READ(IARCH) (ABSKTR(KLAM),KLAM=1,NLP)
      READ(IARCH) (SPRTR(KLAM),KLAM=1,NLP)
      WRITE(7,531) KT,TAU(KT)
      WRITE(7,532) (ABSKTR(KLAM),KLAM=1,NLP)
      WRITE(7,532) (SPRTR(KLAM),KLAM=1,NLP)
81    CONTINUE
      WRITE(7,*) ' '
!      do k=1,ntau
!      vert_scale_height=BOLTZK*T(k)/(emu(k)*XMH*log10(grav)) 
!      write(*,*) V(k),palfa,BOLTZK,T(k),XMH,emu(k)
!      write(*,*) vert_scale_height,palfa*vert_scale_height
!      enddo
C
C first, identify the depth points where logtau_ross is -4, -3, ..., 2
C for printing of opacities
       kdp(1) = 1
       kn = 2
      DO 233 K=1,JTAU
       taul = log10(tau(k))
      if (taul.gt.xdp(kn)) then
       k1 = max(k-1,1)
       kdp(kn) = k
       if ( abs(taul-xdp(kn)) .gt. abs( xdp(kn)-log10(tau(k1)) )  )
     *     kdp(kn) = k1
       if (kn.eq.8) go to 234
       kn = kn + 1
      end if
233   continue
234   continue
       kdp(8) = jtau

      HCK=143922240.

      DO 22 KTAU=1,JTAU
      DO 20 IE=1,NEL
      NJP=NJ(IE)
      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),abmarcs(IE,1),
     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
   20 CONTINUE
      DO 21 KLAM=1,NLP
      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE

C write only for selected layers:
      do 236 kd=1,8
      if (ktau.eq.kdp(kd)) go to 238
236   continue
      go to 22
238   continue

      IF(KTAU.LE.1) WRITE(7,218)
      IF(KTAU.GT.1) WRITE(7,219)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,220) TAUI
      WRITE(7,221) TI,PEI,kr,log10(tau(kr))
      WRITE(7,222)
      PESUM=0.
      DO 32 IE=1,NEL
      PEP(IE)=0.
      NJP=NJ(IE)
      IF(NJP.LT.2) GOTO 32
      DO 33 JJ=2,NJP
   33 PEP(IE)=PEP(IE)+abmarcs(IE,1)*ANJON(IE,JJ)*(JJ-1)
   32 PESUM=PESUM+PEP(IE)
      DO 19 IE=1,NEL
      NJP=NJ(IE)
      PEP(IE)=PEP(IE)/PESUM
      WRITE(7,223)IEL(IE),abmarcs(IE,1),PEP(IE),(ANJON(IE,JJ),JJ=1,NJP)

   19 CONTINUE
      WRITE(7,2250)
      WRITE(7,225) (ABNAME(KP),KP=1,NPROV)

      HCKT=HCK/TI
      DO 18 KLAM=1,NLP
      ABKLA=ABSKA(KLAM)
      STIM=1.-EXP(-HCKT/XLR(KLAM))
C      if (lin_cia.ne.1) 
C     *   abkla = abkla+(prov(nprova,klam)+prov(nprova-1,klam))*stim
C      WRITE(iwrop,2261)XLR(KLAM)/1.d4
C     *    ,log10( max(1.0d-99,ABSKA(KLAM)) )
C     *    ,log10( max(1.0d-99,SPRIDA(KLAM)) )
C     *    ,(log10( max(1.0d-99,PROV(J,KLAM)) ),J=1,NPROVA-2)
C     *    ,(log10( max(1.0d-99,PROV(J,KLAM)) ),J=NPROVA+1,NPROVA+NPROVS)
      DO 180 J=1,NPROVA-2
  180 PROV(J,KLAM)=PROV(J,KLAM)/ABKLA*STIM
      if (lin_cia.ne.1) 
     *   abkla = abkla+(prov(nprova,klam)+prov(nprova-1,klam))*stim
      DO 1801 J=NPROVA-1,NPROVA
 1801 PROV(J,KLAM)=PROV(J,KLAM)/ABKLA*STIM
      DO 181 J=1,NPROVS
  181 PROV(NPROVA+J,KLAM)=PROV(NPROVA+J,KLAM)/SPRIDA(KLAM)
      WRITE(7,226) XLR(KLAM),ABSKA(KLAM),SPRIDA(KLAM),
     &             (PROV(J,KLAM),J=1,NPROV)
   18 CONTINUE
   22 CONTINUE

C - skip output of colours:
      NOCOL=1
      IF(NOCOL.EQ.1) GO TO 1
C
C
      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C        CONVERT TO 'PHYSICAL' FLUXES
      DO 24 J=1,NLB
   24 FLUXME(J)=3.14159*FLUXME(J)
      DO 25 J=1,NLB
      JNORM=J
      IF(XLB(J).GT.5000.) GOTO 26
   25 CONTINUE
   26 STMAGN=-2.5*log10(FLUXME(JNORM))
C
C        WRITE(7,FLUXES)
C
      WRITE(7,521)
      N=0
      DO 52 I=1,NLB,4
      WAVFLX(N+1)=(XLB(I)+XLB(I+1)+XLB(I+2)+XLB(I+3))/4.
      WAVFLX(N+2)=100.*(A*(FLUXME(I)+FLUXME(I+3))+B*(FLUXME(I+1)
     1+FLUXME(I+2)))
      N=N+2
      IF(N.LT.10) GO TO 52
      WRITE(7,520) WAVFLX
      N=0
   52 CONTINUE
      IF(N.NE.0) WRITE(7,520)(WAVFLX(I),I=1,N)
  520 FORMAT(1X,5(F12.0,E13.5))
  521 FORMAT('1  BELL''S FLUX              (CENTER WAVE & FLUX)'//)
C
C
      WRITE(7,219)
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,227)
      LLB=NLB
   60 IF(LLB.LT.NLB) WRITE(7,219)
C      IF(LLB.LT.NLB) WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,
C     &                            IDRAB4,IDRAB5,IDRAB6
      IF(LLB.LT.NLB) WRITE(7,227)
      WRITE(7,228)
      IRAD = 0
      IF(LLB.LT.NLB) GOTO 65
      I=0
      I2=50
      I3=100
      I4=150
      GOTO 650
   65 CONTINUE
      I = I + 149
      I2 = I2 + 149
      I3 = I3 + 149
      I4 = I4 + 149
  650 CONTINUE
      ILI=2
      IF(LLB.GE.51) ILI=3
      IF(LLB.GE.101)ILI=4
      IF(LLB.GE.151)ILI=5
   66 I4=I4+1
   67 I3=I3+1
   68 I2=I2+1
   69 I=I+1
      IRAD = IRAD + 1
      IF(IRAD.GT.50) GOTO 70
      GOTO(70,61,62,63,64),ILI
   64 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3),
     &             XLB(I4),FLUXME(I4),FLUMAG(I4)
      IF(I4.GE.NLB) ILI=4
      GOTO 66
   63 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3)
      IF(I3.GE.NLB) ILI=3
      GOTO 67
   62 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2)
      IF(I2.GE.NLB) ILI=2
      GOTO 68
   61 WRITE(7,229) XLB(I),FLUXME(I),FLUMAG(I)
      IF(I.GE.NLB) GOTO 70
      GOTO 69
   70 LLB = LLB - 200
      IF(LLB.GT.0) GOTO 60
C
C        COMPUTE U, B, V, R, I AND COLOURS
      SUMCOL=0.
      SUMNOR=0.
      WRITE(7,270)
      DO 40 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3000..OR.XXLB.GT.4200.) GOTO 40
      INDEXx=INT(XXLB/100.)-29
      VIKT=W(I)*UW(INDEXx)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   40 CONTINUE
      UFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 41 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3500..OR.XXLB.GT.5600.) GOTO 41
      INDEXx=INT(XXLB/100.)-34
      VIKT=W(I)*BW(INDEXx)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   41 CONTINUE
      BFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 42 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.4700..OR.XXLB.GT.7200.) GOTO 42
      INDEXx=INT(XXLB/100.)-46
      VIKT=W(I)*VW(INDEXx)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   42 CONTINUE
      VFLUX=SUMCOL/SUMNOR
      DO 43 I=1,NLB
      JNORM=I
      IF(XLB(I).GE.7000.) GOTO 44
   43 CONTINUE
   44 RFLUX=FLUXME(JNORM)
      XINTCALL = 9000.0D+0
      CALL TINT(NLB,XLB,FLUXME,XINTCALL,XIFLUX)
      UMAG=-2.5*log10(UFLUX)-STMAGN
      BMAG=-2.5*log10(BFLUX)-STMAGN
      VMAG=-2.5*log10(VFLUX)-STMAGN
      RMAG=-2.5*log10(RFLUX)-STMAGN
      XIMAG=-2.5*log10(XIFLUX)-STMAGN
      UB=UMAG-BMAG
      BV=BMAG-VMAG
      UV=UMAG-VMAG
      RI=RMAG-XIMAG
      VR=VMAG-RMAG
      VI=VMAG-XIMAG
C      WRITE(7,300) TEFF,GLOG,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
      WRITE(7,2698)
      WRITE(7,271) UB,BV,UV
      WRITE(7,272) XLB(JNORM),RI,VR,VI
      WRITE(7,273) UMAG,BMAG,VMAG,RMAG,XIMAG


      do k=1,ntau
      pall(k) = pg(k)+pe(k)+prad(k)+pturb(k)
       if(k.ge.2 .and. k.le.ntau-1) then
         dz=(z(k+1)-z(k-1))/2.
         rho=ro(k)                        !(ro(k+1)+ro(k-1))/2.
         pall1 = pg(k+1)+pe(k+1)+prad(k+1)+pturb(k+1)
         pall2 = pg(k-1)+pe(k-1)+prad(k-1)+pturb(k-1)
         dp=(pall1-pall2)/2.
       else if (k.eq.1) then
         dz = z(2)-z(1)
         rho = (ro(2)+ro(1))/2.
         dp = pg(2)+pe(2)+prad(2)+pturb(2)
     &    -(pg(1)+pe(1)+prad(1)+pturb(1))
       else if (k.eq.ntau) then
         dz = z(ntau)-z(ntau-1)
         rho = (ro(ntau)+ro(ntau-1))/2.
         dp = pall(ntau)-pall(ntau-1)
       end if
      
      end do

   1  CONTINUE !read next archiv-file and write it in a formatted file.


5191  FORMAT(I3,f7.2,1P5E11.3,2x,4e11.3)
5190  FORMAT(//' H Y D R O S T A T I C    E Q U I L I B R I U M')
5192  FORMAT('  K log(tau)    z(k)     pall(k)    ro(k)     ',
     & 'dp/dz       ro*g      pall = pg +   pe     +  prad    +  pturb')

******* FORMAT SATS: *******************

  200 FORMAT(//' M O D E L  P A R A M E T E R S  E T C .'//)
  282 FORMAT(' * THE FOLLOWING MODEL WAS COMPUTED BY  S C M A R C S  ',
     & 'with Teff,log(g) = ',f6.0,f5.2)
  283 FORMAT(1X,81('*'))
  284 FORMAT(1X,81('*')////)
  201 FORMAT(' EFFECTIVE TEMPERATURE  ',F12.0,'  KELVIN'/' TOTAL ',
     +'FLUX (= SIGMA*TEFF**4)',1PE11.3,'  ERGS/S/CM**2'/' ACCELERATION',
     +' OF GRAVITY',0PF12.1,'  CM/S**2 (i.e., log(g) = ',f5.2,')'
     +/'Z/Zo (i.e. Oxygen/Oxygen_o) =',1PE8.1,'  C/O =',0PF6.2
     +//' CONVECTION PARAMETERS'/
     +' PALFA (L/HP)=',F5.2,',  PNY (NY)=',F5.2,',  PY (Y)=',F6.3)
 2017 FORMAT(' EFFECTIVE TEMPERATURE  ',F12.0,'  KELVIN'/' TOTAL ',
     +'FLUX (= SIGMA*TEFF**4)',1PE11.3,'  ERGS/S/CM**2'/' ACCELERATION',
     +' OF GRAVITY',0PF12.5,'  CM/S**2 (i.e., log(g) = ',f5.2,')'
     +/' Z/Zo (i.e. Oxygen/Oxygen_o) =',1PE8.1,'  C/O =',E8.1
     +//' CONVECTION PARAMETERS'/
     +' PALFA (L/HP)=',0PF5.2,',  PNY (NY)=',0PF5.2,',  PY (Y)=',0PF6.3)
 2011 FORMAT(' Convection was, however, excluded in the uppermost'
     +    'NOCONV=',I4,' layers')
 2010 FORMAT(' the Mihalas parameter was MIHAL=',I4,' layers')
 2012 FORMAT(' XMAX, TAUM, facply, moltsuji = ',1pe8.1,0p2f8.1,i3)
 2013 FORMAT(' Molecular equilibrium was treated by Tsujis routine')
 2014 FORMAT(' Molecular equilibrium was treated by the Marcs routin')
 2015 FORMAT(3x,A60)
 2018 FORMAT(' Molecular equilibrium was treated by JF fit to JANAF 99')
20181 FORMAT(' Molecular equilibrium was treated by GGChem')
20182 FORMAT(' Electronpressure for continuum opacities was',
     &          ' from GGChem')
 2019 FORMAT(' Molecular equilibrium was by Gibbs Energy Minimisation')
  256 FORMAT(' TURBULENCE PRESSURE IS NEGLECTED (PBETA<=0.)')
  257 FORMAT(' TURBULENCE PRESSURE IS INCLUDED AND PBETA=',F5.2)
  250 FORMAT(' CONVECTION HAS BEEN INCLUDED IN THIS MODEL; ISTRAL=',I3)
  251 FORMAT(' CONVECTION HAS  N O T  BEEN INCLUDED IN THIS MODEL; '
     *  ,'ISTRAL=',I3)
  252 FORMAT(' LINE BLANKETING HAS  N O T  BEEN INCLUDED IN THIS MODEL;'
     * ,' ILINE=',I3)
 2529 FORMAT(' CIA described by Linskys formulas was included')
  253 FORMAT(' LINE BLANKETING HAS BEEN INCLUDED IN THIS MODEL;'
     * ,' ILINE=',I3)
 2531 FORMAT(' THE MASS WAS ASSUMED TO BE RELM =',F5.1)
 2539 FORMAT(' -- but this computation is plane parallel ')
 2541 FORMAT(' For',f5.0'<T<',f5.0,'K',f6.1,
     * '-',f6.1,' % of C2H2 was assumed to be grains, and'/,
     *  ' for',f5.0'<T<',f5.0,'K',f6.1,
     * '-',f6.1,' % was assumed to be molcular C2H2')
 2534 FORMAT(' THE OS SAMPLING HAVE BEEN PERFORMED IN',I3,' INTERVALS')
 2535 FORMAT(' FROM',F6.1,' CM-1 TO',F9.1,' CM-1')
 2536 FORMAT(' (..which is identical to from',F9.1,'  to',F6.1,' mu)')
 2537 FORMAT(' THE CHOSEN INTERVAL-BEGINNING AND STEPLENGTH IN CM-1:')
 2538 FORMAT(8F9.1)
  254 FORMAT(' THE STROEMGREN EQUATION HAS BEEN USED FOR THE UPPERMOST',
     *I4,' POINTS')
  255 FORMAT(////' THE FOLLOWING MODEL WAS OBTAINED AFTER',I3,
     *' ITERATION(S)')
  202 FORMAT(//'  LOG. ABUNDANCES USED IN MODEL CALCULATIONS'/2X,'H'
     *,5X,'HE',4X,'C',5X,'N',5X,'O',5X,'NE',4X,'NA',4X,'MG',4X,'AL',4X,
     *'SI',4X,'S',5X,'K',5X,'CA',4X,'CR',4X,'FE',4X,'NI',4X,'TI')
  203 FORMAT (20F6.2)
  204 FORMAT(//' C O R R E C T I O N S  I N  T H E  L A S T ',
     *' I T E R A T I O N')
  205 FORMAT(' K',' LOG(TauR)',2X,'TAUROSS',6X,'T',7X,'DELTA(T)',5X,
     *'FCONV',4X,'DELTA(FCONV)',2X,'K')
  206 FORMAT(I3,F7.2,1PE11.2,0P2F11.2,1PE12.3,E12.2,I5)
  207 FORMAT(//' M O D E L  A T M O S P H E R E   (CGS UNITS)')
  208 FORMAT(3X,' K',4X,'TAUROSS',3X,'TAU(5000)',2X,'GEOM. DEPTH',5X,
     *'T',8X,'PE',10X,'PG',9X,'PRAD',8X,'PTURB',5X,'KAPPAROSS',5X,'K')
  209 FORMAT(I5,1P2E12.4,E13.4,0PF9.2,1P5E12.4,I4)
 2061 FORMAT(//' S U M   O F   P A R T I A L    P R E S S U R E S')
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
 2063 FORMAT('  K  P_MOL, P_NEU_HCNO, P_ION_HCNO, P_NEU_HE'
     &,' P_ION_HE, P_NON_HHECNO, P6_ION PG_JON, P6 PHe')
 2071 FORMAT(//' M A S S  L O S S AND S P H E R E AND  L E V I T A T',
     *' I O N')
 2081 FORMAT(2X,' K',2X,'TauRoss',2X,'log(tau)',2X,'density,g/cm3',2X,
     *'Z,cm',2x,'T',7x,'arad',5X,'dpr/dz',4X,'V-loss,km/s',3x,'rossos'
     * ,4x,'p-rossos',6x,'cross')
 2091 FORMAT(I5,1PE11.2,0PF6.2,1P2E11.2,0PF7.0,1P2E11.2,0PF8.4,
     * 1P3E12.3)
  210 FORMAT(//' T H E R M O D Y N A M I C A L  Q U A N T I T I E S ',
     *' AND  C O N V E C T I O N  (CGS UNITS)')
  211 FORMAT(2X,' K',3X,'TAUROSS',3X,'DENSITY',5X,'MU',9X,'CP',10X,
     *'CV',7X,'ADGRAD',9X,'Q',7X,'SOUND VEL.',2X,'CONV. VEL.',1X,
     *'FCONV/F',3X,'K')
  212 FORMAT(I4,1P2E11.3,0PF8.3,6(1PE12.3),0PF9.5,I3)
  213 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (CGS UNITS)')
 2131 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (Tsuji)')
 2132 FORMAT(//' L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S   (GGCHEM)')
 2133 FORMAT(//' SELECTED L O G A R I T H M I C  M O L E C U L A R ',
     *' P A R T I A L  P R E S S U R E S  from GGCHEM')
  214 FORMAT(' K',3X,'P(H)',2X,'P(H-)',2X,'P(H2)',2X,'P(H2+)',1X,
     *'P(H2O)',1X,'P(OH)',2X,'P(CH)',2X,'P(CO)',2X,'P(CN)',2X,'P(C2)',2X
     *,'P(N2)',2X,'P(O2)',2X,'P(NO)',2X,'P(NH)',1X,'P(TIO)',4X,'K')
  215 FORMAT(I3,15f8.3,2X,I3)
 2141 FORMAT(' K',2X,'TiO',5X,'TiO2',4X,'TiC2',4X,'TiO-KE',5X,
     *    'P(He)',3x,'P(HI)',3x,'P(Hx2)',3X,'P6',6X,'PG',5X,' K')
 2151 FORMAT(I3,9F8.3,2X,I3)
  216 FORMAT(' K',4X,'C2H2',4X,'HCN',5X,'C2H',5X,'HS',6X,'SIH',5X
     * ,'C3H',5X,'C3',6X,'CS',6X,'SIC',5X,'SIC2',4X,'NS',6X,'SIN',5X
     * ,'SIO',5X,'SO',6X,'S2',6X,'SIS',6X,'K')

 2161 FORMAT(' K',3x,' VO ',3x,'ZrO',4x,'CaH',4x,'LaO',5x,'CH4',4x
     * ,'CH2',4x,'CH3',4x,'NH3',4x,'CO2',4x,'Si2C',3x,'SiO2',3x
     * ,'H2S',4x,'CS2',4x,'CaOH',3x,'NO2',4x,'AlH',5x,'K ')

  217 FORMAT(I3,16F8.3,2X,I3)
 2171 FORMAT(' K',3x'AlO ',3x,'SiH4',3x,'SO2',4x,'YO',5x,'C5H',4x
     * ,'C4H',4x,'C4',5x,'C5 ',4x,'C6H ',3x,'C4H6',3x,'C6H2',3x,'C6H4'
     * ,3x,'C10H7',2x,'HC5N',3x,'HC7N',3x,'HC9N',3x,'K')

 2181 FORMAT(' K',3x,'HC11N',2x,'C4H4S',2x,'C4H4O',2x,'C5H5N',2x
     * ,'C6H4',3x,'C6H5O',2x,'C6H6O',2x,'C10H8',2x,'C10H16',1x,
     * 'C14H10a','C14H10p','C18H12t','C18H12b','C22H14c','C24H12'
     * ,1x,'C106H56','K')

  218 FORMAT(//' I O N I Z A T I O N  C O N D I T I O N S  AND ',
     *' A B S O R P T I O N  M E C H A N I S M S')
  219 FORMAT(1H1)
  220 FORMAT(/' ****************'/' * TAU=',F8.4,' *'/
     *' ****************')
  221 FORMAT(/' T=',F7.0,'  PE=',1PE9.2,' lag:',I3,' lg10_tau(kr)='
     #    ,0pf7.2)
  222 FORMAT(' ELEMENT  ABUNDANCE'
     *                      ,' ELCONT IONIZATION FRACTIONS',17X,
     *'PARTITION FUNCTIONS'/30X,'I',7X,'II',6X,'III',5X,'IV',12X,
     *'I',9X,'II',8X,'III',7X,'IV')
  223 FORMAT(6X,A2,1PE12.3,0PF6.3,4(F8.4))
  224 FORMAT(1H+,62X,4(1PE10.2))
 2250 FORMAT(////16X,'F R A C T I O N S  O F  C O N T I N U O U S  A B',
     *' S O R P T I O N  A N D  S C A T T E R I N G')
  225 FORMAT(' WAVELENGTH   ABS       SCAT    ',16(1X,A5))
2251  FORMAT('  W-mu log->ABS SCAT',3x,16A7)
2252  FORMAT(' k-depth-point in atmosphere:',i3,' ~log(tau_ross)=',f8.2)
2253  FORMAT(' Units: log_10 of values of opacity in cm^2/g*' )
  226 FORMAT(F11.0,1p2E10.3,0p16F6.3)
2261  FORMAT(F6.2,18f7.2)
2871  format(' wn_mu cont.op line-op. log_partp*akapmol for the',
     *        i4,' molecules:',15(2x,a4,1x))

288   format(///' O P A C I T Y   OF   O S - M O L E C U L E S',
     *  '  in log-cm^2/g-*')
287   format('  log(tau)  cont.op  line-op. partp*akapmol for:'
     *        ,15(2x,a4,1x))
286   format(i3,f7.2,15f7.2)
2861   format(i3,2f7.2,1p11e10.3,/10e10.3)
2862   format(1p13e10.3)
2863   format(f7.2,16f7.2)
285   format(' os-wn',i5,':',f8.1,' (=',f8.4,' mu)')
C285   format(' os-wn',i5,':',f8.1)
  227 FORMAT(//' F L U X E S (PHYSICAL FLUXES IN ERGS/S/CM**2/ANGSTROM)'
     */)
  228 FORMAT(4(6X,'LAMBDA',5X,'FLUX',7X,'MAGN'))
  229 FORMAT(4(F12.0,1PE11.2,0PF9.3))
  260 FORMAT(////' MOLECULES OTHER THAN H2 AND H2+ HAVE NOT BEEN',
     *' CONSIDERED IN THE TEMPERATURE - ELECTRON PRESSURE - PRESSURE',
     *' BALANCE')
 2698 FORMAT(1H0/' (UBV TRANSMISSION FUNCTIONS AFTER MATTHEWS +',
     *' SANDAGE.  AIR MASS = 1)')
  270 FORMAT(' G I A N T  L I N E  C O L O U R S')
  271 FORMAT(8X,'U - B =',F8.3//9X,'B - V =',F8.3//9X,'U - V =',F8.3
     *////)
  272 FORMAT(' TENTATIVE CONTINUUM COLOURS ',
     *' (R AT',F6.0,'A AND I AT 9000A)'//9X,
     *'R - I =',F8.3,//9X,'V - R =',F8.3//9X,'V - I =',F8.3)
  273 FORMAT(////' U =',F10.3,'  B =',F10.3,'  V =',F10.3,'  R =',F10.3,
     *'  I =',F10.3)
  300 FORMAT(/' TEFF=',F6.0,' LOG G=',F5.1,2X,6A4/)
      RETURN
         END
C
      SUBROUTINE MAINB
      implicit real*8 (a-h,o-z)
C     
      include 'parameter.inc'
C     
      real*8 krome_tmax,krome_photo_scale
      COMMON /UTPUT/IREAD,IWRIT
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /CSTYR/MIHAL,NOCONV 
      COMMON /CXMAX/XMAX /CTAUM/TAUM
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX                          
      COMMON /CANGLE/XMY(6),XMY2(6),H(6),MMY
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CLINE4/ILINE
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE 
      COMMON /CPOLY/FACPLY,MOLTSUJI
      COMMON /CMETBL/METBL
      COMMON /MASSE/RELM
      COMMON /CMETPE/ PPEL(NDP), METPE
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLIN/lin_cia
      COMMON /CORRECT/TDIFF,TCONV,KORT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /STATEC/DUM1(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /Cspec/spec(nwl,3),ispec
      common /cosexp/ lops,nops
      common /dpeset/ dpein,dtin
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom
      common /ch4/ nch4
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      common /noneq_time/ dt_start,dt_max,dt_inc,krome_tmax
      common /noneq_output/ krome_output,krome_debug,krome_return
      common /starspec/ stellar_spectrum(nwreal),index_wlambda
      DATA TSUN,GSUN,RSUN/5800.,4.44,7E10/
      character(len=200) :: spectrum_file


C INITIATIONS
      HPLNCK=6.62554E-27
      BOLTZK=1.38054E-16
      CLIGHT=2.997925E10
      ECHARG=4.80298E-10
      RYDBRG=1.097373E5
      STEFAN=5.675E-5
      CLIGHT=2.99793E10
      PI=3.14159265
      PI4C=PI*4./CLIGHT
C
C LOGICAL UNITS
      IREAD=5
      IWRIT=7
C
C TEMPERATURE, GRAVITATION
      READ(5,62) TEFF,G,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,RELM
      
      !print*, teff
      
C
C KONV, MIHAL
CUGJ      ILINE=0
      ISTRAL=0
      READ(5,631) NOCONV,XMAX,TAUM,FACPLY,METBL
     
      IF(METBL.GT.1) METBL=1
      IF(TAUM.EQ.0.) TAUM=50.
      IF(XMAX.EQ.0.) XMAX=1.E10
      READ(5,63) ILINE,AMLOSS,MOLTSUJI,LIN_CIA,ISPEC
      READ(5,51) MIHAL,KONSG,KORT,TDIFF,TCONV
      
      read(5,635) lops,nops,dpein,dtin,metpe

C CONVECTION PARAMETERS
      READ(5,50) PALFA,PBETA,PNY,PY,VFIX
      
C
C MY POINTS  (warning: there are other variables with the same names xmy,h...)
C                      
C
      READ(5,512) MMY,NCORE,KDIFF
      read(5, 1234) irrin,steff,rstar,semimajor,f_irrad,
     > input_star_spec,spectrum_scale

      if (irrin == 1) then
      print*, "Irradiation is turned on"
       if (input_star_spec==1) then
        print*, 'Using input spectrum for irradiation'

        spectrum_file = trim("./data/stellar_spectrum.dat")
        stellar_spectrum = 0.0d0
        print*, "Stellar specturm scaling set to", spectrum_scale    
        open(unit=20, file=spectrum_file, status='old', action='read')
        do j=1,nwreal
         read(20,*) stellar_spectrum(j) 
        end do
        close(20)

       else
      print*, 'Using blackbody spectrum for irradiation'
      print*, "Stellar effective temperature of ", steff, " Kelvin "
      print*, "Stellar radius of  ", rstar, " solar radii"
      print*, "Planet at ",semimajor, "AU from star" 
      
      endif
      
      else

      print*, "Irradiation is turned off" 
      end if   
      read(5, 1235) irrinp, tbottom, reflect, nch4
      if (irrinp == 1) then
      print*, "Modelling a rocky planet with a surface temperature of ", 
     * tbottom, " and a surface albedo of ", reflect
      end if
      read(5, 1236) krome_on, krome_photo_on, krome_photo_scale
      read(5, 1237) dt_start, dt_max, krome_tmax, dt_inc
      read(5, 1238) krome_output,krome_debug,krome_return
      if (krome_on.eq.1) then
      print*, "Non-equilibrium Chemistry is turned on"
        if (krome_photo_on.eq.1) then
         print*, "Photorates included"
         print*, "with photo chemistry scaled by", krome_photo_scale
        endif
      else
      print*, "Non-equilibrium Chemistry is turned off"      
      endif

      !print*, irrinp
      FLUX=STEFAN*TEFF**4/PI
      GRAV=10.**G
      RELRAD=0.5*(log10(RELM)+GSUN-G)
      RELLUM=4.*log10(TEFF/TSUN)+2.*RELRAD
      RELATM=-3.+log10(TSUN/TEFF/RELM)+0.5*RELLUM
      BOLMAG=4.72-2.5*RELLUM
      RADIUS=RSUN*10.**RELRAD
      DO 100 K=1,NDP
100   RR(K)=RADIUS

      TAURAT=10.**(0.1*KDIFF-0.01)
C
      ACALL = 0.0D+0
      BCALL = 1.0D+0
      CALL GAUSI(MMY,ACALL,BCALL,H,XMY)
      DO 110 I=1,MMY
110   XMY2(I)=XMY(I)**2
C
C TAU SCALE
      CALL TAUSCA
C
C NTAU is number of depth points in input-model. JTAU is number of depth
C points demanded in the input-file (input tauscale). In most model
C computations, these two numbers will be identical, but if more (or
C less) depth points are wanted to be computed compared to the input model,
C the values from the ntau layers of the input model will be scaled
C (inter/extra-polated)to the jtau layers of the new model to estimate a
C new model as a starting model for the iterations.
C
      
      IF (NTAU.NE.JTAU) CALL SCALEMOD
      
      DO 200 K=1,NTAU
      TAULN(K)=log(TAU(K))
200   CONTINUE
C
     
C      WRITE(7,59)NOCONV
C      WRITE(7,57)MIHAL
C      WRITE(7,571)KONSG,KORT,TDIFF
      
      RETURN
50    FORMAT(5(7X,F8.0))
51    FORMAT(3(7X,I3,5X),7X,F8.0,7X,F8.0)
512   FORMAT(3(7X,I3,5X))
1234  format(1(7X,I4,4X), 4(7X, F8.0), 1(7X,I4,4X), 1(7X,E8.1))
1235  format(1(8X,I4,3X), 2(7X, F8.0), 1(7X, I2))
52    FORMAT(20X,'LOG G  =',F10.2,10X,'LOG (ATM/R) =',F5.2,10X,
     & 'LOG (R/RSUN)=',F5.2)
53    FORMAT(/20X,'PALFA  =',F10.2)
54    FORMAT(20X,'PBETA  =',F10.2)
55    FORMAT(20X,'PNY    =',F10.0)
56    FORMAT(20X,'PY     =',F10.3)
57    FORMAT(20X,'MIHAL  =',I5)
571   FORMAT(20X,'KONSG  =',I5,'  KORT  =',I5,'  TDIFF  =',F8.0)
58    FORMAT(20X,'MYPNTS =',I10,'(for spherical only)')
59    FORMAT(/20X,'NOCONV =',I10)
60    FORMAT('0* MODEL PARAMETERS')
61    FORMAT(/20X,'TEFF   =',F10.0,10X,'LOG (L/LSUN)=',F5.2,10X,
     & 'BOLOM. MAGN.=',F6.2)
62    FORMAT(2(7X,F8.0),7X,6A4,6X,F8.0)
621   FORMAT(2(7X,F8.1),7X,6A4,6X,F8.1)
631   FORMAT(7X,I3,12X,E8.1,7X,f8.0,7x,f8.0,7x,I3)
632   FORMAT(7X,I3,12X,1pE8.1,7X,0pf8.1,7x,f8.1,7x,I3)
63    FORMAT(7X,I3,12X,E8.1,7X,I3,2(12X,I3))
633   FORMAT(7X,I3,12X,1pE8.1,3(7X,I3))
634   FORMAT(2(7X,F8.1),2(7X,f8.4))
635   FORMAT(7X,I3,12x,i3,12x,f7.5,8x,f7.5,8x,i3)
64    FORMAT(/20X,'NCORE  =',I10/20X,'KDIFF  =',I10)
65    FORMAT(20X,'M/MSUN =',F10.1)
66    FORMAT(20X,'XMAX   =',1PE10.2)
67    FORMAT(20X,'TAUM   =',F10.2)
68    FORMAT(20X,'METBL  =',I10)
69    FORMAT(20X,'FACPLY =',F10.3)
1236  FORMAT(7X,I3,12X,I3,12X,E8.1)
1237  FORMAT(7X,E8.1,7X,E8.1,7X,E8.1,7X,F8.2)
1238  FORMAT(7X,I3,12X,I3,12X,I3)
      END
C
      SUBROUTINE MATINV(A,N)
      implicit real*8 (a-h,o-z)
C
C 'MATINV' IS A STANDARD ROUTINE FOR MATRIX INVERSION (USED IN THE
C MIHALAS CODE). 'MATINV' ASSUMES THAT A(J,J) IS NONZERO AND MAKES
C NO PIVOTING. THIS IS SOMETIMES ADVANTAGEOUS FOR THE NUMERICAL
C ACCURACY (IN THE MARCS CODE FOR INSTANCE).
C
      include 'parameter.inc'
      DIMENSION A(NDP,NDP)
C
      IF(N.EQ.1)GOTO 25
      DO 5 I=2,N
      IM1=I-1
      DO 2 J=1,IM1
      JM1=J-1
      if (A(J,J) == 0.0) A(J,J) = 1.0e-30
      DIV=A(J,J)
      !if (div == 0.0) print*, "J at div == 0 " , J
      SUM=.0
      IF(JM1.LT.1)GOTO 2
      DO 1 L=1,JM1
    1 SUM=SUM+A(I,L)*A(L,J)
    2 A(I,J)=(A(I,J)-SUM)/DIV
      DO 4 J=I,N
      SUM=.0
      DO 3 L=1,IM1
    3 SUM=SUM+A(I,L)*A(L,J)
      A(I,J)=A(I,J)-SUM
    4 CONTINUE
    5 CONTINUE
      DO 13  II=2,N
      I=N+2-II
      IM1=I-1
      IF(IM1.LT.1)GOTO 13
      DO 12 JJ=1,IM1
      J=I-JJ
      JP1=J+1
      SUM=.0
      IF(JP1.GT.IM1)GOTO 12
      DO 11 K=JP1,IM1
   11 SUM=SUM+A(I,K)*A(K,J)
   12 A(I,J)=-A(I,J)-SUM
   13 CONTINUE
      DO 17 II=1,N
      I=N+1-II
      DIV=A(I,I)
      IP1=I+1
      IF(IP1.GT.N)GOTO 17
      DO 16 JJ=IP1,N
      J=N+IP1-JJ
      SUM=.0
      DO 15 K=IP1,J
   15 SUM=SUM+A(I,K)*A(K,J)
      A(I,J)=-SUM/DIV
   16 CONTINUE
   17 A(I,I)=1./A(I,I)
      DO 24 I=1,N
      DO 23 J=1,N
      K0=MAX0(I,J)
      IF(K0.EQ.J)GOTO 22
      SUM=.0
   20 DO 21 K=K0,N
   21 SUM=SUM+A(I,K)*A(K,J)
      GOTO 23
   22 SUM=A(I,K0)
      IF(K0.EQ.N)GOTO 23
      K0=K0+1
      GOTO 20
   23 A(I,J)=SUM
   24 CONTINUE
      RETURN
   25 A(1,1)=1./A(1,1)
      RETURN
      END
C
      SUBROUTINE MODJON(JONTYP,IOUTS)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'

C
C 'JONTYP' DETERMINES THE QUALITY OF THE IONIZATION EQUILIBRIUM COMPUTED
C IN JON.
C
C JONTYP.LE.0  A NUMBER OF ELEMENTS AND IONIZATION STAGES ARE DISREGARDE
C        EQ.1  PARTITION FUNCTIONS ARE TAKEN CONSTANT OR TEMPERATURE DEP
C              ONLY, AS SPECIFIED IN INJON-DATA. ELEMENTS AND IONIZATION
C              STAGES AS SPECIFIED IN INJON-DATA.
C       .EQ.2  FULL PARTITION FUNCTIONS, ASYMTOTIC PART ACCORDING TO
C              BASHEK ET AL.
C       .GE.3  FULL PARTITION FUNCTIONS, ASYMTOTIC PART ACCORDING TO
C              FISCHEL ET AL. ALL POSSIBLE ELEMENTS AND IONIZATION STAGE
C              ARE TAKEN INTO ACCOUNT.
C
C EACH INCREASE OF JONTYP CORRESPONDS TO AN ADDING OF THE QUALITY CITED
C ON TOP OF PREVIOS ONES.
C
C IOUTS .EQ.0  NO PRINT OUT OF NEW JON-PARAMETERS
C       .GT.0  PRINTOUT OF NEW JON-PARAMETERS
C
C ENTRY POINT 'MODMOL' ADDED 76.03.23  *NORD*
C
C MOLTYP.LE.0  GIVES FULL HANDLING OF MOLECULES
C        GT.0  GIVES H-MOLECULES ONLY
C

      COMMON /CI3/DUM(885),IDUM(215),IFISH
      COMMON /CI4/TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      COMMON /CI6/TP,IQFIX(16,5),NQTEMP
      DIMENSION IELS(16),IONS(16,5),IQFS(16,5)
      DATA JONOLD/1/,IMAX/16/,JMAX/5/
      INTEGER MOLH, JUMP
C
C CHECK IF READY
      IF(JONTYP.EQ.JONOLD) RETURN
      IF(JONTYP.GT.3.AND.JONOLD.EQ.3) RETURN
      IF(JONTYP.LT.0.AND.JONOLD.EQ.0) RETURN
      TP=0.
      IF(JONTYP.LT.JONOLD) GO TO 60
C
C ZERO TO ONE
      IF(JONOLD.GT.0) GO TO 20
      DO 10 I=1,IMAX
      IELEM(I)=IELS(I)
      DO 10 J=1,JMAX
      ION(I,J)=IONS(I,J)
10    CONTINUE
      JONOLD=1
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C ONE TO TWO
20    IF(JONOLD.GT.1) GO TO 30
      DO 21 I=1,IMAX
      DO 21 J=1,JMAX
      IQFS(I,J)=IQFIX(I,J)
      IQFIX(I,J)=2
21    CONTINUE
      IFISH=0
      JONOLD=2
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C TWO TO THREE
30    IF(JONOLD.GT.2) GO TO 90
      IFISH=1
      DO 31 I=1,IMAX
      IELS(I)=IELEM(I)
      IELEM(I)=1
      DO 31 J=1,JMAX
      IONS(I,J)=ION(I,J)
      ION(I,J)=1
31    CONTINUE
      JONOLD=3
      GO TO 90
C
C THREE TO TWO
60    IF(JONOLD.LT.3) GO TO 70
      IFISH=0
      DO 61 I=1,IMAX
      IELEM(I)=IELS(I)
      DO 61 J=1,JMAX
      ION(I,J)=IONS(I,J)
61    CONTINUE
      JONOLD=2
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C TWO TO ONE
70    IF(JONOLD.LT.2) GO TO 80
      DO 71 I=1,IMAX
      DO 71 J=1,JMAX
      IQFIX(I,J)=IQFS(I,J)
71    CONTINUE
      JONOLD=1
      IF(JONOLD.EQ.JONTYP) GO TO 90
C
C ONE TO ZERO
80    IF(JONOLD.LT.1) GO TO 90
      DO 81 I=1,IMAX
      IELS(I)=IELEM(I)
      DO 81 J=1,JMAX
      IONS(I,J)=ION(I,J)
81    CONTINUE
      ION(2,2)=0
      ION(3,3)=0
      IELEM(4)=0
      IELEM(6)=0
      IELEM(7)=0
      ION(8,3)=0
      ION(9,3)=0
      ION(9,4)=0
      ION(10,3)=0
      IELEM(11)=0
      IELEM(12)=0
      IELEM(13)=0
      IELEM(14)=0
      IELEM(16)=0
      ION(15,3)=0
      JONOLD=0
C
C RETURN
90    CONTINUE
      IF(IOUTS.GT.0) CALL INJON(2)
      RETURN
C
C ENTRY MODMOL
      ENTRY MODMOL(MOLTYP,IOUTS)
      MOLH=0
      !IF(MOLTYP.LE.0) MOLH=0
      !IF(MOLTYP.GT.0) MOLH=1
      IF(IOUTS.GT.0) CALL INJON(2)
      WRITE(7,101) MOLTYP
101   FORMAT('0MOLTYP=',I1)
      RETURN
      END
C
      SUBROUTINE MOL(T,PE,G2,GC,GN,GO,ABUC,ABUO,ABUN,
     &               XIH,XKHM,XIHM,XNEN,F1,F2,F3,F4,F5)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
C        THIS ROUTINE COMPUTES DISSOCIATION EQUILIBRIA FOR THE MOLECULES H2,H2+
C        H2O,OH,CH,CO,CN,C2,O2,N2,NH AND NO WITH H,H-,H+,C,C+,O,O+,N,N+ CON-
C        SIDERED, USING A NEWTON-RAPHSON SCHEME. SOME FEATURES COME FROM THE
C        MONSTER AND FROM MIHALAS, METH. COMP. PHYS. 7,1.
C
C        G2=N(HII)/N(HI), GC=N(CII)/N(CI) ETC.
C        ABUC= THE NUMBER OF CARBON NUCLEI PER HYDROGEN NUCLEUS, ABUO AND ABUN
C        ARE THE CORRESPONDING VALUES FOR OXYGEN AND NITROGEN.
C        XIH = THE IONIZATION ENERGY OF HYDROGEN
C        XIHM= THE 'DISSOCIATION ENERGY' OF H-
C        XKHM= THE 'DISSOCIATION CONSTANT' OF H-
C        XNEN= THE NUMBER OF ELECTRONS PER UNIT VOLUME FROM ELEMENTS OTHER THAN
C        HYDROGEN, CARBON, OXYGEN AND NITROGEN.
C
C        THE SUBSCRIPT IN AKD(I),PK(I)  HAS THE FOLLOWING MEANING
C        I=1 H-, 2 H2, 3 H2+, 4 H2O, 5 OH, 6 CH, 7 CO, 8 CN, 9 C2, 10 N2,
C         11 O2, 12 NO, 13 NH, 14 C2H2, 15 HCN, 16 C2H, 17 -, 18 HS
C         19 SIH, 20 C3H, 21 C3, 22 CS, 23 SIC, 24 SIC2, 25 NS
C         26 SIN, 27 SIO, 28 SO, 29 S2, 30 SIS, 31 TiO, 32 TiO2, 33 TiC2
C        WHEREAS AKA, AK1 - AK4  HAVE  I=I-1  (CF. LOOP 100)
C
C        NMOL IS THE NUMBER OF MOLECULES CONSIDERED
C
C        THIS ROUTINE CALLS MOLMAT AND AINV3
C        THE DATA FOR COMPUTING THE DISSOCIATION CONSTANTS ARE FROM TSUJI
C        (ASTRON. ASTROPHYS.  23,411 (1973))
C        WITH CORRECTIONS OF NH(TO D0=3.46 EV),
C        BG 831114.
C
C     841210 KE:  CHANGED DISS.EN. FOR CS  WITH -0.5 EV ( IN AK1 )
C     950625 UGJ: CHANGED DISS.EN. FOR TIO TO 6.89 EV (Costes 95) BY CHANGING
C                 AK1(30) FROM 8.5956 (Tsuji_73_D0=7.26eV) to 8.2256
C     950625 UGJ: CHANGED DISS.EN. FOR CN TO 7.77 EV (Costes 94) BY CHANGING
C                 AK1(6) FROM 8.2793 (Tsuji_73_D0=7.9eV) to 8.1493
C
C        THE ROUTINE GIVES FH, FC, FO, FN, FE. FH=P(HI)/PH, FC=P(CI)/PH ETC.,
C        WHERE PH=NH*KT (NH IS THE NUMBER OF HYDROGEN NUCLEI PER CM3).
C
      DOUBLE PRECISION AKA(32),AK1(32),AK2(32),AK3(32),AK4(32),AKD(33),
     & PK(33),TH,WKH,WKC,WKN,WKO,FH,FC,FN,FO,FE,FS,FK,FT,FF(8),
     & CAM,ROOT,EMAX,DIFF,A(8,8),F(8),D(8),S,R,PELLE(4),PALLE(4),
     & R1,R2,R3,R4,R5,ROT(4),ZIM(4),X02AAF
C     & CAM,ROOT,EMAX,DIFF,A(7,7),F(7),D(7),S,R,PELLE(4),PALLE(4),
C ARGUMENTS TO STATEMENT FUNCTION 880916, apollo
      DOUBLE PRECISION X,DX
      EQUIVALENCE (FF(1),FH),(FF(2),FC),(FF(3),FO),(FF(4),FN),
     & (FF(5),FE),(FF(6),FS),(FF(7),FK),(FF(8),FT)
      dimension qr(22),tqr(22),ptio(ndp)
     &          ,xnti1(ndp),xnti2(ndp)
      DIMENSION B1(5),B2(5),DIS(10)
      COMMON /CMOL1/EH,FFE,FFH,FHE,FFC,FCE,FFN,FNE,FFO,FOE,FFK,FKE,
     &              FFS,FSE,FFT,FTE
      COMMON /CMOL2/PPK(33),NMOL
      COMMON /CI5/ abmarcs(18,ndp),ANJON(18,5),DUMT(99)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /ggchembool/ iggcall
      data tqr/1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500.,
     &         5000., 5500., 6000., 6500., 7000., 8000., 9000., 10000.,
     &        12000., 14000., 16000., 18000., 20000., 50000./
      data qr/1.708, 1.911, 2.035, 2.097, 2.112, 2.091, 2.042, 1.972,
     &        1.887, 1.789, 1.684, 1.574, 1.462, 1.241, 1.032, 0.884,
     &        0.541, 0.330, 0.194, 0.111, 0.062, 0.001/
      DATA B1/2.6757,1.4772,0.60602,0.12427,0.00975/,
     *B2/2.9216,2.0036,1.7231,0.82685,0.15253/,
     *DIS/9.50,4.38,3.47,11.11,7.90,6.12,9.76,5.12,6.51,3.21/
      DATA AKA/12.739D0,11.206998D0,25.420D0,12.371D0,12.135D0,13.820D0,
     *12.805D0,12.804D0,13.590D0,13.228D0,12.831D0,12.033D0,
     +38.184D0,25.635D0,25.063D0,
     +0.0D0,12.019D0,11.852D0,40.791D0,25.230D0,13.436D0,12.327D0,
     +25.623D0,12.543D0,12.399D0,13.413D0,12.929D0,12.960D0,
     +13.182D0,13.398D0,27.901D0,27.018D0/
      DATA AK1/5.1172D0,2.7942767D0,10.522D0,5.0578D0,4.0760D0,11.795D0,
     *8.1493D0,6.5178D0,10.585D0,5.5181D0,7.1964D0,4.0935D0,
     +17.365D0,13.833D0,12.291D0,
     +0.0D0,4.2922D0,3.7418D0,21.762D0,14.445D0,8.0574D0,5.0419D0,
     +13.085D0,5.9563D0,5.4876D0,8.8710D0,6.0100D0,5.0952D0,
     +7.1147D0,8.2256D0,14.031D0,13.534D0/
      DATA AK2/1.2572D-1,-7.9196803D-2,1.6939D-1,1.3822D-1,1.2768D-1,
     *1.7217D-1,6.4162D-2,9.7719D-2,2.2067D-1,6.9935D-2,1.7349D-1,
     +1.3629D-1,2.1512D-2,1.3827D-1,-1.9036D-2,
     +0.0000D+0,1.4913D-1,1.5999D-1,9.3377D-1,1.2547D-1,1.8754D-1,
     +1.3941D-1,-5.5227D-2,2.0901D-1,9.5301D-2,1.5042D-1,1.6253D-1,
     +1.8027D-1,1.9300D-1,0.40873D0,0.42156D0,0.45875D0/
      DATA AK3/1.4149D-2,-2.4790744D-2,1.8368D-2,1.6547D-2,1.5473D-2,
     *2.2888D-2,7.3627D-3,1.2739D-2,2.9997D-2,8.1511D-3,2.3065D-2,
     +1.6643D-2,8.8961D-5,1.8122D-2,-4.4498D-3,
     +0.0000D+0,1.8666D-2,2.0629D-2,1.3863D-1,1.7390D-2,2.5507D-2,
     +1.9363D-2,-9.3363D-3,2.8986D-2,1.3369D-2,1.9581D-2,2.1665D-2,
     +2.4324D-2,2.5826D-2,5.7937D-2,6.1271D-2,6.6158D-2/
      DATA AK4/6.3021D-4,0.D0,8.1730D-4,7.7224D-4,7.2661D-4,1.1349D-3,
     *3.4666D-4,6.2603D-4,1.4993D-3,3.7970D-4,1.1380D-3,7.8691D-4,
     +-2.8720D-5,9.1645D-4,-2.3073D-4,
     +0.0000D+0,8.9438D-4,9.9897D-4,7.4549D-3,8.8394D-4,1.2735D-3,
     +9.6202D-4,-4.9876D-4,1.4621D-3,6.9396D-4,9.4828D-4,1.0676D-3,
     +1.2049D-3,1.2648D-3,2.9287D-3,3.1476D-3,3.3834D-3/
C
C STATEMENT FUNCTION FOR CORRECTIONS LIMITED TO FACTOR TWO
      ASUM(X,DX)=DMIN1(2.D0*X,DMAX1(0.5D0*X,X+DX))
C COMPUTATION OF DISSOCIATION CONSTANTS K(AB) (AKD) AND PE/K(AB) (PK).
      TMEM=T

      TH=5040.D0/T
      AKD(1)=XKHM
      PK(1)=PE/XKHM
      PELOG=LOG10(PE)
      DO 100 J=2,33
        M=J-1
        
        AKD(J)=AKA(M)-(AK1(M)-(AK2(M)-(AK3(M)-AK4(M)*TH)*TH)*TH)*TH
        PK(J)=PELOG-AKD(J)
100   CONTINUE

      PK(4)=PELOG+PK(4)
      PK(14)=2.0*PELOG+PK(14)
      PK(15)=PELOG+PK(15)
      PK(16)=PELOG+PK(16)
      PK(20)=2.0*PELOG+PK(20)
      PK(21)=PELOG+PK(21)
      PK(24)=PELOG+PK(24)
      PK(32)=PELOG+PK(32)
      PK(33)=PELOG+PK(33)
      DO 101 J=2,33
101   PK(J)=10.D0**PK(J)

        thta=5040./t
*
* find the ratio  u(Ti II)/u(Ti I)  by interpolating in array qr
*
      
        j=1
        do 662 jj=1,22-1
          if(t.le.tqr(jj)) goto 663
          j=jj
662     continue
663     uratio=qr(j) + (qr(j+1)-qr(j))/(tqr(j+1)-tqr(j))*(t-tqr(j))
*
* q1 is the ratio  n(Ti II)/n(Ti I)
*
       
        q1=.6667*t**2.5*uratio*exp(-6.82*11605./t)/pe
      
*
*
* q2 is the ratio  n(TiO)/n(Ti I)
*
         tiokp = 10.**akd(31)

*
* xnti1 is n(Ti I)    (i.e. per cm**3)
*

         ptio(kl) =  0.
         xnti1(kl) = 0.
         xnti2(kl) = 0.
C
      ANJON(17,2) = 1./(1.+1./q1)
      ANJON(17,1) = 1./(1.+q1)

C COMPUTATION OF STARTING VALUES FOR FH,FC ETC.
      XNENSK=XNEN-abmarcs(11,kl)*ANJON(11,2)-abmarcs(10,kl)*ANJON(10,2)
      GS=ANJON(11,2)/ANJON(11,1)
      GK=ANJON(10,2)/ANJON(10,1)
      GT=ANJON(17,2)/ANJON(17,1)
      AS=abmarcs(11,kl)/abmarcs(1,kl)
      AK=abmarcs(10,kl)/abmarcs(1,kl)
      AT=abmarcs(17,kl)/abmarcs(1,kl)
C
C START WITH VALUES OF FH,FC,FN,FO,FS,FK,FE FROM LAST ITTERATION
      IF (MOLOLD.EQ.1) THEN
         FH = FOLD(KL,1)
         FC = FOLD(KL,2)
         FN = FOLD(KL,3)
         FO = FOLD(KL,4)
         FS = FOLD(KL,5)
         FK = FOLD(KL,6)
         FE = FOLD(KL,7)

         FT = FOLD(KL,8)

         GO TO 159
      END IF
C
C ...OR ESTIMATE THESE VALUES...
      WKH=1.D0+G2+PK(1)
      WKC=1.D0+GC
      WKO=1.D0+GO
      WKN=1.D0+GN
C
C FH, FROM H=HI+HII+H2+H-
      FE=XNENSK

      CAM=FE*WKH/(2.D0*PK(2))
      ROOT=DSQRT(DABS(CAM*CAM+FE/PK(2)))
      FH=-CAM+ROOT
      FE=FH*(G2-PK(1))+XNENSK

      IF (FE.GT.0.D0) GO TO 110
      R=PK(2)/PK(1)/PK(1)
      CAM=(WKH*XNENSK/PK(1)-1.D0-2.D0*XNENSK*R)
      S=R-WKH/PK(1)
      CAM=CAM/2.D0/S
      ROOT=DSQRT(CAM*CAM-R*XNENSK*XNENSK/S)
      FE=-CAM-ROOT
      IF (FE.LT.0.) FE=-CAM+ROOT
      FH=(XNENSK-FE)/PK(1)
110   CONTINUE
C
C FN, FROM N=NI+NII+N2
      R1=PK(10)/FE
      R2=WKN
      R3=-abmarcs(4,kl)
      FN=(-R2+DSQRT(R2*R2-4.D0*R1*R3))/(2.D0*R1)
C
C FC, FROM C=CI+CII+CO+C2+C2H2+HCN+CN AND O=OI+OII+CO+H2O
      NDEG =3
      IF (abmarcs(5,kl).GE.abmarcs(3,kl)) GO TO 130
      R1=PK(7)/(FE*WKO)
      R2=2.D0*FH*FH*PK(14)/(FE*FE*FE)+2.D0*PK(9)/FE
      R3=WKC+FN*PK(8)/FE+FH*FN*PK(15)/(FE*FE)
      R4=WKO+FH*FH*PK(4)/(FE*FE)
C                                  R4 = WATER CORRECTION, 4-JAN-82/NORDLUND
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      SUMLOG=R1LOG+R2LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/2.

         R1=R1*10.**RADD
         R2=R2*10.**RADD
      END IF
      PELLE(1)=R1*R2
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      R3LOG=LOG10(R3)
      R4LOG=LOG10(R4)
      SUMLOG=R1LOG+R2LOG+R3LOG+R4LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/4.

         R1=R1*10.**RADD
         R2=R2*10.**RADD
         R3=R3*10.**RADD
         R4=R4*10.**RADD
      END IF
      PELLE(2)=R2*R4+R1*R3
      PELLE(3)=R3*R4+R1*(abmarcs(5,kl)-abmarcs(3,kl))
      PELLE(4)=-abmarcs(3,kl)*R4

120   FORMAT(1X,5D10.3)
      DO 121 IP=1,4
        PALLE(IP)=PELLE(IP)
121   CONTINUE
      NPOL=NDEG+1
      IFAIL=0
      CALL C02AEF(PELLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 122 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(3,kl)) GO TO 122
        FC=ROT(III)
        IFLAG=IFLAG+1
122   CONTINUE
      FO=abmarcs(5,kl)/(R4+FC*R1)
C.... PRINT 120,FC**3*PALLE(1),FC**2*PALLE(2),FC*PALLE(3),PALLE(4)
      IF (IFLAG.EQ.1) GO TO 135
      IF (IFLAG.EQ.0) PRINT 123
      IF (IFLAG.GE.2) PRINT 124
123   FORMAT('0MOL: NO ROOT FOUND FOR CARBON INITIAL RATIO FC')
124   FORMAT('0MOL: SEVERAL ROOTS FOR CARBON INITIAL RATIO FC')
      STOP ' stop at loop 122 in mol: root problem in mol '
C
C FO, FROM O=OI+OII+CO+H2O+SIO, C=CI+CII+CO+CN+HCN, AND SI=SI(I)+SI(II)+SIO
130   CONTINUE
      R1=PK(7)/FE
      R2=PK(27)/FE
      R3=WKC+FN*PK(8)/FE+FH*FN*PK(15)/(FE*FE)
      R4=1.D0+GK
      R5=WKO+FH*FH*PK(4)/(FE*FE)
      R1LOG=LOG10(R1)
      R2LOG=LOG10(R2)
      R3LOG=LOG10(R3)
      R4LOG=LOG10(R4)
      R5LOG=LOG10(R5)
      SUMLOG=R1LOG+R2LOG+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.

         R1=R1*10.**RADD
         R2=R2*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(1)=R1*R2*R5
      SUMLOG=MAX(R1LOG+R4LOG,R2LOG+R3LOG)+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.

         R1=R1*10.**RADD
         R4=R4*10.**RADD
         R2=R2*10.**RADD
         R3=R3*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(2)=(R1*R4+R2*R3)*R5+R1*R2*
     *    (abmarcs(3,kl)+abmarcs(10,kl)-abmarcs(5,kl))
      SUMLOG=R3LOG+R4LOG+R5LOG
      IF (SUMLOG.GE.34.) THEN
         RADD=(34.-SUMLOG)/3.
         
         R3=R3*10.**RADD
         R4=R4*10.**RADD
         R5=R5*10.**RADD
      END IF
      PALLE(3)=R5*R4*R3-(R1*R4+R2*R3)*abmarcs(5,kl)+R1*R4*abmarcs(3,kl)
     & +R2*R3*abmarcs(10,kl)
      PALLE(4)=-R3*R4*abmarcs(5,kl)
      DO 131 IP=1,4
        PELLE(IP)=PALLE(IP)
131   CONTINUE
C.... PRINT 120,PELLE
      NPOL=NDEG+1
      IFAIL=0
      CALL C02AEF(PELLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 132 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(5,kl)) GO TO 132
        FO=ROT(III)
        IFLAG=IFLAG+1
132   CONTINUE
      FC=abmarcs(3,kl)/(WKC+R1*FO+FN*PK(8)/FE+FN*FH*PK(15)/(FE*FE))
C.... PRINT 120,FO**3*PALLE(1),FO**2*PALLE(2),FO*PALLE(3),PALLE(4)
      IF (IFLAG.EQ.1) GO TO 135
      IF (IFLAG.EQ.0) PRINT 133
      IF (IFLAG.GE.2) PRINT 134
133   FORMAT('0MOL: NO ROOT FOUND FOR OXYGEN INITIAL RATIO FO')
134   FORMAT('0MOL: SEVERAL ROOTS FOR OXYGEN INITIAL RATIO FO')
      STOP ' stop at loop 132 in mol;  root problems in mol '
135   CONTINUE
C
C FK, FROM SI=SI(I)+SI(II)+SIS+SIH+SIC2+SIO, AND S=S(I)+S(II)+CS+SIS
      R1=1.0D0+GS+FC*PK(22)/FE
      R2=PK(30)/FE
      R3=1.D0+GK+FH*PK(19)/FE+FC*FC*PK(24)/FE
      R4=PK(27)/FE
      R5=1.0D0+GO+FC*PK(7)/FE
      ALFA=abmarcs(10,kl)-abmarcs(11,kl)-abmarcs(5,kl)
      PALLE(1)=R2*R3*R4
      PALLE(2)=R3*(R1*R4+R2*R5)-ALFA*R2*R4
      PALLE(3)=R1*R3*R5-ALFA*(R1*R4+R2*R5)-abmarcs(11,kl)*R1*R4-
     &  abmarcs(5,kl)*R2*R5
      PALLE(4)=-R1*R5*abmarcs(10,kl)
      NDEG=3
      IFAIL=0
      NPOL=NDEG+1
      CALL C02AEF(PALLE,NPOL,ROT,ZIM,X02AAF(DUM),IFAIL)
      IFLAG=0
      DO 140 III=1,3
        IF (ROT(III).LE.0.D0.OR.ROT(III).GT.abmarcs(10,kl)) GO TO 140
        IFLAG=IFLAG+1
        FK=ROT(III)
140   CONTINUE
C
C FS, FROM S=S(I)+S(II)+SIS+CS
      FS=abmarcs(11,kl)/(R1+FK*R2)
      IF (IFLAG.EQ.1) GO TO 152
      IF (IFLAG.LE.0) PRINT 150
      IF (IFLAG.GT.1) PRINT 151
150   FORMAT('0MOL: NO ROOT FOUND FOR SILICON INITIAL RATIO FK')
151   FORMAT('0MOL: SEVERAL ROOTS FOR SILICON INITIAL RATIO FK')
      STOP ' stop at loop 140;  root problems in mol '
152   CONTINUE
C
159   CONTINUE  
C  if use of old FH,FC,.. go directly here
C
C NEWTON-RAPHSON IMPROVEMNT USING ALL RELEVANT MOLECULES, ATOMS AND IONS.
C DIFF GIVES THE APPROX. ACCURACY TO WHICH FH ETC. HAVE TO CONVERGE
C BEFORE THE ITERATIONS ARE STOPPED.
C
      DIFF=1.D-3
       M=8


      DO 163 J=1,500
c        
        CALL MOLMAT(PK,G2,GC,GN,GO,GS,GK,GT,ABUC,ABUN,ABUO
     &             ,AS,AK,AT,FH,FC,FN,FO,FS,FK,FT,FE,XNENSK,F,A)


        CALL AINV3(A,M)

        EMAX=0.D0
        DO 162 L=1,M
          D(L)=0.D0
          DO 161 LL=1,M
            D(L)=D(L)-A(L,LL)*F(LL)
161       CONTINUE
          FF(L)=ASUM(FF(L),D(L))
          FF(L) = max(FF(L),1.d-99)
          if (D(L).le.1.d-99) go to 162
          EMAX=DMAX1(DABS(D(L)/FF(L)),EMAX)

162     CONTINUE
Ctemp UGJ 020308:   IF (EMAX.LT.DIFF) GO TO 170
163   CONTINUE
        GO TO 170
      write(7,*)' MOL: THE DESIRED ACCURACY WAS NOT ACHIEVED AFTER'
      WRITE(7,*)' 500 ITERATIONS. LAST VALUES AND DIFFERENCES WERE'
      WRITE(7,164) FF
      WRITE(7,164) D
164   FORMAT(1P8D10.3)
      STOP ' stop at molecular equilibrium iteration in mol '
C
C COMPUTATION OF THE INNER ENERGY. DEH2 AND DEH2P ARE THE SUM OF
C DISSOCIATION, ROTATION AND VIBRATION ENERGIES (IN EV PER MOLECULE) FOR
C H2 AND H2+. DIS(I) IS THE DISSOCIATION ENERGY FOR THE MOLECULE (I+3)
C IN THE LIST OF MOLECULES (VALUES ARE FROM TSUJI). FOR THESE MOLECULES
C THE ROTATION AND VIBRATION ENERGIES ARE NEGLECTED.
C
170   CONTINUE
c
C      write(7,168) FH,FC,FN,FO,FS,FK,FT,FE
C168   FORMAT(' After molmath call: FH,C,N,O,S,Si,Ti,E='/,1P8E9.2)
c      PRINT 160,FH,FC,FN,FO,FS,FK,FT,FE,J
c
   98 TETA=TH
      DEH2=(B1(1)-(B1(2)-(B1(3)-(B1(4)-B1(5)*TETA)*TETA)*TETA)*TETA)*
     *8.617E-5*T-4.476
      DEH2P=(B2(1)-(B2(2)-(B2(3)-(B2(4)-B2(5)*TETA)*TETA)*TETA)*TETA)*
     *8.617E-5*T-2.648
      FHE=FH/FE
      FCE=FC/FE
      FOE=FO/FE
      FNE=FN/FE
      FKE=FK/FE
      FSE=FS/FE
      FTE=FT/FE
      EH2=(-2.*XIH+DEH2)*FHE*FH*PK(2)
      EH2P=(DEH2P-XIH)*FHE*FH*G2*PK(3)
      EHM=-(XIHM+XIH)*FH*PK(1)
      EHJ=-XIH*FH
      EH2O=-(2.*XIH+DIS(1))*FHE*FH*FO*PK(4)
      EOH=-(XIH+DIS(2))*FOE*FH*PK(5)
      ECH=-(XIH+DIS(3))*FHE*FC*PK(6)
      ECO=-DIS(4)*FCE*FO*PK(7)
      ECN=-DIS(5)*FCE*FN*PK(8)
      EC2=-DIS(6)*FCE*FC*PK(9)
      EN2=-DIS(7)*FNE*FN*PK(10)
      EO2=-DIS(8)*FOE*FO*PK(11)
      ENO=-DIS(9)*FNE*FO*PK(12)
      ENH=-(DIS(10)+XIH)*FNE*FH*PK(13)
      EH=EH2+EH2P+EHM+EHJ+EH2O+EOH+ECH+ECO+ECN+EC2+EN2+EO2+ENO+ENH
C                            NOTE THAT ENERGIES INCLUDE ONLY MOLECULES 1-13
C
C PICK UP SINGLE PRECISION VALUES
      FFH=FH
      FFC=FC
      FFN=FN
      FFO=FO
      FFK=FK
      FFS=FS
      FFT=FT
      FFE=FE
c      PRINT*,' FFO,FFE,FOE ',FFO,FFE,FOE
      NMOL=33
      DO 180 I=1,NMOL
        PPK(I)=PK(I)
180   CONTINUE
      F1=FH
      F2=G2*FH
      F3=FH*PK(1)
      F4=FH*FHE*G2*PK(3)
      F5=FH*FHE*PK(2)
      T=TMEM
C
      FOLD(KL,1) = FH
      FOLD(KL,2) = FC
      FOLD(KL,3) = FN
      FOLD(KL,4) = FO
      FOLD(KL,5) = FS
      FOLD(KL,6) = FK
      FOLD(KL,7) = FE
      FOLD(KL,8) = FT


      RETURN
      END
C
      SUBROUTINE MOLEQ(T,PE,G2,XIH,XKHM,XIHM,XNENH,F1,F2,F3,
     &                    F4,F5,FE,FSUM,EH)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES DISSOCIATION EQUILIBRIA FOR THE MOLECULES
C        H2 AND H2+ WITH H+, H AND H- CONSIDERED. IT MAINLY FOLLOWS MIHA
C        METH. COMP. PHYS. 7, 1 (1967).
C
C        THE INNER ENERGY OF THE HYDROGEN GAS, EH, IS ALSO EVALUATED.
C
C        XIH=THE IONIZATION ENERGY OF HYDROGEN
C        XKHM=THE 'DISSOCIATION CONSTANT' OF H-
C        XIHM=THE 'DISSOCIATION ENERGY' OF H-.
C        XNENH=THE NUMBER OF ELECTRONS PER UNIT VOLUME FROM ELEMENTS OTH
C        HYDROGEN (Q IN MIHALAS'S ARTICLE)
C        G2,F1,F2 ETC. SEE REF.
C
C
      COMMON/UTPUT/IREAD,IWRIT
C
C        CALL MOLFYS FOR PHYSICAL DATA
      CALL MOLFYS(T,XKH2,XKH2P,DEH2,DEH2P)
C
C        CALCULATION OF THE EQUILIBRIUM
      G3=PE/XKHM
      G4=PE/XKH2P
      G5=PE/XKH2
      A=1.+G2+G3
      E=G2*G4/G5
      B=2.*(1.+E)
      C=G5
      D=G2-G3
      C1=C*B*B+A*D*B-E*A*A
      C2=2.*A*E-D*B+A*B*XNENH
      C3=-(E+B*XNENH)
      CAM=C2/(2.*C1)
      ROOT=SQRT(CAM*CAM-C3/C1)
      F1D=-CAM+ROOT
      IF(F1D.GT.1.0)F1D=-CAM-ROOT
      F5D=(1.0-A*F1D)/B
      F4D=E*F5D
      F3D=G3*F1D
      F2D=G2*F1D
      FED=F2D-F3D+F4D+XNENH
      FSUMD=F1D+F2D+F3D+F4D+F5D
      F1=F1D
      F2=F2D
      F3=F3D
      F4=F4D
      F5=F5D

      FE=FED

      FSUM=FSUMD
C
C        CALCULATION OF THE ENERGIES
      EH2=(-2.*XIH+DEH2)*F5
      EH2P=(-XIH+DEH2P)*F4
      EHM=-(XIHM+XIH)*F3
      EHJ=-XIH*F1
      EH=EHJ+EHM+EH2+EH2P
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE MOLFYS(T,XKH2,XKH2P,DEH2,DEH2P)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE GIVES DISSOCIATION CONSTANTS XKH2 (=N(H I)*N(H I)/(NH2))
C        AND XKH2P (=N(H I)*N(H II)/N(H2+)), EXPRESSED IN NUMBER PER CM3 AND
C        THE SUM OF DISSOCIATION, ROTATION AND VIBRATION ENERGIES, DEH2 AND
C        DEH2P FOR H2 AND H2+, RESPECTIVELY (EXPRESSED IN ERGS PER MOLECULE)
C        THE DATA ARE FROM VARDYA, M.N.R.A.S. 129, 205 (1965) AND EARLIER
C        REFERENCES. THE DISSOCIATION CONSTANT FOR H2 IS FROM TSUJI,
C        ASTRON. ASTROPHYS. 1973.
C
      DIMENSION A1(5),A2(4),B1(5),B2(5),TE(5)
      DATA A1/12.739,-5.1172,1.2572E-1,-1.4149E-2,6.3021E-4/,
     *A2/11.20699 ,-2.794276 ,-0.079196   ,0.024790   /,
     *B1/2.6757,-1.4772,0.60602,-0.12427,0.009750/,
     *B2/2.9216,-2.0036,1.7231,-0.82685,0.15253/

      TEX=5040./T
      TE(1)=1.
      DO1 K=1,4
    1 TE(K+1)=TE(K)*TEX
      XKH2=0.
      XKH2P=0.
      DEH2=0.
      DEH2P=0.
      DO2 K=1,4
      XKH2=A1(K)*TE(K)+XKH2
      XKH2P=A2(K)*TE(K)+XKH2P
      DEH2=B1(K)*TE(K)+DEH2
    2 DEH2P=B2(K)*TE(K)+DEH2P
      XKH2=A1(5)*TE(5)+XKH2
      DEH2=(B1(5)*TE(5)+DEH2)*8.617E-5*T-4.476
      DEH2P=(B2(5)*TE(5)+DEH2P)*8.617E-5*T-2.648
      XKH2=10.**XKH2
      XKH2P=10.**XKH2P
      RETURN
      END
C
      SUBROUTINE MOLMAT(PK,GGH,GGC,GGN,GGO,GGS,GGK,GGT,AAC,AAN,AAO,
     &               AAS,AAK,AAT,FH,FC,FN,FO,FS,FK,FT,FE,XYNEN,F,A)
      implicit real*8 (a-h,o-z)

      DIMENSION PK(33),F(8),A(8,8)


      GH=GGH
      GC=GGC
      GN=GGN
      GO=GGO
      GS=GGS
      GK=GGK
      GT=GGT
      XNEN=XYNEN
      AC=AAC
      AN=AAN
      AO=AAO
      AS=AAS
      AK=AAK
      AT=AAT
CCC
      FHE=FH/FE
      FCE=FC/FE
      FNE=FN/FE
      FOE=FO/FE
      FSE=FS/FE
      FKE=FK/FE
      FTE=FT/FE
      H=1.0+GH+PK(1)+FCE*PK(6)+FNE*PK(13)+FOE*PK(5)+FSE*PK(18)
     *  +FCE*FNE*PK(15)+FCE*FCE*(PK(16)+FCE*PK(20))+FKE*PK(19)
      HH=2.0*FHE*(PK(2)+GH*PK(3)+FCE*FCE*PK(14)+FOE*PK(4))
      C=1.0+GC+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+FSE*PK(22)+FKE*PK(23)+FHE*
     *  FNE*PK(15)
      CC=2.0*FCE*(PK(9)+FHE*FHE*PK(14)+FHE*PK(16)+FKE*PK(24)+FTE*PK(33))
      CCC=3.0*FCE*FCE*(PK(21)+FHE*PK(20))
      XN=1.0+GN+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+FSE*PK(25)+FKE*PK(26)+
     *   FHE*FCE*PK(15)
      XNN=2.0*FNE*PK(10)
      O=1.0+GO+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE*PK(28)+FKE*PK(27)+
     *  FTE*PK(31)+FHE*FHE*PK(4)
      OO=2.0*FOE*(PK(11)+FTE*PK(32))
      S=1.0+GS+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+FOE*PK(28)+FKE*PK(30)
      SS=2.0*FSE*PK(29)
      XK=1.0+GK+FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *   PK(30)+FCE*FCE*PK(24)
      T=1.0+GT+FOE*PK(31)+FOE*FOE*PK(32)+FCE*FCE*PK(33)
      F(1)=FH*(H+HH)-1.0
      F(2)=FC*(C+CC+CCC)-AC
      F(3)=FO*(O+OO)-AO
      F(4)=FN*(XN+XNN)-AN
      F(5)=FH*(GH-PK(1)+FHE*GH*PK(3))+FC*GC+FN*GN+FO*GO+FS*GS+FK*GK
     *     +XNEN-FE
      F(6)=FS*(S+SS)-AS
      F(7)=FK*XK-AK
      F(8)=FT*T-AT
      A(1,1)=H+2.0*HH
      A(1,2)=FHE*(PK(6)+FNE*PK(15)+2.0*FCE*PK(16)+3.0*FCE*FCE*PK(
     *       20)+4.0*FCE*FHE*PK(14))
      A(1,3)=FHE*(PK(5)+2.0*FHE*PK(4))
      A(1,4)=FHE*(PK(13)+FCE*PK(15))
      A(1,5)=-FHE*(2.0*FHE*(PK(2)+GH*PK(3))+FCE*PK(6)+FNE*PK(13)+
     *       FOE*PK(5)+FSE*PK(18)+FKE*PK(19))-2.0*FHE*(FCE*(FNE*
     *       PK(15)+FCE*PK(16))+2.0*FHE*FOE*PK(4))-3.0*FCE*FCE*
     *       FHE*(FCE*PK(20)+2.0*FHE*PK(14))
      A(1,6)=FHE*PK(18)
      A(1,7)=FHE*PK(19)
      A(1,8)=0.
      A(2,1)=FCE*(PK(6)+FNE*PK(15)+FCE*(4.0*FHE*PK(14)+2.0*PK(16)+
     *       3.0*FCE*PK(20)))
      A(2,2)=C+2.0*CC+3.0*CCC
      A(2,3)=FCE*PK(7)
      A(2,4)=FCE*(PK(8)+FHE*PK(15))
      A(2,5)=-FCE*(2.0*FCE*PK(9)+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+
     *       FSE*PK(22)+FKE*PK(23))-2.0*FCE*(FHE*FNE*PK(15)+FCE*(3.0
     *       *FCE*PK(21)+2.0*FHE*PK(16)+2.0*FKE*PK(24)))-3.0*FCE
     *       *FCE*FHE*(2.0*FHE*PK(14)+3.0*FCE*PK(20))
      A(2,6)=FCE*PK(22)
      A(2,7)=FCE*(PK(23)+2.0*FCE*PK(24))
      A(2,8)=4.*FCE*PK(33)
      A(3,1)=FOE*(PK(5)+FHE*PK(4)*2.0)
      A(3,2)=FOE*PK(7)
      A(3,3)=O+2.0*OO
      A(3,4)=FOE*PK(12)
      A(3,5)=-FOE*(2.0*FOE*PK(11)+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE
     *       *PK(28)+FKE*PK(27)+2.0*FHE*FHE*PK(4))
      A(3,6)=FOE*PK(28)
      A(3,7)=FOE*PK(27)
      A(3,8)=FOE*PK(31)+4.*FOE*PK(32)
      A(4,1)=FNE*(PK(13)+FCE*PK(15))
      A(4,2)=FNE*(PK(8)+FHE*PK(15))
      A(4,3)=FNE*PK(12)
      A(4,4)=XN+2.0*XNN
      A(4,5)=-FNE*(2.0*FNE*PK(10)+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+
     *       FSE*PK(25)+FKE*PK(26)+2.0*FHE*FCE*PK(15))
      A(4,6)=FNE*PK(25)
      A(4,7)=FNE*PK(26)
      A(4,8)=0.
      A(5,1)=GH*(1.0+2.0*FHE*PK(3))-PK(1)
      A(5,2)=GC
      A(5,3)=GO
      A(5,4)=GN
      A(5,5)=-GH*FHE*FHE*PK(3)-1.0
      A(5,6)=GS
      A(5,7)=GK
      A(5,8)=GT
      A(6,1)=FSE*PK(18)
      A(6,2)=FSE*PK(22)
      A(6,3)=FSE*PK(28)
      A(6,4)=FSE*PK(25)
      A(6,5)=-FSE*(2.0*FSE*PK(29)+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+
     *       FOE*PK(28)+FKE*PK(30))
      A(6,6)=S+2.0*SS
      A(6,7)=FSE*PK(30)
      A(6,8)=0.
      A(7,1)=FKE*PK(19)
      A(7,2)=FKE*(PK(23)+2.0*FCE*PK(24))
      A(7,3)=FKE*PK(27)
      A(7,4)=FKE*PK(26)
      A(7,5)=-FKE*(FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *       PK(30)+2.0*FCE*FCE*PK(24))
      A(7,6)=FKE*PK(30)
      A(7,7)=XK
      A(7,8)=0.
      A(8,1)=0.
      A(8,2)=FTE*FCE*PK(33)
      A(8,3)=FTE*(PK(31)+FOE*PK(32))
      A(8,4)=0.
      A(8,5)=-FTE*(FOE*PK(31)+2.*FOE*FOE*PK(32)+2.*FCE*FCE*PK(33))
      A(8,6)=0.
      A(8,7)=0.
      A(8,8)=T

CCC
      RETURN
      END
C
      SUBROUTINE MONTON(XI,N)
      implicit real*8 (a-h,o-z)
C
C      PARAMETER(NDIM=100)
      DIMENSION XI(N)
C
      DO 50 I=2,N
      IF(XI(I).LT.XI(I-1)) GOTO 50
      XI(I)=XI(I-1)-.05
50    CONTINUE
C
C
      RETURN
C
C
      E  N  D
C
      SUBROUTINE MULT(A,B,C,D,N,M)
      implicit real*8 (a-h,o-z)
C
C MULT MULTIPLIES MATRIX C (N*M) WITH SQUARE MATRIX B (N*N) AND PLACES
C THE RESULT IN MATRIX A. MATRIX D (N*M) IS USED FOR SCRATCH.
C FIRST DIMENSION MUST BE NDP.
C TIMING 850 MS FOR 40*40 MATRICES AND FORH.
C
      include 'parameter.inc'
      DIMENSION A(NDP,M),B(NDP,N),C(NDP,M),D(NDP,M)
C
C ZEROSET
      DO 100 J=1,M
      DO 100 I=1,N
100   D(I,J)=0.
C
C MULTIPLY
      DO 200 J=1,M
      DO 200 I=1,N
      DO 200 K=1,N
200   D(I,J)=B(I,K)*C(K,J)+D(I,J)
C
C RESTORE
      DO 300 J=1,M
      DO 300 I=1,N
300   A(I,J)=D(I,J)
C
      RETURN
      END
C
      SUBROUTINE OLDARC(IARCH)
      implicit real*8 (a-h,o-z)
C
C OLDARC RESTARTS FROM OLD ARCHIV FILE ON LUN 'IARCH'
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &   VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /TAUC/TAU(NDP),DUMTAU(NDP),JTAU
      COMMON /CTEFF/TEFF,FLUX /CG/G,KONSG
      DIMENSION ABUND(20)
C
C READ OLD ARCHIV FILE

      READ(IARCH) DUM
      READ(IARCH) TEF,FLX,GD,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITER,NEL,(ABUND(I),I=1,NEL)
      ITER=0
      READ(IARCH) DUM
      READ(IARCH) NTAU
      PI=3.14159
      DO 100 K=1,NTAU
      READ(IARCH) A,A,A,ZZ(K),TT(K),PPE(K),PG,PPR(K),PPT(K),ROSS(K)
      PP(K)=PG+PPR(K)+PPT(K)
      GG(K)=0.
      ZZ(K)=0.
      READ(IARCH) RRO,EMU,CP,CV,AGRAD,Q,U,V,ANCONV
      FFC(K)=FLUX*ANCONV
      VV(K)=V
      DD(K)=0.
      IF(VV(K).GT.0.) DD(K)=2.*PI*FFC(K)/(PALFA*RRO*CP*TT(K)*VV(K))
      READ(IARCH) DUM
      READ(IARCH) DUM
100   CONTINUE
C
C INTERPOLATE SURFACE TEMPERATURE
      YB=TAU(1)/TAU(2)
      YA=1.-YB
      TT(1)=YA*TT(1)+YB*TT(2)
C
C PRINT MODEL INFORMATION
      WRITE(7,51) IARCH
51    FORMAT('0MODEL FROM OLD ARCHIV FILE ON LUN',I3)
      GRAV=log10(G)
      WRITE(7,50) TEFF,GRAV,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6
50    FORMAT(' TEFF=',F6.0,'  log10(G)=',F5.2,'  IDENTIFICATION ',6A4)
C
C MOVE INTEGERS TO HALF INTEGERS
      DO 110 K=2,NTAU
      FFC(NTAU+2-K)=0.5*(FFC(NTAU+2-K)+FFC(NTAU+1-K))
      VV(NTAU+2-K)=0.5*(VV(NTAU+2-K)+VV(NTAU+1-K))
      DD(NTAU+2-K)=0.5*(DD(NTAU+2-K)+DD(NTAU+1-K))
      ZZ(NTAU+2-K)=0.5*(ZZ(NTAU+2-K)+ZZ(NTAU+1-K))
110   CONTINUE
C
      REWIND IARCH
      RETURN
      END
C
      SUBROUTINE OLDSTA
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
C 'OLDSTA/NEWSTA' SAVES MODEL DATA ON AN ASCII FILE.
C
c      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
c     &      VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /STATEC/B(11*NDP),RO(NDP),NTAU,ITER
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C READ OLD STATE FROM unit 16

      READ(16,*) NTAU,ITER
      ITER=0
      DO 100 I=1,11
      IST=(I-1)*NDP

      READ(16,*) (B(IK),IK=IST+1,IST+NTAU)

100   CONTINUE
      
      MINTAU=10.*NDP+1
      MAXTAU=10.*NDP+NTAU
      DO 104 IK=MINTAU,MAXTAU
104   B(IK)=2.3025851*B(IK)    
C transform tauscale to ln for use in SOLVE
      IF (MOLOLD.EQ.1) READ(16,*) ((FOLD(K,I),I=1,8),K=1,NTAU)

      RETURN
      END
C
C WRITE NEW STATE
      SUBROUTINE NEWSTA
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
C
C 'OLDSTA/NEWSTA' SAVES MODEL DATA ON AN ASCII FILE.
C
      COMMON /STATEC/B(11*NDP),RO(NDP),NTAU,ITER
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CIT/IT,ITMAX
C
C      REWIND 17

      OPEN(UNIT=17,FILE='arcivaab.dat',STATUS='UNKNOWN')
      WRITE(17,*) NTAU,IT
      MINTAU=10.*NDP+1
      MAXTAU=10.*NDP+NTAU
      DO 104 IK=MINTAU,MAXTAU
      K=IK-MINTAU+1
      B(IK)=log(TAU(K))
      B(IK)=0.4342945*B(IK)    
C transform tauscale to log10 for store
104   CONTINUE
      DO 100 I=1,11
      IST=(I-1)*NDP
      WRITE(17,*) (B(IK),IK=IST+1,IST+NTAU)
100   CONTINUE
      do 120 K=1,NTAU
120   WRITE(17,122) (FOLD(K,I),I=1,8)
122   FORMAT (1P7E11.3)
      CLOSE(17)
      RETURN
      END
C
C
      SUBROUTINE ONFROM(LUN,LR)
      implicit real*8 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /STATEC/A(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /SPACE1/APP(13*NDP),TAULNP(NDP),SP1DUM((2*NDP-5)*NDP)
     &     ,space1dum2(3*ndp),space1dum3(2*ndp*ndp)
      DIMENSION AP(10*NDP)
      EQUIVALENCE (AP(1),APP(1))
      DATA LREC/0/
C
C RESUME PERMITS RESTART WITH NEW TAU-SCALE.
C
      REWIND LUN
      LREC=0
1     READ(LUN,END=2,ERR=99) APP,TAULNP,NTAUP
C
C OLD RECORD FORMAT
      GO TO 98
C
C NEW RECORD FORMAT
99    BACKSPACE LUN
      READ(LUN,END=2) AP,TAULNP,NTAUP,ITER
98    CONTINUE
      LREC=LREC+1
      GO TO 1
2     BACKSPACE LUN
C
C PRINT HEADING
      WRITE(7,52) LUN,LREC
52    FORMAT('0SAVED VALUES FROM LUN',I3,', RECORD',I2)
C
201   CONTINUE
      KP=2
      DO 100 K=1,NTAU
C
101   IF(TAULN(K).LE.TAULNP(KP)) GO TO 102
      IF(KP.EQ.NTAUP) GO TO 102
      KP=KP+1
      GO TO 101
C
102   CONTINUE
      P=(TAULN(K)-TAULNP(KP-1))/(TAULNP(KP)-TAULNP(KP-1))
      Q=1.-P
C
      DO 100 I=1,10
      J=(I-1)*NDP
      B=AP(KP+J)
      C=AP(KP+J-1)
      IF(B.LE.0..OR.C.LE.0.) GO TO 103
      A(K+J)=EXP(P*log(B)+Q*log(C))
      GO TO 100
103   A(K+J)=P*B+Q*C
100   CONTINUE
C
      RETURN
C
      ENTRY SAVEON(LUN,LR)
C
      LREC=LREC+1
      WRITE(LUN) A,TAULN,NTAU,ITER
C
C PRINT MESSAGE
      WRITE(7,51) LUN,LREC
51    FORMAT('0SAVED ON LUN',I3,', RECORD',I2)
C
      RETURN
C
C RESUME
      ENTRY RESUME(LUN,LR)
      DO 200 I=1,LR
200   READ(LUN) AP,TAULNP,NTAUP,ITER
      LREC=LR
      WRITE(7,52) LUN,LREC
      GO TO  201
C
      END
C
      SUBROUTINE OPAC(J,X,S)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE ADMINISTERS COMPUTATION OF OPACITIES.
C        J IS THE WAVELENGTH NUMBER (IN THE XL ARRAY), X AND S ARE THE
C        NORMALIZED ABSORPTION AND SCATTERING COEFFICIENTS, RESPECTIVELY
C
      include 'parameter.inc'
C
      DIMENSION PRXC(NDP),PRXO(NDP),PRXW(NDP),PRXT(NDP),PRSC(NDP)
      DIMENSION ABSK(NDP),SPRID(NDP),ABSK1(NDP),SPRID1(NDP)
      DIMENSION X(NDP),S(NDP)
      DIMENSION V(NDP),CON(NDP)
      EQUIVALENCE (V(1),CON(1))
      common/carciv/ larciv    !=1 if called from arciv, otherwise = 0
      COMMON /TAUC/TAU(NDP),DUMT(NDP),JJTAU
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     &      VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),
     & NTAU,ITER
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/CLINE1/XLINLO,XLINUP,TSKAL(30),
     &             PESKAL(30),IPEBEG(30),IPEEND(30),LINUN,NTSKAL,NPSKAL
      COMMON/CLINE2/FACTOR(NDP,2,2),IT(NDP),IPE1(NDP),IPE2(NDP)
       COMMON/CLINE4/ILINE
      COMMON /CXMAX/XMAX
      COMMON /ODFCD/  CONODF(NDP,8)
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
C
      COMMON /COPINF/ SUMOP(maxosmol,NDP),SUMKAP(maxosmol,NDP)
      COMMON /CMOL1/DMUDMU(9),FOE,XMUDMUD(6)
      COMMON /CI4/dumdum,IDUMDUM(96),MOLH,JUMP
      COMMON /CARC3/ F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON /DENSTY/ BPZ(NDP),PRH2O(NDP)
      common/ci5/abmarcs(18,ndp),anjon(18,5),h(5),part(18,5),
     *dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      INTEGER MOLH, JUMP
      COMMON /CMETPE/ PPEL(NDP), METPE
      
      common /ggchemmu/ggmu(NDP),ggrho(NDP),ppsum(ndp),ppappsum(ndp),
     &   ppnonappsum(ndp),tg(ndp),pges(ndp)
     &  ,ppat1sum(ndp),ppat2sum(ndp),ppmolsum(ndp),ppgs(ndp)
      character atnames*2, molnames*8,molnames2*4
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)

      common /ggchemdetabs / f1_dt(ndp), f5_dt(ndp), 
     >                       rC(ndp), rMg(ndp), rAl(ndp), 
     >                       rSi(ndp), rHe(ndp), ro_dt(ndp)
C      
! Dust
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common /dustplot/ x_gas(ndp,nwl), s_gas(ndp,nwl), gas_opac(ndp)
C
      dimension x_cloud(ndp,nwl), s_cloud(ndp,nwl),
     *  x_cont(ndp,nwl), s_cont(ndp,nwl), x_line(ndp,nwl)
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP,FIRST
C                                                              
      DIMENSION PELOG(NDP)
c
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
C      COMMON/COPPRR/ xconop(120,10),xlineop(120,10)    !100wn,10dpt
C      COMMON/CONLIN/rconop(nwl,ndp),rlineop(nwl,ndp)
      COMMON/COPsum/ SSUM(NDP),XSUM(NDP),CONSUM(NDP)
c
       dimension tau_lambda(5)
       DATA JJ/4/
       DATA FIRST/.TRUE./
       data tau_lambda/ 8.0e+4, 9.0e+4,1.0e+5, 1.1e+5, 1.2e+5 /
C      EXP10(A)=EXP(2.302585*A)
C
      IF (J.EQ.1) THEN
         IF (WLOS(1).LT.XL(1,2)) THEN
           PRINT*,' WLOS(1),XL(1,2),XL(2,2) = ',WLOS(1),XL(1,2),XL(2,2)
           STOP ' WLOS(1) (==1.E8/WNEND) < first continuums point'
         END IF
         if (idust == 1) then
            call dust_opac_eps_interp
         end if
         IF (NOSMOL.GT.0) THEN
            CALL OSTABLOOK
         END IF       !if nosmol > 0
      END IF          !if j=1 (i.e. calling OPAC(j,x,s) for first oswn)

c
C ROSSELAND OPACITY OR NOT
       JTAU=JJTAU
C
C MURIEL EST GENTILLE
C
C***********************************************************************
C SPECIAL
      IF (J.EQ.1) THEN
      DO 1234 K=1,JTAU
      PRXC(K)=0.0
      PRXO(K)=0.0
      PRXT(K)=0.0
      PRXW(K)=0.0
      PRSC(K)=0.0
 1234 CONTINUE
      ENDIF
C END OF SPECIAL
C
      IF(J.EQ.1) JJ=4
      IF(J.NE.1)GO TO 9
      ISWITCH=1

C        COMPUTATION OF CONTINUOUS ABSORPTION COEFFICIENTS AND INTERPOLATION
      NEWT=1
      JMEM=1
      IMEM=2
C   ******** HERE WE ASSUME THAT THE FIRST SET IS USED FOR ROSSELAND MEAN
      JMEM1=2
      IMEM1=2
      CALL ABSKO(NEWT,JTAU,T,PPEL,IMEM,JMEM,ABSK,SPRID,-1)
      NEWT=0
      do k=1,jtau
      if(metpe.eq.1) then
       CALL ABSKO(NEWT,JTAU,T(k),PE(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      else if(metpe.eq.2) then
      if (T(k) .gt. 2000.) then 
       CALL ABSKO(NEWT,JTAU,T(k),PE(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      else            
       CALL ABSKO(NEWT,JTAU,T(k),PPEL(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      endif
      end if
      enddo
    9 IF(WLOS(J).LE.XL(JMEM1,IMEM1))GO TO 11

      JMEM=JMEM1
      IMEM=IMEM1
      JMEM1=JMEM1+1
      IF (IMEM.GT.10 .OR. JMEM1.GE.20) 
     *                STOP ' WNOS(1) < last continuums point '
      DO8 K=1,JTAU
      ABSK(K)=ABSK1(K)
    8 SPRID(K)=SPRID1(K)
      IF(JMEM1.LE.NL(IMEM1))GO TO 10
      JMEM1=1
      IMEM1=IMEM1+1

10    continue
      do k=1,jtau
      if(metpe.eq.1) then
       CALL ABSKO(NEWT,JTAU,T(k),PE(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      else if(metpe.eq.2) then
      if (T(k) .gt. 2000.) then
       CALL ABSKO(NEWT,JTAU,T(k),PE(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      else    
       CALL ABSKO(NEWT,JTAU,T(k),PPEL(k),IMEM1,JMEM1,ABSK1,SPRID1,-1)
      endif
      end if
      enddo
      

      GO TO 9
C        INTERPOLATION
11    CONTINUE

      DO12 K=1,JTAU
      
      DIFXL=(WLOS(J)-XL(JMEM,IMEM))/(XL(JMEM1,IMEM1)-XL(JMEM,IMEM))

      X(K)=(ABSK(K)+(ABSK1(K)-ABSK(K))*DIFXL)/XKAPR(K)
      x_cont(k,j) = x(k)*XKAPR(K)
      PRXC(K)=X(K)
      S(K)=(SPRID(K)+(SPRID1(K)-SPRID(K))*DIFXL)/XKAPR(K)
      s_cont(k,j) = s(k)*XKAPR(k)
      PRSC(K)=S(K)
   12 CONTINUE  
         
C
C        COMPUTATION OF LINE-ABSORPTION COEFFICIENTS
C
      IF (J.EQ.1) THEN
      DO 19 K=1,JTAU
          CONSUM(K) = 0.
          XSUM(K) = 0.
          SSUM(K) = 0.
19    CONTINUE
      END IF
C
      DO 468 K=1,JTAU
          CONSUM(K) = CONSUM(K) + CONOS(K,J)/XKAPR(K)
          XSUM(K) = XSUM(K) + X(K)
          SSUM(K) = SSUM(K) + S(K)
          x_line(k,j) = CONOS(k,j)
          X(K)=X(K)+CONOS(K,J)/XKAPR(K)
468   CONTINUE
C
C
      do k=1, jtau
      x_gas(k,j) = x(k)*xkapr(k)
      s_gas(k,j) = s(k)*xkapr(k)
      end do
! Addition of dust absorption & scattering
      if(idust == 1) then
      do k=1,jtau   
      x(k) = x(k)+(dust_abs(k,j)/xkapr(k))
      s(k) = s(k)+(dust_sca(k,j)/xkapr(k))
      end do
      end if


   25 CONTINUE
      DO  K=1,JTAU
      X(K)=X(K)*XMAX/(X(K)+XMAX)
      end do
      if ( idust==1) then
      
         do l=1, 5
            tau_c = 0.0
            if ((wlos(j)<= tau_lambda(l)) 
     >      .and. (wlos(j+1)> tau_lambda(l))) then
            if (l==1) then
            open(unit=4297, file='tau_cloud.dat', status="replace", 
     *       position="append", action="write")
            else 
            open(unit=4297, file='tau_cloud.dat', status="old", 
     *       position="append", action="write")
            end if
            write(4297,*) 'table for wavelength (in aa):'
            write(4297,*) wlos(j)
            do k=1, jtau
                  x_cloud(k,j) = dust_abs(k,j)
                  s_cloud(k,j) = dust_sca(k,j)
                  if (k>1) then
                        opac_k_b = opac_k
                  end if
                  opac_k = (x_gas(k,j) + s_gas(k,j) + 
     *             x_cloud(k,j) + s_cloud(k,j))/xkapr(k)
                  if (k==1) then
                  tau_c = opac_k * tau(1)
                  else
                  tau_c = tau_c + 
     *             0.5*(tau(k)-tau(k-1))*(opac_k_b+opac_k)
                  end if
                  write(4297,'(6e24.15)') x_cont(k,j), s_cont(k,j), 
     *            x_line(k,j),x_cloud(k,j), s_cloud(k,j), tau_c
            end do
            write(4297,*)
            write(4297,*)
            close(4297)
            end if

         end do
      end if
50    FORMAT(' X',6(' ***',1PE12.5))
C
C
      FIRST = .FALSE.

      RETURN
      END
C
C
      SUBROUTINE PEMAKE(T,PE,PG,PEX)
      implicit real*8 (a-h,o-z)
C
C 'PEMAKE-R' IS COMPATIBLE WITH 'PEMAKE' AND SOMEWHAT FASTER. A MODIFIED
C REGULA-FALSI PROCEDURE IS USED ON THE LOG-LOG PG-PE RELATION.
C 76.03.08  *NORD*
C
      DATA IT,N,EPS/0,20,1.E-3/
C
C START

      A=log(PE)
      PEX=PE
      CALL JON(T,PE,1,FA,RO,E,0)
C******WRITE(7,40) T,PG,PE,FA
40    FORMAT(' T,PG,PE,PGP=',4E11.4)
      IT=IT+1
      FA=log(FA/PG)
      IF(ABS(FA).LT.EPS) GOTO 101
      B=A-0.69*FA
      PEX=EXP(B)
C ONE PEMAKE ITERATION, CF. PEMAKE
      CALL JON(T,PEX,1,FB,RO,E,0)
C******WRITE(7,40) T,PG,PE,FB
      IT=IT+1
      FB=log(FB/PG)
      IF(ABS(FB).LT.EPS) GOTO 101
      X=B
C
C LOOP
      DO 100 I=1,N
      XOLD=X
C
C INTERPOLATE TO FIND NEW X
      X=A-(B-A)/(FB-FA)*FA
      PEX=EXP(X)
      IF(ABS(X-XOLD).LT.EPS) GOTO 101
      CALL JON(T,PEX,1,FX,RO,E,0)
C******WRITE(7,40) T,PG,PEX,FX
      IT=IT+1
      FX=log(FX/PG)
C
C CHECK IF A OR B CLOSEST TO X
      IF(ABS(A-X).LT.ABS(B-X)) GOTO 102
      A=X
      FA=FX
      GOTO 100
102   B=X
      FB=FX
C
C END OF LOOP
100   CONTINUE
      WRITE(7,51) N,T,PE,PG,A,B,FA,FB,EPS
51    FORMAT('0***PEMAKE, MAX ITER.: N,T,PE,PG,A,B,FA,FB,EPS=',
     * /,1X,I2,8E11.4)
      RETURN
C
C NORMAL END
101   CONTINUE
      RETURN
C
C COUNT ENTRY
      ENTRY PECNT
      WRITE(7,52) IT
52    FORMAT('0TOTAL NUMBER OF CALLS TO JON FROM PEMAKE-R =',I5)
      RETURN

      END
C
      FUNCTION QAS(H,XL,A,Z,PFISH,IFISH)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES THE ASYMPTOTIC PARTS OF THE PARTITION
C        FUNCTIONS FOLLOWING
C           BASCHEK ET AL., ABH. HAMB. VIII, 26 (1966) IF IFISH = 0
C           FISCHEL AND SPARKS, AP. J. 164, 359 (1971) IF IFISH = 1
C           (APPROXIMATING THE ZETA FUNCTIONS BY INTEGRALS).
C
C        XL=QUANTUM NUMBER FOR THE FIRST LEVEL OF THE ASYMPTOTIC PART
C        H=QUANTUM NUMBER OF THE CUT (FOR IFISH=0)
C        A=DZ(FISCHEL AND SPARKS)=ALFA(BASCHEK ET AL.)
C        PFISH=P(FISCHEL AND SPARKS), ONLY NECESSARY IF IFISH = 1
C
C
      COMMON/UTPUT/IREAD,IWRIT
C
C        WHICH TYPE
      IF(IFISH.GT.0)GO TO 1
C
C        BASCHEK ET AL.
      QAS=0.333333*(H*(H+1.)*(H+0.5)-XL*(XL+1.)*(XL+0.5)) +
     *                       A*(H-XL)+0.5*A*A*(H-XL)/(H*XL)
      RETURN
C
C        FISCHEL AND SPARKS
    1 P=PFISH*Z
C
C        FISCHEL AND SPARKS, EQ. (26)
      P2=P*P
      P3=P2*P
      IF(P.LE.XL)GO TO 2
      XLM1=XL-1.
      R2=XLM1*XLM1
      R3=R2*XLM1
      QAS=1.3333333*P3+0.5*P2+0.16666667*P+1.33333333*A*P-0.4*A*A/P-
     *0.33333333*R3-0.5*R2-0.16666667*XLM1-A*XLM1+0.5*A*A/XL
      RETURN
C
C        FISCHEL AND SPARKS, EQ. (27)
    2 AXL2=A/(XL*XL)
      QAS=P3*P/XL*(1.+AXL2*(0.33333333+0.1*AXL2))
      RETURN
      END
C
      FUNCTION QTRAV(TETA,HP,J,JA)
      implicit real*8 (a-h,o-z)
C
C        HERE THE PARTITION FUNCTIONS ACCORDING TO TRAVING ET AL., ABH. HAMB.
C        STERNW. VIII, 1 (1966) ARE COMPUTED. THE SYMBOLS ARE GIVEN
C        IN THE COMMENTS AT THE BEGINNING OF SUBROUTINE INJON.
C        FUNCTION QAS IS CALLED.
C
C        DIMENSIONS NECESSARY
C        A(5),ASDE(KMAX),H(5),QPRIM(KMAX)
C        KMAX IS THE TOTAL NUMBER OF ELECTRON CONFIGURATIONS.
C        DIMENSIONS OF ARRAYS IN COMMON /CI3/ ARE COMMENTED ON IN SUBROUTINE
C        INJON.
C
      DIMENSION ASDE(80),H(5),QPRIM(80)
      COMMON/CI3/ALFA(300),GAMMA(300),G0(45),G2(80),XION(80),XL(80),
     *JBBEG(45),JCBEG(45),NK(45),NL(80),IFISH
      COMMON/CI7/A(5),PFISH,ITP
C
C
C        STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)
C
      FLJ=J
      JB=JBBEG(JA)
      JC1=JCBEG(JA)
      NKP=NK(JA)
      QSUM=0.
C
C        WE START THE LOOP OVER DIFFERENT ELECTRON CONFIGURATIONS ('THE K-LOOP')
      DO5 K=1,NKP
      JC2=NL(JB)+JC1-1
C
C        IS TETA=PRECEDING TETA
      IF(ITP.GT.0)GO TO 4
      PRA=XION(JB)*TETA
      IF(PRA.LT.12.)GO TO 1
      ASDE(JB)=0.
      GO TO 2
    1 ASDE(JB)=G2(JB)*EXP10(-PRA)
C
    2 QPRIM(JB)=0.
      IF(NL(JB).LE.0)GO TO 4
      DO3 L=JC1,JC2
      PRE=GAMMA(L)*TETA
      IF(PRE.GT.12.)GO TO 3
      QPRIM(JB)=QPRIM(JB)+ALFA(L)*EXP10(-PRE)
    3 CONTINUE
    4 JC1=JC2+1
      QSUM=QPRIM(JB)+ASDE(JB)*QAS(HP,XL(JB),A(J),FLJ,PFISH,IFISH)
     *             +QSUM
    5 JB=JB+1
C        END OF 'THE K-LOOP'
      QTRAV=G0(JA)+QSUM
C
      RETURN
      END
C
      FUNCTION ROSSOP(T,PE, nlayer)
      implicit real*8 (a-h,o-z)
C
C 'ROSSOP' CALCULATES THE ROSSELAND MEAN OPACIY AS DEFINED ON THE
C FIRST WAVELENGTH SCALE IN 'INITAB' INPUT. THIS VERSION WITH
C SCATTERING.                     (13/02/89)
C
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
      DATA NEWT/2/
C
      CALL ABSKO(NEWT,1,T,PE,1,0,RSP,DUM,nlayer)
      NEWT=1
      ROSSOP=RSP
      RETURN
C
      END
C
C
      SUBROUTINE ROSSOS
      implicit real*8 (a-h,o-z)
C
      include 'parameter.inc'
C
      DIMENSION X(NDP),S(NDP),SUMW(NDP),sumabs(ndp)
     *   ,sumwy(ndp),sumwyxs(ndp),dummy(ndp)
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
C      REAL*8 Y,YA,SUMW,ROSSO,PTAUO
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common /cdrift/ idust, ieps, idustopac, icloud_conv
C      COMMON/CONLIN/rconop(nwl,ndp),rlineop(nwl,ndp)

C CALCULATE DETAILED ROSSELAND MEAN
C      KL=1
C      DUMMY(1)=ROSSOP(TT(1),PPE(1))


      DO 116 K=1,NTAU
C        KL=K
        SUMW(K)=0.
        ROSSO(K)=0.
        dummy(k) = 1.
C        sumwy(k) = 0.
C        sumwyxs(k) = 0.
C        sumabs(k) = 0.
C        DUMMY(k)=ROSSOP(TT(K),PPE(K))
116   CONTINUE
C          write(7,*) 
C     *     ' j,k,wlos(j),ya,xkapr(k),xkapr(k)*x(k),xkapr(k)*s(k)',
C     *     ' rosso(k),rconop(j,k),rlineop(j,k),tt(k) '

      !print*, "opac call in rossop"
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSO(K)=ROSSO(K)+WLSTEP(J)*YA/(xkapr(k)*(X(K)+S(K)))

117   CONTINUE
1171  format(i5,i3,1p8e12.3,0pf8.0)
C
C      write(7,*) ' k,sum(x(k)+s(k)),rosso(k),sumw(k),sumwy,sumwyxs = '
      DO 111 K=1,NTAU
C        write(7,*) k,sumabs(k),rosso(k),sumw(k),sumwy(k),sumwyxs(k)
        ROSSO(K)=SUMW(K)/ROSSO(K)
        PTAUO(K)=GRAV*TAU(K)/ROSSO(K)
111   CONTINUE
C

C      write(7,*)' k,pp,ptauo,rosso,tt ='
C      do 112 k=1,ntau
C      write(7,1172)k,pp(k),ptauo(k),rosso(k),tt(k)
C112   continue
1172  format(i5,1p3e12.3,0pf8.0)
C      write(7,*) ' ******** end from rossos '

      RETURN
      END
C
      SUBROUTINE SCALE(LUN)
      implicit real*8 (a-h,o-z)
C
C THIS ROUTINE SCALES THE MODEL PRESENTLY IN THE COMMON /STATEC/ (AND
C IN THE MODEL FILE) TO THE NEW EFFECTIVE TEMPERATURE AND GRAVITY GIVEN
C IN THE ORDINARY COMMON'S /CTEFF/ AND /CG/.
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     *VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),ROSTAT(NDP),
     *NTAU,ITER
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CTEFF/TEFF,FLUX /CG/GRAV,KONSG
      COMMON /CMETPE/ PPEL(NDP), METPE
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      COMMON /MIXC/PALFA,PBETA,PNY,PY    /CSTYR/ MIHAL,NOCONV
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CISPH/ISPH
C
C READ OLD MODEL

      CALL OLDSTA
      ITER=0
C
C TRY TO READ OLD MODEL FILE ARCHIV RECORDS
      READ(LUN,END=100,ERR=100) IDUM
      READ(LUN,END=100,ERR=100) TEFOLD,FDUM,GRVOLD
      GO TO 110
C
C UNSUCCESSFUL, NO OLD ARCHIV DATA
100   CONTINUE
      DO 101 K=1,NTAU
      IF(TAU(K).GE.0.95) GO TO 102
101   CONTINUE
102   TEFOLD=TT(K)/1.08
      GRVOLD=1.
C EMPIRICAL FITTING
C
C PRINT
110   G1=log10(GRAV)
      G2=log10(GRVOLD)
      WRITE(7,50) TEFF,G1,TEFOLD,G2
50    FORMAT('0SCALING TO TEFF=',F6.0,' log10(G)=',F5.2,' FROM TEFF='
     &,F6.0,' log10(G)=',F5.2)
C
C SCALE THE TEMPERATURE
      DO 111 K=1,NTAU
      PPR(K)=PPR(K)*(TEFF/TEFOLD)**4
      FFC(K)=FFC(K)*(TEFF/TEFOLD)**4
111   TT(K)=TT(K)*TEFF/TEFOLD
C
C TRY TO AVOID DIVIDE BY ZERO IN TRYCK
      DO 112 K=1,NTAU
        CROSS(K)=1.
  112 CONTINUE
C
C INTEGRATE PRESSURE EQUATION
        if (isph.eq.1) then 
          CALL TRYCK_sph
        else 
          CALL TRYCK
        end if
C
C SCALE DD AND VV
      DO 120 K=NOCONV,NTAU
      KL=K
      IF(FFC(K).EQ.0.) GO TO 120
      TMEAN=.5*(TT(K)+TT(K-1))

      PEMEAN=.5*(PPE(K)+PPE(K-1))
      PRMEAN=.5*(PPR(K)+PPR(K-1))
      ROSSMN=.5*(ROSS(K)+ROSS(K-1))
      CALL TERMON(K,TMEAN,PEMEAN,PRMEAN,PG,PGT,PGPE,RO,ROT,ROPE,
     &      CP,ADIA,Q)
      HSCALE=(PG+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      THETA=OMEGA/(1.+PY*OMEGA**2)
      GAMMA=CP*RO/(8.*STEFAN*TMEAN**3*THETA)
      DD(K)=(GRAV*HSCALE*Q/PNY*(PALFA**2*RO*CP*TMEAN/(2.*PI*FFC(K)))**2
     &)**(-.333333)
      VV(K)=PALFA*SQRT(GRAV*HSCALE*Q/PNY*DD(K))
      GG(K)=GAMMA*VV(K)
120   CONTINUE
      RETURN
      END
C 
      SUBROUTINE SETDIS(NLB,XLB,LMAX,IFIRST,ILAST)
      implicit real*8 (a-h,o-z)
C
C        THIS SUBROUTINE DISTRIBUTES THE NLB WAVELENGTHS GIVEN IN XLB IN
C        WAVELENGTH SETS, WITH MAX. LMAX WAVENLENGTHS IN EACH. THE FIRST
C        SET NUMBER IS IFIRST, THE LAST ILAST. IF MORE SETS ARE NECESSAR
C        EXECUTION IS STOPPED WITH A PRINT-OUT.
C
      DIMENSION XLB(500)
      COMMON/CXLSET/XL(20,10),NSET,NL(10)
      COMMON/CLINE1/XLINLO,XLINUP,TSKAL(30),
     &             PESKAL(30),IPEBEG(30),IPEEND(30),LINUN,NTSKAL,NPSKAL
      COMMON/CLINE3/GLAMD(100),JLBDS
      COMMON/UTPUT/IREAD,IWRIT
        COMMON/CLINE4/ILINE
C
      NLP=1
      NSET=IFIRST
      DO2 J=1,NLB
       IF(ILINE.LE.0)GOTO1
      IF(XLB(J).LT.XLINLO.OR.XLB(J).GT.XLINUP)GO TO 1
      DO11 K=1,JLBDS
      KMEM=K
      IF(XLB(J).LE.GLAMD(K))GO TO 12
   11 CONTINUE
   12 IF(J.EQ.NLB)GO TO 1
      IF(XLB(J+1).LE.GLAMD(KMEM))GO TO 2
    1 CONTINUE
      XL(NLP,NSET)=XLB(J)
      IF(J.EQ.NLB)GO TO 3
      NLP=NLP+1
      IF(NLP.LE.LMAX)GO TO 2
      NL(NSET)=NLP-1
      NLP=1
      NSET=NSET+1
      IF(NSET.GT.ILAST)GO TO 4
    2 CONTINUE
C
    4 WRITE(IWRIT,200)
      STOP 'SETDIS 1'
C
    3 NL(NSET)=NLP
      RETURN
C
  200 FORMAT(1H ,'TOO FEW SETS ALLOWED OR TOO MANY WAVELENGTH POINTS',
     &' WANTED IN SETDIS')
      END
C
      SUBROUTINE STARTM
      implicit real*8 (a-h,o-z)
C
C 'STARTM' FINDS A STARTING MODEL, WITH CONVECTIVE FLUX, AND WITH
C APPROXIMATELY THE CORRECT EFFECTIVE TEMPERATURE.
C
C DIMENSIONS
      include 'parameter.inc'
C
      DIMENSION PTAU(NDP)
C
C STATE VARIABLES
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),
     &               ZZ(NDP),DD(NDP),VV(NDP),FFC(NDP),
     &               PPE(NDP),TT(NDP),TAULN(NDP),ROSTAT(NDP),NTAU,ITER
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
C OWN COMMONS
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,
     &              IDRAB6,IARCH
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA IEDIT/8/
C
C STATEMENT FUNCTIONS
      TF(TAUX)=TEFF*EFF*(.75*(TAUCNV+TAUX))**.25-
     &          EFF*DTBLNK*(1.-2.*TAUX/(TAUBLN+TAUX))
C
C TIME

      CALL CLOCK
C
C START UP
      READ(5,66) TAUCNV,DTBLNK,TAUBLN
      EFF=1.09*(0.75*(TAUCNV+1.))**(-0.25)
      TSURF=TF(TAU(1))
      IF(TSURF/TEFF.LT.0.60) EFF=EFF*0.60*TEFF/TSURF
      DT1=1.E10
      TEOLD=-DT1
94    CONTINUE
      
C
C THE PRESSURE BOUNDARY CONDITION IS PP(1)/PP(2) = PGFACT. TO FIND
C APPROXIMATE STARTING VALUES WE START WITH AN ARBITRARY ELECTRON PRESSU
C AND INTEGRATE THE PRESSURE EQUATION FROM 1 TO 2. THIS GIVES NEW GAS PR
C AND NEW ELECTRON PRESSURES. ITERATION IS PERFORMED UNTIL SUFFICIENT AC
C IS OBTAINED.
      TT(1)=TF(TAU(1))
      TT(2)=TF(TAU(2))
      PPR(1)=4./3.*STEFAN*TT(1)**4/CLIGHT
      PPR(2)=4./3.*STEFAN*TT(2)**4/CLIGHT
      PPE(1)=1.E-3
      PGFACT=(TAU(1)/TAU(2))**0.667
      PPE(2)=PPE(1)/SQRT(PGFACT)
      DPEL=0.
      DFDPEL=3.
      DO 90 K=1,20
      KL=1
      ROSS(1)=ROSSOP(TT(1),PPE(1),1)
      PG1=PGC
      PG2=PG1/PGFACT
      KL=2
      CALL PEMAKE(TT(2),PPE(2),PG2,PPE(2))
      ROSS(2)=ROSSOP(TT(2),PPE(2),2)
      DPG=GRAV*DTAULN(2)*.5*(TAU(1)/ROSS(1)+TAU(2)/ROSS(2))
      CALL ZEROF(log((PG2-PG1)/DPG),DPEL,DFDPEL)
      PPE(1)=PPE(1)*EXP(DPEL)
      IF(ABS(DPEL).LT.0.001) GOTO 88
90    CONTINUE
88    CONTINUE
      PP(1)=PG1+PPR(1)
      PP(2)=PG2+PPR(2)
      KL=1
      ROSS(1)=ROSSOP(TT(1),PPE(1),1)
      KL=2
      ROSS(2)=ROSSOP(TT(2),PPE(2),2)
      PTAU(1)=GRAV*TAU(1)/ROSS(1)
      PTAU(2)=GRAV*TAU(2)/ROSS(2)
C
C TAU LOOP
      RO=0.
      DO 99 K=1,NTAU
      KL=K
      PPT(K)=0.
      ZZ(K)=0.
      GG(K)=0.
      DD(K)=0.
      VV(K)=0.
      FFC(K)=0.
      IF(K.LE.2) GO TO 99
C
C NEW TEMPERATURE
      TT(K)=TF(TAU(K))
      PPR(K)=1.33*STEFAN*TT(K)**4/CLIGHT
      PPE(K)=PPE(K-1)*TAU(K)/TAU(K-1)
      PPT(K)=PPT(K-1)
C
C ITERATE THREE TIMES
      DO 92 IT=1,3
C
C COMPUTE ROSSELAND MEAN AND NEW PRESSURE
      KL=K
      ROSS(K)=ROSSOP(TT(K),PPE(K),k)
      PTAU(K)=GRAV*TAU(K)/ROSS(K)
      PP(K)=PP(K-1)+.5*DTAULN(K)*(PTAU(K)+PTAU(K-1))
      DLNP=log(PP(K)/PP(K-1))
C
C NEW ELECTRON PRESSURE
      PPT(K)=MIN(0.5*PP(K),PPT(K))
      PG=PP(K)-PPR(K)-PPT(K)
      PPEK=PPE(K)
      KL=K
      CALL PEMAKE(TT(K),PPEK,PG,PPE(K))
      IF(K.LE.NOCONV) GO TO 97
C
C CONSIDER CONVECTION
      TMEAN=.5*(TT(K)+TT(K-1))
      PEMEAN=.5*(PPE(K)+PPE(K-1))
      PRMEAN=.5*(PPR(K)+PPR(K-1))
      KL=K
      CALL TERMON(K,TMEAN,PEMEAN,PRMEAN,PG,PGT,PGPE,RO,ROT,
     &      ROPE,CP,ADIA,Q)
      IF(IT.EQ.1) GO TO 95
C
C USE OLD GRAD-ADIA AT SECOND IT-LOOP
      IF(FFC(K).EQ.0.) GO TO 97
      GO TO 96
95    CONTINUE
      GRAD=log(TT(K)/TT(K-1))/DLNP
      IF((GRAD.LE.ADIA.AND.VFIX.EQ.0.).OR.K.LE.NOCONV.OR.PALFA.EQ.0.)
     * GOTO 97
C
C IF CONVECTION ACCORDING TO THESE CONDITIONS AND ESTIMATED CONVECTIVE
C FLUX GREATER THAN TOTAL FLUX, ADJUST GRADIENT UNTIL ESTIMATED CONVECTI
C FLUX EQUALS TOTAL FLUX.
      HSCALE=(PG+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSS(K-1)
      THETA=OMEGA/(1.+PY*OMEGA**2)
      GAMMA=CP*RO/(8.*STEFAN*TT(K-1)**3*THETA)
      VV(K)=VVMLT(GRAD-ADIA,GRAV*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      GG(K)=GAMMA*VV(K)
      DD(K)=GG(K)/(1.+GG(K))*(GRAD-ADIA)
      FFC(K)=PALFA*RO*CP*TT(K-1)*VV(K)*DD(K)/2./PI
      FFCM=FLUX*TAU(K)/10.
      IF(FFCM.GT.FLUX) FFCM=FLUX
      IF(FFC(K).LT.FFCM.AND.VFIX.EQ.0.) GO TO 97
      DD(K)=(GRAV*HSCALE*Q/PNY*(PALFA**2*RO*CP*TT(K-1)/(2.*PI*FFCM))**2
     &)**(-.333333)
      FFC(K)=FFCM
      VV(K)=PALFA*SQRT(GRAV*HSCALE*Q/PNY*DD(K))
      GG(K)=GAMMA*VV(K)
C
C ADJUST TT(K) TO GIVE GRAD-ADIA
96    CONTINUE
      GRAD=ADIA+(1.+GG(K))/GG(K)*DD(K)
      TT(K)=TT(K-1)*EXP(GRAD*DLNP)
      EFF=EFF*TT(K)/TF(TAU(K))
      PPR(K)=1.33*STEFAN*TT(K)**4/CLIGHT
      PPT(K)=MIN(0.5*PP(K),PBETA*RO*VV(K)**2)
      PG=PP(K)-PPR(K)-PPT(K)
      PPEK=PPE(K)
      KL=K
      CALL PEMAKE(TT(K),PPEK,PG,PPE(K))
C
97    CONTINUE
      IF(ABS(TAU(K)-1.).GT.0.01) GOTO 93
C
C IMPROVE FIT TO EFFECTIVE TEMPERATURE
      TE=TT(K)-1.07*TEFF-DTBLNK
      IF(ABS(TE).LT.100.) GOTO 93
      DT1=DT1*TE/(TEOLD-TE)
      TT(1)=TT(1)+DT1
      EFF=EFF*TT(1)/TF(TAU(1))
      TEOLD=TE
      GO TO 94
93    CONTINUE
C
C END OF IT-LOOP
92    CONTINUE
C
C END OF TAU LOOP, PRINT.
99    continue
C
C TIME
      CALL CLOCK
C
      ITER=0
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1',74X,'STARTM(',I1,')',5X,'ITERATION',I3,5X,6A4)
49    FORMAT(13('1234567890'),'123')
51    FORMAT(I8,1P10E12.4,I4)
52    FORMAT(T12,'TAU',T24,'PRAD',T36,'PTURB',T48,'PTOT',T60,'GAMMA',
     *T72,'DELTA',T84,'VCONV',T96,'FCONV',T108,'PE',T120,'TEMP')
61    FORMAT(' STARTING VALUES, TAUCNV,DTBLNK,TAUBLN,EFF=',4E10.3)
66    FORMAT(7(7X,F8.0))
      END
C
      SUBROUTINE TABS(NT,T)
      implicit real*8 (a-h,o-z)
C
C        THIS ROUTINE COMPUTES  FACTORS FOR INTERPOLATION IN T (TETA IF
C        ITETA(KOMP) IS GREATER THAN ZERO) IN THE ABKOF TABLE, INITIATED BY
C        SUBROUTINE INABS. CONCERNING THE OTHER CONTROL INTEGERS, SEE INABS.
C        THE RESULTING FACTORS ARE PUT IN AFAK. THE NUMBER OF FACTORS FOR
C        THE COMPONENT KOMP AT TEMPERATURE T(NTP) IS GIVEN IN
C        NOFAK((NKOMP-KOMPR)*(NTP-1)+KOMP-KOMPR). HERE KOMPR IS THE NUMBER
C        OF COMPONENTS WITH T-INDEP. COEFFICIENTS. NOFAK=0 MEANS THAT THE
C        ABSORPTION COEFFICIENT SHOULD BE =0. NPLATS (INDEX AS FOR NOFAK)
C        GIVES THE ARRAY INDEX OF THE TEMPERATURE POINT AT WHICH THE
C        INTERPOLATION IN ABKOF SHOULD START.
C
C        NT=NUMBER OF TEMPERATURES
C        T= ARRAY OF TEMPERATURES
C
C        DIMENSIONS NECESSARY
C        AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM),T(1)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS OF ARRAYS IN COMMON /CA1/
C        AND /CA2/ ARE COMMENTED ON IN SUBROUTINE INABS.
C        IFADIM SHOULD BE AT LEAST =(NKOMP-KOMPR)*NT, WHERE NKOMP IS THE NUMBER
C               OF COMPONENTS, KOMPR THE NUMBER OF TEMPERATURE-INDEPENDENT
C               COMPONENTS AND NT THE NUMBER OF TEMPERATURE POINTS (IN THE PARA-
C               METER LIST).
C        KFADIM SHOULD BE AT LEAST =KOMPR*NT+(NKOMP-KOMPR)*NT*NUM, WHERE NUM IS
C               BETWEEN 2 AND 3 AND DEPENDENT ON THE TYPE OF TEMPERATURE
C               INTERPOLATION USED.
C
C
      include 'parameter.inc'
C
C      PARAMETER (IFADIM=1000,KFADIM=4000)
      DIMENSION T(NT)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(30,2),TBOT(30,2),IDEL(30),ISVIT(30),ITETA(30),
     *KVADT(30),MAXET(30),MINET(30),NTM(30,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(4000),KOMPLA(600),KOMPR,KOMPS,NKOMP
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
C
      IFAK=1
      KFAK=1
      NSVIT=1
C        THIS IS JUST A DUMMY STATEMENT TO GIVE NSVIT A FORMAL VALUE
C
      DO81 NTP=1,NT
      TP=T(NTP)
     
      KFAK=KFAK+KOMPR
      DO81 KOMP=KOMPS,NKOMP

     
      IF(ISVIT(KOMP).GT.0)GO TO (51,61,70),NSVIT
      IF(ITETA(KOMP).LE.0)GO TO 2
    1 TS=5040./T(NTP)
      
      GO TO 3
    2 TS=T(NTP)
C
C        SEARCHING
      
    
    3 IF((TS-TBOT(KOMP,1)).GE.0.)GO TO 10
      IF(MINET(KOMP).LE.0)GO TO 70
C
C        EXTRAPOLATION DOWNWARDS
      IF(NEXTT.GT.0)WRITE(IWRIT,200)TS,KOMP
      INTA=1
      AP=(TS-TBOT(KOMP,1))/DELT(KOMP,1)
      IP=0
      GO TO 60
C
C        SEARCHING CONTINUES

   10 INTAP=1
      IDP=IDEL(KOMP)
      
      DO11 I=1,IDP
      
      AP=(TS-TBOT(KOMP,I))/DELT(KOMP,I)
      IP=INT(AP)
      INTA=IP+INTAP
      INAP=NTM(KOMP,I)-1+INTAP
      IF(INTA.LE.INAP) GO TO 20
   11 INTAP=INAP+1
      IF(MAXET(KOMP).LE.0)GO TO 70
C
C        EXTRAPOLATION DOWNWARDS
      
      IF(NEXTT.GT.0)WRITE(IWRIT,200)TS,KOMP
      INTA=INAP
      IP=NTM(KOMP,IDP)-1
      GO TO 60
C
   20 IF(KVADT(KOMP).LE.0)GO TO 60
C
C        QUADRATIC INTERPOLATION
   21 IF(INTA.LT.INAP)GO TO 50
      INTA=INTA-1
      IP=IP-1
C
   50 DXX1=AP-DFLOAT(IP)
      DXX2=DXX1-1.
      DXX3=DXX1-2.
      A1=DXX2*DXX3*0.5
      A2=-DXX1*DXX3
      A3=DXX1*DXX2*0.5
   51  AFAK(KFAK)=A1
      AFAK(KFAK+1)=A2
      AFAK(KFAK+2)=A3
      NPLATS(IFAK)=INTA
      NOFAK(IFAK)=3
      IFAK=IFAK+1
      KFAK=KFAK+3
      NSVIT=1
      GO TO 80
C
C        LINEAR INTER/EXTRAPOLATION
   60 A2=AP-DFLOAT(IP)
      A1=1.-A2
   61 AFAK(KFAK)=A1
      AFAK(KFAK+1)=A2
      NPLATS(IFAK)=INTA
      NOFAK(IFAK)=2
      IFAK=IFAK+1
      KFAK=KFAK+2
      NSVIT=2
      GO TO 80
C
C        OUTSIDE TABLE. ABS.COEFF. SHOULD BE = 0
   70 IF(NUTZT.GT.0)WRITE(IWRIT,201)TS,KOMP
      NOFAK(IFAK)=0
      IFAK=IFAK+1
      NSVIT=3
C
   80 CONTINUE
      IF(KFAK.GT.KFADIM)GO TO 90
      IF(IFAK.GT.IFADIM+1)GO TO 91
   81 CONTINUE
C
      GO TO 92
   90 WRITE(IWRIT,202)KFAK,KFADIM,NT
      STOP 'TABS 1'
   91 WRITE(IWRIT,203)IFAK,IFADIM,NT
      STOP 'TABS 2'
   92 CONTINUE
C
  200 FORMAT(33H EXTRAPOLATION IN TABS, T (TETA)=,E12.5,5X,
     *12HCOMPONENT NO,I5)
  201 FORMAT(24H ZERO IN TABS, T (TETA)=,E12.5,5X,12HCOMPONENT NO,I5)
  202 FORMAT(6H KFAK=,I5,5X,11H GT KFADIM=,I5,5X,12HIN TABS, NT=,I5)
  203 FORMAT(6H IFAK=,I5,5X,11H GT IFADIM=,I5,5X,12HIN TABS, NT=,I5)
      RETURN
      END
C
      SUBROUTINE TAET(T,PE,PG,RO,E)
      implicit real*8 (a-h,o-z)
C
C TAET SIMULATION
      
      CALL JON(T,PE,1,PG,RO,E,0)

      RETURN
      END
C
      SUBROUTINE TAUSCA
      implicit real*8 (a-h,o-z)
C
C 'TAUSCA' INITIATES A TAU SCALE FROM INPUT LOGTAU AND LOGTAU-DIFFERENCE
C *NORD*
      
C
      include 'parameter.inc'
C
      DIMENSION TAULNX(NDP)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /STATEC/DUM1(10*NDP),TAULN(NDP),RO(NDP),NTAU,ITER
     
      K=1
C
C READ LOGTAU AND LOGTAU-DIFFERENCE
      READ(5,50) T1,D1
      !print*, T1, D1
C
1     CONTINUE
      READ(5,50) T2,D2
      !print*, T2, D2
      TLIM=T2-.5*D1
      
      T=T1
2     CONTINUE
C
C NEW TAU POINT
      TAU(K)=10.**T
      K=K+1
      T=T+D1
      IF(T.LT.TLIM) GO TO 2
C
C NEW LOGTAU-DIFFERENCE
      T1=T2
      D1=D2
      IF(D1.GT.0.) GO TO 1
C
C END OF TAU SCALE
      TAU(K)=10.**T
      JTAU=K
CUGJ      NTAU=JTAU  commented 12.6.90 to allow NTAU(input-model) .ne. JTAU
C
      DO 3 K=1,JTAU
3     TAULNX(K)=log(TAU(K))
      DO 4 K=2,JTAU
4     DTAULN(K)=TAULNX(K)-TAULNX(K-1)
C

      RETURN
50    FORMAT(5(7X,F8.0))
      END
C
      SUBROUTINE TERMO(k,T,PE,PRAD,P,RO,CP,CV,TGRAD,Q,U2)
      implicit real*8 (a-h,o-z)
       include 'parameter.inc'
C
C
C        RUTINEN BERAEKNAR  OVANSTAAENDE STORHETER UTGAAENDE
C        FRAAN T, PE OCH PRAD (STRAALNINGSTRYCKET). INGEN AV STOR-
C        HETERNA GES ELLER ERHAALLES LOGARITMERAD. METODEN MED DIF-
C        FERENSFORMLER AER VARDAYAS
C        VI ANVAENDER RUTINEN TAET.
C        *** OBSERVERA. OCKSAA TAET MAASTE ARBETA MED OCH GE
C        I C K E  L O G A R I T M E R A D E  S T O R H E T E R . *******
C
C VERSION OF 73.02.05. DOUBLE PRECISION ADDED. *NORD*
C
C 13-OCT-1998 14:19:53.82:
C vi fick negativa Cv i flera modeller. Det beror
C(troligen) paa att naer du beraeknar de partiella derivatorna du
C behoever foer att beraekna Cv saa tar man ett alltfoer kort steg
C ( 0.001 av T resp P ). Vi oekade till 0.005 eller 0.01 och daa gick
C det mycket baettre! (Antagligen fluktuerar part. deriv. starkt eller
C saa aer det helt enkelt numeriskt brus vid laaga temperaturer).
C Kjell

C
      DIMENSION PGH(4),ROH(4),TH(4),EH(4),HH(4),PH(4),EP(4),PRAH(4)
      COMMON /CMETPE/ PPEL(NDP), METPE
      common /cu2warning/nu2warning
      nu2warning = 0

C
      xk_boltz=1.380649d-16
      xm_p=1.66053906660d-24 

      DEREP=0.01
      DERET=0.001
      DELPE=DEREP*PE
      DELT=DERET*T
      PINV=1./(2.*DELPE)
      TINV=1./(2.*DELT)
C     
      
      CALL JON(T,PE,1,PG,RO,EPP,0)
      P=PRAD+PG
C
      PEP=PE-DELPE
      CALL JON(T,PEP,1,PGH(1),ROH(1),EP(1),0)
      PRAH(1)=PRAD
      TH(1)=T
      PEP=PE+DELPE
      CALL JON(T,PEP,1,PGH(2),ROH(2),EP(2),0)
      PRAH(2)=PRAD
      TH(2)=T
      TP=T-DELT
      CALL JON(TP,PE,1,PGH(3),ROH(3),EP(3),0)
      PRAH(3)=PRAD*(1.-4.*DERET)
      TH(3)=TP
      TP=T+DELT
      CALL JON(TP,PE,1,PGH(4),ROH(4),EP(4),0)
      PRAH(4)=PRAD*(1.+4.*DERET)
      TH(4)=TP

      if (METPE .EQ. 2 .and. T.lt.2000.) then
        call calc_adiaindex(k, PG, T, gamma, xmmw)
        cp = gamma/(gamma-1)*xk_boltz/(xmmw*xm_p)
        cv = cp/gamma

C       Squared soundspeed
        xcs2inv = xmmw*xm_p/xk_boltz
        RO=pg/T*xcs2inv

        ROH(1)=PGH(1)/TH(1)*xcs2inv 
        ROH(2)=PGH(2)/TH(2)*xcs2inv  
        ROH(3)=PGH(3)/TH(3)*xcs2inv  
        ROH(4)=PGH(4)/TH(4)*xcs2inv   
C
        DO I=1,4
          PH(I)=PRAH(I)+PGH(I)
        end do
C
        DROT=(ROH(4)-ROH(3))*TINV
        DROP=(ROH(2)-ROH(1))*PINV
        DHP=(HH(2)-HH(1))*PINV
        DPT=(PH(4)-PH(3))*TINV
        DPP=(PH(2)-PH(1))*PINV
        DPGT=(PGH(4)-PGH(3))*TINV
        DPGP=(PGH(2)-PGH(1))*PINV
C
        HJALP=DROT-DROP*DPT/DPP
        TGRAD=-P*HJALP/(CP*RO*RO)
        Q=-T/RO*(DROT-DROP*DPGT/DPGP)
        U2=CP*DPP/(CV*DROP)

      else
C     Case, where METPE not 2, normal stellar marcs
C
        DO1 I=1,4
        EH(I)=3.*PRAH(I)/ROH(I)+EP(I)
        PH(I)=PRAH(I)+PGH(I)
    1   HH(I)=EH(I)+PH(I)/ROH(I)
C
        DET=(EH(4)-EH(3))*TINV
        DEP=(EH(2)-EH(1))*PINV
        DROT=(ROH(4)-ROH(3))*TINV
        DROP=(ROH(2)-ROH(1))*PINV
        DHT=(HH(4)-HH(3))*TINV
        DHP=(HH(2)-HH(1))*PINV
        DPT=(PH(4)-PH(3))*TINV
        DPP=(PH(2)-PH(1))*PINV
        DPGT=(PGH(4)-PGH(3))*TINV
        DPGP=(PGH(2)-PGH(1))*PINV
C
        CV=DET-DEP*DROT/DROP
        CP=DHT-DHP*DPT/DPP
        HJALP=DROT-DROP*DPT/DPP
        TGRAD=-P*HJALP/(CP*RO*RO)
        Q=-T/RO*(DROT-DROP*DPGT/DPGP)
        U2=CP*DPP/(CV*DROP)
      end if  
      if (U2.lt.0.) then
           nu2warning = 1
           write(6,*) ' U2 becomes negative in TERMO '
           write(6,*) ' T,RO,PE,PRAD,P,CV,CP,Q,DPP,DROP,TGRAD,U2:'
           write(6,777) T,RO,PE,PRAD,P,CV,CP,Q,DPP,DROP,TGRAD,U2
777       format(1p6e12.3)
      end if
C
      RETURN
      END
C
      SUBROUTINE TERMON(k,T,PE,PRAD,PG,DPGT,DPGP,RO,DROT,
     &      DROP,CP,TGRAD,Q)
      implicit real*8 (a-h,o-z)
       include 'parameter.inc'
C
C
C        RUTINEN BERAEKNAR  OVANSTAAENDE STORHETER UTGAAENDE
C        FRAAN T, PE OCH PRAD (STRAALNINGSTRYCKET). INGEN AV STOR-
C        HETERNA GES ELLER ERHAALLES LOGARITMERAD. METODEN MED DIF-
C        FERENSFORMLER AER VARDAYAS
C        VI ANVAENDER RUTINEN TAET.
C        *** OBSERVERA. OCKSAA TAET MAASTE ARBETA MED OCH GE
C        I C K E  L O G A R I T M E R A D E  S T O R H E T E R . *******
C
C VERSION OF 73.02.05. DOUBLE PRECISION ADDED. *NORD*
C
C
      DIMENSION PGH(4),ROH(4),TH(4),EH(4),HH(4),PH(4),EP(4),PRAH(4)
      COMMON /CMETPE/ PPEL(NDP), METPE
C
      xk_boltz=1.380649d-16
      xm_p=1.66053906660d-24 

      DEREP=0.01
      DERET=0.001
      DELPE=DEREP*PE
      DELT=DERET*T
      PINV=1./(2.*DELPE)
      TINV=1./(2.*DELT)
C     
      
      CALL JON(T,PE,1,PG,RO,EPP,0)
      P=PRAD+PG
C
      PEP=PE-DELPE
      CALL JON(T,PEP,1,PGH(1),ROH(1),EP(1),0)
      PRAH(1)=PRAD
      TH(1)=T
      PEP=PE+DELPE
      CALL JON(T,PEP,1,PGH(2),ROH(2),EP(2),0)
      PRAH(2)=PRAD
      TH(2)=T
      TP=T-DELT
      CALL JON(TP,PE,1,PGH(3),ROH(3),EP(3),0)
      PRAH(3)=PRAD*(1.-4.*DERET)
      TH(3)=TP
      TP=T+DELT
      CALL JON(TP,PE,1,PGH(4),ROH(4),EP(4),0)
      PRAH(4)=PRAD*(1.+4.*DERET)
      TH(4)=TP

      if (METPE.EQ.2 .and. T.lt.2000.) then
        call calc_adiaindex(k, PG, T, gamma, xmmw)
        cp = gamma/(gamma-1)*xk_boltz/(xmmw*xm_p)
        cv = cp/gamma


C       Squared soundspeed
        xcs2inv = xmmw*xm_p/xk_boltz
        RO=pg/T*xcs2inv

        ROH(1)=PGH(1)/TH(1)*xcs2inv 
        ROH(2)=PGH(2)/TH(2)*xcs2inv  
        ROH(3)=PGH(3)/TH(3)*xcs2inv  
        ROH(4)=PGH(4)/TH(4)*xcs2inv  
C
        DO I=1,4
          PH(I)=PRAH(I)+PGH(I)
        end do
C
        DROT=(ROH(4)-ROH(3))*TINV
        DROP=(ROH(2)-ROH(1))*PINV
        DHP=(HH(2)-HH(1))*PINV
        DPT=(PH(4)-PH(3))*TINV
        DPP=(PH(2)-PH(1))*PINV
        DPGT=(PGH(4)-PGH(3))*TINV
        DPGP=(PGH(2)-PGH(1))*PINV
      else
C     Case, where METPE=2, normal stellar marcs
     
        DO1 I=1,4
        EH(I)=3.*PRAH(I)/ROH(I)+EP(I)
        PH(I)=PRAH(I)+PGH(I)
    1   HH(I)=EH(I)+PH(I)/ROH(I)
C
        DET=(EH(4)-EH(3))*TINV
        DEP=(EH(2)-EH(1))*PINV
        DROT=(ROH(4)-ROH(3))*TINV
        DROP=(ROH(2)-ROH(1))*PINV
        DHT=(HH(4)-HH(3))*TINV
        DHP=(HH(2)-HH(1))*PINV
        DPT=(PH(4)-PH(3))*TINV
        DPP=(PH(2)-PH(1))*PINV
        DPGT=(PGH(4)-PGH(3))*TINV
        DPGP=(PGH(2)-PGH(1))*PINV
C
        CV=DET-DEP*DROT/DROP
        CP=DHT-DHP*DPT/DPP
      end if      
      HJALP=DROT-DROP*DPT/DPP
      TGRAD=-P*HJALP/(CP*RO*RO)
      Q=-T/RO*(DROT-DROP*DPGT/DPGP)
      U2=CP*DPP/(CV*DROP)
C
      RETURN
      END

      SUBROUTINE calc_adiaindex(k, p, T, gamma, xmmw)
      implicit real*8 (a-h,o-z)
        include 'parameter.inc'
        integer i,j,k
        character atnames*2, molnames*8
        character diagempta*2(22), diagemptm*8(543)
        logical thermo_initialized, molactive, atactive
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
        common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
        
        DATA thermo_initialized / .FALSE. /
C     from init_thermo      
        common / thermodata / xmoltha(543,3,7),xmolthTlim(543,4),
     >               xattha(22,3,7),xatthTlim(22,4), thexp(7),
     >               xtemplimlow,xtemplimup,molactive(543),
     >               atactive(22)
        
        cpRsum=0.0d0
        psum=0.0d0
        diagempta(:)=''
        diagemptm(:)=''

        call ggchem(k, T, p)   ! sets ggmuk for here and for termo and termon
        xmmw = ggmuk

C       ADS: initialize thermodynamic data (NASA polynomials)      
        if (.NOT. thermo_initialized) then
          call init_thermo()
          thermo_initialized = .TRUE.
        end if
        
C       The nasa polynomials are only defined between 200K and 20000K        
C       We will use the boundary values if we cross this temperature
        TN = min(max(xtemplimlow+1.0d-60,T),xtemplimup-1.0d-60)

C       Molecules        
        do i=1,543
C         Determine the temperature regime of the nasa polynomial        
          if ((TN.ge.xmolthTlim(i,1)).and.(TN.lt.xmolthTlim(i,2))) 
     &    then 
             treg = 1
          elseif ((TN.ge.xmolthTlim(i,2)).and.(TN.lt.xmolthTlim(i,3))) 
     &    then
             treg = 2
          elseif ((TN.ge.xmolthTlim(i,3)).and.(TN.le.xmolthTlim(i,4))) 
     &    then
             treg = 3
          else 
             ! If the temperature is outside the defining range of the polynomial,
             ! we will not use it
             cycle 
          endif 
C         Calculate the polynomial
          if (molactive(i)) then
            polsum = 0.0d0
            do j=1,7
              polsum=polsum+xmoltha(i,treg,j)*TN**thexp(j)
            end do
            cpRsum=cpRsum+polsum*ppallmol(k,i)
            psum=psum+ppallmol(k,i)
          end if
        end do

C       ATOMS
        do i=1,22
C         Determine the temperature regime of the nasa polynomial        
          if ((TN.ge.xatthTlim(i,1)).and.(TN.lt.xatthTlim(i,2))) 
     &    then 
             treg = 1
          elseif ((TN.ge.xatthTlim(i,2)).and.(TN.lt.xatthTlim(i,3))) 
     &    then
             treg = 2
          elseif ((TN.ge.xatthTlim(i,3)).and.(TN.le.xatthTlim(i,4))) 
     &    then
             treg = 3
          else 
             ! If the temperature is outside the defining range of the polynomial,
             ! we will not use it
             cycle 
          endif 
C         Calculate the polynomial
          if (atactive(i)) then
            polsum = 0.0d0
            do j=1,7
              polsum=polsum+xattha(i,treg,j)*TN**thexp(j)
            end do
            cpRsum=cpRsum+polsum*ppallat(k,i)
            psum=psum+ppallat(k,i)
          end if
        end do
        if (psum/pgesk .LT. 0.97d0) then
C         NOTE: if the code fails here, we need to change the matching between ggchem and nasa polynomials.
C           This can be easily done by diagnosing the here generated file and improving the matching in
C           init_thermo.        

          print*, 'ISSUE IN CALCULATING THERMODYNAMIC'
          print*, 'QUANTITIES (HEATCAPACITIES)!'
          print*, 'The matching between ggchem indexes and nasa9'
          print*, 'could not match a representative fraction.'
          print*, 'nasa9 file: data/nasa9.dat'
          print*, 'ggchem file: pp.dat'
          print*, 'Find the matched AND NOT MATCHED molecules/atoms'
          print*, 'here: nasa2ggchem.dat'          
          
          open(unit=8895, file='nasa2ggchem.dat', status='replace')
          write(8895, *) '# Matched elements:'
          write(8895, *) merge(atnames, diagempta, atactive)
          write(8895, *) '# UNmatched elements:'
          write(8895, *) merge(diagempta, atnames, atactive)
          write(8895, *) ''
          write(8895, *) '# Matched molecules:'
          write(8895, *) merge(molnames, diagemptm, molactive)
          write(8895, *) '# UNmatched molecules:'
          write(8895, *) merge(diagemptm, molnames, molactive)
          close(8895)
          stop
        endif
        
C       Final heat capacities        
        cpR = cpRsum/psum
        cvR = cpR - 1
        gamma = cpR/cvR

      END SUBROUTINE calc_adiaindex

C
      SUBROUTINE TINT(N,X,Y,XINT,YINT)
      implicit real*8 (a-h,o-z)
C
C        DENNA ENDIMENSIONELLA INTERPOLATIONSRUTIN ARBETAR MED SUCCESIV
C        HALVERING. I PRINCIP BAADE INTER- OCH EXTRAPOLERAR DEN VILLIGT.
C        VARNING SKRIVS DOCK UT VID EXTRAPOLATION. INTERPOLATIONERNA SKE
C        MED TREPUNKTSFORMEL, SUBR.  I N P 3  .
C        ***** OBSERVERA. TABELLERNA SKALL VARA V A E X A N D E ******
C
      DIMENSION X(N),Y(N),ARG(3),FUNK(3)
      COMMON/UTPUT/IREAD,IWRIT
C
      NOEV=N
      NED=1
    1 NP=(NED+NOEV)/2
      IF(XINT-X(NP))2,2,3
    2 NOEV=NP
      GO TO 4
    3 NED=NP
    4 IF(NOEV-NED-2)5,5,1
    5 IF(NOEV-2)7,7,6
    6 J=NOEV-3
      GO TO 8
    7 J=NOEV-2
    8 DO9 K=1,3
      JP=K+J
      ARG(K)=X(JP)
    9 FUNK(K)=Y(JP)
      IF(ARG(1)-XINT)11,11,12
   11 IF(ARG(3)-XINT)12,13,13
   12 WRITE(IWRIT,200)XINT,ARG
  200 FORMAT(37H VARNING. EXTRAPOLATION I TINT. XINT=,E12.5,3X,4HARG=,
     *3E12.5)
   13 CALL INP3(ARG,FUNK,XINT,YINT)
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE SOLVE(NEW)
      implicit real*8 (a-h,o-z)
C
C SOLVE PERFORMS ONE NEWTON-RAPSON ITERATION ON THE MODELATMOSPHERE PROB
C INCLUDING LOCAL CONVECTION.THE STATE OF THE ATMOSPHERE IS DESCRIBED BY
C NUMBER OF VARIABLES SUCH AS TEMPERATURE,ELECTRON PRESSURE,TOTAL PRESSU
C CONVECTIVE FLUX ETC..TO EACH VARIABLE CORRESPONDS A CERTAIN CONDITIONA
C EQUATION WICH DETERMINES THAT VARIABLE,ASSUMING THE OTHER VARIABLES BE
C KNOWN.
C
C NAMING CONVENTION.THE VARIABLES HAVE NAMES WITH A DOUBLE OCCURANCE OF
C FIRST LETTER,CORRECTIONS (TO BEE COMPUTED IN THIS ITERATION) HAVE SING
C OCCURANCE OF FIRST LETTER.RIGHTHANDSIDENAMES BEGIN WITH R.
C
C VARIABLES ARE CENTERED ON INTEGER AND HALFINTEGER TAU-POINTS AS INDICA
C BY 'I' OR 'H' ON THE FOLLOWING COMMENT CARDS.
C
C VARIABLE CORRECTION                                           HALF-INT
C PPR      PR         RADIATION PRESSURE                             I
C PPT      PT         TURBULENT PRESSURE                        H
C PP       P          TOTAL PRESSURE                                 I
C GG       G          CONVECTIVE EFFICIENCY,GAMMA               H
C ZZ       Z          GEOMETRIC HEIGTH                          H
C DD       D          GRADIENT DIFFERENCE,DELTA-DELTAPRIME      H
C VV       V          CONVECTIVE VELOCITY                       H
C FFC      FC         CONVECTIVE FLUX                           H
C PPE      PE         ELECTRON PRESSURE                              I
C TT       T          TEMPERATURE                                    I
C XJ                  MEAN INTENSITY                                 I
C
C NAMES OF THE PARTIAL DERIVATIVES ARE FORMED WITH A FIRST PART FROM THE
C EQUATION TO WICH IT BELONGS,AND A SECOND PART WICH IS THE NAME OF THE
C VARIABLE WITH RESPECT TO WICH THE DERIVATIVE IS TAKEN.(DOUBLE OCCURANC
C IF NECESSARY TO AVOID CONFUSION).
C
C PRPR,PRT
C PTPT,PTV,PTPE,PTT
C PPP,PPPE,PPTT
C GGG,GV,GPE,GT
C ZZZ,ZPE,ZT
C DDD,DP,DG,DPE,DT
C VVV,VZ,VD,VPE,VT
C FCFC,FCP,FCD,FCV,FCPE,FCT
C PEPE,PEPR,PEPT,PEP,PET
C TXJ,TFC,TTT
C
C STATE VARIABLES
      include 'parameter.inc'
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     *VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),stro(ndp),
     *NTAU,ITER
      common /ckdtpe/dpex,kdtpe
      common /dpeset/ dpein,dtin, pe_corr(ndp)
C
C DIMENSIONS
      DIMENSION PTAU(NDP),ROSSP(NDP),SUMW(NDP),ROSST(NDP),ROSSPE(NDP)
     *,XL(500),W(500)
     *,XJ1(NDP),XJ2(NDP),XJ3(NDP),XJT1(NDP),XJT2(NDP),XJT3(NDP)
     *,XJPE1(NDP),XJPE2(NDP),XJPE3(NDP)
     *,PR(NDP),PRT(NDP,NDP),PRJ(NDP)
     *,PT(NDP),PTV(NDP),PTPE(NDP),PTT(NDP)
     *,P(NDP),PPPE(2*NDP),PPTT(2*NDP)
     *,GV(NDP),GPE(NDP),GT(NDP)
     *,DP(2*NDP),DG(NDP),DPE(NDP,NDP),DT(NDP,NDP),DV(NDP)
      DIMENSION D(NDP),DTS(2*NDP),DPS(2*NDP),DPES(2*NDP)
     *,V(NDP),VD(NDP),VPE(NDP,NDP),VT(NDP,NDP)
     *,FC(NDP),FCD(NDP),FCV(NDP),FCPE(NDP,NDP),FCT(NDP,NDP)
     *,PE(NDP),PEPE(NDP,NDP),PET(NDP,NDP)
     *,T(NDP),TTT(NDP,NDP),TPE(NDP,NDP),TJ1(NDP),TJ2(NDP),TTTS(NDP)
     *,PRPE(NDP,NDP),RTS(NDP)
     *,SCRATC(NDP,NDP),DBPL(NDP)
     *,XT(NDP),ST(NDP),DLNX(NDP),XLOG(NDP)
     *,XPE(NDP),SPE(NDP)
     *,RPR(NDP),RP(NDP),RD(NDP),RV(NDP),RFC(NDP),RPE(NDP),RT(NDP)
     *,TAUTAU(NDP) 
      LOGICAL NEWV, exist
      real*8 a,b,c,aa,bb,cc,aaa,ccc,STBZ,IR,RS,R,TIR, bpl_var
      character*24 idmodl
      logical:: first_call_rad = .True.
      logical:: file_exists
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C NLAM=NUMBER OF LAMBDAPOINTS, XL=LAMBDAPOINTS, W=INTEGRATIONWEIGHTS.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /CPRINT/NPRINT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CVAAGL/XL,W,NLAM
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /CMETPE/ PPEL(NDP), METPE

      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     * STEFAN
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
CUGJ FFR in excess     COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE,FFR(NDP)
      dimension ffr(ndp)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE 
C OWN COMMONS
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,dumtran(4*ndp),idumtran(3)
      COMMON /CTRAN2/EJ(NDP),TOTEJ(NDP),TOTIR(NDP),E(NDP),TOTE(NDP),
     & E_P(NDP),EJ_P(NDP), EJ_PLANET(NDP), E_PLANET(NDP), f(nwl)
      COMMON /CANGLE/XMU(6),XMU2(6),H(6),MMU_PP
      COMMON /CSURF/HSURF,Y1(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP) /RHOC/RHO(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CARC2/T,FC,FLUXME(NWL),TAU5(NDP),INORD
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CI8/PGC,ROC,EC
      COMMON /NEWMO/NEWMOD
      COMMON /MASSE/RELM
      COMMON /CORRECT/TDIFF,TCONV,KORT
      COMMON /CIT/IT,ITMAX
C
C SPACE ALLOCATION
      COMMON /SPACE1/XJ1,XJ2,XJ3,TJ1,TJ2,XJT1,XJT2,XJT3,PRJ,TTT,PRT
     &  ,XJPE1,XJPE2,XJPE3,TPE,PRPE
      COMMON /SPACE2_PP/FCT,FCPE,PET,PEPE
! Dust
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
! IRRADIATION
      COMMON /CROSSIR/ROSSIR(NDP),ROSSPIR(NDP),SUMWIR(NDP),TAURIR(NDP)
      COMMON /CPLANCKIR/PLANCKIR(NDP),PLANCKPIR(NDP),SUMWPIR(NDP),
     &     TAUPIR(NDP)
      COMMON /CIR/TAUIR(NDP),XIR(NDP,NWL),SIR(NDP,NWL),synspec(nwl),
     * DTAUIR(NDP), DTAUPLANET(NDP), DTAUP(NDP)
     
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom         !irrin=1~comp.irrad,steff=rad*
      common /dustplot/ x_gas(ndp,nwl), s_gas(ndp,nwl), gas_opac(ndp)
      
      common /ch4/ nch4
      character*24 file_name
      character*8 file_id
      common /ctcorlast/tcorlast 

      DATA IVERS,IEDIT/21,1/
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      common /noneq_output/ krome_output,krome_debug,krome_return
      common /photochem/ FLUX_RAD(ndp,nwreal) !second dimension should be nwtot, in most cases 7949
      common /starspec/ stellar_spectrum(nwreal),index_wlambda

      
      if (krome_on.EQ.1) then
       if (krome_photo_on.EQ.1) then
        if (first_call_rad.eq..True.) then
         if (nwtot.ne.nwreal) then
         write(*,*) "Number of wavelengths after OS calc", nwtot, 
     >   "not the same as number of wavelengths for FLUX_RAD", nwreal,
     >   "please adjust nwreal to match nwtot"
         stop
         endif
        FLUX_RAD(:,:)=0 !initialize radiative flux array with zeros in case somethign goes wrong in the intensity calculations        

        !write(*,*) "Open krome_flux_rad"
        inquire(file="krome_flux_rad.dat",exist=file_exists)

        if (file_exists) then
        write(*,*) "Found krome_flux_rad.dat, write data into FLUX_RAD"
        open(unit=7373,file="krome_flux_rad.dat",status='old',readonly)
        read(7373,*) !read first line before actually writing the data into FLUX_RAD
        do k=1,ntau 
         do j=1,nwreal
          read(7373,*) dummy_k,dummy_wl,FLUX_RAD(k,j) !first two entries dont matter just dummies
         enddo
        enddo
        close(7373)
        else
        write(*,*) "Did not find krome_flux_rad.dat, consider using
     > such a file for better convergence" 
        endif 
         if (krome_debug.eq.1) then
            open(unit=7676,file='BPL_sun.dat')
            open(unit=7777,file='BPL_upper.dat')
            open(unit=7878,file='XJ_upper.dat')
            open(unit=7979,file='XJ_lower.dat')
          endif        
C       
        first_call_rad=.False.
        endif
       endif
      endif       



C     
C IN THIS SECTION THE MEAN INTENSITY IS ELIMINATED IN THE TRANSPORT EQUA
C LEAVING THE EXPLICIT TEMPERATURE DEPENDANCE OF FLUX AND RADIATION PRES
C IN THE MATRICES TTT AND PRT.
C
      IF (mmu_pp.GT.nrays) STOP  
     & ' solve: increase nrays in parameter.inc'
      IF (mmu_pp.GT.6) STOP ' solve: increase dimension to mmu_pp'
      ITER=ITER+1
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
C
C ZEROSET
      !if (idust==1) call dust_eps
      DO 110 I=1,NTAU
      PPRG(I) = 0.
      RT(I)=0.
      RPR(I)=0.
      ROSSP(I)=0.
      SUMW(I)=0.
      DO 110 J=1,NTAU
      TTT(I,J)=0.
      TPE(I,J)=0.
      PRT(I,J)=0.
      PRPE(I,J)=0.
110   CONTINUE
      kdtpe = 0
C
C 
C CALCULATE DETAILED ROSSELAND MEAN
      REWIND 11
      KL=1
      !print*, "rossop line 8118"
      DUMMY=ROSSOP(TT(1),PPE(1), 1)
      PGA=PGC

      DO 116 K=1,NTAU
        KL=K
        !print*, "rossop line 8126"
        !print*, "gas pressure ", pp(k)
        DUMMY=ROSSOP(TT(K),PPE(K), k)
        ROSS(K)=1.0
        RHO(K)=ROC
        SUMW(K)=0.
        ROSSP(K)=0.
116   CONTINUE
      !print*, "opac call for ross "
      pe_corr(1:ntau) = 0.0
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(11) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSP(K)=ROSSP(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
C         if (j/500*500.eq.j  .and.  k.eq.10)
C     *         write(7,1171) j,wlos(j),wlstep(j),ya,y,x(k),s(k)
117   CONTINUE
1171  format(i5,1p8e12.3)
      REWIND 11
      
C
C TEMPERATURE AND ELECTRON PRESSURE PERTURBATIONS.
C KEEP THEM SMALL, TO STAY ON THE LINEAR PART.
C      DTX=0.001
C      DPEX=0.001
      DTX=dtin
      DPEX=dpein
      DO 111 K=1,NTAU
        KL=K
        T(K)=TT(K)*DTX
        PE(K)=PPE(K)*DPEX
        TT(K)=TT(K)+T(K)
        ROSSP(K)=SUMW(K)/ROSSP(K)
        SUMW(K)=0.
        ROSST(K)=0.
111   CONTINUE
      kdtpe = 1
C

C FIRST WAVELENGTH LOOP, TO CALCULATE XT,ST AND SAVE.
      REWIND 12
      pe_corr(1:ntau) = 0.0
      DO 112 J=1,NWTOT
       ! print*, "opac call for XT and ST "
        CALL OPAC(J,X,S)
        WRITE(12) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 112 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
          XIR(K,J)=X(K)
          SIR(K,J)=S(K)
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 112
          ROSST(K)=ROSST(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
112   CONTINUE
c        gas_opac(1:ntau) = 0.0
c       !integrate the gas opacities for plots
c       do k=1, nwtot
c             gas_opac(1:ntau) = gas_opac(1:ntau) 
c      *       + x_gas(1:ntau, k)*wlstep(k)*(10**-8.)
c      *       + s_gas(1:ntau, k)*wlstep(k)*(10**-8.)
c       end do
c       open(unit=2396, file='gas_opac.dat', status='unknown')
c       do n=1, ntau
c        write(2396, '(1E16.8,2x,1E15.6)') tt(n), gas_opac(n)
c       end do
c       close(2396)
      DO 113 K=1,NTAU
        ROSST(K)=SUMW(K)/ROSST(K)
        SUMW(K)=0.
        ROSSPE(K)=0.
        TT(K)=TT(K)-T(K)
        PPE(K)=PPE(K)+PE(K)
        pe_corr(k) = pe(k)
113   CONTINUE
        kdtpe = 2            !information to tstgem about computing dpg/dpe
      REWIND 12
      CALL CLOCK

C
C SECOND WAVELENGTH LOOP, TO CALCULATE XPE,SPE AND SAVE.
      REWIND 14
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1),1)
      PGPE=PGC
      !print*, "opac call for XPE and SPE "
      DO 114 J=1,NWTOT
        
        CALL OPAC(J,X,S)
        WRITE(14) X,S
        Y=((WLOS(J)/1.E4)**2)**3
        DO 114 K=1,NTAU
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 114
          ROSSPE(K)=ROSSPE(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
114   CONTINUE
      kdtpe = 3
      REWIND 14
      CALL CLOCK
C
C FROM THIS POINT ON, ROSS() HOLDS THE TRUE ROSSELAND MEAN.  CROSS HOLDS
C THE RATIO OF THE TRUE TO APPROXIMATE MEANS, WHICH ARE NEEDED IN TRYCK.
C      write (6,*) 'rosseland values'
C      write(7,*)' k,tt(k),ppe(k),ross(k),cross(k) in solve ='
      DO 123 K=1,NTAU
        KL=K
        ROSSPE(K)=SUMW(K)/ROSSPE(K)
        PPE(K)=PPE(K)-PE(K)
        !print*, "rossop line 8126"
        !print*, "gas pressure ", pp(k)
        CROSS(K)=ROSSP(K)/ROSSOP(TT(K),PPE(K),k)
        ROSS(K)=ROSSP(K)
        PTAU(K)=GRAV*TAU(K)/ROSS(K)
C        write (6,'(1x,i3,4(1pe12.3))')
C     &     k,ross(k),cross(k),rosst(k),rosspe(k)
C      write(7,1171)k,tt(k),ppe(k),ross(k),cross(k)
123   CONTINUE
      CALL CLOCK
      

C
C RIGHT HAND SIDE IN PRESSURE EQUATION
        KL=1
      CALL TAET(TT(1),PPE(1),PG,RO,DUM)
      DLNP=1./(1.+(ROSSPE(1)-ROSS(1))/ROSS(1)*PGA/(PGPE-PGA))
      RP(1)=GRAV*TAU(1)/(ROSS(1)*DLNP)+PPR(1)-PP(1)
C SIMPSONS RULE
      DO 101 K=2,NTAU
      F0=PTAU(K-1)
      F1=FOUR(PTAU,TAULN,K,NTAU)
      F2=PTAU(K)
      RP(K)=(F0+4.*F1+F2)*DTAULN(K)/6.-(PP(K)-PP(K-1))
101   CONTINUE
C
C TIME
      CALL CLOCK
      MSA=0
C     CALL MSLEFT(MSA)
      open(unit=960,file='wlos.dat',status='replace')
     

C WAVELENGTH LOOP
      
      FTOT=0.
      write(960,*) '  J,  wlos[mu],   wn[cm-1]'
      DO 150 J=1,NWTOT
      wlambda = wlos(j)
      wlam_mu = wlos(j)/1.e4
      wn_cm1 = 1.e8/wlos(j)
      if(J/10*10.eq.J) write(960,996) J, wlam_mu, wn_cm1
      index_wlambda = j
      write(960,996) K, J, wlos(J), WLSTEP(J), Y
996   FORMAT(i5, i6, f10.2, f10.2, f10.2) 
C
      DO K=1,NTAU
            BPLAN(K)=BPL(TT(K),WLOS(J))
            DBPL(K)=DIVBP(TT(K),WLOS(J))  
      END DO

C
C CALCULATE OPACITY DERIVATIVES AT CONSTANT GAS PRESSURE
      READ (11) X,S
      READ (12) XT,ST
      READ (14) XPE,SPE
      DO 131 K=1,NTAU
        X(K)=X(K)/ROSS(K)
        S(K)=S(K)/ROSS(K)
        XT(K)=XT(K)/ROSST(K)
        ST(K)=ST(K)/ROSST(K)
        XPE(K)=XPE(K)/ROSSPE(K)
        SPE(K)=SPE(K)/ROSSPE(K)
        XLOG(K)=log10(X(K))
        XT(K)=log(XT(K)/X(K))*X(K)/T(K)
        ST(K)=log(ST(K)/S(K))*S(K)/T(K)
        XPE(K)=log(XPE(K)/X(K))*X(K)/PE(K)
        SPE(K)=log(SPE(K)/S(K))*S(K)/PE(K)
        DLNX(K)=XT(K)*(TT(K)/2.3)/X(K)
131   CONTINUE
C
C TIME
      MS=MSA

      MSOPAC=MS-MSA
C
C SOLVE TRANSPORTEQUATION WITH OLD STRATIFICATION.
      CALL TRANEQ
      DO 132 K=1,NTAU
      IF(XT(K)*(XJ(K)-BPLAN(K)).GT.X(K)*DBPL(K))
     & XT(K)=X(K)*DBPL(K)/(XJ(K)-BPLAN(K))
132   CONTINUE
      MS=MSA
C     CALL MSLEFT(MSA)
      MSTRAN=MS-MSA
C
C FLUX TO PRINT
      HSURF=MAX(HSURF,1.D-99)
      HFLUX1=4.*PI*HSURF
      HFLUX2=4.*PI*HFLUX(NTAU-1)
      FLUXME(J)=HFLUX1/PI
      GFLUX1=4.*HSURF*WLOS(J)**2/CLIGHT
      GFLUX2=4.*HFLUX(NTAU-1)*WLOS(J)**2/CLIGHT
      FFLUX1=-2.5*log10(MAX(1.D-99,GFLUX1))
      FFLUX2=-2.5*log10(MAX(1.D-99,GFLUX2))
      spec(j,1) = wlos(j)
      spec(j,2) = hsurf
      spec(j,3) = fluxme(j)

C
C INITIATE MATRICES, TAU LOOP.
      DO 140 K=1,NTAU
C EQUATION OF RADIATIVE TRANSPORT
      IF(K.GT.1) GO TO 142
C UPPER BOUNDARY
      FKB=XK(1)/XJ(1)
      FKC=XK(2)/XJ(2)
      DC=(X(1)+X(2)+S(1)+S(2))*(TAU(2)-TAU(1))
      SQ3=SQRT(3.)
      SOURCE=(X(1)*BPLAN(1)+S(1)*XJ(1))/(X(1)+S(1))
      YA=0.
      YB=0.
      YC=0.
      DO 141 I=1,MMU_PP
      TAU1=TAU(1)*(X(1)+S(1))/XMU(I)
      EXP1=EXP(-TAU1)
      YA=YA+H(I)*XMU(I)*(1.-EXP1)
      YC=YC+H(I)*XMU(I)*EXP1*TAU1*(XPE(1)+SPE(1))/(X(1)+S(1))
141   YB=YB+H(I)*XMU(I)*EXP1*TAU1*(XT(1)+ST(1))/(X(1)+S(1))
      XG=HFLUX(1)+YA*SOURCE
      Y=DC**2*X(1)/(8.*(X(1)+S(1)))
      XJ1(K)=0.
      FH=HFLUX(1)/XJ(1)
      FG=XG/XJ(1)
      XJ2(K)=(FKC-FKB)-(.5*DC*(FG-YA*S(1)/(X(1)+S(1)))+Y)
      XJ3(K)=FKC
      XSC=X(2)+S(2)+X(1)+S(1)
      XJT2(1)=-(HFLUX(1)*.5*DC+Y*(XJ(K)-BPLAN(K))*2.)*(XT(1)+ST(1))/XSC
     * -Y*(XJ(1)-BPLAN(1))*(S(1)*XT(1)-X(1)*ST(1))/X(1)/(X(1)+S(1))
     * +(YA*0.5*DC*X(1)/(X(1)+S(1))+Y)*DBPL(1)
     * +YB*0.5*DC*SOURCE
     * -YA*.5*DC*(XJ(1)-BPLAN(1))*(S(1)*XT(1)-X(1)*ST(1))/(X(1)+S(1))**2
      XJT3(1)=-(HFLUX(1)*.5*DC+Y*(XJ(1)-BPLAN(1))*2.)/XSC*(XT(2)+ST(2))
c ??
      XJPE2(1)=-(HFLUX(1)*.5*DC+Y*(XJ(K)-BPLAN(K))*2.)*(XPE(1)+SPE(1))/
     * XSC
     * -Y*(XJ(1)-BPLAN(1))*(S(1)*XPE(1)-X(1)*SPE(1))/X(1)/(X(1)+S(1))
     * +YC*0.5*DC*SOURCE
     * -YA*.5*DC*(XJ(1)-BPLAN(1))*(S(1)*XPE(1)-X(1)*SPE(1))/(X(1)+S(1))
     * **2
      XJPE3(1)=-(HFLUX(1)*.5*DC+Y*(XJ(1)-BPLAN(1))*2.)/XSC*(XPE(2)+SPE(2
     *  ))
      GO TO 144
142   DA=DC
      FKA=FKB
      FKB=FKC
      IF(K.EQ.NTAU) GO TO 143
C INNER PARTS
      FKC=XK(K+1)/XJ(K+1)
      DC=(X(K+1)+S(K+1)+X(K)+S(K))*(TAU(K+1)-TAU(K))
      DB=DA+DC
C NOTE THAT DA AND DC ARE 2* AND DB 4* THE NORMAL DTAU'S
      A=8./(DA*DB)
      C=8./(DC*DB)
      B=-A-C
      AA=A*XK(K-1)
      BB=B*XK(K)
      CC=C*XK(K+1)
      AAA=-((AA+BB*A/(A+C))*(1.+DA/DB)+(CC+BB*C/(A+C))*DA/DB)
      CCC=-((CC+BB*C/(A+C))*(1.+DC/DB)+(AA+BB*A/(A+C))*DC/DB)
      XSA=X(K-1)+S(K-1)+X(K)+S(K)
      XSC=X(K+1)+S(K+1)+X(K)+S(K)
      ABK=X(K)/(X(K)+S(K))
      XJT1(K)=AAA/XSA*(XT(K-1)+ST(K-1))
      XJT2(K)=(AAA/XSA+CCC/XSC)*(XT(K)+ST(K))+ABK*DBPL(K)
     & +(XT(K)*S(K)-ST(K)*X(K))/(X(K)+S(K))**2*(BPLAN(K)-XJ(K))
      XJT3(K)=CCC/XSC*(XT(K+1)+ST(K+1))
c ??
      XJPE1(K)=AAA/XSA*(XPE(K-1)+SPE(K-1))
      XJPE2(K)=(AAA/XSA+CCC/XSC)*(XPE(K)+SPE(K))
     & +(XPE(K)*S(K)-SPE(K)*X(K))/(X(K)+S(K))**2*(BPLAN(K)-XJ(K))
      XJPE3(K)=CCC/XSC*(XPE(K+1)+SPE(K+1))
      XJ1(K)=FKA*8./DA/DB
      XJ2(K)=(FKA-FKB)*8./(DA*DB)+(FKC-FKB)*8./(DC*DB)-ABK
      XJ3(K)=FKC*8./DB/DC
      GO TO 144
143   CONTINUE
C LOWER BOUNDARY
      XJ1(K)=0.
      XJ2(K)=-1.
      XJ3(K)=0.
      XJT1(K)=0.
      XJT2(K)=DBPL(K)
      XJPE1(K)=0.
      XJPE2(K)=0.
      XJPE3(K)=0.
144   CONTINUE
C     
C TEMPERATURE EQUATION
      IF (K.GT.MIHAL) GO TO 145

C RADIATIVE EQUILIBRIUM
      Y=-WLSTEP(J)*X(K)
      IF(K.GT.2) then
             Y=Y*DB/(X(K)+S(K))
      end if
      if((irrinp>0) .and. (k.eq.ntau) ) then
            RT(K)=-(TT(K) - tbottom)
            TTT(K,K)=1.0
            TPE(K,K)=1.0
        else
            RT(K)=RT(K)-Y*(XJ(K)-BPLAN(K))
            TJ2(K)=Y
            TJ1(K)=0.
            TTT(K,K)=TTT(K,K)-Y*DBPL(K)+Y*(XJ(K)-BPLAN(K))*XT(K)/X(K)
            if (x(k)==0.) then
                  print*, "X is zero at k ", k
            end if
            TPE(K,K)=TPE(K,K)+Y*(XJ(K)-BPLAN(K))*XPE(K)/X(K)  
      end if
      
      GO TO 146
145   CONTINUE
C FLUX CONSTANCY

      Y=8.*WLSTEP(J)/DA
      RT(K)=RT(K)-Y*(XK(K)-XK(K-1))
      TJ1(K)=-FKA*Y
      TJ2(K)=FKB*Y
      TTT(K,K-1)=TTT(K,K-1)-Y*(XK(K)-XK(K-1))*(XT(K-1)+ST(K-1))/XSA
      TTT(K,K)=TTT(K,K)-Y*(XK(K)-XK(K-1))*(XT(K)+ST(K))/XSA
      TPE(K,K-1)=TPE(K,K-1)-Y*(XK(K)-XK(K-1))*(XPE(K-1)+SPE(K-1))/XSA
      TPE(K,K)=TPE(K,K)-Y*(XK(K)-XK(K-1))*(XPE(K)+SPE(K))/XSA

146   CONTINUE
C

C EQUATION OF RADIATION PRESSURE
      Y=PI4C*WLSTEP(J)
      RPR(K)=RPR(K)+Y*XK(K)
      PRJ(K)=-Y*FKB
      PPRG(K) = PPRG(K) + Y*HFLUX(K)*(X(K)+S(K))*ROSS(K)
C
C END OF TAU LOOP
140   CONTINUE
C
C ELIMINATE THIS WAVELENGTH
      IF(NEW.EQ.1) CALL ALGEBN(NTAU)
C
C TIME
      MS=MSA
C     CALL MSLEFT(MSA)
      MS=MS-MSA
C
C END OF WAVELENGTH LOOP
      HW1=HFLUX1*WLSTEP(J)
      FTOT=FTOT+HW1
      HW2=HFLUX2*WLSTEP(J)
      WAVEN=1.E4/WLOS(J)
      HSURF=MAX(1.D-99,HSURF)
      TRAD1=1.438E8/(WLOS(J)*
     #   log(1.0D+0+1.191D7*0.25D+0/HSURF*WAVEN**5))
      X01=log10(X(01))
      X25=log10(X(25))
      S01=log10(S(01))
      S25=log10(S(25))
30    FORMAT(1X,26F5.2)


      if (krome_debug.eq.1) then
            write(7676,'(2(999E17.8e3))') WLOS(J), BPL(steff,WLOS(J))
            Rsun_au=0.00465047
            Rstar_au= rstar*Rsun_au
            delta_omega = (Rstar_au/(semimajor))**2.0 /
     &        (4.0*(f_irrad))
            bstar_upper=BPL(steff,WLOS(J))*delta_omega
            write(7777,'(2(999E17.8e3))') WLOS(J), bstar_upper
            write(7878,'(2(999E17.8e3))') WLOS(J), XJ(1)
            write(7979,'(2(999E17.8e3))') WLOS(J), XJ(ntau)  
      endif

      if (krome_on.eq.1) then
       if (krome_photo_on.eq.1) then
            do K=1, ntau
            !calculate the radiative flux in eV cm-2 s-1 Hz-1 sr-1 from XJ with units erg s-1 cm-2 Å-1 for krome
            aa_to_cm_conv=1E-8 !converts Angstrom to centimeter
            ergs_to_eV_conv=6.242E11 !converts ergs to eV 
            FLUX_RAD(K,J)=max(0.0,XJ(K)*WLOS(J)*
     >       (WLOS(j)*aa_to_cm_conv/CLIGHT)*ergs_to_eV_conv/(4*pi))            
            enddo
       end if
      end if

150   CONTINUE
      if (krome_debug.eq.1) then
       close(7676)
       close(7777)
       close(7878)
       close(7979)
      endif
      close(960)
      ftot_a = 0.0
      do j=1, nwtot
            ftot_a = f(j) + ftot_a
      end do
      !if (model_type > 1) then 
      !print*, "Bond albedo estimate ", 1.0 - ftot_a/nwtot
      !end if
      if (newmod .eq. 2) go to 900
      
      TEFFP=TEFF*(FTOT/FLUX/PI)**.25
      print*, "Current Teff... ", teffp
      !WRITE(7,65) FTOT,TEFFP
      DO 154 K=1,NTAU
      Y=0.
      DO 155 L=1,NTAU
155   Y=Y+TTT(K,L)
      IF(PFD) WRITE(7,151) K,Y,(TTT(K,L),L=1,NTAU)
      IF(K.LE.MIHAL.AND.Y.LE.0.)TTT(K,K)=TTT(K,K)-Y
  154 CONTINUE
151   FORMAT(' TTT='/(I3,E11.4,4(/10E11.4)))
156   CONTINUE
C
C TIME
      CALL CLOCK
C
C PRINT PRESSURE EQUATION
      IF(PF) WRITE(7,48) FNORD,ITER,idmodl
      IF(PF) WRITE(7,62)
      IF(PF) WRITE(7,63)
      DO 161 K=1,NTAU
        KL=K
      ROSST(K)=(ROSST(K)-ROSS(K))/T(K)
      ROSSPE(K)=(ROSSPE(K)-ROSS(K))/PE(K)
C
C TAU SCALES
      IF(K.GT.1) GO TO 162
      CALL KAP5(TT(1),PPE(1),ABSK,1)
      TAU5(1)=TAU(1)*ABSK/ROSS(1)
      TAUP=TAU(1)/CROSS(1)
      YC=ABSK/ROSS(K)
      YD=1./CROSS(K)
      GO TO 163
162   CONTINUE
      CALL KAP5(TT(K),PPE(K),ABSK,k)
      YA=ABSK/ROSS(K)
      YB=1./CROSS(K)
      TAU5(K)=TAU5(K-1)+.5*(YA+YC)*(TAU(K)-TAU(K-1))
      TAUP=TAUP+.5*(YB+YD)*(TAU(K)-TAU(K-1))
      YC=YA
      YD=YB
163   CONTINUE
      ROSPE=ROSSPE(K)*PPE(K)/ROSS(K)
      ROST=ROSST(K)*TT(K)/ROSS(K)
      IF(PF) WRITE(7,51) K,TAU(K),PTAU(K),ROSS(K),ROSPE,ROST
     * ,YC,YD,TAU5(K),TAUP,PT(K),K
161   CONTINUE
900   CONTINUE
      CALL CLOCK
      RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C+DEF,D=MATRIX
C
C IN THIS SECTION WE COMPUTE MATRIX ELEMENTS FOR THE REST OF THE PROBLEM
C
C
      ENTRY MATRIX
C TIME
C     CALL MSLEFT(MSA)
C
C ZEROSET
      IF(PF) WRITE(7,48) FNORD,ITER,idmodl
      IF(PF) WRITE(7,57)
      IF(PF) WRITE(7,55)
      GRAD=0.
      DO 201 I=1,NTAU
      DO 201 J=1,NTAU
      VT(I,J)=0.
      VPE(I,J)=0.
      FCT(I,J)=0.
      FCPE(I,J)=0.
      DT(I,J)=0.
      DPE(I,J)=0.
      PET(I,J)=0.
      PEPE(I,J)=0.
201   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 230
      DO 231 KK=1,NTAU
      K=(1+NTAU)-KK
      IF(K.LE.NOCONV) GOTO 230
231   VV(K)=VFIX*1.E5
230   CONTINUE
C

C TAU LOOP
      DO 200 K=1,NTAU
        KL=K
C
C TERMODYNAMICAL QUANTITIES WITH PARTIAL DERIVATIVES
      K1=MAX0(1,K-1)
      TMEAN=.5*(TT(K)+TT(K1))
      PEMEAN=.5*(PPE(K)+PPE(K1))
      PRMEAN=.5*(PPR(K)+PPR(K1))
      ROSSMN=.5*(ROSS(K)+ROSS(K1))
      IF(K.GT.NOCONV) GO TO 213
C NO CONVECTION
      CALL TERMON(K,TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO,ROT,
     &                  ROPE,CP,ADIA,Q)
      RO1=RO
      PG1=PG
      PG1T=0.
      PG1PE=0.
      CPT=0.
      CPPE=0.
      ADIAT=0.
      ADIAPE=0.
      QT=0.
      QPE=0.
      GO TO 212
C CONVECTION
213   CONTINUE
      DERET=.01
      DEREP=.15
      TDELT=TMEAN*DERET
      PEDELT=PEMEAN*DEREP
      T1=TMEAN+TDELT
      PE1=PEMEAN+PEDELT
      CALL TERMON(K,T1,PEMEAN,PRMEAN,Y,YA,YB,YC,YD,YE,CP1,
     &        ADIA1,Q1)
      CALL TERMON(K,TMEAN,PE1,PRMEAN,Y,YA,YB,YC,YD,YE,CP2,
     &        ADIA2,Q2)
      CALL TERMON(K,TMEAN,PEMEAN,PRMEAN,PG1,PG1T,PG1PE,RO,
     &        ROT,ROPE,CP,ADIA,Q)
      CALL TERMON(K,TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO1,Y,
     &        YA,YB,YC,YD)
      CPPE=(CP2-CP)/PEDELT
      ADIAPE=(ADIA2-ADIA)/PEDELT
      QPE=(Q2-Q)/PEDELT
      CPT=(CP1-CP)/TDELT
      ADIAT=(ADIA1-ADIA)/TDELT
      QT=(Q1-Q)/TDELT
212   CONTINUE
      RHO(K)=RO1
C
C DEPTH SCALE
      IF(K.EQ.1) ZZ(K)=0.
      IF(K.GT.1) ZZ(K)=ZZ(K-1)-.5*DTAULN(K)*
     & (TAU(K)/ROSS(K)/RHO(K)+TAU(K-1)/ROSS(K-1)/RHO(K-1))
      IF(TAULN(K).LT.0.0) KK0=K
C
C RADIATION PRESSURE
      RPR(K)=RPR(K)-PPR(K)
CC    PRJ IS NOW FREE
CC    PRT IS ALREADY INITIATED
CC    PRPR IS UNITY
C
C TOTAL PRESSURE
      IF(K.GT.1) GO TO 202
      Y=GRAV*TAU(1)/(ROSS(1)*DLNP)
      PPTT(1)=Y*ROSST(1)/ROSS(1)
      PPPE(1)=Y*ROSSPE(1)/ROSS(1)
      GO TO 203
202   Y=GRAV*DTAULN(K)*.5
      YY=Y*TAU(K)/ROSS(K)**2
      Y=Y*TAU(K-1)/ROSS(K-1)**2
      PPPE(K)=YY*ROSSPE(K)
      PPTT(K)=YY*ROSST(K)
      PPPE(K+NTAU-1)=Y*ROSSPE(K-1)
      PPTT(K+NTAU-1)=Y*ROSST(K-1)
203   CONTINUE
C
C CONVECTION EFFICIENCY GAMMA
      HSCALE=(PG1+PRMEAN)/GRAV/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      IF(PALFA.EQ.0.) OMEGA=HSCALE*RO*ROSSMN
      Y=PY*OMEGA**2
      YY=(Y-1.)/(Y+1.)
      THETA=OMEGA/(1.+Y)
      GAMMA=-CP*RO/(8.*STEFAN*TMEAN**3*THETA)
CC    GGG IS UNITY
      GV(K)=GAMMA
      IF(PBETA.GT.0.) VV(K)=MIN(VV(K),SQRT(0.5*PP(K)/PBETA/RO))
      GG(K)=-GAMMA*VV(K)
      ROSPEM=.5*(ROSSPE(K)+ROSSPE(K1))
      ROSSTM=.5*(ROSST(K)+ROSST(K1))
      GPE(K)=-GG(K)*(CPPE/CP+ROPE/RO+YY*(ROSPEM/ROSSMN+PG1PE/PG1))
      GT(K)=-GG(K)*(CPT/CP+ROT/RO-3./TMEAN+YY*(ROSSTM/ROSSMN+PG1T/PG1))
CC    RG IS ZERO
C
C GRADIENT DIFFERENCE
      IF(K.LE.NOCONV) GO TO 206
      DELP=PP(K)-PP(K-1)
      DELT=TT(K)-TT(K-1)
      PM=PP(K)+PP(K-1)
      TM=TT(K)+TT(K-1)
      Y=1.+GG(K)
      YY=-GG(K)/Y
      GRAD=log(TT(K)/TT(K-1))/log(PP(K)/PP(K-1))
      NEWV=DD(K).GT.0..AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 263
      VV(K)=SQRT(GRAV*HSCALE*Q*PALFA**2*DD(K)/PNY)
      IF(PF) WRITE(7,64) K,VV(K)
      GO TO 203
263   CONTINUE
      NEWV=GRAD.GE.ADIA.AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 204
      VV(K)=VVMLT(GRAD-ADIA,GRAV*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      IF(PF) WRITE(7,64) K,VV(K)
      GO TO 203
204   CONTINUE
      YYY=GRAD-ADIA
C DDD IS UNITY
C ******* NEXT STATEMENT FIXES T80G4M0 BUT NOT TESTED FOR ALL MODELS
      IF(DD(K).LE.0..AND.VV(K).GT.0.) DD(K)=-YY*YYY
      RD(K)=-YY*YYY-DD(K)
      DG(K)=-YYY/Y**2
      DT(K,K)=YY*(GRAD*(1./DELT-1./TM)-.5*ADIAT)
      DT(K,K-1)=YY*(GRAD*(-1./DELT-1./TM)-.5*ADIAT)
      DP(K)=YY*GRAD*(1./PM-1./DELP)
      DP(K+NTAU-1)=YY*GRAD*(1./PM+1./DELP)
      DPE(K,K)=-.5*YY*ADIAPE
      DPE(K,K-1)=-.5*YY*ADIAPE
      GO TO 205
206   DD(K)=0.
      RD(K)=0.
      DG(K)=0.
      DP(K)=0.
      DP(K+NTAU-1)=0.
205   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 280
      Y=0.
      IF(K.GT.NOCONV) Y=-VV(K)
      GOTO 207
280   CONTINUE
C
C CONVECTIVE VELOCITY
CC    VVV IS UNITY
      Y=0.
      IF(DD(K).LE.0.) GOTO 207
      Y=-SQRT(GRAV*HSCALE*Q*PALFA**2*DD(K)/PNY)
      VD(K)=Y*.5/DD(K)
      IF(-Y.GT.2.*VV(K)) VD(K)=VD(K)*2.
      VT(K,K)=.25*Y*(QT/Q+PGT/PG-ROT/RO)
      VT(K,K-1)=VT(K,K)
      VPE(K,K)=.25*Y*(QPE/Q+PGPE/PG-ROPE/RO)
      VPE(K,K-1)=VPE(K,K)
      GO TO 208
207   VD(K)=0.
208   CONTINUE
      RV(K)=-Y-VV(K)
C
C TURBULENT PRESSURE.
CC    PTPT IS UNITY
      Y=-PBETA*VV(K)**2
      PPT(K)=-RO*Y
      PTT(K)=Y*ROT
      PTPE(K)=Y*ROPE
      PTV(K)=-PBETA*2.*VV(K)*RO
C
C CONVECTIVE FLUX
      Y=-CP*RO*PALFA*TMEAN/2./PI
CC    FCFC IS UNITY
      YY=Y*VV(K)*DD(K)
      RFC(K)=-YY-FFC(K)
      FCD(K)=Y*VV(K)
      FCV(K)=Y*DD(K)
      IF(K.LE.NOCONV) GO TO 217
      FCT(K,K)=.5*YY*(CPT/CP+ROT/RO+1./TMEAN)
      FCT(K,K-1)=FCT(K,K)
      FCPE(K,K)=.5*YY*(CPPE/CP+ROPE/RO)
      FCPE(K,K-1)=FCPE(K,K)
217   CONTINUE
C
C ELECTRON PRESSURE
209   RPE(K)=PP(K)-PG-PPR(K)
      PET(K,K)=PGT
      PEPE(K,K)=PGPE
CC    PEPR=PEPT=1.  PEP=-1.
210   CONTINUE
C
C TEMPERATURE
CC    TTT IS ALREADY INITIATED
      IF(K.GT.MIHAL) GO TO 261
C STRMGREN CONDITION
      RT(K)=RT(K)+(FFC(K+1)-FFC(K))
      GO TO 262
261   CONTINUE
C FLUXCONSTANCY
      RT(K)=RT(K)+FLUX-FFC(K)
262   CONTINUE

C END OF TAU LOOP
      PGT=PGT*TT(K)/PG
      PGPE=PGPE*PPE(K)/PG
      IF(PF) WRITE(7,51) K,TAU(K),HSCALE,ADIA,GRAD,CP,Q,PG,RO,PGPE
     * ,PGT,K
200   CONTINUE

C
C SUBTRACT CENTERED TURBULENT PRESSURE
      DO 216 K=2,NTAU
      K1=MIN0(K+1,NTAU)
216   RPE(K)=RPE(K)-.5*(PPT(K)+PPT(K1))
C
C SUBTRACT ZZ(TAU=1) FROM ZZ
      KK0 = MAX(1,KK0)
      ZZ0=ZZ(KK0)
      DO 283 K=1,NTAU
283   ZZ(K)=ZZ(K)-ZZ0
C
C TIME
      CALL CLOCK
C
C PRINT
      IF(PFE.OR.PFD) WRITE(7,48) FNORD,ITER,idmodl
      IF(PFE.OR.PFD) WRITE(7,54)
      IF(PFE.OR.PFD) WRITE(7,52)
      IF(PFE.OR.PFD)WRITE(7,51)(I,TAU(I),RPR(I),PPT(I),RP(I),GG(I),RD(I)
     * ,RV(I),RFC(I),RPE(I),RT(I),I,I=1,NTAU)
C
C SAVE DD-MATRICES
      DO 270 K=1,NTAU
      DPS(K)=DP(K)
      DPES(K)=DPE(K,K)
      DTS(K)=DT(K,K)
      D(K)=RD(K)
      IF(K.EQ.1) GO TO 270
      DPS(K+NTAU-1)=DP(K+NTAU-1)
      DPES(K+NTAU-1)=DPE(K,K-1)
      DTS(K+NTAU-1)=DT(K,K-1)
270   CONTINUE
C
C+DEF,D=ELIMIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GAUSS ELIMINATION TO UPPER TRIANGULAR FORM.
C
C TIME
      CALL CLOCK
C
C RADIATION PRESSURE
      RP(1)=RP(1)+RPR(1)
      DO 320 I=2,NTAU
      RPE(I)=RPE(I)-RPR(I)
      DO 301 J=1,NTAU
      PET(I,J)=PET(I,J)-PRT(I,J)
301   CONTINUE
C
C TURBULENT PRESSURE
      PET(I,I)=PET(I,I)-PTT(I)
      PEPE(I,I)=PEPE(I,I)-PTPE(I)
310   CONTINUE
C
C CONVECTIVE EFFICIENCY
      DT(I,I)=DT(I,I)-.5*DG(I)*GT(I)
      DT(I,I-1)=DT(I,I-1)-.5*DG(I)*GT(I)
      DPE(I,I)=DPE(I,I)-.5*DG(I)*GPE(I)
      DPE(I,I-1)=DPE(I,I-1)-.5*DG(I)*GPE(I)
      DV(I)=-DG(I)*GV(I)
C
C TOTAL PRESSURE
      DPI=-DP(I+NTAU-1)
      DPE(I,I-1)=DPE(I,I-1)-DPI*PPPE(I+NTAU-1)
      DT(I,I-1)=DT(I,I-1)-DPI*PPTT(I+NTAU-1)
      DPE(I,I)=DPE(I,I)-DPI*PPPE(I)
      DT(I,I)=DT(I,I)-DPI*PPTT(I)
      RD(I)=RD(I)-DPI*RP(I)
      DP(I)=DP(I)-DPI
320   CONTINUE
C
      RPI=0.
      DV(1)=0.
      DO 321 I=1,NTAU
C WE HAVE THE MATRICES       PPP       PPPE  PPTT
C AND                        PEPP      PEPE  PET
C WE WANT TO SUBTRACT FROM PEPE AND PET THE PRODUCTS OF PEPP*PPP-INVERS
C PPPE AND PET RESPECTIVELY. PPP IS BIDIAGONAL WITH UNITY ON THE DIAGONA
C MINUS UNITY ON THE SUBDIAGONAL. ITS INVERS IS A MATRIX WITH UNITY EVER
C UNDER AND ON THE DIAGONAL. PEPP IS MINUS UNITY. THUS THE PRODUCT OF PE
C INVERS WITH PPPE IS MATRIX OF THE FOLLOWING TYPE. IN EVERY COLUMN EACH
C ELEMENT IS THE SUM OF ALL ELEMENTS ABOVE THAT POINT (AND INCLUDING) IN
C PPPE MATRIX. SIMILARILY FOR THE PPTT-MATRIX.
      RPI=RPI+RP(I)
      RPE(I)=RPE(I)+RPI
      RD(I)=RD(I)-DP(I)*RPI
      Y=PPPE(I)
      YY=PPTT(I)
      YY=YY+PRT(1,I)
      JMIN=MAX0(I,1)
      DO 321 J=JMIN,NTAU
      IF(J.NE.I+1) GO TO 322
      Y=Y+PPPE(I+NTAU)
      YY=YY+PPTT(I+NTAU)
322   CONTINUE
      

      PEPE(J,I)=PEPE(J,I)+Y

      PET(J,I)=PET(J,I)+YY
      DPE(J,I)=DPE(J,I)-DP(J)*Y
      DT(J,I)=DT(J,I)-DP(J)*YY
321   CONTINUE
C
C GRADIENT DIFFERENCE
      DO 350 I=2,NTAU
      DVI=DV(I)
      VDI=VD(I)
      FCDI=FCD(I)
      RFC(I)=RFC(I)-FCDI*RD(I)
      RV(I)=RV(I)-VDI*RD(I)
      VVVI=MAX(0.5d+0,1.0D+0-VDI*DVI)
      FCV(I)=FCV(I)-FCDI*DVI
      DO 340 J=1,NTAU
      VT(I,J)=VT(I,J)-VDI*DT(I,J)
      VPE(I,J)=VPE(I,J)-VDI*DPE(I,J)
      FCT(I,J)=FCT(I,J)-FCDI*DT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCDI*DPE(I,J)
340   CONTINUE
C
C TURBULENT VELOCITY
      PEVI=-PTV(I)
      RV(I)=RV(I)/VVVI
      FCVI=FCV(I)
      RFC(I)=RFC(I)-FCVI*RV(I)
      RPE(I)=RPE(I)-PEVI*RV(I)
      DO 350 J=1,NTAU
      VT(I,J)=VT(I,J)/VVVI
      VPE(I,J)=VPE(I,J)/VVVI
      FCT(I,J)=FCT(I,J)-FCVI*VT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCVI*VPE(I,J)
      PET(I,J)=PET(I,J)-PEVI*VT(I,J)
      PEPE(I,J)=PEPE(I,J)-PEVI*VPE(I,J)
350   CONTINUE
C
C TIME
      CALL CLOCK
C
C CONVECTIVE FLUX
      !print*, "convective flux "
      DO 360 I=1,NTAU
      RT(I)=RT(I)-RFC(I)
      DO 361 J=1,NTAU
      !if (j==2) print*, TTT(I,J)
      TTT(I,J)=TTT(I,J)-FCT(I,J)
      !if (j==2) print*, TTT(I,J)
      TPE(I,J)=TPE(I,J)-FCPE(I,J)
361   CONTINUE
      IF(I.GT.MIHAL) GO TO 360
      RT(I)=RT(I)+RFC(I+1)
      DO 362 J=1,NTAU
      TTT(I,J)=TTT(I,J)+FCT(I+1,J)
      TPE(I,J)=TPE(I,J)+FCPE(I+1,J)
362   CONTINUE
360   CONTINUE
C
C ELECTRON PRESSURE
      
      CALL MATINV(PEPE,NTAU)
      !print*, "PET before MULT ", pet(1,1)
      !print*, "PEPE before MULT ", PEPE(1,1)
      !print*, "PET before MULT ", PET(1,1)
      !print*, "SCRATC before MULT ", SCRATC(1,1)
      CALL MULT(PET,PEPE,PET,SCRATC,NTAU,NTAU)
      !print*, "pet after ", pet(1,1)
      CALL MULT(RPE,PEPE,RPE,SCRATC,NTAU,1)
      !print*, "electron pressure "
      DO 374 I=1,NTAU
      DO 374 J=1,NTAU
      SUMA=0.
      DO 375 L=1,NTAU
      !print*, "tpe(i,l) ", tpe(i,l)
      !print*, "pet(l,j) ", pet(l,j)
      SUMA=SUMA+TPE(I,L)*PET(L,J)
      !if (j==2) print*, SUMA
375   CONTINUE
      
      TTT(I,J)=TTT(I,J)-SUMA
      ! print*, "i ", i 
      ! print*, "j ", 
      RT(I)=RT(I)-TPE(I,J)*RPE(J)
374   CONTINUE
C
C TIME
      CALL CLOCK
C      
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C BACKSUBSTITUTION IN GAUSS ELIMINATION SCHEME.

C
C INITIATE
      CALL MATINV(TTT,NTAU)
C
      DO 400 I=1,NTAU
      PE(I)=RPE(I)
      FC(I)=RFC(I)
      V(I)=RV(I)
      PT(I)=0.
      PR(I)=RPR(I)
      P(I)=RP(I)
C
C
C SOLVE FOR TEMPERATURE CORRECTION
      T(I)=0.
      TTTS(I)=0.
      RTS(I)=0.
      DO 399 J=1,NTAU
      T(I)=T(I)+TTT(I,J)*RT(J)
      
   
399   CONTINUE
      ! print*, "i ", i
      ! print*, "t corr",t(i) 
400   CONTINUE

C
C
999   FORMAT(I3,1P5E10.3)
C
C---

C If the maximum temperature correction the program suggest to any model layer
C is smaller than TCONV (from input), then we will finish with computing one 
C more iteration where we will write the output etc. 
C UGJ 25/11/2004:
C In order to save time in interpolating in two-dimensional (i.e. atomic) 
C absorption coefficient files and in order not to compute a new atomic
C absorption coefficient file for each iteration, we start with a one-dimesional
C atomic file for a model with near-by parameters, and then recompute a new 
C (possible temporary) 1D atomic file for the almost converged model and use that 
C in the final iteration (if this cause a total sum of TCORMX of the following 
C iterations to be larger than 20.*TCONV then we repeat and make one more 
C atomic file for the more coverged model). In this way the opacity can be 
C calculated in the same way for the molecular (Doppler broadening only, and
C hence 1D absorption coefficients, dependent on temperature only) and the 
C atomic (T,Pg dependent Vogit profiles) absorption coefficient.

      TCORMX=ABS(T(1))
      DO 405 I=2,NTAU
      PM = T(I)/ABS(T(I))
      TCORMXM=PM*MAX(ABS(TCORMX),ABS(T(I)))
 405  TCORMX=MAX(TCORMX,ABS(T(I)))
         if(newosatom.eq.2) then 
                     tcormxend = tcormxend + tcormx
         end if
      IF(TCORMX.LE. 5.*TCONV )  THEN 
         if(newosatom.eq.1) then 
C                    call osatom
                     newosatom = 2
                     tcormxend = 0.
         end if
      END IF
      IF(TCORMX.LE. TCONV )  THEN 
         
         ITSTOP=.TRUE.
C        if(newosatom.eq.2 .and. tcormxend.gt.20.0*tconv) call osatom
      END IF
      IF(TCORMX.GT.TCONV )  THEN !catch cases where in the onemor iteration TCORMX can become larger than TCONV
               ITSTOP=.FALSE.
      ENDIF
      PRINT406, TCORMXM,IT
      inquire(file="tcormx.dat", exist=exist)
      if (exist) then
      open(unit=987,file= 'tcormx.dat',status='replace')
      write(987, *) TCORMXM
      close(987)
      else 
      open(unit=987,file= 'tcormx.dat',status='new')
      write(987, *) TCORMXM
      close(987)
      end if
      
406   FORMAT(' Max corr. to T wanted was ',F8.1,' K for iteration ',I3)
C
C
C     CALL DISPLA(' TCORMX ',TCORMX)
C---
C
C CHECK T CORR
      IF(KORT.EQ.1) THEN
      DO 401 I=1,NTAU
      T(I)=T(I)/SQRT(1.+25.*(T(I)/TT(I))**2)
401   CONTINUE
      TCORMX=ABS(T(1))
      DO 407 I=2,NTAU
      PM = T(I)/ABS(T(I))
 407  TCORMX=MAX(TCORMX,ABS(T(I)))
      
      
      PRINT4061, TCORMX
4061  FORMAT(' Max corr. to T wanted for kort=1 was',F6.1)
      END IF
C
      IF (KORT.EQ.4.AND.IT.GE.4  .OR. KORT.EQ.5) THEN
      PPK = TDIFF
      DO 4011 I=1,NTAU
      if (abs(T(I)).GT.abs(PPK))
     - T(I)=PPK*t(i)/abs(t(i))
4011  CONTINUE
C
      TCORMX=ABS(T(1))
      DO 4071 I=2,NTAU
      PM = T(I)/ABS(T(I))
      TCORMX=MAX(TCORMX,ABS(T(I)))
4071  TCORMXM=PM*TCORMX
      
      !print*, "maximum temperature correction LINE 9129"
      write(file_id,'(i3)') iter
      file_name = 'PT_struct_' // trim(adjustl(file_id)) // '.dat'
      !print*, file_name
      open(unit=73972, file=file_name, status='new')
      do i=1, ntau
      write(73972,*) 
     * TT(i), ",", PP(i), ",", PPE(i), "," ,ROSS(i)
      end do
      close(73972)
      PRINT408, TCORMX,KORT
408   FORMAT(' Maximum correction applied was',F6.1,
     *     ' for applied kort =',I2)
      END IF
      
      IF (KORT.EQ.6) THEN
      IF(IT.EQ.1) TCORPRV = 1.D3
      IF(IT.EQ.1) TCORMIN = 1.D3
      PPK = TDIFF
      PPKK = 2.D0*TCORPRV
      PPK = MIN(PPKK,PPK)
      PPPK = MIN(TCORMIN,PPK)
      TCORMX=ABS(T(1))
      DO 4081 I=2,NTAU
      PM = T(I)/ABS(T(I))
      TCORMX=MAX(TCORMX,ABS(T(I)))
4081  TCORMXM=PM*TCORMX
        TCORPRV=TCORMX
        TCORMIN=min(TCORMX,tcormin)
C       IF (TCORMX.GT.ppk) THEN
        IF (TCORMX.GT.pppk) THEN
        DO 4082 I=1,NTAU
        if (abs(T(I)).GT.PPK) then
C         T(I)=PPK*t(i)/abs(t(i))
          T(I)=PPPK*t(i)/abs(t(i))
          TCORPRV=ABS(T(I))
        end if
4082    CONTINUE
        END IF

      PRINT4083, TCORPRV,KORT
4083  FORMAT(' Maximum correction applied was',F6.1,
     *     ' for applied kort =',I2)
      END IF

C
C SUBTRACT TEMPERATURE
      DO 410 I=1,NTAU
      PT(I)=PT(I)-PTT(I)*T(I)
      P(I)=P(I)-PPTT(I)*T(I)
      IF(I.GT.1) P(I)=P(I)-PPTT(I+NTAU-1)*T(I-1)
      DO 410 J=1,NTAU
      V(I)=V(I)-VT(I,J)*T(J)
      PE(I)=PE(I)-PET(I,J)*T(J)
      FC(I)=FC(I)-FCT(I,J)*T(J)
      PR(I)=PR(I)-PRT(I,J)*T(J)
410   CONTINUE
C
C SUBTRACT ELECTRON PRESSURE
      DO 420 I=1,NTAU
      PT(I)=PT(I)-PTPE(I)*PE(I)
      P(I)=P(I)-PPPE(I)*PE(I)
      IF(I.GT.1) P(I)=P(I)-PPPE(I+NTAU-1)*PE(I-1)
      DO 420 J=1,NTAU
      V(I)=V(I)-VPE(I,J)*PE(J)
      FC(I)=FC(I)-FCPE(I,J)*PE(J)
420   CONTINUE
C
C TOTAL PRESSURE
      P(1)=P(1)+PR(1)-RPR(1)
      DO 425 I=2,NTAU
      P(I)=P(I)+P(I-1)
      !print*, "i ", i , " p ", p(i)
425   CONTINUE
C
C SUBTRACT VELOCITY
      DO 430 I=1,NTAU
      PT(I)=PT(I)-PTV(I)*V(I)
430   CONTINUE
C
C PRINT
      ! IF(PFE.OR.PFD) WRITE(7,48) FNORD,ITER,idmodl
      ! IF(PFE.OR.PFD) WRITE(7,56)
      ! IF(PFE.OR.PFD) WRITE(7,52)
      DO 440 I=1,NTAU
C
C SOLVE FOR CONVECTIVE EFFICIENCY AND GRADIENT DIFFERENCE
      I1=MAX0(I-1,1)
      GI=-.5*(GT(I)*(T(I)+T(I1))+GPE(I)*(PE(I)+PE(I1)))-GV(I)*V(I)
      GG(I)=GG(I)+GI
      D(I)=D(I)-DG(I)*GI-DPS(I)*P(I)-DPES(I)*PE(I)-DTS(I)*T(I)
      IF(I.GT.1) D(I)=D(I)-DPS(I+NTAU-1)*P(I-1)-DPES(I+NTAU-1)
     &*PE(I-1)-DTS(I+NTAU-1)*T(I-1)
c       IF(PFE.OR.PFD)WRITE(7,51)I,TAU(I),PR(I),PT(I),P(I),GI,D(I),V(I)
c      * ,FC(I),PE(I),T(I),I
440   CONTINUE
C
C TIME
      CALL CLOCK
C
C APPLY CORRECTIONS
      IPRESS=0
      DO 450 I=1,NTAU
      TT(I)=TT(I)+T(I)
      VV(I)=MAX(VV(I)+V(I),0.0D+0)
      FFC(I)=FFC(I)+FC(I)
      PPT(I)=MAX(PPT(I)+PT(I),0.0D+0)
      PPT(I)=MIN(PPT(I),0.5*PP(I))
      PPR(I)=PPR(I)+PR(I)
      DD(I)=DD(I)+D(I)
C
C IF TOO VIOLENT CHANGES TO PPE OR PP, SET IPRESS FOR AN EXTRA
C PRESSURE INTEGRATION AFTER CORRECTIONS HAVE BEEN APPLIED.
      IF(ABS(P(I)/PP(I)).LT.0.5.AND.ABS(PE(I)/PPE(I)).LT.0.5) GOTO 451
      IPRESS=1
      GOTO 450
C
451   PPE(I)=PPE(I)+PE(I)
      PP(I)=PP(I)+P(I)
450   CONTINUE
      IF(IPRESS.EQ.1) CALL TRYCK

C
C PRINT PRESENT STATE OF THE ATMOSPHERE
      ENTRY PRESNT
      IF(.NOT.PFE) RETURN
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
      WRITE(7,48) FNORD,ITER,idmodl
      WRITE(7,53)
      WRITE(7,52)
      I=0
      Y=0.
      TT0=(TAU(2)*TT(1)-TAU(1)*TT(2))/(TAU(2)-TAU(1))
      WRITE(7,51) I,Y,PPR(1),Y,PPR(1),Y,Y,Y,Y,Y,TT0,I
      WRITE(7,51) (I,TAU(I),PPR(I),PPT(I),PP(I),GG(I),DD(I)
     * ,VV(I),FFC(I),PPE(I),TT(I),I,I=1,NTAU)
C
C TIME
      CALL CLOCK
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1MARCS',F5.1,5X,'SOLVE(61) 12-NOV-80',5X,
     & '                   ',5X,'ITERATION',I3,5X,A24/)
49    FORMAT(13('1234567890'),'123')
50    FORMAT(4(/2X,1P10E13.5))
51    FORMAT(I3,1P10E12.4,I4)
52    FORMAT(7X,'TAU',9X,'PRAD',8X,'PTURB',7X,'PTOT',8X,'GAMMA',
     & 7X,'DELTA',7X,'VCONV',7X,'FCONV',7X,'PE',10X,'TEMP')
53    FORMAT(' STATE OF MODEL ATMOSPHERE')
54    FORMAT(' *=R.H. SIDES     *',23X,'*',23X,'*',
     ,11X,'*',11X,'*',11X,'*',11X,'*')
55    FORMAT(6X,'TAU',9X,'HSCALE',6X,'ADIA',8X,'GRAD',8X,'CP',10X,
     *'Q',11X,'PG',10X,'RO',10X,'PGPE',8X,'PGT')
56    FORMAT(' CORRECTIONS')
57    FORMAT(' THERMODYNAMICALS')
58    FORMAT(1X,I3,F9.1,F8.1,1P2E10.3,0PF6.3,1PE10.3,0PF8.3,
     * F7.0,4F5.1,1PE10.3,2(0PF6.1),2I5,I4)
59    FORMAT('   J',4X,'WAVEL',2X,'WEIGHT',2X,'F(WAVEL)',3X,
     * 'CONTRIB',1X,'WAVEN',2X,'F(WAVEN)',4X,'MAGN',2X,'TRAD',3X,
     * 'X01',2X,'X25',2X,'S01',2X,'S25',3X,'CONTRIB',2X,
     * 'DX01  DX25 OPAC TRAN ALGB')
60    FORMAT(' SAVED VALUES FROM LOGICAL UNIT',I3)
61    FORMAT(' STARTING VALUES, QTEMP=',F5.2)
62    FORMAT(' PRESSURE EQUATION')
63    FORMAT(6X,'TAU',9X,'PTAU',8X,'ROSS',8X,'ROSSPE',6X,'ROSST',7X,
     &'KAP5',8X,'ROSSP',7X,'TAU5',8X,'TAUP')
64    FORMAT(' K=',I2,4X,'NEW V =',1PE10.3)
65    FORMAT(' TOTAL FLUX=',1PE11.4,' ERGS/CM**2/S  TEFF=',0PF6.0,' K')
66    FORMAT(7X,F8.0)
67    FORMAT(1X,10F7.1)
68    FORMAT(' ITERATION',I3)
      
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRANEQ
      implicit real*8 (a-h,o-z)
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU_PP ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. AITKEN EXTRAPOLATION ACCELLERATES THE CONVERGENCE.
C 3. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 4. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 5. THE SCATTERING IS TREATED AS DIPOLE SCATTERING INSTEAD OF THE NORMA
C    USED ISOTROPIC APPROXIMATION. THIS CAN BE DONE VERY SIMPLY IN THE
C    ITERATING CANNON SCHEME.
C 6. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQS.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CANGLE/XMU(6),XMU2(6),H(6),MMU_PP
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,dumtran(4*ndp),idumtran(3)
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),DUM(3*NDP),P(NDP)
     & ,SP1(NDP,6),SP2(NDP,6),SP3(NDP,6),AN(NDP),AD(NDP),BD(NDP)
     & ,FACT(NDP),DSO(NDP),SP2DUM((4*NDP-29)*NDP)
      DIMENSION A(7)
      COMMON /CSURF/HSURF,Y1(NRAYS)
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom  
      common /starspec/ stellar_spectrum(nwreal),index_wlambda

C
C INITIATE
      
      if (irrin>0) then
        if (input_star_spec==1) then
        bstar = spectrum_scale*stellar_spectrum(index_wlambda)
        else
        bstar = bpl(steff, wlambda)
        endif
      end if
      DO  K=1,JTAU
      FACT(K)=1.
      DSO(K)=0.
      XJ(K)=0.
      XK(K)=0.
      ERROR(K)=BPLAN(K)*X(K)/(X(K)+S(K))
      SOURCE(K)=0.

      end do
C
C CALCULATE THE MATRIX ELEMENTS
      CALL TRANFR
      CALL TRANSC
C
C ITERATION LOOP
      ITMAX=7
      DO 110 IT=1,ITMAX
110   A(IT)=0.
      DO 140 IT=1,ITMAX
      ITM=IT
C
C SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
      h_irrad = 1.0
      CALL SCATTR

      DO 120 K=1,JTAU
      if (irrin>0) then 
       XJ(K)=XJ(K)+P(K) + Pstar(K)
      else
       XJ(K)=XJ(K)+P(K)
      end if 

      if (irrin>0) then
      XK(K)=XK(K)+.333333*(P(K)+Pstar(k))
      else
      XK(K)=XK(K)+.333333*P(K)
      end if
C
C AITKEN EXTRAPOLATION USED FOR CONVERGENCE ACCELLERATION
      if (irrin>0) then
      DS=ERROR(K)+(P(K)+Pstar(k))*S(K)/(X(K)+S(K))
      else
      DS=ERROR(K)+P(K)*S(K)/(X(K)+S(K))
      end if
      IF(DSO(K).NE.0.) 
     #   FACT(K)=MIN(1.25D+0,MAX(0.8D+0,FACT(K)-DS/DSO(K)))
      DS=DS/FACT(K)
      IF(IT.GE.2) DSO(K)=DS
120   SOURCE(K)=SOURCE(K)+DS

      !different call for taking into account the gaussian weights 
      !in the irradiation routine
      
C
C SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION

      CALL FORMAL

C
C CHECK ERROR IN SOURCE FUNCTION
      DO 130 K=1,JTAU
130   A(IT)=MAX(A(IT),ABS(ERROR(K)/SOURCE(K)))
C
C END OF ITERATION LOOP
      IF(A(IT).LT.0.0001) GO TO 141
140   CONTINUE
      
50    FORMAT(' MAXFEL =',12F10.7)
141   CONTINUE
C     
      if (wlambda == wlos(1)) then
            open(unit=51298, file='flux_toa.dat', status='replace')
            write(51298, '(2e24.15)') wlambda, HSURF
            close(51298)
      else
            open(unit=51298, file='flux_toa.dat', status="old", 
     *       position="append", action="write")
            write(51298,'(2e24.15)') wlambda, HSURF
            close(51298)
      end if
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C
      SUBROUTINE TRANFR
      implicit real*8 (a-h,o-z)
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION SOURCE
C ERROR IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. FLUX AND INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /CANGLE/XMU(6),XMU2(6),H(6),MMU_PP
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,dumtran(4*ndp),idumtran(3)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),DUM(3*NDP),P(NDP)
     *,SP1(NDP,6),SP2(NDP,6),SP3(NDP,6),AN(NDP),AD(NDP),BD(NDP)
     *,FACT(NDP),DSO(NDP),C(6),T(6),EX(6),SP2DUM((4*NDP-29)*NDP-18)
      COMMON /CSURF/HSURF,Y1(NRAYS)
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom    
      real*8, parameter :: pi = 3.14159265
C
C MU LOOP
      JTAU1=JTAU-1
      JTAU2=JTAU-2
      IF (mmu_pp.GT.nrays) STOP
     &    ' tranfr:increase nrays in parameter.inc'
      IF (mmu_pp.GT.6) STOP ' tranfr: increase dimensions incl mmu_pp'
      
      DO 110 I=1,MMU_PP
      
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU(I)
      A=1./DTAUB
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
C LET P BE THE EVEN PART OF THE INTENSITY, THEN
C
C         P(2)= P(1) + D*P'(1) + .5*D2*P''(1)
C OR      P(2)= P(1) + D*(P(1)-I(1,-MU)) + .5*D2*(P(1)-S(1)) .
C WHERE   I(1,-MU) = S(1)*(1.-EXP(-T))
C
C THE DIFFERENCE AS COMPARED TO THE USUAL SECOND ORDER BOUNDARY CONDITIO
C IS THE ADDITIONAL TERM   I(1,-MU)=S(1)*(1.-EXP(-T)). THUS THE COEFFICI
C FOR S(1) IN THE FIRST EQUATION SHOULD BE CHANGED AS FOLLOWS
C         S(1)=S(1)*(1.+C*(1.-EXP(-T))
C WHERE   C=2./D
C *NORD* 751009
      C(I)=2.*A
      T(I)=TAU(1)*(X(1)+S(1))/XMU(I)
      IF(T(I).LE.0.1) THEN
        EX(I)=T(I)*(1.-.5*T(I)*(1.-.3333*T(I)))
      ELSE IF(T(I).LE.675.0) THEN
        EX(I)=1.0-EXP(-T(I))
      ELSE IF(T(I).GT.675.0) THEN
        EX(I)=1.0
      END IF
C
C K=2,JTAU-1
      DO 100 K=2,JTAU1
      DTAUA=DTAUB
      DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU(I)
      DTAUC=.5*(DTAUA+DTAUB)
      AD(K)=.166667*DTAUA/DTAUC
      BD(K)=.166667*DTAUB/DTAUC
      SP1(K,I)=-1./(DTAUA*DTAUC)+AD(K)
      SP2(K,I)=1.
100   SP3(K,I)=-1./(DTAUB*DTAUC)+BD(K)
C
C K=JTAU
      SP2(JTAU,I)=1.
C
C END OF MU LOOP
110   CONTINUE
C
C ELIMINATE SUBDIAGONAL, SAVE FACTORS IN SP1
      DO 121 I=1,MMU_PP
      DO 120 K=1,JTAU2
      SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
      SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
120   SP2(K,I)=SP2(K,I)-SP3(K,I)
121   SP2(JTAU-1,I)=SP2(JTAU-1,I)-SP3(JTAU-1,I)
C
      RETURN
C
      ENTRY FORMAL
C
C ZEROSET
      DO 130 K=1,JTAU
      AN(K)=(3.*XK(K)-XJ(K))/8.*S(K)/(X(K)+S(K))
      XK(K)=0.
130   XJ(K)=0.
C
C MU LOOP
      XH(1)=0.
      HSURF=0.
      DO 170 I=1,MMU_PP
      
C
C INITIATE APPROXIMATIVE SOURCE FUNCTION
      P(1)=SOURCE(1)+AN(1)*(3.*XMU2(I)-1.)
C NOTE THE ANISOTROPIC SCATTERING CORRECTION
      S0=P(1)
      if (irrin>0) then
      Rsun_au=0.00465047
      Rstar_au= rstar*Rsun_au
      delta_omega = (Rstar_au/(semimajor))**2.0 /
     &  (4.0*(f_irrad))
      P_star=(delta_omega/H(i)) * bstar * 
     &  (C(i)*(1-EX(i)))
      P(1)=P(1) *(1.+C(I)*EX(I))*(1.-delta_omega)
      P(1) = P(1)+P_star
      else
      P(1)=P(1) *(1.+C(I)*EX(I))
      end if

      DO K=2,JTAU1
      P(K)=(1.-AD(K)-BD(K))*(SOURCE(K)+AN(K)*(3.*XMU2(I)-1.))
     & +AD(K)*(SOURCE(K-1)+AN(K-1)*(3.*XMU2(I)-1.))
     & +BD(K)*(SOURCE(K+1)+AN(K+1)*(3.*XMU2(I)-1.))
      end do 
      
      P(JTAU)=SOURCE(JTAU)
C
C ACCUMULATE RIGHT HAND SIDE
      DO 150 K=1,JTAU2
150   P(K+1)=P(K+1)+SP1(K,I)*P(K)
C
      
C BACKSUBSTITUTE
      DO 160 K=1,JTAU1
      P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K,I)*P(JTAU-K+1))/SP2(JTAU-K,I)
      XK(JTAU-K)=XK(JTAU-K)+H(I)*P(JTAU-K)*XMU2(I)
160   XJ(JTAU-K)=XJ(JTAU-K)+H(I)*P(JTAU-K)
C     
C END OF MU LOOP
      XK(JTAU)=XK(JTAU)+H(I)*P(JTAU)*XMU2(I)
      R1=P(1)-S0*EX(I)
      XH(1)=XH(1)+H(I)*XMU(I)*R1
      P0=P(1)*(1.-EX(I))+.5*S0*EX(I)**2
      HSURF=HSURF+H(I)*XMU(I)*P0
      Y1(I)=2.*P0
C HSURF AND Y1(6) ARE THE FLUX AND INTENSITIES AT THE SURFACE
170   CONTINUE
      XJ(JTAU)=P(JTAU)
C
C 'XJ' IS THE NEW MEAN INTENSITY
      DO 180 K=1,JTAU
180   ERROR(K)=(X(K)*BPLAN(K)+S(K)*XJ(K))/(X(K)+S(K))-SOURCE(K)
C
C FLUX AND SECOND MOMENT
      DO 190 K=2,JTAU
190   XH(K)=2.*(XK(K)-XK(K-1))/(X(K)+S(K)+X(K-1)+S(K-1))/
     /(TAU(K)-TAU(K-1))
C
      RETURN

      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRANSC
      implicit real*8 (a-h,o-z)
C
C SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
C IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE MU POINT.
C 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' CONTAINS
C THE ESTIMATED MEAN INTENSITY ON EXIT. TRANSC CALCULATES THE MATRIX
C ELEMENTS FOR SCATTR.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,dumtran(4*ndp),idumtran(3)
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2_PP/SOURCE(NDP),ERROR(NDP),SP1(NDP),SP2(NDP),
     +            SP3(NDP),P(NDP),SP2DUM((4*NDP-6)*NDP)
      common /cirinp/steff,reflect,f_irrad,h_irrad,
     > wlambda,bstar,spectrum_scale,irrinp,irrin,input_star_spec
      common /irradcs/Pstar(ndp),rstar, semimajor,tbottom
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     * STEFAN
      DATA XMU,XMU2/0.5773503,0.3333333/
C
C K=1
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU
      A=1./DTAUB
      B=A**2
      SP2(1)=2.*A+X(1)/(S(1)+X(1))
      SP3(1)=-2.*B
      C=2.*A
      T=TAU(1)*(X(1)+S(1))/XMU
      EX=T*(1.-.5*T*(1.-.33333*T))
      IF(T.GT.0.1) EX=1.-EXP(-T)
      SP2(1)=SP2(1)-C*S(1)/(X(1)+S(1))*EX
C
C K=2,JTAU-1
      JTAU1=JTAU-1
      DO 100 K=2,JTAU1
      DTAUA=DTAUB
      DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU
      DTAUC=.5*(DTAUA+DTAUB)
      A=1./(DTAUA*DTAUC)
      B=1./(DTAUB*DTAUC)
      SP1(K)=-A
      SP2(K)=X(K)/(S(K)+X(K))
      SP3(K)=-B
100   CONTINUE
C
C K=JTAU
      SP2(JTAU)=X(JTAU)/(X(JTAU)+S(JTAU))
C
C ELIMINATE SUBDIAGONAL
      JTAU2=JTAU-2
      DO 110K=1,JTAU2
      SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
      SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
110   SP2(K)=SP2(K)-SP3(K)
      SP2(JTAU-1)=SP2(JTAU-1)-SP3(JTAU-1)
      RETURN
C
      ENTRY SCATTR
C
C INITIATE INHOMOGENOUS TERMS
      DO 120 K=1,JTAU
120   P(K)=ERROR(K)
C PRELIM
      if (irrin<=0) then
            P(1)=P(1)*(1.+C*EX)
            
      else if (irrin>0) then
      !IRRADIATION
      Rsun_au=0.00465047
      Rstar_au= rstar*Rsun_au
      delta_omega = (Rstar_au/(semimajor))**2.0 / 
     & (4.0*(f_irrad))
      Pstar(1) = (C*(1-EX))* 
     & (delta_omega/h_irrad) *bstar
      P(1)=P(1) *(1.+C*EX)*(1.-delta_omega) 
      end if
C
C ACCUMULATE INHOMOGENOUS TERMS
      DO 130 K=1,JTAU2
130   P(K+1)=P(K+1)+SP1(K)*P(K)
C
      if (irrin>0) then
      do k=1, JTAU1
       Pstar(k+1) = SP1(k)*Pstar(k)
      end do
      if (irrinp >0) then
       Pstar(JTAU) = Pstar(jtau) + reflect*Pstar(jtau)
      end if 
      end if

C BACKSUBSTITUTE
      P(JTAU)=P(JTAU)/SP2(JTAU)
      Pstar(jtau)=Pstar(jtau)/SP2(JTAU)
      DO 140 K=1,JTAU1
      Pstar(JTAU-K)=(Pstar(JTAU-K)
     & -SP3(JTAU-K)*Pstar(JTAU-K+1))/SP2(JTAU-K)
140   P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K)*P(JTAU-K+1))/SP2(JTAU-K)
C
      RETURN
      END
C*
C*NEW PDS MEMBER FOLLOWS
C*
      SUBROUTINE TRYCK
      implicit real*8 (a-h,o-z)
C
C TRYCK IS A FAST PRESSURE INTEGRATION ROUITINE. IT IS FAST BECAUSE OF
C TWO REASONS: 1) IT INTEGRATES THE DIFFFERENTIAL EQUATION FOR LN(P)
C AS A FUNCTION OF LN(TAU). 2) IT ITERATES DIRECTLY ON THE ELECTRON
C PRESSURE, KEEPING THE NUMBER OF CALLS TO ABSKO TO A MINIMUM.
C ASSUMING A POWER LAW BEHAVIOUR OF PP,TT,PPE,ETC.: PP=C*TAU**DP,ETC., O
C CAN SHOW THAT DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST)/(1.+ROSSPE/PGPE).
C THE ANSATZ FOR PP IMPLIES PP(1)=TAU(1)*GRAV/(ROSS(1)*DP), WHICH
C SERVES AS A BOUNDARY CONDITION.
C 790516 *NORD*
C
C TRYCK/KOL IS A VERSION THAT USES CROSS(K) AS A FACTOR ON ROSSOP.
C 801105 *NORD*
C
      include 'parameter.inc'
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),
     & NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CROSSOS/ ROSSO(NDP),PTAUO(NDP)
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE 
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./
CUGJ98DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./

C
C START
C     CALL MSLEFT(MSA)
      MSA=0
      EPS = 1.0e-3;
      DT=0.
C USE 'DLNT/DLNTAU'=DT=0. TO BE COMPATIBLE WITH SOLVE. OTHERWISE
C DT=(TT(2)/TT(1)-1.)/DLNTAU(2)
      NABSKO=0
      KK=1
      IF(PPE(1).LE.0.) PPE(1)=PEDEF

C
C ITERATE ON BOUNDARY CONDITION, USING PARTIAL DERIVATIVES
      i =0
100   CONTINUE
      KL=1
C      call rossos

      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1),1)
      PP(1)=PGC+PPT(1)+PPR(1)
      PG=PGC
      ROSST=CROSS(1)*ROSSOP(TT(1)*(1.+RELT),PPE(1),1)
      PGT=PGC
      ROSSPE=CROSS(1)*ROSSOP(TT(1),PPE(1)*(1.+RELPE),1)
      PGPE=PGC
      PGT=(PGT/PG-1.)/RELT
      PGPE=(PGPE/PG-1.)/RELPE
      ROSST=(ROSST/ROSS(1)-1.)/RELT
      ROSSPE=(ROSSPE/ROSS(1)-1.)/RELPE
      NABSKO=NABSKO+3
      DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST))/(1.+ROSSPE/PGPE)
      DP=MAX(DP,0.1D+0)
      

      DLNPE=log(GRAV*TAU(1)/(PG*ROSS(1)*DP))/(PGPE+ROSSPE)

      PPE(1)=PPE(1)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) then
      i = i+1
       if (i>200) then
        eps = eps*2.0
        i = 0
       end if
       GOTO 100
      else 
        eps = 1.0e-3
      end if
C
C END BOUNDARY CONDITION
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1),1)
      NABSKO=NABSKO+1
      PP(1)=PGC+PPT(1)+PPR(1)

C
C TAU LOOP
      DPE=(DP-DT*PGT)/PGPE
      DEDLNP=-(PGPE*PG/PP(1)+.5*DLNTAU(2)*GRAV*TAU(1)/(PP(1)*ROSS(1))*
     & (PGPE*PG/PP(1)+ROSSPE))

      DO 110 K=2,NTAU
      PPE(K)=PPE(K-1)*EXP(DPE*DLNTAU(K))

      NABSKO=0
C
C ITERATION LOOP
      DLNPE=0.
      i = 0 !counter for eps
111   CONTINUE
      KL=K
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K),k)
      PP(K)=PGC+PPT(K)+PPR(K)
      NABSKO=NABSKO+1
      ERROR=(.5*DLNTAU(K)*GRAV*(TAU(K-1)/(PP(K-1)*ROSS(K-1))+
     & TAU(K)/(PP(K)*ROSS(K)))-log(PP(K)/PP(K-1)))
      CALL ZEROF(ERROR,DLNPE,DEDLNP)   
      PPE(K)=PPE(K)*EXP(DLNPE)

C2023  format(2i4,1p7e12.3,0pf8.3)
      IF(ABS(DLNPE).GT.EPS) then
       i = i+1
       if (i>200) then
        eps = eps*2.0
        i = 0
       end if
       GOTO 111
      else 
        eps = 1.0e-3
      end if
C
C END TAU LOOP
C      call rossos
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K),k)
      NABSKO=NABSKO+1
      PP(K)=PGC+PPT(K)+PPR(K)
      DP=GRAV*TAU(K)/(PGC*ROSS(K))
      DPE=log(PPE(K)/PPE(K-1))/DLNTAU(K)

110   CONTINUE

      MSB=0
      MSB=MSA-MSB

      RETURN
      END
C
      FUNCTION TRQUAD(N,X,F,W)
      implicit real*8 (a-h,o-z)
C
      DIMENSION X(N),F(N),W(2*N)
C was : dim x(1) etc...
C
C TRAPEZOIDAL QUADRATURE PLUS NEXT ORDER CORRECTION FOR NON-
C -EQUIDISTANT GRID.
      N1=N-1
      Q=0.
      DO 100 K=2,N
      W(K)=X(K)-X(K-1)
      W(N+K)=(F(K)-F(K-1))/W(K)
100   Q=Q+W(K)*(F(K-1)+F(K))
      Q=Q*6.
      DO 101 K=2,N1
101   Q=Q+(W(K+1)-W(K))*(W(K)*W(N+K+1)+W(K+1)*W(N+K))
      W1=((W(2)+0.5*W(3))*W(N+2)-0.5*W(2)*W(N+3))*2.0/(W(2)+W(3))
      WN=((W(N)+0.5*W(N1))*W(N+N)-0.5*W(N)*W(N+N1))*2.0/(W(N)+W(N1))
      Q=0.083333333*(Q+W(2)**2*W1-W(N)**2*WN)
      TRQUAD=Q
      RETURN
      END
C
C
C
C ugj950523:  here the sphereical part begins @@@@@
C
C
      SUBROUTINE SOLVE_sph(NEW)
      implicit real*8 (a-h,o-z)
C
C SOLVE PERFORMS ONE NEWTON-RAPSON ITERATION ON THE MODELATMOSPHERE PROBLEM
C INCLUDING LOCAL CONVECTION.THE STATE OF THE ATMOSPHERE IS DESCRIBED BY A
C NUMBER OF VARIABLES SUCH AS TEMPERATURE,ELECTRON PRESSURE,TOTAL PRESSURE,
C CONVECTIVE FLUX ETC..TO EACH VARIABLE CORRESPONDS A CERTAIN CONDITIONAL
C EQUATION WICH DETERMINES THAT VARIABLE,ASSUMING THE OTHER VARIABLES BEEING
C KNOWN.
C
C NAMING CONVENTION.THE VARIABLES HAVE NAMES WITH A DOUBLE OCCURANCE OF THE
C FIRST LETTER,CORRECTIONS (TO BEE COMPUTED IN THIS ITERATION) HAVE SINGLE
C OCCURANCE OF FIRST LETTER.RIGHTHANDSIDENAMES BEGIN WITH R.
C
C VARIABLES ARE CENTERED ON INTEGER AND HALFINTEGER TAU-POINTS AS INDICATED
C BY 'I' OR 'H' ON THE FOLLOWING COMMENT CARDS.
C
C VARIABLE CORRECTION                                           HALF-INTEGER
C PPR      PR         RADIATION PRESSURE                             I
C PPT      PT         TURBULENT PRESSURE                        H
C PP       P          TOTAL PRESSURE                                 I
C GG       G          CONVECTIVE EFFICIENCY,GAMMA               H
C ZZ       Z          GEOMETRIC HEIGTH                          H
C DD       D          GRADIENT DIFFERENCE,DELTA-DELTAPRIME      H
C VV       V          CONVECTIVE VELOCITY                       H
C FFC      FC         CONVECTIVE FLUX                           H
C PPE      PE         ELECTRON PRESSURE                              I
C TT       T          TEMPERATURE                                    I
C XJ                  MEAN INTENSITY                                 I
C
C NAMES OF THE PARTIAL DERIVATIVES ARE FORMED WITH A FIRST PART FROM THE
C EQUATION TO WICH IT BELONGS,AND A SECOND PART WICH IS THE NAME OF THE
C VARIABLE WITH RESPECT TO WICH THE DERIVATIVE IS TAKEN.(DOUBLE OCCURANCE
C IF NECESSARY TO AVOID CONFUSION).
C
C PRPR,PRT
C PTPT,PTV,PTPE,PTT
C PPP,PPPE,PPTT
C GGG,GV,GPE,GT
C ZZZ,ZPE,ZT
C DDD,DP,DG,DPE,DT
C VVV,VZ,VD,VPE,VT
C FCFC,FCP,FCD,FCV,FCPE,FCT
C PEPE,PEPR,PEPT,PEP,PET
C TXJ,TFC,TTT
C
C  PARAMETRIZED VERSION (in NDP and NRAYS).     B.PLEZ 20-NOV-88. 
C
      include 'parameter.inc'
C
C STATE VARIABLES
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),ROSTATEC(NDP),
     & NTAU,ITER
      common /ckdtpe/dpex,kdtpe
      common /dpeset/ dpein,dtin

C
C DIMENSIONS
      DIMENSION PTAU(NDP),ROSSP(NDP),SUMW(NDP),ROSST(NDP),ROSSPE(NDP)
     *,XL(500),W(500)
     *,XJ1(NDP),XJ2(NDP),XJ3(NDP),XJT1(NDP),XJT2(NDP),XJT3(NDP)
     *,XJPE1(NDP),XJPE2(NDP),XJPE3(NDP)
     *,PR(NDP),PRT(NDP,NDP),PRJ(NDP)
     *,PRPE(NDP,NDP)
     *,PT(NDP),PTV(NDP),PTPE(NDP),PTT(NDP)
     *,P(NDP),PPPE(2*NDP),PPTT(2*NDP)
     *,GV(NDP),GPE(NDP),GT(NDP)
     *,DP(2*NDP),DG(NDP),DPE(NDP,NDP),DT(NDP,NDP),DV(NDP)
      DIMENSION D(NDP),DTS(2*NDP),DPS(2*NDP),DPES(2*NDP)
     *,V(NDP),VD(NDP),VPE(NDP,NDP),VT(NDP,NDP)
     *,FC(NDP),FCD(NDP),FCV(NDP),FCPE(NDP,NDP),FCT(NDP,NDP)
     *,PE(NDP),PEPE(NDP,NDP),PET(NDP,NDP)
     *,T(NDP),TTT(NDP,NDP),TPE(NDP,NDP),TJ1(NDP),TJ2(NDP)
     *,SCRATC(NDP,NDP),DBPL(NDP)
     *,XT(NDP),ST(NDP),DLNX(NDP),XLOG(NDP)
     *,XPE(NDP),SPE(NDP)
     *,RPR(NDP),RP(NDP),RD(NDP),RV(NDP),RFC(NDP),RPE(NDP),RT(NDP)
     *,TLAST(NDP),TVD(NDP),IA(5)
      LOGICAL NEWV, exist
C
C CONNECTIONS VIA COMMON.
C THE COMMENTED COMMONS MUST BE INITIATED OUTSIDE THIS ROUTINE BEFORE IT
C IS CALLED.
C JTAU=NUMBER OF TAUPOINTS, TAU=TAUSCALE.
C NLAM=NUMBER OF LAMBDAPOINTS, XL=LAMBDAPOINTS, W=INTEGRATIONWEIGHTS.
C MIHAL=LOWER LIMIT OF RADIATIVE EQUILIBRIUM CONDITION, TAUMAX NOT USED.
C PALFA,PBETA,PNY,PNY = MIXING LENGTH THEORY COEFFICIENTS.
C GRAV=SURFACE GRAVITY, TEFF=EFFECTIVE TEMPERATURE, FLUX=STEFAN*TEFF**4/PI
      CHARACTER MOLNAME*4,OSFIL*60,SAMPLING*3
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      COMMON /CLEVETAT/GEFF(NDP),PPRG(NDP),AMLOSS
      COMMON /CLEVPRINT/ PRJ2(NDP),masslinf
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      common /CPRINT/NPRINT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /CVAAGL/XL,W,NLAM
      COMMON /CSTYR/MIHAL,NOCONV /DEBUG/KDEBUG
      COMMON /MIXC/PALFA,PBETA,PNY,PY /CVFIX/VFIX
      COMMON /CG/GRAV,KONSG /CTEFF/TEFF,FLUX
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     * STEFAN
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
CUGJ FFR in excess     COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP),FFR(NDP)
      dimension ffr(ndp)
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE 
C OWN COMMONS
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /CSURF/HSURF,Y1(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP) /RHOC/RHO(NDP)
      COMMON /CARC1/ISTRAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &              IARCH
      COMMON /CARC2/T,FC,FLUXME(NWL),TAU5(NDP),INORD
      COMMON /Cspec/spec(nwl,3),ispec
      COMMON /CI8/PGC,ROC,EC
      COMMON /NEWMO/NEWMOD
      COMMON /MASSE/RELM
      COMMON /CORRECT/TDIFF,TCONV,KORT
      COMMON /CIT/IT,ITMAX
C
C SPACE ALLOCATION
      COMMON /SPACE1/XJ1,XJ2,XJ3,TJ1,TJ2,XJT1,XJT2,XJT3,PRJ,TTT,PRT
     &  ,XJPE1,XJPE2,XJPE3,TPE,PRPE
      COMMON /SPACE2/ SPACEDUM1(NDP*7+NDP*NRAYS*5+NRAYS*2),
     &       SPACEDUM2(NDP*2),PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),
     &       MMU(NDP),KSPACE2DUM(NRAYS+1)


C
C IN THIS SECTION THE MEAN INTENSITY IS ELIMINATED IN THE TRANSPORT EQUATIONS
C LEAVING THE EXPLICIT TEMPERATURE DEPENDANCE OF FLUX AND RADIATION PRESSURE
C IN THE MATRICES TTT AND PRT.
C
      ITER=ITER+1
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
C
      
C
C ZEROSET
      DO 110 I=1,NTAU
      PPRG(I)=0.
      RT(I)=0.
      RPR(I)=0.
      FFR(I)=0.
      ROSSP(I)=0.
      SUMW(I)=0.
      DO 110 J=1,NTAU
        TTT(I,J)=0.
        TPE(I,J)=0.
        PRT(I,J)=0.
        PRPE(I,J)=0.
110   CONTINUE
      kdtpe = 0
C
C CALCULATE DETAILED ROSSELAND MEAN
      REWIND 11
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1),1)
      PGA=PGC
      DO 116 K=1,NTAU
        KL=K
	DUMMY=ROSSOP(TT(K),PPE(K),k)
	ROSS(K)=1.0
        RHO(K)=ROC
	SUMW(K)=0.
	ROSSP(K)=0.
116   CONTINUE
C      DO 117 J=1,NLAM
      !print*, "opac call for ROSS"
      DO 117 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(11) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 117 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSSP(K)=ROSSP(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 117
          ROSSP(K)=ROSSP(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
117   CONTINUE
      REWIND 11
C
C TEMPERATURE AND ELECTRON PRESSURE PERTURBATIONS.
C KEEP THEM SMALL, TO STAY ON THE LINEAR PART.
C      DTX=0.001
C      DPEX=0.001
      DTX=dtin
      DPEX=dpein
      DO 111 K=1,NTAU
        T(K)=TT(K)*DTX
        PE(K)=PPE(K)*DPEX
        TT(K)=TT(K)+T(K)
	ROSSP(K)=SUMW(K)/ROSSP(K)
	SUMW(K)=0.
	ROSST(K)=0.
111   CONTINUE
      kdtpe = 1
C
C FIRST WAVELENGTH LOOP, TO CALCULATE XT,ST AND SAVE.
      REWIND 12
C      DO 112 J=1,NLAM
      DO 112 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(12) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 112 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSST(K)=ROSST(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 112
          ROSST(K)=ROSST(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
112   CONTINUE
      DO 113 K=1,NTAU
	ROSST(K)=SUMW(K)/ROSST(K)
	SUMW(K)=0.
	ROSSPE(K)=0.
        TT(K)=TT(K)-T(K)
        PPE(K)=PPE(K)+PE(K)
113   CONTINUE
        kdtpe = 2            !information to tstgem about computing dpg/dpe
      REWIND 12
      CALL CLOCK
C
C SECOND WAVELENGTH LOOP, TO CALCULATE XPE,SPE AND SAVE.
      REWIND 14
      KL=1
      DUMMY=ROSSOP(TT(1),PPE(1),1)
      PGPE=PGC
C      DO 114 J=1,NLAM
      DO 114 J=1,NWTOT
        CALL OPAC(J,X,S)
        WRITE(14) X,S
C        Y=((XL(J)/1.E4)**2)**3
        Y=((WLOS(J)/1.E4)**2)**3
        DO 114 K=1,NTAU
C          YA=EXP(-1.438E8/(TT(K)*XL(J)))
          YA=EXP(-1.438E8/(TT(K)*WLOS(J)))
          YA=YA/(1.-YA)**2/Y
C          SUMW(K)=SUMW(K)+W(J)*YA
C          ROSSPE(K)=ROSSPE(K)+W(J)*YA/(ROSS(K)*(X(K)+S(K)))
          SUMW(K)=SUMW(K)+WLSTEP(J)*YA
        if (wlos(j).le.5000. .or. wlos(j).ge.1.e5) go to 114
          ROSSPE(K)=ROSSPE(K)+WLSTEP(J)*YA/(ROSS(K)*(X(K)+S(K)))
114   CONTINUE
      kdtpe = 3
      REWIND 14
      CALL CLOCK
C
C FROM THIS POINT ON, ROSS() HOLDS THE TRUE ROSSELAND MEAN.  CROSS HOLDS
C THE RATIO OF THE TRUE TO APPROXIMATE MEANS, WHICH ARE NEEDED IN TRYCK.
      
      DO 123 K=1,NTAU
        KL=K
	ROSSPE(K)=SUMW(K)/ROSSPE(K)
        PPE(K)=PPE(K)-PE(K)
	CROSS(K)=ROSSP(K)/ROSSOP(TT(K),PPE(K),k)
        ROSS(K)=ROSSP(K)
        PTAU(K)=GRAV*TAU(K)/ROSS(K)
	
123   CONTINUE
      CALL CLOCK
C
C RIGHT HAND SIDE IN PRESSURE EQUATION
      KL=1
      CALL TAET(TT(1),PPE(1),PG,RO,DUM)
      DLNP=1./(1.+(ROSSPE(1)-ROSS(1))/ROSS(1)*PGA/(PGPE-PGA))
      GRVR=GRAV*(RADIUS/RR(1))**2
      IF (KONSG.EQ.1) GRVR=GRAV   !test for effect of varying g
      RP(1)=GRVR*TAU(1)/(ROSS(1)*DLNP)+PPR(1)-PP(1)
C SIMPSONS RULE
      DO 101 K=2,NTAU
      F0=PTAU(K-1)
      F1=FOUR(PTAU,TAULN,K,NTAU)
      F2=PTAU(K)
      RP(K)=(F0+4.*F1+F2)*DTAULN(K)/6.-(PP(K)-PP(K-1))
101   CONTINUE
C
C CALCULATE RADII
      RR(1)=0.
      DO 102 K=2,NTAU
      IF (TAU(K).LT.0.67) K0=K
102   RR(K)=RR(K-1)-0.5*DTAULN(K)
     & *(TAU(K)/(ROSS(K)*RHO(K))+TAU(K-1)/(ROSS(K-1)*RHO(K-1)))

      Y=RADIUS-RR(K0)
      DO 103 K=1,NTAU
      RR(K)=RR(K)+Y
      IF (KONSG.EQ.1) RR(K)=RADIUS      !study effect of const R and g
103   PTAU(K)=PTAU(K)*(RADIUS/RR(K))**2
C
      
104   CONTINUE
C
C TIME
      CALL CLOCK
      MSA=0
C
C WAVELENGTH LOOP
      
      FTOT=0.

      DO 150 J=1,NWTOT
      DO 130 K=1,NTAU
      BPLAN(K)=BPL(TT(K),WLOS(J))
130   DBPL(K)=DIVBP(TT(K),WLOS(J))
      
C
C CALCULATE OPACITY DERIVATIVES AT CONSTANT GAS PRESSURE
      READ (11) X,S
      READ (12) XT,ST
      READ (14) XPE,SPE
      DO 131 K=1,NTAU
        X(K)=X(K)/ROSS(K)
        S(K)=S(K)/ROSS(K)
        XT(K)=XT(K)/ROSST(K)
        ST(K)=ST(K)/ROSST(K)
        XPE(K)=XPE(K)/ROSSPE(K)
        SPE(K)=SPE(K)/ROSSPE(K)
C	if (j.eq.100.and.k.eq.1) write (7,'(1x,7(1pe10.2))')
C     &   x(k),t(k),xt(k),log(xt(k)/x(k)),log(xt(k)/x(k))*tt(k)/t(k)
        XLOG(K)=log10(X(K))
        XT(K)=log(XT(K)/X(K))*X(K)/T(K)
        ST(K)=log(ST(K)/S(K))*S(K)/T(K)
        XPE(K)=log(XPE(K)/X(K))*X(K)/PE(K)
        SPE(K)=log(SPE(K)/S(K))*S(K)/PE(K)
        DLNX(K)=XT(K)*(TT(K)/2.3)/X(K)
C	if (j.eq.100) write (7,'(1x,7(1pe10.2))')
C     &    x(k),xt(k)/x(k)*tt(k),xpe(k)/x(k)*ppe(k)
131   CONTINUE
C TIME
      MS=MSA
      MSA=0
      MSOPAC=MS-MSA
C
C SOLVE TRANSPORTEQUATION WITH OLD STRATIFICATION.
      CALL TRANEQ_sph
C
C ??
      DO 132 K=1,NTAU
      IF(XT(K).LT.0.0.AND.XT(K)*(XJ(K)-BPLAN(K)).GT.X(K)*DBPL(K))
     & XT(K)=X(K)*DBPL(K)/(XJ(K)-BPLAN(K))
132   CONTINUE
      MS=MSA
      MSA=0
      MSTRAN=MS-MSA
C
C FLUX TO PRINT
C      HFLUX1=4.*PI*HSURF*(RR(1)/RADIUS)**2
      HFLUX1=4.*PI*HSURF
      HFLUX1=MAX(1.0D-99,HFLUX1)
      HFLUX2=4.*PI*HFLUX(NTAU)*(RR(NTAU)/RADIUS)**2
      FLUXME(J)=HFLUX1/PI
C      GFLUX1=HFLUX1/PI*XL(J)**2/CLIGHT
C      GFLUX2=HFLUX2/PI*XL(J)**2/CLIGHT
      GFLUX1=HFLUX1/PI*WLOS(J)**2/CLIGHT
      GFLUX2=HFLUX2/PI*WLOS(J)**2/CLIGHT
      FFLUX1=-2.5*log10(MAX(GFLUX1,1.0D-99))
      FFLUX2=-2.5*log10(MAX(GFLUX2,1.0D-99))
      spec(j,1) = wlos(j)
      spec(j,2) = hsurf
      spec(j,3) = fluxme(j)
C

C
C SUM UP RADIATIVE FLUXES
      DO 133 K=1,NTAU
      FFR(K)=FFR(K)+WLSTEP(J)*HFLUX(K)    
133   CONTINUE
C
C UPPER BOUNDARY
      DO 140 K=1,NTAU
      IF (K.GT.1) GO TO 143
      PB=XJ(1)/FJ(1)
      PC=XJ(2)/FJ(2)
      IF (TAUS(1).LT.0.1) GO TO 141
      EX=EXP(-TAUS(1))
      EX1=1.-EX
      GO TO 142
141   EX1=TAUS(1)*(1.-0.5*TAUS(1)*(1.-0.333333*TAUS(1)))
      EX=1.-EX1
142   YA=DTAUS(2)*(EX1+0.5*DTAUS(2))/(X(1)+S(1))
      YB=((1.+DTAUS(2))*(PB-SOURCE(1))+SOURCE(1)*EX)*DTAUS(2)
     & /(X(1)+S(1)+X(2)+S(2))
      XJ2(1)=DTAUS(2)+0.5*DTAUS(2)**2-YA*S(1)*FJ(1)
      XJ3(1)=-1.
      XJT2(1)=-YA*X(1)*DBPL(1)+YB*(XT(1)+ST(1))
     & -YA*(BPLAN(1)-XJ(1))/(X(1)+S(1))*(S(1)*XT(1)-X(1)*ST(1))
     & -DTAUS(2)*SOURCE(1)*EX*TAUS(1)/(X(1)+S(1))*(XT(1)+ST(1))
      XJT3(1)=YB*(XT(2)+ST(2))
      XJPE2(1)=YB*(XPE(1)+SPE(1))
     & -YA*(BPLAN(1)-XJ(1))/(X(1)+S(1))*(S(1)*XPE(1)-X(1)*SPE(1))
     & -DTAUS(2)*SOURCE(1)*EX*TAUS(1)/(X(1)+S(1))*(XPE(1)+SPE(1))
      XJPE3(1)=YB*(XPE(2)+SPE(2))
      GO TO 170
C
C INTERNAL POINTS
143   PA=PB
      PB=PC
      IF (K.EQ.JTAU) GO TO 144
      PC=XJ(K+1)/FJ(K+1)
      DTAUSK=0.5*(DTAUS(K+1)+DTAUS(K))
      XJ1(K)=1./DTAUS(K)
      XJ3(K)=1./DTAUS(K+1)
      XJ2(K)=-DTAUSK*(1.-FJ(K)*S(K)/(X(K)+S(K)))
      YA=((PB-PA)/DTAUS(K)-(PB-SOURCE(K))*0.5*DTAUS(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      YB=(-(PC-PB)/DTAUS(K+1)-(PB-SOURCE(K))*0.5*DTAUS(K+1))
     & /(X(K)+S(K)+X(K+1)+S(K+1))
      XJT1(K)=YA*(XT(K-1)+ST(K-1))
      XJT3(K)=YB*(XT(K+1)+ST(K+1))
      XJT2(K)=DTAUSK*X(K)/(X(K)+S(K))*DBPL(K)+(YA+YB)*(XT(K)+ST(K))
     & +DTAUSK*(BPLAN(K)-XJ(K))/(X(K)+S(K))**2*(S(K)*XT(K)-X(K)*ST(K))
      XJPE1(K)=YA*(XPE(K-1)+SPE(K-1))
      XJPE3(K)=YB*(XPE(K+1)+SPE(K+1))
      XJPE2(K)=(YA+YB)*(XPE(K)+SPE(K))
     & +DTAUSK*(BPLAN(K)-XJ(K))/(X(K)+S(K))**2*(S(K)*XPE(K)-X(K)*SPE(K))
      GO TO 170
C
C OPTICALLY THICK POINTS
144   XJ1(K)=0.
      XJ2(K)=-1.
      XJ3(K)=0.
      XJT1(K)=0.
      XJT2(K)=DBPL(K)
      XJT3(K)=0.
      XJPE1(K)=0.
      XJPE2(K)=0.
      XJPE3(K)=0.
C
C TEMPERATURE EQUATION
170   IF (K.GT.MIHAL) GO TO 171
C      Y=W(J)*X(K)
      Y=WLSTEP(J)*X(K)
      IF (K.GT.2) Y=Y*((TAU(K+1)-TAU(K))*(X(K)+S(K)+X(K+1)+S(K+1))
     & +(TAU(K)-TAU(K-1))*(X(K)+S(K)+X(K-1)+S(K-1)))/(X(K)+S(K))
      RT(K)=RT(K)+Y*(XJ(K)-BPLAN(K))
      TJ2(K)=-Y*FJ(K)
      TJ1(K)=0.
c ??
***      TTT(K,K)=TTT(K,K)+AMAX1(0.,Y*DBPL(K)+Y*(BPLAN(K)-XJ(K))
***     & *XT(K)/X(K))
      TTT(K,K)=TTT(K,K)+Y*DBPL(K)+Y*(BPLAN(K)-XJ(K))*XT(K)/X(K)
      TPE(K,K)=TPE(K,K)+Y*(BPLAN(K)-XJ(K))*XPE(K)/X(K)
      GO TO 172
C171   Y=4.*W(J)
171   Y=4.*WLSTEP(J)
      RT(K)=RT(K)-Y*HFLUX(K)
      XHK=0.557*(PB-PA)/DTAUS(K)
      FH=HFLUX(K)/XHK
C
C DEBUG
      IF (FH.GT.0.7.AND.FH.LT.1.4) GO TO 174
      
      IF (K.GT.JTAU1) GO TO 174
      NMU=MMU(K-1)
      NMU=MMU(K)

174   CONTINUE
C
      TJ1(K)=-Y*FH*0.557/DTAUS(K)
      TJ2(K)=Y*FH*0.557/DTAUS(K)
      TTT(K,K-1)=TTT(K,K-1)-Y*HFLUX(K)*(XT(K-1)+ST(K-1))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TTT(K,K)=TTT(K,K)-Y*HFLUX(K)*(XT(K)+ST(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TPE(K,K-1)=TPE(K,K-1)-Y*HFLUX(K)*(XPE(K-1)+SPE(K-1))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
      TPE(K,K)=TPE(K,K)-Y*HFLUX(K)*(XPE(K)+SPE(K))
     & /(X(K)+S(K)+X(K-1)+S(K-1))
172   CONTINUE
C
C EQUATION OF RADIATIVE PRESSURE
C      Y=PI4C*W(J)
      Y=PI4C*WLSTEP(J)
      RPR(K)=RPR(K)+Y*XK(K)
      PRJ(K)=-Y*XK(K)*FJ(K)/XJ(K)
      PRJ2(K)=PRJ(K)
CUGJ 900503: The gradient of the radiative pressure
      PPRG(K) = PI4C*WLSTEP(J)*HFLUX(K)*(X(K)+S(K))*ROSS(K) + PPRG(K)
C
C END OF TAU LOOP
140   CONTINUE

C
C ELIMINATE THIS WAVELENGTH
      IF(NEW.EQ.1) CALL ALGEBN(NTAU)
C
C TIME
      MS=MSA
      MSA=0
      MS=MS-MSA
C
C END OF WAVELENGTH LOOP
C      HW1=HFLUX1*W(J)
      HW1=HFLUX1*WLSTEP(J)
      FTOT=FTOT+HW1
      HW2=HFLUX2*WLSTEP(J)
C      HW2=HFLUX2*W(J)
C      WAVEN=1.E4/XL(J)
C      TRAD1=1.438E8/(XL(J)*log(1.+1.191E7*PI/HFLUX1*WAVEN**5))
      WAVEN=1.E4/WLOS(J)
      TRAD1=1.438E8/(WLOS(J)*
     #  log(1.0D+0+1.191D7*PI/HFLUX1*WAVEN**5))
      !X01=log10(X(01))
      !X25=log10(X(25))
      !S01=log10(S(01))
      !S25=log10(S(25))

150   CONTINUE
C
   
      TEFFP=TEFF*(FTOT/(FLUX*PI))**.25
      print*, "Current Teff... ", teffp
      Y=FTOT/PI
      
      DO 154 K=1,NTAU
      Y=0.
      DO 155 L=1,NTAU
155   Y=Y+TTT(K,L)
      
      IF (K.LE.MIHAL.AND.Y.LT.0.0) TTT(K,K)=TTT(K,K)-Y
154   CONTINUE
      DO 153 K=1,NTAU
153   FFR(K)=FFR(K)*4./FLUX*(RR(K)/RADIUS)**2

C
C TIME
      CALL CLOCK
C
C PRINT PRESSURE EQUATION
      
      DO 161 K=1,NTAU
      KL=K
      ROSST(K)=(ROSST(K)-ROSS(K))/T(K)
      ROSSPE(K)=(ROSSPE(K)-ROSS(K))/PE(K)
C
C TAU SCALES
      IF(K.GT.1) GO TO 162

      CALL KAP5(TT(1),PPE(1),ABSK,1)
      TAU5(1)=TAU(1)*ABSK/ROSS(1)
      TAUP=TAU(1)*CROSS(1)
      YC=ABSK/ROSS(1)
      YD=CROSS(1)
      GO TO 163
162   CONTINUE
      CALL KAP5(TT(K),PPE(K),ABSK,k)
      YA=ABSK/ROSS(K)
      YB=CROSS(K)
      TAU5(K)=TAU5(K-1)+.5*(YA+YC)*(TAU(K)-TAU(K-1))
      TAUP=TAUP+.5*(YB+YD)*(TAU(K)-TAU(K-1))
      YC=YA
      YD=YB
163   CONTINUE
      ROSPE=ROSSPE(K)*PPE(K)/ROSS(K)
      ROST=ROSST(K)*TT(K)/ROSS(K)
     
161   CONTINUE
      CALL CLOCK
900   CONTINUE

      RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C+DEF,D=MATRIX
C
C IN THIS SECTION WE COMPUTE MATRIX ELEMENTS FOR THE REST OF THE PROBLEM
C
C
      ENTRY MATRIX_sph
C TIME
      MSA=0
C
C ZEROSET
      
      GRAD=0.
      DO 201 I=1,NTAU
      DO 201 J=1,NTAU
      VT(I,J)=0.
      VPE(I,J)=0.
      FCT(I,J)=0.
      FCPE(I,J)=0.
      DT(I,J)=0.
      DPE(I,J)=0.
      PET(I,J)=0.
      PEPE(I,J)=0.
201   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 230
      DO 231 KK=2,NTAU
      K=(1+NTAU)-KK
      IF(K.LE.NOCONV) GOTO 230
231   VV(K)=MAX(VV(K+1)*EXP(-DTAULN(K+1)/VFIX),VV(K))
230   CONTINUE
C
C TAU LOOP
      DO 200 K=1,NTAU
      KL=K
C
C TERMODYNAMICAL QUANTITIES WITH PARTIAL DERIVATIVES
      K1=MAX0(1,K-1)
      TMEAN=.5*(TT(K)+TT(K1))
      PEMEAN=.5*(PPE(K)+PPE(K1))
      PRMEAN=.5*(PPR(K)+PPR(K1))
      ROSSMN=.5*(ROSS(K)+ROSS(K1))
      IF(K.GT.NOCONV) GO TO 213
C NO CONVECTION
      CALL TERMON(K,TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO,ROT,
     &                  ROPE,CP,ADIA,Q)
      RO1=RO
      PG1=PG
      PG1T=0.
      PG1PE=0.
      CPT=0.
      CPPE=0.
      ADIAT=0.
      ADIAPE=0.
      QT=0.
      QPE=0.
      GO TO 212
C CONVECTION
213   CONTINUE
      DERET=.01
      DEREP=.15
      TDELT=TMEAN*DERET
      PEDELT=PEMEAN*DEREP
      T1=TMEAN+TDELT
      PE1=PEMEAN+PEDELT
      CALL TERMON(K,T1,PEMEAN,PRMEAN,Y,YA,YB,YC,YD,YE,CP1,
     &        ADIA1,Q1)
      CALL TERMON(K,TMEAN,PE1,PRMEAN,Y,YA,YB,YC,YD,YE,CP2,
     &        ADIA2,Q2)
      CALL TERMON(K,TMEAN,PEMEAN,PRMEAN,PG1,PG1T,PG1PE,RO,
     &        ROT,ROPE,CP,ADIA,Q)
      CALL TERMON(K,TT(K),PPE(K),PPR(K),PG,PGT,PGPE,RO1,Y,
     &        YA,YB,YC,YD)   
      CPPE=(CP2-CP)/PEDELT
      ADIAPE=(ADIA2-ADIA)/PEDELT
      QPE=(Q2-Q)/PEDELT
      CPT=(CP1-CP)/TDELT
      ADIAT=(ADIA1-ADIA)/TDELT
      QT=(Q1-Q)/TDELT
212   CONTINUE
      RHO(K)=RO1
C
C DEPTH SCALE
      IF(K.EQ.1) ZZ(K)=0.
      IF(K.GT.1) ZZ(K)=ZZ(K-1)-.5*DTAULN(K)*
     & (TAU(K)/ROSS(K)/RHO(K)+TAU(K-1)/ROSS(K-1)/RHO(K-1))
      IF(TAULN(K).LT.0.0) KK0=K
C
C RADIATION PRESSURE
      RPR(K)=RPR(K)-PPR(K)
CC    PRJ IS NOW FREE
CC    PRT IS ALREADY INITIATED
CC    PRPR IS UNITY
C
C TOTAL PRESSURE
      GRVR=GRAV*(RADIUS/RR(K))**2
      GEFF(K)=GRVR
      IF (KONSG.EQ.1) GRVR=GRAV   !test for effect of varying g
      IF(K.GT.1) GO TO 202
      Y=GRVR*TAU(1)/(ROSS(1)*DLNP)
      PPTT(1)=Y*ROSST(1)/ROSS(1)
      PPPE(1)=Y*ROSSPE(1)/ROSS(1)
      GO TO 203
202   Y=GRVR*DTAULN(K)*.5
      YY=Y*TAU(K)/ROSS(K)**2
      Y=Y*TAU(K-1)/ROSS(K-1)**2
      PPPE(K)=YY*ROSSPE(K)
      PPTT(K)=YY*ROSST(K)
      PPPE(K+NTAU-1)=Y*ROSSPE(K-1)
      PPTT(K+NTAU-1)=Y*ROSST(K-1)
203   CONTINUE
C
C CONVECTION EFFICIENCY GAMMA
      HSCALE=(PG1+PRMEAN)/GRVR/RO
      OMEGA=PALFA*HSCALE*RO*ROSSMN
      IF(PALFA.EQ.0.) OMEGA=HSCALE*RO*ROSSMN
      Y=PY*OMEGA**2
      YY=(Y-1.)/(Y+1.)
      THETA=OMEGA/(1.+Y)
      GAMMA=-CP*RO/(8.*STEFAN*TMEAN**3*THETA)
CC    GGG IS UNITY
      GV(K)=GAMMA
      IF(PBETA.GT.0.) VV(K)=MIN(VV(K),SQRT(0.5*PP(K)/PBETA/RO))
      GG(K)=-GAMMA*VV(K)
      ROSPEM=.5*(ROSSPE(K)+ROSSPE(K1))
      ROSSTM=.5*(ROSST(K)+ROSST(K1))
      GPE(K)=-GG(K)*(CPPE/CP+ROPE/RO+YY*(ROSPEM/ROSSMN+PG1PE/PG1))
      GT(K)=-GG(K)*(CPT/CP+ROT/RO-3./TMEAN+YY*(ROSSTM/ROSSMN+PG1T/PG1))
CC    RG IS ZERO
C
C GRADIENT DIFFERENCE
      IF(K.LE.NOCONV) GO TO 206
      DELP=PP(K)-PP(K-1)
      DELT=TT(K)-TT(K-1)
      PM=PP(K)+PP(K-1)
      TM=TT(K)+TT(K-1)
      Y=1.+GG(K)
      YY=-GG(K)/Y
      GRAD=log(TT(K)/TT(K-1))/log(PP(K)/PP(K-1))
      NEWV=DD(K).GT.0..AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 263
      VV(K)=SQRT(GRVR*HSCALE*Q*PALFA**2*DD(K)/PNY)
      
      GO TO 203
263   CONTINUE
      NEWV=GRAD.GE.ADIA.AND.VV(K).EQ.0..AND.PALFA.GT.0..AND.K.GT.2
      IF(.NOT.NEWV) GO TO 204
      VV(K)=VVMLT(GRAD-ADIA,GRVR*HSCALE*Q*PALFA**2/PNY,GAMMA**2)
      
      GO TO 203
204   CONTINUE
      YYY=GRAD-ADIA
C DDD IS UNITY
C ******* NEXT STATEMENT FIXES T80G4M0 BUT NOT TESTED FOR ALL MODELS
      IF(DD(K).EQ.0..AND.VV(K).GT.0.) DD(K)=-YY*YYY
      RD(K)=-YY*YYY-DD(K)
      DG(K)=-YYY/Y**2
      DT(K,K)=YY*(GRAD*(1./DELT-1./TM)-.5*ADIAT)
      DT(K,K-1)=YY*(GRAD*(-1./DELT-1./TM)-.5*ADIAT)
      DP(K)=YY*GRAD*(1./PM-1./DELP)
      DP(K+NTAU-1)=YY*GRAD*(1./PM+1./DELP)
      DPE(K,K)=-.5*YY*ADIAPE
      DPE(K,K-1)=-.5*YY*ADIAPE
      GO TO 205
206   DD(K)=0.
      RD(K)=0.
      DG(K)=0.
      DP(K)=0.
      DP(K+NTAU-1)=0.
205   CONTINUE
C
C VFIX OPTION
      IF(VFIX.EQ.0.) GOTO 280
      Y=0.
      IF(K.GT.NOCONV) Y=-VV(K)
      GOTO 207
280   CONTINUE
C
C CONVECTIVE VELOCITY
CC    VVV IS UNITY
      Y=0.
      IF(DD(K).LE.0.) GOTO 207
      Y=-SQRT(GRVR*HSCALE*Q*PALFA**2*DD(K)/PNY)
      VD(K)=Y*.5/DD(K)
      IF(-Y.GT.2.*VV(K)) VD(K)=VD(K)*2.
      VT(K,K)=.25*Y*(QT/Q+PGT/PG-ROT/RO)
      VT(K,K-1)=VT(K,K)
      VPE(K,K)=.25*Y*(QPE/Q+PGPE/PG-ROPE/RO)
      VPE(K,K-1)=VPE(K,K)
      GO TO 208
207   VD(K)=0.
208   CONTINUE
      RV(K)=-Y-VV(K)
C
C TURBULENT PRESSURE.
CC    PTPT IS UNITY
      Y=-PBETA*VV(K)**2
      PPT(K)=-RO*Y
      PTT(K)=Y*ROT
      PTPE(K)=Y*ROPE
      PTV(K)=-PBETA*2.*VV(K)*RO
C
C CONVECTIVE FLUX
      Y=-CP*RO*PALFA*TMEAN/2./PI
CC    FCFC IS UNITY
      YY=Y*VV(K)*DD(K)
      RFC(K)=-YY-FFC(K)
      FCD(K)=Y*VV(K)
      FCV(K)=Y*DD(K)
      IF(K.LE.NOCONV) GO TO 217
      FCT(K,K)=.5*YY*(CPT/CP+ROT/RO+1./TMEAN)
      FCT(K,K-1)=FCT(K,K)
      FCPE(K,K)=.5*YY*(CPPE/CP+ROPE/RO)
      FCPE(K,K-1)=FCPE(K,K)
217   CONTINUE
C
C ELECTRON PRESSURE
209   RPE(K)=PP(K)-PG-PPR(K)
      PET(K,K)=PGT
      PEPE(K,K)=PGPE
CC    PEPR=PEPT=1.  PEP=-1.
210   CONTINUE
C
C TEMPERATURE
CC    TTT IS ALREADY INITIATED
      IF(K.GT.MIHAL) GO TO 261
C STRMGREN CONDITION
      
      RT(K)=RT(K)+(FFC(K+1)-FFC(K))
      GO TO 262
261   CONTINUE
C FLUXCONSTANCY
      RT(K)=RT(K)+FLUX*(RADIUS/RR(K))**2-FFC(K)
262   CONTINUE
C
C END OF TAU LOOP
      PGT=PGT*TT(K)/PG
      PGPE=PGPE*PPE(K)/PG
      
200   CONTINUE
C
C SUBTRACT CENTERED TURBULENT PRESSURE
      DO 216 K=2,NTAU
      K1=MIN0(K+1,NTAU)
216   RPE(K)=RPE(K)-.5*(PPT(K)+PPT(K1))
C
C SUBTRACT ZZ(TAU=1) FROM ZZ
      KK0 = MAX(1,KK0)
      ZZ0=ZZ(KK0)
      DO 283 K=1,NTAU
283   ZZ(K)=ZZ(K)-ZZ0
C
C TIME
      CALL CLOCK

C SAVE DD-MATRICES
      DO 270 K=1,NTAU
      DPS(K)=DP(K)
      DPES(K)=DPE(K,K)
      DTS(K)=DT(K,K)
      D(K)=RD(K)
      IF(K.EQ.1) GO TO 270
      DPS(K+NTAU-1)=DP(K+NTAU-1)
      DPES(K+NTAU-1)=DPE(K,K-1)
      DTS(K+NTAU-1)=DT(K,K-1)
270   CONTINUE
C
C+DEF,D=ELIMIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GAUSS ELIMINATION TO UPPER TRIANGULAR FORM.
C
C TIME
      CALL CLOCK
C
C RADIATION PRESSURE
      RP(1)=RP(1)+RPR(1)
      DO 320 I=2,NTAU
      RPE(I)=RPE(I)-RPR(I)
      DO 301 J=1,NTAU
      PET(I,J)=PET(I,J)-PRT(I,J)
301   CONTINUE
C
C TURBULENT PRESSURE
      PET(I,I)=PET(I,I)-PTT(I)
      PEPE(I,I)=PEPE(I,I)-PTPE(I)
310   CONTINUE
C
C CONVECTIVE EFFICIENCY
      DT(I,I)=DT(I,I)-.5*DG(I)*GT(I)
      DT(I,I-1)=DT(I,I-1)-.5*DG(I)*GT(I)
      DPE(I,I)=DPE(I,I)-.5*DG(I)*GPE(I)
      DPE(I,I-1)=DPE(I,I-1)-.5*DG(I)*GPE(I)
      DV(I)=-DG(I)*GV(I)
C
C TOTAL PRESSURE
      DPI=-DP(I+NTAU-1)
      DPE(I,I-1)=DPE(I,I-1)-DPI*PPPE(I+NTAU-1)
      DT(I,I-1)=DT(I,I-1)-DPI*PPTT(I+NTAU-1)
      DPE(I,I)=DPE(I,I)-DPI*PPPE(I)
      DT(I,I)=DT(I,I)-DPI*PPTT(I)
      RD(I)=RD(I)-DPI*RP(I)
      DP(I)=DP(I)-DPI
320   CONTINUE
C
      RPI=0.
      DV(1)=0.
      DO 321 I=1,NTAU
C WE HAVE THE MATRICES       PPP       PPPE  PPTT
C AND                        PEPP      PEPE  PET
C WE WANT TO SUBTRACT FROM PEPE AND PET THE PRODUCTS OF PEPP*PPP-INVERS WITH
C PPPE AND PET RESPECTIVELY. PPP IS BIDIAGONAL WITH UNITY ON THE DIAGONAL AND
C MINUS UNITY ON THE SUBDIAGONAL. ITS INVERS IS A MATRIX WITH UNITY EVERYWHERE
C UNDER AND ON THE DIAGONAL. PEPP IS MINUS UNITY. THUS THE PRODUCT OF PEPP*PPP-
C INVERS WITH PPPE IS MATRIX OF THE FOLLOWING TYPE. IN EVERY COLUMN EACH
C ELEMENT IS THE SUM OF ALL ELEMENTS ABOVE THAT POINT (AND INCLUDING) IN THE
C PPPE MATRIX. SIMILARILY FOR THE PPTT-MATRIX.
      RPI=RPI+RP(I)
      RPE(I)=RPE(I)+RPI
      RD(I)=RD(I)-DP(I)*RPI
      Y=PPPE(I)
      YY=PPTT(I)
      YY=YY+PRT(1,I)
      JMIN=MAX0(I,1)
      DO 321 J=JMIN,NTAU
      IF(J.NE.I+1) GO TO 322
      Y=Y+PPPE(I+NTAU)
      YY=YY+PPTT(I+NTAU)
322   CONTINUE
      PEPE(J,I)=PEPE(J,I)+Y
      PET(J,I)=PET(J,I)+YY
      DPE(J,I)=DPE(J,I)-DP(J)*Y
      DT(J,I)=DT(J,I)-DP(J)*YY
321   CONTINUE
C
C GRADIENT DIFFERENCE
      DO 350 I=2,NTAU
      DVI=DV(I)
      VDI=VD(I)
      FCDI=FCD(I)
      RFC(I)=RFC(I)-FCDI*RD(I)
      RV(I)=RV(I)-VDI*RD(I)
      VVVI=MAX(0.5D+0,1.0D+0-VDI*DVI)
      FCV(I)=FCV(I)-FCDI*DVI
      DO 340 J=1,NTAU
      VT(I,J)=VT(I,J)-VDI*DT(I,J)
      VPE(I,J)=VPE(I,J)-VDI*DPE(I,J)
      FCT(I,J)=FCT(I,J)-FCDI*DT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCDI*DPE(I,J)
340   CONTINUE
C
C TURBULENT VELOCITY
      PEVI=-PTV(I)
      RV(I)=RV(I)/VVVI
      FCVI=FCV(I)
      RFC(I)=RFC(I)-FCVI*RV(I)
      RPE(I)=RPE(I)-PEVI*RV(I)
      DO 350 J=1,NTAU
      VT(I,J)=VT(I,J)/VVVI
      VPE(I,J)=VPE(I,J)/VVVI
      FCT(I,J)=FCT(I,J)-FCVI*VT(I,J)
      FCPE(I,J)=FCPE(I,J)-FCVI*VPE(I,J)
      PET(I,J)=PET(I,J)-PEVI*VT(I,J)
      PEPE(I,J)=PEPE(I,J)-PEVI*VPE(I,J)
350   CONTINUE
C
C TIME
      CALL CLOCK
C
C CONVECTIVE FLUX
      DO 360 I=1,NTAU
      RT(I)=RT(I)-RFC(I)
      DO 361 J=1,NTAU
      TTT(I,J)=TTT(I,J)-FCT(I,J)
      TPE(I,J)=TPE(I,J)-FCPE(I,J)
361   CONTINUE
      IF(I.GT.MIHAL) GO TO 360
      RT(I)=RT(I)+RFC(I+1)
      DO 362 J=1,NTAU
      TTT(I,J)=TTT(I,J)+FCT(I+1,J)
      TPE(I,J)=TPE(I,J)+FCPE(I+1,J)
362   CONTINUE
360   CONTINUE
C
C ELECTRON PRESSURE
      CALL MATINV(PEPE,NTAU)
      CALL MULT(PET,PEPE,PET,SCRATC,NTAU,NTAU)
      CALL MULT(RPE,PEPE,RPE,SCRATC,NTAU,1)
      DO 374 I=1,NTAU
      DO 374 J=1,NTAU
      SUMA=0.
      DO 375 L=1,NTAU
      SUMA=SUMA+TPE(I,L)*PET(L,J)
375   CONTINUE
      TTT(I,J)=TTT(I,J)-SUMA
      RT(I)=RT(I)-TPE(I,J)*RPE(J)
374   CONTINUE
C
C TIME
      CALL CLOCK
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C BACKSUBSTITUTION IN GAUSS ELIMINATION SCHEME.
C
C INITIATE
      CALL MATINV(TTT,NTAU)
      DO 400 I=1,NTAU
      PE(I)=RPE(I)
      FC(I)=RFC(I)
      V(I)=RV(I)
      PT(I)=0.
      PR(I)=RPR(I)
      P(I)=RP(I)
C
C SOLVE FOR TEMPERATURE CORRECTION
      T(I)=0.
      DO 400 J=1,NTAU
      T(I)=T(I)+TTT(I,J)*RT(J)
400   CONTINUE
      
!      WRITE(7,68) ITER
!      WRITE(7,67) (T(I),I=1,NTAU)
C---
      TCORMX=ABS(T(1))
      DO 405 I=2,NTAU
 405  TCORMX=MAX(TCORMX,ABS(T(I)))
      write(*,*) TCORMX,TCONV
      IF(TCORMX.LE. tconv)  ITSTOP=.TRUE.
      
    
      !PRINT406, TCORMX,ITER
      inquire(file="tcormx.dat", exist=exist)
      if (exist) then
      open(unit=987,file= 'tcormx.dat',status='replace')
      write(987, *) TCORMXM
      close(987)
      else 
      open(unit=987,file= 'tcormx.dat',status='new')
      write(987, *) TCORMXM
      close(987)
      end if
      !print*, "max corr wanted ", tcormx
406   FORMAT(' Max corr. to T wanted was',F6.1,' K for iteration',I2)
C
C CHECK T CORR
C
      IF(KORT.LE.0) GO TO 412 !no limiting in temperature correction
C
C LIMITING T CORRECTION TO ONE TENTH OF THE LOCAL TEMPERATURE
C

      IF(KORT.EQ.1) THEN
      DO 401 I=1,NTAU
      T(I)=T(I)/SQRT(1.+25.*(T(I)/TT(I))**2)
401   CONTINUE
      TCORMX=ABS(T(1))
      DO 407 I=2,NTAU
      PM = T(I)/ABS(T(I))
 407  TCORMX=MAX(TCORMX,ABS(T(I)))
      PRINT4061, TCORMX
4061  FORMAT(' Max corr. to T wanted for kort=1 was',F6.1)
      END IF


C      DO 401 I=1,NTAU
C       T(I)=T(I)/SQRT(1.+100.*(T(I)/TT(I))**2)
C 401  CONTINUE
C
      IF(KORT.LE.1) GO TO 412 !no temperature inversion fix in corr
C
C Make sure correction doesnt impose a temperature inversion 
C (UGJ/900614)
C
      if(kort.eq.2) then
      DO 402 I=2,NTAU
       IF (TT(I)+T(I).LE.TT(I-1)+T(I-1)) THEN
         IF (I.EQ.2) THEN
           T(I)=TT(I-1)+T(I-1)-TT(I)+2.
          ELSE
           T(I)=1.1*(TT(I-1)+T(I-1))-0.1*(TT(I-2)+T(I-2))-TT(I)
         END IF
       END IF

 402  CONTINUE
      
C
C      IF(KORT.LE.2) GO TO 412 !no .le.-limitits in temperature correction
C
C Avoid growing oscillations in temperature correction
C (UGJ/900614)
C

c
      IF (NFIRST.eq.0) THEN
         DO 403 I=1,NTAU
403      TLAST(I)=T(I)
         NFIRST=1
      END IF
      DO 404 I=1,NTAU
      TAVDIF=0.
      IA(1)=MAX(1,I-2)
      IA(2)=MAX(2,I-1)
      IA(3)=MAX(3,I)
      IA(4)=MIN(MAX(4,I+1),NTAU)
      IA(5)=MIN(MAX(5,I+2),NTAU)
      DO 409 KIA=1,5
409   TAVDIF=TAVDIF+ABS( T(IA(KIA)) )
      TVD(I)=TAVDIF*0.2
      IF( T(I)/TLAST(I).LE.-1.0 .AND. ABS(T(I)).GT.TAVDIF ) THEN
        T(I)=-TLAST(I)
      END IF
404   CONTINUE
      endif
C
C
C      IF(KORT.LE.3) GO TO 412 ! damping in temperature correction
C
C Damp oscillations in temperature correction
C (UGJ/900719)
C
C      DO 407 I=1,NTAU
C      IF ( T(I)/TLAST(I).LE.-0.6  .AND. ABS(T(I)).GT.TVD(I) ) 
C     -     T(I)=-0.3*TLAST(I)
C407   CONTINUE
C
C
C      IF(KORT.LE.4) GO TO 412 ! damping in temperature correction
C

      IF (KORT.EQ.4.AND.IT.GE.4  .OR. KORT.EQ.5) THEN
      PPK = TDIFF
      DO 4011 I=1,NTAU
      if (abs(T(I)).GT.abs(PPK))
     - T(I)=PPK*t(i)/abs(t(i))
4011  CONTINUE
C
      TCORMX=ABS(T(1))
      DO 4071 I=2,NTAU
      PM = T(I)/ABS(T(I))
      TCORMX=MAX(TCORMX,ABS(T(I)))
4071  TCORMXM=PM*TCORMX
      PRINT408, TCORMX,KORT
      

408   FORMAT(' Maximum correction applied was',F6.1,
     *     ' for applied kort =',I2)
      END IF



C Limit temperature correction to TDIFF
C (AB/950519)
C

      
      PPK = TDIFF
      DO  I=1,NTAU
      IF (ABS(T(I)).GT.ABS(PPK)) 
     -     T(I)=PPK*T(I)/ABS(T(I))
      end do
      
C 
C
412   continue
C
      TCORMX=ABS(T(1))
      DO 414 I=2,NTAU
414   TCORMX=MAX(TCORMX,ABS(T(I)))
      !PRINT415, TCORMX,ITER
415   FORMAT(' Max corr. applied to T was',F6.1,' K for iteration ',I3)
C
C
      DO 4081 I=1,NTAU
4081   TLAST(I)=T(I)
C
C
     
!      WRITE(7,*) ' applied corrections: '
!      WRITE(7,67) (T(I),I=1,NTAU)
C
C SUBTRACT TEMPERATURE
      DO 410 I=1,NTAU
      PT(I)=PT(I)-PTT(I)*T(I)
      P(I)=P(I)-PPTT(I)*T(I)
      IF(I.GT.1) P(I)=P(I)-PPTT(I+NTAU-1)*T(I-1)
      DO 410 J=1,NTAU
      V(I)=V(I)-VT(I,J)*T(J)
      PE(I)=PE(I)-PET(I,J)*T(J)
      FC(I)=FC(I)-FCT(I,J)*T(J)
      PR(I)=PR(I)-PRT(I,J)*T(J)
410   CONTINUE
C
C SUBTRACT ELECTRON PRESSURE
      DO 420 I=1,NTAU
      PT(I)=PT(I)-PTPE(I)*PE(I)
      P(I)=P(I)-PPPE(I)*PE(I)
      IF(I.GT.1) P(I)=P(I)-PPPE(I+NTAU-1)*PE(I-1)
      DO 420 J=1,NTAU
      V(I)=V(I)-VPE(I,J)*PE(J)
      FC(I)=FC(I)-FCPE(I,J)*PE(J)
420   CONTINUE
C
C TOTAL PRESSURE
      P(1)=P(1)+PR(1)-RPR(1)
      DO 425 I=2,NTAU
      P(I)=P(I)+P(I-1)
425   CONTINUE
C
C SUBTRACT VELOCITY
      DO 430 I=1,NTAU
      PT(I)=PT(I)-PTV(I)*V(I)
430   CONTINUE
C

      DO 440 I=1,NTAU
C
C SOLVE FOR CONVECTIVE EFFICIENCY AND GRADIENT DIFFERENCE
      I1=MAX0(I-1,1)
      GI=-.5*(GT(I)*(T(I)+T(I1))+GPE(I)*(PE(I)+PE(I1)))-GV(I)*V(I)
      GG(I)=GG(I)+GI
      D(I)=D(I)-DG(I)*GI-DPS(I)*P(I)-DPES(I)*PE(I)-DTS(I)*T(I)
      IF(I.GT.1) D(I)=D(I)-DPS(I+NTAU-1)*P(I-1)-DPES(I+NTAU-1)
     &*PE(I-1)-DTS(I+NTAU-1)*T(I-1)
      
440   CONTINUE
C
C TIME
      CALL CLOCK
C
C APPLY CORRECTIONS
      IPRESS=0
      DO 450 I=1,NTAU
      TT(I)=TT(I)+T(I)
      VV(I)=MAX(VV(I)+V(I),0.0D+0)
      FFC(I)=FFC(I)+FC(I)
      PPT(I)=MAX(PPT(I)+PT(I),0.0D+0)
      PPT(I)=MIN(PPT(I),0.5*PP(I))
      PPR(I)=PPR(I)+PR(I)
      DD(I)=DD(I)+D(I)
C
C IF TOO VIOLENT CHANGES TO PPE OR PP, SET IPRESS FOR AN EXTRA
C PRESSURE INTEGRATION AFTER CORRECTIONS HAVE BEEN APPLIED.
      IF(ABS(P(I)/PP(I)).LT.0.5.AND.ABS(PE(I)/PPE(I)).LT.0.5) GOTO 451
      IPRESS=1
      GOTO 450
C
451   PPE(I)=PPE(I)+PE(I)
      PP(I)=PP(I)+P(I)
450   CONTINUE
C
C WE MUST NOT USE TRYCK TOWARDS THE END OF ITERATIONS, BECAUSE ROSS(K)
C IS DEFINED BY OPAC CALLS, RATHER THAN ROSSOP CALLS.
      IF(IPRESS.EQ.1) CALL TRYCK_sph
C
C PRINT PRESENT STATE OF THE ATMOSPHERE
*************************************************************************
      ENTRY PRESNT_sph
*
      IF(.NOT.PFE) RETURN
      INORD=IEDIT+10*IVERS
      FNORD=.1*INORD
      
      I=0
      Y=0.
      TT0=(TAU(2)*TT(1)-TAU(1)*TT(2))/(TAU(2)-TAU(1))
      
C
C TIME
      CALL CLOCK
      RETURN
C
C FORMATS
45    FORMAT(' TIME',I6,' MSEC')
48    FORMAT('1SCMARCS',F5.1,5X,'SOLVE/SPH(53)  7-NOV-80',5X,
     & '.....................................ITERATION',I3,
     & 5X,6A4/)
49    FORMAT(13('1234567890'),'123')
50    FORMAT(4(/2X,1P10E13.5))
51    FORMAT(I3,1P10E12.4,I4)
52    FORMAT(T7,'TAU',T19,'PRAD',T31,'PTURB',T43,'PTOT',T55,'GAMMA',
     *T67,'DELTA',T79,'VCONV',T91,'FCONV',T103,'PE',T115,'TEMP')
53    FORMAT(' STATE OF MODEL ATMOSPHERE')
54    FORMAT(' *=R.H. SIDES     *',23X,'*',23X,'*',
     ,11X,'*',11X,'*',11X,'*',11X,'*')
55    FORMAT(6X,'TAU',9X,'HSCALE',6X,'ADIA',8X,'GRAD',8X,'CP',10X,
     *'Q',11X,'PG',10X,'RO',10X,'PGPE',8X,'PGT')
56    FORMAT(' CORRECTIONS')
57    FORMAT(' THERMODYNAMICALS')
58    FORMAT(I4,F9.1,F8.1,1P2E10.3,0PF6.3,1PE10.3,0PF8.3,
     * F7.0,4F5.1,1PE10.3,2(0PF6.1),2I4,F6.3)
59    FORMAT('   J',4X,'WAVEL',2X,'WEIGHT',2X,'F(WAVEL)',3X,
     * 'CONTRIB',1X,'WAVEN',2X,'F(WAVEN)',4X,'MAGN',2X,'TRAD',3X,
     * 'X01',2X,'X25',2X,'S01',2X,'S25',3X,'CONTRIB',2X,
     * 'DX01  DX25 OPAC TRAN ALGB')
60    FORMAT('0SAVED VALUES FROM LOGICAL UNIT',I3)
61    FORMAT(' STARTING VALUES, QTEMP=',F5.2)
62    FORMAT(' PRESSURE EQUATION')
63    FORMAT(6X,'TAU',9X,'PTAU',8X,'ROSS',8X,'ROSSPE',6X,'ROSST',7X,
     &'KAP5',8X,'ROSSP',7X,'TAU5',8X,'TAUP',8X,'RADIUS')
64    FORMAT(' K=',I2,4X,'NEW V =',1PE10.3)
65    FORMAT('0TOTAL FLUX=',1PE11.4,' ERGS/CM**2/S    FLUX/PI=',1PE11.4,
     & ' ERGS/CM**2/S    TEFF=',0PF6.0,' K')
67    FORMAT(1X,10F7.1)
68    FORMAT(' ITERATION',I3,'     COMPUTED CORRECTION')
69    FORMAT(' ITERATION',I3,'     APPLIED  CORRECTION')
      END
C
      SUBROUTINE TRANEQ_sph
      implicit real*8 (a-h,o-z)
C
C TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
C FEATURES:
C
C 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
C    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
C    A COMPLICATED (MMU ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
C    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
C    THE ERROR.
C 2. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
C    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
C    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
C    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
C    MENT.
C 3. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
C    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
C    TERMS AD AND BD IN SUBROUTINE TRANFR.
C 4. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
C    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
C    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
C    SHOULD BE USED.
C
C THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQS.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
C
      PARAMETER (ITMAX=12)
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU   
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) 
      COMMON /RHOC/RHO(NDP)
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE 
      COMMON /TRDBUG/IDEBUG
      LOGICAL DEBUG
      DIMENSION A(ITMAX)
      DATA DEBUG/.FALSE./
C
C INITIATE, XJ IS SET TO THE DIFFUSION LIMIT VALUE
109   IDEBUG=0
      DO 100 K=1,JTAU
      IF (K.GT.1) GO TO 101
      DTAUB=(TAU(2)-TAU(1))*0.5*(X(2)+S(2)+X(1)+S(1))
      DBPLB=(BPLAN(2)-BPLAN(1))/DTAUB
      D2BPL=0.
      GO TO 102
101   IF (K.EQ.JTAU) GO TO 102
      DTAUA=DTAUB
      DTAUB=(TAU(K+1)-TAU(K))*0.5*(X(K)+S(K)+X(K+1)+S(K+1))
      DBPLA=DBPLB
      DBPLB=(BPLAN(K+1)-BPLAN(K))/DTAUB
      DTAUC=0.5*(DTAUA+DTAUB)
      D2BPL=(DBPLB-DBPLA)/DTAUC
102   XH(K)=D2BPL
      XJ(K)=BPLAN(K)+0.333333*(X(K)+S(K))/X(K)*D2BPL
      XK(K)=0.333333*BPLAN(K)+(0.2+0.111111*S(K)/X(K))*D2BPL
      FJ(K)=1.
100   SOURCE(K)=BPLAN(K)

C CALCULATE THE MATRIX ELEMENTS
      CALL TRRAYS_sph
      CALL TRANFR_sph
      CALL FORMAL_sph
      NIMP1=NIMPAC+1
C_temp      IF (DEBUG) PRINT 132,XJ,SOURCE,ERROR,FJ
C_temp     & ,((PFEAU(I,K),K=1,NDP),I=1,NIMP1)
      IF (IDEBUG.GT.1) GO TO 150
C
C ITERATION LOOP
      DO 110 IT=1,ITMAX
110   A(IT)=0.
      DO 140 IT=1,ITMAX
      ITM=IT
C
C SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
      CALL TRANSC_sph
      CALL SCATTR_sph
C_temp      IF (DEBUG) PRINT 122,EX(ISCAT),DUM,P,DTAUS
C_temp 122   FORMAT(' EX,SP1,SP2,SP3,P,DTAUS=',E10.4/(/4(1X,1P,10E12.4/)))
C
C CORRECTION TO THE SOURCE FUNCTION
      DO 120 K=1,JTAU1
      P(K)=ERROR(K)+P(K)*FJ(K)*S(K)/(X(K)+S(K))
      A(IT)=MAX(A(IT),ABS(P(K)/SOURCE(K)))
120   CONTINUE
C
C CHECK ERROR IN SOURCE FUNCTION
      IF (A(IT).LT.0.001) GO TO 141
      DO 130 K=1,JTAU1
130   SOURCE(K)=SOURCE(K)+P(K)
C
C SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
      CALL FORMAL_sph
      NTAU=KIMPAC(ISCAT)
C
C NOTE THAT FJ() SHOULD ONLY BE PICKED UP ABOVE JTAU0.  THE ISCAT
C BECOMES TO INCLINED BELOW JTAU0.
      DO 131 K=1,JTAU0
131   FJ(K)=XJ(K)/PFEAU(ISCAT,K)

132   FORMAT(' XJ,SO,ERR,FJ,PF='/(/4(1X,1P,10E12.4/)))
      IF (IDEBUG.GT.1) GO TO 150
C
C END OF ITERATION LOOP
140   CONTINUE
C
C NOT CONVERGED
      IDEBUG=1

C
C CONVERGED, IF IN FIRST ITERATION, HAVE TO CALCULATE FJ().
141   IF (ITM.GT.1) GO TO 143
      NTAU=KIMPAC(ISCAT)
      DO 144 K=1,NTAU
144   FJ(K)=XJ(K)/PFEAU(ISCAT,K)
143   CONTINUE
C
C CALCULATE MOMENTS, AND CHECK DEBUG CONTROL

      CALL TRMOM_sph

      IF (DEBUG.AND.IDEBUG.GT.1) STOP ' stop in traneq_sph at 150 '
150   continue
      IF (DEBUG.AND.IDEBUG.EQ.1) IDEBUG=0
      DEBUG=IDEBUG.GT.1
      IF (DEBUG) GO TO 109
C
      RETURN
      END
C
      SUBROUTINE TRANFR_sph
      implicit real*8 (a-h,o-z)
C
C FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE'.
C 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
C SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
C OF THE PROBLEM. INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
C 80.08.05 *NORD*
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP),
     &  FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE 
      COMMON /CSURF/HSURF,YSURF(NRAYS)
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) /RHOC/RHO(NDP)
      COMMON /TRDBUG/IDEBUG
C
      DIMENSION FUN(NRAYS),DER(NRAYS*2),DMU(NRAYS+1)
      DIMENSION TAULOG(NDP)
C
C MU LOOP
      DO 131 I=1,NIMPAC

      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C CALCULATE DTAUT ALONG THE RAY
      ZOLD=0.0
      DO 100 K=1,NTAU
      Z=SQRT(RR(NTAU-K+1)**2-PIMPAC(I)**2)
      XMU(I,NTAU-K+1)=-Z/RR(NTAU-K+1)
      IF (K.EQ.1) GO TO 100
      DZ=Z-ZOLD
      DZDR=DZ/(RR(NTAU-K+1)-RR(NTAU-K+2))
      DTAUT(NTAU-K+2)=DZDR*0.5*(X(NTAU-K+1)+S(NTAU-K+1)
     & +X(NTAU-K+2)+S(NTAU-K+2))*(TAU(NTAU-K+2)-TAU(NTAU-K+1))
100   ZOLD=Z
      TAUT(1)=DZDR*(X(1)+S(1))*TAU(1)
      DO 101 K=2,NTAU
101   TAUT(K)=TAUT(K-1)+DTAUT(K)
C

C SAVE THE TAU SCALE FOR THE RADIAL RAY (I=1).  THIS IS USED IN
C THE EXTRA/INTERPOLATION OF THE PFEAU FOR MU=0., FURTHER DOWN.
C
      IF (I.EQ.1) THEN
        DO 102 K=1,NTAU
102       TAULOG(K)=log(TAUT(K))
      ENDIF
C
C K=1
      A=1./DTAUT(2)
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
      EX(I)=TAUT(1)*(1.-0.5*TAUT(1)*(1.-0.333333*TAUT(1)))
      IF (TAUT(1).GT.0.1) EX(I)=1.-EXP(-TAUT(1))
      SP2(1,I)=SP2(1,I)/(1.+2.*A*EX(I))
      SP3(1,I)=SP3(1,I)/(1.+2.*A*EX(I))
C
C K=2,NTAU-1
      DO 110 K=2,NTAU1
      DTAUC=0.5*(DTAUT(K)+DTAUT(K+1))
C...      AD(K,I)=0.1666667*DTAUT(K)/DTAUC
C...      BD(K,I)=0.1666667*DTAUT(K+1)/DTAUC
      AD(K,I)=0.
      BD(K,I)=0.
      SP1(K,I)=-1./(DTAUT(K)*DTAUC)+AD(K,I)
      SP2(K,I)=1.
110   SP3(K,I)=-1./(DTAUT(K+1)*DTAUC)+BD(K,I)
C
C K=NTAU
      AD(NTAU,I)=0.0
      BD(NTAU,I)=0.0
      SP1(NTAU,I)=0.0
      SP2(NTAU,I)=1.0
      SP3(NTAU,I)=0.0
      IF (I.LE.NCORE) GO TO 120
      AD(NTAU,I)=0.3333333
      SP1(NTAU,I)=0.3333333-2./DTAUT(NTAU)**2
120   CONTINUE
C
C ELIMINATE SUBDIAGONAL
      DO 130 K=1,NTAU1
      SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
      SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
130   SP2(K,I)=SP2(K,I)-SP3(K,I)
131   CONTINUE
C
C FIND A GOOD RAY FOR TRANSC.  THE XMU VALUES ARE INCREASING FROM
C -1.0 TOWARDS 0.0, AND WE THUS ALWAYS FIND A VALUE FOR ISCAT.
C
C      IMAX=MIN0(MMU(JTAU0),NCORE)
C
C NCORE MUST BE SMALLER THAN MMU(K) FOR ALL K, SO EFFECTIVELY THE
C MIN0(MMU(),NCORE) IS EQUAL TO NCORE.  SHOULD ONE REQUIRE ISCAT TO
C BE A CORE RAY?  PROBABLY NOT, SINCE WITH A LARGE TAUM (ALLOWED),
C THE ISCAT RAY IS FORCED TOWARDS SMALLER IMPACT PARAMETERS, WHICH
C DEGRADES PERFORMANCE.
C
      IMAX=MMU(JTAU0)
      TMP1=1.
      DO 132 I=1,IMAX
C
C      IF (XMU(I,JTAU0).LT.-0.577) ISCAT=I
C
C  Better to look for the ray which is closest to the Eddington angle
C
	TMP2=ABS(XMU(I,JTAU0)+0.577)
	IF (TMP2.LT.TMP1) THEN
	  TMP1=TMP2
	  ISCAT=I
	ENDIF
132   CONTINUE

      RETURN
C---------------------------------------------------------------------
C
      ENTRY FORMAL_sph
C
C MU LOOP
      DO 170 I=1,NIMPAC
      NTAU=KIMPAC(I)
      NTAU1=NTAU-1
C
C INITIATE
      P(1)=SOURCE(1)
      DO 140 K=2,NTAU
140   P(K)=(1.-AD(K,I)-BD(K,I))*SOURCE(K)+AD(K,I)*SOURCE(K-1)+
     & BD(K,I)*SOURCE(K+1)
      IF(I.LE.NCORE) P(JTAU1)=SOURCE(JTAU1)+XMU(I,JTAU1)**2*XH(JTAU1)
C
C ACCUMULATE RIGHT HAND SIDE
      DO 150 K=1,NTAU1
150   P(K+1)=P(K+1)+SP1(K,I)*P(K)
C
C BACKSUBSTITUTE
      PFEAU(I,NTAU)=P(NTAU)/SP2(NTAU,I)
      DO 160 K=1,NTAU1
      PFEAU(I,NTAU-K)=(P(NTAU-K)-
     & SP3(NTAU-K,I)*PFEAU(I,NTAU-K+1))/SP2(NTAU-K,I)
      IF (PFEAU(I,NTAU-K).LE.0.0) GO TO 230
160   CONTINUE
C
C END MU LOOP
      YSURF(I)=2.*(1.-EX(I))*PFEAU(I,1)+EX(I)**2*SOURCE(1)
      FUN(I)=-XMU(I,1)*(PFEAU(I,1)-SOURCE(1)*EX(I))
      IF (YSURF(I).LE.0.0) GO TO 231
170   CONTINUE

C
      IF (NCORE.EQ.NIMPAC) THEN
C
C  TEMPORARY SECURITY TO AVOID DIVIDE BY ZERO ERRORS IN COMPUTATION OF PX,
C  WHEN THE OPACITY IS SO LARGE THAT THERE ARE ONLY NCORE RAYS.
C  IN THAT CASE, ONE USES THE OLD INTERPOLATION ROUTINE
C

      DO 1810 K=1,JTAU1
        II=MMU(K)
        IF (KIMPAC(II).EQ.K) GO TO 1810
        PX=-XMU(II-2,K)/(XMU(II-1,K)-XMU(II-2,K))
        QX=1.-PX
        PFEAU(II,K)=EXP(log(PFEAU(II-2,K))*QX+log(PFEAU(II-1,K))*PX)
      
1810   CONTINUE

      ELSE
C
C HERE IS THE NORMAL ONE.  EXTRAPOLATE FOR K<4, THEN INTERPOLATE
C
      I=NIMPAC
      DO 181 K=1,JTAU1
        IF (K.NE.KIMPAC(I)) THEN
          PX=(TAULOG(K)-TAULOG(KIMPAC(I)))/
     &       (TAULOG(KIMPAC(I-1))-TAULOG(KIMPAC(I)))
          QX=1.-PX

          IF ( (PFEAU(MMU(KIMPAC(I)),KIMPAC(I)).LE.1.E-20) .OR.
     &    (PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)).LE.1.E-20) ) THEN

             PFEAU(MMU(KIMPAC(I)),KIMPAC(I))=
     &                MAX(1.0D-99,PFEAU(MMU(KIMPAC(I)),KIMPAC(I)))
             PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1))=
     &                MAX(1.0D-99,PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)))
          END IF
          PFEAU(MMU(K),K)=log( PFEAU(MMU(KIMPAC(I)),KIMPAC(I)) )
     &       *QX + log( PFEAU(MMU(KIMPAC(I-1)),KIMPAC(I-1)) )*PX
          IF(PFEAU(MMU(K),K).GE.85.) THEN

1812         FORMAT(5I3,1P7E9.2)
             PFEAU(MMU(K),K)=85.
          END IF
          PFEAU(MMU(K),K)=EXP( PFEAU(MMU(K),K) )
        ENDIF
        IF (K.EQ.(KIMPAC(I-1)-1)) THEN
          I=I-1
        ENDIF
181   CONTINUE

      ENDIF

C
C CALCULATE MEAN INTENSITY
      DO 190 K=1,JTAU1
      XJ(K)=0.
      NMU=MMU(K)
      DO 191 I=2,NMU
      DMU(I)=XMU(I,K)-XMU(I-1,K)
      DER(I)=(PFEAU(I,K)-PFEAU(I-1,K))/DMU(I)
191   XJ(K)=XJ(K)+DMU(I)*(PFEAU(I,K)+PFEAU(I-1,K))
      XJ(K)=XJ(K)*6.
      NMU1=NMU-1
      DO 192 I=2,NMU1
192   XJ(K)=XJ(K)+(DMU(I+1)-DMU(I))*(DMU(I)*DER(I+1)+DMU(I+1)*DER(I))
      XJ(K)=XJ(K)+DMU(2)**2*DER(2)-DMU(NMU)**2*DER(NMU)
      XJ(K)=XJ(K)*0.083333333
      ERROR(K)=(XJ(K)*S(K)+BPLAN(K)*X(K))/(X(K)+S(K))-SOURCE(K)
190   CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY TRMOM_sph
C
C FLUX AT TAU(1)

      XH(1)=TRQUAD_sph(MMU(1),XMU,FUN,DER)
C
C  CALCULATE SECOND MOMENT XK
C
      DO 201 K=1,JTAU1
      NMU=MMU(K)
      DO 200 I=1,NMU
200   FUN(I)=PFEAU(I,K)*XMU(I,K)**2

201   XK(K)=TRQUAD_sph(NMU,XMU(1,K),FUN,DER)
C
C CALCULATE FIRST MOMENT, XH, FROM MOMENT RELATION
      DO 211 K=JTAU1+1,JTAU
211   XH(K)=( XK(K)-XK(K-1)+
     &       (XJ(K)+XJ(K-1)-3.*(XK(K)+XK(K-1)))*
     &       (RR(K-1)-RR(K))/(RR(K)+RR(K-1))    )*2.0/
     &          ( (TAU(K)-TAU(K-1)) * (X(K)+S(K)+X(K-1)+S(K-1)) )
C
C CALCULATE FIRST MOMENT BY USING R=DP/DTAU.  THIS IS MORE ACCURATE
C IN THE OPTICALLY THIN PARTS.
      ZOLD=0.0
      DO 212 K=2,JTAU1
        NMU=MMU(K)
        DO 213 I=1,NMU-1
          DZDR=(XMU(I,K)*RR(K)-XMU(I,K-1)*RR(K-1))/(RR(K-1)-RR(K))
          DTAU=DZDR*0.5*(X(K-1)+S(K-1)+X(K)+S(K))*(TAU(K)-TAU(K-1))
          DMU(I)=
     &      -(XMU(I,K)*RR(K)+XMU(I,K-1)*RR(K-1))/(RR(K)+RR(K-1))
          FUN(I)=DMU(I)*(PFEAU(I,K)-PFEAU(I,K-1))/DTAU
213     CONTINUE
        FUN(NMU)=0.
        DMU(NMU)=0.

        XH(K)=-TRQUAD_sph(NMU,DMU,FUN,DER)
212   CONTINUE
C
C SURFACE FLUX
      NMU=MMU(1)
      PX=-XMU(NMU-2,1)/(XMU(NMU-1,1)-XMU(NMU-2,1))
      QX=1.-PX
      YSURF(NMU)=EXP(log(YSURF(NMU-2))*QX+log(YSURF(NMU-1))*PX)
      DO 220 I=1,NMU
220   FUN(I)=-XMU(I,1)*YSURF(I)

      HSURF=0.5*TRQUAD_sph(NMU,XMU,FUN,DER)

      RETURN
C------------------------------------------------------------------
C
C EMERGENCY EXIT
230   KK=NTAU-K
      IDEBUG=2
C      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
C232   FORMAT('0NON-POSITIVE RESULT AT I,K,J0,J1,N=',5I3)
      GO TO 233
231   KK=0
      IDEBUG=3
C      PRINT 232,I,KK,JTAU0,JTAU1,NTAU
C233   PRINT 237,NCORE,ISCAT,DIFLOG,RADIUS,EX(I)
233   continue
237   FORMAT(' NCORE,ISCAT,DIFLOG,RADIUS,EX(I)=',2I3,3G12.3)
C      WRITE (13,*) JTAU,TAU,X,S,BPLAN,RADIUS,RR,RHO,ROSS
C      PRINT 236,TAU,X,S,BPLAN,RR,RHO,ROSS,SOURCE
C     & ,(YSURF(I),(PFEAU(I,K),K=1,39),I=1,NMU)
236   FORMAT('0TAU=',4(/10E12.4)/'0X=',4(/10E12.4)/'0S=',4(/10E12.4)
     & /'0BPLAN=',4(/10E12.4)/'0RR=',4(/10E12.4)/'0RHO=',4(/10E12.4)
     & /'0ROSS=',4(/10E12.4)/'0SOURCE=',4(/10E12.4)
     & /'0YSURF,PFEAU='/(10E12.4))
      
      RETURN
      END
C
      SUBROUTINE TRANSC_sph
      implicit real*8 (a-h,o-z)
C
C SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
C IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE RAY.
C 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' GIVES THE
C ESTIMATED MEAN INTENSITY CORRECTION (FJ*P). TRANSC CALCULATES THE MATRIX
C ELEMENTS FOR SCATTR.
C 79.06.21 *NORD*
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),SP1(NDP),SP2(NDP),SP3(NDP),
     &  DUM(NDP,NRAYS,3),AD(NDP,NRAYS),BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/DIFLOG,RADIUS,RR(NDP),NCORE 
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C CALCULATE TAUS ALONG THE RAY, SPIRAL AT EDDINGTON ANGLE AT DEPTH.
      Z=SQRT(RR(JTAU0)**2-PIMPAC(ISCAT)**2)
      ZOLD=Z
      DO 101 K=2,JTAU
      IF (JTAU-K+2.GT.JTAU0) GO TO 103
      Z=SQRT(RR(JTAU-K+1)**2-PIMPAC(ISCAT)**2)
      DZ=Z-ZOLD
      DZDR=DZ/(RR(JTAU-K+1)-RR(JTAU-K+2))
      GO TO 104
103   DZDR=1.732
104   DTAUS(JTAU-K+2)=DZDR*0.5*(X(JTAU-K+1)+S(JTAU-K+1)
     & +X(JTAU-K+2)+S(JTAU-K+2))*(TAU(JTAU-K+2)-TAU(JTAU-K+1))
101   ZOLD=Z
C
      TAUS(1)=DZDR*(X(1)+S(1))*TAU(1)
C
C K=1
      A=1./DTAUS(2)
      B=A**2
      SP2(1)=1.-FJ(1)*S(1)/(X(1)+S(1))+2.*A
      SP3(1)=-2.*B
      T=TAUS(1)
      EX(ISCAT)=TAUS(1)*(1.-0.5*TAUS(1)*(1.-0.333333*TAUS(1)))
      IF (TAUS(1).GT.0.1) EX(ISCAT)=1.-EXP(-TAUS(1))
      SP2(1)=SP2(1)-2.*A*EX(ISCAT)*FJ(1)*S(1)/(X(1)+S(1))
      SP2(1)=SP2(1)/(1.+2.*A*EX(ISCAT))
      SP3(1)=SP3(1)/(1.+2.*A*EX(ISCAT))
C
C K=2,NTAU-1
      DO 100 K=2,NTAU1
      DTAUC=0.5*(DTAUS(K)+DTAUS(K+1))
      SP1(K)=-1./(DTAUS(K)*DTAUC)
      SP2(K)=1.-FJ(K)*S(K)/(X(K)+S(K))
100   SP3(K)=-1./(DTAUS(K+1)*DTAUC)
C
C K=NTAU
      SP1(NTAU)=0.0
      SP2(NTAU)=X(NTAU)/(X(NTAU)+S(NTAU))
      SP3(NTAU)=0.0
C
C ELIMINATE SUBDIAGONAL
      DO 120 K=1,NTAU1
      SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
      SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
120   SP2(K)=SP2(K)-SP3(K)
121   CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY SCATTR_sph
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C ACCUMULATE RIGHT HAND SIDE
      P(1)=ERROR(1)
      DO 150 K=1,NTAU1
150   P(K+1)=ERROR(K+1)+SP1(K)*P(K)
C
C BACKSUBSTITUTE
      P(NTAU)=P(NTAU)/SP2(NTAU)
      DO 160 K=1,NTAU1
160   P(NTAU-K)=(P(NTAU-K)-SP3(NTAU-K)*P(NTAU-K+1))/SP2(NTAU-K)
C
      RETURN
      END
C
      SUBROUTINE TRYCK_sph
      implicit real*8 (a-h,o-z)
C
C TRYCK IS A FAST PRESSURE INTEGRATION ROUITINE. IT IS FAST BECAUSE OF
C TWO REASONS: 1) IT INTEGRATES THE DIFFFERENTIAL EQUATION FOR LN(P)
C AS A FUNCTION OF LN(TAU). 2) IT ITERATES DIRECTLY ON THE ELECTRON
C PRESSURE, KEEPING THE NUMBER OF CALLS TO ABSKO TO A MINIMUM.
C ASSUMING A POWER LAW BEHAVIOUR OF PP,TT,PPE,ETC.: PP=C*TAU**DP,ETC., ONE
C CAN SHOW THAT DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST)/(1.+ROSSPE/PGPE).
C THE ANSATZ FOR PP IMPLIES PP(1)=TAU(1)*GRVR/(ROSS(1)*DP), WHICH
C SERVES AS A BOUNDARY CONDITION.
C 790516 *NORD*
C
      include 'parameter.inc'
C
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PPE(NDP),TT(NDP),TAULN(NDP),RO(NDP),
     &NTAU,ITER
      COMMON /TAUC/TAU(NDP),DLNTAU(NDP),JTAU 
      COMMON /CG/GRAV,KONSG
      COMMON /ROSSC/ROSS(NDP),CROSS(NDP)
      COMMON /CI8/PGC,RHOC,EC
      COMMON /CSPHER/TAURAT,RADIUS,RR(NDP),NCORE 
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      DATA EPS,RELT,RELPE,PEDEF/1.E-3,1.E-3,1.E-3,1./
C
C START
      MSA=0
C      WRITE(7,101)
C101   FORMAT('1PRESSURE INTEGRATION'/'  K',6X,'TAU',10X,'TT',
C     & 9X,'PPE',8X,'PTOT',8X,'ROSS',10X,'DP',9X,'NABSKO')
      DT=0.
C USE 'DLNT/DLNTAU'=DT=0. TO BE COMPATIBLE WITH SOLVE. OTHERWISE
C DT=(TT(2)/TT(1)-1.)/DLNTAU(2)
      NABSKO=0
      KK=1
      IF(PPE(1).LE.0.) PPE(1)=PEDEF
C
C ITERATE ON BOUNDARY CONDITION, USING PARTIAL DERIVATIVES
      GRVR=GRAV*(RADIUS/RR(1))**2
C test for effect of constant gravity...
      IF (KONSG.EQ.1) GRVR=GRAV   
100   CONTINUE
      KL=1
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1),1)
      PP(1)=PGC+PPT(1)+PPR(1)
      PG=PGC
      ROSST=CROSS(1)*ROSSOP(TT(1)*(1.+RELT),PPE(1),1)
      PGT=PGC
      ROSSPE=CROSS(1)*ROSSOP(TT(1),PPE(1)*(1.+RELPE),1)
      PGPE=PGC
      PGT=(PGT/PG-1.)/RELT
      PGPE=(PGPE/PG-1.)/RELPE
      ROSST=(ROSST/ROSS(1)-1.)/RELT
      ROSSPE=(ROSSPE/ROSS(1)-1.)/RELPE
      NABSKO=NABSKO+3
      DP=(1.+DT*(ROSSPE*PGT/PGPE-ROSST))/(1.+ROSSPE/PGPE)
      DP=MAX(DP,0.1D+0)
      DLNPE=log(GRVR*TAU(1)/(PG*ROSS(1)*DP))/(PGPE+ROSSPE)
      PPE(1)=PPE(1)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) GOTO 100
C
C END BOUNDARY CONDITION
      ROSS(1)=CROSS(1)*ROSSOP(TT(1),PPE(1),1)
      NABSKO=NABSKO+1
      PP(1)=PGC+PPT(1)+PPR(1)
C      WRITE(7,102) KK,TAU(1),TT(1),PPE(1),PP(1),ROSS(1),DP,NABSKO
C102   FORMAT(I3,6E12.5,I12)
C
C TAU LOOP
      DPE=(DP-DT*PGT)/PGPE
      DEDLNP=-(PGPE*PG/PP(1)+.5*DLNTAU(2)*GRVR*TAU(1)/(PP(1)*ROSS(1))*
     & (PGPE*PG/PP(1)+ROSSPE))
      DO 110 K=2,NTAU
      GRVR=GRAV*(RADIUS/RR(K))**2
C test for effect of constant gravity....
      IF (KONSG.EQ.1) GRVR=GRAV   
      PPE(K)=PPE(K-1)*EXP(DPE*DLNTAU(K))
      NABSKO=0
C
C ITERATION LOOP
      DLNPE=0.
111   CONTINUE
      KL=K
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K),k)
      PP(K)=PGC+PPT(K)+PPR(K)
      NABSKO=NABSKO+1
      ERROR=(.5*DLNTAU(K)*GRVR*(TAU(K-1)/(PP(K-1)*ROSS(K-1))+
     & TAU(K)/(PP(K)*ROSS(K)))-log(PP(K)/PP(K-1)))
C       print*,'error,dedlnp ', error,dedlnp
      CALL ZEROF(ERROR,DLNPE,DEDLNP)
      PPE(K)=PPE(K)*EXP(DLNPE)
      IF(ABS(DLNPE).GT.EPS) GOTO 111
C
C END TAU LOOP
      ROSS(K)=CROSS(K)*ROSSOP(TT(K),PPE(K),k)
      NABSKO=NABSKO+1
      PP(K)=PGC+PPT(K)+PPR(K)
      DP=GRVR*TAU(K)/(PGC*ROSS(K))
      DPE=log(PPE(K)/PPE(K-1))/DLNTAU(K)
C      WRITE(7,102) K,TAU(K),TT(K),PPE(K),PP(K),ROSS(K),DP,NABSKO
110   CONTINUE
C
C END
      MSB=0
      MSB=MSA-MSB
C     PRINT 120,MSB
120   FORMAT(' TRYCK_sph TIME=',I5,' MS')
      RETURN
      END
C
C
      FUNCTION TRQUAD_sph(N,X,F,W)
      implicit real*8 (a-h,o-z)
C
      DIMENSION X(N),F(N),W(2*N)
C was : dim x(1) etc...
C
C TRAPEZOIDAL QUADRATURE PLUS NEXT ORDER CORRECTION FOR NON-
C -EQUIDISTANT GRID
      N1=N-1
      Q=0.
      DO 100 K=2,N
      W(K)=X(K)-X(K-1)
      W(N+K)=(F(K)-F(K-1))/W(K)
100   Q=Q+W(K)*(F(K-1)+F(K))
      Q=Q*6.
      DO 101 K=2,N1

101   Q=Q+(W(K+1)-W(K))*(W(K)*W(N+K+1)+W(K+1)*W(N+K))
      W1=((W(2)+0.5*W(3))*W(N+2)-0.5*W(2)*W(N+3))*2.0/(W(2)+W(3))
      WN=((W(N)+0.5*W(N1))*W(N+N)-0.5*W(N)*W(N+N1))*2.0/(W(N)+W(N1))
      Q=0.083333333*(Q+W(2)**2*W1-W(N)**2*WN)

      TRQUAD_sph=Q
      RETURN
      END
C
C
      SUBROUTINE TRRAYS_sph
      implicit real*8 (a-h,o-z)
C
      include 'parameter.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  PIMPAC(NRAYS),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP),MMU(NDP),KIMPAC(NRAYS),NIMPAC
      COMMON /CSPHER/TLIM,RADIUS,RR(NDP),NCORE 
      COMMON /CSTYR/MIHAL,NOCONV
      COMMON /CTAUM/TAUM
*
*  Distribute rays, based on two depth indices, JTAU0 and JTAU1.
*
*  JTAU0 is the largest depth index for which TAU < SQRT((X+S)/X),
*  and represents the "surface", where radiation is released.
*  Above JTAU0, linearized perturbations in the radiation field are
*  represented by a single ray, which is chosen in TRANFR as the ray
*  which is closest to the Eddington angle 1./sqrt(3.) at JTAU0.
*
*  JTAU1 is the largest depth index for which TAU < TAUM*SQRT((X+S)/X),
*  and represents the border to the "core", where diffusion is a good
*  approximation.  The detailed formal solution of the radiative transfer
*  on a set of parallel rays is only performed above JTAU1.  TAUM should
*  be at least 50 to 100.
*
      TAUT(1)=TAU(1)*(X(1)+S(1))
      JTAU0=1
      JTAU1=1
      ICASE=1
      DO 100 K=2,JTAU
      
      TAUT(K)=TAUT(K-1)+0.5*(X(K)+S(K)+X(K-1)+S(K-1))*(TAU(K)-TAU(K-1))
      GO TO (101,102,103,104,105),ICASE
101   IF (TAUT(K).LT.1.0) GO TO 105
      ICASE=2
102   IF (TAUT(K).LT.SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=3
103   IF (TAUT(K).LT.TAUM) GO TO 105
      ICASE=4
104   IF (TAUT(K).LT.TAUM*SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=5
105   IF (ABS(XH(K))/BPLAN(K).GT.0.10) ICASE=MIN0(ICASE,4)
      GO TO (106,106,107,107,100),ICASE
106   JTAU0=K
107   JTAU1=K
100   CONTINUE
      JTAU1=MAX0(JTAU1,3)
*
*  Distribute rays in the core.  This is done in such a way that the
*  mu values at depth JTAU0 are equidistant.  The rationale for this
*  is to get a good representation of the the "core" part of the
*  radiation as a function of mu, in the optically thin regions.
*
      RRK=RR(JTAU0)
      NCORE1=NCORE-1
      IF (NCORE1.EQ.0) PRINT*,' NCORE1= ',NCORE1
      DR=(RR(JTAU0)-SQRT(RR(JTAU0)**2-RR(JTAU1)**2))/NCORE1
      DO 110 I=1,NCORE1
      PIMPAC(I)=SQRT(RR(JTAU0)**2-RRK**2)
      RRK=RRK-DR
110   KIMPAC(I)=JTAU1
C
C RAYS IN ATMOSPHERE
      I=NCORE
      KI=JTAU1
120   KIMPAC(I)=KI
      PIMPAC(I)=RR(KI)
      KIP=KI
121   KI=KI-1
      
      IF (TAU(KIP)/TAU(KI).LT.TLIM.AND.KI.GT.1) GO TO 121
      I=I+1
      IF (KI.GE.3) GO TO 120
      KIMPAC(I)=0
      PIMPAC(I)=RR(1)
      NIMPAC=I-1
      IF (NIMPAC.LT.NRAYS) GO TO 131
      PRINT 122,NIMPAC,NRAYS,NCORE,TLIM
122   FORMAT(' ** NMBR OF RAYS (NIMPAC) TOO LARGE =',I3
     &  ,' parameter NRAYS= ',I3,'  NCORE,TLIM =',
     & I3,F6.3)
      STOP ' stop in trrays_sph at 122 '
C
C FIND THE NUMBER OF MU-PNTS FOR EACH K, PLUS ONE EXTRA FOR MU=0.0
131   II=NIMPAC+1
      DO 130 K=1,JTAU1
      MMU(K)=II
      XMU(II,K)=0.0
      IF (K+1.EQ.KIMPAC(II-1)) II=II-1
130   CONTINUE
      RETURN
      END
C
C
C_ugj950523:  here the sphereical part ends@@@@@
C
      SUBROUTINE VAAGL(NLB,XL,W)
      implicit real*8 (a-h,o-z)
C
C        DENNA RUTIN BERAEKNAR VAAGLAENGDSPUNKTER OCH MOTSVARANDE
C        VIKTER.  GLAMD(I),I=1,JLBDS  AER DISKONTINUITETER ELLER DEL-
C        NINGSPUNKTER I VAAGLAENGDSLED, INKLUSIVE VAAGLAENGDSSKALANS
C        AENDPUNKTER.  MLD(I), I=1,JLBDS-1   AER DET OENSKADE ANTALET
C        VAAGLAENGDSPUNKTER I RESP. INTERVALL.   DESSA STORHETER
C        INLAESES I RUTINEN.
C        VI ANVAENDER SUBR.  G A U S I .
C        ****   OBSERVERA. JLBDS FAAR HOEGST VARA 15 MED HAER BRUKADE DI
C        OCH FORMATSATSEN 100 ******
C
      DIMENSION XL(500),W(500),GLAMD(100),MLD(100),XLP(10),WP(10)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CLINE3/GLAMD,JLBDS
      COMMON/LDOPAC/ ALES,BLES
C

      READ(IREAD,100)JLBDS,ALESX,BLESX

      IF(JLBDS.GT.1) GOTO 31
C ALLOW ONE POINT STANDARD OPACITY
      READ(IREAD,102)XL(1)

      W(1)=1.
      NLB=1
      RETURN
31    CONTINUE
      IF(JLBDS.GT.10)  GO TO 77
C SET UV OPACITY CONSTANTS
      ALES=ALESX
      BLES=BLESX
   77 CONTINUE
      JP=JLBDS-1
      DO3 K=1,JP 
   3   READ(IREAD,102)GLAMD(K),MLD(K)
      READ(IREAD,102)GLAMD(JLBDS)
      I=0
      DO2 K=1,JP
      IF(MLD(K).GT.0)GO TO 21
      JIP=2
      XLP(1)=GLAMD(K)
      WP(1)=(GLAMD(K+1)-GLAMD(K))*0.5
      XLP(2)=GLAMD(K+1)
      WP(2)=WP(1)
      MLD(K)=1
      IF(K.EQ.JP)GO TO 22
      IF(MLD(K+1).EQ.0)GO TO 23
   22 MLD(K)=2
      GO TO 23
   21 CONTINUE
      CALL GAUSI(MLD(K),GLAMD(K),GLAMD(K+1),WP,XLP)
      JIP=MLD(K)
   23 CONTINUE
      DO1 J=1,JIP
      IP=I+J
      XL(IP)=XLP(J)
      W(IP)=WP(J)+W(IP)
    1 CONTINUE
    2 I=I+MLD(K)
      NLB=I

  100 FORMAT(I5,5X,2F10.0)
  102 FORMAT(F10.5,I5)
      RETURN
      END
C
      FUNCTION VVMLT(A,B,C)
      implicit real*8 (a-h,o-z)
C
C COMPUTE SQRT(B*(A+.5/B/C-SQRT(.5/B/C*(2.*A+.5/B/C)))) WITH SPECIAL
C CONSIDERATION ON SMALL C-VALUES. 73.12.02  *NORD*
C
      D=.5/B/C
      IF(D.GT.20.*A) GO TO 1
      E=A+D-SQRT(D*(2.*A+D))
      GO TO 2
1     E=.5*MAX(A,0.0D+0)**2/D
2     VVMLT=SQRT(B*E)
      RETURN
      END
C
      DOUBLE PRECISION   FUNCTION X02AAF(X)
      implicit real*8 (a-h,o-z)
C
C  RETURNS THE SMALLEST VALUE X02AAF SUCH THAT 1+X02AAF>1
C
C CYBER ?
C     X02AAF=7.2E-15
C IBM DOUBLE PREC.
C      X02AAF = 2.2205D-16
C VAX DOUBLE PREC.
C     X02AAF=1.E-16
C
C APOLLO/DOMAIN WORKSTATION DN 3000/4000  (TRY CALCULATING IT . . . )
C
      X02AAF = 1.0D-15
C
      RETURN
          END
C
      SUBROUTINE X1MAKE(N,XXI,W,XI)
      implicit real*8 (a-h,o-z)
C
C      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NDIM=100)
      DIMENSION XI(NDIM),W(NDIM),XXI(NDIM)
C
      I=N
      DO 5 J=1,N-1
        XI(I)=XXI(I-1)
        I=I-1
    5 CONTINUE
C
      F=.5
      XI(1)=XI(2)+F*(XI(2)-XI(3))/(W(3)/W(2)-1.)
C
C
      RETURN
C
C
C
      E  N  D
C
      SUBROUTINE XIINIT(WVAL,K,XIVAL,N,XI,W)
      implicit real*8 (a-h,o-z)
C
C      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NDIM=100)
      DIMENSION XI(NDIM),W(NDIM)
      DIMENSION A(NDIM),B(NDIM),C(NDIM),ITYPE(NDIM)
C
      D1=(XI(2)-XI(1))/(W(2)-W(1))
      D2=(XI(3)-XI(2))/(W(3)-W(2))
      YPI=D1*D1/D2
      DO 105 I=2,N
        HI=W(I)-W(I-1)
        DI=(XI(I)-XI(I-1))/HI
        IF(YPI/DI.GT.1.) THEN
          ITYPE(I)=0
          C(I)=(1./DI-1./YPI)/DI/HI
          BB=1./YPI-2.*XI(I-1)*C(I)
          A(I)=W(I-1)-XI(I-1)*(BB+C(I)*XI(I-1))
          B(I)=-.5*BB/C(I)
          YPI=1./(BB+2.*C(I)*XI(I))
        ELSE
          ITYPE(I)=1
          C(I)=(DI-YPI)/HI
          B(I)=YPI-2.*W(I-1)*C(I)
          A(I)=XI(I-1)-W(I-1)*(B(I)+W(I-1)*C(I))
          YPI=B(I)+2.*C(I)*W(I)
        END IF
  105 CONTINUE
C
C
      RETURN
C
      ENTRY XIMAKE(WVAL,K,XIVAL,N,XI,W)
C
      IF(WVAL.LE.W(K).OR.K.EQ.N) GOTO 12
C
      DO 10 I=1,N-K-1
        K=K+1
        IF(WVAL.LE.W(K)) GO TO 12
   10 CONTINUE
C
      K=N

   12 IF(ITYPE(K).EQ.1) THEN
      XIVAL=A(K)+WVAL*(B(K)+C(K)*WVAL)
      ELSE
        YYY=B(K)*B(K) - (A(K)-WVAL)/C(K)
        IF(YYY.LT.0.) THEN
C        PRINT 1234,K,B(K),A(K),WVAL,C(K),YYY
C      DO 2727 II7=1,N
C      PRINT 1235,II7,W(II7),XI(II7),A(II7),B(II7),C(II7)
C2727  CONTINUE
C
1234  FORMAT('XIMAKE: K,B(K),A(K),WVAL,C(K),YYY '/I5,1P5E13.5)
1235  FORMAT(I5,F8.5,F10.5,1P3E15.5)
         STOP 'XIMAKE-KRASCH'
      ENDIF
        XIVAL=B(K)-DSQRT(B(K)*B(K)-(A(K)-WVAL)/C(K))
      END IF
C
C
      RETURN
C
C
      E  N  D
C
      SUBROUTINE XMETAL(ID,N,W,XI)
      implicit real*8 (a-h,o-z)
C
C      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'parameter.inc'
C
      DIMENSION W(6),XI(6)
      COMMON /ODFAD/ V(NDP,4)
C
      DO 5 I=1,N-2
        XI(I)=V(ID,I)
    5 CONTINUE
C
      FN=.5
      XI(N-1)=XI(N-2)+FN*((XI(N-3)-XI(N-2))*(1.-W(N-2))/
     &                          (W(N-3)-W(N-2)))
C
      RETURN
C
C
C
      E  N  D
C
      SUBROUTINE ZEROF(F,DX,DFDX)
      implicit real*8 (a-h,o-z)
C
C FIND DX=-F/DFDX, TO MAKE F ZERO. IF DX=0 AT ENTRY, THEN DFDX IS A
C START APPROXIMATION, AND IT IS THE FIRST CALL. OTHERWISE USE OLD
C INFO.  780926/NORDLUND.
C
      !write(*,*) 'F',F
      !write(*,*) 'FOLD',FOLD
      !write(*,*) 'DXOLD',DXOLD
      if (F.eq.FOLD) then
      write(*,*) "F equals F_OLD, which would make the derivative zero,
     >   currently no fix in place Sorry :( Try a different setup"
       stop
      endif  
      IF(DX.NE.0.) DFDX=(F-FOLD)/DXOLD
      !write(*,*) 'DFDX',DFDX
      DX=-F/DFDX
      !write(*,*) 'DX',DX
C
C TRY TO AVOID OVERFLOW WHEN A JUMP IN TEMPERATURE IS PRESENT 
      IF (DX.GT.2.) DX=2.
C
      FOLD=F
      DXOLD=DX
      !write(*,*) 'FOLD',FOLD
      !write(*,*) 'DXOLD',DXOLD
      RETURN
      END

      SUBROUTINE OSTABLOOK

      implicit real*8 (a-h,o-z)
       !use omp_lib
       include 'parameter.inc'
       character atnames*2, molnames*8, molnames2*4, shn*8
       character(len=100) :: mol_file
       !BCE November 2022 - routine is being revamped for correct consideration
       !of molecular opacities.


      DIMENSION OPJV(NWL,NDP),AKAPMOL(NDP)
     *      ,opav(20,nwl),wnav(nwl)
     *      ,TMOL(MTEMP),OPLN(MTEMP),FXLN(MTEMP)
     *      ,reliso(15),OPT(MTEMP),DADT(MTEMP),WT(3,MTEMP)
      dimension xx(ndp), yy (ndp)

      dimension fxi(ndp), fyi(ndp), fxyi(ndp)
      CHARACTER MOLNAME*4,OSFIL*60,MOLID*4,SAMPLING*3
      LOGICAL FIRST

      NAMELIST /INPUTOSMOL/ MOLID, KTEMP, TMOL, NWNOS
     &   ,VKMS, KISO, RELISO, RATIS, JDERIV, L_PER_STELLAR, lchrom

      COMMON /COPINF/ SUMOP(maxosmol,NDP),SUMKAP(maxosmol,NDP)
      common/carciv/ larciv    !=1 if called from arciv, otherwise = 0
      common/eostab/ xmin, ymin, dx, dy, f(mtemp,mpe)
     *       , fx(mtemp,mpe), fy(mtemp,mpe), fxy(mtemp,mpe), nx, ny
      COMMON /STATEC/PPR(NDP),PPT(NDP),PP(NDP),GG(NDP),ZZ(NDP),DD(NDP),
     & VV(NDP),FFC(NDP),PE(NDP),T(NDP),TAULN(NDP),RO(NDP),
     & NTAU,ITER
      COMMON /CMETPE/ PPEL(NDP), METPE
      !COMMON /ROSSC/XKAPR(NDP),CROSS(NDP)
      !COMMON/CXLSET/XL(20,10),NSET,NL(10)
      !COMMON /CVAAGL/XLB(500),W(500),NLB
      COMMON/COS/WNOS(NWL),CONOS(NDP,NWL),WLOS(NWL),WLSTEP(NWL)
     *    ,KOS_STEP,NWTOT,NOSMOL,NEWOSATOM,NEWOSATOMLIST
     *    ,nchrom,OSFIL(maxosmol),MOLNAME(maxosmol),SAMPLING
      !COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(33)
      COMMON/CI4/ TMOLIM,IELEM(16),ION(16,5),MOLH,JUMP
      common/ci5/abmarcs(18,ndp),anjon(18,5),h(5),part(18,5),
     *dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
      COMMON /CPOLY/FACPLY,MOLTSUJI
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),partp(ndp,0:maxmol),
     &  partpp(ndp,0:maxmol)
      common /cmasabs/ masabs(3)
      common /cosexp/ lops,nops
      common /catoms_head/ p6(9),t_at(17)
      COMMON /PJONINF/ P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP),
     & P_NEU_HE(NDP),P_ION_HE(NDP), P_NON_HHECNO(NDP), PG_JON(NDP), 
     & HN_JON(NDP), RO_JON(NDP), P6_JON(NDP)
      COMMON /CKMOL/KMOL(MAXOSMOL)

C     parameter(nspec=892)
      common /cgem/pres_gem(ndp,nspec)
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      common /cdrift/ idust, ieps, idustopac
      character name_gem*8
      common /cabink/abink(ndp,nspec)
      common /ch4/ nch4
      dimension ptot(ndp),pp_sum(ndp),pg(ndp)
      dimension trpe(ndp), trphe(ndp),akapmax(ndp)
      dimension pe_gem(ndp),ptot1(ndp),dptot(ndp)
     &  ,pe1(ndp),dptot2(ndp),dpe2(ndp),dpe(ndp)
      integer krome_photo_on,call_counter
      INTEGER MOLH, JUMP
      DATA FIRST/.TRUE./

      common /ggchemmu/ggmu(NDP),ggrho(NDP),ppsum(ndp),ppappsum(ndp),
     &   ppnonappsum(ndp),tg(ndp),pges(ndp)
     &  ,ppat1sum(ndp),ppat2sum(ndp),ppmolsum(ndp),ppgs(ndp)
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      !integer, dimension(75) :: idmarcspart, idggchempart

      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
      common /ggchemdetabs / f1_dt(ndp), f5_dt(ndp), 
     >                       rC(ndp), rMg(ndp), rAl(ndp), 
     >                       rSi(ndp), rHe(ndp), ro_dt(ndp)
      logical ggchem_mol(maxosmol), ggchem_index_read, species_found
      integer ggchem_index(maxosmol), molno
      character(len=5) molnames_new(maxosmol)
      common /molupdate/ molnames_new, 
     * ggchem_index,
     * molno, ggchem_mol, ggchem_index_read
      common /dpeset/ dpein,dtin, pe_corr(ndp)
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      character*20 file_id, file_name
      dimension nmid(20)
      data nmid/3,4,16,29,33,34,37,39,44,53,59,62,8*1/
  

      if (first) then
      TOSREAD = 0.
      TPART = 0.
      call_counter = 0
      end if

      call timex1
      do it=1,ntau
         do jv=1,nwtot
            conos(it,jv) = 0.
         end do
      end do
C
C CALCULATE THE PARTIAL PRESSURE
      TBPART=SECOND()
      MOLM=MOLH
      MOLH=0
      div = 1./(xmy(1)*xmh)
      
      do k=1,ntau
        ptot(k)=pp(k)-ppr(k)-ppt(k)+pe_corr(k)
       
      end do
      DO  K=1,NTAU
        KL=K         

C ------------ USING GGCHEM TO COMPUTE PARTIALPRESSURES ------------  
          
! use here ggrho(k) from GGchem instead ro(k) from jon and skip calling
! jon here
            call GGCHEM(k,t(k),ptot(k))
            
            ppel(k) = ppelGG
            ggmu(k) = ggmuk
            ggrho(k) = ggrhok
            ppsum(k) = ppsumk
            ppappsum(k) = ppappsumk
            ppnonappsum(k) = ppnonappsumk
            ppat1sum(k) = ppat1sumk
            ppat2sum(k) = ppat2sumk
            ppmolsum(k) = ppmolsumk
            ppgs(k) = ppgsk
            tg(k) = tgk
            pges(k) = pgesk
            ro(k) = ggrho(k)
            f1_dt(k) = f1gg
            f5_dt(k) = f5gg
            rC(k) = rCgg
            rMg(k) = rMggg
            rAl(k) = rAlgg
            rSi(k) = rSigg
            rHe(k) = rHegg
            ro_dt(k) = ggrhok
      enddo

      
C     Updating stuff for ggchem
      DO K=1,NTAU
        partpp(k,0) = ppallat(k,1)
         partpp(k,2) = (ppallmol(k,1) /3.7095D3/T(K) )**2 /RO(K)
         XNHE = abmarcs(2,k) / (XMH*XMY(k))

         partpp(k,17) = 1.38053e-16 * T(K) * XNHE
         partpp(k,17)=
     *     partpp(k,2)*partpp(k,17)/( 3.7095D3 * T(K) )**2
         partpp(k,17) = partpp(k,17) * RO(K)
         PHE = partpp(K,17)
         P6_JON(K) = xmettryck(k,1)+0.42*PHE+0.85*ppallmol(k,1)     
      end do
C ------------ USING KROME TO COMPUTE NON-EQ CHEMISTRY------------

      if (krome_on.eq.1) then
       call krome_solve(ntau,T,ptot)
      endif
C ------------ KROME DONE------------

      !MOLH=MOLM
      molh=0


C     ADS: Code that matches the ggchem index to the read in molnames
C     Idea: call ggchem once to generate ggchem index of molnames and atnames for partial pressures
C     in ggchem. We will then try to match these to those in mol_names.dat, where we could have 
C     molnames or atomnames inside. The order defined there is the order in which the opacities are
C     processed and read in. The order in mol_names.dat can be whatever.

C     Lets now search for the ggchem index of our species      
      if (.not. ggchem_index_read) then
C           Note: we only need to do this once!!       
      do nm=1, molno
          species_found = .FALSE.
C         We first try to match the species to a molecule
C         Species can be in mixed uppercase and lowercase to be matched to GGchem
          molloop: do imol = 1, 543
             shn = trim(upper(molnames_new(nm)))
             if (trim(molnames(imol)) == shn) then
                  ggchem_index(nm) = imol
                  ggchem_mol(nm) = .TRUE.
                  species_found = .TRUE.
                  exit molloop
             endif
          end do molloop
C         We then check if our species is an atom          
          atmloop: do iatm = 1, 22
             shn = trim(molnames_new(nm))
             if (trim(atnames(iatm)) == shn) then
                  ggchem_index(nm) = iatm
                  if (.not. species_found) then
                        ggchem_mol(nm) = .FALSE.
                        species_found = .TRUE.
                  else
                        print*,' Species could be ',
     &                     'atom or molecule! Check: '
     &                     ,nm,molnames_new(nm)                    
                        STOP 'Check species name in mol_names.dat'
                  endif
                  exit atmloop
             endif
          end do atmloop
          if (.not. species_found) then
            print*,' species is not an atom or molecule from ggchem: '
     &         ,nm,molnames_new(nm)
            STOP 'Check species name in mol_names.dat'
          endif
      end do
      ggchem_index_read = .TRUE.
      endif
      write(file_id,'(i2)') iter
      file_name = 'ppress_' // trim(adjustl(file_id)) // '.dat'
      open(unit=23962, file=file_name, status='replace')
      file_name = 'sumop_' // trim(adjustl(file_id)) // '.dat'
      open(unit=23963, file=file_name, status='replace')
C      open(unit=23975, file='osplots.dat', status='replace')

      if(larciv.eq.1)then
               open(unit=424,file='kappa_plot.dat',status='replace')
               open(unit=430,file='opplot.dat',status='replace')
               ltop = 44
               write(424,426)ltop, t(ltop), nwtot
426            format(' ltop =',i3,' temp = ',f8.2,' nwtot =',i9)
               write(424,*)' molecule number and name etc in OS list:'
            write(424,428)
428         format
     *      ('name nm jvmax wl_mu{op_first/last/max}',
     *       ' wnmax_cm-1 akapmax sumop:')
      end if

      do nm=1, molno
       write(23962,*) molnames_new(nm)
       write(23963,*) molnames_new(nm)
C           if(larciv.eq.1) write(424,425) nm, molnames_new(nm)
C425        format(i4,2x,a4)
          do it=1, ntau
            sumop(nm,it) = 0.
            sumkap(nm,it) =0.
            akapmax(it) = 0.
          end do
          
          opjv(:,:) = 0.
          call opac_wrapper_interp(ptot, t, nm, opjv, ntau)
          if (molnames_new(nm)=='CH4') then
          if (nch4 == 1) then
            open(unit=423, file='fch4.dat', readonly)
            read(423, *) fch4
            close(423)
            opjv(:,:) = fch4 * opjv(:,:)
          end if
          end if
          jvfirst = 0
          jvlast = 0
          locfirst = 0
          lop10 = 0
          do  jv=1,nwtot
            jvn = min(nwtot-1,jv)
            jvpr = max(1,jv-1)
            dw = (wnos(jvn)-wnos(jvpr))/2.

C ....printout for understanding and plotting the opacities
            if(larciv.eq.1)then
              if(opjv(jv,ltop).gt.0.0 .and. locfirst.eq.0) then   ! first opacity gt 0
                       jvfirst = jv
                       locfirst = 1
              end if
C                       if(locfirst.eq.1 .and. lop10.eq.10) then !average over say 10 frequencies
                        do i=1,11
                        nmi = i
                        if(nm.eq.nmid(i)) go to 3401
                        go to 3402
                        end do
3401                    continue
                       if(lop10.eq.50) then !average over say 10 frequencies
                       write(430,*) 'nmi,iop wnav(iop),opav(nmi,iop) = '
                       iop = iop+1
                       opav(nmi,iop) = 0.
                       wnav(iop) = 0.
                       lop10 = 0
                       do jop = 0,49
                       opav(nmi,iop) = 
     &                         opav(nmi,iop) + opjv(jvn-jop,ltop)
                       wnav(iop) = wnav(iop) + wnos(jvn-jop)
                       end do
                       opav(nmi,iop) = opav(nmi,iop)/50.
                       wnav(iop) = wnav(iop)/50.
                       write(430,3410) 
     *                    nmi,iop,wnav(iop),1.e4/wnav(iop),opav(nmi,iop)
3410                   format('1.loop',2i4,f10.2,f10.6,1pe12.3)
                       end if
                       lop10 = lop10 + 1
3402                    continue
              if(opjv(jv,ltop) .eq. 0.0) jvlast = jv    ! when the rest of opacities are 0
            end if
C ....the plot            

            jvs=nwtot
            jvl = nwtot-jv + 1
            if (jvl/nops*nops == jvl) jvs = jvl

            do it=1, ntau
                  TRIX=T(it)*RO(it)
                  if (ggchem_mol(nm)) then
                        if (jv==1) then
                        write(23962,*) ppallmol(it,ggchem_index(nm))
                        end if
                        ppx_os = ppallmol(it,ggchem_index(nm))
                  else
                        ppx_os = ppallat(it,ggchem_index(nm))
                  endif
                  px_os = ppx_os / TRIX*1.20274D-8    

                  akapmol(it) = opjv(jv,it)
                  if(it.eq.ltop .and. akapmol(it) .gt. akapmax(it)) then
                        akapmax(it) = akapmol(it)
                        jvmax = jvn
                        wnmax = wnos(jvn)
                  end if
                  sumop(nm,it) = sumop(nm,it) + akapmol(it)*dw
                  conos(it,jvl) = 
     &            px_os*akapmol(it) + conos(it,jvl)
                  sumkap(nm,it) = 
     &            sumkap(nm,it) + akapmol(it)*px_os*dw
            end do
          end do 

          write(23962, *) 
          do it=1, ntau
          write(23963, *) sumop(nm, it)
          end do
           if (larciv.eq.1) then
            write(424,427)molnames_new(nm),nm,jvmax
     &            ,1.e4/wnos(jvfirst),1.e4/wnos(jvlast),1.e4/wnmax
     &            ,wnmax,akapmax(ltop),sumop(nm,ltop)
           end if
      end do
           if (larciv.eq.1) then
              do i=1,11
              write(430,432)molnames_new(nmid(i)),nmid(i)
     &            ,1.e4/wnos(jvfirst),1.e4/wnos(jvlast)
     &            ,opav(i,1)
              write(430,*) 'i,iop = ',i,iop
              do j=1,50
              if (j.le.50)
     &          write(430,431)i,nmid(i),1.e4/wnav(j),wnav(j),opav(i,j)    !(opav(i,j),in=1,11)
C     *                    nmi,iop,wnav(iop),1.e4/wnav(iop),opav(nmi,iop)
              end do
              end do
           end if  
        close(424)
        close(430)
        write(23963,*) 
      close(23962)
      close(23963)
427           format(2x,a4,1x,i3,i6,3f10.6,f10.2,1p2e12.3)
431           format(2i3,f10.6,f10.2,1p11e12.3)
432           format(2x,a4,1x,i3,3f10.6,f10.2,1p2e12.3)

      FIRST = .FALSE.

      call timex1

      contains
      function upper(strIn) result(strOut)
      ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
      ! Original author: Clive Page

      implicit none

      character(len=*), intent(in) :: strIn
      character(len=len(strIn)):: strOut
      integer :: i,j

      do i = 1, len(strIn)
            j = iachar(strIn(i:i))
            if (j>= iachar("a") .and. j<=iachar("z") ) then
                  strOut(i:i) = achar(iachar(strIn(i:i))-32)
            else
                  strOut(i:i) = strIn(i:i)
            end if
      end do

      end function upper

      END
C
C
      FUNCTION TG01B(II,ND,N,X,F,D,XX)
      implicit real*8 (a-h,o-z)
c      IMPLICIT REAL*8 (A-H),(O-Z)
      DIMENSION X(ND),F(ND),D(ND)
      COMMON /TG01BA/I1,IN,KK
      DATA I1,IN/1,1/
C ROUTINE TO CALCULATE VALUE FXX OF SPLINE IN POINT XX WHEN N KNOTS
C (XI,FI) WITH DERIVATIVE DI ARE GIVEN.
C II<0 => SEARCH THE WHOLE RANGE. II >= 0 => FUNCTION HAS PREVIOUSLY 
C BEEN ENNTERED WITH A SMALLER VALUE OF XX.
C COMMON VALUES I1 AND IN CONTROLS WHAT TO DO IF XX IS OUTSIDE X INTERVAL.
C
C II NEGATIVE, RESET
      IF(II.LT.0) KK=2  
C   
C CHECK IF OUTSIDE  
      IF(XX.LT.X(1)) GOTO 110   
      IF(XX.GT.X(N)) GOTO 120   
      DO 100 K=KK,N 
      IF(XX.LT.X(K)) GOTO 101   
100   CONTINUE  
      KK=N  
      GOTO 102  
101   KK=K  
C   
C CALCULATE FUNCTION
102   DX=X(KK)-X(KK-1)  
      DF=F(KK)-F(KK-1)  
      P=(XX-X(KK-1))/DX 
      Q=1.-P
      TG01B=Q*F(KK-1)+P*F(KK)+P*Q*  
     & (Q*(D(KK-1)*DX-DF)-P*(D(KK)*DX-DF))
1021  FORMAT(I3,1P5E13.5)
      RETURN
C   
C BEFORE X(1)   
110   TG01B=0.  
      IF(I1.LE.0) RETURN
      TG01B=F(1)
      IF(I1.EQ.1) RETURN
      TG01B=TG01B+(XX-X(1))*D(1)
      IF(I1.EQ.2) RETURN
      DX=X(2)-X(1)  
      D2=2.*(3.*(F(2)-F(1))/DX**2-(2.*D(1)+D(2))/DX)
      TG01B=TG01B+.5*(XX-X(1))**2*D2
      IF(I1.EQ.3) RETURN
      D36=(D(1)+D(2)-2.*(F(2)-F(1))/DX)/DX**2   
      TG01B=TG01B+(XX-X(1))*(XX-X(1))**2*D36
      RETURN
C   
C AFTER X(N)
120   TG01B=0.  
      IF(IN.LE.0) RETURN
      TG01B=F(N)
      IF(IN.EQ.1) RETURN
      TG01B=TG01B+(XX-X(N))*D(N)
      IF(IN.EQ.2) RETURN
      DX=X(N)-X(N-1)
      D2=2.*(-3.*(F(N)-F(N-1))/DX**2+(2.*D(N)+D(N-1))/DX)   
      TG01B=TG01B+.5*(XX-X(N))**2*D2
      IF(IN.EQ.3) RETURN
      D36=(D(N)+D(N-1)-2.*(F(N)-F(N-1))/DX)/DX**2   
      TG01B=TG01B+(XX-X(N))*(XX-X(N))**2*D36
      END
C
C
C
      SUBROUTINE DAYTIM(ADATE,ATIME)
      implicit real*8 (a-h,o-z)
C   This routine returns date and time in 9 character format.
C   The routine is highly machine dependent !
C   This version is for VAX.
C
      CHARACTER*9 ADATE,ATIME
      ADATE=' ' 
      ATIME=' ' 
C_ursa      CALL DATE_AND_TIME(ADATE)  
C_ursa      CALL TIME(ATIME)
C      ADATE(5:5)=CHAR(ICHAR(ADATE(5:5))+32)
C      ADATE(6:6)=CHAR(ICHAR(ADATE(6:6))+32)
      RETURN
      END
C
C
      function second()
      implicit real*8 (a-h,o-z)
c  Added by Bjorn S. Nilsson, nbi, on 19-Feb-1991;  for DEC-stations
c      real times(2)
c      second = etime(times)
      second=0
      return
      end
C
C
C      FUNCTION SECOND(TIME0)
C  This routine returns seconds used in floating format.
C  The routine is highly machine dependent !
C  This version is for VAX.
C
C      DATA ICALL/0/ 
C      IF(ICALL.EQ.0) CALL LIB$INIT_TIMER()  
C      ICALL=1
C      CALL LIB$STAT_TIMER(2,LCSEC)  
C      SECOND=0.01*DFLOAT(LCSEC)  
C      RETURN
C      END
C
      SUBROUTINE TIMEX
      implicit real*4 (a-h,o-z)
C   This routine prints total accumulated time and time spent since
C  last call.
C
      CHARACTER*20 FORM/'(A,F05.2,A,F05.2,A)'/
      REAL(KIND=4), SAVE, DIMENSION(2):: time_last=(/0.,0./), time_0
      REAL(KIND=4), SAVE, DIMENSION(2):: time_now=(/0.,0./)
      ENTRY TIMEX0
C      TIME_LAST=SECOND(DUMTIM)  
      call etime(time_last)
      time_0 = time_last
C      TIME_LAST=SECOND()
      RETURN
      ENTRY TIMEX1
C      TIME_NOW=SECOND(DUMTIM)
C      TIME_NOW=SECOND()
      call etime(time_now)
      L1=5
      IF(TIME_NOW(1).GT.99.) L1=6
      IF(TIME_NOW(1).GT.999.) L1=7
      IF(TIME_NOW(1).GT.9999.) L1=8
      IF(TIME_NOW(1).GT.99999.) L1=9
      IF(TIME_NOW(1).GT.999999.) L1=10
      WRITE(FORM(5:6),'(I2.2)') L1
      DELTA_TIME=TIME_NOW(1)-TIME_LAST(1)
      tot_time = time_now(1) - time_0(1)
      L1=5
      IF(DELTA_TIME.GT.99.) L1=6
      IF(DELTA_TIME.GT.999.) L1=7
      IF(DELTA_TIME.GT.9999.) L1=8
      IF(DELTA_TIME.GT.99999.) L1=9
      IF(DELTA_TIME.GT.999999.) L1=10
      WRITE(FORM(13:14),'(I2.2)') L1
      
      TIME_LAST=TIME_NOW
      RETURN
      END
C
C
      SUBROUTINE TIMEF
      implicit real*4 (a-h,o-z)
C   This routine prints total accumulated time and time spent since
C  last call.
C
      REAL(KIND=4), SAVE, DIMENSION(2):: time_last=(/0.,0./), time_0
      REAL(KIND=4), SAVE, DIMENSION(2):: time_now=(/0.,0./)
      CHARACTER*20 FORM/'(A,F05.2,A,F05.2,A)'/
      ENTRY TIME0
C      TIME_LAST=SECOND(DUMTIM)  
      call etime(time_last)
      time_0 = time_last
C      TIME_LAST=SECOND()
      RETURN
      ENTRY TIME1
C      TIME_NOW=SECOND(DUMTIM)
C      TIME_NOW=SECOND()
      call etime(time_now)
      L1=5
      IF(TIME_NOW(1).GT.99.) L1=6
      IF(TIME_NOW(1).GT.999.) L1=7
      IF(TIME_NOW(1).GT.9999.) L1=8
      IF(TIME_NOW(1).GT.99999.) L1=9
      IF(TIME_NOW(1).GT.999999.) L1=10
      WRITE(FORM(5:6),'(I2.2)') L1
      DELTA_TIME=TIME_NOW(1)-TIME_LAST(1)
      tot_time = time_now(1) - time_0(1)
      L1=5
      IF(DELTA_TIME.GT.99.) L1=6
      IF(DELTA_TIME.GT.999.) L1=7
      IF(DELTA_TIME.GT.9999.) L1=8
      IF(DELTA_TIME.GT.99999.) L1=9
      IF(DELTA_TIME.GT.999999.) L1=10
      WRITE(FORM(13:14),'(I2.2)') L1
      PRINT FORM, ' Total time spent: ',tot_time,' sec. Added time: ',
     1 DELTA_TIME,' sec.'
      TIME_LAST=TIME_NOW
      RETURN
      END
      
!-----------------------------------------------------------------------
! GETTIME: Prints total accumulated time and time spent since last call
! Juncher 2015
! UGJ, 2019: etime is a real*4 routine and variables have to be defined
! as such implicit as well as explicit in order to work.
!-----------------------------------------------------------------------
      subroutine gettime(io)
      
      implicit real*4 (a-h,o-z)
      include 'parameter.inc'
      integer :: tot_hour, tot_min, tot_sec
      real*4, save, dimension(2):: time_0, time_last, time_now
      
      if(io .eq. 0) then
        call etime(time_last)
        time_0 = time_last
      else
        call etime(time_now)

        tot_hour = int((time_now(1)-time_0(1))/3600.)
        tot_min  = int((time_now(1)-time_0(1)-tot_hour*3600.)/60.)
        tot_sec  = int(time_now(1)-time_0(1)-tot_hour*3600.-tot_min*60.)
        
        time_last = time_now
        
        if(io .eq. 1) then
          write(*,'(a19,i2,2(a1,i2.2))')' Total time spent: ',
     *      tot_hour, ':', tot_min, ':', tot_sec
        end if
      end if

      return
      end


      SUBROUTINE SCALEMOD
      implicit real*8 (a-h,o-z)
C
C This routine scales all parameters in one model to a model with another
C number of depth points. The (new) tau-scale from the input file is in
C common TAUC, and the variables to be scaled are those in common STATEC
C
      include 'parameter.inc'
C
      DIMENSION W(3,NDP),SX(NDP),SXDER(NDP)
C
      COMMON /STATEC/PR(NDP),PT(NDP),P(NDP),G(NDP),Z(NDP),D(NDP),
     &V(NDP),FC(NDP),PE(NDP),T(NDP),TAULN(NDP),ROSTAT(NDP),
     &NTAU,ITER
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /TG01BA/I1,IN,KK
      COMMON /CMOLRAT/ FOLD(NDP,8),MOLOLD,KL
C
C      DATA I1,IN/2,2/ !0~F==0. outside interpolation interval
C                      1~F==value in first (last) point -||-
C                      2~F from linear extrapolation -||-
C ---------------- ***************** ----------------------------  
C CALCULATE log-log SPLINE DERIVATIVES of variables in STATEC
C
C SxDER IS THE DERIVATIVE IN ALL NTAU NODES WHEN
C TAULN AND SxK ARE GIVEN 
C Calculate interpolated values at depths of new tauscale. Save
C in STATEC and TAUC

      
      if (ntau.gt.ndp) STOP  ' increase dimension ndp '
      I1 = 2
      IN = 2
          !print*, "ntau ", ntau
          DO 1 K=1,NTAU
1         Sx(K) = LOG(P(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 11 K=1,JTAU
            !print*, k
            !print*, "TAU(K)", tau(k)
          X = LOG(TAU(K))
11        P(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C          DO 2 K=1,NTAU
C2         Sx(K) = LOG(G(K))
C          CALL TB04A(NTAU,tauln,SX,SXDER,W)
C          DO 12 K=1,JTAU
C          X = LOG(TAU(K))
C12        G(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C          DO 3 K=1,NTAU
C3         Sx(K) = LOG(Z(K))
C          CALL TB04A(NTAU,tauln,SX,SXDER,W)
C          DO 13 K=1,JTAU
C          X = LOG(TAU(K))
C13        Z(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
          DO 4 K=1,NTAU
4         Sx(K) = LOG(PE(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 14 K=1,JTAU
          X = LOG(TAU(K))
14        PE(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
          DO 5 K=1,NTAU
5         Sx(K) = LOG(T(K))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 15 K=1,JTAU
          X = LOG(TAU(K))
15        T(K)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
C
C
      IF (MOLOLD.EQ.1) THEN
         DO 26 M=1,8
          DO 6 K=1,NTAU
6         Sx(K) = LOG(FOLD(K,M))
          CALL TB04A(NTAU,tauln,SX,SXDER,W)
          DO 16 K=1,JTAU
          X = LOG(TAU(K))
16        FOLD(K,M)=EXP(TG01B(-1,NDP,NTAU,TAULN,SX,SXDER,X))
26       CONTINUE
      END IF
C
      NTAU=JTAU
      DO 20 K=1,NTAU
20    G(K)=G(1)
C20    TAULN(K)=LOG(TAU(K))
C
      
      
301   FORMAT(I3,1P5E11.3,0PF8.0)
C
      RETURN
      END
c
c
      SUBROUTINE TB04A(N,X,F,D,W)
      implicit real*8 (a-h,o-z)
C
      DIMENSION X(N),F(N),D(N),W(3,N)   
C   
C THIS VERSION OF TB04A IS COMPATIBLE WITH RKU*HARWELL.TB04A, BUT   
C FOR UNKNOWN REASONS 30 % FASTER.  770815
C INPUT ARE X(I), FUNCTION VALUE IN N KNOTS. THEN THE DERIVATIVE IS 
C CALCULATED IN THE KNOTS. W=0 AT SUCCESFULL RETURN, OTHERWISE W=1.
C X(I) SHOULD BE IN STRICTH INCREASING ORDER, X1<X2<...<XN.
C         
C FIRST POINT   
      CXB=1./(X(2)-X(1))
      CXC=1./(X(3)-X(2))
      DFB=F(2)-F(1) 
      DFC=F(3)-F(2) 
      W(1,1)=CXB*CXB
      W(3,1)=-CXC*CXC   
      W(2,1)=W(1,1)+W(3,1)  
      D(1)=2.*(DFB*CXB*CXB*CXB-DFC*CXC*CXC*CXC) 
C   
C INTERIOR POINTS   
      N1=N-1
      DO 100 K=2,N1 
      CXA=CXB   
      CXB=1./(X(K+1)-X(K))  
      DFA=DFB   
      DFB=F(K+1)-F(K)   
      W(1,K)=CXA
      W(3,K)=CXB
      W(2,K)=2.*(CXA+CXB)   
      D(K)=3.*(DFB*CXB*CXB+DFA*CXA*CXA) 
100   CONTINUE  
C   
C LAST POINT
      W(1,N)=CXA*CXA
      W(3,N)=-CXB*CXB   
      W(2,N)=W(1,N)+W(3,N)  
      D(N)=2.*(DFA*CXA*CXA*CXA-DFB*CXB*CXB*CXB) 
C   
C ELIMINATE AT FIRST POINT  
      C=-W(3,1)/W(3,2)  
      W(1,1)=W(1,1)+C*W(1,2)
      W(2,1)=W(2,1)+C*W(2,2)
      D(1)=D(1)+C*D(2)  
      W(3,1)=W(2,1) 
      W(2,1)=W(1,1) 
C   
C ELIMINATE AT LAST POINT   
      C=-W(1,N)/W(1,N-1)
      W(2,N)=W(2,N)+C*W(2,N-1)  
      W(3,N)=W(3,N)+C*W(3,N-1)  
      D(N)=D(N)+C*D(N-1)
      W(1,N)=W(2,N) 
      W(2,N)=W(3,N) 
C   
C ELIMINATE SUBDIAGONAL 
      DO 110 K=2,N  
      C=-W(1,K)/W(2,K-1)
      W(2,K)=W(2,K)+C*W(3,K-1)  
      D(K)=D(K)+C*D(K-1)
110   CONTINUE  
C   
C BACKSUBSTITUTE
      D(N)=D(N)/W(2,N)  
      DO 120 KK=2,N 
      K=(N+1)-KK
      D(K)=(D(K)-W(3,K)*D(K+1))/W(2,K)  
120   CONTINUE  
C   
      RETURN
      END



      FUNCTION TG01BT(II,ND,N,X,F,D,XX) 
      implicit real*8 (a-h,o-z)
      DIMENSION X(ND),F(ND),D(ND)
      COMMON /TG01BA/I1,IN,KK
      I1 = 1     ! we fix extrapolated values to be == end value
      IN = 1
C ROUTINE TO CALCULATE VALUE FXX OF SPLINE IN POINT XX WHEN N KNOTS
C (XI,FI) WITH DERIVATIVE DI ARE GIVEN.
C II<0 => SEARCH THE WHOLE RANGE. II >= 0 => FUNCTION HAS PREVIOUSLY 
C BEEN ENTERED WITH A SMALLER VALUE OF XX.
C COMMON VALUES I1 AND IN CONTROLS WHAT TO DO IF XX IS OUTSIDE X INTERVAL.
C
C II NEGATIVE, RESET
      IF(II.LT.0) KK=2  
C   
C CHECK IF OUTSIDE  
      IF(XX.LT.X(1)) GOTO 110   
      IF(XX.GT.X(N)) GOTO 120   
      DO 100 K=KK,N 
      IF(XX.LT.X(K)) GOTO 101   
100   CONTINUE  
      KK=N  
      GOTO 102  
101   KK=K  
C   
C CALCULATE FUNCTION
102   DX=X(KK)-X(KK-1)  
      DF=F(KK)-F(KK-1)  
      P=(XX-X(KK-1))/DX 
      Q=1.-P
      TG01BT=Q*F(KK-1)+P*F(KK)+P*Q*  
     & (Q*(D(KK-1)*DX-DF)-P*(D(KK)*DX-DF))  
      RETURN
C   
C BEFORE X(1)   
110   TG01BT=0.  
      IF(I1.LE.0) RETURN
      TG01BT=F(1)
      IF(I1.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(1))*D(1)
      IF(I1.EQ.2) RETURN
      DX=X(2)-X(1)  
      D2=2.*(3.*(F(2)-F(1))/DX**2-(2.*D(1)+D(2))/DX)
      TG01BT=TG01BT+.5*(XX-X(1))**2*D2
      IF(I1.EQ.3) RETURN
      D36=(D(1)+D(2)-2.*(F(2)-F(1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(1))*(XX-X(1))**2*D36
      RETURN
C   
C AFTER X(N)
120   TG01BT=0.  
      IF(IN.LE.0) RETURN
      TG01BT=F(N)
      IF(IN.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(N))*D(N)
      IF(IN.EQ.2) RETURN
      DX=X(N)-X(N-1)
      D2=2.*(-3.*(F(N)-F(N-1))/DX**2+(2.*D(N)+D(N-1))/DX)   
      TG01BT=TG01BT+.5*(XX-X(N))**2*D2
      IF(IN.EQ.3) RETURN
      D36=(D(N)+D(N-1)-2.*(F(N)-F(N-1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(N))*(XX-X(N))**2*D36
      END
C
C
C
      SUBROUTINE TB04AT(ND,N,X,F,D,W)
      implicit real*8 (a-h,o-z)
C
      DIMENSION X(ND),F(ND),D(ND),W(3,ND)   
C   
C THIS VERSION OF TB04A IS identical to TB04A, except that it can be
C called with arrays which might be dimensioned larger than the part
C of it which is used for the interpolation.
C INPUT ARE X(I), FUNCTION VALUE IN N KNOTS. THEN THE DERIVATIVE IS 
C CALCULATED IN THE KNOTS. W=0 AT SUCCESFULL RETURN, OTHERWISE W=1.
C X(I) SHOULD BE IN STRICTH INCREASING ORDER, X1<X2<...<XN.
C In the call from SUBROUTINE PROFILE, X(1), X(2),...,X(NTAU) are the 
C temperature values at the NTAU optical depth values, and the derivative
C dPg/dT (=D(N))is calculated in the NTAU points, for later use to 
C calculate Pg in the points where the temperature of the
C absorption coefficient is known.  
C         
C FIRST POINT   
      CXB=1./(X(2)-X(1))
      CXC=1./(X(3)-X(2))
      DFB=F(2)-F(1) 
      DFC=F(3)-F(2) 
      W(1,1)=CXB*CXB
      W(3,1)=-CXC*CXC   
      W(2,1)=W(1,1)+W(3,1)  
      D(1)=2.*(DFB*CXB*CXB*CXB-DFC*CXC*CXC*CXC) 
C   
C INTERIOR POINTS   
      N1=N-1
      DO 100 K=2,N1 
      CXA=CXB   
      CXB=1./(X(K+1)-X(K))  
      DFA=DFB   
      DFB=F(K+1)-F(K)   
      W(1,K)=CXA
      W(3,K)=CXB
      W(2,K)=2.*(CXA+CXB)   
      D(K)=3.*(DFB*CXB*CXB+DFA*CXA*CXA) 
100   CONTINUE  
C   
C LAST POINT
      W(1,N)=CXA*CXA
      W(3,N)=-CXB*CXB   
      W(2,N)=W(1,N)+W(3,N)  
      D(N)=2.*(DFA*CXA*CXA*CXA-DFB*CXB*CXB*CXB) 
C   
C ELIMINATE AT FIRST POINT  
      C=-W(3,1)/W(3,2)  
      W(1,1)=W(1,1)+C*W(1,2)
      W(2,1)=W(2,1)+C*W(2,2)
      D(1)=D(1)+C*D(2)  
      W(3,1)=W(2,1) 
      W(2,1)=W(1,1) 
C   
C ELIMINATE AT LAST POINT   
      C=-W(1,N)/W(1,N-1)
      W(2,N)=W(2,N)+C*W(2,N-1)  
      W(3,N)=W(3,N)+C*W(3,N-1)  
      D(N)=D(N)+C*D(N-1)
      W(1,N)=W(2,N) 
      W(2,N)=W(3,N) 
C   
C ELIMINATE SUBDIAGONAL 
      DO 110 K=2,N  
      C=-W(1,K)/W(2,K-1)
      W(2,K)=W(2,K)+C*W(3,K-1)  
      D(K)=D(K)+C*D(K-1)
110   CONTINUE  
C   
C BACKSUBSTITUTE
      D(N)=D(N)/W(2,N)  
      DO 120 KK=2,N 
      K=(N+1)-KK
      D(K)=(D(K)-W(3,K)*D(K+1))/W(2,K)  
120   CONTINUE  
C   
      RETURN
      END
C
   
      subroutine atoms_head
* reads integrated metal opacity file in ascii format
      implicit real*8 (a-h,o-z)
C     implicit none
*
      integer np6,ntemp,maxpnt,nel
      parameter (np6 = 9)
      parameter (ntemp = 17)
      parameter (maxpnt = 153910)
      parameter (nel = 92)
*
      integer ispec(184),itemp,ip6,nwave,i,n,nwav
      real*8 waven
C     real p6(np6),temp(ntemp),abund(nel)
      real abund(nel)
      real opl(np6,ntemp)
      real xite
      character*60 filename
      character*80 string
*
c     data filename /'/ste1/uffegj/atoms/metals_sun_ascii.x03'/
      common /catoms_head/ p6(9),temp(17)
c     print *,'name of the input file ?'
c     read(*,'(a)') filename(1)
c     open(11,file=filename(1),form='unformatted',status='old')
c     print *,'name of the output file ?'
c     read(*,'(a)') filename(2)
c     open(41,file=filename,form='formatted',recl=1246,
c    &     status='unknown',readonly)
      open(13,file='metals.output',status='unknown')
*
      read(41,1080) string
      write(13,1080) string
* 1) info record
 1080 format(a80)
      
      read(41,1080) string
      write(13,1080) string
* 2) info record
      
      read(41,1000) ispec
      write(13,1000) ispec
* 3) ispec:  integers identifying all species included:
*            ex: 9201 = U II, 600 = C I    (max number=2*92)
 1000 format(184i5)
      
      read(41,1010) xite
      write(13,1010) xite
* 4) xite: [km/s]  microturbulence velocity (for line broadening)
 1010 format(f5.2)
      
      read(41,1020) ip6     !1020 or *
      if(ip6.ne.np6) then
       
        stop 'Data file not as expected'
      endif
      write(13,1020) np6
* 5) np6: number of P6 values (must always=9)
 1020 format(i3)
      
      read(41,1030) (p6(i),i=1,np6)
      write(13,1030) (p6(i),i=1,np6)
* 6) P6: [dyn] Damping pressures for "van der Waals" line broadening
*        to be computed as P(HI) + 0.42*P(HeI) + 0.85*P(H2)
 1030 format(1p,9e10.2,0p)
      
      read(41,1020) itemp
      if(itemp.ne.ntemp) then
        
        stop 'Data file not as expected'
      endif
      write(13,1020) ntemp
* 7) ntemp: number of T values (must always=17)
     
      read(41,1070) (temp(i),i=1,ntemp)
      write(13,1070) (temp(i),i=1,ntemp)
* 8) temp: [K] temperature values for excitation, ionization, chemical
*       equilibrium, line broadening (LTE)
 1070 format(17f8.0)
      
      read(41,1040) nwave
      write(13,1040) nwave
* 9) nwave will probably be 153910, the expected number of wavelength points
 1040 format(i10)
      
      read(41,1080)
      write(13,*) 'dummy record'
* 10) dummy record: unly used in single species files
      
      read(41,1080)
      write(13,*) 'dummy record'
* 11) dummy record: unly used in single species files
      
      read(41,1050,err=3879) abund
      write(13,1050) abund
* 12) abund: logarithmic number abundances used in the compilation of the file
*            on a scale where H = 12.00, 92 values, -99. means not used/known
 1050 format(92f7.2)
      
      goto 3880
 3879 continue
      
 3880 continue
*
* 12 header lines ready, write 153910 opacity data recods:
*
      nwav=0
     
      write(13,*) ' now starting big read; maxpnt,np6,ntemp=',
     &     maxpnt,np6,ntemp
C     do n=1,maxpnt+1
C     read(41,1060,end=98,err=99) 
C    &   waven,((opl(ip6,itemp),ip6=1,np6),itemp=1,ntemp)
* waven: [cm-1] vacuum wave number
* opl: log(10) opacity, opacity in cm2/g stellar matter
*      opl = -30.0 or -40.0 means no significant metal line opacity found
* acurracy: Line data: b-b transition line data of neutral and singly-ionized
*           atoms were all adopted from VALD, January 1998. (cf A&AS 112, 525)
*           For 36 rare species no line data is there or the solar abundance=0
*           For the species where continuous opacities are considered,
*           the line opacity data [cm2/g] is complete to better than 1% of
*           the continuous opacity (at the relevant T and wavel). Based on the
*           experience from these species similar completness levels were
*           adopted for all other species.

      return
      end


      SUBROUTINE atomw
C
C Mean atomic weights are stored into the vector watom, AB2001
C Species 1:92 are included, as in irwin.dat
C Deuterium is added as number 93, because it is treated as a seperate
C atom in gfits.data of molecules. Actually W(H)=1.007825, while Earth
C mixture of H and D has W(H+D)=1.0079 as given for watom(1).
C

      IMPLICIT real*8 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)


      DO i=1,200
         watom(i) = 0.0
      ENDDO

      watom(1)  =   1.0079 !H  Hydrogen
      watom(2)  =   4.0026 !He Helium
      watom(3)  =   6.941  !Li Lithium
      watom(4)  =   9.0121 !Be Beryllium
      watom(5)  =  10.81   !B  Boron
      watom(6)  =  12.011  !C  Carbon
      watom(7)  =  14.0067 !N  Nitrogen
      watom(8)  =  15.9994 !O  Oxygen
      watom(9)  =  18.9984 !F  Fluorine
      watom(10) =  20.179  !Ne Neon
      watom(11) =  22.9897 !Na Sodium
      watom(12) =  24.305  !Mg Magnesium
      watom(13) =  26.9814 !Al Aluminum
      watom(14) =  28.0855 !Si Silicon
      watom(15) =  30.9737 !P  Phosphorus
      watom(16) =  32.06   !S  Sulfur
      watom(17) =  35.453  !Cl Chlorine
      watom(18) =  39.948  !Ar Argon
      watom(19) =  39.0983 !K  Potassium
      watom(20) =  40.08   !Ca Calcium
      watom(21) =  44.9559 !Sc Scandium
      watom(22) =  47.88   !Ti Titanium
      watom(23) =  50.9415 !V  Vanadium
      watom(24) =  51.996  !Cr Chromium
      watom(25) =  54.9380 !Mn Manganese
      watom(26) =  55.847  !Fe Iron
      watom(27) =  58.9332 !Co Cobalt
      watom(28) =  58.96   !Ni Nickel
      watom(29) =  63.546  !Cu Copper
      watom(30) =  65.38   !Zn Zinc
      watom(31) =  69.72   !Ga Gallium
      watom(32) =  72.59   !Ge Germanium
      watom(33) =  74.9216 !As Arsenic
      watom(34) =  78.96   !Se Selenium
      watom(35) =  79.904  !Br Bromine
      watom(36) =  83.80   !Kr Krypton
      watom(37) =  85.4678 !Rb Rubidium
      watom(38) =  87.62   !Sr Strontium
      watom(39) =  88.9059 !Y  Yttrium
      watom(40) =  91.22   !Zr Zirconium
      watom(41) =  92.9064 !Nb Niobium
      watom(42) =  95.94   !Mo Molybdenum
      watom(43) =  97.907  !Tc Technetium
      watom(44) = 101.07   !Ru Ruthenium
      watom(45) = 102.9055 !Rh Rhodium
      watom(46) = 106.42   !Pd Palladium
      watom(47) = 107.868  !Ag Silver
      watom(48) = 112.41   !Cd Cadmium
      watom(49) = 114.82   !In Indium
      watom(50) = 118.69   !Sn Tin
      watom(51) = 121.75   !Sb Antimony
      watom(52) = 127.60   !Te Tellurium
      watom(53) = 126.9045 !I  Iodine
      watom(54) = 131.29   !Xe Xenon
      watom(55) = 132.9054 !Cs Cesium
      watom(56) = 137.33   !Ba Barium
      watom(57) = 138.9055 !La Lanthanum
      watom(58) = 140.12   !Ce Cerium
      watom(59) = 140.9077 !Pr Praseodymium
      watom(60) = 144.24   !Nd Neodymium
      watom(61) = 144.913  !Pm Promethium
      watom(62) = 150.36   !Sm Samarium
      watom(63) = 151.96   !Eu Europium
      watom(64) = 157.25   !Gd Gadolinium
      watom(65) = 158.9254 !Tb Terbium
      watom(66) = 162.50   !Dy Dysprosium
      watom(67) = 164.9304 !Ho Holmium
      watom(68) = 167.26   !Er Erbium
      watom(69) = 168.9342 !Tm Thulium
      watom(70) = 173.04   !Yb Ytterbium
      watom(71) = 174.967  !Lu Lutetium
      watom(72) = 178.49   !Hf Hafnium
      watom(73) = 180.9479 !Ta Tantalum
      watom(74) = 183.85   !W  Tungsten
      watom(75) = 186.207  !Re Rhenium
      watom(76) = 190.2    !Os Osmium
      watom(77) = 192.22   !Ir Iridium
      watom(78) = 195.08   !Pt Platinum
      watom(79) = 196.9665 !Au Gold
      watom(80) = 200.59   !Hg Mercury
      watom(81) = 204.383  !Tl Thallium
      watom(82) = 207.2    !Pb Lead
      watom(83) = 208.9804 !Bi Bismuth
      watom(84) = 208.982  !Po Polonium
      watom(85) = 209.987  !At Astatine
      watom(86) = 222.018  !Rn Radon
      watom(87) = 223.020  !Fr Francium
      watom(88) = 226.0254 !Ra Radium
      watom(89) = 227.0278 !Ac Actinium
      watom(90) = 232.0381 !Th Thorium
      watom(91) = 231.0359 !Pa Protactinium
      watom(92) = 238.051  !U  Uranium
      watom(93) = 2.01410  !D  Deuterium

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE atomnam
C
C 1 or 2 character names of atoms.
C Species 1:92 are included, as in irwin.dat
C


      implicit REAL*8 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)

      DO i=1,200
         atomname(i) = '  '
      ENDDO

      atomname(1)  = 'H ' ! Hydrogen
      atomname(2)  = 'He' ! Helium
      atomname(3)  = 'Li' ! Lithium
      atomname(4)  = 'Be' ! Beryllium
      atomname(5)  = 'B ' ! Boron
      atomname(6)  = 'C ' ! Carbon
      atomname(7)  = 'N ' ! Nitrogen
      atomname(8)  = 'O ' ! Oxygen
      atomname(9)  = 'F ' ! Fluorine
      atomname(10) = 'Ne' ! Neon
      atomname(11) = 'Na' ! Sodium
      atomname(12) = 'Mg' ! Magnesium
      atomname(13) = 'Al' ! Aluminum
      atomname(14) = 'Si' ! Silicon
      atomname(15) = 'P ' ! Phosphorus
      atomname(16) = 'S ' ! Sulfur
      atomname(17) = 'Cl' ! Chlorine
      atomname(18) = 'Ar' ! Argon
      atomname(19) = 'K ' ! Potassium
      atomname(20) = 'Ca' ! Calcium
      atomname(21) = 'Sc' ! Scandium
      atomname(22) = 'Ti' ! Titanium
      atomname(23) = 'V ' ! Vanadium
      atomname(24) = 'Cr' ! Chromium
      atomname(25) = 'Mn' ! Manganese
      atomname(26) = 'Fe' ! Iron
      atomname(27) = 'Co' ! Cobalt
      atomname(28) = 'Ni' ! Nickel
      atomname(29) = 'Cu' ! Copper
      atomname(30) = 'Zn' ! Zinc
      atomname(31) = 'Ga' ! Gallium
      atomname(32) = 'Ge' ! Germanium
      atomname(33) = 'As' ! Arsenic
      atomname(34) = 'Se' ! Selenium
      atomname(35) = 'Br' ! Bromine
      atomname(36) = 'Kr' ! Krypton
      atomname(37) = 'Rb' ! Rubidium
      atomname(38) = 'Sr' ! Strontium
      atomname(39) = 'Y ' ! Yttrium
      atomname(40) = 'Zr' ! Zirconium
      atomname(41) = 'Nb' ! Niobium
      atomname(42) = 'Mo' ! Molybdenum
      atomname(43) = 'Tc' ! Technetium
      atomname(44) = 'Ru' ! Ruthenium
      atomname(45) = 'Rh' ! Rhodium
      atomname(46) = 'Pd' ! Palladium
      atomname(47) = 'Ag' ! Silver
      atomname(48) = 'Cd' ! Cadmium
      atomname(49) = 'In' ! Indium
      atomname(50) = 'Sn' ! Tin
      atomname(51) = 'Sb' ! Antimony
      atomname(52) = 'Te' ! Tellurium
      atomname(53) = 'I ' ! Iodine
      atomname(54) = 'Xe' ! Xenon
      atomname(55) = 'Cs' ! Cesium
      atomname(56) = 'Ba' ! Barium
      atomname(57) = 'La' ! Lanthanum
      atomname(58) = 'Ce' ! Cerium
      atomname(59) = 'Pr' ! Praseodymium
      atomname(60) = 'Nd' ! Neodymium
      atomname(61) = 'Pm' ! Promethium
      atomname(62) = 'Sm' ! Samarium
      atomname(63) = 'Eu' ! Europium
      atomname(64) = 'Gd' ! Gadolinium
      atomname(65) = 'Tb' ! Terbium
      atomname(66) = 'Dy' ! Dysprosium
      atomname(67) = 'Ho' ! Holmium
      atomname(68) = 'Er' ! Erbium
      atomname(69) = 'Tm' ! Thulium
      atomname(70) = 'Yb' ! Ytterbium
      atomname(71) = 'Lu' ! Lutetium
      atomname(72) = 'Hf' ! Hafnium
      atomname(73) = 'Ta' ! Tantalum
      atomname(74) = 'W ' ! Tungsten
      atomname(75) = 'Re' ! Rhenium
      atomname(76) = 'Os' ! Osmium
      atomname(77) = 'Ir' ! Iridium
      atomname(78) = 'Pt' ! Platinum
      atomname(79) = 'Au' ! Gold
      atomname(80) = 'Hg' ! Mercury
      atomname(81) = 'Tl' ! Thallium
      atomname(82) = 'Pb' ! Lead
      atomname(83) = 'Bi' ! Bismuth
      atomname(84) = 'Po' ! Polonium
      atomname(85) = 'At' ! Astatine
      atomname(86) = 'Rn' ! Radon
      atomname(87) = 'Fr' ! Francium
      atomname(88) = 'Ra' ! Radium
      atomname(89) = 'Ac' ! Actinium
      atomname(90) = 'Th' ! Thorium
      atomname(91) = 'Pa' ! Protactinium
      atomname(92) = 'U ' ! Uranium
      atomname(93) = 'D ' ! Deuterium

      RETURN
      END


      module mie_precision
      use, intrinsic :: iso_fortran_env ! Requires fortran 2008
      implicit none

      !!!
      ! Different sets of single, double and quad precision availible for DIHRT
      ! Try different sets should one fail to compile/give errors for any reason
      ! This module should be compiled first in the DIHRT chain and used in every
      ! module/subroutine when required
      !!!

      private
      public :: sp, dp, qp


      ! Fortran 2008 intrinsic precisions - reccomonded if possible
      integer, parameter :: sp = REAL32
      integer, parameter :: dp = REAL64
      integer, parameter :: qp = REAL128

      end module mie_precision

      module mie_data
        
      use mie_precision,ONLY: sp,dp,qp
      !use drift_data,ONLY: n_dust => NDUST
      implicit none
        
        
        
      ! Constants
      real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp) 
      real(kind=dp), parameter :: pi2 = 8.0_dp * atan(1.0_dp)
      real(kind=dp), parameter :: onethird = 1.0_dp/3.0_dp

      ! Parameters
      !integer, parameter :: n_wl = 31     ! Number of wavelengths
      logical, parameter :: Brug = .True. ! Use Bruggeman method? .False. = LLL method
      integer, parameter :: a_type = 0    ! 0 = mean grain size, 1 = effecitve grain size


      end module mie_data
!      ********************************************************************
!New Mie routine for cloud opacity calculation (June 2023)

! ******************************************************************
! MIEX: MIE SCATTERING CODE FOR LARGE GRAINS
! BCE: Mie code by Sebastian Wolf - see contact below!
!  _____________________________________________________
!  Contact information:   swolf@mpia.de (Sebastian Wolf)
! ==================================================================

      module datatype
      implicit none
      integer, parameter, public ::  r1=selected_real_kind(1)  ! real*4
      integer, parameter, public ::  r2=selected_real_kind(9)  ! real*8 (double precision)
      end module datatype


! ====================================================================================================
! Collection of subroutines
! ====================================================================================================
      module mie_routines
      private :: aa2
      public  :: shexqnn2
      contains
  ! ==================================================================================================
  ! Subroutine for calculations of the ratio of derivative to the function for Bessel functions
  ! of half order with complex argument: J'(n)/J(n). The calculations are given by the recursive
  ! expression ``from top to bottom'' beginning from n=num.
  ! *  a=1/x (a=2*pi*a(particle radius)/lambda - size parameter).
  ! *  ri - complex refractive index.
  ! *  ru-array of results.
  ! - this routine is based on the routine 'aa' published by
  !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
  !                           Astrophysics and Space Physics Review 12,  1 (2002)
  ! ==================================================================================================
      pure subroutine aa2( a, ri, num, ru )
      use datatype

      implicit none

      ! variables for data exchange.....................................................................
      real(kind=r2), intent(in)                   :: a
      complex(kind=r2), intent(in)                :: ri
      integer, intent(in)                         :: num
      complex(kind=r2), dimension(:), intent(out) :: ru

      ! local variables.................................................................................
      integer :: i, i1, j, num1
      complex(kind=r2) :: s, s1
      !-------------------------------------------------------------------------------------------------
      ! initialisierung: not necessary (+ slowes the code down remarkably)
      ! ru(:) = (0.0, 0.0)

      s       = a / ri
      ru(num) = real(num+1,kind=r2) * s
      num1    = num - 1
      do j=1, num1
            i     = num - j
            i1    = i + 1
            s1    = i1 * s
            ru(i) = s1 - 1.0_r2 / (ru(i1) + s1)
      end do
      end subroutine aa2


      !===================================================================================================
      ! shexqnn2
      ! --------
      ! - for a given size parameter 'x' and (complex) refractive index 'ri' the following quantities
      !   are determined:
      !   * Qext     - extinction effiency
      !   * Qsca     - scattering effiency
      !   * Qabs     - absorption effiency
      !   * Qbk      - backscattering effiency
      !   * Qpr      - radiation pressure effiency
      !   * albedo   - Albedo
      !   * g        - g scattering assymetry factor
      !   * SA1, SA2 - scattering amplitude function
      ! - further input parameters
      !   * doSA = .true.  ->  calculation of the scattering amplitudes
      !   * nang ... half number of scattering angles theta in the intervall 0...PI/2
      !              (equidistantly distributed)
      ! - this routine is based on the routine 'shexqnn' published by
      !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
      !                           Astrophysics and Space Physics Review 12,  1 (2002)
      !===================================================================================================
      subroutine shexqnn2( ri, x, Qext, Qsca, Qabs, Qbk, 
     * Qpr, albedo, g, ier, SA1, SA2, doSA, nang )
      use datatype

      implicit none

      ! variables for data exchange.....................................................................
      complex(kind=r2), intent(in)  :: ri
      real(kind=r2), intent(in)     :: x
      real(kind=r2), intent(out)    :: Qext, Qsca, Qabs, 
     > Qbk, Qpr, albedo, g
      integer, intent(out)          :: ier
      complex(kind=r2), dimension(:), intent(out) :: SA1, SA2
      logical, intent(in)           :: doSA
      integer, intent(in)           :: nang

      ! local variables.................................................................................
      integer       :: iterm, nterms, num, iu0, 
     * iu1, iu2, iang2, iang
      real(kind=r2) :: r_iterm, factor, eps, pi, ax, 
     * besJ0, besJ1, besJ2, besY0, besY1, besY2, b, an, 
     * y, ass, w1, qq, fac, an2, P, T, Si, Co, z, xmin

      complex(kind=r2) :: ra0, rb0, ra1, rb1, r, 
     * ss, s1, s2, s3, s, rr

      real(kind=r2),dimension(0:1)   :: fact
      real(kind=r2),dimension(:),allocatable,save :: mu, fpi, 
     *                                     fpi0, fpi1, ftau
      complex(kind=r2),dimension(:),allocatable,save :: ru
      !$omp threadprivate(ru,mu,fpi,fpi0,fpi1,ftau)
      !-------------------------------------------------------------------------------------------------
      ! Maximum number of terms to be considered
      nterms = 100000000
      ! this works for x up to 1.d9, but needs a hell more of memory!!
      !nterms= 550000000

      ! Accuracy to be achieved
      eps    = 1.0e-15_r2

      ! Minimum size parameter
      xmin   = 1.0e-6_r2

      !-------------------------------------------------------------------------------------------------
      ! initialization
      if (.not.allocated(ru)) then
            allocate( ru(1:nterms), mu(1:nang), fpi(1:nang), 
     >              fpi0(1:nang), fpi1(1:nang), ftau(1:nang))
      endif
      ier     = 0
      Qext    = 0.0_r2
      Qsca    = 0.0_r2
      Qabs    = 0.0_r2
      Qbk     = 0.0_r2
      Qpr     = 0.0_r2
      albedo  = 0.0_r2
      g       = 0.0_r2
      fact(0) = 1.0_r2
      fact(1) = 1.0e+250_r2
      factor  = 1.0e+250_r2

      ! null argument
      if (x <= xmin) then
            ier = 1
            print *, "<!> Error in subroutine shexqnn2:"
            print *, "    - Mie scattering limit exceeded:"
            print *, "      current size parameter: ", x
      else
            pi = 4.0_r2 * atan(1.0_r2) ! PI = 3.14...
            ax = 1.0_r2 / x
            b  = 2.0_r2 * ax**2
            ss = (0.0_r2, 0.0_r2)
            s3 = (0.0_r2,-1.0_r2)
            an = 3.0_r2

            ! define the number for subroutine aa2 [Loskutov (1971)]
            y   = sqrt( RI * conjg(ri) )  *  x
            num = 1.25 * y + 15.5

            if      ( y<1.0_r2 ) then
            num = 7.5 * y + 9.0
            else if ( (y>100.0_r2) .and. (y<50000.0_r2) ) then
            num = 1.0625 * y + 28.5
            else if ( y>=50000.0_r2 ) then
            num=1.005*y+50.5
            end if

            if(num > nterms) then
            ier = 2
            print *, "<!> Error in subroutine shexqnn2:"
            print *, "    - Maximum number of terms  : ", nterms
            print *, "    - Number of terms required : ", num
      print *, 
     *  "** Solution: Increase default value of the variable 'nterm' **"
            else
            ! logarithmic derivative to Bessel function (complex argument)
            call aa2(ax,ri,num,ru)

            ! ------------------------------------------------------------------------------------------
            ! FIRST TERM
            ! ------------------------------------------------------------------------------------------
            ! initialize term counter
            iterm = 1

            ! Bessel functions
            ass = sqrt( pi / 2.0_r2 * ax )
            w1  = 2.0_r2/pi * ax
            Si  = sin(x)/x
            Co  = cos(x)/x

            ! n=0
            besJ0 =  Si / ass
            besY0 = -Co / ass
            iu0   = 0

            ! n=1
            besJ1 = ( Si * ax - Co) / ass
            besY1 = (-Co * ax - Si) / ass
            iu1   = 0
            iu2   = 0

            ! Mie coefficients
            s   = ru(1) / ri + ax
            s1  = s * besJ1 - besJ0
            s2  = s * besY1 - besY0
            ra0 = s1 / (s1 - s3 * s2)   ! coefficient a_1

            s   = ru(1) * ri + ax
            s1  = s * besJ1 - besJ0
            s2  = s * besY1 - besY0
            rb0 = s1 / (s1 - s3 * s2)   ! coefficient b_1

            ! efficiency factors
            r    = -1.5_r2 * (ra0-rb0)
            Qext = an * (ra0 + rb0)
            Qsca = an * (ra0 * conjg(ra0)  +  rb0 * conjg(rb0))

            ! scattering amplitude functions
            if (doSA) then
                  do iang=1, nang
                  mu(iang) = cos( (real(iang,kind=r2)-1.0_r2) * 
     >             (pi/2.0_r2)/real(nang-1,kind=r2) )
                  end do

                  fpi0(:) = 0.0_r2
                  fpi1(:) = 1.0_r2
                  SA1(:)  = cmplx( 0.0_r2, 0.0_r2 )
                  SA2(:)  = cmplx( 0.0_r2, 0.0_r2 )

                  r_iterm = real(iterm,kind=r2)  ! double precision
                  fac     = (2.0*r_iterm + 1.0_r2) / 
     >             (r_iterm * (r_iterm+1.0_r2))

                  do iang=1, nang
                  iang2      = 2 * nang - iang

                  fpi(iang)  = fpi1(iang)
                  ftau(iang) = r_iterm * mu(iang) * fpi(iang)  - 
     >              (r_iterm+1.0) * fpi0(iang)

                  P          = (-1.0)**(iterm-1)
                  SA1(iang)  = SA1(iang)   +   fac * 
     >             (ra0*fpi(iang)  + rb0*ftau(iang))

                  T          = (-1.0)**iterm
                  SA2(iang)  = SA2(iang)   +   fac * 
     >             (ra0*ftau(iang) + rb0*fpi(iang) )

                  if  ( iang /= iang2 )  then
                        SA1(iang2) = SA1(iang2)   +   
     >                   fac * (ra0*fpi( iang)*P + rb0*ftau(iang)*T)
                        SA2(iang2) = SA2(iang2)   +   
     >                   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                  end if
                  end do

                  iterm   = iterm + 1
                  r_iterm = real(iterm, kind=r2)

                  do iang=1, nang
                  fpi1(iang) = ((2.0*r_iterm-1.0) / 
     >             (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                  fpi1(iang) = fpi1(iang)   -   
     >             r_iterm * fpi0(iang)/(r_iterm-1.0)
                  fpi0(iang) = fpi(iang)
                  end do
            else
                  ! start value for the next terms
                  iterm = 2
            end if

            ! ------------------------------------------------------------------------------------------
            ! 2., 3., ... num
            ! ------------------------------------------------------------------------------------------
            z = -1.0_r2

            do
                  an  = an + 2.0_r2
                  an2 = an - 2.0_r2

                  ! Bessel functions
                  if(iu1 == iu0) then
                  besY2 = an2 * ax * besY1 - besY0
                  else
                  besY2 = an2 * ax * besY1 - besY0 / factor
                  end if
                  if(dabs(besY2) > 1.0e+300_r2) then
                  besY2 = besY2 / factor
                  iu2   = iu1 + 1
                  end if
                  besJ2 = (w1 + besY2 * besJ1) / besY1

                  ! Mie coefficients
                  r_iterm = real(iterm,kind=r2)

                  s   = ru(iterm) / ri + r_iterm * ax
                  if(iu1>1) then
                  ier=1
                  return
                  endif
                  if(iu2>1) then
                  ier=1
                  return
                  endif
                  s1  = s * besJ2 / fact(iu2) - 
     >             besJ1 / fact(iu1) ! Subscript #1 of the array FACT has value 2 which is greater than the upper bound of 1
                  s2  = s * besY2 * fact(iu2) - 
     >             besY1 * fact(iu1)
                  ra1 = s1 / (s1 - s3 * s2)                        ! coefficient a_n, (n=iterm)

                  s   = ru(iterm) * ri + r_iterm * ax
                  s1  = s * besJ2 / fact(iu2) - 
     >             besJ1 / fact(iu1)
                  s2  = s * besY2 * fact(iu2) - 
     >             besY1 * fact(iu1)
                  rb1 = s1 / (s1 - s3 * s2)                        ! coefficient b_n, (n=iterm)

                  ! efficiency factors
                  z  = -z
                  rr = z * (r_iterm + 0.5_r2) * (ra1 - rb1)
                  r  = r + rr
                  ss = ss + (r_iterm - 1.0_r2) * 
     >             (r_iterm + 1.0_r2) / r_iterm * (ra0 * conjg(ra1)  
     >                   + rb0 * conjg(rb1)) 
     >            + an2 / r_iterm / (r_iterm - 1.0_r2) * 
     >            (ra0 * conjg(rb0))
                  qq   = an * (ra1 + rb1)
                  Qext = Qext + qq
                  Qsca = Qsca + an * (ra1 * conjg(ra1) 
     >             + rb1 * conjg(rb1))

                  ! leaving-the-loop criterion
                  if ( dabs(qq / qext) < eps ) then
                  exit
                  end if

                  ! Bessel functions
                  besJ0 = besJ1
                  besJ1 = besJ2
                  besY0 = besY1
                  besY1 = besY2
                  iu0   = iu1
                  iu1   = iu2
                  ra0   = ra1
                  rb0   = rb1

                  ! scattering amplitude functions
                  if (doSA) then
                  r_iterm = real(iterm,kind=r2)
                  fac      = (2.0 * r_iterm+1.0) / 
     >             (r_iterm * (r_iterm+1.0))

                  do iang=1, nang
                        iang2      = 2 * nang - iang

                        fpi(iang)  = fpi1(iang)
                        ftau(iang) = r_iterm * mu(iang) * fpi(iang) 
     >                    -  (r_iterm+1.0) * fpi0(iang)

                        P          = (-1.0)**(iterm-1)
                        SA1(iang)  = SA1(iang)   +   fac * 
     >                   (ra0*fpi(iang) + rb0*ftau(iang))

                        T          = (-1.0)**iterm
                        SA2(iang)  = SA2(iang)   +   fac * 
     >                  (ra0*ftau(iang) + rb0*fpi(iang))

                        if  ( iang /= iang2 ) then
                        SA1(iang2) = SA1(iang2)   +   
     >                   fac * (ra0*fpi(iang)*P  + rb0*ftau(iang)*T)
                        SA2(iang2) = SA2(iang2)   +   
     >                   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                        end if
                  end do

                  iterm   = iterm + 1
                  r_iterm = real(iterm,kind=r2)

                  do iang=1, nang
                        fpi1(iang) = ((2.0*r_iterm-1.0) / 
     >                   (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                        fpi1(iang) = fpi1(iang)   -   
     >                   r_iterm * fpi0(iang)/(r_iterm-1.0)
                        fpi0(iang) = fpi(iang)
                  end do
                  else
                  iterm = iterm + 1
                  endif

                  if ( iterm==num ) then
                  exit
                  else
                  cycle
                  end if
            end do

            ! efficiency factors (final calculations)
            Qext   = b * Qext
            Qsca   = b * Qsca
            Qbk    = 2.0_r2 * b * r * conjg(r)
            Qpr    = Qext - 2.0_r2 * b * ss
            Qabs   = Qext - Qsca
            albedo = Qsca / Qext
            g      = (Qext - Qpr) / Qsca
            end if
      end if
      !deallocate( ru, mu, fpi, fpi0, fpi1, ftau )
      ier=0
      return
      end subroutine shexqnn2
      end module mie_routines

      module mie_opacity
      use mie_precision
      use mie_data
      use mie_routines
      implicit none

      complex(kind=dp), allocatable, 
     > dimension(:,:) :: e_inc, N_inc
      real(kind=dp), parameter :: 
     > vol2rad = (3.0_dp/(4.0_dp*pi))**(1.0_dp/3.0_dp)

      contains

      subroutine calc_Mie(LL,N_eff,kext_dust,g_dust,
     > a_dust,kabs_dust,ksca_dust)
    !! Uses MieX large size parameter routines: Wolf & Voshchinnikov (2004)
    !! References: Helling et al. (2008), Lee et al. (2015b), Lee et al. (2016)
      
      implicit none
      include 'parameter.inc'
      
      real(kind=dp), dimension(4), intent(in) :: LL
      complex(kind=dp), dimension(nwl), intent(in) ::  N_eff
      real(kind=dp), dimension(nwl), intent(out) :: 
     > kext_dust, g_dust, a_dust, kabs_dust, ksca_dust
      integer, parameter :: rnang = 1
      integer :: rier, l, nwtot, kos_step
      real(kind=dp) :: x, rQext, rQsca, rQabs, rQbk, rQpr, 
     > ralbedo, rg, a, cross_sec_g
      complex(kind=dp), dimension(rnang) :: rSA1, rSA2
      complex(kind=dp) :: N_efft
      logical, parameter :: rdoSA = .False.
      real(kind=dp) :: wnos, conos, wlos, wlstep
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot

      
      ! Use mean (0) or effective (1) cloud particle radius [cm]
      if (a_type == 0) then
            a = (LL(2)/LL(1)) * vol2rad
      else if (a_type == 1) then
            a = (LL(4)/LL(3)) * vol2rad
      end if

      ! Cross sectional area * g-1 [cm2 g-1]
      cross_sec_g = (LL(1))*pi*a**2

      ! Mie Theory Step----------------------------
      do l = 1, nwtot

            ! Size parameter - limit to 1e-6 and 10000 (micron)(2000 for low memory)
            x = (pi2*(a*1.0000e+8_dp))/wlos(l)
            x = max(1.00001e-2_dp, x)
            !x = min(1.0000e+8_dp, x)
            ! if (x<=1.0e-2) then
            ! kabs_dust(l) = 0.0
            ! ksca_dust(l) = 0.0
            ! else 
            x = min(1.0000e+8_dp, x)
            N_efft = N_eff(l)

            !print*, N_efft, x
            ! Mie theory routine - careful with memory in parallel and large n_wl
            call shexqnn2(N_efft, x, rQext, rQsca, rQabs, 
     >      rQbk, rQpr, ralbedo, rg, 
     >       rier, rSA1, rSA2, rdoSA, rnang)

            kext_dust(l) = cross_sec_g * rQext
            a_dust(l) = ralbedo
            g_dust(l) = rg
            kabs_dust(l) = cross_sec_g * rQabs
            ksca_dust(l) = cross_sec_g * rQsca
            !end if

      end do

      end subroutine calc_Mie

      subroutine calc_emt(Vol,N_eff)
      !! Effective medium theory routines
      !! Calculates (n,k) constants of mixed material cloud particle
      !! Uses the Bruggeman method with LLL method backup
      !! References: Helling et al. (2008), Lee et al. (2015b), Lee et al. (2016)
      
      implicit none
      include 'parameter.inc'
      integer :: n, l, j, n_dust, nwtot, kos_step
      integer :: n_lay, n_eps
      real(kind=dp) :: wnos, conos, wlos, wlstep
      real(kind=dp), dimension(max_inc),intent(in) :: Vol
      complex(kind=dp), dimension(nwl) :: e_eff
      complex(kind=dp), dimension(nwl), intent(out) ::  N_eff
      complex(kind=dp) :: e_eff0, N_eff0
      logical :: errflag = .False.
      logical :: first_call = .True.
      complex(kind=dp) :: cri_inc
      complex(kind=dp), allocatable, dimension(:,:) :: e_inc, N_inc
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      real(kind=dp) :: dabstable, dscatable, eps, temp, pgas,
     * rhod, rho_sw, rho_dtg, z_sw, z_marcs_init, 
     * z_marcs, pg_read, tt_init, eps_init, eps_new,
     * r_null, f_eps, f_opac

      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps

      common /cdustnew/ cri_inc(max_inc, nwl), n_dust

      ! On first call, allocate local global arrays
      if (first_call .eqv. .True.) then
            allocate(e_inc(n_dust, nwtot), N_inc(n_dust, nwtot))
            do l = 1, nwtot
            do n = 1, n_dust
            N_inc(n,l) = cri_inc(n,l)
            e_inc(n,l) = m2e(N_inc(n,l)) ! N_inc**2
            end do
            end do

            first_call = .False.
      end if

      ! Main wavelength loops
      
      do l = 1, nwtot
            ! Effective Medium Theory (EMT) step
            N_eff0 = (0.0_dp,0.0_dp)
            do n = 1, n_dust
            N_eff0 = N_eff0 + Vol(n) * N_inc(n,l)
            end do

            ! Try Bruggeman approach
            if (Brug .eqv. .True.) then
            ! Call Newton-Raphson minimization
            call NR(N_inc(1:n_dust,l),Vol(1:n_dust),
     >       N_eff0,N_eff(l),errflag)

            ! if fails (errflag) use LLL method
            if (errflag .eqv. .True.) then
            e_eff0 = (0.0_dp,0.0_dp)
            do n = 1, n_dust
                  e_eff0 = e_eff0 + Vol(n)*e_inc(n,l)**(onethird)
            end do
            e_eff(l) = e_eff0**3
            N_eff(l) = e2m(e_eff(l))
            errflag = .False.
            end if
            ! Try LLL method
            else if (Brug .eqv. .False.) then
            e_eff0 = (0.0_dp,0.0_dp)
            do n = 1, n_dust
                  e_eff0 = e_eff0 + Vol(n)*e_inc(n,l)**(onethird)
            end do
            e_eff(l) = e_eff0**3
            N_eff(l) = e2m(e_eff(l))
            end if

      end do


      end subroutine calc_emt

      ! ------------- Functions for LLL theory ------------- !!
      pure complex(kind=dp) function e2m(e)
      
      implicit none
      
      real(kind=dp) :: ereal, eimag, n, k
      real(kind=dp) :: sqrte2
      complex(kind=dp), intent(in) :: e

      ereal = real(e,kind=dp)
      eimag = aimag(e)
      sqrte2 = sqrt(ereal*ereal + eimag*eimag)
      n = sqrt(0.5_dp * ( ereal + sqrte2))
      k = sqrt(0.5_dp * (-ereal + sqrte2))
      e2m = cmplx(n,k,kind=dp)
      end function e2m


      pure complex(kind=dp) function m2e(m)
      
      implicit none

      real(kind=dp) :: ereal, eimag, n, k
      complex(kind=dp), intent(in) :: m

      n = real(m,kind=dp)
      k = aimag(m)
      ereal = n*n - k*k
      eimag = 2.0_dp * n * k
      m2e = cmplx(ereal,eimag,kind=dp)
      end function m2e

      ! -------------------------------------------------------------------------------------------------
      ! A program for function minimization using the Newton-Raphson method.
      ! -------------------------------------------------------------------------------------------------
      pure subroutine NR(M_inc,V_inc,M_eff0,M_eff,unphysical)
      implicit none

      integer, parameter :: itmax = 30
      integer :: it
      real(kind=dp), dimension(:), intent(in)  :: V_inc
      real(kind=dp), dimension(2,2) :: DF
      real(kind=dp), dimension(2) ::  Fold, FF, FF1, FF2, FF3, FF4
      real(kind=dp), dimension(2) :: corr, xx, xnew
      real(kind=dp) :: qual, de1, de2
      complex(kind=dp), dimension(:), intent(in) :: M_inc
      complex(kind=dp), intent(in)  :: M_eff0
      complex(kind=dp), intent(out) :: M_eff
      logical, intent(inout) :: unphysical

      M_eff = M_eff0

      do it = 1, itmax
            xx(1) = real(M_eff,kind=dp)
            xx(2) = aimag(M_eff)
            call Bruggeman(M_eff,M_inc,V_inc,FF)
            qual = FF(1)*FF(1) + FF(2)*FF(2)
            de1 = xx(1)*1.0e-5_dp
            de2 = xx(2)*1.0e-5_dp
            call Bruggeman(M_eff
     >       +cmplx(de1,0.0_dp,kind=dp),M_inc,V_inc,FF1)
            call Bruggeman(M_eff
     >       -cmplx(de1,0.0_dp,kind=dp),M_inc,V_inc,FF2)
            call Bruggeman(M_eff
     >       +cmplx(0.0_dp,de2,kind=dp),M_inc,V_inc,FF3)
            call Bruggeman(M_eff
     >       -cmplx(0.0_dp,de2,kind=dp),M_inc,V_inc,FF4)
            DF(1,1) = (FF1(1)-FF2(1)) / (2.0_dp*de1)
            DF(1,2) = (FF3(1)-FF4(1)) / (2.0_dp*de2)
            DF(2,1) = (FF1(2)-FF2(2)) / (2.0_dp*de1)
            DF(2,2) = (FF3(2)-FF4(2)) / (2.0_dp*de2)
            Fold = FF
            call gauss(2,2,DF,corr,FF)
            corr = -corr
            call eff_pullback(2,xx,corr,Fold,
     >       xnew,unphysical,M_inc,V_inc)
            if (unphysical) then
      !        print*, qual, 'qual'
            exit
            end if
            M_eff = cmplx(xnew(1),xnew(2),kind=dp)
            if (abs(qual) < 1.0e-13_dp) then
            exit
            end if
      end do

      end subroutine NR


      !!! Combine using Bruggeman formula
      ! define function to be minimized
      pure subroutine Bruggeman(M_eff,M_inc,V_inc,FF)
      implicit none

      integer :: i
      real(kind=dp), dimension(:), intent(in) :: V_inc
      real(kind=dp), intent(out) :: FF(2)
      complex(kind=dp), intent(in) :: M_eff
      complex(kind=dp), dimension(:), intent(in) :: M_inc
      complex(kind=dp) :: fun, mm2, mmi2

      mm2 = M_eff**2
      fun = cmplx(0.0_dp,0.0_dp,kind=dp)
      do i = 1, size(V_inc)
            mmi2 = M_inc(i)**2
            fun = fun + V_inc(i)*(mmi2 - mm2)/(mmi2 + 2.0_dp*mm2)
      end do

      FF(1) = real(fun,kind=dp)
      FF(2) = aimag(fun)

      end subroutine Bruggeman

      ! Pull Back Eff
      pure subroutine eff_pullback(N,xx,dx,Fold,
     > xnew,unphysical,M_inc,V_inc)
      implicit none

      integer, intent(in) :: N
      integer, parameter :: itmax = 20
      integer :: it
      real(kind=dp), dimension(:),intent(in) :: V_inc
      real(kind=dp), dimension(N), intent(in) :: xx,Fold
      real(kind=dp), dimension(N), intent(inout) :: dx
      real(kind=dp), dimension(N), intent(out) :: xnew
      real(kind=dp), dimension(2) :: Fnew(2)
      real(kind=dp) :: fac,qold,qnew
      complex(kind=dp), dimension(:),intent(in) :: M_inc
      complex(kind=dp) :: mm
      logical, intent(out) :: unphysical

      qold = Fold(1)*Fold(1)+Fold(2)*Fold(2)
      fac = 1.0_dp

      do it = 1, itmax
            xnew = xx + fac*dx

            if ((xnew(1) > 0.0_dp).and.(xnew(2) > 0.0_dp)) then
            mm = cmplx(xnew(1),xnew(2),kind=dp)
            call Bruggeman(mm,M_inc,V_inc,Fnew)
            qnew = Fnew(1)*Fnew(1)+Fnew(2)*Fnew(2)
            !write(*,*) it,qold,qnew
            unphysical = .False.
            if (qnew < qold) then
            exit
            end if
            else
            !write(*,*) it,"negative (n,k)",xnew
            unphysical = .True.
            endif

            fac = fac*0.7_dp
      enddo

      end subroutine eff_pullback

      !**********************************************************************
      pure subroutine gauss(Nd,N,a,x,b)
      implicit none
      !**********************************************************************
      !*****                                                            *****
      !*****   Diese Routine loesst ein lineares Gleichungssystem       *****
      !*****   der Form    (( a )) * ( x ) = ( b )     nach x auf.      *****
      !*****   Der Algorithmus funktioniert, indem die Matrix a         *****
      !*****   auf Dreiecksform gebracht wird.                          *****
      !*****                                                            *****
      !*****   EINGABE:  Nd = Dimension der Vektoren, der Matrix        *****
      !*****              N = Dimension der Gl-Systems (N<=Nd)          *****
      !*****              a = (N x N)-Matrix                            *****
      !*****              b = (N)-Vektor                                *****
      !*****   AUSGABE:   x = (N)-Vektor                                *****
      !*****                                                            *****
      !**********************************************************************
      !*
      integer, intent(in) :: Nd, N
      integer :: i, j, k, kmax
      real(kind=dp), dimension(Nd,Nd), intent(inout) :: a
      real(kind=dp), dimension(Nd), intent(inout) :: b
      real(kind=dp), dimension(Nd), intent(out) :: x
      real(kind=dp) :: c, amax

      do i = 1, N-1
      !*       ------------------------------------------
      !*       ***  MAX-Zeilentausch der i-ten Zeile  ***
      !*       ------------------------------------------
            kmax = i
            amax = abs(a(i,i))
            do k = i+1, N
            if (abs(a(k,i)) > amax) then
            amax = abs(a(k,i))
            kmax = k
            endif
            end do

            if (kmax /= i) then
            do j = 1, N
            c = a(i,j)
            a(i,j) = a(kmax,j)
            a(kmax,j) = c
            end do
            c = b(i)
            b(i) = b(kmax)
            b(kmax) = c
            end if
      !*
      !*       ---------------------------------
      !*       ***  bringe auf Dreiecksform  ***
      !*       ---------------------------------
            do k = i+1, N
            c = a(k,i) / a(i,i)
            a(k,i) = 0.0_dp
            do j = i+1, N
            a(k,j) = a(k,j) - c * a(i,j)
            end do
            b(k) = b(k) - c * b(i)
            end do
      !*
      end do
      !*
      !*     --------------------------
      !*     ***  loese nach x auf  ***
      !*     --------------------------
      do i = N, 1, -1
            c = 0.0_dp
            if (i < N) then
            do j = i+1, N
            c = c + a(i,j) * x(j)
            end do
            end if
            x(i) = (b(i) - c) / a(i,i)
      end do

      end subroutine gauss
      end module mie_opacity

!-----------------------------------------------------------------------
! Calculates the dust opacity for each wavelength at each depth layer 
! by interpolating the dust opacity table. 
! Using interpolation scheme from petitRadtrans. 
! BCE 2022
!-----------------------------------------------------------------------
      subroutine dust_opac_eps_interp

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'

      ! internal variables)
      DOUBLE PRECISION:: kappa_abs_0(nwl), 
     * kappa_abs_1(nwl), kappa_sca_0(nwl),
     * kappa_sca_1(nwl), abs_r(nwl),
     * sca_r(nwl), eps_0(max_eps),
     * eps_1(max_eps), eps_r(max_eps),
     * rhod_0(max_eps), rhod_1(max_eps)
      double precision:: kappa_cloud, kappa_cloud_old, 
     * edepletion, kappa_cloud_diff, epsilon_cloud,
     * kappa_cloud_abs, kappa_cloud_int
      dimension edepletion(max_eps, ndp)
      dimension kappa_cloud_diff(ndp), kappa_cloud_old(ndp)
      dimension kappa_cloud_abs(ndp)
      dimension kappa_cloud_int(ndp)
      dimension pgmarcs(ndp), ptot(ndp)
      dimension a1(nwl), a2(nwl), s1(nwl), s2(nwl)
      dimension rho(ndp)
      dimension icloud_diff(ndp)
      COMMON /CSTYR/MIHAL,NOCONV
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     * vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     * ntau,iter
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      common/ci9/ai(16)
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      common/cabinit/abinit(natms),kelem(natms),nelem
      common/ci5/abmarcs(18,ndp),anjon(18,5),h(5),part(18,5),dxi,
     *           f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      common/ci1/fl2(5),parco(45),parq(180),shxij(5),tparf(4),
     *  xiong(16,5),eev,enamn(ndp),sumh(ndp),xkbol,nj(16),iel(16),
     *  summ(ndp),nel
      common /cdusteps/ ielnr(max_eps)
      common /cdustindex/ iabinit(max_eps), iabmarcs(max_eps), 
     * isw(max_eps)
      common /dpeset/ dpein,dtin, pe_corr(ndp)
      dimension sum(ndp)
      dimension pg(ndp),diff(ndp,max_lay),imin(ndp),itemp(12),
     *      imin2(ndp) 
      dimension f_abs(nwl), f_sca(nwl)
      dimension diff_Ps_vals(max_lay), diff_Ts_vals(max_lay)
      dimension mx_elm(18)
      data eev/1.602095e-12/,xmh/1.67339e-24/ 
      data mx_elm /1, 2,6,7,8,10,11,12,13,14,16,19,20,23,26,27,22,17/
            
      dust_abs(1:ntau,1:nwtot) = 0.
      dust_sca(1:ntau,1:nwtot) = 0.
      
      if (idustopac==1) then
      open(unit=976, file='f_cloud.in')
      read(976,*) f_opac
      print*, "f is ", f_opac
      close(976)

      open(unit=389, file='./data/abinit_abmarcs.dat')
      read(389,*)
      do n=1, max_eps !dependent on amount of elements in abinit_abmarcs.dat
        read(389, '(12x,i2,8x,i2,4x,i2)') 
     >       iabinit(n), iabmarcs(n), isw(n) 
      
      end do
      close(389)  
      end if
      do i=1,18
            abmarcs(i,1:ntau) = abinit(mx_elm(i))
      end do
      istruc_len = ntau
      if (n_lay > 0) then
      p_min = MINVAL(pgas(1:n_lay))
      p_max = MAXVAL(pgas(1:n_lay))

      
      if (idustopac ==1) then
            pgmarcs(1:istruc_len) = pg_read(1:istruc_len) 
      else 
            do k=1,istruc_len
                  ptot(k)=pp(k)-ppr(k)-ppt(k)+pe_corr(k)   
            end do
            pgmarcs(1:istruc_len) = ptot(1:istruc_len)
      end if
      icloud_count = 0

      do i_str=1, istruc_len
        call search_intp_ind(pgas, n_lay, pgmarcs(i_str), 1, m)
        if (pgmarcs(i_str) <= pgas(m) .and. (m==1)) then
           
            dust_abs(i_str, 1:nwtot) = 
     >      dabstable(m, 1:nwtot)  

            dust_sca(i_str, 1:nwtot) =   
     >      dscatable(m, 1:nwtot) 


            eps_new(1:n_eps, i_str) = 
     >      eps(1:n_eps, m) 


        else if ((pgmarcs(i_str)>= pgas(m)) .and. 
     >         (pgmarcs(i_str)<= pgas(m+1))) then

            kappa_abs_0(1:nwtot) = dabstable(m, 1:nwtot)
            kappa_abs_1(1:nwtot) = dabstable(m+1, 1:nwtot)

            kappa_sca_0(1:nwtot) = dscatable(m, 1:nwtot)
            kappa_sca_1(1:nwtot)= dscatable(m+1, 1:nwtot)

            eps_0(1:n_eps) = eps(1:n_eps, m)
            eps_1(1:n_eps) = eps(1:n_eps, m+1)

            p0 = pgas(m)
            p1 = pgas(m+1)
            pmarcs = pgmarcs(i_str)

            p1_p0 = p1/p0 
            p1_p = p1/pmarcs
            p_p0 = pmarcs/p0

            dust_abs(i_str, 1:nwtot) = 
     >       kappa_abs_0(1:nwtot)**(log10(p1_p)/log10(p1_p0)) * 
     >       kappa_abs_1(1:nwtot)**(log10(p_p0)/log10(p1_p0))

            dust_sca(i_str, 1:nwtot) =   
     >       kappa_sca_0(1:nwtot)**(log10(p1_p)/log10(p1_p0)) * 
     >       kappa_sca_1(1:nwtot)**(log10(p_p0)/log10(p1_p0))

            eps_new(1:n_eps, i_str) = 
     >       eps_0(1:n_eps)**(log10(p1_p)/log10(p1_p0)) * 
     >       eps_1(1:n_eps)**(log10(p_p0)/log10(p1_p0))

            else if (pgmarcs(i_str)>pgas(m+1)) then
            if (icloud_count==0) then
            ncloud = i_str-1
            icloud_count = 1
            end if
            dust_abs(i_str,1:nwtot) = 0.0
            dust_sca(i_str,1:nwtot) = 0.0
            eps_new(1:n_eps, i_str) = -1.0

            end if

        
      end do
      else if (n_lay <= 0) then
      print*, "0 layers in DRIFT - no cloud!"
      do i_str=1, istruc_len
            dust_abs(i_str,1:nwtot) = 0.0
            dust_sca(i_str,1:nwtot) = 0.0
            eps_new(1:n_eps, i_str) = -1.0
      end do
      ncloud = 0
      !idust =0
      end if

      do i=1, ntau
        do m=1, n_eps
         do n=1, max_eps
            if (ielnr(m)==isw(n)) then
              eps_init_el = 10.0**(abinit(iabinit(n))-12.0)
              if (eps_new(m,i) <0.0) then
                  eps_new(m,i) = eps_init_el
                 
              end if             
            exit
            end if   
          end do
        end do
      end do

      do m=1, n_eps
      epsilon_cloud(m, 1:ntau) = eps_new(m,1:ntau)
      end do
      
      
      if (idustopac==1) then
      do i=1, ntau
      kappa_cloud(i, 1:nwtot) = 0.0
      do j=1, nwtot
            kappa_cloud(i,j) = 
     *       dust_abs(i,j) + dust_sca(i,j)
      !print*, dust_abs(i,j) + dust_sca(i,j)
     
      end do
      end do
      kappa_cloud_int(1:ntau) = 0.0
      do i=1, nwtot
      kappa_cloud_int(1:ntau) = kappa_cloud_int(1:ntau)+
     *  kappa_cloud(1:ntau, i)*wlstep(i)*1.0e-4
      end do
c       open(unit=0302, file='kappa_cloud_lay.dat', status="replace", 
c      *       position="append", action="write")
c       !print*, "Integrated cloud opacity"
c       do i=1, ntau
      
c       write(0302, *) kappa_cloud_int(i)
c       end do
c       close(0302)

      open(unit=880, file='driftmarcs_it.in')
      read(880,*) it_driftmarcs
      close(880)

      end if
      
      !if (icloud_conv == 0) then

      call cloud_opac(f_opac)

      do i=1, ntau
        do m=1, n_eps
         do n=1, max_eps
            if (ielnr(m)==isw(n)) then
              abmarcs(iabmarcs(n),i) = 
     *         log10(epsilon_cloud(m,i)) + 12.0  
            exit
            end if   
          end do
        end do
      end do
      
      
      sum(1:ntau)=0.0
      xmy(1:ntau)=0.0
      summ(1:ntau)=0.0
      sumh(1:ntau)=0.0
      enamn(1:ntau)=0.0
      do i=1,16
        abmarcs(i,1:ntau)=10.**abmarcs(i,1:ntau)
        sum(1:ntau)=sum(1:ntau)+abmarcs(I,1:ntau)
      end do
      abmarcs(17,1:ntau)=10.**abmarcs(17,1:ntau)
      abmarcs(18,1:ntau)=10.**abmarcs(18,1:ntau)
      aha=abmarcs(1,1)
      do i=1,16
        abmarcs(i,1:ntau)=abmarcs(i,1:ntau)/aha
        summ(1:ntau)=summ(1:ntau)+abmarcs(i,1:ntau)
        xmy(1:ntau)=xmy(1:ntau)+abmarcs(i,1:ntau)*ai(i)
      end do
      abmarcs(17,1:ntau)=abmarcs(17,1:ntau)/aha
      abmarcs(18,1:ntau)=abmarcs(18,1:ntau)/aha
      xmy(1:ntau)=xmy(1:ntau)/ai(1)
      sumh(1:ntau)=sum(1:ntau)/aha-1.
      summ(1:ntau)=summ(1:ntau)-abmarcs(1,1:ntau)-abmarcs(3,1:ntau)-
     *  abmarcs(4,1:ntau)-abmarcs(5,1:ntau)
     
      enamn(1:ntau) = eev/(xmh*xmy(1:ntau))
      
      idustopac = 0

      do i=1, ntau
      kappa_cloud(i, 1:nwtot) = 0.0
      do j=1, nwtot
            kappa_cloud(i,j) =  
     *       dust_abs(i,j) + dust_sca(i,j)
      end do
      end do
      kappa_cloud_int(1:ntau) = 0.0
      do i=1, nwtot
      kappa_cloud_int(1:ntau) = kappa_cloud_int(1:ntau)+
     *  kappa_cloud(1:ntau, i)*wlstep(i)*1.0e-4
      end do
      open(unit=0307, file='kappa_cloud_lay.dat', status="replace", 
     *       position="append", action="write")

      do i=1, ntau
      write(0307, *) kappa_cloud_int(i)
      end do
      close(0307)

      open(unit=0407, file='dust_abs.dat', status="replace", 
     *       position="append", action="write")

      do i=1, ntau
      write(0407, *) dust_abs(i, 1:nwtot)
      end do
      close(0407)
      open(unit=0507, file='dust_sca.dat', status="replace", 
     *       position="append", action="write")

      do i=1, ntau
      write(0507, *) dust_sca(i, 1:nwtot)
      end do
      close(0507)
      open(unit=0707, file='epsilon_cloud.dat', status='replace',
     * position='append', action='write' )
      do m=1, n_eps
      write(0707, *) ielnr(m), epsilon_cloud(m,1:ntau)
      end do
      close(0707)    

      if (it_driftmarcs>1) then
      print*, 'Checking for opacity differences.'
      
      open(0302,file='kappa_cloud_lay_old.dat',
     &  status='old',readonly)
      do i=1, ntau
            read(0302,*) kappa_cloud_old(i)
      end do
      close(0302)
      !do i=1, ncloud
      !print*, i
      !print*, 'new ', kappa_cloud_int(i)
      !print*, 'old ', kappa_cloud_old(i)
      !end do
      do i=1, ncloud
      
       if (kappa_cloud_int(i)>0.0) then
       mag1 = nint(log10(kappa_cloud_int(i)))
       else
       mag1 = 0
       end if
      if (kappa_cloud_old(i)>0.0) then
      mag2 = nint(log10(kappa_cloud_old(i)))
      else
      mag2=0
      end if
      
      if (kappa_cloud_old(i)<=1.0e-20)  then
        if (mag1 == mag2 .or. mag1 == mag2-1 .or. 
     *   mag1 == mag2+1) then
          kappa_cloud_diff(i) = 0.0
          icloud_diff(i) = 0
        else
          if(kappa_cloud_old(i)==0.0) then
            kappa_cloud_old(i)=1.0e-30
            kappa_cloud_diff(i) = 
     * ((kappa_cloud_old(i) - kappa_cloud_int(i)))/
     * kappa_cloud_old(i)
          else
          kappa_cloud_diff(i) = 
     * ((kappa_cloud_old(i) - kappa_cloud_int(i)))/
     * kappa_cloud_old(i)
          end if
      !     if (kappa_cloud_diff(i)<=1.0e-1) then
      !       icloud_diff(i) = 0
      !       else
      !       icloud_diff(i) = 1
      !     end if
        end if 
      else 
      !if (i==1) then
      kappa_cloud_diff(i) = 
     * ((kappa_cloud_old(i) - kappa_cloud_int(i)))/
     * kappa_cloud_old(i)
      !if (mag1==mag2) then
      !icloud_diff(i) = 0
      !else
      !icloud_diff(i) = 1
      !end if

      !else
      !kappa_cloud_diff(i) = 
      !* (abs(kappa_cloud_int(i) - kappa_cloud_old(i)))/
      !* kappa_cloud_old(i)
      !if (kappa_cloud_diff(i)<=1.0e-1) then
      !icloud_diff(i) = 0
      !else
      !icloud_diff(i) = 1
      !end if
      end if

      !else
      !kappa_cloud_diff(i) = 
      !* (abs(kappa_cloud_int(i) - kappa_cloud_old(i)))/
      !* kappa_cloud_old(i)
      !if (kappa_cloud_diff(i)<=5.0e-2) then
      !icloud_diff(i) = 0
      !else
      !icloud_diff(i) = 1
      !end if
      !end if
      end do
      
      ! do i=1, ntau
      ! print*, i, kappa_cloud_diff(i)
      ! end do
      do i=1, ntau
      kappa_cloud_abs(i) = abs(kappa_cloud_diff(i))
      end do

      !idiff_kappa = maxval(abs(icloud_diff(1:ntau)))
      index_kappa = maxloc(kappa_cloud_abs(1:ntau), dim=1)
      print*, "Max cloud opacity difference was ", 
     * kappa_cloud_diff(index_kappa)
      open(unit=990, file='kappa_diff.dat', status='replace')
      write(990,*) kappa_cloud_diff(index_kappa)
      close(990)
      end if
      return
      end
!-----------------------------------------------------------------------
      subroutine cloud_opac(f_dust)
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      integer:: i, ielnr, n_eps, ntau, nwtot
      real*8 :: kappa_cloud, epsilon_cloud,f_dust
      real*8 :: epsilon_cloud_old, 
     * dust_abs, dust_sca

      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     & vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     & ntau,iter
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot

      
      open(0207,file='dust_abs_old.dat',
     &  status='old',readonly)
      do i=1, ntau
            read(0207,*) dust_abs_old(i, 1:nwtot)
      end do
      close(0207)
      open(0107,file='dust_sca_old.dat',
     &  status='old',readonly)
      do i=1, ntau
            read(0107,*) dust_sca_old(i, 1:nwtot)
      end do
      close(0107)

      open(0507,file='epsilon_cloud_old.dat',
     &  status='old',readonly)
      do i=1, n_eps
            read(0507,*) ielnr, epsilon_cloud_old(i, 1:ntau)
      end do
      close(0507)
      do i=1, ntau

      dust_abs(i, 1:nwtot) = 
     *   f_dust* dust_abs(i, 1:nwtot) + 
     *   (1.0-f_dust)*dust_abs_old(i,1:nwtot)

      dust_sca(i, 1:nwtot) = 
     *   f_dust* dust_sca(i, 1:nwtot) + 
     *   (1.0-f_dust)*dust_sca_old(i,1:nwtot)
      end do
      
      do i=1, n_eps
      epsilon_cloud(i, 1:ntau) = 
     *       f_dust*epsilon_cloud(i,1:ntau) +
     *       (1.0-f_dust)*epsilon_cloud_old(i, 1:ntau)
      end do

      return
      end 
!-----------------------------------------------------------------------
      subroutine search_intp_ind(binbord,binbordlen,arr,arrlen,intpint)

      implicit none

      INTEGER            :: binbordlen, arrlen, intpint(arrlen)
      DOUBLE PRECISION   :: binbord(binbordlen),arr(arrlen)
      INTEGER            :: i_arr
      INTEGER            :: pivot, k0, km

      ! carry out a binary search for the interpolation bin borders

      do i_arr = 1, arrlen

      if (arr(i_arr) >= binbord(binbordlen)) then
            intpint(i_arr) = binbordlen - 1
      else if (arr(i_arr) <= binbord(1)) then
            intpint(i_arr) = 1
      else

            k0 = 1
            km = binbordlen
            pivot = (km+k0)/2

            do while(km-k0>1)

            if (arr(i_arr) >= binbord(pivot)) then
                  k0 = pivot
                  pivot = (km+k0)/2
            else
                  km = pivot
                  pivot = (km+k0)/2
            end if

            end do

            intpint(i_arr) = k0

      end if

      end do
      return
      end


!-----------------------------------------------------------------------
! Writes a MARCS output file to be read by DRIFT
! Juncher 2015
!-----------------------------------------------------------------------
      subroutine marcs2drift
      
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      character flag(ndp)
      dimension pg(ndp),surfgrav(ndp),v(ndp),emu(ndp),rad(ndp)
      dimension flip_rad(ndp), abundances(ndp,100), pg2(ndp)
      dimension kappa_cloud_int(ndp)
      real*8 :: kappa_cloud, kappa_cloud_int
      common /tauc/tau(ndp),dtauln(ndp),jtau
      common /cg/grav,konsg /cteff/teff,flux
      common /masse/relm
      common /mixc/palfa,pbeta,pny,py /cvfix/vfix
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     & vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     & ntau,iter
      common /cstyr/mihal,noconv
      common /rossc/xkapr(ndp),cross(ndp)
      common /cabinit/abinit(natms),kelem(natms),nelem
      common /tsuji/ nattsuji,nmotsuji,parptsuji(500)
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdrift/ idust, ieps, idustopac
      common /cmolrat/ fold(ndp,8),molold,kl
      common /cmetpe/ ppel(ndp), metpe
      common /cdustopac/ dust_abs(ndp,nwl), dust_sca(ndp,nwl),
     *      dust_abs_old(ndp,nwl), dust_sca_old(ndp,nwl),
     *      kappa_cloud(ndp,nwl),epsilon_cloud(max_eps,ndp),
     *      epsilon_cloud_old(max_eps,ndp)
      common /cdustindex/ iabinit(max_eps), iabmarcs(max_eps), 
     * isw(max_eps)
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      open(unit=389, file='./data/abinit_abmarcs.dat')
      read(389,*)
      do n=1, max_eps !dependent on amount of elements in abinit_abmarcs.dat
        read(389, '(12x,i2,8x,i2,4x,i2)') 
     >       iabinit(n), iabmarcs(n), isw(n) 
      
      end do
      close(389)

      open(unit=2873, file='nlay_nwtot.in', status="replace", 
     *       position="append", action="write")
      write(2873,*) ntau
      write(2873,*) nwtot
      close(2873)

      !print*, "at marcs2drift"

      sun_rad = 6.96342e10   ! cm
      
! Radius
      relr = sqrt(relm/grav*10**(4.44))

      flip_rad(1:ntau) = 0.
      do k=1,ntau
        kl = k
        furem = fure
        !call termo(k,tt(k),ppel(k),ppr(k),ptot,rro,cp,cv,agrad,q,u2)
        fure = 1./(xkapr(k)*ro(k))
        if(k .eq. 1) cycle
        flip_rad(k) = flip_rad(k-1) + (tau(k)-tau(k-1))*(fure+furem)*0.5
      end do
      rad(1:ntau) = 0.
      do k=1,ntau
        rad(k) = flip_rad(ntau-k+1)+relr*sun_rad
      end do
      
! Gas pressure
      do k=1,ntau
        kl = k
      if(metpe.eq.1) then
        call jon(tt(k),ppe(k),1,pgx,rox,dumx,0)
      else if(metpe.eq.2) then
        pgx =  pp(k)-ppr(k)-ppt(k)
        !call jon(tt(k),ppel(k),1,pgx,rox,dumx,0)
      end if
        pg(k) = PGx
      end do

! Surface gravity
      do k=1,ntau
        surfgrav(k) = (relr*sun_rad/rad(k))*grav
      end do 

! Convective velocity
      do k=1,ntau
        if(k .gt. 1) go to 13
        v(k) = 0.
        go to 15
        
13      if(k .eq. ntau) go to 14
        ya=(tau(k)-tau(k-1))/(tau(k+1)-tau(k-1))
        yb=1.-ya
        v(k)=ya*vv(k+1)+yb*vv(k)
        if(vv(k).gt.0..and.vv(k+1).gt.0.) v(k)=
     &    exp(ya*log(vv(k+1))+yb*log(vv(k)))
        go to 15
        
14      continue
        ya=(2.*tau(k)-tau(k-1)-tau(k-2))/(tau(k)-tau(k-2))
        yb=1.-ya
        v(k)=ya*vv(k)+yb*vv(k-1)
15      continue
      end do

! Convection flag
      do k=1,ntau
        if(v(k) .eq. 0.) then
          flag(k) = 'F'
        else
          flag(k) = 'T'
        end if
      end do   
      
! Mean molecular mass
      do k=1,ntau
        kl = k
        !call termo(k,tt(k),ppel(k),ppr(k),ptot,rro,cp,cv,agrad,q,u2)
        emu(k) = (1.38*ro(k)*tt(k))/(1.67e-8*pg(k))

      end do  
      
! Save to file
      
      !print*, "writting marcs2drift file"
      open(unit=33, file='marcs2drift.dat')
      
      write(33,'(a2)') ' !'
      write(33,'(a43)')' ! MARCS output file to be read in by DRIFT'
      if(idust .eq. 0) then
        write(33,'(a2)') ' ! Dust not included'
      else
        write(33,'(a2)') ' ! Dust included'
      end if
      write(33,'(a2)') ' !'
      write(33,'(a32,a19)') ' ! Model parameters: Teff, logg,',
     *  ' mixing, overshoot:'
      write(33,'(f12.3,f13.3,f13.3,f13.3)') teff, log10(grav), 
     *  palfa, 2.200
      write(33,'(a2)') ' !'
      write(33,'(a31)') ' ! Number of atmosphere layers:'
      write(33,'(i5)') ntau
      write(33,'(a2)') ' !'
      write(33,'(a42)') ' ! Number of elements in abundances table:'
      write(33,'(i5)') nelem
      write(33,'(a2)') ' !'
      write(33,'(a32)') ' ! Z of the considered elements:'
      do i=1,nelem,8
        if(i .gt. nelem-8) then
          write(33,'(8(i5.2,1x))') (kelem(j), j=i,nelem)
        else
          write(33,'(8(i5,1x))') (kelem(j), j=i,i+7)
        end if
      end do      
      write(33,'(a2)') ' !'
      write(33,'(a5,6a16,a6,a16)') ' !  #', 'Rad [cm]', 
     *  'Temp [K]', 'Pgas [dyn cm-2]', 'Ro [g cm-3]', 'g [cm s-2]',
     *  'v_conv [cm s-1]', 'Flag', 'mu [amu]'
      do k=1,ntau
        write(33,'(i5,6e16.8,a6,e16.8)') k, rad(k), tt(k), pg(k),
     *  ro(k), surfgrav(k), v(k), flag(k), emu(k)
      end do
      write(33,'(a2)') ' !'
      write(33,'(a41)') ' ! Initial Element abundances for each Z:'
      do i=1,nelem-1,8
        if(i .gt. nelem-9) then
          write(33,'(8f6.2)') (abinit(j), j=i,nelem-1)
        else
          write(33,'(8f6.2)') (abinit(j), j=i,i+7)
        end if
      end do
      close(33)
      if (idust == 0) then
      open(unit=0307, file='kappa_cloud_lay.dat', status="replace", 
     *       position="append", action="write")

      do i=1, ntau
      kappa_cloud_int(i) = 0.0
      write(0307, *) kappa_cloud_int(i)
      end do
      close(0307)

      open(unit=0407, file='dust_abs.dat', status="replace", 
     *       position="append", action="write")

      
      do i=1, ntau
      dust_abs(i, 1:nwtot) = 0.0
      write(0407, *) dust_abs(i, 1:nwtot)
      end do
      close(0407)
      open(unit=0507, file='dust_sca.dat', status="replace", 
     *       position="append", action="write")

      do i=1, ntau
      dust_sca(i, 1:nwtot) = 0.0
      write(0507, *) dust_sca(i, 1:nwtot)
      end do
      close(0507)

      do m=1, max_eps
              eps_init_el = 10.0**(abinit(iabinit(m))-12.0)
              epsilon_cloud(m,1:ntau) = eps_init_el            
      end do


      open(unit=0707, file='epsilon_cloud.dat', status='replace',
     * position='append', action='write' )
      do m=1, max_eps
      write(0707, *) isw(m), epsilon_cloud(m,1:ntau)
      end do
      close(0707)    
      end if
      return 
      end

    

************************************************************************
c routine that calls GGChem
c updated jan 2021
c BCE
c Takes as arguments:
c k - the atmosphere layer number
c temp - temperature of the atmospheric layer (Kelvin)
c pgas - gas pressures at the layer
c
c Input for GGchem is written in marcs2ggchem.in file.
c GGchem is called with model dimension 0, with a single temperature
c and pressure. For now it is also called with a typical solar
c abundance initially.
c 
c marcs2gg_presmo.dat and marcs2gg_part.dat have correspondence
c between molecular ids in MARCS and GGchem for the PRESMO array
c and the PARTP/PARTPP arrays.
c 
c GGchem_ppel contains overall output from GGchem (it is written in
c ggchems main.f ) such as total electron pressure, total atom pressure
c and total molecular pressure. 
!
! Updated by BCE oct 22 to compute atomic fractions for continuum
! correctly
!
************************************************************************

      subroutine GGCHEM(k,temp,pgas)
      
      implicit real*8 (a-h, o-z)
      include 'parameter.inc'     
      character*2, dimension(17) :: atmarcsnames
      character*2, dimension(83) :: atnames_gg
      character*4 abcname
      real*8, dimension(48) :: abundforgg
      real*8,intent(in) :: temp,pgas
      real*8 :: pphsum
      character atnames*2, molnames*8, molnames2*4
      character(len=2) :: Al, He, Si, C, Mg
      integer, dimension(500) :: Al_id, He_id, Si_id, 
     > C_id, Mg_id,H_id
      real*8, dimension(500) :: Al_x, He_x, Si_x, C_x, 
     > Mg_x, H_x 
      common /ggchemresults/
     > tgk,pgesk,ppelGG,ggmuk,ggrhok,ppsumk,ppappsumk,ppnonappsumk,
     > ppat1sumk,ppat2sumk,ppmolsumk,ppgsk,rhon_total, f1gg, f5gg,
     > rCgg, rMggg, rAlgg, rSigg, rHegg
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)     
     >                ,atnames(22),molnames(543),molnames2(75)
      common /ggchembool/ iggcall
      
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      common/ci5/abmarcs(18,ndp),anjon(18,5),h(5),part(18,5),dxi,
     *           f1,f2,f3,f4,f5,xkhm,xmh,xmy(ndp)
      common/ci1/fl2(5),parco(45),parq(180),shxij(5),tparf(4),
     *  xiong(16,5),eev,enamn(ndp),sumh(ndp),xkbol,nj(16),iel(16),
     *  summ(ndp),nel
      common/cabinit/abinit(natms),kelem(natms),nelem
      common/cabnames/abcname(natms)
      logical ggchem_mol(maxosmol), ggchem_index_read
      integer ggchem_index(maxosmol), molno
      character(len=5) molnames_new(maxosmol)
      common /molupdate/ molnames_new, 
     * ggchem_index,
     * molno, ggchem_mol, ggchem_index_read
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp_sw(max_lay),pgas_sw(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdusteps/ ielnr(max_eps)
      common /cdustindex/ iabinit(max_eps), iabmarcs(max_eps), 
     * isw(max_eps)

      bar=1.Q+6
      
      
      do n=1, 48
            abundforgg(n) = abinit(n)
      end do
      
      
      if (idust ==1) then
      do m=1, n_eps
         do n=1, max_eps
            if (ielnr(m)==isw(n)) then
            abundforgg(iabinit(n)) = 
     >       log10(abmarcs(iabmarcs(n), k))+12.
             
              exit
            end if
          end do
      end do
      end if

      open(unit=546, file='abund_drift.in', status='replace')
      do n=1, 48
            if (n <36  .or. n==44 ) then
            write(546,'(a5f6.2)') (abcname(n),abundforgg(n))
            
            end if
      end do
      close(546)
      
      open(unit=70, file='marcs2ggchem.in', status='replace')
      write(70, '(19a)') '# selected elements'
      write(70, '(99a)') 'H He C N O Na Mg Si F Fe Al Ca Cr Ti S Cl K Li
     & V Zr Be P el'
      write(70, '(38a)') '# name of files with molecular kp-data'
      write(70,'(a54)') 'dispol_BarklemCollet.dat             
     &   ! dispol_file '
      write(70,'(a54)')
     & 'dispol_StockKitzmann_withoutTsuji.dat   ! dispol_file2'
      write(70,'(a54)') 'dispol_WoitkeRefit.dat         
     &         ! dispol_file3'
      write(70,'(a64)') '# abundance options 1=EarthCrust, 2=Ocean,
     & 3=Solar, 4=Meteorites'
      write(70,'(a34)') '0                     ! abund_pick'
      write(70,'(a34)') 'abund_drift.in                    '
      write(70,'(a30)') '0                    ! verbose'
      write(70,'(a27)') '# equilibrium condensation?'
      write(70,'(a36)') '.false.               ! model_eqcond'   
      write(70,'(a42)') '0                     ! model_dim  (0,1,2)'
      write(70,'(a36)') '.true.                ! model_pconst'
      write(70,*) temp, '                ! Tmax [K]'
      write(70,*) pgas/bar,'                   ! pmax [bar]' 
      close(70)

      call system('./GGchem/ggchem marcs2ggchem.in > ggchem_out.txt')

      open(unit=990,file='GGchem_ppel')
        read(990,*) Tg,pgesk,ppelGG,ggmuk,ggrhok,ggrhodust,ppsumk
     &     ,ppappsumk,ppnonappsumk,ppat1sumk,ppat2sumk,ppmolsumk
     &     ,ppgsk,rhon_total
      close(990)
      open(unit=809, file='./marcs2gg_partpp.dat')
             read(809,124) 
     * ((idmarcspart(m), idggchempart(m),molnames2(m)), m=1, 75)
124          format(i4,3x,i3,4x,a4)
      close(809)      
        open(unit=321, file='ndensity.dat')
        
            read(321,*) (rhonallat(k,m),m=1,22)
            read(321,*) (rhonallmol(k,m), m=1,543)
        close(321)

      
        open(unit=707,file='pp.dat')
              read(707,*) (ppallat(k,m),m=1,22)
              read(707,*) (ppallmol(k,m),m=1,543)
              
                read(707,708) ((km,atnames(m)),m=1,22)
                read(707,*)
                read(707,709) (molnames(m),m=1,543)


708             format(i4,18x,a2)
709             format(10a8)
              
        close(707)
        
        

        Al = 'AL'
        He = 'HE'
        Si = 'SI'
        Mg = 'MG'
        C = 'C'

        iAl = 1
        iHe = 1
        iSi = 1
        iMg = 1
        iC = 1
        iH = 1
        

        do m=1, 543
            if (index(molnames(m), 'H') /= 0 ) then
              if (index(molnames(m), "HE") == 0 ) then
              indH = index(molnames(m), 'H')
              if (molnames(m) == 'SI(CH3)4') then
                  H_x(iH) = 12.0
              else if (molnames(m)(indH:indH+1) == 'H2') then
                  H_x(iH) = 2.0
              else if (molnames(m)(indH:indH+1) == 'H3') then
                  H_x(iH) = 3.0
              else if (molnames(m)(indH:indH+1) == 'H4') then
                  H_x(iH) = 4.0
              else if (molnames(m)(indH:indH+1) == 'H5') then
                  H_x(iH) = 5.0
              else if (molnames(m)(indH:indH+1) == 'H6') then
                  H_x(iH) = 6.0
              else if (molnames(m)(indH:indH+1) == 'H7') then
                  H_x(iH) = 7.0
              else if (molnames(m)(indH:indH+2) == 'H)2') then
                  H_x(iH) = 2.0
              else
                  H_x(iH) = 1.0
              end if

               H_id(iH) = m  
               iH = iH +1
              end if
            end if

            if (index(molnames(m), Al) /= 0 ) then
              if (index(molnames(m), "AL2") /= 0 ) then
                  Al_x(iAl) = 2.0
              else if (index(molnames(m), ")2") /= 0 ) then 
                  Al_x(iAl) = 2.0
              else 
                  Al_x(iAl) = 1.0
              end if
              Al_id(iAl) = m  
              iAl = iAl +1
               
            end if
            if (index(molnames(m), He) /= 0 ) then
              if (index(molnames(m), "HE2") /= 0 ) then
                  He_x(iHe) = 2.0
              else if (index(molnames(m), ")2") /= 0 ) then 
                  He_x(iHe) = 2.0
              else 
                  He_x(iHe) = 1.0
              end if
               He_id(iHe) = m
              iHe = iHe +1
              
            end if
            if (index(molnames(m), Mg) /= 0 ) then
              if (index(molnames(m), "MG2") /= 0 ) then
                  Mg_x(iMg) = 2.0
              else 
                  Mg_x(iMg) = 1.0
              end if
                  
              Mg_id(iMg) = m 
              iMg = iMg +1
            end if
            if (index(molnames(m), Si) /= 0 ) then
              if (index(molnames(m), "SI2") /= 0 ) then
                  Si_x(iSi) = 2.0
              else 
                  Si_x(iSi) = 1.0
              end if

              Si_id(iSi) = m    
              iSi = iSi +1
               
            end if
            if (index(molnames(m), 'C') /= 0 ) then
              indC = index(molnames(m), 'C')
              if (molnames(m) == 'CRC2') then
                  C_x(iC) = 2.0
                  C_id(iC) = m
                  iC = iC +1
              else if ((molnames(m)(indC:indC+1) /= 'CL') 
     &         .and. (molnames(m)(indC:indC+1) /= 'CA')          
     &         .and. (molnames(m)(indC:indC+1) /= 'CR')) 
     &          then
              if (molnames(m)(indC:indC+1) == 'C5') then
                  C_x(iC) = 5.0
              else if (molnames(m)(indC:indC+3) == 'CO)5') then
                  C_x(iC) = 5.0
              else if (molnames(m)(indC:indC+1) == 'C4') then 
                  C_x(iC) = 4.0
              else if (molnames(m)(indC:indC+4) == 'CH3)4') then 
                  C_x(iC) = 4.0
              else if (molnames(m)(indC:indC+1) == 'C3') then 
                  C_x(iC) = 3.0
              else if (molnames(m)(indC:indC+1) == 'C2') then
                  C_x(iC) = 2.0
              else if (molnames(m)(indC:indC+3) == 'CN)2') then
                  C_x(iC) = 2.0
              else
                  C_x(iC) = 1.0
              end if
              C_id(iC) = m
              iC = iC +1
              end if
            end if
        end do 
        
        
        pp_C_tot = ppallat(k,3)
        do m =1, iC-1
            pp_C_tot = pp_C_tot + C_x(m)*ppallmol(k,C_id(m))
        end do        
        rCgg = ppallat(k,3) / pp_C_tot

        pp_Mg_tot = ppallat(k,7)
        do m =1, iMg-1
            pp_Mg_tot = pp_Mg_tot + Mg_x(m)*ppallmol(k,Mg_id(m))
        end do       
        rMggg = ppallat(k,7) / pp_Mg_tot 

        pp_Si_tot = ppallat(k,8)
        do m =1, iSi-1
            pp_Si_tot = pp_Si_tot + Si_x(m)*ppallmol(k,Si_id(m))
        end do        
        rSigg = ppallat(k,8) / pp_Si_tot
        
        pp_Al_tot = ppallat(k,11)
        do m =1, iAl-1
            pp_Al_tot = pp_Al_tot + Al_x(m)*ppallmol(k,Al_id(m))
        end do        
        rAlgg = ppallat(k,11) / pp_Al_tot

        pp_He_tot = ppallat(k, 2)
        do m =1, iHe-1
            pp_He_tot = pp_He_tot + He_x(m)*ppallmol(k,He_id(m))
        end do        

        rHegg = ppallat(k,2) / pp_He_tot

        phydrototal = ppallat(k,1) 
        do m=1, iH-1
         phydrototal = phydrototal + H_x(m)*ppallmol(k,H_id(m))
        end do
        f1gg = ppallat(k,1) / phydrototal
        f2gg = ppallmol(k,446) / phydrototal
        f3gg = ppallmol(k,447) / phydrototal
        f4gg = ppallmol(k,15) / phydrototal
        f5gg = ppallmol(k,1) / phydrototal

C       do n =1,imarcs2gg
C              gg_partpp(k,idmarcspart(n)) = ppallmol(k,idggchempart(n))
C       enddo
      
        iggcall = 1
      return
      end

!-----------------------------------------------------------------------
! Reads in the dust data from the DRIFT output file and makes a table
! of dust opacities
! Juncher 2015
! BCE 2023 updats to use new mie routine
!-----------------------------------------------------------------------
      subroutine drift2marcs
      use mie_precision,only: sp,dp,qp
      use mie_data
      use mie_opacity
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      real*8    :: L0(max_lay), L1(max_lay), 
     * L2(max_lay), L3(max_lay), object
      !real*8, parameter :: pi = 3.14159265
      character (len=50) :: dust_names, dust_files
      character (len=50) :: find_dust, dname
      character(99) get_path
      character(:), allocatable :: path
      complex(kind=dp) :: cri_inc
      logical   :: first
      complex(kind=dp), dimension(nwl) :: N_eff
      real(kind=dp), dimension(4) :: LL
      real(kind=dp),dimension(nwl) :: k_ext_lay, g_lay, a_lay,
     > kabs_dust, ksca_dust

      dimension :: V_inc(max_inc,max_lay),
     *             a(max_lay)
      dimension :: tt_read(ndp),rad_read(ndp),
     *  ro_read(ndp), surfgrav_read(ndp), v_read(ndp), 
     *  flag_read(ndp), emu_read(ndp)


      dimension :: wn(2000),p(2),step(2),var(2)
      common /statec/ppr(ndp),ppt(ndp),pp(ndp),gg(ndp),zz(ndp),dd(ndp),
     * vv(ndp),ffc(ndp),ppe(ndp),tt(ndp),tauln(ndp),ro(ndp),
     * ntau,iter
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      common /cdustdata/ dabstable(max_lay,nwl),dscatable(max_lay,nwl),
     *    eps(max_eps,max_lay),temp(max_lay),pgas(max_lay), 
     *    rhod(max_lay), rho_sw(max_lay), rho_dtg(ndp),
     *    z_sw(max_lay),z_marcs_init(ndp), 
     *    z_marcs(ndp), pg_read(ndp), tt_init(ndp),
     *    eps_init(max_eps, ndp),  eps_new(max_eps, ndp), 
     *      r_null,f_eps,f_opac, n_lay, n_eps
      common /cdustdata2/ dust_names(max_inc), dust_files(max_inc)
      common /cdusteps/ ielnr(max_eps)
      common /cdrift/ idust, ieps, idustopac, icloud_conv
      common /cdustnew/ cri_inc(max_inc,nwl), n_dust

! Read in dust data from DRIFT
      open(110,file='./out_default/drift2marcs.dat',
     &  status='old',readonly)
      do i=1,50
        read(110,'(A)') find_dust
        if (find_dust(5:17) == ' dust species') then
            read (find_dust(:4),'(i4)') n_dust
            exit
        end if
      end do
      print*, "Number of dust species: ", n_dust
      if(n_dust .gt. max_inc) then
        print *, 'Error: increase max_inc.'
        stop
      end if
      do i=1,n_dust
        read(110,*) dust_names(i)
        !print*, dust_names(i)
      end do
      read(110,'(i4)') n_eps
      if (n_eps>max_eps) then 
      print*, 'Error: increase max_eps'
      stop
      end if
      print*, "Number of affected elements by cloud formation: ", n_eps
      !if(n_eps .gt. max_eps) then
       ! print *, 'Error: increase max_eps.'
       ! stop
      !end if
      do i=1,n_eps
        read(110,*) ielnr(i)
      end do
      read(110,*)
      read(110,*)
      i = 1
      do
        read(110,'(8e20.12,99e20.12)',iostat=io) 
     *    z_sw(i),temp(i),rho_sw(i),pgas(i),L0(i),L1(i),
     *    L2(i), L3(i), rhod(i), V_inc(1:n_dust,i),
     *    eps(1:n_eps,i)
        if(io .lt. 0) exit
        if(L0(i) .le. 0.) then
            exit
        else
          a(i) = (3./4./pi)**(1./3.)*L1(i)/L0(i) 
        end if
        i = i + 1
      end do
      close(110)
      n_lay = i-1
      print*, "DRIFT layers #", n_lay
      open(unit=189, file='dust_file.dat', readonly)
      do n=1, max_inc
            read(189,'(A)') dname
            !print*, 'dust in dust file ', trim(dname)
            n_on = 0
            do i=1, n_dust
              if (trim(dname)== dust_names(i)) then
                  !print*, "dust present ", dust_names(i)
                  read(189, '(A)') get_path
                  path= trim(get_path)
                  dust_files(i) = path
                  n_on =1
                  exit
              end if
            end do
            if (n_on == 0) read(189, *) 
      end do
      close(189)

      open(unit=33, file='marcs2drift.dat')
      do i=1,22
            read(33,*)
      end do
      do k=1,ntau
        
        read(33,'(i5,6e16.8,a6,e16.8)') k_read, rad_read(k), 
     *  tt_read(k), pg_read(k), ro_read(k),
     *  surfgrav_read(k), v_read(k), flag_read(k), emu_read(k)
      if (k==1) r_null = rad_read(k)
      
      end do
      close(33)
      if(n_lay .gt. max_lay) then
        print *, 'Warning: only used the first 1000 layers from DRIFT.'
      end if
      
! Read in optical constants
      
      call optical_data(n_dust, cri_inc)

! Make opacity table
      do j=1,n_lay
      !calculate effective medium theory
      call calc_emt(V_inc(:,j), N_eff)
      LL(1) = L0(j)
      LL(2) = L1(j)
      LL(3) = L2(j)
      LL(4) = L3(j)
      call calc_Mie(LL,N_eff,k_ext_lay,g_lay,
     > a_lay, kabs_dust, ksca_dust)
      dabstable(j, :) = kabs_dust
      dscatable(j, :) = ksca_dust
      end do
     

      return
      end

      !-----------------------------------------------------------------------
! Reads in optical data for solids
! Juncher 2015
! BCE 2023 updated for different interpolation/extrapolation
!-----------------------------------------------------------------------
      subroutine optical_data(n_dust, cri_inc)
      use mie_precision,ONLY: sp,dp,qp
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      integer :: i, j, k, io, u, nlines, counter, n_dust, n_lines
      real*8, allocatable, dimension(:) :: wl_work, n_work, k_work
      real*8  :: a, b, nkdata(3,3000)
      real*8  :: ndata(max_inc,nwl), kdata(max_inc,nwl)
      complex(kind=dp) :: cri_inc(max_inc, nwl)
      character :: filename(16)*50
      character(len=200) :: info
      logical :: conducting,ex
      character (len=50) :: dust_names, dust_files
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      common /cdustdata2/ dust_names(max_inc), dust_files(max_inc)
      !common /cdustnew/ cri_inc(max_inc,nwl), n_dust
      
      do i=1, n_dust
        u = 6785
        !print*, dust_files(i)
        open(unit=u,file=dust_files(i))
        
        read(u,*) n_lines, conducting
        read(u,*) ; read(u,'(A200)') info; 
     *  read(u,*) ; read(u,*)
        allocate(wl_work(n_lines),n_work(n_lines),k_work(n_lines))

        do l = 1, n_lines
            read(u,*) wl_work(l),n_work(l),k_work(l)
            !print*, wl_work(l),n_work(l),k_work(l)
            n_work(l) = max(0.0D+0, n_work(l))
            k_work(l) = max(0.0D+0, k_work(l))
            if (dust_names(i) == 'CaSiO3[s]') then
            n_work(l) = 1.0
            k_work(l) = 0.0
            end if
        enddo
        close(u)

        wl_work(:) = wl_work(:)*1.0e+4
        do l = 1, nwtot
            !if (i==4) print*, wlos(l)
            ! If required wavelength is less than available data - keep constant
         if (wlos(l) < wl_work(1)) then
            !if (i==4) print*, "first wl smaller for SiO2"
            ndata(i,l) = n_work(1)
            kdata(i,l) = k_work(1)
            ! If required wavelength is greater than available data - extrapolate
            ! Non conducting: n is constant - k is linear decreasing
            ! Conducting: n and k are log-log extrapolated
         else if (wlos(l) > wl_work(n_lines)) then
            !if (i==4) print*, "extrapolation!"
            if (conducting .eqv. .False.) then
                ndata(i,l) = n_work(n_lines)
                kdata(i,l) = k_work(n_lines)*wl_work(n_lines)/wlos(l)
            else if (conducting .eqv. .True.) then
                do l1 = n_lines,1,-1
                  if (wl_work(l1) < 0.7D+0*wl_work(n_lines)) then  ! data can be noisy, so
                    exit                                           ! it's safer to use larger
                  endif                                            ! region to get the slope
                enddo
                fac = log(wlos(l)/wl_work(n_lines))/
     *           log(wl_work(l1)/wl_work(n_lines))
                ndata(i,l) = exp(log(n_work(n_lines)) 
     &            + fac*log(n_work(l1)/n_work(n_lines)))
                kdata(i,l) = exp(log(k_work(n_lines)) 
     &            + fac*log(k_work(l1)/k_work(n_lines)))
            endif
              ! Data is availible in the required wavelength range - log-log interpolation
              
       else
              ! Loop across work arrays untill straddle point is point then interpolate
            do l1 = 1, n_lines - 1
              if (wlos(l) >= wl_work(l1) .and. 
     *         wlos(l) <= wl_work(l1+1)) then
                  fac = log(wlos(l)/wl_work(l1))/
     *             log(wl_work(l1+1)/wl_work(l1))
                  ndata(i,l) = exp(log(n_work(l1)) 
     *             + fac*log(n_work(l1+1)/n_work(l1)))
                  if (k_work(l1) <= 0.0D+0 .or. 
     *             k_work(l1+1) <= 0.0) then
                    kdata(i,l) = 0.0D+0
                  else
                    kdata(i,l) = exp(log(k_work(l1)) 
     *                      + fac*log(k_work(l1+1)/k_work(l1)))
                  endif
                  exit
              endif
            enddo
       endif
      enddo
      deallocate(wl_work,n_work, k_work)
      
      end do
      do i=1,nwtot
        cri_inc(1:n_dust,i) = 
     >   dcmplx(ndata(1:n_dust,i),kdata(1:n_dust,i))
      end do
      
      end
!--------------------------------------------
      subroutine krome_solve(ntau,T,ptot)
!--------------------------------------------
      !TODO/TO CONSIDER
      !what if some speciesnames have capitals in them
      !what if speciesnames are ordered differently

      use krome_main !use krome (mandatory)
      use krome_user !use utility (for krome_idx_* constants and others)

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'

C      implicit none
      integer,parameter::nsp=krome_nmols !number of species (common) 
      integer::ntau, i,j,k, istep,header_size,buffer
      integer::ss_istep(ntau)
      integer,dimension(22,2) :: index_at  !indices of atomic specs, first index marcs, second index krome
      integer,dimension(543,2):: index_mol !indices of mol specs, first index marcs, second index krome   
      integer,dimension(75,2) :: index_mol2!indices of mol specs, first index marcs, second index krome        
      integer,dimension(1000) :: index_not !indices of species NOT found in marcs mols and atoms, only krome index needed
      integer:: M_index
      integer:: atom_counter,mol_counter,mol2_counter,not_counter !counters for different kinds of species
      logical:: not_found
      logical:: first_call = .True.
      logical:: first_not_call = .True.
      integer:: krome_photo_on
      real*8::Tgas,dt,num_den(ntau,nsp)
      real*8::R, R_cgs,Na, Pcon(ntau), T(ntau), ptot(ndp)
      real*8::num_den_mol(ndp,543), num_den_at(ndp,22)
      real*8::mix_rat(ntau,nsp)
      real*8::num_den_cont(100000,nsp),time_cont(100000)
      real*8::time,dtmax,dt_inc
      real*8::ss_time(ntau)
      real*8::krome_tmax
      real*8::conv_crit
      logical::use_conv = .False.
      logical::is_conv
      real*8,dimension(nwreal)::FLUX_RAD_eV
      real*8::photo_bins_high,photo_bins_low,photo_bins_nominator
      real*8::aa_to_m_conv,J_to_eV_conv,HC_to_SI_conv
      character(len=100)::spec_name
      character(len=8),dimension(nsp)::chem_spec
      character atnames*2, molnames*8, molnames2*4    
      real*8 :: krome_photo_scale
      integer ,parameter:: output_freq=100
      common/cos/wnos(nwl),conos(ndp,nwl),wlos(nwl),wlstep(nwl),
     *    kos_step,nwtot
      COMMON/COSWR/osresl
      COMMON /CPF/PF,PFE,PFD,FIXROS,ITSTOP
      LOGICAL PF,PFE,PFD,FIXROS,ITSTOP
      COMMON /NATURE/BOLTZK,CLIGHT,ECHARG,HPLNCK,PI,PI4C,RYDBRG,
     *STEFAN
      common /cit/it,itmax
      common /ggchempp/ppallat(ndp,22),ppallmol(ndp,543)
     >                ,rhonallat(ndp,22),rhonallmol(ndp,543)
     >                ,gg_partpp(ndp,400)
     >                ,ppat(22),ppmol(543)
     >                ,idmarcspres(32),idggchempres(32)
     >                ,idmarcspart(75),idggchempart(75)
     >                ,atnames(22),molnames(543),molnames2(75)
      common /noneq/ krome_on,krome_photo_on,krome_photo_scale
      common /noneq_time/ dt_start,dt_max,dt_inc,krome_tmax
      common /noneq_output/ krome_output,krome_debug,krome_return
      common /photochem/ FLUX_RAD(ndp,nwreal) !second dimension should be nwtot, in most cases 7949

      call krome_init() !init krome (mandatory)

      R = 8.31446261815324 !Gas constant in m^3 Pa K^-1 mol^-1
      R_cgs = 8.31446261815324E-3
      Na = 6.02214076D23 !Avogadros number in mol^-1
      if (use_conv .eqv. .True.) then
            conv_crit=1e-16
      endif

      !Convert pressures in dyne/cm^2 to number densities in molecules/cm^3 
      !write relevant species out from info.log
      if (first_call.eq..True.) then !initialization of krome. Checking for molecules in MARCS and finding their indices, setting Photobins if photochem is needed
      header_size=5 !info.log header size
      open(unit=12,file='./krome/MARCS_build/info.log',status = 'old')     
      do i=1,nsp+header_size
        read(unit=12,fmt='(A100)') spec_name
        if (i-header_size.lt.10) buffer=3
        if ((i-header_size.ge.10).and.(i-header_size.lt.100)) buffer=4
        if ((i-header_size.ge.100)) buffer=5
        if (i.gt.header_size) then      
         do j=buffer,len(chem_spec)+buffer                                      
          if (spec_name(j:j)=='k') exit !k starts the id part of the species in info.log                        
          write(chem_spec(i-header_size)(j-(buffer-1):j-(buffer-1)),
     >       '(A1)') spec_name(j:j)
         enddo       
         chem_spec(i-header_size)=chem_spec(i-header_size)(1:j-buffer-1)!cut down empty ends of the string 
        endif 
      enddo
      close(12)
      !if (first_call.eq..True.) then
        !write(*,*) nsp,'Species are found'
        write(*,*) 'The following',nsp,
     >   'species are found in your krome build'
        write(*,*) chem_spec(1:nsp)
        !first_call = .False.
      !endif      
      !check if species is in atomnames, molnames, molnames2 and safe indices
      !check also if speices not in MARCS database
      atom_counter=0
      mol_counter=0
      mol2_counter=0
      not_counter=0
      
      !write(*,*) atomnames
      !write(*,*) molnames
      !write(*,*) molnames2
      do i=1,krome_nmols
       not_found=.True.
       do k=1,22 !number of atoms    
        if (trim(atnames(k)).eq.trim(chem_spec(i))) then
         atom_counter=atom_counter+1        
         index_at(atom_counter,1)=k
         index_at(atom_counter,2)=i
         not_found=.False.
         cycle
        endif
       enddo
       do k=1,543 !number of molecues in molnames
        if (trim(molnames(k)).eq.trim(chem_spec(i))) then
         mol_counter=mol_counter+1    
         index_mol(mol_counter,1)=k
         index_mol(mol_counter,2)=i
         not_found=.False.     
         cycle    
        endif
       enddo 
       do k=1,75 !number of molecues in molnames2
        if (trim(molnames2(k)).eq.trim(chem_spec(i))) then
         mol2_counter=mol2_counter+1      
         index_mol2(mol2_counter,1)=k
         index_mol2(mol2_counter,2)=i
         not_found=.False.    
         cycle     
        endif        
       enddo
       if (trim(chem_spec(i)).eq.'M') then
        !write(*,*) 'M identified as', i
        M_counter=1
        M_index = i
        not_found=.False.
        cycle
       endif  
       if (not_found==.True.) then
        not_counter=not_counter+1
        index_not(not_counter)= i
        if (first_not_call.eq. .True.) then
        write(*,*) trim(chem_spec(i)),
     >  ' has not been found in MARCS' 
        endif
       endif  
      enddo

      if (not_counter.gt.0) then
       if (first_not_call.eq..True.) then
       write(*,*) 'Species not found in MARCS will be set to ',
     > 'default abundance of 1E-20'  
       first_not_call=.False.
       endif
      endif
      if (krome_photo_on.eq.1) then !initializing of photobins
         J_to_eV_conv=6.242E18!Joule to eV
         HC_to_SI_conv=1E-9 !convert hplanck and clight to SI from CGS
         aa_to_m_conv=1E-10 !Angstrom to meters for wavelenths
         photo_bins_nominator=HPLNCK*CLIGHT*HC_to_SI_conv
         photo_bins_low=photo_bins_nominator/(WLOS(nwtot)*aa_to_m_conv)
         photo_bins_high=photo_bins_nominator/(WLOS(1)*aa_to_m_conv)
         photo_bins_low=photo_bins_low*J_to_eV_conv
         photo_bins_high=photo_bins_high*J_to_eV_conv
         call krome_set_photobinE_log(photo_bins_low,photo_bins_high)
         
       if (krome_debug.eq.1) then
         open(unit=5656,file='krome_bins_mid.dat')
         open(unit=5959,file='krome_bins_delta.dat')
         open(unit=4242,file='krome_bins_photoJ.dat')
         open(unit=3535,file='krome_bins_rates.dat')
       endif
      endif
      first_call=.False. 
      endif
      !end of initilazation routine
      do k = 1,ntau
        Pcon(k) = 0.1*Na/(R*T(k))*1E-6
        num_den_mol(k,:)=ppallmol(k,:)*Pcon(k)
        num_den_at(k,:)=ppallat(k,:)*Pcon(k)
      enddo                        

      do k=1,ntau
        do i=1,atom_counter
         num_den(k,index_at(i,2)) = num_den_at(k,index_at(i,1))
        enddo    
        do i=1,mol_counter
         num_den(k,index_mol(i,2)) = num_den_mol(k,index_mol(i,1))    
        enddo
        do i=1,mol2_counter
         num_den(k,index_mol2(i,2)) = num_den_mol(k,index_mol2(i,1))
        enddo
        do i=1,not_counter
         num_den(k,index_not(i)) = 1d-20
        enddo
        if (M_counter==1) then
         num_den(k,M_index) = ptot(k)*Pcon(k)
        endif
        mix_rat(k,:) = num_den(k,:)/(ptot(k)*Pcon(k))

      enddo
      !write(*,*) 'before krome num_den',num_den(1,:)
      !write header of full output file
      if (krome_debug.eq.1) then
        open(unit=13,file='krome_full_output.dat')
        write(13,'(A6,A9,A13,A16)',Advance = 'No') 'Layer ','Time [s] ',
     >'Timestep [s] ','Temperature [K] '
        do i=1,nsp
         write(13,'(A8,A9)',Advance = 'No') chem_spec(i),' [cm^-3] '
        enddo
         write(13,'(A,/)') ' '
      endif 

      !main loop

      do k=1,ntau
        istep = 1
        dt = dt_start
        time=0.0

        FLUX_RAD_eV=FLUX_RAD(k,nwreal:1:-1) !invert FLUX_RAD to align with energies as KROME requires it

        if (krome_photo_on.eq.1) then
         call krome_set_photoBinJ(FLUX_RAD_eV(:))
         call krome_photoBin_scale(krome_photo_scale)
        write(3535,'(I3)') k
         call krome_explore_flux(num_den(k,:),T(k),3535,FLUX_RAD_eV(k))        
         if (krome_debug.eq.1) then         
          if (k.eq.1) then !only print values for first layer
           write(4242,*) krome_get_photoBinJ()
           write(5656,*) krome_get_photobinE_mid()
           write(5959,*) krome_get_photobinE_delta()   
          endif
         endif 
        endif

        do
         call krome(num_den(k,:), T(k), dt) !call KROME

         num_den_cont(istep,:) = num_den(k,:)
         time_cont(istep) = time
         if (krome_debug.eq.1) then
          if(istep==1 .or. istep==2 .or. 
     & mod(istep,output_freq)==0) then
           write(13,'(I3,4(999E17.8e3))') k,time,dt,T(k),
     >     num_den(k,:)     
          end if
         end if
         dt = min(dt*dt_inc,dt_max)
       if (use_conv.eqv..True.) then
         if (dt.ge.dt_max) then !break loop if convergence is reached and timestep on max timestep length
            is_conv=.True.
            do j=1,nsp
 
            if (abs(num_den_cont(istep,j)-num_den_cont(istep-1,j))
     >                     /num_den_cont(istep,j)            
     >          .ge.conv_crit) then
                is_conv=.False.
                !write(*,*) k,j
                !write(*,*) time,istep,istep-1
                !write(*,*) num_den_cont(istep,j),num_den_cont(istep-1,j)
                !write(*,*) (abs(num_den_cont(istep,j)
     >          !           -num_den_cont(istep-1,j))
     >          !           /num_den_cont(istep,j))
                !write(*,*) abs(num_den_cont(istep,j)
     >          !           -num_den_cont(istep-1,j))
                exit
            endif

            enddo
            if (is_conv.eq..True.) then
              if (krome_debug.eq.1) then
                  write(13,'(I3,4(999E17.8e3))') k,time,dt,T(k),
     >            num_den(k,:)        
              endif
              exit   
            endif
          endif
         endif
         if(time>krome_tmax) then !break loop if maximum time is reached
           if (krome_debug.eq.1) then
            write(13,'(I3,4(999E17.8e3))') k,time,dt,T(k),
     >      num_den(k,:)        

           endif
          if (use_conv.eqv..True.) then
           write(*,*) 'Layer', k, 
     >      'did not converge within given time'
           write(*,*) 'relative change between timesteps'
           write(*,*) abs(num_den_cont(istep,j)-num_den_cont(istep-1,j))
     >                     /num_den_cont(istep,j) 
           write(*,*) 'with given conv. criteria', conv_crit 
           write(*,*) ' '    
            
          endif
          exit
         endif
         time = time + dt !increase time         
         istep = istep + 1 !increase timestep
        end do
      end do
      if (krome_debug.eq.1) then         
       close(3535)
      endif
      !final output
      if (krome_output.eq.1) then
        open(unit=77,file='krome_final_output.dat')
        write(77,'(A6,A9,A16)',Advance = 'No') 'Layer ','Time [s] '
     >,'Temperature [K] '
        do i=1,nsp
         write(77,'(A8,A9)',Advance = 'No') chem_spec(i),' [cm^-3] '
        enddo
         write(77,'(A,/)') ' '
        do k = 1,ntau
          write(77,'(I3,3(999E17.8e3))') k,time,T(k),
     >    num_den(k,:)
        enddo
      endif


C Returning the krome values to MARCS
      if (krome_return == 1) then
        do k=1,ntau
          do i=1,atom_counter
           num_den_at(k,index_at(i,1)) = num_den(k,index_at(i,2))
          enddo
          do i=1,mol_counter
           num_den_mol(k,index_mol(i,1)) = num_den(k,index_mol(i,2))
          enddo
          ppallmol(k,:) = num_den_mol(k,:)/Pcon(k)
          ppallat(k,:) = num_den_at(k,:)/Pcon(k)
        enddo
      endif


      if (krome_photo_on.eq.1) then
       if (krome_output.eq.1) then
        if ((ITSTOP.eq..True.).or.(it.eq.ITMAX)) then !write out all of fluxrad at the end of the iteration
         open(unit=7070,file='krome_flux_rad.dat')
         write(7070,'(A6,A17,A15,A24)') 'Layer ','Wavelength Index '
     >          ,'Wavelength [A] ','Fluxrad [eV/s/hz/cm2/sr]'       
         do k=1,ntau
          do j=1,nwreal
              write(7070,'(I3,2(999E17.8e3))') k,WLOS(j)
     >         ,FLUX_RAD(k,j)
         enddo
        enddo
        close(7070) 
        endif
       endif
      endif
      return      
      end
