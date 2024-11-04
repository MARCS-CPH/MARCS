
!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-09-23 13:10:57
  !  Changeset ed8d3da
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,r3,p1,p2

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    !O2
    !O2
    dn(idx_O2) = &
        -k(1)*n(idx_O2)*n(idx_O)*n(idx_M) &
        +2.d0*k(2)*n(idx_O3)*n(idx_O) &
        -k(3)*n(idx_O2) &
        +k(4)*n(idx_O3)

    !O
    !O
    dn(idx_O) = &
        -k(1)*n(idx_O2)*n(idx_O)*n(idx_M) &
        -k(2)*n(idx_O3)*n(idx_O) &
        +2.d0*k(3)*n(idx_O2) &
        +k(4)*n(idx_O3)

    !M
    !M
    dn(idx_M) = &
        -k(1)*n(idx_O2)*n(idx_O)*n(idx_M) &
        +k(1)*n(idx_O2)*n(idx_O)*n(idx_M)

    !O3
    !O3
    dn(idx_O3) = &
        +k(1)*n(idx_O2)*n(idx_O)*n(idx_M) &
        -k(2)*n(idx_O3)*n(idx_O) &
        -k(4)*n(idx_O3)

    !CR

    !CR
    dn(idx_CR) = 0.d0

    !g

    !g
    dn(idx_g) = 0.d0

    !Tgas

    !Tgas
    dn(idx_Tgas) = 0.d0

    !dummy

    !dummy
    dn(idx_dummy) = 0.d0

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(1)*n(idx_O)*n(idx_M)  &
          -k(3)
      pdj(2) =  &
          -k(1)*n(idx_O)*n(idx_M)  &
          +2.d0*k(3)
      pdj(3) =  &
          -k(1)*n(idx_O)*n(idx_M)  &
          +k(1)*n(idx_O)*n(idx_M)
      pdj(4) =  &
          +k(1)*n(idx_O)*n(idx_M)
    elseif(j==2) then
      pdj(1) =  &
          -k(1)*n(idx_O2)*n(idx_M)  &
          +2.d0*k(2)*n(idx_O3)
      pdj(2) =  &
          -k(1)*n(idx_O2)*n(idx_M)  &
          -k(2)*n(idx_O3)
      pdj(3) =  &
          -k(1)*n(idx_O2)*n(idx_M)  &
          +k(1)*n(idx_O2)*n(idx_M)
      pdj(4) =  &
          +k(1)*n(idx_O2)*n(idx_M)  &
          -k(2)*n(idx_O3)
    elseif(j==3) then
      pdj(1) =  &
          -k(1)*n(idx_O2)*n(idx_O)
      pdj(2) =  &
          -k(1)*n(idx_O2)*n(idx_O)
      pdj(3) =  &
          -k(1)*n(idx_O2)*n(idx_O)  &
          +k(1)*n(idx_O2)*n(idx_O)
      pdj(4) =  &
          +k(1)*n(idx_O2)*n(idx_O)
    elseif(j==4) then
      pdj(1) =  &
          +2.d0*k(2)*n(idx_O)  &
          +k(4)
      pdj(2) =  &
          -k(2)*n(idx_O)  &
          +k(4)
      pdj(4) =  &
          -k(2)*n(idx_O)  &
          -k(4)
    elseif(j==5) then
    elseif(j==6) then
    elseif(j==7) then

    elseif(j==8) then
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[O2_dot]/d[O2]
    pd(1,1) =  &
        -k(1)*n(idx_O)*n(idx_M)  &
        -k(3)

    !d[O_dot]/d[O2]
    pd(2,1) =  &
        -k(1)*n(idx_O)*n(idx_M)  &
        +2.d0*k(3)

    !d[M_dot]/d[O2]
    pd(3,1) =  &
        -k(1)*n(idx_O)*n(idx_M)  &
        +k(1)*n(idx_O)*n(idx_M)

    !d[O3_dot]/d[O2]
    pd(4,1) =  &
        +k(1)*n(idx_O)*n(idx_M)

    !d[O2_dot]/d[O]
    pd(1,2) =  &
        -k(1)*n(idx_O2)*n(idx_M)  &
        +2.d0*k(2)*n(idx_O3)

    !d[O_dot]/d[O]
    pd(2,2) =  &
        -k(1)*n(idx_O2)*n(idx_M)  &
        -k(2)*n(idx_O3)

    !d[M_dot]/d[O]
    pd(3,2) =  &
        -k(1)*n(idx_O2)*n(idx_M)  &
        +k(1)*n(idx_O2)*n(idx_M)

    !d[O3_dot]/d[O]
    pd(4,2) =  &
        +k(1)*n(idx_O2)*n(idx_M)  &
        -k(2)*n(idx_O3)

    !d[O2_dot]/d[M]
    pd(1,3) =  &
        -k(1)*n(idx_O2)*n(idx_O)

    !d[O_dot]/d[M]
    pd(2,3) =  &
        -k(1)*n(idx_O2)*n(idx_O)

    !d[M_dot]/d[M]
    pd(3,3) =  &
        -k(1)*n(idx_O2)*n(idx_O)  &
        +k(1)*n(idx_O2)*n(idx_O)

    !d[O3_dot]/d[M]
    pd(4,3) =  &
        +k(1)*n(idx_O2)*n(idx_O)

    !d[O2_dot]/d[O3]
    pd(1,4) =  &
        +2.d0*k(2)*n(idx_O)  &
        +k(4)

    !d[O_dot]/d[O3]
    pd(2,4) =  &
        -k(2)*n(idx_O)  &
        +k(4)

    !d[O3_dot]/d[O3]
    pd(4,4) =  &
        -k(2)*n(idx_O)  &
        -k(4)

  end subroutine jex

end module krome_ode
