MODULE EOS_F_CONTRIBUTIONS

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: F_IDEAL_GAS, F_HARD_SPHERE, F_CHAIN_TPT1, F_CHAIN_TPT_D,  &
            F_CHAIN_HU_LIU, F_CHAIN_HU_LIU_RC, F_SPT, F_DISP_PCSAFT,  &
            F_ASSOCIATION, F_ION_DIPOLE_TBH, F_ION_ION_PrimMSA,  &
            F_ION_ION_nonPrimMSA, F_LC_MayerSaupe, F_pert_theory

CONTAINS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_IDEAL_GAS ( fid )

  USE EOS_VARIABLES, ONLY: nc, ncomp, x, rho, PI, KBOL, NAv

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   ::  fid
  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL, DIMENSION(nc)                    :: rhoi
  !-----------------------------------------------------------------------------

  !h_Planck = 6.62606896E-34 ! Js
  DO i = 1, ncomp
     rhoi(i) = x(i) * rho
     ! debroglie(i) = h_Planck *1d10  &       ! in units Angstrom
     !               *SQRT( 1.0 / (2.0*PI *1.0 / NAv / 1000.0 * KBOL*T) )
     !               ! *SQRT( 1.0 / (2.0*PI *mm(i) /NAv/1000.0 * KBOL*T) )
     ! fid = fid + x(i) * ( LOG(rhoi(i)*debroglie(i)**3) - 1.0 )
     fid = fid + x(i) * ( LOG(rhoi(i)) - 1.0 )
  END DO

END SUBROUTINE F_IDEAL_GAS

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_HARD_SPHERE ( m_mean2, fhs )

  USE EOS_VARIABLES, ONLY: z0t, z1t, z2t, z3t, eta, rho

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       ::  m_mean2
  REAL, INTENT(IN OUT)                   ::  fhs
  !-----------------------------------------------------------------------------
  REAL                                   :: z0, z1, z2, z3, zms
  !-----------------------------------------------------------------------------

  rho = eta / z3t
  z0 = z0t * rho
  z1 = z1t * rho
  z2 = z2t * rho
  z3 = z3t * rho
  zms = 1.0 - z3

  fhs= m_mean2*( 3.0*z1*z2/zms + z2**3 /z3/zms/zms + (z2**3 /z3/z3-z0)*LOG(zms) )/z0


END SUBROUTINE F_HARD_SPHERE

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_CHAIN_TPT1 ( fhc )

  USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, z0t, z1t, z2t, z3t,  &
       rho, eta, dij_ab, gij

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      ::  fhc
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j
  REAL                                   :: z0, z1, z2, z3, zms
  !-----------------------------------------------------------------------------

  rho = eta / z3t
  z0 = z0t * rho
  z1 = z1t * rho
  z2 = z2t * rho
  z3 = z3t * rho
  zms = 1.0 - z3

  DO i = 1, ncomp
     DO j = 1, ncomp
        gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 / zms**3
     END DO
  END DO

  fhc = 0.0
  DO i = 1, ncomp
     fhc = fhc + x(i) *(1.0- mseg(i)) *LOG(gij(i,i))
  END DO

END SUBROUTINE F_CHAIN_TPT1


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_CHAIN_TPT_D ( fhc )

  USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, z0t, z1t, z2t, z3t, rho, eta,  &
       mseg, dij_ab, gij

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      ::  fhc
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j
  REAL, DIMENSION(nc)                    :: gij_hd
  REAL                                   :: z0, z1, z2, z3, zms
  !-----------------------------------------------------------------------------

  rho = eta / z3t
  z0 = z0t * rho
  z1 = z1t * rho
  z2 = z2t * rho
  z3 = z3t * rho
  zms = 1.0 - z3

  DO i = 1, ncomp
     DO j = 1, ncomp
        gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 / zms**3
     END DO
  END DO

  DO i = 1, ncomp
     gij_hd(i) = 1.0/(2.0*zms) + 3.0*dij_ab(i,i)*z2 / zms**2
  END DO

  fhc = 0.0
  DO i = 1, ncomp
     IF ( mseg(i) >= 2.0 ) THEN
        fhc = fhc - x(i) * ( mseg(i)/2.0 * LOG( gij(i,i) ) + ( mseg(i)/2.0 - 1.0 ) * LOG( gij_hd(i)) )
     ELSE
        fhc = fhc + x(i) * ( 1.0 - mseg(i) ) * LOG( gij(i,i) )
     END IF
  END DO

END SUBROUTINE F_CHAIN_TPT_D


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_CHAIN_HU_LIU ( fhc )

  USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, eta

  ! This subroutine calculates the hard chain contribution of the TPT-Liu-Hu Eos.
  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      ::  fhc
  !-----------------------------------------------------------------------------
  REAL                                   :: a2, b2, c2, a3, b3, c3
  REAL                                   :: a20, b20, c20, a30, b30, c30
  REAL                                   :: sum1, sum2, am, bm, cm
  REAL                                   :: zms
  !-----------------------------------------------------------------------------

  zms = 1.0 - eta

  sum1 = SUM( x(1:ncomp)*(mseg(1:ncomp)-1.0) )
  sum2 = SUM( x(1:ncomp)/mseg(1:ncomp)*(mseg(1:ncomp)-1.0)*(mseg(1:ncomp)-2.0) )

  a2  =  0.45696
  a3  = -0.74745
  b2  =  2.10386
  b3  =  3.49695
  c2  =  1.75503
  c3  =  4.83207
  a20 = - a2 + b2 - 3.0*c2
  b20 = - a2 - b2 + c2
  c20 =   c2
  a30 = - a3 + b3 - 3.0*c3
  b30 = - a3 - b3 + c3
  c30 =   c3
  am  = (3.0 + a20) * sum1 + a30 * sum2
  bm  = (1.0 + b20) * sum1 + b30 * sum2
  cm  = (1.0 + c20) * sum1 + c30 * sum2

  fhc = - ( (am*eta - bm) / (2.0*zms) + bm/2.0/zms**2 - cm *LOG(ZMS) )


END SUBROUTINE F_CHAIN_HU_LIU

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_CHAIN_HU_LIU_RC ( fhs, fhc )

  USE EOS_VARIABLES, ONLY: mseg, chiR, eta

  ! This subroutine calculates the hard chain contribution of the TPT-Liu-Hu Eos.
  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       ::  fhs
  REAL, INTENT(OUT)                      ::  fhc
  !-----------------------------------------------------------------------------
  REAL                                   :: a2, b2, c2, a3, b3, c3
  REAL                                   :: para1,para2,para3,para4
  REAL                                   :: aLH,bLH,cLH
  !-----------------------------------------------------------------------------

  ! This routine is only for pure components

  a2 = 0.45696
  b2 = 2.10386
  c2 = 1.75503

  para1 = -0.74745
  para2 = 0.299154629727814
  para3 = 1.087271036653154
  para4 = -0.708979110326831
  a3 = para1 + para2*chiR(1) + para3*chiR(1)**2 + para4*chiR(1)**3
  b3 = 3.49695 - (3.49695 + 0.317719074806190)*chiR(1)
  c3 = 4.83207 - (4.83207 - 3.480163780334421)*chiR(1)

  aLH = mseg(1)*(1.0 + ((mseg(1)-1.0)/mseg(1))*a2 + ((mseg(1)-1.0)/mseg(1))*((mseg(1)-2.0)/mseg(1))*a3 )
  bLH = mseg(1)*(1.0 + ((mseg(1)-1.0)/mseg(1))*b2 + ((mseg(1)-1.0)/mseg(1))*((mseg(1)-2.0)/mseg(1))*b3 )
  cLH = mseg(1)*(1.0 + ((mseg(1)-1.0)/mseg(1))*c2 + ((mseg(1)-1.0)/mseg(1))*((mseg(1)-2.0)/mseg(1))*c3 )

  fhc = ((3.0 + aLH - bLH + 3.0*cLH)*eta - (1.0 + aLH + bLH - cLH)) / (2.0*(1.0-eta)) + &
       (1.0 + aLH + bLH - cLH) / ( 2.0*(1.0-eta)**2 ) + (cLH - 1.0)*LOG(1.0-eta)

  fhc = fhc - fhs

END SUBROUTINE F_CHAIN_HU_LIU_RC


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_SPT
!
! This routine calculates the residual Helmholtz energy contribution
! obtained from Scaled Particle Theory (SPT).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_SPT ( fhs, fhc )

  USE EOS_VARIABLES, ONLY: nc, ncomp, PI, mseg, chiR, eta


  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                   ::  fhc
  REAL, INTENT(IN)                    ::  fhs
  !-----------------------------------------------------------------------------
  real, dimension(nc,nc)              :: C0, C1, C2

  INTEGER                             :: i,j
  REAL                                :: msegav,aexv_1,bexv_1,chiRav,aexv_2,aexv_3
  REAL                                :: bexv_2,Vm,B2,a_spt,Sm,psi_spt
  REAL, DIMENSION(nc,nc)              :: B2iso
  LOGICAL                             :: Vex_williamson, Vex_correlation
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Define which model to use for calculating the excluded volume
  !-----------------------------------------------------------------------------
  Vex_williamson = .False.
  Vex_correlation = .True.

  ! ---------------------------------------------------------------------
  ! Calculate isotropic 2nd virial coefficient B2iso
  ! ---------------------------------------------------------------------

  IF ( Vex_williamson ) THEN
     DO i = 1, ncomp
        DO j = 1, ncomp
           msegav = ( mseg(i) + mseg(j) ) / 2.0
           C0(i,j) = (pi/6.0)*(11.0*msegav - 3.0)
           C1(i,j) = (pi/6.0)*3.53390 *(mseg(i)-1.0) * (mseg(j)-1.0)
           C2(i,j) = 0.0
           B2iso(i,j) = C0(i,j)/2.0 + (pi/8.0)*C1(i,j) + (1.0/3.0)*C2(i,j)
        END DO
     END DO
  ELSEIF ( Vex_correlation ) THEN
     aexv_1 = 4.63
     bexv_1 = 0.305
     DO i = 1, ncomp
        DO j = 1, ncomp
           msegav = ( mseg(i) + mseg(j) ) / 2.0
           chiRav = ( chiR(i) + chiR(j) ) / 2.0
           aexv_2 = -4.70 + 7.84/msegav
           aexv_3 = 1.31 - 6.18/msegav
           bexv_2 = -0.171 + 3.32/msegav

           C0(i,j) = (pi/6.0)* (  (11.0*msegav - 3.0) + (mseg(i) - 1.0)*  &
                (mseg(j) -1.0)*( aexv_1*(1.0-chiRav) + aexv_2*(1.0-chiRav)**2  &
                + aexv_3*(1.0-chiRav)**3 ) )
           C1(i,j) = (pi/6.0)*3.53390*(mseg(i) - 1.0)*(mseg(j) -1.0)*chiRav**2
           C2(i,j) = (pi/6.0)*(mseg(i) - 1.0)*(mseg(j) -1.0)*( bexv_1*(1.0-chiRav)  &
                + bexv_2*(1.0-chiRav)**2 )
           B2iso(i,j) = C0(i,j)/2.0 + (pi/8.0)*C1(i,j) + (1.0/3.0)*C2(i,j)
        END DO
     END DO
  ELSE
     write(*,*) 'no model for excluded volume defined'
     Stop
  END IF

  !write(*,*)'f_spt'
  !write(*,*)'C0, C1, C2 :',C0(1,1),C1(1,1),C2(1,1)

  ! This only holds for pure components!!!!!
  Vm = (PI/6.0) * mseg(1)
  B2 = B2iso(1,1) / Vm;
  a_spt = (B2 - 1.0) / 3.0
  Sm = PI*mseg(1)
  psi_spt = (Sm/(9.0*Vm)) * (3.0*a_spt - (Sm/(4.0*Vm))*(1.0 - (mseg(1)-1.0)*(1.0/4.0) ) )
  fhc = (psi_spt -1.0)*LOG(1.0-eta) + 3.0*a_spt*eta/(1.0-eta) + psi_spt*eta/( (1.0-eta)**2 )

  fhc = fhc - fhs

END SUBROUTINE F_SPT

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_DISP_PCSAFT ( fdsp )

  USE EOS_VARIABLES, ONLY: PI, rho, eta, z0t, apar, bpar, order1, order2

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdsp
  !-----------------------------------------------------------------------------
  INTEGER                                :: m
  REAL                                   :: I1, I2, c1_con, z3, zms, m_mean
  !-----------------------------------------------------------------------------

  z3  = eta
  zms = 1.0 - eta
  m_mean = z0t / ( PI / 6.0 )

  I1   = 0.0
  I2   = 0.0
  DO m = 0, 6
     I1 = I1 + apar(m) * z3**m
     I2 = I2 + bpar(m) * z3**m
  END DO

  c1_con= 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3**2 )/zms**4  &
       + (1.0-m_mean)*( 20.0*z3-27.0*z3**2 +12.0*z3**3 -2.0*z3**4 )/(zms*(2.0-z3))**2 )

  fdsp  = -2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2

END SUBROUTINE F_DISP_PCSAFT

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!!$SUBROUTINE F_DISP_CKSAFT ( fdsp )
!!$
!!$  USE EOS_VARIABLES, ONLY: nc, ncomp, PI, TAU, t, rho, eta, x, z0t, mseg, vij, uij, parame, um
!!$  USE EOS_CONSTANTS, ONLY: DNM
!!$  IMPLICIT NONE
!!$
!!$  !-----------------------------------------------------------------------------
!!$  REAL, INTENT(IN OUT)                   :: fdsp
!!$  !-----------------------------------------------------------------------------
!!$  INTEGER                                :: i, j, n, m
!!$  REAL                                   :: zmr, nmr, m_mean
!!$  !-----------------------------------------------------------------------------
!!$
!!$  m_mean = z0t / ( PI / 6.0 )
!!$
!!$  DO i = 1, ncomp
!!$     DO j = 1, ncomp
!!$        vij(i,j)=(0.5*((parame(i,2)*(1.0-0.12 *EXP(-3.0*parame(i,3)/t))**3 )**(1.0/3.0)  &
!!$             +(parame(j,2)*(1.0-0.12 *EXP(-3.0*parame(j,3)/t))**3 )**(1.0/3.0)))**3
!!$     END DO
!!$  END DO
!!$  zmr = 0.0
!!$  nmr = 0.0
!!$  DO i = 1, ncomp
!!$     DO j = 1, ncomp
!!$        zmr = zmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)*uij(i,j)
!!$        nmr = nmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)
!!$     END DO
!!$  END DO
!!$  um = zmr / nmr
!!$  fdsp  = 0.0
!!$  DO n = 1, 4
!!$     DO m = 1, 9
!!$        fdsp =  fdsp + DNM(n,m) * (um/t)**n *(eta/TAU)**m
!!$     END DO
!!$  END DO
!!$  fdsp = m_mean * fdsp
!!$
!!$
!!$END SUBROUTINE F_DISP_CKSAFT

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_ASSOCIATION ( eps_kij, k_kij, fhb )

  USE EOS_VARIABLES, ONLY: nc, nsite, ncomp, t, z0t, z1t, z2t, z3t, rho, eta, x,  &
       parame, sig_ij, dij_ab, gij, nhb_typ, mx, nhb_no
  USE utilities

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: eps_kij, k_kij
  REAL, INTENT(IN OUT)                   :: fhb
  !-----------------------------------------------------------------------------
  LOGICAL                                :: assoc
  INTEGER                                :: i, j, k, l, no, ass_cnt, max_eval
  REAL, DIMENSION(nc,nc)                 :: kap_hb
  REAL, DIMENSION(nc,nc,nsite,nsite)     :: eps_hb
  REAL, DIMENSION(nc,nc,nsite,nsite)     :: delta
  REAL, DIMENSION(nc,nsite)              :: mx_itr
  REAL                                   :: err_sum, sum0, attenu, tol, ass_s1, ass_s2
  REAL                                   :: z0, z1, z2, z3, zms
  !-----------------------------------------------------------------------------

  assoc = .false.
  DO i = 1,ncomp
     IF (NINT(parame(i,12)) /= 0) assoc = .true.
  END DO
  IF (assoc) THEN

     rho = eta / z3t
     z0 = z0t * rho
     z1 = z1t * rho
     z2 = z2t * rho
     z3 = z3t * rho
     zms = 1.0 - z3

     DO i = 1, ncomp
        DO j = 1, ncomp
           gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 / zms**3
        END DO
     END DO


     DO  i = 1, ncomp
        IF ( NINT(parame(i,12)) /= 0 ) THEN
           nhb_typ(i) = NINT( parame(i,12) )
           kap_hb(i,i) = parame(i,13)
           no = 0
           DO k = 1,nhb_typ(i)
              DO l = 1,nhb_typ(i)
                 eps_hb(i,i,k,l) = parame(i,(14+no))
                 no = no + 1
              END DO
           END DO
           DO k = 1,nhb_typ(i)
              nhb_no(i,k) = parame(i,(14+no))
              no = no + 1
           END DO
        ELSE
           nhb_typ(i) = 0
           kap_hb(i,i)= 0.0
           DO k = 1,nsite
              DO l = 1,nsite
                 eps_hb(i,i,k,l) = 0.0
              END DO
           END DO
        END IF
     END DO

     DO i = 1,ncomp
        DO j = 1,ncomp
           IF ( i /= j .AND. (nhb_typ(i) /= 0.AND.nhb_typ(j) /= 0) ) THEN
              ! kap_hb(i,j)= (kap_hb(i,i)+kap_hb(j,j))/2.0
              ! kap_hb(i,j)= ( ( kap_hb(i,i)**(1.0/3.0) + kap_hb(j,j)**(1.0/3.0) )/2.0 )**3
              kap_hb(i,j) = (kap_hb(i,i)*kap_hb(j,j))**0.5  &
                   *((parame(i,2)*parame(j,2))**3 )**0.5  &
                   / (0.5*(parame(i,2)+parame(j,2)))**3
              kap_hb(i,j)= kap_hb(i,j)*(1.0-k_kij)
              DO k = 1,nhb_typ(i)
                 DO l = 1,nhb_typ(j)
                    IF ( k /= l .AND. nhb_typ(i) >= 2 .AND. nhb_typ(j) >= 2 ) THEN
                       eps_hb(i,j,k,l) = (eps_hb(i,i,k,l)+eps_hb(j,j,l,k))/2.0
                       ! eps_hb(i,j,k,l) = (eps_hb(i,i,k,l)*eps_hb(j,j,l,k))**0.5
                       eps_hb(i,j,k,l) = eps_hb(i,j,k,l)*(1.0-eps_kij)
                    ELSE IF ( nhb_typ(i) == 1 .AND. l > k ) THEN
                       eps_hb(i,j,k,l) = (eps_hb(i,i,k,k)+eps_hb(j,j,l,k))/2.0
                       eps_hb(j,i,l,k) = (eps_hb(i,i,k,k)+eps_hb(j,j,l,k))/2.0
                       eps_hb(i,j,k,l) = eps_hb(i,j,k,l)*(1.0-eps_kij)
                       eps_hb(j,i,l,k) = eps_hb(j,i,l,k)*(1.0-eps_kij)
                    END IF
                 END DO
              END DO
           END IF
        END DO
     END DO

     !--------------------------------------------------------------------------
     ! setting the self-association to zero for ionic compounds
     !--------------------------------------------------------------------------
     DO i = 1,ncomp
        IF ( parame(i,10) /= 0) kap_hb(i,i)=0.0
        DO j = 1,ncomp
           IF ( parame(i,10) /= 0 .AND. parame(j,10) /= 0 ) kap_hb(i,j) = 0.0
        END DO
     END DO
     ! kap_hb(1,2)=0.050
     ! kap_hb(2,1)=0.050
     ! eps_hb(2,1,1,1)=465.0
     ! eps_hb(1,2,1,1)=465.0
     ! nhb_typ(1) = 1
     ! nhb_typ(2) = 1
     ! nhb_no(1,1)= 1.0
     ! nhb_no(2,1)= 1.0


     DO i = 1, ncomp
        DO k = 1, nhb_typ(i)
           DO j = 1, ncomp
              DO l = 1, nhb_typ(j)
                 delta(i,j,k,l) =gij(i,j) *kap_hb(i,j) *(EXP(eps_hb(i,j,k,l)/t) - 1.0) *sig_ij(i,j)**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !              IF ((i+j).EQ.3) delta(i,j,k,l)=94.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              END DO
           END DO
           IF ( mx(i,k) == 0.0 ) mx(i,k) = 1.0
        END DO
     END DO

     !--------------------------------------------------------------------------
     ! constants for Iteration
     !--------------------------------------------------------------------------
     attenu = 0.7
     tol = 1.E-10
     IF (eta < 0.2) tol = 1.E-11
     IF (eta < 0.01) tol = 1.E-14
     max_eval = 200
     err_sum = 2.0 * tol

     !--------------------------------------------------------------------------
     ! Iterate over all components and all sites
     !--------------------------------------------------------------------------
     ass_cnt = 0
     DO WHILE ( err_sum > tol .AND. ass_cnt <= max_eval )

        ass_cnt = ass_cnt + 1

        DO i = 1, ncomp
           DO k = 1, nhb_typ(i)
              sum0 = 0.0
              DO j = 1, ncomp
                 DO l = 1, nhb_typ(j)
                    sum0 = sum0 +  x(j) * mx(j,l) * nhb_no(j,l) * delta(i,j,k,l)
                 END DO
              END DO
              mx_itr(i,k) = 1.0 / (1.0 + sum0 * rho)
           END DO
        END DO

        err_sum = 0.0
        DO i = 1, ncomp
           DO k = 1, nhb_typ(i)
              err_sum = err_sum + ABS( mx_itr(i,k) - mx(i,k) )    ! / ABS(mx_itr(i,k))
              mx(i,k) = mx_itr(i,k) * attenu + mx(i,k) * (1.0 - attenu)
              IF ( mx(i,k) <= 0.0 ) mx(i,k) = 1.E-50
              IF ( mx(i,k) > 1.0 )  mx(i,k) = 1.0
           END DO
        END DO

     END DO

     IF ( err_sum /= err_sum ) write (*,*) 'F_EOS: association NaN',ass_cnt, rho, sum0
     IF ( err_sum /= err_sum ) read (*,*)
     IF ( ass_cnt >= max_eval .AND. err_sum > SQRT(tol) ) THEN
        WRITE(*,'(a,2G15.7)') 'F_EOS: Max_eval violated (mx) Err_Sum = ',err_sum,tol
     END IF

     DO i = 1, ncomp
        ass_s1 = 0.0
        ass_s2 = 0.0
        DO k = 1, nhb_typ(i)
           ass_s1 = ass_s1 + nhb_no(i,k) * ( 1.0 - mx(i,k) )
           ass_s2 = ass_s2 + nhb_no(i,k) * LOG(mx(i,k))
        END DO
        fhb = fhb + x(i) * ( ass_s2 + ass_s1 / 2.0 )
     END DO

  END IF

END SUBROUTINE F_ASSOCIATION


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_ION_DIPOLE_TBH ( fhend )

  USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, t, rho, eta, x, z0t, parame, uij, sig_ij

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fhend
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, dipole, ions
  REAL                                   :: m_mean
  REAL                                   :: fh32, fh2, fh52, fh3
  REAL                                   :: e_elem, eps_cc0, rho_sol, dielec
  REAL                                   :: polabil, ydd, kappa, x_dipol, x_ions
  REAL, DIMENSION(nc)                    :: my2dd, z_ii, e_cd, x_dd, x_ii
  REAL                                   :: sig_c, sig_d, sig_cd, r_s
  REAL                                   :: I0cc, I1cc, I2cc, Icd, Idd
  REAL                                   :: Iccc, Iccd, Icdd, Iddd
  !-----------------------------------------------------------------------------

  m_mean = z0t / ( PI / 6.0 )

  !-----------------------------------------------------------------------------
  ! Dieletric Constant of Water
  !-----------------------------------------------------------------------------
  e_elem   = 1.602189246E-19   ! in Unit [Coulomb]
  eps_cc0  = 8.854187818E-22   ! in Unit [Coulomb**2 /(J*Angstrom)]
  ! Correlation of M. Uematsu and E. U. Frank
  ! (Static Dieletric Constant of Water and Steam)
  ! valid range of conditions 273,15 K <=T<= 823,15 K
  ! and density <= 1150 kg/m3 (i.e. 0 <= p <= 500 MPa)
  rho_sol = rho * 18.015 * 1.E27/ NAv
  rho_sol = rho_sol/1000.0
  dielec = 1.0+(7.62571/(t/293.15))*rho_sol +(2.44E2/(t/293.15)-1.41E2  &
       + 2.78E1*(t/293.15))*rho_sol**2   &
       + (-9.63E1/(t/293.15)+4.18E1*(t/293.15)  &
       - 1.02E1*(t/293.15)**2 )*rho_sol**3  +(-4.52E1/(t/293.15)**2   &
       + 8.46E1/(t/293.15)-3.59E1)*rho_sol**4


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dielec = 1.0

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !-----------------------------------------------------------------------------
  ! Dipole-Ion Term
  !-----------------------------------------------------------------------------
  dipole = 0
  ions   = 0
  fhend   = 0.0
  DO i = 1, ncomp
     IF ( parame(i,6) /= 0.0 .AND. uij(i,i) /= 0.0 .AND. x(i) > 0.0 ) THEN
        my2dd(i) = (parame(i,6))**2  *1.E-49 / (uij(i,i)*KBOL*sig_ij(i,i)**3 *1.E-30)
        dipole = 1
     ELSE
        my2dd(i) = 0.0
     END IF

     z_ii(i) = parame(i,10)
     IF ( z_ii(i) /= 0.0 .AND. uij(i,i) /= 0.0 .AND. x(i) > 0.0 ) THEN
        e_cd(i) = ( parame(i,10)*e_elem* 1.E5 / SQRT(1.11265005) )**2   &
             / ( uij(i,i)*kbol*sig_ij(i,i)*1.E-10 )
        ions = 1
     ELSE
        e_cd(i) = 0.0
     END IF
  END DO


  IF ( dipole == 1 .AND. ions == 1 ) THEN

     ydd     = 0.0
     kappa   = 0.0
     x_dipol = 0.0
     x_ions  = 0.0
     polabil = 0.0
     DO i = 1, ncomp
        ydd = ydd + x(i)*(parame(i,6))**2  *1.E-49/ (kbol*t*1.E-30)
        kappa = kappa + x(i)  &
             *(parame(i,10)*e_elem* 1.E5/SQRT(1.11265005))**2  /(KBOL*t*1.E-10)
        IF (parame(i,10) /= 0.0) THEN
           x_ions = x_ions + x(i)
        ELSE
           polabil = polabil + 4.0*PI*x(i)*rho*1.4573 *1.E-30  &
                / (sig_ij(3,3)**3 *1.E-30)
        END IF
        IF (parame(i,6) /= 0.0) x_dipol= x_dipol+ x(i)
     END DO
     ydd   = ydd * 4.0/9.0 * PI * rho
     kappa = SQRT( 4.0 * PI * rho * kappa )

     fh2 = 0.0
     sig_c = 0.0
     sig_d = 0.0
     DO i=1,ncomp
        x_ii(i) = 0.0
        x_dd(i) = 0.0
        IF(parame(i,10) /= 0.0 .AND. x_ions /= 0.0) x_ii(i) = x(i)/x_ions
        IF(parame(i,6) /= 0.0 .AND. x_dipol /= 0.0) x_dd(i) = x(i)/x_dipol
        sig_c  = sig_c + x_ii(i)*parame(i,2)
        sig_d  = sig_d + x_dd(i)*parame(i,2)
     END DO
     sig_cd = 0.5 * (sig_c + sig_d)

     r_s = 0.0
     ! DO i=1,ncomp
     !   r_s=r_s + rho * x(i) * dhs(i)**3
     ! END DO
     r_s = eta*6.0 / PI / m_mean

     I0cc =  - (1.0 + 0.97743 * r_s + 0.05257*r_s*r_s)  &
          /(1.0 + 1.43613 * r_s + 0.41580*r_s*r_s)
     ! I1cc =  - (10.0 - 2.0*z3 + z3*z3) /20.0/(1.0 + 2.0*z3)
     I1cc =  - (10.0 - 2.0*r_s*pi/6.0 + r_s*r_s*pi/6.0*pi/6.0)  &
          /20.0/(1.0 + 2.0*r_s*pi/6.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! I2cc =  + (z3-4.0)*(z3*z3+2.0) /24.0/(1.0+2.0*z3)
     !  relation of Stell and Lebowitz
     I2cc = -0.33331+0.7418*r_s - 1.2047*r_s*r_s  &
          + 1.6139*r_s**3 - 1.5487*r_s**4 + 0.6626*r_s**5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Icd =  (1.0 + 0.79576 *r_s + 0.104556 *r_s*r_s)  &
          /(1.0 + 0.486704*r_s - 0.0222903*r_s*r_s)
     Idd =  (1.0 + 0.18158*r_s - 0.11467*r_s*r_s)  &
          /3.0/(1.0 - 0.49303*r_s + 0.06293*r_s*r_s)
     Iccc=  3.0*(1.0 - 1.05560*r_s + 0.26591*r_s*r_s)  &
          /2.0/(1.0 + 0.53892*r_s - 0.94236*r_s*r_s)
     Iccd=  11.0*(1.0 + 2.25642 *r_s + 0.05679  *r_s*r_s)  &
          /6.0/(1.0 + 2.64178 *r_s + 0.79783  *r_s*r_s)
     Icdd=  0.94685*(1.0 + 2.97323 *r_s + 3.11931  *r_s*r_s)  &
          /(1.0 + 2.70186 *r_s + 1.22989  *r_s*r_s)
     Iddd=  5.0*(1.0   + 1.12754*r_s + 0.56192*r_s*r_s)  &
          /24.0/(1.0 - 0.05495*r_s + 0.13332*r_s*r_s)

     IF ( sig_c <= 0.0 ) WRITE (*,*) 'error in Henderson ion term'

     fh32= - kappa**3 /(12.0*pi*rho)
     fh2 = - 3.0* kappa**2  * ydd*Icd /(8.0*pi*rho) / sig_cd  &
          - kappa**4 *sig_c/(16.0*pi*rho)*I0cc
     IF (sig_d /= 0.0) fh2 = fh2 - 27.0* ydd * ydd*Idd  &
          /(8.0*pi*rho) / sig_d**3
     fh52= (3.0*kappa**3  * ydd +  kappa**5 *sig_c**2 *I1cc)  &
          /(8.0*pi*rho)
     fh3 =  - kappa**6 * sig_c**3 /(8.0*pi*rho) *(I2cc-Iccc/6.0)  &
          + kappa**4  * ydd *sig_c/(16.0*pi*rho)  &
          *( (6.0+5.0/3.0*sig_d/sig_c)*I0cc + 3.0*sig_d/sig_c*Iccd )  &
          + 3.0*kappa**2  * ydd*ydd /(8.0*pi*rho) / sig_cd  &
          *( (2.0-3.21555*sig_d/sig_cd)*Icd +3.0*sig_d/sig_cd*Icdd )
     IF (sig_d /= 0.0) fh3 = fh3 + 27.0*ydd**3   &
          /(16.0*pi*rho)/sig_d**3 *Iddd

     fhend = ( fh32 + (fh32*fh32*fh3-2.0*fh32*fh2*fh52+fh2**3 )  &
          /(fh2*fh2-fh32*fh52)  )  &
          /   ( 1.0 + (fh32*fh3-fh2*fh52) /(fh2*fh2-fh32*fh52)  &
          - (fh2*fh3-fh52*fh52) /(fh2*fh2-fh32*fh52)  )
     !----------
     ! fH32= - kappa**3 /(12.0*PI*rho)
     ! fH2 = - 3.0* kappa**2  * ydd*Icd /(8.0*PI*rho) / sig_cd
     ! fH52= (3.0*kappa**3  * ydd)/(8.0*PI*rho)
     ! fH3 =  + kappa**4  * ydd *sig_c/(16.0*PI*rho)  &
     !         *( (6.0+5.0/3.0*sig_d/sig_c)*0.0*I0cc + 3.0*sig_d/sig_c*Iccd)  &
     !             + 3.0*kappa**2  * ydd*ydd /(8.0*PI*rho) / sig_cd  &
     !         *( (2.0-3.215550*sig_d/sig_cd)*Icd +3.0*sig_d/sig_cd*Icdd )

     ! fHcd = (  + (fH32*fH32*fH3-2.0*fH32*fH2*fH52+fH2**3 )  &
     !                                                /(fH2*fH2-fH32*fH52)  )  &
     !          /   ( 1.0 + (fH32*fH3-fH2*fH52) /(fH2*fH2-fH32*fH52)  &
     !                     - (fH2*fH3-fH52*fH52) /(fH2*fH2-fH32*fH52)  )

  END IF

END SUBROUTINE F_ION_DIPOLE_TBH


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_ION_ION_PrimMSA ( fcc )

  USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, t, rho, x, parame, mx

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fcc
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, cc_it, ions
  REAL                                   :: e_elem, eps_cc0, rho_sol, dielec
  REAL                                   :: x_ions
  REAL                                   :: cc_sig1, cc_sig2, cc_sig3
  REAL, DIMENSION(nc)                    :: z_ii, x_ii, sigm_i, my2dd
  REAL                                   :: alpha_2, kappa, ii_par
  REAL                                   :: cc_omeg, p_n, q2_i, cc_q2, cc_gam
  REAL                                   :: cc_error(2), cc_delt
  REAL                                   :: rhs, lambda, lam_s
  !-----------------------------------------------------------------------------

  !----------------Dieletric Constant of Water--------------------------
  e_elem   = 1.602189246E-19   ! in Unit [Coulomb]
  eps_cc0  = 8.854187818E-22   ! in Unit [Coulomb**2 /(J*Angstrom)]
  ! Correlation of M. Uematsu and E. U. Frank
  ! (Static Dieletric Constant of Water and Steam)
  ! valid range of conditions 273,15 K <=T<= 823,15 K
  ! and density <= 1150 kg/m3 (i.e. 0 <= p <= 500 MPa)
  rho_sol = rho * 18.015 * 1.E27/ NAv
  rho_sol = rho_sol/1000.0
  dielec = 1.0+(7.62571/(t/293.15))*rho_sol +(2.44E2/(t/293.15)-1.41E2  &
       +2.78E1*(t/293.15))*rho_sol**2   &
       +(-9.63E1/(t/293.15)+4.18E1*(t/293.15)  &
       -1.02E1*(t/293.15)**2 )*rho_sol**3  +(-4.52E1/(t/293.15)**2   &
       +8.46E1/(t/293.15)-3.59E1)*rho_sol**4


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dielec = 1.0

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !----------------Ion-Ion: primitive MSA -------------------------------
  ! the (dipole moment)**2 [my**2] corresponds to an attraction from
  ! point charges of [ SUM(xi * zi**2 * e_elem**2) * 3 * di**2 ]

  ! parame(ion,6))**2  * 1.E-49 / (kbol*T)
  !                      = (e_elem* 1.E5/SQRT(1.112650050))**2
  !                         *x(i)*zi**2  *3.0*sig_ij(1,1)**2  *1.E-20

  ! parame(ion,6))**2  = (e_elem* 1.E5/SQRT(1.112650050))**2  /1.E-49
  !                         *x(i)*zi**2  *3.0*sig_ij(i,i)**2  *1.E-20

  ! with the units
  ! my**2       [=] D**2 = 1.E-49 J*m3
  ! e_elem **2  [=] C**2 = 1.E5 / SQRT(1.112650050) J*m


  ions   = 0
  x_ions = 0.0
  fcc = 0.0
  DO i = 1, ncomp
     z_ii(i)   = parame(i,10)
     IF (z_ii(i) /= 0.0) THEN
        sigm_i(i) = parame(i,2)
     ELSE
        sigm_i(i) = 0.0
     END IF
     IF (z_ii(i) /= 0.0) ions = 1
     IF (z_ii(i) /= 0.0) x_ions = x_ions + x(i)
  END DO

  IF (ions == 1 .AND. x_ions > 0.0) THEN

     cc_sig1 = 0.0
     cc_sig2 = 0.0
     cc_sig3 = 0.0
     DO i=1,ncomp
        IF (z_ii(i) /= 0.0) THEN
           x_ii(i) = x(i)/x_ions
        ELSE
           x_ii(i) =0.0
        END IF
        cc_sig1 = cc_sig1 +x_ii(i)*sigm_i(i)
        cc_sig2 = cc_sig2 +x_ii(i)*sigm_i(i)**2
        cc_sig3 = cc_sig3 +x_ii(i)*sigm_i(i)**3
     END DO


     ! alpha_2 = 4.0*PI*e_elem**2 /eps_cc0/dielec/kbol/T
     alpha_2 = e_elem**2 /eps_cc0 / dielec / KBOL/t
     kappa   = 0.0
     DO i = 1, ncomp
        kappa = kappa + x(i)*z_ii(i)*z_ii(i)*mx(i,1)
     END DO
     kappa = SQRT( rho * alpha_2 * kappa )
     ii_par= kappa * cc_sig1

     ! Temporaer: nach der Arbeit von Krienke verifiziert
     ! noch nicht fuer Mischungen mit unterschiedl. Ladung erweitert
     ! ii_par = SQRT( e_elem**2 /eps_cc0/dielec/kbol/T  &
     !          *rho*(x(1)*Z_ii(1)**2 + x(2)*Z_ii(2)**2 ) )*cc_sig1


     cc_gam = kappa/2.0

     ! noch offen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     cc_delt = 0.0
     DO i = 1, ncomp
        cc_delt = cc_delt + x(i)*mx(i,1)*rho*sigm_i(i)**3
     END DO
     cc_delt= 1.0 - PI/6.0*cc_delt

     cc_it = 0
13   CONTINUE
     j = 0
     cc_it = cc_it + 1
131  CONTINUE
     j = j + 1
     cc_omeg = 0.0
     DO i = 1, ncomp
        cc_omeg = cc_omeg +x(i)*mx(i,1)*sigm_i(i)**3 /(1.0+cc_gam*sigm_i(i))
     END DO
     cc_omeg = 1.0 + PI/2.0 / cc_delt * rho * cc_omeg
     p_n = 0.0
     DO i = 1, ncomp
        p_n = p_n + x(i)*mx(i,1)*rho / cc_omeg*sigm_i(i)*z_ii(i) / (1.0+cc_gam*sigm_i(i))
     END DO
     q2_i = 0.0
     cc_q2= 0.0
     DO i = 1, ncomp
        q2_i = q2_i + rho*x(i)*mx(i,1)*( (z_ii(i)-pi/2.0/cc_delt*sigm_i(i)**2 *p_n)  &
             /(1.0+cc_gam*sigm_i(i)) )**2
        cc_q2 = cc_q2 + x(i)*mx(i,1)*rho*z_ii(i)**2 / (1.0+cc_gam*sigm_i(i))
     END DO
     q2_i = q2_i*alpha_2 / 4.0

     cc_error(j) = cc_gam - SQRT(q2_i)
     IF (j == 1) cc_gam = cc_gam*1.000001
     IF (j == 2) cc_gam = cc_gam - cc_error(2)* (cc_gam-cc_gam/1.000001)/(cc_error(2)-cc_error(1))

     IF ( j == 1 .AND. ABS(cc_error(1)) > 1.E-15 ) GO TO 131
     IF ( cc_it >= 10 ) THEN
        WRITE (*,*) ' cc error'
        STOP
     END IF
     IF ( j /= 1 ) GO TO 13

     fcc= - alpha_2 / PI/4.0 /rho* (cc_gam*cc_q2  &
          + pi/2.0/cc_delt *cc_omeg*p_n**2 ) + cc_gam**3 /pi/3.0/rho
     !  Restricted Primitive Model
     !      fcc=-(3.0*ii_par*ii_par+6.0*ii_par+2.0  &
     !                             -2.0*(1.0+2.0*ii_par)**1.50)  &
     !               /(12.0*PI*rho *cc_sig1**3 )

     !      fcc = x_ions * fcc

     my2dd(3) = (parame(3,6))**2  *1.E-19 /(KBOL*t)
     my2dd(3) = (1.84)**2  *1.E-19 /(kbol*t)

     rhs   = 12.0 * PI * rho * x(3) * my2dd(3)
     lam_s = 1.0
12   CONTINUE
     lambda = (rhs/((lam_s+2.0)**2 ) + 16.0/((1.0+lam_s)**4 ) )**0.5
     IF ( ABS(lam_s-lambda) > 1.E-10 )THEN
        lam_s = ( lambda + lam_s ) / 2.0
        GO TO 12
     END IF

     ! f_cd = -(ii_par*ii_par)/(4.0*PI*rho*m_mean *cc_sig1**3 )  &
     !         *(dielec-1.0)/(1.0 + parame(3,2)/cc_sig1/lambda)
     ! write (*,*) ' ',f_cd,fcc,x_ions
     ! f_cd = f_cd/(1.0 - fcc/f_cd)
     ! fcc = 0.0

  END IF


END SUBROUTINE F_ION_ION_PrimMSA


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_ION_ION_nonPrimMSA ( fdd, fqq, fdq, fcc )

  USE EOS_VARIABLES, ONLY: nc, ncomp, x, parame, mseg

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdd, fqq, fdq, fcc
  !-----------------------------------------------------------------------------
  INTEGER                                :: dipole
  !REAL                                   :: A_MSA !, A_CC, A_CD, A_DD, U_MSA, chempot
  REAL, DIMENSION(nc)                    :: x_export, msegm
  !-----------------------------------------------------------------------------

  dipole = 0
  IF ( SUM( parame(1:ncomp,6) ) > 1.E-5 ) dipole = 1

  IF ( dipole /= 0 ) THEN      ! alternatively ions and dipoles = 1
     fdd = 0.0
     fqq = 0.0
     fdq = 0.0
     fcc = 0.0
     msegm(:)    = mseg(:)     ! the entries of the vector mseg and x are changed
     x_export(:) = x(:)        ! in SEMIRESTRICTED because the ions should be positioned first
     ! that is why dummy vectors msegm and x_export are defined
     !CALL SEMIRESTRICTED (A_MSA,A_CC,A_CD,A_DD,U_MSA,  &
     !                    chempot,ncomp,parame,t,eta,x_export,msegm,0)
     !fdd = A_MSA
     write (*,*) 'why are individual contrib. A_CC,A_CD,A_DD not used'
     stop
  END IF

END SUBROUTINE F_ION_ION_nonPrimMSA


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_LC_MayerSaupe ( flc )

  USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, phas, t, rho, eta,  &
       x, mseg, E_lc, S_lc, dhs

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: flc
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k
  INTEGER                                :: liq_crystal, count_lc, steps_lc
  REAL                                   :: alpha_lc, tolerance, deltay
  REAL                                   :: integrand1, integrand2, accel_lc
  REAL                                   :: error_lc, u_term, sphase
  REAL, DIMENSION(nc)                    :: z_lc, S_lc1, S_lc2, sumu
  REAL, DIMENSION(nc,nc)                 :: u_lc, klc
  !-----------------------------------------------------------------------------
  INTEGER :: stabil
  COMMON /stabil / stabil
  !-----------------------------------------------------------------------------


  klc(1,2) = 0.0
  klc(2,1) = klc(1,2)

  alpha_lc = 1.0
  accel_lc = 4.0
  IF ( eta < 0.35 ) accel_lc = 1.3
  IF ( eta < 0.15 ) accel_lc = 1.0

  liq_crystal = 0
  DO i = 1, ncomp
     DO j = 1, ncomp
        E_lc(i,j) = (E_lc(i,i)*E_lc(j,j))**0.5 *(1.0-klc(i,j)) !combining rule
        ! E_LC(i,j)= ( E_LC(i,i)+E_LC(j,j) ) * 0.5   !combining rule
        ! S_LC(i,j)= ( S_LC(i,i)+S_LC(j,j) ) * 0.5   !combining rule
        IF (E_lc(i,j) /= 0.0) liq_crystal = 1
     END DO
  END DO
  ! S_LC(1,2) = 0.0
  ! S_LC(2,1) = S_LC(1,2)
  ! E_LC(1,2) = 60.0
  ! E_LC(2,1) = E_LC(1,2)

  IF ( liq_crystal == 1 .AND. phas == 1 .AND. stabil == 0 ) THEN

     count_lc = 0
     tolerance = 1.E-6

     steps_lc = 200
     deltay = 1.0 / REAL(steps_lc)

     ! --- dimensionless function U_LC repres. anisotr. intermolecular interactions in l.c.

     DO i = 1, ncomp
        DO j = 1, ncomp
           u_lc(i,j) = 2.0/3.0*pi*mseg(i)*mseg(j) *(0.5*(dhs(i)+dhs(j)))**3  &   ! sig_ij(i,j)**3
                *(E_lc(i,j)/t+S_lc(i,j))*rho
        END DO
     END DO


     DO i=1,ncomp
        ! S_lc2(i) = 0.0         !for isotropic
        S_lc2(i) = 0.5           !for nematic
        S_lc1(i) = S_lc2(i)
     END DO

1    CONTINUE

     DO i = 1, ncomp
        IF (S_lc2(i) <= 0.3) S_lc1(i) = S_lc2(i)
        IF (S_lc2(i) > 0.3) S_lc1(i) = S_lc1(i) + (S_lc2(i)-S_lc1(i))*accel_lc
     END DO

     count_lc = count_lc + 1

     !--------------------------------------------------------------------------
     ! single-particle orientation partition function Z_LC in liquid crystals
     !--------------------------------------------------------------------------

     DO i = 1, ncomp
        sumu(i) = 0.0
        DO j = 1, ncomp
           sumu(i) = sumu(i) + x(j)*u_lc(i,j)*S_lc1(j)
        END DO
     END DO

     DO i = 1, ncomp
        z_lc(i) = 0.0
        integrand1 = EXP(-0.5*sumu(i))           !eq. for Z_LC with y=0
        DO k = 1, steps_lc
           integrand2 = EXP(0.5*sumu(i)*(3.0*(deltay*REAL(k)) **2 -1.0))
           z_lc(i) = z_lc(i) + (integrand1 + integrand2)/2.0*deltay
           integrand1 = integrand2
        END DO    !k-index integration
     END DO    !i-index Z_LC(i) calculation

     !--------------------------------------------------------------------------
     ! order parameter S_lc in liquid crystals
     !--------------------------------------------------------------------------

     error_lc = 0.0
     DO i = 1, ncomp
        S_lc2(i) = 0.0
        integrand1 = -1.0/z_lc(i)*0.5*EXP(-0.5*sumu(i))  !for S_lc with y=0
        DO k = 1, steps_lc
           integrand2 = 1.0/z_lc(i)*0.5*(3.0*(deltay*REAL(k))  &
                **2 -1.0)*EXP(0.5*sumu(i)*(3.0 *(deltay*REAL(k))**2 -1.0))
           S_lc2(i) = S_lc2(i) + (integrand1 + integrand2)/2.0*deltay
           integrand1 = integrand2
        END DO    !k-index integration
        error_lc = error_lc + ABS(S_lc2(i)-S_lc1(i))
     END DO    !i-index Z_LC(i) calculation

     sphase = 0.0
     DO i = 1, ncomp
        sphase = sphase + S_lc2(i)
     END DO
     IF (sphase < 1.E-4) THEN
        error_lc = 0.0
        DO i = 1, ncomp
           S_lc2(i)= 0.0
           z_lc(i) = 1.0
        END DO
     END IF


     ! write (*,*) count_LC,S_lc2(1)-S_lc1(1),S_lc2(2)-S_lc1(2)
     IF (error_lc > tolerance .AND. count_lc < 400) GO TO 1
     ! write (*,*) 'done',eta,S_lc2(1),S_lc2(2)

     IF (count_lc == 400) WRITE (*,*) 'LC iteration not converg.'

     !--------------------------------------------------------------------------
     ! the anisotropic contribution to the Helmholtz energy
     !--------------------------------------------------------------------------

     u_term = 0.0
     DO i = 1, ncomp
        DO j = 1, ncomp
           u_term = u_term + 0.5*x(i)*x(j)*S_lc2(i) *S_lc2(j)*u_lc(i,j)
        END DO
     END DO

     flc = 0.0
     DO i = 1, ncomp
        IF (z_lc(i) /= 0.0) flc = flc - x(i) * LOG(z_lc(i))
     END DO
     flc = flc + u_term

  END IF
  ! write (*,'(i2,i2,4(f15.8))') phas,stabil,flc,eta,S_lc2(1),x(1)


END SUBROUTINE F_LC_MayerSaupe


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_pert_theory ( fdsp )

  USE EOS_VARIABLES, ONLY: nc, PI, rho, eta, z0t, order1, order2
  USE EOS_NUMERICAL_DERIVATIVES, ONLY: disp_term
  USE DFT_MODULE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdsp
  !-----------------------------------------------------------------------------
  REAL                                   :: I1, I2
  REAL                                   :: z3, zms, c1_con, m_mean
  !-----------------------------------------------------------------------------

  ! caution: positive sign of correlation integral is used here !
  ! (the Helmholtz energy terms are written with a negative sign, while I1 and I2 are positive)

  IF (disp_term == 'PT1') THEN

     CALL f_dft ( I1, I2)
     c1_con = 0.0
     I2 = 0.0
     fdsp  = + ( - 2.0*PI*rho*I1*order1 )

  ELSEIF (disp_term == 'PT2') THEN

     CALL f_dft ( I1, I2)
     z3 = eta
     zms = 1.0 - z3
     m_mean = z0t / ( PI / 6.0 )
     c1_con = 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3**2 )/zms**4  &
          + (1.0 - m_mean)*( 20.0*z3 -27.0*z3**2 +12.0*z3**3 -2.0*z3**4 )  &
          /(zms*(2.0-z3))**2 )
     fdsp  = + ( - 2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2 )

  ELSEIF (disp_term == 'PT_MIX') THEN

     CALL f_pert_theory_mix ( fdsp )

  ELSEIF (disp_term == 'PT_MF') THEN

     ! mean field theory
     I1 = - ( - 8.0/9.0 - 4.0/9.0*(rc**-9 -3.0*rc**-3 ) - tau_cut/3.0*(rc**3 -1.0) )
     fdsp  = + ( - 2.0*PI*rho*I1*order1 )
     write (*,*) 'caution: not thoroughly checked and tested'

  ELSE
     write (*,*) 'define the type of perturbation theory'
     stop
  END IF

  ! I1 = I1 + 4.0/9.0*(2.5**-9 -3.0*2.5**-3 )
  ! fdsp  = + ( - 2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2 )

END SUBROUTINE F_pert_theory

END MODULE EOS_F_CONTRIBUTIONS
