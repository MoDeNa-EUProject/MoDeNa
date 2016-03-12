!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module DFT_MODULE
!
! This module contains parameters and variables for DFT calculations.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module DFT_MODULE

  use PARAMETERS, only: nc
  implicit none
  save

  INTEGER, PARAMETER                    :: NDFT = 4000
  INTEGER, PARAMETER                    :: r_grid = 200
  INTEGER                               :: kmax, den_step
  LOGICAL                               :: shift, WCA, MFT
  REAL                                  :: rc, rg, dzr, dzp, tau_cut
  REAL                                  :: d_hs, dhs_st, z3t_st
  REAL                                  :: z_ges
  REAL, DIMENSION(r_grid)               :: x1a
  REAL, DIMENSION(NDFT)                 :: x2a
  REAL, DIMENSION(r_grid,NDFT)          :: ya, y1a, y2a, y12a
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub
  REAL                                  :: fres_temp

  REAL, DIMENSION(r_grid)               :: x1a_11
  REAL, DIMENSION(NDFT)                 :: x2a_11
  REAL, DIMENSION(r_grid,NDFT)          :: ya_11, y1a_11, y2a_11, y12a_11
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_11
  REAL, DIMENSION(r_grid)               :: x1a_12
  REAL, DIMENSION(NDFT)                 :: x2a_12
  REAL, DIMENSION(r_grid,NDFT)          :: ya_12, y1a_12, y2a_12, y12a_12
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_12
  REAL, DIMENSION(r_grid)               :: x1a_22
  REAL, DIMENSION(NDFT)                 :: x2a_22
  REAL, DIMENSION(r_grid,NDFT)          :: ya_22, y1a_22, y2a_22, y12a_22
  REAL, DIMENSION(r_grid,NDFT,4,4)      :: c_bicub_22

End Module DFT_MODULE



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DFT_RAD_INT
!
! This subroutine integrates the kernel of the perturbation theory in
! radial direction (more precisely in distance coordinate r^ hat). A
! non-mean field approach is taken, using a radial distribution function
! (currently at a Percus-Yevick level).
!
! The first and second order contributions to the perturbation theory
! are calculated - although the second order contribution is rarely used.
!
! ToDo: comment the subroutine and remove the GOTO constructs
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE DFT_RAD_INT (i,j,ih,zz1,rhop,f_int_r,my_int_r,  &
     f_int2_r,my_int2_r,my_int3_r)

  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: j
  INTEGER, INTENT(IN OUT)                :: ih
  REAL, INTENT(IN)                       :: zz1
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: f_int_r
  REAL, INTENT(OUT)                      :: my_int_r
  REAL, INTENT(OUT)                      :: f_int2_r
  REAL, INTENT(OUT)                      :: my_int2_r
  REAL, INTENT(OUT)                      :: my_int3_r

  !-----------------------------------------------------------------------------
  INTEGER                                :: k

  REAL                                   :: fint0, fint0_2, fint1, fint1_2
  REAL                                   :: myint0, myint0_2, myint0_3, myint1
  REAL                                   :: myint1_2, myint1_3
  REAL                                   :: dg_drho, dg_dr
  REAL                                   :: rad, xg, rdf, rho_bar, ua, rs
  REAL                                   :: analytic1, analytic2, tau_rs
  LOGICAL                                :: shortcut
  !-----------------------------------------------------------------------------

  shortcut = .true.
  fint0    = rc * tau_cut                 ! first order term
  fint0_2  = rc * tau_cut*tau_cut         ! 2nd order
  myint0   = rc * tau_cut                 ! first order term
  myint0_2 = rc * tau_cut*tau_cut         ! 2nd order
  myint0_3 = 0.0                          ! 2nd order

  f_int_r  = 0.0
  f_int2_r = 0.0
  my_int_r = 0.0
  my_int2_r= 0.0
  my_int3_r= 0.0
  !-------------this block only speeds up the integration
  IF ( shortcut ) THEN
     rs = MAX( rg, zz1 ) ! +dzr
     IF ( rs > rc ) WRITE (*,*) 'error !!!!'
     analytic1 = 0.4*rs**-10 - 0.4*rc**-10 - rs**-4  + rc**-4
     analytic2 = 16.0/22.0 * (rs**-22 - rc**-22 )  &
          -2.0*rs**-16 +2.0*rc**-16  +1.6*rs**-10 - 1.6*rc**-10
     f_int_r  = f_int_r  + analytic1
     f_int2_r = f_int2_r + analytic2
     my_int_r = my_int_r + analytic1
     my_int2_r= my_int2_r+ analytic2
     IF ( rs == zz1 ) GO TO 10
     tau_rs  = 4.0 * ( rs**-12 - rs**-6 )
     fint0   = rs * tau_rs
     fint0_2 = rs * tau_rs*tau_rs
     myint0  = rs * tau_rs
     myint0_2= rs * tau_rs*tau_rs
     rad = rs   ! the simple integration scheme: set to rc
     k = 0 + NINT( (rc-rs)/dzr )   ! in simple scheme: set to 0
  ELSE
     rad = rc
     k = 0
  END IF
  !-------------
  !rad = rc
  !k = 0

4 CONTINUE
  rad = rad - dzr
  k = k + 1
  ua = 4.0 * ( rad**-12 -rad**-6 )
  xg = rad / d_hs
  rho_bar = ( rhop(i) + rhop(j) )/2.0
  rdf = 1.0
  dg_drho = 0.0
  IF ( rad <= rg ) THEN
     CALL BI_CUB_SPLINE (rho_bar,xg,ya,x1a,x2a,y1a,y2a,y12a,  &
          c_bicub,rdf,dg_drho,dg_dr,den_step,ih,k)
     !              write (*,*) 'vor ',rdf,dg_drho
     !              CALL BI_CUB_SPLINE2 (rho_bar,xg,ya,x1a,x2a,y1a,
     !     &                      rdf,dg_drho,den_step,ih,k)
     !              write (*,*) 'nach',rdf,dg_drho
     !              read (*,*)
  END IF

  fint1    = rdf * rad * ua
  fint1_2  = rdf * rad * ua * ua
  myint1   = rad * (rdf + 0.5*rhop(i)*dg_drho) * ua
  myint1_2 = rdf * rad * ua * ua
  myint1_3 = dg_drho * rad * ua * ua
  ! intf(k) = fint1
  f_int_r  = f_int_r  + dzr * (fint1 + fint0)/2.0
  f_int2_r = f_int2_r + dzr * (fint1_2 + fint0_2)/2.0
  my_int_r = my_int_r + dzr * (myint1 + myint0)/2.0
  my_int2_r= my_int2_r+ dzr * (myint1_2 + myint0_2)/2.0
  my_int3_r= my_int3_r+ dzr * (myint1_3 + myint0_3)/2.0

  fint0   = fint1
  fint0_2 = fint1_2
  myint0  = myint1
  myint0_2= myint1_2
  myint0_3= myint1_3
  IF ( zz1 >= 1.0 .AND. rad-dzr+1.E-8 >= zz1 ) GO TO 4  ! integration down to ABS(zz1)
  IF ( zz1 < 1.0  .AND. rad-dzr+1.E-8 >= 1.0 ) GO TO 4  ! integration down to ABS(zz1) but for r^hat<1, g(r)=0 (so stop at r^hat=1)

  ! IF (k.GT.30) THEN
  !   stepno = k
  !   CALL SPLINE_PARA (dzr,intf,utri,stepno)
  !   CALL SPLINE_INT (f_int_r,dzr,intf,utri,stepno)
  ! ENDIF

10 CONTINUE

  analytic1 = 4.0/10.0*rc**-10  - rc**-4
  analytic2 = 16.0/22.0*rc**-22 - 2.0*rc**-16 + 1.6*rc**-10
  f_int_r  = f_int_r  + analytic1
  f_int2_r = f_int2_r + analytic2
  my_int_r = my_int_r + analytic1
  my_int2_r= my_int2_r+ analytic2

END SUBROUTINE DFT_RAD_INT


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE f_dft ( I1_dft, I2_dft )

  USE EOS_VARIABLES, ONLY: nc, PI, rho, parame
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: I1_dft
  REAL, INTENT(OUT)                      :: I2_dft

  !-----------------------------------------------------------------------------
  INTEGER                                :: k,ih
  ! REAL                                   :: z3
  REAL                                   :: ua, ua_c, ua_2, ua_c_2, rm
  REAL                                   :: int10, int11, int20, int21
  REAL                                   :: dg_drho
  REAL                                   :: rad, xg, rdf, rho_st, msegm
  REAL                                   :: sig_ij
  REAL                                   :: dg_dr, dzr_org !,rdf_d
  ! REAL                                 :: intgrid(0:NDFT),intgri2(0:NDFT)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! constants
  !-----------------------------------------------------------------------------
  msegm = parame(1,1)
  rho_st = rho * parame(1,2)**3

  ua_c = 4.0 * ( rc**-12 - rc**-6 )
  ua_c_2 = ua_c * ua_c
  rm = 2.0**(1.0/6.0)

  int10     = rc*rc* ua_c
  int20     = rc*rc* ua_c_2
  ! intgrid(0)= int10
  ! intgri2(0)= int20


  sig_ij = parame(1,2)


  I1_dft = 0.0
  I2_dft = 0.0
  rad    = rc
  !dzr    = dzp / 2.0    ! this line is obsolete. dzr is defined in DFT-nMF2 (dimensionless)
  dzr_org= dzr
  k = 0
  ih = 85

  DO WHILE ( rad-dzr+1.E-9 >= 1.0 )

     rad = rad - dzr
     ! IF (rad <= 8.0) dzr = dzp
     ! IF (rad <= rg) dzr = dzp/2.0
     k = k + 1
     xg = rad / dhs_st
     ua = 4.0 * ( rad**-12 - rad**-6 )
     ua_2 = ua * ua
     rdf  = 1.0
     dg_drho = 0.0
     IF ( rad <= rg ) THEN
        CALL BI_CUB_SPLINE (rho_st,xg,ya,x1a,x2a,y1a,y2a,y12a,  &
             c_bicub,rdf,dg_drho,dg_dr,den_step,ih,k)
     END IF

     int11 = rdf*rad*rad* ua
     int21 = rdf*rad*rad* ua_2
     I1_dft= I1_dft + dzr*(int11+int10)/2.0
     I2_dft= I2_dft + dzr*(int21+int20)/2.0
     int10 = int11
     int20 = int21

  END DO

  dzr = dzr_org

  ! stepno = k
  ! CALL SPLINE_PARA (dzr,intgrid,utri,stepno)
  ! CALL SPLINE_INT (I1,dzr,intgrid,utri,stepno)

  !     caution: 1st order integral is in F_EOS.f defined with negative sign
  I1_dft= - I1_dft - ( 4.0/9.0 * rc**-9 - 4.0/3.0 * rc**-3 )

  ! CALL SPLINE_PARA (dzr,intgri2,utri,stepno)
  ! CALL SPLINE_INT (I2,dzr,intgri2,utri,stepno)

  I2_dft = I2_dft + 16.0/21.0 * rc**-21 - 32.0/15.0 * rc**-15 + 16.0/9.0 * rc**-9


END SUBROUTINE f_dft


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE rdf_matrix (msegm)

  USE PARAMETERS, ONLY: PI
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                        :: msegm

  !-----------------------------------------------------------------------------
  INTEGER                                 :: i, k = 0
  REAL                                    :: rdf, rad, xg, rho_rdf, z3
  !-----------------------------------------------------------------------------


  DO i = 1, den_step

     ! rho_rdf= rhob(1) + (rhob(2)-rhob(1))*REAL(i)/den_step
     rho_rdf = 1.E-5 + (1.0) * REAL(i-1) / REAL(den_step-1)   ! segment density, rho_s*sigma**3
     rho_rdf = rho_rdf / msegm                                ! molar density, rho_m*sigma**3
     rad = rc
     k = 0

     do while ( rad-dzr+1.E-9 >= 0.95 )

        rad = rad - dzr
        k = k + 1
        xg = rad / d_hs
        z3 = rho_rdf * msegm * PI/6.0  * d_hs**3
        ya(i,k) = 1.0
        IF ( xg <= rg .AND. z3 > 0.0 ) CALL rdf_int ( z3, msegm, xg, rdf )
        IF ( xg <= rg ) ya(i,k) = rdf
        ! ya(i,k) = y(x1a(i), x2a(k))  with x1a: density-vector, x2a: r-vector
        x1a(i) = rho_rdf
        x2a(k) = xg

     end do

  END DO

  if ( xg > 1.0 ) stop 'rdf_matrix: 0.95*sigma is too high for lower bound'

  WRITE (*,*) ' done with calculating g(r)',d_hs

  kmax = k
  CALL bicub_derivative ( ya, x1a, x2a, y1a, y2a, y12a, den_step, kmax )
  CALL bicub_c ( ya, x1a, x2a, y1a, y2a, y12a, c_bicub, den_step, kmax )

END SUBROUTINE rdf_matrix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BI_CUB_SPLINE2 ( rho_rdf, xg, ya, x1a, x2a, y1a,  &
     rdf, dg_drho, i_max, ih, k )

  USE DFT_MODULE, ONLY: NDFT, r_grid
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rho_rdf
  REAL, INTENT(IN OUT)                   :: xg
  REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
  REAL, INTENT(IN)                       :: x1a(r_grid)
  REAL, INTENT(IN OUT)                   :: x2a(NDFT)
  REAL, INTENT(IN)                       :: y1a(r_grid,NDFT)
  REAL, INTENT(OUT)                      :: rdf
  REAL, INTENT(OUT)                      :: dg_drho
  INTEGER, INTENT(IN OUT)                :: i_max
  !INTEGER, INTENT(IN OUT)                :: k_max
  INTEGER, INTENT(OUT)                   :: ih
  INTEGER, INTENT(IN OUT)                :: k

  !-----------------------------------------------------------------------------
  REAL                                   :: as,bs,cs,ds,rho_t,del_rho
  !REAL :: y2a(r_grid,NDFT),y12a(r_grid,NDFT), c_bicub(r_grid,NDFT,4,4),c(4,4)
  !-----------------------------------------------------------------------------

  IF ( rho_rdf < x1a(1) ) THEN
     dg_drho = 0.0
     rdf = 1.0
     RETURN
  END IF
  IF ( x1a(ih) <= rho_rdf.AND.rho_rdf < x1a(ih+1) ) GO TO 10
  IF ( ih > 2 ) THEN
     IF ( x1a(ih-1) <= rho_rdf.AND.rho_rdf < x1a(ih) ) THEN
        ih = ih - 1
        GO TO 10
     END IF
  END IF
  !      write (*,*) 'in ',ih
  CALL hunt(x1a,i_max,rho_rdf,ih)
  !      write (*,*) 'out',ih
10 CONTINUE
  IF ( (x2a(k) > (xg + 1.E-10)) .OR. (x2a(k) < (xg - 1.E-10)) )  &
       WRITE (*,*) 'error in bi-cubic-spline',x2a(k),xg


  del_rho = x1a(ih+1) - x1a(ih)

  as =  ya(ih,k)
  bs =  y1a(ih,k)
  cs =  3.0*(ya(ih+1,k)-as)/del_rho**2 - (2.0*bs+y1a(ih+1,k))/del_rho
  ds = -2.0*(ya(ih+1,k)-as)/del_rho**3 + (bs+y1a(ih+1,k))/del_rho**2

  rho_t = rho_rdf - x1a(ih)
  rdf = as + bs*rho_t + cs*rho_t*rho_t + ds*rho_t**3
  dg_drho = bs + 2.0*cs*rho_t + 3.0*ds*rho_t*rho_t

END SUBROUTINE BI_CUB_SPLINE2



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BI_CUB_SPLINE ( rho_rdf, xg, ya, x1a, x2a, y1a, y2a, y12a,  &
     c_bicub, rdf, dg_drho, dg_dr, i_max, ih, k )

  USE DFT_MODULE, ONLY: NDFT, r_grid
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: rho_rdf
  REAL, INTENT(IN OUT)                   :: xg
  REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
  REAL, INTENT(IN)                       :: x1a(r_grid)
  REAL, INTENT(IN)                       :: x2a(NDFT)
  REAL, INTENT(IN)                       :: y1a(r_grid,NDFT)
  REAL, INTENT(IN)                       :: y2a(r_grid,NDFT)
  REAL, INTENT(IN)                       :: y12a(r_grid,NDFT)
  REAL, INTENT(IN)                       :: c_bicub(r_grid,NDFT,4,4)
  REAL, INTENT(OUT)                      :: rdf
  REAL, INTENT(OUT)                      :: dg_drho
  REAL, INTENT(OUT)                      :: dg_dr
  INTEGER, INTENT(IN OUT)                :: i_max
  !INTEGER, INTENT(IN OUT)                :: k_max
  INTEGER, INTENT(OUT)                   :: ih
  INTEGER, INTENT(IN)                    :: k

  !-----------------------------------------------------------------------------
  INTEGER                                :: m, n

  REAL                                   :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
  REAL                                   :: c(4,4)
  !-----------------------------------------------------------------------------

  IF ( rho_rdf < x1a(1) ) THEN
     dg_drho = 0.0
     dg_dr = 0.0
     rdf = 1.0
     RETURN
  END IF
  IF ( x1a(ih) <= rho_rdf .AND. rho_rdf < x1a(ih+1) ) GO TO 10
  IF ( ih > 2 ) THEN
     IF ( x1a(ih-1) <= rho_rdf .AND. rho_rdf < x1a(ih) ) THEN
        ih = ih - 1
        GO TO 10
     END IF
  END IF
  ! write (*,*) 'in ',ih
  CALL hunt ( x1a, i_max, rho_rdf, ih )
  ! write (*,*) 'out',ih
10 CONTINUE
  IF ( x2a(k) /= xg ) THEN
     ! write (*,*) 'error bi-cubic-spline',k,x2a(k),xg
     ! DO k=1,NDFT
     !   write (*,*) k,x2a(k)
     ! ENDDO
     ! stop
  END IF



  y(1) = ya(ih,k)
  y(2) = ya(ih+1,k)
  y(3) = ya(ih+1,k+1)
  y(4) = ya(ih,k+1)

  y1(1) = y1a(ih,k)
  y1(2) = y1a(ih+1,k)
  y1(3) = y1a(ih+1,k+1)
  y1(4) = y1a(ih,k+1)

  y2(1) = y2a(ih,k)
  y2(2) = y2a(ih+1,k)
  y2(3) = y2a(ih+1,k+1)
  y2(4) = y2a(ih,k+1)

  y12(1) = y12a(ih,k)
  y12(2) = y12a(ih+1,k)
  y12(3) = y12a(ih+1,k+1)
  y12(4) = y12a(ih,k+1)

  x1l = x1a(ih)
  x1u = x1a(ih+1)
  x2l = x2a(k)
  x2u = x2a(k+1)

  DO m = 1, 4
     DO n = 1, 4
        c(m,n) = c_bicub( ih, k, m, n )
     END DO
  END DO
  CALL bcuint ( x1l, x1u, x2l, x2u, rho_rdf, xg, c, rdf, dg_drho, dg_dr )

END SUBROUTINE BI_CUB_SPLINE



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bicub_derivative ( ya, x1a, x2a, y1a, y2a, y12a, i_max, k_max )

  USE DFT_MODULE, ONLY: NDFT, r_grid
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
  REAL, INTENT(IN)                       :: x1a(r_grid)
  REAL, INTENT(IN)                       :: x2a(NDFT)
  REAL, INTENT(OUT)                      :: y1a(r_grid,NDFT)
  REAL, INTENT(OUT)                      :: y2a(r_grid,NDFT)
  REAL, INTENT(OUT)                      :: y12a(r_grid,NDFT)
  INTEGER, INTENT(IN)                    :: i_max
  INTEGER, INTENT(IN)                    :: k_max

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, k
  !-----------------------------------------------------------------------------


  DO i = 2, i_max-1
     DO k = 2, k_max-1
        y1a(i,k) = (ya(i+1,k)-ya(i-1,k)) / (x1a(i+1)-x1a(i-1))
        y2a(i,k) = (ya(i,k+1)-ya(i,k-1)) / (x2a(k+1)-x2a(k-1))
        y12a(i,k)= (ya(i+1,k+1)-ya(i+1,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
             /((x1a(i+1)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
     END DO
  END DO

  i = 1
  DO k = 1, k_max
     y1a(i,k) = 0.0
     y2a(i,k) = 0.0
     y12a(i,k)= 0.0
  END DO

  k = 1
  DO i = 1, i_max
     y1a(i,k) = 0.0
     y2a(i,k) = 0.0
     y12a(i,k)= 0.0
  END DO


  i = i_max
  DO k = 2, k_max-1
     y1a(i,k) = (ya(i,k)-ya(i-1,k)) / (x1a(i)-x1a(i-1))
     y2a(i,k) = (ya(i,k+1)-ya(i,k-1)) / (x2a(k+1)-x2a(k-1))
     y12a(i,k)= (ya(i,k+1)-ya(i,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
          /((x1a(i)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
  END DO


  k = k_max
  DO i = 2, i_max-1
     y1a(i,k) = (ya(i+1,k)-ya(i-1,k)) / (x1a(i+1)-x1a(i-1))
     y2a(i,k) = (ya(i,k)-ya(i,k-1)) / (x2a(k)-x2a(k-1))
     y12a(i,k)= (ya(i+1,k)-ya(i+1,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
          /((x1a(i+1)-x1a(i-1))*(x2a(k)-x2a(k-1)))
  END DO

  k = k_max
  i = i_max
  y1a(i,k) = (ya(i,k)-ya(i-1,k)) / (x1a(i)-x1a(i-1))
  y2a(i,k) = (ya(i,k)-ya(i,k-1)) / (x2a(k)-x2a(k-1))
  y12a(i,k)= (ya(i,k)-ya(i,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
       /((x1a(i)-x1a(i-1))*(x2a(k)-x2a(k-1)))

END SUBROUTINE bicub_derivative



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bicub_c ( ya, x1a, x2a, y1a, y2a, y12a, c_bicub, i_max, k_max )

  USE DFT_MODULE, ONLY: NDFT, r_grid
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
  REAL, INTENT(IN)                       :: x1a(r_grid)
  REAL, INTENT(IN)                       :: x2a(NDFT)
  REAL, INTENT(IN)                       :: y1a(r_grid,NDFT)
  REAL, INTENT(IN)                       :: y2a(r_grid,NDFT)
  REAL, INTENT(IN)                       :: y12a(r_grid,NDFT)
  REAL, INTENT(OUT)                      :: c_bicub(r_grid,NDFT,4,4)
  INTEGER, INTENT(IN)                    :: i_max
  INTEGER, INTENT(IN)                    :: k_max

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, k, m, n
  REAL                                   :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
  REAL                                   :: c(4,4)
  !-----------------------------------------------------------------------------

  DO i = 1, i_max-1
     DO k = 1, k_max-1
        y(1)=ya(i,k)
        y(2)=ya(i+1,k)
        y(3)=ya(i+1,k+1)
        y(4)=ya(i,k+1)

        y1(1)=y1a(i,k)
        y1(2)=y1a(i+1,k)
        y1(3)=y1a(i+1,k+1)
        y1(4)=y1a(i,k+1)

        y2(1)=y2a(i,k)
        y2(2)=y2a(i+1,k)
        y2(3)=y2a(i+1,k+1)
        y2(4)=y2a(i,k+1)

        y12(1)=y12a(i,k)
        y12(2)=y12a(i+1,k)
        y12(3)=y12a(i+1,k+1)
        y12(4)=y12a(i,k+1)

        x1l=x1a(i)
        x1u=x1a(i+1)
        x2l=x2a(k)
        x2u=x2a(k+1)

        CALL bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
        DO m=1,4
           DO n=1,4
              c_bicub(i,k,m,n)=c(m,n)
           END DO
        END DO

     END DO
  END DO

END SUBROUTINE bicub_c



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!SUBROUTINE bcuint ( y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, c, ansy, ansy1, ansy2 )
SUBROUTINE bcuint ( x1l, x1u, x2l, x2u, x1, x2, c, ansy, ansy1, ansy2 )

  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !REAL, INTENT(IN OUT)                   :: y(4)
  !REAL, INTENT(IN OUT)                   :: y1(4)
  !REAL, INTENT(IN OUT)                   :: y2(4)
  !REAL, INTENT(IN OUT)                   :: y12(4)
  REAL, INTENT(IN OUT)                   :: x1l
  REAL, INTENT(IN OUT)                   :: x1u
  REAL, INTENT(IN OUT)                   :: x2l
  REAL, INTENT(IN OUT)                   :: x2u
  REAL, INTENT(IN OUT)                   :: x1
  REAL, INTENT(IN OUT)                   :: x2
  REAL, INTENT(IN)                       :: c(4,4)
  REAL, INTENT(OUT)                      :: ansy
  REAL, INTENT(OUT)                      :: ansy1
  REAL, INTENT(OUT)                      :: ansy2

  !-----------------------------------------------------------------------------
  !U    USES bcucof
  INTEGER                                :: i
  REAL                                   :: t, u
  !-----------------------------------------------------------------------------

  ! call bcucof ( y, y1, y2, y12, x1u-x1l, x2u-x2l, c )

  IF ( x1u == x1l .OR. x2u == x2l ) call paus ('bad input in bcuint')
  t = (x1-x1l) / (x1u-x1l)
  u = (x2-x2l) / (x2u-x2l)
  ansy  = 0.0
  ansy2 = 0.0
  ansy1 = 0.0
  DO  i = 4, 1, -1
     ansy  = t *ansy  + ( (c(i,4)*u + c(i,3))*u+c(i,2) )*u + c(i,1)
     ansy2 = t *ansy2 + ( 3.*c(i,4)*u+2.*c(i,3) )*u + c(i,2)
     ansy1 = u *ansy1 + ( 3.*c(4,i)*t+2.*c(3,i) )*t + c(2,i)
  END DO
  ansy1 = ansy1 / (x1u-x1l)
  ansy2 = ansy2 / (x2u-x2l)

END SUBROUTINE bcuint



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: y(4)
  REAL, INTENT(IN)                       :: y1(4)
  REAL, INTENT(IN)                       :: y2(4)
  REAL, INTENT(IN)                       :: y12(4)
  REAL, INTENT(IN)                       :: d1
  REAL, INTENT(IN)                       :: d2
  REAL, INTENT(OUT)                      :: c(4,4)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i,j,k,l
  REAL                                   :: d1d2,xx,cl(16),wt(16,16),x(16)
  SAVE wt
  DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*  &
       0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,  &
       1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,  &
       -6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,  &
       10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,  &
       -2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,  &
       2,-2,2*0,-1,1/
  !-----------------------------------------------------------------------------

  d1d2 = d1 * d2
  DO  i = 1, 4
     x(i)    = y(i)
     x(i+4)  = y1(i)*d1
     x(i+8)  = y2(i)*d2
     x(i+12) = y12(i)*d1d2
  END DO
  DO i = 1, 16
     xx = 0.0
     DO  k = 1, 16
        xx = xx + wt(i,k) * x(k)
     END DO
     cl(i) = xx
  END DO
  l = 0
  DO i = 1, 4
     DO j = 1, 4
        l = l + 1
        c(i,j) = cl(l)
     END DO
  END DO

END SUBROUTINE bcucof


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE hunt
!
! Given an array xx(1:n), and given a value x, returns a value jlo
! such that x is between xx(jlo) and xx(jlo+1). xx(1:n) must be
! monotonic, either increasing or decreasing. jlo=0 or jlo=n is
! returned to indicate that x is out of range. jlo on input is taken
! as the initial guess for jlo on output.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE hunt ( xx, n, x, jlo )

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n
  INTEGER, INTENT(OUT)                   :: jlo
  REAL, INTENT(IN)                       :: xx(n)
  REAL                                   :: x

  !-----------------------------------------------------------------------------
  INTEGER                                :: inc,jhi,jm
  LOGICAL                                :: ascnd
  !-----------------------------------------------------------------------------

  ascnd = xx(n) >= xx(1)
  IF( jlo <= 0 .OR. jlo > n ) THEN
     jlo = 0
     jhi = n + 1
     GO TO 3
  END IF
  inc = 1
  IF( x >= xx(jlo) .EQV. ascnd ) THEN
1    jhi = jlo + inc
     IF ( jhi > n ) THEN
        jhi = n + 1
     ELSE IF ( x >= xx(jhi) .EQV. ascnd ) THEN
        jlo = jhi
        inc = inc + inc
        GO TO 1
     END IF
  ELSE
     jhi = jlo
2    jlo = jhi - inc
     IF ( jlo < 1 ) THEN
        jlo = 0
     ELSE IF ( x < xx(jlo) .EQV. ascnd ) THEN
        jhi = jlo
        inc = inc + inc
        GO TO 2
     END IF
  END IF
3 IF (jhi-jlo == 1 ) THEN
     IF ( x == xx(n)) jlo = n - 1
     IF ( x == xx(1) ) jlo = 1
     RETURN
  END IF
  jm = ( jhi + jlo ) / 2
  IF ( x >= xx(jm) .EQV. ascnd ) THEN
     jlo = jm
  ELSE
     jhi = jm
  END IF
  GO TO 3
END SUBROUTINE hunt




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE dens_inv_coeff (rhob)

  USE BASIC_VARIABLES
  USE DFT_MODULE
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rhob(2)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, no_step_l, no_step_h
  REAL                                   :: my_eta_div,den_st,den_min,den_max, den_div
  REAL                                   :: lnphi(np,nc)
  !-----------------------------------------------------------------------------
  REAL :: rho_array1(200),rho_array2(200),  &
       my_rho_array1(200),my_rho_array2(200)
  COMMON /deninv/ rho_array1,rho_array2,my_rho_array1,  &
       my_rho_array2,my_eta_div,den_min,den_max, no_step_l,no_step_h
  !-----------------------------------------------------------------------------

  den_min = rhob(2) * z3t_st * 0.8
  den_max = rhob(1) * z3t_st * 1.1
  den_div = 0.04

  no_step_l = 80
  no_step_h = 200

  my_eta_div = 0.0
  IF (den_min < den_div) THEN
     DO i=1,no_step_l
        dense(1) = den_min + (den_div-den_min) * REAL(i-1)/REAL(no_step_l-1)
        densta(1)= dense(1)
        rho_array1(i) = dense(1)
        CALL fugacity (lnphi)
        my_rho_array1(i) = EXP( lnphi(1,1) ) *dense(1)/z3t_st / parame(1,2)**3    ! myrho = my_res(T,rho) + ln(rho)
     END DO
     my_eta_div = my_rho_array1(no_step_l)
  END IF

  den_st = MAX( den_min, den_div )
  DO i = 1, no_step_h
     dense(1) = den_st + (den_max-den_st) * REAL(i-1)/REAL(no_step_h-1)
     densta(1)= dense(1)
     rho_array2(i) = dense(1)
     CALL fugacity (lnphi)
     my_rho_array2(i) = lnphi(1,1) + LOG( dense(1)/z3t_st /parame(1,2)**3 )  ! myrho = my_res(T,rho) + ln(rho)
     IF ( i >= 2 ) THEN
        IF ( (my_rho_array2(i)- my_rho_array2(i-1)) < 0.0 ) THEN
           call paus ('DENS_INV_COEFF: derivative negative')
        END IF
     END IF
  END DO

END SUBROUTINE dens_inv_coeff


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION DENS_INVERS ( rhob, my_rho )

  USE BASIC_VARIABLES
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rhob(2)
  REAL, INTENT(IN OUT)                   :: my_rho

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, no_step_l, no_step_h
  REAL                                   :: my_eta_div, den_min, den_max, eta
  !-----------------------------------------------------------------------------
  REAL :: rho_array1(200),rho_array2(200),  &
       my_rho_array1(200),my_rho_array2(200)
  COMMON /deninv/ rho_array1,rho_array2,my_rho_array1,  &
       my_rho_array2,my_eta_div,den_min,den_max, no_step_l,no_step_h
  !-----------------------------------------------------------------------------

  ! den_min = rhob(2)*z3t_st * 0.8
  ! den_max = rhob(1)*z3t_st * 1.1
  ! den_div = 0.04

  ! no_step_l = 40
  ! no_step_h = 100

  IF ( EXP(my_rho) < my_eta_div ) THEN
     DO i = 2, no_step_l
        IF ( my_rho_array1(i) > EXP(my_rho) ) THEN
           eta = rho_array1(i-1) + (rho_array1(i)-rho_array1(i-1))  &
                / (my_rho_array1(i)-my_rho_array1(i-1))  &
                * (EXP(my_rho)-my_rho_array1(i-1))
           DENS_INVERS = eta/z3t_st
           GO TO 1
        END IF
     END DO
     WRITE (*,*) 'error 1',EXP(my_rho),my_eta_div
     ! call paus (' ')
  END IF
1 CONTINUE

  IF ( EXP(my_rho) >= my_eta_div ) THEN
     DO i = 2, no_step_h
        IF ( my_rho_array2(i) > my_rho ) THEN
           eta = rho_array2(i-1)+(rho_array2(i)-rho_array2(i-1))  &
                / (my_rho_array2(i)-my_rho_array2(i-1)) *(my_rho-my_rho_array2(i-1))
           DENS_INVERS = eta/z3t_st
           GO TO 2
        END IF
     END DO
     WRITE (*,*) 'error 2',my_rho,my_eta_div
     ! call paus (' ')
  END IF
2 CONTINUE

  IF ( DENS_INVERS < (rhob(2)*0.8) ) THEN
     WRITE (*,*) 'lower bound',DENS_INVERS,my_rho
     DENS_INVERS = rhob(2)
  END IF
  IF ( DENS_INVERS > (rhob(1)*1.1) ) THEN
     WRITE (*,*) 'upper bound',DENS_INVERS,my_rho
     DENS_INVERS = rhob(1)
  END IF

END FUNCTION DENS_INVERS
