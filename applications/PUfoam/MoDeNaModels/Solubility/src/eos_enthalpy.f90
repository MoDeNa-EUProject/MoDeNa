!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine enthalpy_ig_QSPR ( cpig, hig, sig )

  use parameters, only: RGAS
  use basic_variables, only: nc, np, ncomp, nphas, t, p, xi, parame
  implicit none

  !-----------------------------------------------------------------------------
  real, dimension(np), intent(out)       :: cpig
  real, dimension(np), intent(out)       :: hig
  real, dimension(np), intent(out)       :: sig

  !-----------------------------------------------------------------------------
  integer                                :: i, ph
  real                                   :: p0
  real                                   :: t0
  real, dimension(nc)                    :: cp_ig, h_ig, s_ig
  real                                   :: m_par, sig_par, eps_par
  real                                   :: dip_par, quad_par, eps_ab_par, kap_ab_par
  real                                   :: p1, p2, p3, p4, p5, p6
  real                                   :: c1a, c2a, c3a, c4a, c5a, c6a
  real                                   :: c1b, c2b, c3b, c4b, c5b, c6b
  real                                   :: a_coeff, b_coeff
  real                                   :: ICPC300, ICPC400
  real                                   :: T1, T2
  !-----------------------------------------------------------------------------

  t0  = 298.0

  do i = 1, ncomp

     m_par   = parame(i, 1) 
     sig_par = parame(i, 2) 
     eps_par = parame(i, 3) 

     eps_ab_par = parame(i,14)
     kap_ab_par = parame(i,13)
     dip_par  = parame(i, 6)
     quad_par = parame(i, 7)

     p1 = m_par * eps_par / t
     p2 = m_par * sig_par**3
     p3 = m_par * sig_par**3 * eps_par / t
     p4 = m_par * sig_par**6 * kap_ab_par *( EXP( eps_ab_par/t ) - 1.0 )
     p5 = m_par * sig_par**3 * quad_par
     p6 = 1.0

     T1 = 300.0
     T2 = 400.0

     !--------------------------------------------------------------------------
     ! nAnP
     !--------------------------------------------------------------------------

     ! ICPC @  300 K

     c1a = -5763.04893
     c2a = 1232.30607
     c3a = -239.3513996
     c4a = 0.0
     c5a = 0.0
     c6a = -15174.28321

     ! ICPC @  400 K

     c1b = -8171.26676935062
     c2b = 1498.01217504596
     c3b = -315.515836223387
     c4b = 0.0
     c5b = 0.0
     c6b = -19389.5468655708

     !--------------------------------------------------------------------------
     ! nAP
     !--------------------------------------------------------------------------

     !  
     !!    ICPC @  300K:
     !    c1a = 5177.19095226181
     !    c2a = 919.565206504576
     !    c3a = -108.829105648889
     !    c4a = 0
     !    c5a = -3.93917830677682
     !    c6a = -13504.5671858292
     !    
     !!    ICPC @  400K:
     !    c1b = 10656.1018362315
     !    c2b = 1146.10782703748
     !    c3b = -131.023645998081
     !    c4b = 0
     !    c5b = -9.93789225413177
     !    c6b = -24430.12952497

     !--------------------------------------------------------------------------
     ! AP      
     !--------------------------------------------------------------------------

     !!    ICPC @  300K:
     !    c1a = 3600.32322462175
     !    c2a = 1006.20461224949
     !    c3a = -151.688378113974
     !    c4a = 7.81876773647109d-07
     !    c5a = 8.01001754473385
     !    c6a = -8959.37140957179

     !!    ICPC @  400K:
     !    c1b = 7248.0697641199
     !    c2b = 1267.44346171358
     !    c3b = -208.738557800023
     !    c4b = 0.000170238690157906
     !    c5b = -6.7841792685616
     !    c6b = -12669.4196622924

     ICPC300 = c1a*P1/300.0 + c2a*P2 + c3a*P3/300.0 + c4a*P4/300.0 + c5a*P5 + c6a*P6 

     ICPC300 = ICPC300 /1000.0    

     ICPC400 = c1b*P1/400.0 + c2b*P2 + c3b*P3/400.0 + c4b*P4/400.0 + c5b*P5 + c6b*P6 

     ICPC400 = ICPC400 / 1000.0 

     !--------------------------------------------------------------------------
     ! linear approximation
     !--------------------------------------------------------------------------

     b_coeff = ( icpc400 - icpc300 ) / ( T2 - T1 )
     a_coeff = icpc300 - b_coeff * T1

     t0 = 298.0
     p0 = 1.E5

     cp_ig(i)= a_coeff + b_coeff * t      
     h_ig(i) = a_coeff * (t-t0) + 0.5 * b_coeff * ( t**2 - t0**2 )
     s_ig(i) = a_coeff * log( t/t0 ) + b_coeff * (t-t0) - RGAS * log( p/p0 )

  end do

  !-----------------------------------------------------------------------------
  ! Ideal gas properties
  !-----------------------------------------------------------------------------
  do ph = 1, nphas

     hig(ph)  = sum( xi( ph, 1:ncomp ) * h_ig( 1:ncomp ) )
     cpig(ph) = sum( xi( ph, 1:ncomp ) *cp_ig( 1:ncomp ) )

     sig(ph)  = 0.0
     do i = 1, ncomp
        if ( xi(ph,i) > 0.0 ) sig(ph) = sig(ph) + xi(ph,i) * ( s_ig(i) - RGAS * LOG( xi(ph,i) ) )
     end do

  end do

end subroutine enthalpy_ig_QSPR



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine enthalpy_ig ( coef_ig, cpig, hig, sig )

  use parameters, only: RGAS
  use basic_variables, only: nc, np, ncomp, nphas, t, p, xi
  implicit none

  !-----------------------------------------------------------------------------
  real, dimension(nc,5), intent(in)      :: coef_ig
  real, dimension(np), intent(out)       :: cpig
  real, dimension(np), intent(out)       :: hig
  real, dimension(np), intent(out)       :: sig

  !-----------------------------------------------------------------------------
  integer                                :: i, ph
  real                                   :: p0, t0, t02, t03, t04
  real, dimension(nc)                    :: cp_ig, h_ig, s_ig
  !-----------------------------------------------------------------------------

  p0 = 1.E5

  t0  = 298.15
  t02 = t0 * t0
  t03 = t02 * t0
  t04 = t02 * t02
  do i = 1, ncomp

     cp_ig(i)= coef_ig(i,1) + coef_ig(i,2)*t + coef_ig(i,3)*t*t + coef_ig(i,4)*t**3
     h_ig(i) = coef_ig(i,1)*(t-t0) + coef_ig(i,2)/2.0*(t*t-t02)  &
             + coef_ig(i,3)/3.0*(t**3-t03) + coef_ig(i,4)/4.0*(t**4 - t04)
     s_ig(i) = coef_ig(i,1) * log( t/t0 ) + coef_ig(i,2) * (t-t0)  &
             + coef_ig(i,3)/2.0*(t*t-t02) + coef_ig(i,4)/3.0*(t**3-t03) - RGAS * log( p/p0 )

  end do

  do ph = 1, nphas

     hig(ph)  = sum( xi( ph, 1:ncomp ) * h_ig( 1:ncomp ) )
     cpig(ph) = sum( xi( ph, 1:ncomp ) *cp_ig( 1:ncomp ) )
     sig(ph)  = 0.0
     do i = 1, ncomp
        if ( xi(ph,i) > 0.0 ) sig(ph)  = sig(ph) + xi(ph,i) * ( s_ig(i) - RGAS * LOG( xi(ph,i) ) )
     end do

  end do

END SUBROUTINE enthalpy_ig


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE H_EOS

  USE PARAMETERS, ONLY: RGAS
  ! USE BASIC_VARIABLES, ONLY: mm           ! only for speed of sound 
  USE EOS_CONSTANTS
  USE EOS_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL                                   :: zges, df_dt, dfdr, ddfdrdr
  REAL                                   :: cv_res, df_dt2, df_drdt
  REAL                                   :: fact, dist, t_tmp, rho_0
  REAL                                   :: fdr1, fdr2, fdr3, fdr4
  ! REAL                                   :: cp_ig_mix, w_square, mm_mean   ! only for speed of sound


  INTEGER                                :: i, m
  REAL                                   :: dhsdt(nc), dhsdt2(nc)
  REAL                                   :: z0, z1, z2, z3, z1tdt, z2tdt, z3tdt
  REAL                                   :: z1dt, z2dt, z3dt, zms, gii
  REAL                                   :: fhsdt, fhsdt2
  REAL                                   :: fchdt, fchdt2
  REAL                                   :: fdspdt, fdspdt2
  REAL                                   :: fhbdt, fhbdt2
  REAL                                   :: sumseg, I1, I2, I1dt, I2dt, I1dt2, I2dt2
  REAL                                   :: c1_con, c2_con, c3_con, c1_dt, c1_dt2
  REAL                                   :: z1tdt2, z2tdt2, z3tdt2
  REAL                                   :: z1dt2, z2dt2, z3dt2

  INTEGER                                :: j, k, l, no, ass_cnt, max_eval
  LOGICAL                                :: assoc
  REAL                                   :: dij, dijdt, dijdt2
  REAL                                   :: gij1dt, gij2dt, gij3dt
  REAL                                   :: gij1dt2, gij2dt2, gij3dt2
  REAL, DIMENSION(nc,nc)                 :: gijdt, gijdt2, kap_hb
  REAL, DIMENSION(nc,nc,nsite,nsite)     :: ass_d_dt, ass_d_dt2, eps_hb, delta, deltadt, deltadt2
  REAL, DIMENSION(nc,nsite)              :: mxdt, mxdt2, mx_itr, mx_itrdt, mx_itrdt2
  REAL                                   :: attenu, tol, suma, sumdt, sumdt2, err_sum

  INTEGER                                :: dipole
  REAL                                   :: fdddt, fdddt2
  REAL, DIMENSION(nc)                    :: my2dd, my0
  REAL, DIMENSION(nc,nc)                 :: idd2, idd2dt, idd2dt2, idd4, idd4dt, idd4dt2
  REAL, DIMENSION(nc,nc,nc)              :: idd3, idd3dt, idd3dt2
  REAL                                   :: factor2, factor3
  REAL                                   :: fdd2, fdd3, fdd2dt, fdd3dt, fdd2dt2, fdd3dt2
  REAL                                   :: eij, xijmt, xijkmt

  INTEGER                                :: qudpole
  REAL                                   :: fqqdt, fqqdt2
  REAL, DIMENSION(nc)                    :: qq2
  REAL, DIMENSION(nc,nc)                 :: iqq2, iqq2dt, iqq2dt2, iqq4, iqq4dt, iqq4dt2
  REAL, DIMENSION(nc,nc,nc)              :: iqq3, iqq3dt, iqq3dt2
  REAL                                   :: fqq2, fqq2dt, fqq2dt2, fqq3, fqq3dt, fqq3dt2

  INTEGER                                :: dip_quad
  REAL                                   :: fdqdt, fdqdt2
  REAL, DIMENSION(nc)                    :: myfac, q_fac
  REAL, DIMENSION(nc,nc)                 :: idq2, idq2dt, idq2dt2, idq4, idq4dt, idq4dt2
  REAL, DIMENSION(nc,nc,nc)              :: idq3, idq3dt, idq3dt2
  REAL                                   :: fdq2, fdq2dt, fdq2dt2, fdq3, fdq3dt, fdq3dt2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! Initializing
  !-----------------------------------------------------------------------------

  CALL PERTURBATION_PARAMETER

  rho = eta/z3t
  z0 = z0t*rho
  z1 = z1t*rho
  z2 = z2t*rho
  z3 = z3t*rho

  sumseg = z0t / (PI/6.0)
  zms    = 1.0 - z3


  !-----------------------------------------------------------------------------
  ! first and second derivative of f to density (dfdr,ddfdrdr)
  !-----------------------------------------------------------------------------

  CALL P_EOS

  zges = (pges * 1.E-30)/(kbol*t*rho)

  dfdr    = pges/(eta*rho*(kbol*t)/1.E-30)
  ddfdrdr = pgesdz/(eta*rho*(kbol*t)/1.E-30) - dfdr*2.0/eta - 1.0/eta**2 


  !-----------------------------------------------------------------------------
  ! Helmholtz Energy f/kT = fres
  !-----------------------------------------------------------------------------

  CALL F_EOS


  !-----------------------------------------------------------------------------
  ! derivative of some auxilliary properties to temperature
  !-----------------------------------------------------------------------------

  DO i = 1, ncomp
     dhsdt(i)=parame(i,2) *(-3.0*parame(i,3)/t/t)*0.12*EXP(-3.0*parame(i,3)/t)
     dhsdt2(i) = dhsdt(i)*3.0*parame(i,3)/t/t  &
          + 6.0*parame(i,2)*parame(i,3)/t**3  *0.12*EXP(-3.0*parame(i,3)/t)
  END DO

  z1tdt = 0.0
  z2tdt = 0.0
  z3tdt = 0.0
  DO i = 1, ncomp
     z1tdt = z1tdt + x(i) * mseg(i) * dhsdt(i)
     z2tdt = z2tdt + x(i) * mseg(i) * 2.0*dhs(i)*dhsdt(i)
     z3tdt = z3tdt + x(i) * mseg(i) * 3.0*dhs(i)*dhs(i)*dhsdt(i)
  END DO
  z1dt  = PI / 6.0*z1tdt *rho
  z2dt  = PI / 6.0*z2tdt *rho
  z3dt  = PI / 6.0*z3tdt *rho


  z1tdt2 = 0.0
  z2tdt2 = 0.0
  z3tdt2 = 0.0
  DO i = 1, ncomp
     z1tdt2 = z1tdt2 + x(i)*mseg(i)*dhsdt2(i)
     z2tdt2 = z2tdt2 + x(i)*mseg(i)*2.0 *( dhsdt(i)*dhsdt(i) +dhs(i)*dhsdt2(i) )
     z3tdt2 = z3tdt2 + x(i)*mseg(i)*3.0 *( 2.0*dhs(i)*dhsdt(i)*  &
          dhsdt(i) +dhs(i)*dhs(i)*dhsdt2(i) )
  END DO
  z1dt2  = PI / 6.0*z1tdt2 *rho
  z2dt2  = PI / 6.0*z2tdt2 *rho
  z3dt2  = PI / 6.0*z3tdt2 *rho


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT hard spheres to temp. (fhsdt)
  !-----------------------------------------------------------------------------

  fhsdt = 6.0/PI/rho*(  3.0*(z1dt*z2+z1*z2dt)/zms + 3.0*z1*z2*z3dt/zms/zms  &
       + 3.0*z2*z2*z2dt/z3/zms/zms  &
       + z2**3 *(2.0*z3*z3dt-z3dt*zms)/(z3*z3*zms**3 )  &
       + (3.0*z2*z2*z2dt*z3-2.0*z2**3 *z3dt)/z3**3 *LOG(zms)  &
       + (z0-z2**3 /z3/z3)*z3dt/zms  )

  fhsdt2= 6.0/PI/rho*( 3.0*(z1dt2*z2+2.0*z1dt*z2dt+z1*z2dt2)/zms  &
       + 6.0*(z1dt*z2+z1*z2dt)*z3dt/zms/zms  &
       + 3.0*z1*z2*z3dt2/zms/zms + 6.0*z1*z2*z3dt*z3dt/zms**3   &
       + 3.0*z2*(2.0*z2dt*z2dt+z2*z2dt2)/z3/zms/zms  &
       - z2*z2*(6.0*z2dt*z3dt+z2*z3dt2)/(z3*z3*zms*zms)  &
       + 2.0*z2**3 *z3dt*z3dt/(z3**3  *zms*zms)  &
       - 4.0*z2**3 *z3dt*z3dt/(z3*z3 *zms**3 )  &
       + (12.0*z2*z2*z2dt*z3dt+2.0*z2**3 *z3dt2)/(z3*zms**3 )  &
       + 6.0*z2**3 *z3dt*z3dt/(z3*zms**4 )  &
       - 2.0*(3.0*z2*z2*z2dt/z3/z3-2.0*z2**3 *z3dt/z3**3 ) *z3dt/zms  &
       -(z2**3 /z3/z3-z0)*(z3dt2/zms+z3dt*z3dt/zms/zms)  &
       +  ( (6.0*z2*z2dt*z2dt+3.0*z2*z2*z2dt2)/z3/z3  &
       - (12.0*z2*z2*z2dt*z3dt+2.0*z2**3 *z3dt2)/z3**3   &
       + 6.0*z2**3 *z3dt*z3dt/z3**4  )* LOG(zms)    )


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT of chain term to T (fchdt)
  !-----------------------------------------------------------------------------

  fchdt  = 0.0
  fchdt2 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        dij=dhs(i)*dhs(j)/(dhs(i)+dhs(j))
        dijdt =(dhsdt(i)*dhs(j) + dhs(i)*dhsdt(j)) / (dhs(i)+dhs(j))  &
             - dhs(i)*dhs(j)/(dhs(i)+dhs(j))**2  *(dhsdt(i)+dhsdt(j))
        dijdt2=(dhsdt2(i)*dhs(j) + 2.0*dhsdt(i)*dhsdt(j)  &
             + dhs(i)*dhsdt2(j)) / (dhs(i)+dhs(j))  &
             - 2.0*(dhsdt(i)*dhs(j) + dhs(i)*dhsdt(j))  &
             / (dhs(i)+dhs(j))**2  *(dhsdt(i)+dhsdt(j))  &
             + 2.0* dhs(i)*dhs(j) / (dhs(i)+dhs(j))**3   &
             * (dhsdt(i)+dhsdt(j))**2   &
             - dhs(i)*dhs(j)/(dhs(i)+dhs(j))**2  *(dhsdt2(i)+dhsdt2(j))
        gij1dt = z3dt/zms/zms
        gij2dt = 3.0*( z2dt*dij+z2*dijdt )/zms/zms +6.0*z2*dij*z3dt/zms**3 
        gij3dt = 4.0*(dij*z2)* (dijdt*z2 + dij*z2dt) /zms**3   &
             + 6.0*(dij*z2)**2  * z3dt /zms**4 
        gij1dt2 = z3dt2/zms/zms +2.0*z3dt*z3dt/zms**3 
        gij2dt2 = 3.0*( z2dt2*dij+2.0*z2dt*dijdt+z2*dijdt2 )/zms/zms  &
             + 6.0*( z2dt*dij+z2*dijdt )/zms**3  * z3dt  &
             + 6.0*(z2dt*dij*z3dt+z2*dijdt*z3dt+z2*dij*z3dt2) /zms**3   &
             + 18.0*z2*dij*z3dt*z3dt/zms**4 
        gij3dt2 = 4.0*(dijdt*z2+dij*z2dt)**2  /zms**3   &
             + 4.0*(dij*z2)* (dijdt2*z2+2.0*dijdt*z2dt+dij*z2dt2) /zms**3   &
             + 24.0*(dij*z2) *(dijdt*z2+dij*z2dt)/zms**4  *z3dt  &
             + 6.0*(dij*z2)**2  * z3dt2 /zms**4   &
             + 24.0*(dij*z2)**2  * z3dt*z3dt /zms**5 
        gijdt(i,j)  = gij1dt + gij2dt + gij3dt
        gijdt2(i,j)  = gij1dt2 + gij2dt2 + gij3dt2
     END DO
  END DO

  DO i = 1, ncomp
     gii = 1.0/zms + 3.0*dhs(i)/2.0*z2/zms/zms + 2.0*dhs(i)*dhs(i)/4.0*z2*z2/zms**3 
     fchdt = fchdt + x(i) * (1.0-mseg(i)) * gijdt(i,i) / gii
     fchdt2= fchdt2+ x(i) * (1.0-mseg(i))  &
          * (gijdt2(i,i) / gii - (gijdt(i,i)/gii)**2 )
  END DO


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT dispersion term to T (fdspdt)
  !-----------------------------------------------------------------------------

  I1 = 0.0
  I2 = 0.0
  I1dt = 0.0
  I2dt = 0.0
  I1dt2= 0.0
  I2dt2= 0.0
  DO m = 0, 6
     I1   = I1   + apar(m)*z3**m
     I2   = I2   + bpar(m)*z3**m
     I1dt = I1dt + apar(m)*z3dt*REAL(m)*z3**(m-1)
     I2dt = I2dt + bpar(m)*z3dt*REAL(m)*z3**(m-1)
     I1dt2= I1dt2+ apar(m)*z3dt2*REAL(m)*z3**(m-1)  &
          + apar(m)*z3dt*z3dt *REAL(m)*REAL(m-1)*z3**(m-2)
     I2dt2= I2dt2+ bpar(m)*z3dt2*REAL(m)*z3**(m-1)  &
          + bpar(m)*z3dt*z3dt *REAL(m)*REAL(m-1)*z3**(m-2)
  END DO

  c1_con= 1.0/ (  1.0 + sumseg*(8.0*z3-2.0*z3**2 )/zms**4   &
       + (1.0 - sumseg)*(20.0*z3-27.0*z3**2  +12.0*z3**3 -2.0*z3**4 )  &
       /(zms*(2.0-z3))**2   )
  c2_con= - c1_con*c1_con *(sumseg*(-4.0*z3**2 +20.0*z3+8.0)/zms**5   &
       + (1.0 - sumseg) *(2.0*z3**3 +12.0*z3**2 -48.0*z3+40.0)  &
       /(zms*(2.0-z3))**3  )
  c3_con= 2.0 * c2_con*c2_con/c1_con - c1_con*c1_con  &
       *( sumseg*(-12.0*z3**2 +72.0*z3+60.0)/zms**6 + (1.0 - sumseg)  &
       *(-6.0*z3**4 -48.0*z3**3 +288.0*z3**2  -480.0*z3+264.0)  &
       /(zms*(2.0-z3))**4  )
  c1_dt = c2_con*z3dt
  c1_dt2 = c3_con*z3dt*z3dt + c2_con*z3dt2

  fdspdt  = - 2.0*PI*rho*(I1dt-I1/t)*order1  &
       - PI*rho*sumseg*(c1_dt*I2+c1_con*I2dt-2.0*c1_con*I2/t)*order2

  fdspdt2 = - 2.0*PI*rho*(I1dt2-2.0*I1dt/t+2.0*I1/t/t)*order1  &
       - PI*rho*sumseg*order2*( c1_dt2*I2 +2.0*c1_dt*I2dt -4.0*c1_dt*I2/t  &
       + 6.0*c1_con*I2/t/t -4.0*c1_con*I2dt/t +c1_con*I2dt2)


  !-----------------------------------------------------------------------------
  ! 1st & 2nd derivative of f/kT association term to T (fhbdt)
  !-----------------------------------------------------------------------------

  fhbdt  = 0.0
  fhbdt2 = 0.0
  assoc = .false.
  DO i = 1, ncomp
     IF ( nhb_typ(i) /= 0 ) assoc = .true.
  END DO
  IF (assoc) THEN

     DO i = 1, ncomp
        IF ( nhb_typ(i) /= 0 ) THEN
           kap_hb(i,i) = parame(i,13)
           no = 0
           DO j = 1,nhb_typ(i)
              DO k = 1,nhb_typ(i)
                 eps_hb(i,i,j,k) = parame(i,(14+no))
                 no = no + 1
              END DO
           END DO
           DO j = 1,nhb_typ(i)
              no = no + 1
           END DO
        ELSE
           kap_hb(i,i) = 0.0
           DO k = 1,nsite
              DO l = 1,nsite
                 eps_hb(i,i,k,l) = 0.0
              END DO
           END DO
        END IF
     END DO

     DO i = 1, ncomp
        DO j = 1, ncomp
           IF ( i /= j .AND. (nhb_typ(i) /= 0 .AND. nhb_typ(j) /= 0) ) THEN
              kap_hb(i,j)= (kap_hb(i,i)*kap_hb(j,j))**0.5  &
                   *((parame(i,2)*parame(j,2))**3 )**0.5  &
                   /(0.5*(parame(i,2)+parame(j,2)))**3 
              ! kap_hb(i,j)= kap_hb(i,j)*(1.0-k_kij)
              DO k = 1,nhb_typ(i)
                 DO l = 1,nhb_typ(j)
                    IF (k /= l) THEN
                       eps_hb(i,j,k,l)=(eps_hb(i,i,k,l)+eps_hb(j,j,l,k))/2.0
                       ! eps_hb(i,j,k,l)=eps_hb(i,j,k,l)*(1.0-eps_kij)
                    END IF
                 END DO
              END DO
           END IF
        END DO
     END DO
     IF (nhb_typ(1) == 3) THEN
        eps_hb(1,2,3,1)=0.5*(eps_hb(1,1,3,1)+eps_hb(2,2,1,2))
        eps_hb(2,1,1,3) = eps_hb(1,2,3,1)
     END IF
     IF (nhb_typ(2) == 3) THEN
        eps_hb(2,1,3,1)=0.5*(eps_hb(2,2,3,1)+eps_hb(1,1,1,2))
        eps_hb(1,2,1,3) = eps_hb(2,1,3,1)
     END IF

     DO i = 1, ncomp
        DO k = 1, nhb_typ(i)
           DO j = 1, ncomp
              DO l = 1, nhb_typ(j)
                 ! ass_d(i,j,k,l)=kap_hb(i,j) *sig_ij(i,j)**3 *(EXP(eps_hb(i,j,k,l)/t)-1.0)
                 ass_d_dt(i,j,k,l)= - kap_hb(i,j) *sig_ij(i,j)**3 *  &
                      eps_hb(i,j,k,l)/t/t*EXP(eps_hb(i,j,k,l)/t)
                 ass_d_dt2(i,j,k,l)= - kap_hb(i,j) *sig_ij(i,j)**3   &
                      * eps_hb(i,j,k,l)/t/t*EXP(eps_hb(i,j,k,l)/t)  &
                      * (-2.0/t - eps_hb(i,j,k,l)/t/t)
              END DO
           END DO
        END DO
     END DO

     DO i = 1, ncomp
        DO k = 1, nhb_typ(i)
           DO j = 1, ncomp
              DO l = 1, nhb_typ(j)
                 delta(i,j,k,l)=gij(i,j)*ass_d(i,j,k,l)
                 deltadt(i,j,k,l) = gijdt(i,j)*ass_d(i,j,k,l) + gij(i,j)*ass_d_dt(i,j,k,l)
                 deltadt2(i,j,k,l)= gijdt2(i,j)*ass_d(i,j,k,l)  &
                      + 2.0*gijdt(i,j)*ass_d_dt(i,j,k,l) +gij(i,j)*ass_d_dt2(i,j,k,l)
              END DO
           END DO
        END DO
     END DO


     !--- constants for iteration ----------------------------------------------
     attenu = 0.7
     tol = 1.E-10
     IF (eta < 0.2) tol = 1.E-11
     IF (eta < 0.01) tol = 1.E-12
     max_eval = 200

     !--- initialize mxdt(i,j) -------------------------------------------------
     DO i = 1, ncomp
        DO k = 1, nhb_typ(i)
           mxdt(i,k) = 0.0
           mxdt2(i,k) = 0.0
        END DO
     END DO


     !--- iterate over all components and all sites ----------------------------
     DO ass_cnt = 1, max_eval

        DO i = 1, ncomp
           DO k = 1, nhb_typ(i)
              suma  = 0.0
              sumdt = 0.0
              sumdt2= 0.0
              DO j = 1, ncomp
                 DO l = 1, nhb_typ(j)
                    suma  = suma  + x(j)*nhb_no(j,l)*  mx(j,l) *delta(i,j,k,l)
                    sumdt = sumdt + x(j)*nhb_no(j,l)*( mx(j,l) *deltadt(i,j,k,l)  &
                         + mxdt(j,l)*delta(i,j,k,l) )
                    sumdt2 = sumdt2 + x(j)*nhb_no(j,l)*( mx(j,l)*deltadt2(i,j,k,l)  &
                         + 2.0*mxdt(j,l)*deltadt(i,j,k,l) + mxdt2(j,l)*delta(i,j,k,l) )
                 END DO
              END DO
              mx_itr(i,k) = 1.0 / (1.0 + suma * rho)
              mx_itrdt(i,k)= - mx_itr(i,k)**2 * sumdt*rho
              mx_itrdt2(i,k)= +2.0*mx_itr(i,k)**3 * (sumdt*rho)**2 - mx_itr(i,k)**2 *sumdt2*rho
           END DO
        END DO

        err_sum = 0.0
        DO i = 1, ncomp
           DO k = 1, nhb_typ(i)
              err_sum = err_sum + ABS(mx_itr(i,k) - mx(i,k))  &
                   + ABS(mx_itrdt(i,k) - mxdt(i,k)) + ABS(mx_itrdt2(i,k) - mxdt2(i,k))
              mx(i,k) = mx_itr(i,k) * attenu + mx(i,k) * (1.0 - attenu)
              mxdt(i,k)=mx_itrdt(i,k)*attenu +mxdt(i,k)* (1.0 - attenu)
              mxdt2(i,k)=mx_itrdt2(i,k)*attenu +mxdt2(i,k)* (1.0 - attenu)
           END DO
        END DO
        IF(err_sum <= tol) GO TO 10

     END DO
     WRITE(6,*) 'CAL_PCSAFT: max_eval violated err_sum = ',err_sum,tol
     STOP
10   CONTINUE

     DO i = 1, ncomp
        DO k = 1, nhb_typ(i)
           ! fhb = fhb + x(i)* nhb_no(i,k)* ( 0.5 * ( 1.0 - mx(i,k) ) + LOG(mx(i,k)) )
           fhbdt = fhbdt  + x(i)*nhb_no(i,k) *mxdt(i,k)*(1.0/mx(i,k)-0.5)
           fhbdt2= fhbdt2 + x(i)*nhb_no(i,k) *(mxdt2(i,k)*(1.0/mx(i,k)-0.5)  &
                -(mxdt(i,k)/mx(i,k))**2 )
        END DO
     END DO

  END IF


  !-----------------------------------------------------------------------------
  ! derivatives of f/kT of dipole-dipole term to temp. (fdddt)
  !-----------------------------------------------------------------------------

  fdddt  = 0.0
  fdddt2 = 0.0
  dipole = 0
  DO i = 1, ncomp
     my2dd(i) = 0.0
     IF ( parame(i,6) /= 0.0 .AND. uij(i,i) /= 0.0 ) THEN
        dipole = 1
        my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**3 *1.E-30)
     END IF
     my0(i) = my2dd(i)        ! needed for dipole-quadrupole-term
  END DO

  IF ( dipole == 1 ) THEN
     DO i = 1, ncomp
        DO j = 1, ncomp
           idd2(i,j)   =0.0
           idd4(i,j)   =0.0
           idd2dt(i,j) =0.0
           idd4dt(i,j) =0.0
           idd2dt2(i,j)=0.0
           idd4dt2(i,j)=0.0
           DO m = 0, 4
              idd2(i,j)  = idd2(i,j)   +ddp2(i,j,m)*z3**m
              idd4(i,j)  = idd4(i,j)   +ddp4(i,j,m)*z3**m
              idd2dt(i,j)= idd2dt(i,j) +ddp2(i,j,m)*z3dt*REAL(m)*z3**(m-1)
              idd4dt(i,j)= idd4dt(i,j) +ddp4(i,j,m)*z3dt*REAL(m)*z3**(m-1)
              idd2dt2(i,j)=idd2dt2(i,j)+ddp2(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                   + ddp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
              idd4dt2(i,j)=idd4dt2(i,j)+ddp4(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                   + ddp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
           END DO
           DO k = 1, ncomp
              idd3(i,j,k)   =0.0
              idd3dt(i,j,k) =0.0
              idd3dt2(i,j,k)=0.0
              DO m = 0, 4
                 idd3(i,j,k)   = idd3(i,j,k)  +ddp3(i,j,k,m)*z3**m
                 idd3dt(i,j,k) = idd3dt(i,j,k)+ddp3(i,j,k,m)*z3dt*REAL(m) *z3**(m-1)
                 idd3dt2(i,j,k)= idd3dt2(i,j,k)+ddp3(i,j,k,m)*z3dt2*REAL(m)  &
                      *z3**(m-1) +ddp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1) *z3**(m-2)
              END DO
           END DO
        END DO
     END DO


     factor2= - PI *rho
     factor3= - 4.0 / 3.0 * PI**2 * rho**2 

     fdd2  =  0.0
     fdd3  =  0.0
     fdd2dt=  0.0
     fdd3dt=  0.0
     fdd2dt2= 0.0
     fdd3dt2= 0.0
     DO i = 1, ncomp
        DO j = 1, ncomp
           xijmt = x(i)*parame(i,3)*parame(i,2)**3 *x(j)*parame(j,3)*parame(j,2)**3   &
                / ((parame(i,2)+parame(j,2))/2.0)**3  *my2dd(i)*my2dd(j)
           eij = (parame(i,3)*parame(j,3))**0.5
           fdd2  = fdd2   +factor2* xijmt/t/t*(idd2(i,j)+eij/t*idd4(i,j))
           fdd2dt= fdd2dt+ factor2* xijmt/t/t*(idd2dt(i,j)-2.0*idd2(i,j)/t  &
                +eij/t*idd4dt(i,j)-3.0*eij/t/t*idd4(i,j))
           fdd2dt2=fdd2dt2+factor2*xijmt/t/t*(idd2dt2(i,j)-4.0*idd2dt(i,j)/t  &
                +6.0*idd2(i,j)/t/t+eij/t*idd4dt2(i,j)  &
                -6.0*eij/t/t*idd4dt(i,j)+12.0*eij/t**3 *idd4(i,j))
           DO k = 1, ncomp
              xijkmt=x(i)*parame(i,3)*parame(i,2)**3  &
                   *x(j)*parame(j,3)*parame(j,2)**3  &
                   *x(k)*parame(k,3)*parame(k,2)**3  &
                   /((parame(i,2)+parame(j,2))/2.0) /((parame(i,2)+parame(k,2))/2.0)  &
                   /((parame(j,2)+parame(k,2))/2.0) *my2dd(i)*my2dd(j)*my2dd(k)
              fdd3   =fdd3   +factor3*xijkmt/t**3 *idd3(i,j,k)
              fdd3dt =fdd3dt +factor3*xijkmt/t**3 * (idd3dt(i,j,k)-3.0*idd3(i,j,k)/t)
              fdd3dt2=fdd3dt2+factor3*xijkmt/t**3  &
                   *( idd3dt2(i,j,k)-6.0*idd3dt(i,j,k)/t+12.0*idd3(i,j,k)/t/t )
           END DO
        END DO
     END DO

     IF ( fdd2 < -1.E-100 .AND. fdd3 /= 0.0 ) THEN
        fdddt = fdd2* (fdd2*fdd2dt - 2.0*fdd3*fdd2dt+fdd2*fdd3dt) / (fdd2-fdd3)**2 
        fdddt2 = ( 2.0*fdd2*fdd2dt*fdd2dt +fdd2*fdd2*fdd2dt2  &
             -2.0*fdd2dt**2 *fdd3  -2.0*fdd2*fdd2dt2*fdd3 +fdd2*fdd2*fdd3dt2 )  &
             /(fdd2-fdd3)**2  + fdddt * 2.0*(fdd3dt-fdd2dt)/(fdd2-fdd3)
     END IF
  END IF


  !-----------------------------------------------------------------------------
  ! derivatives f/kT of quadrupole-quadrup. term to T  (fqqdt)
  !-----------------------------------------------------------------------------

  fqqdt  = 0.0
  fqqdt2 = 0.0
  qudpole = 0
  DO i = 1, ncomp
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     IF (qq2(i) /= 0.0) qudpole = 1
  END DO

  IF (qudpole == 1) THEN

     DO i = 1, ncomp
        DO j = 1, ncomp
           iqq2(i,j)   = 0.0
           iqq4(i,j)   = 0.0
           iqq2dt(i,j) = 0.0
           iqq4dt(i,j) = 0.0
           iqq2dt2(i,j)= 0.0
           iqq4dt2(i,j)= 0.0
           DO m = 0, 4
              iqq2(i,j)   = iqq2(i,j)  + qqp2(i,j,m)*z3**m
              iqq4(i,j)   = iqq4(i,j)  + qqp4(i,j,m)*z3**m
              iqq2dt(i,j) = iqq2dt(i,j)+ qqp2(i,j,m)*z3dt*REAL(m)*z3**(m-1)
              iqq4dt(i,j) = iqq4dt(i,j)+ qqp4(i,j,m)*z3dt*REAL(m)*z3**(m-1)
              iqq2dt2(i,j)= iqq2dt2(i,j)+qqp2(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                   + qqp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
              iqq4dt2(i,j)= iqq4dt2(i,j)+qqp4(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                   + qqp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
           END DO
           DO k = 1, ncomp
              iqq3(i,j,k)   =0.0
              iqq3dt(i,j,k) =0.0
              iqq3dt2(i,j,k)=0.0
              DO m = 0, 4
                 iqq3(i,j,k)   = iqq3(i,j,k)  + qqp3(i,j,k,m)*z3**m
                 iqq3dt(i,j,k) = iqq3dt(i,j,k)+ qqp3(i,j,k,m)*z3dt*REAL(m) * z3**(m-1)
                 iqq3dt2(i,j,k)= iqq3dt2(i,j,k)+qqp3(i,j,k,m)*z3dt2*REAL(m) * z3**(m-1)  &
                      + qqp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
              END DO
           END DO
        END DO
     END DO

     factor2 = -9.0/16.0 * PI *rho
     factor3 =  9.0/16.0 * PI**2 * rho**2 

     fqq2   = 0.0
     fqq3   = 0.0
     fqq2dt = 0.0
     fqq3dt = 0.0
     fqq2dt2= 0.0
     fqq3dt2= 0.0
     DO i = 1, ncomp
        DO j = 1, ncomp
           xijmt = x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5   &
                * x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /sig_ij(i,j)**7.0
           eij = (parame(i,3)*parame(j,3))**0.5
           fqq2  = fqq2   +factor2* xijmt/t/t*(iqq2(i,j)+eij/t*iqq4(i,j))
           fqq2dt= fqq2dt +factor2* xijmt/t/t*(iqq2dt(i,j)-2.0*iqq2(i,j)/t  &
                + eij/t*iqq4dt(i,j)-3.0*eij/t/t*iqq4(i,j))
           fqq2dt2=fqq2dt2+factor2*xijmt/t/t*(iqq2dt2(i,j)-4.0*iqq2dt(i,j)/t  &
                + 6.0*iqq2(i,j)/t/t+eij/t*iqq4dt2(i,j)  &
                - 6.0*eij/t/t*iqq4dt(i,j)+12.0*eij/t**3 *iqq4(i,j))
           DO k = 1, ncomp
              xijkmt = x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /sig_ij(i,j)**3   &
                   * x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /sig_ij(i,k)**3   &
                   * x(k)*uij(k,k)*qq2(k)*sig_ij(k,k)**5 /sig_ij(j,k)**3 
              fqq3   = fqq3   +factor3*xijkmt/t**3 *iqq3(i,j,k)
              fqq3dt = fqq3dt +factor3*xijkmt/t**3 *(iqq3dt(i,j,k)-3.0*iqq3(i,j,k)/t)
              fqq3dt2= fqq3dt2+factor3*xijkmt/t**3   &
                   * ( iqq3dt2(i,j,k)-6.0*iqq3dt(i,j,k)/t+12.0*iqq3(i,j,k)/t/t )
           END DO
        END DO
     END DO

     IF ( fqq2 /= 0.0 .AND. fqq3 /= 0.0 ) THEN
        fqqdt = fqq2* (fqq2*fqq2dt - 2.0*fqq3*fqq2dt+fqq2*fqq3dt) / (fqq2-fqq3)**2 
        fqqdt2 = ( 2.0*fqq2*fqq2dt*fqq2dt +fqq2*fqq2*fqq2dt2  &
             - 2.0*fqq2dt**2 *fqq3  -2.0*fqq2*fqq2dt2*fqq3 +fqq2*fqq2*fqq3dt2 )  &
             / (fqq2-fqq3)**2  + fqqdt * 2.0*(fqq3dt-fqq2dt)/(fqq2-fqq3)
     END IF

  END IF


  !-----------------------------------------------------------------------------
  ! derivatives f/kT of dipole-quadruppole term to T  (fdqdt)
  !-----------------------------------------------------------------------------

  fdqdt = 0.0
  fdqdt2= 0.0
  dip_quad = 0
  DO i = 1, ncomp
     DO j = 1, ncomp
        IF (parame(i,6) /= 0.0 .AND. parame(j,7) /= 0.0) dip_quad = 1
     END DO
     myfac(i) = parame(i,3) * parame(i,2)**4 * my0(i)
     q_fac(i) = parame(i,3) * parame(i,2)**4 * qq2(i)
  END DO

  IF (dip_quad == 1) THEN

     DO i = 1, ncomp
        DO j = 1, ncomp
           idq2(i,j)   = 0.0
           idq4(i,j)   = 0.0
           idq2dt(i,j) = 0.0
           idq4dt(i,j) = 0.0
           idq2dt2(i,j)= 0.0
           idq4dt2(i,j)= 0.0
           IF ( myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0 ) THEN
              DO m = 0, 4
                 idq2(i,j)   = idq2(i,j)  + dqp2(i,j,m)*z3**m
                 idq4(i,j)   = idq4(i,j)  + dqp4(i,j,m)*z3**m
                 idq2dt(i,j) = idq2dt(i,j)+ dqp2(i,j,m)*z3dt*REAL(m)*z3**(m-1)
                 idq4dt(i,j) = idq4dt(i,j)+ dqp4(i,j,m)*z3dt*REAL(m)*z3**(m-1)
                 idq2dt2(i,j)= idq2dt2(i,j)+dqp2(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                      + dqp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
                 idq4dt2(i,j)= idq4dt2(i,j)+dqp4(i,j,m)*z3dt2*REAL(m)*z3**(m-1)  &
                      + dqp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
              END DO

              DO k = 1, ncomp
                 idq3(i,j,k)   = 0.0
                 idq3dt(i,j,k) = 0.0
                 idq3dt2(i,j,k)= 0.0
                 IF ( myfac(k) /= 0.0 .OR. q_fac(k) /= 0.0 ) THEN
                    DO m = 0, 4
                       idq3(i,j,k)  = idq3(i,j,k)  + dqp3(i,j,k,m)*z3**m
                       idq3dt(i,j,k)= idq3dt(i,j,k)+ dqp3(i,j,k,m)*z3dt*REAL(m) *z3**(m-1)
                       idq3dt2(i,j,k)= idq3dt2(i,j,k)+dqp3(i,j,k,m)*z3dt2*REAL(m) *z3**(m-1)  &
                            + dqp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**(m-2)
                    END DO
                 END IF
              END DO
           END IF
        END DO
     END DO

     factor2= -9.0/4.0 * PI * rho
     factor3=  PI**2 * rho**2 

     fdq2  = 0.0
     fdq3  = 0.0
     fdq2dt= 0.0
     fdq3dt= 0.0
     fdq2dt2=0.0
     fdq3dt2=0.0
     DO i = 1, ncomp
        DO j = 1, ncomp
           IF ( myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0 ) THEN
              xijmt = x(i)*myfac(i) * x(j)*q_fac(j) /sig_ij(i,j)**5 
              eij = (parame(i,3)*parame(j,3))**0.5
              fdq2  = fdq2  + factor2* xijmt/t/t*(idq2(i,j)+eij/t*idq4(i,j))
              fdq2dt= fdq2dt+ factor2* xijmt/t/t*(idq2dt(i,j)-2.0*idq2(i,j)/t  &
                   + eij/t*idq4dt(i,j)-3.0*eij/t/t*idq4(i,j))
              fdq2dt2 = fdq2dt2+factor2*xijmt/t/t*(idq2dt2(i,j)-4.0*idq2dt(i,j)/t  &
                   + 6.0*idq2(i,j)/t/t+eij/t*idq4dt2(i,j)  &
                   - 6.0*eij/t/t*idq4dt(i,j)+12.0*eij/t**3 *idq4(i,j))
              DO k = 1, ncomp
                 IF ( myfac(k) /= 0.0 .OR. q_fac(k) /= 0.0 ) THEN
                    xijkmt= x(i)*x(j)*x(k)/(sig_ij(i,j)*sig_ij(i,k)*sig_ij(j,k))**2   &
                         * ( myfac(i)*q_fac(j)*myfac(k)  &
                         + myfac(i)*q_fac(j)*q_fac(k)*1.193735 )

                    fdq3  =fdq3  + factor3*xijkmt/t**3 *idq3(i,j,k)
                    fdq3dt=fdq3dt+ factor3*xijkmt/t**3 * (idq3dt(i,j,k)-3.0*idq3(i,j,k)/t)
                    fdq3dt2=fdq3dt2+factor3*xijkmt/t**3   &
                         *( idq3dt2(i,j,k)-6.0*idq3dt(i,j,k)/t+12.0*idq3(i,j,k)/t/t )
                 END IF
              END DO
           END IF
        END DO
     END DO

     IF (fdq2 /= 0.0 .AND. fdq3 /= 0.0) THEN
        fdqdt = fdq2* (fdq2*fdq2dt - 2.0*fdq3*fdq2dt+fdq2*fdq3dt) / (fdq2-fdq3)**2 
        fdqdt2 = ( 2.0*fdq2*fdq2dt*fdq2dt +fdq2*fdq2*fdq2dt2  &
             - 2.0*fdq2dt**2 *fdq3  -2.0*fdq2*fdq2dt2*fdq3 +fdq2*fdq2*fdq3dt2 )  &
             / (fdq2-fdq3)**2  + fdqdt * 2.0*(fdq3dt-fdq2dt)/(fdq2-fdq3)
     END IF

  END IF
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
  ! total derivative of fres/kT to temperature
  !-----------------------------------------------------------------------------

  df_dt = fhsdt + fchdt + fdspdt + fhbdt + fdddt + fqqdt + fdqdt



  !-----------------------------------------------------------------------------
  ! second derivative of fres/kT to T
  !-----------------------------------------------------------------------------

  df_dt2 = fhsdt2 +fchdt2 +fdspdt2 +fhbdt2 +fdddt2 +fqqdt2 +fdqdt2



  !-----------------------------------------------------------------------------
  ! derivatives of fres/kt to density and to T
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! the analytic derivative of fres/kT to (density and T) (df_drdt)
  ! is still missing. A numerical differentiation is implemented.
  !-----------------------------------------------------------------------------

  fact = 1.0
  dist = t * 100.E-5 * fact
  t_tmp  = t
  rho_0 = rho


  t  = t_tmp - 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t*rho_0
  CALL P_EOS
  fdr1  = pges / (eta*rho_0*(kbol*t)/1.E-30)
  t  = t_tmp - dist
  CALL PERTURBATION_PARAMETER
  eta = z3t*rho_0
  CALL P_EOS
  fdr2  = pges / (eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp + dist
  CALL PERTURBATION_PARAMETER
  eta = z3t*rho_0
  CALL P_EOS
  fdr3  = pges / (eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp + 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t*rho_0
  CALL P_EOS
  fdr4  = pges / (eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp
  CALL PERTURBATION_PARAMETER
  eta = z3t*rho_0
  CALL P_EOS


  df_drdt   = (-fdr4+8.0*fdr3-8.0*fdr2+fdr1)/(12.0*dist)




  !-----------------------------------------------------------------------------
  ! thermodynamic properties
  !-----------------------------------------------------------------------------

  s_res = ( - df_dt *t - fres )*RGAS + RGAS * LOG(zges)
  h_res = ( - t*df_dt  +  zges-1.0  ) * RGAS *t
  cv_res = - (   t*df_dt2 + 2.0*df_dt  ) * RGAS*t
  cp_res = cv_res - RGAS  + RGAS*(zges +eta*t*df_drdt)**2   &
       / (1.0 + 2.0*eta*dfdr +eta**2  *ddfdrdr)

  !-----------------------------------------------------------------------------
  ! speed of sound
  !-----------------------------------------------------------------------------

!!! cp_ig_mix = 1.57797619E-04 * t * t + 1.76863690E-01 * t + 2.98855179E+01      ! in units J/(mol*K)
  ! mm_mean = SUM( x(1:ncomp) * mm(1:ncomp) )

  ! w_square = 1.0 + 2.0*eta*dfdr + eta**2 *ddfdrdr  &
  !     + 1.0/( cp_ig_mix/RGAS - 1.0 + cv_res/RGAS )  &
  !     * ( zges + eta*t *df_drdt )**2
  ! speed_sound = SQRT( w_square * RGAS * t / (mm_mean/1000.0) )


  !-----------------------------------------------------------------------------
  ! write (*,*) 'df_... ', df_dt,df_dt2
  ! write (*,*) 'kreuz ', zges,eta*t*df_drdt,eta*dfdr, eta**2 *ddfdrdr
  ! write (*,*) 'h,cv,cp', h_res,cv_res,cp_res


END SUBROUTINE H_EOS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE H_EOS_num
!
! This subroutine calculates enthalpies and heat capacities (cp) by
! taking numerical derivatieves.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE H_EOS_num

  USE PARAMETERS, ONLY: RGAS
  USE EOS_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL                                   :: dist, fact, rho_0
  REAL                                   :: fres1, fres2, fres3, fres4, fres5
  REAL                                   :: f_1, f_2, f_3, f_4
  REAL                                   :: cv_res, t_tmp, zges
  REAL                                   :: df_dt, df_dtdt, df_drdt, dfdr, ddfdrdr
  !-----------------------------------------------------------------------------


  CALL PERTURBATION_PARAMETER
  rho_0 = eta/z3t


  fact = 1.0
  dist = t * 100.E-5 * fact

  t_tmp  = t

  t  = t_tmp - 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL F_EOS
  fres1  = fres
  t  = t_tmp - dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL F_EOS
  fres2  = fres
  t  = t_tmp + dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL F_EOS
  fres3  = fres
  t  = t_tmp + 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL F_EOS
  fres4  = fres
  t  = t_tmp
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL F_EOS
  fres5  = fres
  ! *(KBOL*T)/1.E-30

  zges = (p * 1.E-30)/(kbol*t*rho_0)


  df_dt   = (-fres4+8.0*fres3-8.0*fres2+fres1)/(12.0*dist)
  df_dtdt = (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1)  &
       /(12.0*(dist**2 ))


  s_res  =  (- df_dt -fres/t)*RGAS*t + RGAS * LOG(zges)
  h_res  =  ( - t*df_dt  +  zges-1.0  ) * RGAS*t
  cv_res = -(   t*df_dtdt + 2.0*df_dt  ) * RGAS*t



  t  = t_tmp - 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL P_EOS
  f_1  = pges/(eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp - dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL P_EOS
  f_2  = pges/(eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp + dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL P_EOS
  f_3  = pges/(eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp + 2.0*dist
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL P_EOS
  f_4  = pges/(eta*rho_0*(kbol*t)/1.E-30)

  t  = t_tmp
  CALL PERTURBATION_PARAMETER
  eta = z3t * rho_0
  CALL P_EOS

  dfdr    = pges / ( eta * rho_0 * (kbol*t) / 1.E-30 )
  ddfdrdr = pgesdz / ( eta * rho_0*(kbol*t)/1.E-30 ) - dfdr * 2.0 / eta - 1.0 / eta**2 

  df_drdt   = ( -f_4 + 8.0*f_3 - 8.0*f_2 + f_1) / ( 12.0 * dist )

  cp_res = cv_res - RGAS + RGAS *(zges + eta*t*df_drdt)**2 / (1.0 +2.0*eta*dfdr +eta**2 *ddfdrdr)

  ! write (*,*) 'n',df_dt,df_dtdt
  ! write (*,*) 'kreuz ', zges,eta*t*df_drdt,eta*dfdr, eta**2 *ddfdrdr
  ! write (*,*) 'h, cv', h_res, cv_res
  ! write (*,*) h_res - t*s_res
  ! write (*,*) cv_res,cp_res,eta
  ! pause

END SUBROUTINE H_EOS_num

