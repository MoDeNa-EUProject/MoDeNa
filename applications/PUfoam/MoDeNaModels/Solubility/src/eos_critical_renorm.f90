!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_EOS_RN
  USE BASIC_VARIABLES, ONLY: ncomp, RGT_variant
  IMPLICIT NONE

  IF (RGT_variant == 'phase_cell' .AND. ncomp >= 2) THEN
     CALL F_EOS_RN_PHASE_SPACE_CELL
  ELSE
     CALL F_EOS_RN_ISOMORPHIC_DENSITY
  END IF

END SUBROUTINE F_EOS_RN



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_CRITICAL_RENORM
  USE BASIC_VARIABLES, ONLY: ncomp, RGT_variant
  IMPLICIT NONE

  IF (RGT_variant == 'phase_cell' .AND. ncomp >= 2) THEN
     CALL PHI_CRITICAL_RENORM_PHASE_SPACE_CELL
  ELSE
     CALL PHI_CRITICAL_RENORM_ISOMORPHIC_DENSITY
  END IF

END SUBROUTINE PHI_CRITICAL_RENORM





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_EOS_RN_ISOMORPHIC_DENSITY
!
! White's recursive procedure is:
!
!  f = fid + f0 + delta_f1 + delta_f2 + delta_f3 + ...
!
! where f is the Helmholtz energy per unit volume, and where f0 is the
! classical residual repulsive part, f0 = f_res,rep.
!
! PC-SAFT can be considered for the classical residual part, with
!
!  f0 = f_PCSAFT
!
! It is
!
!                 INT_0^rho( exp(-Kn^-1*(fns(rho+x)-2fns(rho)+fns(rho-x))) )
! delta_fn= -Kn*ln----------------------------------------------------------
!                 INT_0^rho( exp(-Kn^-1*(fnl(rho+x)-2fnl(rho)+fnl(rho-x))) )
!
! where x in the integrals INT runs from x=0 to x=rho_x=rho. The counter
! kk in belows code takes care of the integrals INT from 0 to rho.
! An outer loop is established (with index k) in order to calculate f for
! the whole density range. First, the f0 is calculated and approximated
! by a cubic spline. Then delta_fn is calculated for 0<=rho<=rhomax. New
! spline parameters are determined after evey iteration n.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_EOS_RN_ISOMORPHIC_DENSITY

  USE PARAMETERS
  USE EOS_VARIABLES
  USE EOS_CONSTANTS
  USE EOS_NUMERICAL, only: F_NUMERICAL
  USE FITTING_RGT_PARAMETERS
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m, n, kk, k05max
  INTEGER                                :: stepno, niter, pos, step2
  REAL                                   :: rho0, rhomax, kn(10), LL
  REAL                                   :: drho, densav1
  REAL                                   :: alpha_l, alpha_s
  REAL                                   :: integrand_l, integrand_s
  REAL                                   :: integr_old_l, integr_old_s
  REAL                                   :: integral_l, integral_s
  REAL                                   :: gn_l, gn_s
  REAL                                   :: I1, w_ratio, del_f(8,0:10000)
  REAL                                   :: phi_crit
  REAL                                   :: alph(0:8,0:10000)
  REAL                                   :: rhovec(0:10000), xdev, dzp1, fres_v
  REAL                                   :: sig_m2
  REAL                                   :: m_mean
  REAL                                   :: chapm, order_chap
  REAL                                   :: fid, rho_i(nc)
  REAL, DIMENSION(5000)                  :: utri
  REAL, DIMENSION(0:5000)                :: alphn, beta, gamma, delta
  REAL                                   :: rho_c, rho_l, rho_r
  REAL                                   :: alp_c, alph_lc, alph_sc
  REAL                                   :: alp_r, alph_lr, alph_sr
  REAL                                   :: alp_l, alph_ll, alph_sl
  REAL                                   :: dzp2, delta_rho

  INTEGER, SAVE                          :: scan = 0
  REAL, SAVE                             :: tempsav = 0.0
  REAL, DIMENSION(5), SAVE               :: xsav = 0.0
  REAL, DIMENSION(0:10000), SAVE         :: alphnsav, betasav, gammasav, deltasav, rhovecsav
  !SAVE                                      alphnsav, betasav, gammasav, deltasav, rhovecsav

  !SAVE tempsav,xsav,scan
  !DATA scan /0/
  !DATA tempsav /0.0/
  !DATA xsav    /0.0,0.0,0.0,0.0,0.0/
  !-----------------------------------------------------------------------------

  IF ( critfit == 1 ) THEN
     IF (t < 0.2*parame(1,3) ) t = 0.2*parame(1,3)
     LLi(1)       = LLfit*sig_ij(1,1)
     phi_criti(1) = phifit
     chap(1)      = chapfit
  END IF

  CALL PERTURBATION_PARAMETER

  LL = 0.0
  phi_crit = 0.0
  sig_m2   = 0.0
  m_mean   = 0.0
  chapm    = 0.0
  order_chap = 0.0
  DO i = 1,ncomp
     IF (LLi(i) == 0.0 .OR. phi_criti(i) == 0.0 .OR. chap(i) == 0.0) STOP
     LL = LL + x(i)*LLi(i)**3
     phi_crit = phi_crit + x(i)*phi_criti(i)
     sig_m2 = sig_m2 + x(i)*mseg(i)*sig_ij(i,i)**2
     m_mean = m_mean + x(i)*mseg(i)
     chapm = chapm + x(i)*chap(i)
     DO j = 1,ncomp
        order_chap = order_chap + x(i)*x(j)* mseg(i)*mseg(j)  &
             *sig_ij(i,j)**3  * uij(i,j)/t *(chap(i)*chap(j))**0.5
     END DO
  END DO
  LL = LL**(1.0/3.0)
  sig_m2 = sig_m2/m_mean

  stepno = 400
  niter  = 5
  w_ratio = 2.0

  rho0 = eta/z3t
  rhomax = 0.0
  DO i = 1,ncomp
     ! rhomax = rhomax + x(i)*mseg(i)*sig_ij(i,i)**3
     rhomax = rhomax + x(i)*mseg(i)*dhs(i)**3   !jg
  END DO
  rhomax = SQRT(2.0)/rhomax


  densav1 = eta

  xdev = SUM( ABS( xsav(1:ncomp)-x(1:ncomp) ) )
  IF (tempsav == t .AND. xdev <= 1.E-10 .AND. scan == 1) GO TO 5

  tempsav =  t
  xsav(1:ncomp) = x(1:ncomp)
  !-----------------------------------------------------------------------------



  dzp1 = rhomax / REAL(stepno)
  alph(0,0) = 0.0
  rhovec(0) = 0.0
  DO k = 1, stepno
     drho =  REAL(k)*dzp1
     rhovec(k) = drho
     alph(0,k) = 0.0
     eta = drho*z3t
     fid =  0.0
     rho_i  = 0.0
     DO i = 1,ncomp
        ! debroglie(i) =  6.62606896E-34 *1.E10   &     ! in units Angstrom
        !                 *SQRT( 1.0 / (2.0*PI *1.0 /6.022045E23/1000.0*KBOL*T) )
        rho_i(i) = x(i)*drho
        IF (rho_i(i) > 0.0) fid = fid + x(i)*(LOG(rho_i(i))-1.0)
     END DO
     CALL F_NUMERICAL
     alph(0,k) = (fres+fid)*drho
     alphn(k)  = alph(0,k)
  END DO
  ! CALL SPLINE_PARA (dzp1,alphn,utri,stepno)
  ! CALL SPLINE_COEFF(beta,gamma,delta,dzp1,alphn,utri,stepno)
  scan = 1


  ! alpha_l =  2.0/3.0*PI *1.5**3 *order1
  ! alpha_s =  alpha_l*phi_crit* 0.2 *1.5**2 *sig_ij(1,1)**2/LL/LL


  DO n = 1,niter

     k05max = stepno/2
     kn(n) = 1.0 / LL**3 / w_ratio**REAL(3*n)
     alpha_l =  16.0/9.0 * PI*order_chap !* chapm
     alpha_s =  alpha_l*phi_crit*  9.0/7.0 * sig_m2  /LL/LL /w_ratio**REAL(2*n)/2.0

     ! DO k = 1,stepno-1        !!!!!
     DO k = 1,stepno*3/4-1       !!!!!
        integr_old_l = 1.0
        integr_old_s = 1.0
        ! hx(1) = 0.0
        ! hyl(1) = integr_old_l
        ! hys(1) = integr_old_s
        ! nn = 1
        integral_l   = 0.0
        integral_s   = 0.0
        step2 = k                        !step2: step# for integrat.(index kk)
        IF (k > k05max) step2 = stepno-k
        ! IF (k.LE.20)step2 = 2*step2
        ! step2 = 2*step2
        ! IF (step2.GT.stepno) step2 = stepno
        dzp2 = rhovec(k)/REAL(step2)
        IF (rhovec(k) > rhomax/2.0) dzp2 = (rhomax-rhovec(k)) / REAL(step2)
        DO kk = 1, step2
           delta_rho = REAL(kk)*dzp2
           rho_c = rhovec(k)
           rho_l = rhovec(k) - delta_rho
           rho_r = rhovec(k) + delta_rho


           eta = rho_c*PI/6.0*mseg(1)*dhs(1)**3
           I1 = 0.0
           DO m = 0,6
              I1 = I1 + apar(m)*eta**REAL(m)
           END DO
           ! alpha_l = + 2.0*PI*I1*order_chap
           ! alpha_s =  alpha_l*phi_crit*  9.0/7.0*sig_m2  /LL/LL  / w_ratio**REAL(2*n)/2.0
           ! pos = INT(rho_c/rhomax*stepno + 1.E-9)
           ! alp_c = alphn(pos) + beta(pos) *(rho_c-rhovec(pos))
           !                    + gamma(pos)*(rho_c-rhovec(pos))**2
           !                    + delta(pos)*(rho_c-rhovec(pos))**3
           alp_c = alph(n-1,k)
           alph_lc  = alp_c +alpha_l*rho_c*rho_c
           alph_sc  = alp_c +alpha_s*rho_c*rho_c

           eta = rho_l*PI/6.0*mseg(1)*dhs(1)**3
           I1 = 0.0
           DO m = 0,6
              I1 = I1 + apar(m)*eta**REAL(m)
           END DO
           ! alpha_l = + 2.0*PI*I1*order_chap
           ! alpha_s =  alpha_l*phi_crit*  9.0/7.0*sig_m2  /LL/LL / w_ratio**REAL(2*n)/2.0
           ! pos = INT(rho_l/rhomax*stepno + 1.E-9)
           ! alp_l = alphn(pos) + beta(pos) *(rho_l-rhovec(pos))
           !                    + gamma(pos)*(rho_l-rhovec(pos))**2
           !                    + delta(pos)*(rho_l-rhovec(pos))**3
           alp_l = alph(n-1,k-kk)
           alph_ll  = alp_l +alpha_l*rho_l*rho_l
           alph_sl  = alp_l +alpha_s*rho_l*rho_l

           eta = rho_r*PI/6.0*mseg(1)*dhs(1)**3
           I1 = 0.0
           DO m = 0,6
              I1 = I1 + apar(m)*eta**REAL(m)
           END DO
           ! alpha_l = + 2.0*PI*I1*order_chap
           ! alpha_s =  alpha_l*phi_crit*  9.0/7.0*sig_m2  /LL/LL / w_ratio**REAL(2*n)/2.0
           ! pos = INT(rho_r/rhomax*stepno + 1.E-9)
           ! alp_r = alphn(pos) + beta(pos) *(rho_r-rhovec(pos))
           !                    + gamma(pos)*(rho_r-rhovec(pos))**2
           !                    + delta(pos)*(rho_r-rhovec(pos))**3
           alp_r = alph(n-1,k+kk)
           alph_lr  = alp_r +alpha_l*rho_r*rho_r
           alph_sr  = alp_r +alpha_s*rho_r*rho_r

           gn_l =  (alph_lr+alph_LL)/2.0 - alph_lc
           gn_s =  (alph_sr+alph_sl)/2.0 - alph_sc
           IF (gn_l < 0.0) gn_l = 0.0
           IF (gn_s < 0.0) gn_s = 0.0
           ! write (*,*) k,kk,Gn_s
           integrand_l = 0.0
           integrand_s = 0.0
           IF ( -gn_l/kn(n) > -300.0 .AND. -gn_l/kn(n) < 300.0) integrand_l = EXP( -gn_l/kn(n) )
           IF ( -gn_s/kn(n) > -300.0 .AND. -gn_s/kn(n) < 300.0) integrand_s = EXP( -gn_s/kn(n) )
           ! nn = nn+1
           ! hx(nn) = delta_rho
           ! hyl(nn) = integrand_l
           ! hys(nn) = integrand_s
           integral_l = integral_l + dzp2 *(integrand_l+integr_old_l)/2.0
           integral_s = integral_s + dzp2 *(integrand_s+integr_old_s)/2.0
           integr_old_l = integrand_l
           integr_old_s = integrand_s

           IF (integrand_l < 1.E-9 .AND. integrand_s < 1.E-9) THEN
              GO TO 15  ! end loop, because no significant contribution is expected
           END IF

        END DO   ! enddo kk (integration over density)
        ! CALL spline(hx,hyl,nn,0.0,0.0,y2)
        ! CALL splint_integral(hx,hyl,y2,nn,0.0,delta_rho,integral_l)
        ! CALL spline(hx,hys,nn,0.0,0.0,y2)
        ! CALL splint_integral(hx,hys,y2,nn,0.0,delta_rho,integral_s)
15      CONTINUE

        del_f(n,k) = 0.0
        IF (integral_s /= 0.0 .AND.  integral_l /= 0.0) THEN
           del_f(n,k) = -kn(n)*LOG(integral_s/integral_l)
        END IF
     END DO   ! enddo k

     del_f(n,0) = 0.0
     DO k = 0,stepno
        alph(n,k) =  alph(n-1,k) + del_f(n,k)
        alphn(k)  =  alph(n,k)
        ! subract the ideal gas part in order to have less non-linearity in the low
        ! density region. The ideal gas part is then added to the final value below.
        IF(n == niter .AND. k /= 0) THEN
           fid =  0.0
           DO i = 1,ncomp
              rho_i(i) = x(i)*rhovec(k)
              IF(rho_i(i) > 0.0) fid = fid +rho_i(i)*(LOG(rho_i(i))-1.0)
           END DO
           alphn(k) = alphn(k)-fid
        END IF
     END DO
     CALL spline_para (dzp1,alphn,utri,stepno)
     CALL spline_coeff(beta,gamma,delta,dzp1,alphn,utri,stepno)

     DO k = 0,stepno
        alphnsav(k) = alphn(k)
        betasav(k) = beta(k)
        gammasav(k) = gamma(k)
        deltasav(k) = delta(k)
        rhovecsav(k)  = rhovec(k)
     END DO

  END DO   ! loop of n cycles
  ! write (71,*) ' '


5 CONTINUE
  ! CALL SPLINE_PARA (dzp1,alphn,utri,stepno)
  ! CALL SPLINE_COEFF(beta,gamma,delta,dzp1,alphn,utri,stepno)

  DO k = 0,stepno
     alphn(k)  = alphnsav(k)
     beta(k)   = betasav(k)
     gamma(k)  = gammasav(k)
     delta(k)  = deltasav(k)
     rhovec(k) = rhovecsav(k)
  END DO
  ! write (*,*) 'beta',beta(70),gamma(70),t



  ! fid =  0.0
  ! DO i = 1,ncomp
  !  ! in units Angstrom
  !  debroglie(i) =  6.62606896d-34 *1d10 *SQRT( 1.0 / (2.0*PI *1.0 /6.022045d23/1000.0*KBOL*T) )
  !  rho_i(i) = x(i)*rho0
  !  fid =  fid + x(i)*(LOG(rho_i(i))-1.0)
  ! ENDDO
  ! !fid =  0.0


  pos = INT(rho0/rhomax*REAL(stepno))
  fres_v = alphn(pos) + beta(pos) *(rho0-rhovec(pos))  &
       + gamma(pos)*(rho0-rhovec(pos))**2  + delta(pos)*(rho0-rhovec(pos))**3
  fres = fres_v/rho0 !-fid

  pges = ( (beta(pos)+ 2.0*gamma(pos)*(rho0-rhovec(pos))  &
       + 3.0* delta(pos)*(rho0-rhovec(pos))**2  )*rho0 -fres_v )  &
       *(KBOL*t)/1.E-30 +   rho0 * (KBOL*t) / 1.E-30   ! ideal gas contribution
  ! if(pos.GT.1)write (69,*) pos,delta(pos)
  !     &          (beta(pos)+ 2.0*gamma(pos)*(rho0-rhovec(pos))
  !     &        + 3.0* delta(pos)*(rho0-rhovec(pos))**2 )
  pgesdz = ( 2.0*gamma(pos)*rho0  &
       + rho0*6.0* delta(pos)*(rho0-rhovec(pos)) ) *(KBOL*t)/1.E-30 /z3t  &
       + 1.0/z3t*(KBOL*t)/1.E-30   ! ideal gas contribution
  pgesd2 = ( 2.0*gamma(pos) + 6.0* delta(pos)*(2.0*rho0-rhovec(pos)) ) *(KBOL*t)/1.E-30 /z3t /z3t
  pgesd3 = ( 12.0* delta(pos) ) *(KBOL*t)/1.E-30 /z3t /z3t /z3t


  eta = densav1

END SUBROUTINE F_EOS_RN_ISOMORPHIC_DENSITY



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_CRITICAL_RENORM_ISOMORPHIC_DENSITY

  USE parameters, ONLY: nc, nsite, PI, KBOL
  USE EOS_CONSTANTS
  USE EOS_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, k
  REAL                                   :: myresq
  REAL                                   :: zres, zges, sumseg
  REAL                                   :: fres1, fres2 = 0.0, fres3, fres4
  REAL                                   :: xconc, d_org, d_new
  REAL                                   :: dicht, gradi
  REAL                                   :: myres(nc), mpart(nc)
  REAL                                   :: gres !,my_pt1, p_pt1
  REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! density iteration or pressure calculation
  !-----------------------------------------------------------------------------
  IF (ensemble_flag == 'tp') THEN
     CALL DENSITY_ITERATION
  ELSEIF (ensemble_flag == 'tv') THEN
     eta = eta_start
     rho = eta/z3t
     sumseg = z0t/(PI/6.0)   ! mean segment number
     CALL F_EOS_RN
  ELSE
     write (*,*) 'PHI_CRITICAL_RENORM: define ensemble, ensemble_flag == (pv) or (pt)'
     stop
  END IF

!!$!---------------density iteration-------------------------------------
!!$IF (dfti >= 98) THEN
!!$  eta = eta_start
!!$  rho = eta/z3t
!!$  sumseg = z0t/(PI/6.0)   ! mean segment number
!!$  CALL F_EOS_RN
!!$ELSE
!!$  CALL densitr_rn
!!$END IF
!!$
  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------
  zges = (p * 1.E-30) / (KBOL*t*eta/z3t)
  IF ( ensemble_flag == 'tv' ) zges = (pges * 1.E-30) / (kbol*t*eta/z3t)
  zres = zges - 1.0
  !z_ges = zges


  IF (ncomp == 1) THEN
     gres = fres + zres
     IF (ensemble_flag == 'tp') lnphi(1) = gres - LOG(zges)
     IF (ensemble_flag == 'tv' .AND. eta >= 0.0) lnphi(1) = gres !+LOG(rho)
  END IF

!!$!--------------------------
!!$IF (dfti == 99) THEN
!!$   ! z3 = rho*PI/6.0*mseg(1)*dhs(1)**3
!!$   CALL phipt1 ( my_pt1, p_pt1 )
!!$   lnphi(1) = lnphi(1)+my_pt1
!!$   z_ges = z_ges+(p_pt1 * 1.E-30)/(KBOL*t*rho)
!!$END IF
!!$!--------------------------


  IF (ncomp /= 1) THEN
     !----------------------derivations with big D--------------------------------
     !      gradi = 1.E-14
     gradi = 1.E-13

     DO  k = 1,ncomp
        xconc = x(k)
        dicht = eta
        d_org = z3t/PI*6.0

        IF (xconc > 1.E-13 .AND. xconc < 0.9999) THEN
           x(k) =  xconc-gradi*10.0**(0.5*(15.0+LOG10(xconc)))
        ELSE
           IF (xconc < 0.9999) THEN
              x(k) = xconc + 2.E-11
           ELSE
              x(k) = xconc - 2.E-10
           END IF
        END IF
        d_new = 0.0
        DO  i = 1,ncomp
           d_new =  d_new +x(i)*mseg(i)*dhs(i)**3
        END DO
        eta = dicht*d_new/d_org
        CALL F_EOS_RN
        fres1 = fres
        tfr_1 = tfr

        IF (xconc > 1.E-13 .AND. xconc < 0.9999) THEN
           !****************************************
           x(k)  = xconc+gradi*10.0**(0.5*(15.0+LOG10(xconc)))
           d_new = 0.0
           DO  i = 1,ncomp
              d_new =  d_new +x(i)*mseg(i)*dhs(i)**3
           END DO
           eta = dicht*d_new/d_org
           CALL F_EOS_RN
           fres2 = fres
           tfr_2 = tfr
           !****************************************
        END IF

        IF (xconc > 1.E-13 .AND. xconc < 0.9999) THEN
           x(k) =  xconc-0.5*gradi*10.0**(0.5*(15.0+LOG10(xconc)))
        ELSE
           x(k) = xconc
        END IF
        d_new = 0.0
        DO  i = 1,ncomp
           d_new =  d_new +x(i)*mseg(i)*dhs(i)**3
        END DO
        eta = dicht*d_new/d_org
        CALL F_EOS_RN
        fres3 = fres
        tfr_3 = tfr

        IF (xconc > 1.E-13 .AND. xconc < 0.9999) THEN
           x(k) =  xconc+0.5*gradi*10.0**(0.5*(15.0+LOG10(xconc)))
        ELSE
           IF (xconc < 0.9999) THEN
              x(k) = xconc + 1.E-11
           ELSE
              x(k) = xconc - 1.E-10
           END IF
        END IF
        d_new = 0.0
        DO  i = 1,ncomp
           d_new =  d_new +x(i)*mseg(i)*dhs(i)**3
        END DO
        eta = dicht*d_new/d_org
        CALL F_EOS_RN
        fres4 = fres
        tfr_4 = tfr

        x(k)  = xconc
        eta = dicht


        IF (xconc > 1.E-13 .AND. xconc < 0.9999) THEN
           mpart(k) = (fres1-8.0*fres3+8.0*fres4-fres2)  &
                /(6.0*gradi*10.0**(0.5*(15.0+LOG10(xconc))))
           !          tfrdx(k) = (tfr_1-8.0*tfr_3+8.0*tfr_4-tfr_2)
           !     &             /(6.0*gradi*10.0**(0.5d0*(15.0+LOG10(xconc))))
           !          write (*,*) tfrdx,eta
        ELSE
           IF (xconc < 0.9999) THEN
              mpart(k) = (-3.0*fres3+4.0*fres4-1.0*fres1)/2.E-11
              !          MPART(k) = (fres4-fres3)/1.E-12
           ELSE
              mpart(k) =  - (-4.0*fres3+6.0*fres4-2.0*fres1)/2.E-10
              !          MPART(k) = (fres3-fres4)/1.E-5
           END IF
        END IF

     END DO

     !--------------------------------------------------------------------------
     ! f res
     !--------------------------------------------------------------------------
     CALL F_EOS_RN


     myresq = 0.0
     DO   i = 1, ncomp
        myresq = myresq - x(i)* mpart(i)
     END DO

     DO  k = 1, ncomp
        myres(k) = myresq + mpart(k)  + fres + zres
        IF (ensemble_flag == 'tp') lnphi(k) = myres(k) - LOG(zges)
        IF (ensemble_flag == 'tv' .AND. eta >= 0.0) lnphi(1) = myres(k) !+LOG(rho)
        ! IF (DFT.GE.98) write (*,*) dft
        ! write (*,*) 'lnphi',k,LNPHI(k),x(k),MYRES(k), -LOG(ZGES)
     END DO
     ! write(*,*) t,p,eta,x(1)
     ! stop

  END IF    !  ncomp.NE.1

END SUBROUTINE PHI_CRITICAL_RENORM_ISOMORPHIC_DENSITY



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_EOS_RN

  USE EOS_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL :: z0,z1,z2,z3
  REAL :: dzetdv,dicht,dist,fact
  REAL :: pgesdz1,pgesdz4,pgesdz5
  REAL :: pges2d1,pges2d4,pges2d5
  !-----------------------------------------------------------------------


  rho = eta/z3t
  z0 = z0t*rho
  z1 = z1t*rho
  z2 = z2t*rho
  z3 = z3t*rho

  !-----------------density-------------------------------------------
  dzetdv = eta*rho

  IF (eta > 0.1) THEN
     fact = 1.0
  ELSE IF (eta <= 0.1 .AND. eta > 0.01) THEN
     fact = 10.0
  ELSE
     fact = 50.0
  END IF
  dist = eta*1.E-2 *fact

  dicht  = eta
  eta  = dicht - 2.0*dist
  CALL F_EOS_RN
  pgesdz1  = pgesdz
  pges2d1  = pgesd2
  !      eta  = dicht - dist
  !      CALL F_EOS_RN
  !      pgesdz2  = pgesdz
  !      pges2d2  = pgesd2
  !      eta  = dicht + dist
  !      CALL F_EOS_RN
  !      pgesdz3  = pgesdz
  !      pges2d3  = pgesd2
  eta  = dicht + 2.0*dist
  CALL F_EOS_RN
  pgesdz4  = pgesdz
  pges2d4  = pgesd2
  eta  = dicht
  CALL F_EOS_RN
  pgesdz5  = pgesdz
  pges2d5  = pgesd2

  !      pgesd3 = (-pgesdz4+16.0*pgesdz3-30.0*pgesdz5+16.0*pgesdz2-pgesdz1)
  !     &              /(12.0*(dist**2 ))
  !      write(*,*) pgesd3
  !      pgesd3 = (pges2d4-pges2d1)/(4.0*dist)
  !      write(*,*) pgesd3
  pgesd3 = (pgesdz4-2.0*pgesdz5+pgesdz1)/(2.0*dist)**2

END SUBROUTINE P_EOS_RN


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE H_EOS_RN

  USE PARAMETERS, ONLY: RGAS
  USE EOS_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL :: dist,fact
  REAL :: fres1,fres2,fres3,fres4,fres5
  REAL :: cv_res,t_tmp,zges,df_dt, df_dtdt
  !-----------------------------------------------------------------------


  CALL PERTURBATION_PARAMETER
  rho = eta/z3t


  fact = 1.0
  dist = t * 100.E-5 * fact

  t_tmp  = t
  t  = t_tmp - 2.0*dist
  CALL F_EOS_RN
  fres1  = fres
  t  = t_tmp - dist
  CALL F_EOS_RN
  fres2  = fres
  t  = t_tmp + dist
  CALL F_EOS_RN
  fres3  = fres
  t  = t_tmp + 2.0*dist
  CALL F_EOS_RN
  fres4  = fres
  t  = t_tmp
  CALL F_EOS_RN
  fres5  = fres
  ! *(KBOL*T)/1.E-30

  zges = (p * 1.E-30) / (KBOL*t*rho)


  df_dt   = (-fres4+8.0*fres3-8.0*fres2+fres1)/(12.0*dist)
  df_dtdt = (-fres4+16.0*fres3-30.0*fres5+16.0*fres2-fres1) / (12.0*(dist**2 ))


  s_res  =  (- df_dt -fres/t)*RGAS*t + RGAS *LOG(zges)
  h_res  =  ( - t*df_dt  +  zges-1.0  ) * RGAS*t
  cv_res = -(   t*df_dtdt + 2.0*df_dt  ) * RGAS*t

END SUBROUTINE H_EOS_RN
