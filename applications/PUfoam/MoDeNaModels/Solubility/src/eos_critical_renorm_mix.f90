Module RGT_PHASE_SPACE_CELL

  use PARAMETERS, only: nc
  implicit none
  save

  real, dimension(nc)                   :: dfdx
  REAL, DIMENSION(300)                  :: x1a
  REAL, DIMENSION(300)                  :: x2a
  REAL, DIMENSION(300,300)              :: ya, y1a, y2a, y12a
  REAL, DIMENSION(300,300,4,4)          :: c_bicub

End Module RGT_PHASE_SPACE_CELL

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_EOS_RN_PHASE_SPACE_CELL
!
! White's recursive procedure is:
!
!  f = fid + f0 + delta_f1 + delta_f2 + delta_f3 + ...
!
! where f is the Helmholtz energy per unit volume, and where f0 is the
! classical residual repulsive part, f0=f_res,rep.
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
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_EOS_RN_PHASE_SPACE_CELL

  USE PARAMETERS
  USE EOS_VARIABLES
  USE EOS_CONSTANTS
  USE EOS_NUMERICAL, only: F_NUMERICAL
  USE FITTING_RGT_PARAMETERS
  USE RGT_PHASE_SPACE_CELL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: n, i, j, m
  INTEGER                                :: k1, k2, kk1, kk2
  INTEGER                                :: k05max
  INTEGER                                :: stepno, niter, pos1, pos2, step1, step2
  REAL                                   :: rho0, rhomax1, rhomax2, kn(10,0:300,0:300), LL
  REAL                                   :: xsav1, xsav2
  REAL                                   :: int1_l, int1_s, int3_l, int3_s
  REAL                                   :: int4_l, int4_s
  REAL                                   :: integral_l, integral_s
  REAL                                   :: del_f(8,0:300,0:300)
  REAL                                   :: phi_crit, rhoi(nc)
  REAL                                   :: alph(0:8,0:300,0:300), w_ratio
  REAL                                   :: rhovec1(0:300), rhovec2(0:300), dzp1, dzp2
  REAL                                   :: fres_v, sig_m2, m_mean
  REAL                                   :: chapm, order_chap, fid
  REAL                                   :: mfrac(2)
  REAL                                   :: rhot, alphn(0:300,0:300), combi(0:300,0:300)
  REAL                                   :: combi2(10,0:300,0:300), fact, I1
  REAL                                   :: rho1, rho2, fdr1, fdr2, fdr11, fdr12, fdr22
  REAL                                   :: fdr111, fdr112, fdr122, fdr222, f_drho, f_drho2
  REAL                                   :: f_drho3, f_drho4, z3tsav
  REAL                                   :: int3_lx, int3_sx, int_o_l, int_o_s
  REAL                                   :: densav1

  INTEGER, SAVE                          :: scan = 0
  REAL, SAVE                             :: tempsav = 0.0
  !-----------------------------------------------------------------------------

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
             *sig_ij(i,j)**3 * uij(i,j)/t *(chap(i)*chap(j))**0.5
     END DO
  END DO
  LL = LL**(1.0/3.0)
  sig_m2 = sig_m2 / m_mean

  stepno = 200     ! for stepno > 200, expand the array dimension of DFT_MODULE
  niter  = 4
  w_ratio = 2.0

  xsav1 = x(1)
  xsav2 = x(2)
  z3tsav = z3t

  rho0 = eta/z3t
  rho1 = x(1)*rho0
  rho2 = x(2)*rho0
  ! rhomax = 0.0
  ! DO i = 1,ncomp
  !   !rhomax = rhomax + x(i)*mseg(i)*sig_ij(i,i)**3
  !   rhomax = rhomax + x(i)*mseg(i)*dhs(i)**3  !jg
  ! ENDDO
  ! rhomax = SQRT(2.0)/rhomax
  rhomax1 =  mseg(1)*dhs(1)**3  !jg
  rhomax2 =  mseg(2)*dhs(2)**3  !jg
  rhomax1 =  SQRT(2.0)/rhomax1
  rhomax2 =  SQRT(2.0)/rhomax2


  densav1 = eta

  IF (tempsav == t .AND. scan == 1) GO TO 5

  tempsav =  t
  !-----------------------------------------------------------------------------


  dzp1 = rhomax1/REAL(stepno)
  dzp2 = rhomax2/REAL(stepno)
  DO k2 = 0,stepno
     DO k1 = 0,stepno
        rhovec1(k1) = REAL(k1)*dzp1
        rhovec2(k2) = REAL(k2)*dzp2
        rhot = rhovec1(k1)+rhovec2(k2)
        x(1) = rhovec1(k1)/rhot
        x(2) = rhovec2(k2)/rhot
        alph(0,k1,k2) = 0.0
        eta = rhot*PI/6.0*(x(1)*mseg(1)*dhs(1)**3  &
             +x(2)*mseg(2)*dhs(2)**3) !*z3t
        IF (eta > 0.7) eta = 0.8-0.1*EXP((0.7-eta)*10.0)
        rhoi(1) = rhovec1(k1)
        rhoi(2) = rhovec2(k2)

        m_mean =  x(1)*mseg(1) + x(2)*mseg(2)
        ! v_mean =  x(1)*mseg(1)*sig_ij(1,1)**3 +x(2)*mseg(2)*sig_ij(2,2)**3
        ! vfrac(1) = mseg(1)*sig_ij(1,1)**3 / v_mean
        ! vfrac(2) = mseg(2)*sig_ij(2,2)**3 / v_mean
        mfrac(1) = mseg(1) / m_mean
        mfrac(2) = mseg(2) / m_mean

        combi(k1,k2) = 0.0
        DO i = 1,ncomp
           DO j = 1,ncomp
              combi(k1,k2) = combi(k1,k2)+16.0/9.0*PI*(  &
                   rhoi(i)*rhoi(j)*mseg(i)*mseg(j)*sig_ij(i,j)**3  &
                   * uij(i,j)/t *(chap(i)+chap(j))*.5 )
              !combi(k1,k2) = combi(k1,k2)+16.0/9.0*PI*(  &
              !                  rhoi(i)*rhoi(j)*mseg(i)*mseg(j)*sig_ij(i,j)**3  &
              !                * uij(i,j)/t *(chap(i)*mseg(i)+chap(j)*mseg(j))*0.5 )  &
              !                / m_mean
           END DO
        END DO

        I1 = 0.0
        DO m = 0,6
           I1  =  I1 + apar(m)*eta**REAL(m)
        END DO
        order1 = 0.0
        DO i = 1,ncomp
           DO j = 1,ncomp
              order1 = order1 + x(i)*x(j)* mseg(i)*mseg(j)  &
                   *sig_ij(i,j)**3 * uij(i,j)/t *(chap(i)+chap(j))*.5
           END DO
        END DO
        ! combi(k1,k2) = + 2.0*PI*rhot*rhot*I1*order1
        ! write (*,*) combi(k1,k2),+ 2.0*PI*rhot*rhot*I1*order1

        DO n = 1,niter
           combi2(n,k1,k2) = 0.0
           DO i = 1,ncomp
              DO j = 1,ncomp
                 combi2(n,k1,k2) = combi2(n,k1,k2)+16.0/9.0*PI*9.0/7.0*(  &
                      rhoi(i)*rhoi(j)*mseg(i)*mseg(j)*sig_ij(i,j)**2  &
                      *(phi_criti(i)*mfrac(i)  &
                      +phi_criti(j)*mfrac(j))*0.5 *sig_ij(i,j)**3  &
                      * uij(i,j)/t *(chap(i)+chap(j))*0.5  &
                      /(mfrac(i)**0.33333333*LLi(i)) /(mfrac(j)**0.33333333*LLi(j))  &
                      )/ ( w_ratio**REAL(2*n) * 2.0 )
                 !         combi2(n,k1,k2) = combi2(n,k1,k2)+16.0/9.0*PI*9.0/7.0*(
                 !     &             rhoi(i)*rhoi(j)*mseg(i)*mseg(j)*sig_ij(i,j)**2
                 !     &            *(phi_criti(i)/mseg(i) + phi_criti(j)/mseg(j))*0.5
                 !     &            *sig_ij(i,j)**3 * uij(i,j)/t *(chap(i)*mseg(i)+chap(j)*mseg(j))*.5
                 !     & /(mfrac(i)**0.33333333*LLi(i))
                 !     & /(mfrac(j)**0.33333333*LLi(j))
                 !     &              )/ ( w_ratio**REAL(2*n) * 2.0 )
              END DO
           END DO
           LL = ( x(1)*mfrac(1)*LLi(1)**3  &
                + x(2)*mfrac(2)*LLi(2)**3 )**(1.0/3.0)
           ! LL = 0.0
           ! DO i = 1,ncomp
           ! DO j = 1,ncomp
           ! LL =  LL +x(i)*x(j)*mfrac(i)*mfrac(j)*LLi(i)**1.5*LLi(j)**1.5
           ! ENDDO
           ! ENDDO
           ! LL = LL**(1.0/3.0)
           fact =  x(1)*sig_ij(1,1)**2 *phi_criti(1)  &
                +x(2)*sig_ij(2,2)**2 *phi_criti(2)
           ! combi2(n,k1,k2) = combi(k1,k2)*9.0/7.0*fact  /LL/LL / ( 2.0**REAL(2*n) * 2.0 )
           kn(n,k1,k2) = 1.0/LL**3 / w_ratio**REAL(3*n)
        END DO
        I1 = 0.0
        DO m = 0,6
           I1 = I1 + apar(m)*eta**REAL(m)
        END DO
        order1 = 0.0
        DO i = 1,ncomp
           DO j = 1,ncomp
              order1 = order1 + x(i)*x(j)* mseg(i)*mseg(j)  &
                   *sig_ij(i,j)**5 * uij(i,j)/t *(chap(i)+chap(j))*.5  &
                   *(phi_criti(i)*phi_criti(j))**0.5 /LLi(i)/LLi(j)
           END DO
        END DO
        DO n = 1,niter
           ! combi2(n,k1,k2) = + 2.0*PI*rhot*rhot*I1*order1  *9.0/7.0 / ( w_ratio**REAL(2*n) * 2.0 )
        END DO


        fid =  0.0
        IF (rhot /= 0.0) THEN
           DO i = 1,ncomp
              ! debroglie(i) =  6.62606896d-34 *1d10  &       ! in units Angstrom
              !                 * SQRT( 1.0 / (2.0*PI *1.0 /6.022045d23/1000.0*KBOL*T) )
              IF (rhoi(i) > 0.0) fid =  fid + rhoi(i)*(LOG(rhoi(i))-1.0)
           END DO
           CALL f_numerical
           alph(0,k1,k2) = fres*rhot + fid
           alphn(k1,k2)  = alph(0,k1,k2)
           ! ya(k1+1,k2+1) = alphn(k1,k2)  ! this line can be deleted after debugging
        END IF
        ! x1a(k1+1) = rhovec1(k1)  ! this line can be deleted after debugging
     END DO
     ! x2a(k2+1) = rhovec2(k2)  ! this line can be deleted after debugging
  END DO
  ! CALL SPLINE_PARA (dzp,alphn,utri,stepno)
  ! CALL SPLINE_COEFF(beta,gamma,delta,dzp,alphn,utri,stepno)
  scan = 1

  x(1) = xsav1
  x(2) = xsav2
  z3t  = z3tsav

  !-----------------------------------------------------------------------------
  !      alpha_l =  16.0/9.0*PI*order_chap !* chapm
  !c      alpha_s =  alpha_l*phi_crit*  9.0/7.0*sig_m2  /LL/LL
  !      alpha_s =  phi_crit*  9.0/7.0*sig_m2  /LL/LL

  !      n = 1   ! for debugging only
  !      GOTO 33
  DO n = 1,niter

     ! Kn(n) = 1.0/LL**3 / 2.0**REAL(3*n)
     k05max = stepno/2


     ! calculate for rhovec2 = 0
     DO k1 = 1,stepno*3/4-1       !!!!!
        int_o_l = 1.0
        int_o_s = 1.0
        integral_l   = 0.0
        integral_s   = 0.0
        step1 = k1                        !step1: step# for integrat.(index kk1)
        IF(k1 > k05max) step1 = stepno-k1
        DO kk1 = 0,step1-1
           CALL integr(n,k1,0,kk1,0,kn,combi,combi2, alph,int1_l,int1_s)
           ! int1_lv(kk1) = int1_l
           ! int1_sv(kk1) = int1_s
           ! integral_l = integral_l + dzp1 *(int1_l+int_o_l)/2.0
           ! integral_s = integral_s + dzp1 *(int1_s+int_o_s)/2.0
           integral_l = integral_l + dzp1 *int1_l
           integral_s = integral_s + dzp1 *int1_s
           int_o_l = int1_l
           int_o_s = int1_s
        END DO
        del_f(n,k1,0) = 0.0
        IF (integral_s /= 0.0.AND. integral_l /= 0.0) THEN
           del_f(n,k1,0) = -kn(n,k1,0)*LOG(integral_s/integral_l)
        END IF
     END DO   ! enddo k1

     ! calculate for rhovec1 = 0
     DO k2 = 1,stepno*3/4-1       !!!!!
        int_o_l = 1.0
        int_o_s = 1.0
        integral_l   = 0.0
        integral_s   = 0.0
        step2 = k2                        !step2: step# for integrat.(index kk2)
        IF(k2 > k05max) step2 = stepno-k2
        DO kk2 = 0,step2-1
           CALL integr(n,0,k2,0,kk2,kn,combi,combi2, alph,int3_l,int3_s)
           ! int3_lv(kk2) = int3_l
           ! int3_sv(kk2) = int3_s
           ! integral_l = integral_l + dzp2 *(int3_l+int_o_l)/2.0
           ! integral_s = integral_s + dzp2 *(int3_s+int_o_s)/2.0
           integral_l = integral_l + dzp2 *int3_l
           integral_s = integral_s + dzp2 *int3_s
           int_o_l = int3_l
           int_o_s = int3_s
        END DO
        del_f(n,0,k2) = 0.0
        IF (integral_s /= 0.0.AND. integral_l /= 0.0) THEN
           del_f(n,0,k2) = -kn(n,0,k2)*LOG(integral_s/integral_l)
        END IF
     END DO   ! enddo k2



     ! DO k = 1,stepno-1        !!!!!
     DO k2 = 0,stepno*3/4-1       !!!!!
        WRITE (*,*) k2,n
        DO k1 = 0,stepno*3/4-1       !!!!!
           int3_lx = 1.0
           int3_sx = 1.0
           ! hx(1) = 0.0
           ! hyl(1) = int3_l
           ! hys(1) = int3_s
           ! nn = 1
           integral_l   = 0.0
           integral_s   = 0.0
           step1 = k1                        !step1: step# for integrat.(index kk1)
           step2 = k2                        !step2: step# for integrat.(index kk2)
           IF(k1 > k05max) step1 = stepno-k1
           IF(k2 > k05max) step2 = stepno-k2

           DO kk2 = 0,step2-1
              DO kk1 = 0,step1-1

                 ! write (*,*)rhovec1(k1)+rhovec2(k2),k1,k2
                 ! CALL INTEGR(n,k1,k2,kk1-1,kk2-1,Kn,combi,combi2,alph,int1_l,int1_s)
                 ! CALL INTEGR(n,k1,k2,kk1,kk2-1,Kn,combi,combi2,alph,int2_l,int2_s)
                 ! CALL INTEGR(n,k1,k2,kk1-1,kk2,Kn,combi,combi2,alph,int3_l,int3_s)
                 CALL integr(n,k1,k2,kk1,kk2,kn,combi,combi2, alph,int4_l,int4_s)
                 ! if (k2.EQ.5.AND.kk1.eq.1)write (*,*) int3_l,int4_l,k1,k2
                 ! if (k2.EQ.5.AND.kk1.eq.1)write (*,*) int1_l,int2_l,kk1,kk2
                 ! write (*,*) 'pp',int3_l,int3_lx,k1,k2,kk1,kk2
                 ! if (k2.EQ.5.AND.kk1.eq.1) read (*,*)

                 ! integral_l = integral_l+dzp1*dzp2*(int4_l+int3_l+int2_l+int1_l)/4.0
                 ! integral_s = integral_s+dzp1*dzp2*(int4_s+int3_s+int2_l+int1_l)/4.0
                 integral_l = integral_l+dzp1*dzp2*int4_l
                 integral_s = integral_s+dzp1*dzp2*int4_s
                 int3_lx = int4_l
                 int3_sx = int4_s

                 IF (int4_l < 1.E-9.AND.int4_s < 1.E-9) THEN
                    ! write (*,*) kk1,kk2,int4_l,int4_s
                    GO TO 15  ! end loop, because no significant contribution is expected
                 END IF

              END DO   ! enddo kk1 (integration over density)
15            CONTINUE
           END DO   ! enddo kk2 (integration over density)


           IF (k1 /= 0.AND.k2 /= 0) del_f(n,k1,k2) = 0.0
           IF (integral_s /= 0.0.AND. integral_l /= 0.0) THEN
              del_f(n,k1,k2) = -kn(n,k1,k2)*LOG(integral_s/integral_l)
           END IF

        END DO   ! enddo k1
     END DO   ! enddo k2

     ! 33   CONTINUE  ! for debugging only

     del_f(n,0,0) = 0.0
     DO k2 = 0,stepno
        DO k1 = 0,stepno
           alph(n,k1,k2) =  alph(n-1,k1,k2)+del_f(n,k1,k2)
           alphn(k1,k2) = alph(n,k1,k2)
           ! subract the ideal gas part in order to have less non-linearity in the low
           ! density region. The ideal gas part is then added to the final value below.
           IF(n == niter) THEN
              fid =  0.0
              rhoi(1) = rhovec1(k1)
              rhoi(2) = rhovec2(k2)
              DO i = 1,ncomp
                 IF (rhoi(i) > 0.0) fid =  fid +rhoi(i)*(LOG(rhoi(i))-1.0)
              END DO
              alphn(k1,k2) = alphn(k1,k2)-fid
              ya(k1+1,k2+1) = alphn(k1,k2)
           END IF
           x1a(k1+1) = rhovec1(k1)
        END DO
        x2a(k2+1) = rhovec2(k2)
     END DO
     ! GOTO 34

     ! DO k2 = 1,stepno+1
     ! DO k1 = 1,stepno+1
     !   ya(k1,k2)
     !   y1a(k1,k2)
     !   y2a(k1,k2)
     !   y12a(k1,k2)
     !   DO m = 1,4
     !     DO n = 1,4
     !       c_bicub(k1,k2,m,n)
     !     ENDDO
     !   ENDDO
     !   x1a(k1)
     ! ENDDO
     ! x2a(k2)
     ! ENDDO


  END DO   ! loop of n cycles
  ! write (71,*) ' '


  ! k2 = 0
  ! DO k1 = 0,stepno/2,1
  !   write (*,*) k1,k2,alph(0,k1,k2),alph(niter,k1,k2) !,del_f(1,k1,k2)
  ! ENDDO
  ! c ENDDO
  ! read (*,*)

  ! 34   CONTINUE
  CALL bicub_derivative2( stepno+1, stepno+1 )
  CALL bicub_c2 ( stepno+1, stepno+1 )


5 CONTINUE


  pos1 = INT(rho1/rhomax1*stepno) + 1  ! plus 1 because the (0:stepno)array is mapped onto a (1:stepno+1)array
  pos2 = INT(rho2/rhomax2*stepno) + 1  ! plus 1 because the (0:stepno)array is mapped onto a (1:stepno+1)array
  CALL bi_cub_spline_crt (rho1,rho2,fres_v,fdr1,fdr2,fdr11,fdr12,fdr22,  &
       fdr111,fdr112,fdr122,fdr222,stepno+1,stepno+1,pos1,pos2)
  fres = fres_v/rho0
  f_drho = x(1)*fdr1 + x(2)*fdr2  ! derivative of Helholtz energy density to rho, d(f'*rho)/d_rho

  f_drho2 = x(1)*x(1)*fdr11 + 2.0*x(1)*x(2)*fdr12 +x(2)*x(2)*fdr22
  f_drho3 = x(1)**3 *fdr111 + 3.0*x(1)*x(1)*x(2)*fdr112  &
       + 3.0*x(1)*x(2)*x(2)*fdr122 + x(2)**3 *fdr222
  f_drho4 = 0.0 ! not calculated !!!

  pges = ( f_drho*rho0 - fres_v ) *(kbol*t)/1.E-30  &
       +   rho0 * (kbol*t) / 1.E-30   ! ideal gas contribution
  pgesdz = ( f_drho2*rho0 ) *(kbol*t)/1.E-30 /z3t  &
       + 1.0/z3t*(kbol*t)/1.E-30   ! ideal gas contribution
  pgesd2 = ( f_drho3*rho0 + f_drho2 ) *(kbol*t)/1.E-30 /z3t /z3t
  pgesd3 = ( f_drho4*rho0 + 2.0*f_drho3 ) *(kbol*t)/1.E-30 /z3t /z3t /z3t

  dfdx(1) = fdr1
  dfdx(2) = fdr2
  ! write (*,*) 'chem.p.',fdr1,fdr2,rho1,rho2

  eta = densav1
  ! write (*,*) 'fp',fres,pges,pgesdz,pgesd2,pgesd3,eta,rho0*z3t
  ! read (*,*)

END SUBROUTINE F_EOS_RN_PHASE_SPACE_CELL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE integr(n,k1,k2,kk1,kk2,kn,combi,combi2, alph,int_l,int_s)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n
  INTEGER, INTENT(IN)                    :: k1
  INTEGER, INTENT(IN)                    :: k2
  INTEGER, INTENT(IN)                    :: kk1
  INTEGER, INTENT(IN)                    :: kk2
  REAL, INTENT(IN)                       :: kn(10,0:300,0:300)
  REAL, INTENT(IN)                       :: combi(0:300,0:300)
  REAL, INTENT(IN)                       :: combi2(10,0:300,0:300)
  REAL, INTENT(IN)                       :: alph(0:8,0:300,0:300)
  REAL, INTENT(OUT)                      :: int_l
  REAL, INTENT(OUT)                      :: int_s


  REAL :: alp_c,alph_lc,alph_sc,alp_l,alph_ll,alph_sl,  &
       alp_r,alph_lr,alph_sr,gn_l,gn_s
  !-----------------------------------------------------------------------------


  alp_c = alph(n-1,k1,k2)
  alph_lc =alp_c + combi(k1,k2)
  alph_sc =alp_c + combi2(n,k1,k2)
  ! write (*,*) 'c',alph(n-1,k1,k2),alph(n-1,k1,k2-kk2) ! ,combi2(n,k1,k2)

  alp_l = alph(n-1,k1-kk1,k2-kk2)
  alph_ll  = alp_l + combi(k1-kk1,k2-kk2)
  alph_sl  = alp_l + combi2(n,k1-kk1,k2-kk2)

  alp_r = alph(n-1,k1+kk1,k2+kk2)
  alph_lr  = alp_r + combi(k1+kk1,k2+kk2)
  alph_sr  = alp_r + combi2(n,k1+kk1,k2+kk2)

  gn_l =  (alph_lr+alph_ll)/2.0 - alph_lc
  gn_s =  (alph_sr+alph_sl)/2.0 - alph_sc
  IF (gn_l < 0.0) gn_l = 0.0
  IF (gn_s < 0.0) gn_s = 0.0
  int_l = 0.0
  int_s = 0.0
  IF( -gn_l/kn(n,k1,k2) > -300.0.AND.-gn_l/kn(n,k1,k2) < 300.0)  int_l = EXP( -gn_l/kn(n,k1,k2) )
  IF( -gn_s/kn(n,k1,k2) > -300.0.AND.-gn_s/kn(n,k1,k2) < 300.0)  int_s = EXP( -gn_s/kn(n,k1,k2) )

END SUBROUTINE integr


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bicub_derivative2 ( i_max, k_max )

  USE RGT_PHASE_SPACE_CELL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                      :: i_max
  INTEGER, INTENT(IN)                      :: k_max
  !-----------------------------------------------------------------------------
  INTEGER :: i,k, steps,lm
  REAL :: d1,d2,gam1,gam2,   &
       d2f_dx1_2(301,301),d2f_dx2_2(301,301)
  REAL :: yspl(301),xspl(301),y1(301),y2(301),yp1,ypn
  !-----------------------------------------------------------------------------

  d1 = x1a(2)-x1a(1)
  d2 = x2a(2)-x2a(1)

  gam1 = SQRT(1.0+d2*d2/d1/d1)/(1.0+d2*d2/d1/d1)
  gam2 = SQRT(1.0+d1*d1/d2/d2)/(1.0+d1*d1/d2/d2)

  DO i = 2,i_max-1
     DO k = 2,k_max-1
        y1a(i,k) = (ya(i+1,k)-ya(i-1,k))/(x1a(i+1)-x1a(i-1))
        y2a(i,k) = (ya(i,k+1)-ya(i,k-1))/(x2a(k+1)-x2a(k-1))
        y12a(i,k) = (ya(i+1,k+1)-ya(i+1,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
             /((x1a(i+1)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
     END DO
  END DO

  i = 1
  DO k = 2,k_max-1
     y1a(i,k) = (-ya(i+2,k)+4.0*ya(i+1,k)-3.0*ya(i,k)) /(x1a(i+2)-x1a(i))
     y2a(i,k) = (ya(i,k+1)-ya(i,k-1))/(x2a(k+1)-x2a(k-1))
     y12a(i,k) = (ya(i+1,k+1)-ya(i+1,k-1)-ya(i,k+1)+ya(i,k-1))  &
          /((x1a(i+1)-x1a(i))*(x2a(k+1)-x2a(k-1)))
  END DO

  k = 1
  DO i = 2,i_max-1
     y1a(i,k) = (ya(i+1,k)-ya(i-1,k))/(x1a(i+1)-x1a(i-1))
     ! write (*,*) y1a(i,k),ya(i,k),(x1a(i+2)-x1a(i))
     ! read (*,*)
     y2a(i,k) = (-ya(i,k+2)+4.0*ya(i,k+1)-3.0*ya(i,k)) /(x2a(k+2)-x2a(k))
     y12a(i,k) = (ya(i+1,k+1)-ya(i+1,k)-ya(i,k+1)+ya(i,k))  &
          /((x1a(i+1)-x1a(i))*(x2a(k+1)-x2a(k)))
  END DO


  i = i_max
  DO k = 2,k_max-1
     y1a(i,k) = (ya(i,k)-ya(i-1,k))/(x1a(i)-x1a(i-1))
     y2a(i,k) = (ya(i,k+1)-ya(i,k-1))/(x2a(k+1)-x2a(k-1))
     y12a(i,k) = (ya(i,k+1)-ya(i,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
          /((x1a(i)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
  END DO


  k = k_max
  DO i = 2,i_max-1
     y1a(i,k) = (ya(i+1,k)-ya(i-1,k))/(x1a(i+1)-x1a(i-1))
     y2a(i,k) = (ya(i,k)-ya(i,k-1))/(x2a(k)-x2a(k-1))
     y12a(i,k) = (ya(i+1,k)-ya(i+1,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
          /((x1a(i+1)-x1a(i-1))*(x2a(k)-x2a(k-1)))
  END DO

  k = k_max
  i = i_max
  y1a(i,k) = (ya(i,k)-ya(i-1,k))/(x1a(i)-x1a(i-1))
  y2a(i,k) = (ya(i,k)-ya(i,k-1))/(x2a(k)-x2a(k-1))
  y12a(i,k) = (ya(i,k)-ya(i,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
       /((x1a(i)-x1a(i-1))*(x2a(k)-x2a(k-1)))

  k = 1
  i = 1
  y1a(i,k) = (-ya(i+2,k)+4.0*ya(i+1,k)-3.0*ya(i,k)) /(x1a(i+2)-x1a(i))
  y2a(i,k) = (-ya(i,k+2)+4.0*ya(i,k+1)-3.0*ya(i,k)) /(x2a(k+2)-x2a(k))
  y12a(i,k) = (ya(i+1,k+1)-ya(i+1,k)-ya(i,k+1)+ya(i,k))  &
       /((x1a(i+1)-x1a(i))*(x2a(k+1)-x2a(k)))


  DO i = 1,i_max
     DO k = 1,k_max
        xspl(k) = x2a(k)
        yspl(k) = ya(i,k)
     END DO
     yp1 = (-ya(i,1+2)+4.0*ya(i,1+1)-3.0*ya(i,1))/(x2a(1+2)-x2a(1))
     ypn = (-ya(i,k_max-2)+4.0*ya(i,k_max-1)-3.0*ya(i,k_max))  &
          /(x2a(k_max-2)-x2a(k_max))
     CALL spline2(xspl,yspl,k_max,yp1,ypn,y1,y2)
     DO k = 1,k_max
        y2a(i,k) =  y1(k)
        d2f_dx2_2(i,k)  =  y2(k)
     END DO
  END DO

  DO k = 1,k_max
     DO i = 1,i_max
        xspl(i) = x1a(i)
        yspl(i) = ya(i,k)
     END DO
     yp1 = (-ya(1+2,k)+4.0*ya(1+1,k)-3.0*ya(1,k))/(x1a(1+2)-x1a(1))
     ypn = (-ya(i_max-2,k)+4.0*ya(i_max-1,k)-3.0*ya(i_max,k))  &
          /(x1a(i_max-2)-x1a(i_max))
     CALL spline2(xspl,yspl,i_max,yp1,ypn,y1,y2)
     DO i = 1,i_max
        y1a(i,k) =  y1(i)
        d2f_dx1_2(i,k) = y2(i)
     END DO
  END DO

  DO lm = 1,k_max-2
     DO i = lm,i_max
        k = i - (lm-1)
        ! steps = (i_max-1) - (lm-1)
        ! dz = (x1a(2)-x1a(1))**2 + (x2a(2)-x2a(1))**2
        ! dz = dz**0.5
        ! fct(k-1) = ya(i,k)
        ! write (*,*) i-1,k-1,steps
        xspl(k) = ( x1a(k)**2 + x2a(k)**2 )**0.5
        yspl(k) = ya(i,k)
     END DO
     !      read (*,*)
     !      CALL SPLINE_PARA (dz,fct,utri,steps)
     !      CALL SPLINE_COEFF(beta,gamma,delta,dz,fct,utri,steps)
     steps = i_max - (lm-1)
     yp1 = y1a(lm,1)*gam1 + y2a(lm,1)*gam2
     ypn = y1a(i_max,steps)*gam1+y2a(i_max,steps)*gam2
     CALL spline2(xspl,yspl,steps,yp1,ypn,y1,y2)
     DO i = lm,i_max
        k = i - (lm-1)
        ! secderiv = 2.0*gamma(k-1)
        ! c  write (*,*) secderiv,y2(k)
        ! secderiv = y2(k)
        ! write (*,*) y12a(i,k),i,k
        y12a(i,k) =  0.5*y2(k)/gam1/gam2 -0.5*(d2f_dx1_2(i,k)*gam1/gam2  &
             +d2f_dx2_2(i,k)*gam2/gam1)
        ! write (*,*) d2f_dx1_2(i,k),d2f_dx2_2(i,k)
        ! write (*,*) y12a(i,k),0.5*secderiv, -0.5*(d2f_dx1_2(i,k)+d2f_dx2_2(i,k))
        ! read (*,*)
     END DO
  END DO

  DO lm = 1,k_max-2
     DO i = 1,i_max-(lm-1)
        k = i + (lm-1)
        !  steps = (i_max-1) - (lm-1)
        ! dz = (x1a(2)-x1a(1))**2 + (x2a(2)-x2a(1))**2
        ! dz = dz**0.5
        ! fct(i-1) = ya(i,k)
        ! c  write (*,*) i,k,steps
        xspl(i) = ( x1a(i)**2 + x2a(i)**2 )**0.5
        yspl(i) = ya(i,k)
     END DO
     ! read (*,*)
     ! CALL SPLINE_PARA (dz,fct,utri,steps)
     ! CALL SPLINE_COEFF(beta,gamma,delta,dz,fct,utri,steps)
     steps = i_max - (lm-1)
     yp1 =  y1a(1,lm)*gam1 +y2a(1,lm)*gam2
     ypn =  y1a(steps,k_max)*gam1+y2a(steps,k_max)*gam2
     CALL spline2(xspl,yspl,steps,yp1,ypn,y1,y2)
     DO i = 1,i_max-(lm-1)
        k = i + (lm-1)
        ! secderiv = 2.0*gamma(i-1)
        ! write (*,*) secderiv,y2(i)
        ! write (*,*) y12a(i,k),i,k
        ! y12a(i,k) =  -0.5*y2(i) +1.0*(d2f_dx1_2(i,k)+d2f_dx2_2(i,k))
        y12a(i,k) =  0.5*y2(i)/gam1/gam2 - 0.5*(d2f_dx1_2(i,k)*gam1/gam2 + d2f_dx2_2(i,k)*gam2/gam1)
        ! write (*,*) y12a(i,k),0.5*y2(i), -0.5*(d2f_dx1_2(i,k)+d2f_dx2_2(i,k))
        ! read (*,*)
     END DO
  END DO

END SUBROUTINE bicub_derivative2


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE spline2
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
! x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
! function at points 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating function at the tabulated
! points xi. If yp1 and/or ypn are equal to 1 � 1030 or larger, the routine is signaled to set
! the corresponding boundary condition for a natural spline, with zero second derivative on
! that boundary.
! Parameter: NMAX is the largest anticipated value of n.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE spline2(x,y,n,yp1,ypn,y1,y2)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)             :: x(n)
  REAL, INTENT(IN)             :: y(n)
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(IN OUT)         :: yp1
  REAL, INTENT(IN OUT)         :: ypn
  REAL, INTENT(OUT)            :: y1(n)
  REAL, INTENT(OUT)            :: y2(n)
  !-----------------------------------------------------------------------------
  REAL :: dx
  INTEGER, PARAMETER :: nmax = 500
  INTEGER :: i,k
  REAL :: p,qn,sig,un,u(nmax)
  !-----------------------------------------------------------------------------


  IF (yp1 > .99E30) THEN
     y2(1) = 0.0
     u(1) = 0.0
  ELSE
     y2(1) = -0.5
     u(1) = (3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  END IF
  DO  i = 2,n-1
     sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
     p = sig*y2(i-1)+2.0
     y2(i) = (sig-1.)/p
     u(i) = (6.0*((y(i+1)-y(i))/(x(i+1)-x(i))  &
          -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  END DO
  IF (ypn > .99D30) THEN
     qn = 0.0
     un = 0.0
  ELSE
     qn = 0.5
     un = (3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  END IF
  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
  DO k = n-1,1,-1
     y2(k) = y2(k)*y2(k+1)+u(k)
  END DO
  DO i = 1,n-1
     dx = x(i+1)-x(i)
     y1(i) = (y(i+1)-y(i))/dx - 1.0/3.0*dx*y2(i)-y2(i+1)/6.0*dx
  END DO
  i = n-1
  dx = x(i+1)-x(i)
  y1(n) = (y(i+1)-y(i))/dx + 1.0/6.0*dx*y2(i)+1.0/3.0*y2(i+1)*dx
  ! y1(1) = yp1
  ! y1(n) = ypn

END SUBROUTINE spline2



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bicub_c2 ( i_max, k_max )

  USE RGT_PHASE_SPACE_CELL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                      :: i_max
  INTEGER, INTENT(IN)                      :: k_max
  !-----------------------------------------------------------------------------
  INTEGER :: i,k,m,n
  REAL :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
  REAL :: c(4,4)
  !-----------------------------------------------------------------------------

  DO i = 1,i_max-1
     DO k = 1,k_max-1
        y(1) = ya(i,k)
        y(2) = ya(i+1,k)
        y(3) = ya(i+1,k+1)
        y(4) = ya(i,k+1)

        y1(1) = y1a(i,k)
        y1(2) = y1a(i+1,k)
        y1(3) = y1a(i+1,k+1)
        y1(4) = y1a(i,k+1)

        y2(1) = y2a(i,k)
        y2(2) = y2a(i+1,k)
        y2(3) = y2a(i+1,k+1)
        y2(4) = y2a(i,k+1)

        y12(1) = y12a(i,k)
        y12(2) = y12a(i+1,k)
        y12(3) = y12a(i+1,k+1)
        y12(4) = y12a(i,k+1)

        x1l = x1a(i)
        x1u = x1a(i+1)
        x2l = x2a(k)
        x2u = x2a(k+1)

        CALL bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
        DO m = 1,4
           DO n = 1,4
              c_bicub(i,k,m,n) = c(m,n)
           END DO
        END DO

     END DO
  END DO


END SUBROUTINE bicub_c2

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE bcuint2
!
! this is a variant of bcuint, where higher derivatives are calculated.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bcuint2 (y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,c,ansy,ansy1,  &
     ansy2,ansy11,ansy12,ansy22,ansy111,ansy112,ansy122,ansy222)

  use utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)         :: y(4)
  REAL, INTENT(IN OUT)         :: y1(4)
  REAL, INTENT(IN OUT)         :: y2(4)
  REAL, INTENT(IN OUT)         :: y12(4)
  REAL, INTENT(IN OUT)         :: x1l
  REAL, INTENT(IN OUT)         :: x1u
  REAL, INTENT(IN OUT)         :: x2l
  REAL, INTENT(IN OUT)         :: x2u
  REAL, INTENT(IN OUT)         :: x1
  REAL, INTENT(IN OUT)         :: x2
  REAL, INTENT(IN)             :: c(4,4)
  REAL, INTENT(OUT)            :: ansy
  REAL, INTENT(OUT)            :: ansy1
  REAL, INTENT(OUT)            :: ansy2
  REAL, INTENT(OUT)            :: ansy11
  REAL, INTENT(OUT)            :: ansy12
  REAL, INTENT(OUT)            :: ansy22
  REAL, INTENT(OUT)            :: ansy111
  REAL, INTENT(OUT)            :: ansy112
  REAL, INTENT(OUT)            :: ansy122
  REAL, INTENT(OUT)            :: ansy222
  !-----------------------------------------------------------------------------
  REAL :: d1,d2
  !U    USES bcucof
  INTEGER :: i,j
  REAL :: t,u
  !-----------------------------------------------------------------------------

  ! call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
  IF( x1u == x1l .OR. x2u == x2l ) call paus ('bad input in bcuint')
  t = (x1-x1l) / (x1u-x1l)
  u = (x2-x2l) / (x2u-x2l)
  IF (u == 0.0) u = 1.E-15
  IF (t == 0.0) t = 1.E-15
  ansy = 0.0
  ansy1 = 0.0
  ansy2 = 0.0
  ansy11 = 0.0
  ansy12 = 0.0
  ansy22 = 0.0
  ansy111 = 0.0
  ansy112 = 0.0
  ansy122 = 0.0
  ansy222 = 0.0
  DO i = 1,4
     DO j = 1,4
        ansy = ansy + c(i,j)*u**REAL(j-1) *t**REAL(i-1)
        ! write (*,*) ansy,c(i,j),u**REAL(j-1),t**REAL(i-1)
        ansy1 =  ansy1+ c(i,j)*REAL(i-1)*u**REAL(j-1) *t**REAL(i-2)
        ansy2 =  ansy2+ c(i,j)*REAL(j-1)*u**REAL(j-2) *t**REAL(i-1)

        ansy11 = ansy11+c(i,j)*REAL(i-1)*REAL(i-2)*u**REAL(j-1) *t**REAL(i-3)
        ansy22 = ansy22+c(i,j)*REAL(j-1)*REAL(j-2)*u**REAL(j-3) *t**REAL(i-1)
        ansy12 = ansy12+c(i,j)*REAL(i-1)*REAL(j-1)*u**REAL(j-2) *t**REAL(i-2)

        ansy111 = ansy111+c(i,j)*REAL(i-1)*REAL(i-2)*REAL(i-3) *u**REAL(j-1) *t**REAL(i-4)
        ansy222 = ansy222+c(i,j)*REAL(j-1)*REAL(j-2)*REAL(j-3) *u**REAL(j-4) *t**REAL(i-1)
        ansy112 = ansy112+c(i,j)*REAL(i-1)*REAL(i-2)*REAL(j-1) *u**REAL(j-2) *t**REAL(i-3)
        ansy122 = ansy122+c(i,j)*REAL(j-1)*REAL(j-2)*REAL(i-1) *u**REAL(j-3) *t**REAL(i-2)
     END DO
  END DO
  ! stop
  d1 = (x1u-x1l)
  d2 = (x2u-x2l)
  ansy1 = ansy1/d1
  ansy2 = ansy2/d2

  ansy11 = ansy11/d1/d1
  ansy22 = ansy22/d2/d2
  ansy12 = ansy12/d1/d2

  ansy111 = ansy111/d1**3
  ansy222 = ansy222/d2**3
  ansy112 = ansy112/d1/d1/d2
  ansy122 = ansy122/d1/d2/d2

END SUBROUTINE bcuint2



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE bi_cub_spline_crt
!
! ............
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE bi_cub_spline_crt (rho1,rho2, fr,fdr1,fdr2,fdr11,fdr12,fdr22,  &
     fdr111,fdr112,fdr122,fdr222,i_max,k_max,ih,k)

  USE RGT_PHASE_SPACE_CELL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: rho1
  REAL, INTENT(IN OUT)                   :: rho2
  REAL, INTENT(OUT)                      :: fr
  REAL, INTENT(IN OUT)                   :: fdr1
  REAL, INTENT(IN OUT)                   :: fdr2
  REAL, INTENT(IN OUT)                   :: fdr11
  REAL, INTENT(IN OUT)                   :: fdr12
  REAL, INTENT(IN OUT)                   :: fdr22
  REAL, INTENT(IN OUT)                   :: fdr111
  REAL, INTENT(IN OUT)                   :: fdr112
  REAL, INTENT(IN OUT)                   :: fdr122
  REAL, INTENT(IN OUT)                   :: fdr222
  INTEGER, INTENT(IN)                    :: i_max
  INTEGER, INTENT(IN)                    :: k_max
  INTEGER, INTENT(OUT)                   :: ih
  INTEGER, INTENT(IN)                    :: k
  !-----------------------------------------------------------------------------
  INTEGER :: m,n
  REAL :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
  REAL :: c(4,4)
  !-----------------------------------------------------------------------------

  IF (rho1 < x1a(1) .AND. rho2 < x2a(1)) THEN
     fr = 1.0
     WRITE (*,*) 'error? BI_CUB_SPLINE_CRT'
     STOP
     RETURN
  END IF
  !      write (*,*) ih,k
  IF (x1a(ih) <= rho1.AND.rho1 < x1a(ih+1).AND.  &
       x2a(k) <= rho2.AND.rho2 < x2a(k+1)) THEN
     GO TO 10
  ELSE
     WRITE (*,*) 'error in BI_CUB_SPLINE_CRT',ih,k
     WRITE (*,*) rho1,x1a(ih),x1a(ih+1)
     WRITE (*,*) rho2,x2a(k),x2a(k+1)
     STOP
  END IF
  IF (ih > 2) THEN
     IF (x1a(ih-1) <= rho1.AND.rho1 < x1a(ih)) THEN
        ih = ih-1
        GO TO 10
     END IF
  END IF
  WRITE (*,*) 'error in BI_CUB_SPLINE_CRT'
  CALL hunt(x1a,i_max,rho1,ih)
10 CONTINUE

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

  DO m = 1,4
     DO n = 1,4
        c(m,n) = c_bicub(ih,k,m,n)
     END DO
  END DO
  CALL bcuint2(y,y1,y2,y12,x1l,x1u,x2l,x2u,rho1,rho2,c,  &
       fr,fdr1,fdr2,fdr11,fdr12,fdr22,fdr111,fdr112,fdr122,fdr222)

END SUBROUTINE bi_cub_spline_crt



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PHI_CRITICAL_RENORM_PHASE_SPACE_CELL
!
! ............
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_CRITICAL_RENORM_PHASE_SPACE_CELL

  USE parameters, ONLY: KBOL
  USE EOS_CONSTANTS
  USE EOS_VARIABLES
  USE RGT_PHASE_SPACE_CELL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: k
  REAL                                   :: zres, zges
  !-----------------------------------------------------------------------------

  CALL DENSITY_ITERATION

  !-----------------------------------------------------------------------------
  zges = (p * 1.E-30) / (KBOL*t*eta/z3t)
  IF ( ensemble_flag == 'tv' ) zges = (pges * 1.E-30) / (kbol*t*eta/z3t)
  zres = zges - 1.0

  DO  k = 1, ncomp
     IF (ensemble_flag == 'tp') lnphi(k) = dfdx(k) - LOG(zges)
     IF (ensemble_flag == 'tv' .AND. eta >= 0.0) lnphi(k) = dfdx(k) ! +LOG(rho)
  END DO

END SUBROUTINE PHI_CRITICAL_RENORM_PHASE_SPACE_CELL
