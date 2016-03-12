MODULE EOS_POLAR

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: F_POLAR, P_POLAR, PHI_POLAR

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_POLAR ( fdd, fqq, fdq )

  USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: fdd, fqq, fdq

  !-----------------------------------------------------------------------------
  INTEGER                                :: dipole
  INTEGER                                :: quadrupole
  INTEGER                                :: dipole_quad
  !-----------------------------------------------------------------------------

  fdd = 0.0
  fqq = 0.0
  fdq = 0.0

  dipole      = 0
  quadrupole  = 0
  dipole_quad = 0
  IF ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
  IF ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
  IF ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

  !-----------------------------------------------------------------------------
  ! dipole-dipole term
  !-----------------------------------------------------------------------------
  IF (dipole == 1) THEN

     IF (dd_term == 'GV') CALL F_DD_GROSS_VRABEC( fdd )
     ! IF (dd_term == 'SF') CALL F_DD_SAAGER_FISCHER( k )
     ! IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  IF (quadrupole == 1) THEN

     IF (qq_term == 'JG') CALL F_QQ_GROSS( fqq )
     ! IF (qq_term == 'SF') CALL F_QQ_SAAGER_FISCHER( k )
     ! IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! dipole-quadrupole cross term
  !-----------------------------------------------------------------------------
  IF (dipole_quad == 1) THEN

     IF (dq_term == 'VG') CALL F_DQ_VRABEC_GROSS( fdq )
     ! IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

  ENDIF

END SUBROUTINE F_POLAR


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_DD_GROSS_VRABEC( fdd )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdd
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  INTEGER                                :: ddit, ddmax
  REAL                                   :: factor2, factor3
  REAL                                   :: xijfa, xijkfa, xijf_j, xijkf_j, eij
  REAL                                   :: fdd2, fdd3
  REAL, DIMENSION(nc)                    :: my2dd, my0, alph_tst, z1dd, z2dd, dderror
  REAL, DIMENSION(nc)                    :: fdd2m, fdd3m, fdd2m2, fdd3m2, fddm, fddm2
  REAL, DIMENSION(nc,nc)                 :: Idd2, Idd4
  REAL, DIMENSION(nc,nc,nc)              :: Idd3
  !-----------------------------------------------------------------------------

  fdd    = 0.0
  ddit   = 0
  ddmax  = 0   ! value assigned, if polarizable compound is present
  fddm(:) = 0.0
  DO i = 1, ncomp
     IF ( uij(i,i) == 0.0 ) write (*,*) 'F_DD_GROSS_VRABEC: do not use dimensionless units'
     IF ( uij(i,i) == 0.0 ) stop
     my2dd(i) = (parame(i,6))**2 *1.E-49 /(uij(i,i)*kbol* mseg(i)*sig_ij(i,i)**3 *1.E-30)
     alph_tst(i) = parame(i,11) / (mseg(i)*sig_ij(i,i)**3 ) * t/parame(i,3)
     IF ( alph_Tst(i) /= 0.0 ) ddmax = 25     ! set maximum number of polarizable RGT-iterations
     z1dd(i) = my2dd(i) + 3.0*alph_tst(i)
     z2dd(i) = 3.0*alph_tst(i)
     my0(i)  = my2dd(i)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Idd2(i,j) = 0.0
        Idd4(i,j) = 0.0
        ! IF (parame(i,6).NE.0.0 .AND. parame(j,6).NE.0.0) THEN
        DO m = 0, 4
           Idd2(i,j) = Idd2(i,j) + ddp2(i,j,m)*eta**m
           Idd4(i,j) = Idd4(i,j) + ddp4(i,j,m)*eta**m
        END DO
        DO k = 1, ncomp
           Idd3(i,j,k) = 0.0
           ! IF (parame(k,6).NE.0.0) THEN
           DO m = 0, 4
              Idd3(i,j,k) = Idd3(i,j,k) + ddp3(i,j,k,m)*eta**m
           END DO
           ! ENDIF
        END DO
        ! ENDIF
     END DO
  END DO

  factor2 = -PI *rho
  factor3 = -4.0/3.0*PI**2 * rho**2

9 CONTINUE

  fdd2m(:)  = 0.0
  fdd2m2(:) = 0.0
  fdd3m(:)  = 0.0
  fdd3m2(:) = 0.0
  fdd2 = 0.0
  fdd3 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        ! IF (parame(i,6).NE.0.0 .AND. parame(j,6).NE.0.0) THEN
        xijfa =x(i)*parame(i,3)/t*parame(i,2)**3  * x(j)*parame(j,3)/t*parame(j,2)**3   &
             /((parame(i,2)+parame(j,2))/2.0)**3  * (z1dd(i)*z1dd(j)-z2dd(i)*z2dd(j))    ! * (1.0-lij(i,j))
        eij = (parame(i,3)*parame(j,3))**0.5
        fdd2= fdd2 + factor2 * xijfa * ( Idd2(i,j) + eij/t*Idd4(i,j) )
        xijf_j = parame(i,3)/t*parame(i,2)**3  *x(j)*parame(j,3)/t*parame(j,2)**3   &
             /((parame(i,2)+parame(j,2))/2.0)**3       ! * (1.0-lij(i,j))
        fdd2m(i)=fdd2m(i)+4.0*SQRT(my2dd(i))*z1dd(j)*factor2* xijf_j *(Idd2(i,j)+eij/t*Idd4(i,j))
        fdd2m2(i)=fdd2m2(i) + 4.0*z1dd(j)*factor2* xijf_j *(Idd2(i,j)+eij/t*Idd4(i,j))
        IF (j == i) fdd2m2(i) =fdd2m2(i) +8.0*factor2* xijf_j*my2dd(i) *(Idd2(i,j)+eij/t*Idd4(i,j))
        DO k = 1, ncomp
           ! IF (parame(k,6).NE.0.0) THEN
           xijkfa = x(i)*parame(i,3)/t*parame(i,2)**3  *x(j)*parame(j,3)/t*parame(j,2)**3   &
                *x(k)*parame(k,3)/t*parame(k,2)**3  / ((parame(i,2)+parame(j,2))/2.0)  &
                /((parame(i,2)+parame(k,2))/2.0) / ((parame(j,2)+parame(k,2))/2.0)  &
                *(z1dd(i)*z1dd(j)*z1dd(k)-z2dd(i)*z2dd(j)*z2dd(k))
           ! *(1.0-lij(i,j))*(1.0-lij(i,k))*(1.0-lij(j,k))
           fdd3 = fdd3 + factor3 * xijkfa * Idd3(i,j,k)
           xijkf_j = parame(i,3)/t*parame(i,2)**3  *x(j)*parame(j,3)/t*parame(j,2)**3   &
                *x(k)*parame(k,3)/t*parame(k,2)**3  /((parame(i,2)+parame(j,2))/2.0)  &
                /((parame(i,2)+parame(k,2))/2.0) /((parame(j,2)+parame(k,2))/2.0)
           ! *(1.0-lij(i,j))*(1.0-lij(i,k))*(1.0-lij(j,k))
           fdd3m(i)=fdd3m(i)+6.0*factor3*SQRT(my2dd(i))*z1dd(j)*z1dd(k) *xijkf_j*Idd3(i,j,k)
           fdd3m2(i)=fdd3m2(i)+6.0*factor3*z1dd(j)*z1dd(k) *xijkf_j*Idd3(i,j,k)
           IF(j == i) fdd3m2(i) =fdd3m2(i)+24.0*factor3*my2dd(i)*z1dd(k) *xijkf_j*Idd3(i,j,k)
           ! ENDIF
        END DO
        ! ENDIF
     END DO
  END DO

  IF (fdd2 < -1.E-50 .AND. fdd3 /= 0.0) THEN
     fdd = fdd2 / ( 1.0 - fdd3/fdd2 )
     IF ( ddmax /= 0 ) THEN
        DO i = 1, ncomp
           ddit = ddit + 1
           fddm(i) =fdd2*(fdd2*fdd2m(i) -2.0*fdd3*fdd2m(i)+fdd2*fdd3m(i)) /(fdd2-fdd3)**2
           fddm2(i) = fdd2m(i) * (fdd2*fdd2m(i)-2.0*fdd3*fdd2m(i) +fdd2*fdd3m(i)) / (fdd2-fdd3)**2  &
                + fdd2*(fdd2*fdd2m2(i) -2.0*fdd3*fdd2m2(i)+fdd2m(i)**2  &
                -fdd2m(i)*fdd3m(i) +fdd2*fdd3m2(i)) / (fdd2-fdd3)**2  &
                - 2.0*fdd2*(fdd2*fdd2m(i) -2.0*fdd3*fdd2m(i) +fdd2*fdd3m(i)) /(fdd2-fdd3)**3   &
                *(fdd2m(i)-fdd3m(i))
           dderror(i)= SQRT( my2dd(i) ) - SQRT( my0(i) ) + alph_Tst(i)*fddm(i)
           my2dd(i) = ( SQRT( my2dd(i) ) - dderror(i) / (1.0+alph_Tst(i)*fddm2(i)) )**2
           z1dd(i) = my2dd(i) + 3.0 * alph_Tst(i)
        ENDDO
        DO i = 1, ncomp
           IF (ABS(dderror(i)) > 1.E-11 .AND. ddit < ddmax) GOTO 9
        ENDDO
        fdd = fdd + SUM( 0.5*x(1:ncomp)*alph_Tst(1:ncomp)*fddm(1:ncomp)**2 )
     ENDIF
  END IF


END SUBROUTINE F_DD_GROSS_VRABEC



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_QQ_GROSS( fqq )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fqq
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  REAL                                   :: factor2, factor3
  REAL                                   :: xijfa, xijkfa, eij
  REAL                                   :: fqq2, fqq3
  REAL, DIMENSION(nc)                    :: qq2
  REAL, DIMENSION(nc,nc)                 :: Iqq2, Iqq4
  REAL, DIMENSION(nc,nc,nc)              :: Iqq3
  !-----------------------------------------------------------------------------


  fqq = 0.0
  DO i = 1, ncomp
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Iqq2(i,j) = 0.0
        Iqq4(i,j) = 0.0
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           DO m = 0, 4
              Iqq2(i,j) = Iqq2(i,j) + qqp2(i,j,m)*eta**m
              Iqq4(i,j) = Iqq4(i,j) + qqp4(i,j,m)*eta**m
           END DO
           DO k = 1, ncomp
              Iqq3(i,j,k) = 0.0
              IF (parame(k,7) /= 0.0) THEN
                 DO m = 0, 4
                    Iqq3(i,j,k) = Iqq3(i,j,k) + qqp3(i,j,k,m)*eta**m
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO

  factor2 = -9.0/16.0*PI *rho
  factor3 =  9.0/16.0*PI**2 * rho**2

  fqq2 = 0.0
  fqq3 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           xijfa=x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
                *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7
           eij = (parame(i,3)*parame(j,3))**0.5
           fqq2= fqq2 +factor2* xijfa * (Iqq2(i,j)+eij/t*Iqq4(i,j))
           DO k = 1, ncomp
              IF (parame(k,7) /= 0.0) THEN
                 xijkfa=x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t/sig_ij(i,j)**3   &
                      *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,k)**3   &
                      *x(k)*uij(k,k)*qq2(k)*sig_ij(k,k)**5 /t/sig_ij(j,k)**3
                 fqq3 = fqq3 + factor3 * xijkfa * Iqq3(i,j,k)
              END IF
           END DO
        END IF
     END DO
  END DO

  IF ( fqq2 < -1.E-50 .AND. fqq3 /= 0.0 ) THEN
     fqq = fqq2 / ( 1.0 - fqq3/fqq2 )
  END IF



END SUBROUTINE F_QQ_GROSS

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_DQ_VRABEC_GROSS( fdq )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdq
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  REAL                                   :: factor2, factor3
  REAL                                   :: xijfa, xijkfa, eij
  REAL                                   :: fdq2, fdq3
  REAL, DIMENSION(nc)                    :: my2dd, myfac, qq2, q_fac
  REAL, DIMENSION(nc,nc)                 :: Idq2, Idq4
  REAL, DIMENSION(nc,nc,nc)              :: Idq3
  !-----------------------------------------------------------------------------


  fdq = 0.0
  DO i = 1, ncomp
     my2dd(i) = (parame(i,6))**2 *1.E-49 /(uij(i,i)*kbol* mseg(i)*sig_ij(i,i)**3  *1.E-30)
     myfac(i) = parame(i,3)/t*parame(i,2)**4  *my2dd(i)
     ! myfac(i)=parame(i,3)/T*parame(i,2)**4  *my2dd_renormalized(i)
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     q_fac(i) = parame(i,3)/t*parame(i,2)**4  *qq2(i)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Idq2(i,j) = 0.0
        Idq4(i,j) = 0.0
        IF (myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0) THEN
           DO m = 0, 4
              Idq2(i,j) = Idq2(i,j) + dqp2(i,j,m)*eta**m
              Idq4(i,j) = Idq4(i,j) + dqp4(i,j,m)*eta**m
           END DO
           DO k = 1, ncomp
              Idq3(i,j,k) = 0.0
              IF (myfac(k) /= 0.0 .OR. q_fac(k) /= 0.0) THEN
                 DO m = 0, 4
                    Idq3(i,j,k) = Idq3(i,j,k) + dqp3(i,j,k,m)*eta**m
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO

  factor2 = -9.0/4.0 * PI *rho
  factor3 =  PI**2 * rho**2

  fdq2 =  0.0
  fdq3 =  0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        IF (myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0) THEN
           xijfa = x(i)*myfac(i) * x(j)*q_fac(j) /sig_ij(i,j)**5
           eij  = (parame(i,3)*parame(j,3))**0.5
           fdq2 = fdq2 +factor2* xijfa*(Idq2(i,j)+eij/t*Idq4(i,j))
           DO k = 1, ncomp
              IF (myfac(k) /= 0.0 .OR. q_fac(k) /= 0.0) THEN
                 xijkfa=x(i)*x(j)*x(k)/(sig_ij(i,j)*sig_ij(i,k)*sig_ij(j,k))**2  &
                      *( myfac(i)*q_fac(j)*myfac(k) + myfac(i)*q_fac(j)*q_fac(k)*1.1937350 )
                 fdq3 = fdq3 + factor3*xijkfa*Idq3(i,j,k)
              END IF
           END DO
        END IF
     END DO
  END DO

  IF (fdq2 < -1.E-50 .AND. fdq3 /= 0.0) THEN
     fdq = fdq2 / ( 1.0 - fdq3/fdq2 )
  END IF

END SUBROUTINE F_DQ_VRABEC_GROSS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_POLAR ( zdd, zddz, zddz2, zddz3, zqq, zqqz, zqqz2, zqqz3, zdq, zdqz, zdqz2, zdqz3 )

  USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term

  !-----------------------------------------------------------------------------
  REAL, INTENT(OUT)                      :: zdd, zddz, zddz2, zddz3
  REAL, INTENT(OUT)                      :: zqq, zqqz, zqqz2, zqqz3
  REAL, INTENT(OUT)                      :: zdq, zdqz, zdqz2, zdqz3

  !-----------------------------------------------------------------------------
  INTEGER                                :: dipole
  INTEGER                                :: quadrupole
  INTEGER                                :: dipole_quad
  !-----------------------------------------------------------------------------

  zdd   = 0.0
  zddz  = 0.0
  zddz2 = 0.0
  zddz3 = 0.0
  zqq   = 0.0
  zqqz  = 0.0
  zqqz2 = 0.0
  zqqz3 = 0.0
  zdq   = 0.0
  zdqz  = 0.0
  zdqz2 = 0.0
  zdqz3 = 0.0

  dipole      = 0
  quadrupole  = 0
  dipole_quad = 0
  IF ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
  IF ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
  IF ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

  !-----------------------------------------------------------------------------
  ! dipole-dipole term
  !-----------------------------------------------------------------------------
  IF (dipole == 1) THEN

     IF (dd_term == 'GV') CALL P_DD_GROSS_VRABEC( zdd, zddz, zddz2, zddz3 )
     ! IF (dd_term == 'SF') CALL F_DD_SAAGER_FISCHER( k )
     IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  IF (quadrupole == 1) THEN

     !IF (qq_term == 'SF') CALL F_QQ_SAAGER_FISCHER( k )
     IF (qq_term == 'JG') CALL P_QQ_GROSS( zqq, zqqz, zqqz2, zqqz3 )
     IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! dipole-quadrupole cross term
  !-----------------------------------------------------------------------------
  IF (dipole_quad == 1) THEN

     IF (dq_term == 'VG') CALL P_DQ_VRABEC_GROSS( zdq, zdqz, zdqz2, zdqz3 )
     IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

  ENDIF

END SUBROUTINE P_POLAR


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_DD_GROSS_VRABEC( zdd, zddz, zddz2, zddz3 )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: zdd, zddz, zddz2, zddz3
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, eij
  REAL                                   :: fdddr, fddd2, fddd3, fddd4
  REAL                                   :: fdd2, fdd2z, fdd2z2, fdd2z3, fdd2z4
  REAL                                   :: fdd3, fdd3z, fdd3z2, fdd3z3, fdd3z4
  REAL, DIMENSION(nc)                    :: my2dd
  REAL, DIMENSION(nc,nc)                 :: Idd2, Idd2z, Idd2z2, Idd2z3, Idd2z4
  REAL, DIMENSION(nc,nc)                 :: Idd4, Idd4z, Idd4z2, Idd4z3, Idd4z4
  REAL, DIMENSION(nc,nc,nc)              :: Idd3, Idd3z, Idd3z2, Idd3z3, Idd3z4
  !-----------------------------------------------------------------------------


  zdd   = 0.0
  zddz  = 0.0
  zddz2 = 0.0
  zddz3 = 0.0
  z3 = eta
  DO i = 1, ncomp
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Idd2(i,j)   = 0.0
        Idd4(i,j)   = 0.0
        Idd2z(i,j)  = 0.0
        Idd4z(i,j)  = 0.0
        Idd2z2(i,j) = 0.0
        Idd4z2(i,j) = 0.0
        Idd2z3(i,j) = 0.0
        Idd4z3(i,j) = 0.0
        Idd2z4(i,j) = 0.0
        Idd4z4(i,j) = 0.0
        ! IF (paramei,6).NE.0.0 .AND. parame(j,6).NE.0.0) THEN
        DO m = 0, 4
           Idd2(i,j)  =Idd2(i,j) + ddp2(i,j,m) *z3**(m+1)
           Idd4(i,j)  =Idd4(i,j) + ddp4(i,j,m) *z3**(m+1)
           Idd2z(i,j) =Idd2z(i,j) +ddp2(i,j,m)*REAL(m+1) *z3**m
           Idd4z(i,j) =Idd4z(i,j) +ddp4(i,j,m)*REAL(m+1) *z3**m
           Idd2z2(i,j)=Idd2z2(i,j)+ddp2(i,j,m)*REAL((m+1)*m) *z3**(m-1)
           Idd4z2(i,j)=Idd4z2(i,j)+ddp4(i,j,m)*REAL((m+1)*m) *z3**(m-1)
           Idd2z3(i,j)=Idd2z3(i,j)+ddp2(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
           Idd4z3(i,j)=Idd4z3(i,j)+ddp4(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
           Idd2z4(i,j)=Idd2z4(i,j)+ddp2(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
           Idd4z4(i,j)=Idd4z4(i,j)+ddp4(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
        END DO
        DO k = 1, ncomp
           Idd3(i,j,k)   = 0.0
           Idd3z(i,j,k)  = 0.0
           Idd3z2(i,j,k) = 0.0
           Idd3z3(i,j,k) = 0.0
           Idd3z4(i,j,k) = 0.0
           ! IF (parame(k,6).NE.0.0) THEN
           DO m = 0, 4
              Idd3(i,j,k)  =Idd3(i,j,k)  +ddp3(i,j,k,m)*z3**(m+2)
              Idd3z(i,j,k) =Idd3z(i,j,k) +ddp3(i,j,k,m)*REAL(m+2)*z3**(m+1)
              Idd3z2(i,j,k)=Idd3z2(i,j,k)+ddp3(i,j,k,m)*REAL((m+2)*(m+1))*z3**m
              Idd3z3(i,j,k)=Idd3z3(i,j,k)+ddp3(i,j,k,m)*REAL((m+2)*(m+1)*m)*z3**(m-1)
              Idd3z4(i,j,k)=Idd3z4(i,j,k)+ddp3(i,j,k,m)*REAL((m+2)*(m+1)*m*(m-1)) *z3**(m-2)
           END DO
           ! ENDIF
        END DO
        ! ENDIF
     END DO
  END DO

  factor2= -PI *rho/z3
  factor3= -4.0/3.0*PI**2 * (rho/z3)**2

  fdd2   = 0.0
  fdd2z  = 0.0
  fdd2z2 = 0.0
  fdd2z3 = 0.0
  fdd2z4 = 0.0
  fdd3   = 0.0
  fdd3z  = 0.0
  fdd3z2 = 0.0
  fdd3z3 = 0.0
  fdd3z4 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        ! IF (parame(i,6).NE.0.0 .AND. parame(j,6).NE.0.0) THEN
        xijfa = x(i)*parame(i,3)/t*parame(i,2)**3  *x(j)*parame(j,3)/t*parame(j,2)**3   &
             / ((parame(i,2)+parame(j,2))/2.0)**3  *my2dd(i)*my2dd(j)
        eij   = (parame(i,3)*parame(j,3))**0.5
        fdd2   = fdd2  +factor2*xijfa*(Idd2(i,j)  +eij/t*Idd4(i,j))
        fdd2z  = fdd2z +factor2*xijfa*(Idd2z(i,j) +eij/t*Idd4z(i,j))
        fdd2z2 = fdd2z2+factor2*xijfa*(Idd2z2(i,j)+eij/t*Idd4z2(i,j))
        fdd2z3 = fdd2z3+factor2*xijfa*(Idd2z3(i,j)+eij/t*Idd4z3(i,j))
        fdd2z4 = fdd2z4+factor2*xijfa*(Idd2z4(i,j)+eij/t*Idd4z4(i,j))
        DO k = 1, ncomp
           ! IF (parame(k,6).NE.0.0) THEN
           xijkfa= x(i)*parame(i,3)/t*parame(i,2)**3  *x(j)*parame(j,3)/t*parame(j,2)**3   &
                *x(k)*parame(k,3)/t*parame(k,2)**3  /((parame(i,2)+parame(j,2))/2.0)  &
                /((parame(i,2)+parame(k,2))/2.0) /((parame(j,2)+parame(k,2))/2.0)  &
                *my2dd(i)*my2dd(j)*my2dd(k)
           fdd3   = fdd3   + factor3 * xijkfa*Idd3(i,j,k)
           fdd3z  = fdd3z  + factor3 * xijkfa*Idd3z(i,j,k)
           fdd3z2 = fdd3z2 + factor3 * xijkfa*Idd3z2(i,j,k)
           fdd3z3 = fdd3z3 + factor3 * xijkfa*Idd3z3(i,j,k)
           fdd3z4 = fdd3z4 + factor3 * xijkfa*Idd3z4(i,j,k)
           ! ENDIF
        END DO
        ! ENDIF
     END DO
  END DO
  IF (fdd2 < -1.E-50 .AND. fdd3 /= 0.0 .AND. fdd2z /= 0.0 .AND. fdd3z /= 0.0) THEN

     fdddr= fdd2* (fdd2*fdd2z - 2.0*fdd3*fdd2z+fdd2*fdd3z) / (fdd2-fdd3)**2
     fddd2=(2.0*fdd2*fdd2z*fdd2z +fdd2*fdd2*fdd2z2  &
          -2.0*fdd2z**2 *fdd3-2.0*fdd2*fdd2z2*fdd3+fdd2*fdd2*fdd3z2)  &
          /(fdd2-fdd3)**2 + fdddr * 2.0*(fdd3z-fdd2z)/(fdd2-fdd3)
     fddd3=(2.0*fdd2z**3  +6.0*fdd2*fdd2z*fdd2z2+fdd2*fdd2*fdd2z3  &
          -6.0*fdd2z*fdd2z2*fdd3-2.0*fdd2z**2 *fdd3z  &
          -2.0*fdd2*fdd2z3*fdd3 -2.0*fdd2*fdd2z2*fdd3z  &
          +2.0*fdd2*fdd2z*fdd3z2+fdd2*fdd2*fdd3z3) /(fdd2-fdd3)**2  &
          + 2.0/(fdd2-fdd3)* ( 2.0*fddd2*(fdd3z-fdd2z)  &
          +     fdddr*(fdd3z2-fdd2z2)  &
          -     fdddr/(fdd2-fdd3)*(fdd3z-fdd2z)**2 )
     fddd4=( 12.0*fdd2z**2 *fdd2z2+6.0*fdd2*fdd2z2**2  &
          +8.0*fdd2*fdd2z*fdd2z3+fdd2*fdd2*fdd2z4-6.0*fdd2z2**2 *fdd3  &
          -12.0*fdd2z*fdd2z2*fdd3z -8.0*fdd2z*fdd2z3*fdd3  &
          -2.0*fdd2*fdd2z4*fdd3-4.0*fdd2*fdd2z3*fdd3z  &
          +4.0*fdd2*fdd2z*fdd3z3+fdd2**2 *fdd3z4 ) /(fdd2-fdd3)**2  &
          + 6.0/(fdd2-fdd3)* ( fddd3*(fdd3z-fdd2z)  &
          -fddd2/(fdd2-fdd3)*(fdd3z-fdd2z)**2  &
          - fdddr/(fdd2-fdd3)*(fdd3z-fdd2z)*(fdd3z2-fdd2z2)  &
          + fddd2*(fdd3z2-fdd2z2) +1.0/3.0*fdddr*(fdd3z3-fdd2z3) )
     zdd   = fdddr*eta
     zddz  = fddd2*eta + fdddr
     zddz2 = fddd3*eta + 2.0* fddd2
     zddz3 = fddd4*eta + 3.0* fddd3

  END IF


END SUBROUTINE P_DD_GROSS_VRABEC



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_QQ_GROSS( zqq, zqqz, zqqz2, zqqz3 )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: zqq, zqqz, zqqz2, zqqz3
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, eij
  REAL                                   :: fqqdr, fqqd2, fqqd3, fqqd4
  REAL                                   :: fqq2, fqq2z, fqq2z2, fqq2z3, fqq2z4
  REAL                                   :: fqq3, fqq3z, fqq3z2, fqq3z3, fqq3z4
  REAL, DIMENSION(nc)                    :: qq2
  REAL, DIMENSION(nc,nc)                 :: Iqq2, Iqq2z, Iqq2z2, Iqq2z3, Iqq2z4
  REAL, DIMENSION(nc,nc)                 :: Iqq4, Iqq4z, Iqq4z2, Iqq4z3, Iqq4z4
  REAL, DIMENSION(nc,nc,nc)              :: Iqq3, Iqq3z, Iqq3z2, Iqq3z3, Iqq3z4
  !-----------------------------------------------------------------------------

  zqq   = 0.0
  zqqz  = 0.0
  zqqz2 = 0.0
  zqqz3 = 0.0
  z3 = eta
  DO i=1,ncomp
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Iqq2(i,j)   = 0.0
        Iqq4(i,j)   = 0.0
        Iqq2z(i,j)  = 0.0
        Iqq4z(i,j)  = 0.0
        Iqq2z2(i,j) = 0.0
        Iqq4z2(i,j) = 0.0
        Iqq2z3(i,j) = 0.0
        Iqq4z3(i,j) = 0.0
        Iqq2z4(i,j) = 0.0
        Iqq4z4(i,j) = 0.0
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           DO m = 0, 4
              Iqq2(i,j)  =Iqq2(i,j) + qqp2(i,j,m)*z3**(m+1)
              Iqq4(i,j)  =Iqq4(i,j) + qqp4(i,j,m)*z3**(m+1)
              Iqq2z(i,j) =Iqq2z(i,j) +qqp2(i,j,m)*REAL(m+1)*z3**m
              Iqq4z(i,j) =Iqq4z(i,j) +qqp4(i,j,m)*REAL(m+1)*z3**m
              Iqq2z2(i,j)=Iqq2z2(i,j)+qqp2(i,j,m)*REAL((m+1)*m)*z3**(m-1)
              Iqq4z2(i,j)=Iqq4z2(i,j)+qqp4(i,j,m)*REAL((m+1)*m)*z3**(m-1)
              Iqq2z3(i,j)=Iqq2z3(i,j)+qqp2(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
              Iqq4z3(i,j)=Iqq4z3(i,j)+qqp4(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
              Iqq2z4(i,j)=Iqq2z4(i,j)+qqp2(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
              Iqq4z4(i,j)=Iqq4z4(i,j)+qqp4(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
           END DO
           DO k=1,ncomp
              Iqq3(i,j,k)   = 0.0
              Iqq3z(i,j,k)  = 0.0
              Iqq3z2(i,j,k) = 0.0
              Iqq3z3(i,j,k) = 0.0
              Iqq3z4(i,j,k) = 0.0
              IF (parame(k,7) /= 0.0) THEN
                 DO m=0,4
                    Iqq3(i,j,k)  =Iqq3(i,j,k) + qqp3(i,j,k,m)*z3**(m+2)
                    Iqq3z(i,j,k)=Iqq3z(i,j,k)+qqp3(i,j,k,m)*REAL(m+2)*z3**(m+1)
                    Iqq3z2(i,j,k)=Iqq3z2(i,j,k)+qqp3(i,j,k,m)*REAL((m+2)*(m+1)) *z3**m
                    Iqq3z3(i,j,k)=Iqq3z3(i,j,k)+qqp3(i,j,k,m)*REAL((m+2)*(m+1)*m) *z3**(m-1)
                    Iqq3z4(i,j,k)=Iqq3z4(i,j,k)+qqp3(i,j,k,m) *REAL((m+2)*(m+1)*m*(m-1)) *z3**(m-2)
                 END DO
              END IF
           END DO

        END IF
     END DO
  END DO

  factor2= -9.0/16.0*PI *rho/z3
  factor3=  9.0/16.0*PI**2 * (rho/z3)**2

  fqq2   = 0.0
  fqq2z  = 0.0
  fqq2z2 = 0.0
  fqq2z3 = 0.0
  fqq2z4 = 0.0
  fqq3   = 0.0
  fqq3z  = 0.0
  fqq3z2 = 0.0
  fqq3z3 = 0.0
  fqq3z4 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           xijfa =x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
                *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7
           eij = (parame(i,3)*parame(j,3))**0.5
           fqq2=  fqq2  +factor2*xijfa*(Iqq2(i,j)  +eij/t*Iqq4(i,j)  )
           fqq2z =fqq2z +factor2*xijfa*(Iqq2z(i,j) +eij/t*Iqq4z(i,j) )
           fqq2z2=fqq2z2+factor2*xijfa*(Iqq2z2(i,j)+eij/t*Iqq4z2(i,j))
           fqq2z3=fqq2z3+factor2*xijfa*(Iqq2z3(i,j)+eij/t*Iqq4z3(i,j))
           fqq2z4=fqq2z4+factor2*xijfa*(Iqq2z4(i,j)+eij/t*Iqq4z4(i,j))
           DO k = 1, ncomp
              IF (parame(k,7) /= 0.0) THEN
                 xijkfa=x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t/sig_ij(i,j)**3   &
                      *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,k)**3   &
                      *x(k)*uij(k,k)*qq2(k)*sig_ij(k,k)**5 /t/sig_ij(j,k)**3
                 fqq3   = fqq3   + factor3 * xijkfa*Iqq3(i,j,k)
                 fqq3z  = fqq3z  + factor3 * xijkfa*Iqq3z(i,j,k)
                 fqq3z2 = fqq3z2 + factor3 * xijkfa*Iqq3z2(i,j,k)
                 fqq3z3 = fqq3z3 + factor3 * xijkfa*Iqq3z3(i,j,k)
                 fqq3z4 = fqq3z4 + factor3 * xijkfa*Iqq3z4(i,j,k)
              END IF
           END DO
        END IF
     END DO
  END DO
  IF (fqq2 < -1.E-50 .AND. fqq3 /= 0.0 .AND. fqq2z /= 0.0 .AND. fqq3z /= 0.0) THEN
     fqqdr = fqq2* (fqq2*fqq2z - 2.0*fqq3*fqq2z+fqq2*fqq3z) /(fqq2-fqq3)**2
     fqqd2= (2.0*fqq2*fqq2z*fqq2z +fqq2*fqq2*fqq2z2  &
          -2.0*fqq2z**2 *fqq3-2.0*fqq2*fqq2z2*fqq3+fqq2*fqq2*fqq3z2)  &
          /(fqq2-fqq3)**2 + fqqdr * 2.0*(fqq3z-fqq2z)/(fqq2-fqq3)
     fqqd3=(2.0*fqq2z**3  +6.0*fqq2*fqq2z*fqq2z2+fqq2*fqq2*fqq2z3  &
          -6.0*fqq2z*fqq2z2*fqq3-2.0*fqq2z**2 *fqq3z  &
          -2.0*fqq2*fqq2z3*fqq3 -2.0*fqq2*fqq2z2*fqq3z  &
          +2.0*fqq2*fqq2z*fqq3z2+fqq2*fqq2*fqq3z3) /(fqq2-fqq3)**2  &
          + 2.0/(fqq2-fqq3)* ( 2.0*fqqd2*(fqq3z-fqq2z)  &
          +     fqqdr*(fqq3z2-fqq2z2) -     fqqdr/(fqq2-fqq3)*(fqq3z-fqq2z)**2 )
     fqqd4=( 12.0*fqq2z**2 *fqq2z2+6.0*fqq2*fqq2z2**2  &
          +8.0*fqq2*fqq2z*fqq2z3+fqq2*fqq2*fqq2z4-6.0*fqq2z2**2 *fqq3  &
          -12.0*fqq2z*fqq2z2*fqq3z -8.0*fqq2z*fqq2z3*fqq3  &
          -2.0*fqq2*fqq2z4*fqq3-4.0*fqq2*fqq2z3*fqq3z  &
          +4.0*fqq2*fqq2z*fqq3z3+fqq2**2 *fqq3z4 ) /(fqq2-fqq3)**2  &
          + 6.0/(fqq2-fqq3)* ( fqqd3*(fqq3z-fqq2z)  &
          -fqqd2/(fqq2-fqq3)*(fqq3z-fqq2z)**2  &
          - fqqdr/(fqq2-fqq3)*(fqq3z-fqq2z)*(fqq3z2-fqq2z2)  &
          + fqqd2*(fqq3z2-fqq2z2) +1.0/3.0*fqqdr*(fqq3z3-fqq2z3) )
     zqq   = fqqdr*eta
     zqqz  = fqqd2*eta + fqqdr
     zqqz2 = fqqd3*eta + 2.0* fqqd2
     zqqz3 = fqqd4*eta + 3.0* fqqd3
  END IF


END SUBROUTINE P_QQ_GROSS

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_DQ_VRABEC_GROSS( zdq, zdqz, zdqz2, zdqz3 )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: zdq, zdqz, zdqz2, zdqz3
  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, m
  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, eij
  REAL                                   :: fdqdr, fdqd2, fdqd3, fdqd4
  REAL                                   :: fdq2, fdq2z, fdq2z2, fdq2z3, fdq2z4
  REAL                                   :: fdq3, fdq3z, fdq3z2, fdq3z3, fdq3z4
  REAL, DIMENSION(nc)                    :: my2dd, myfac, qq2, q_fac
  REAL, DIMENSION(nc,nc)                 :: Idq2, Idq2z, Idq2z2, Idq2z3, Idq2z4
  REAL, DIMENSION(nc,nc)                 :: Idq4, Idq4z, Idq4z2, Idq4z3, Idq4z4
  REAL, DIMENSION(nc,nc,nc)              :: Idq3, Idq3z, Idq3z2, Idq3z3, Idq3z4
  !-----------------------------------------------------------------------------

  zdq   = 0.0
  zdqz  = 0.0
  zdqz2 = 0.0
  zdqz3 = 0.0
  z3 = eta
  DO i = 1, ncomp
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
     myfac(i) = parame(i,3)/t*parame(i,2)**4  *my2dd(i)
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     q_fac(i) = parame(i,3)/t*parame(i,2)**4  *qq2(i)
  END DO


  DO i = 1, ncomp
     DO j = 1, ncomp
        Idq2(i,j)   = 0.0
        Idq4(i,j)   = 0.0
        Idq2z(i,j)  = 0.0
        Idq4z(i,j)  = 0.0
        Idq2z2(i,j) = 0.0
        Idq4z2(i,j) = 0.0
        Idq2z3(i,j) = 0.0
        Idq4z3(i,j) = 0.0
        Idq2z4(i,j) = 0.0
        Idq4z4(i,j) = 0.0
        IF (myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0) THEN
           DO m = 0, 4
              Idq2(i,j)  =Idq2(i,j) + dqp2(i,j,m)*z3**(m+1)
              Idq4(i,j)  =Idq4(i,j) + dqp4(i,j,m)*z3**(m+1)
              Idq2z(i,j) =Idq2z(i,j) +dqp2(i,j,m)*REAL(m+1)*z3**m
              Idq4z(i,j) =Idq4z(i,j) +dqp4(i,j,m)*REAL(m+1)*z3**m
              Idq2z2(i,j)=Idq2z2(i,j)+dqp2(i,j,m)*REAL((m+1)*m)*z3**(m-1)
              Idq4z2(i,j)=Idq4z2(i,j)+dqp4(i,j,m)*REAL((m+1)*m)*z3**(m-1)
              Idq2z3(i,j)=Idq2z3(i,j)+dqp2(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
              Idq4z3(i,j)=Idq4z3(i,j)+dqp4(i,j,m)*REAL((m+1)*m*(m-1)) *z3**(m-2)
              Idq2z4(i,j)=Idq2z4(i,j)+dqp2(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
              Idq4z4(i,j)=Idq4z4(i,j)+dqp4(i,j,m)*REAL((m+1)*m*(m-1)*(m-2)) *z3**(m-3)
           END DO
           DO k = 1, ncomp
              Idq3(i,j,k)   = 0.0
              Idq3z(i,j,k)  = 0.0
              Idq3z2(i,j,k) = 0.0
              Idq3z3(i,j,k) = 0.0
              Idq3z4(i,j,k) = 0.0
              IF (myfac(k) /= 0.0.OR.q_fac(k) /= 0.0) THEN
                 DO m = 0, 4
                    Idq3(i,j,k)  =Idq3(i,j,k) + dqp3(i,j,k,m)*z3**(m+2)
                    Idq3z(i,j,k)=Idq3z(i,j,k)+dqp3(i,j,k,m)*REAL(m+2)*z3**(m+1)
                    Idq3z2(i,j,k)=Idq3z2(i,j,k)+dqp3(i,j,k,m)*REAL((m+2)*(m+1)) *z3**m
                    Idq3z3(i,j,k)=Idq3z3(i,j,k)+dqp3(i,j,k,m)*REAL((m+2)*(m+1)*m) *z3**(m-1)
                    Idq3z4(i,j,k)=Idq3z4(i,j,k)+dqp3(i,j,k,m)  &
                         *REAL((m+2)*(m+1)*m*(m-1)) *z3**(m-2)
                 END DO
              END IF
           END DO

        END IF
     END DO
  END DO

  factor2= -9.0/4.0*PI *rho/z3
  factor3=  PI**2 * (rho/z3)**2

  fdq2   = 0.0
  fdq2z  = 0.0
  fdq2z2 = 0.0
  fdq2z3 = 0.0
  fdq2z4 = 0.0
  fdq3   = 0.0
  fdq3z  = 0.0
  fdq3z2 = 0.0
  fdq3z3 = 0.0
  fdq3z4 = 0.0
  DO i = 1, ncomp
     DO j = 1, ncomp
        IF (myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0) THEN
           xijfa =x(i)*myfac(i) * x(j)*q_fac(j) /sig_ij(i,j)**5
           eij = (parame(i,3)*parame(j,3))**0.5
           fdq2=  fdq2  +factor2*xijfa*(Idq2(i,j)  +eij/t*Idq4(i,j)  )
           fdq2z =fdq2z +factor2*xijfa*(Idq2z(i,j) +eij/t*Idq4z(i,j) )
           fdq2z2=fdq2z2+factor2*xijfa*(Idq2z2(i,j)+eij/t*Idq4z2(i,j))
           fdq2z3=fdq2z3+factor2*xijfa*(Idq2z3(i,j)+eij/t*Idq4z3(i,j))
           fdq2z4=fdq2z4+factor2*xijfa*(Idq2z4(i,j)+eij/t*Idq4z4(i,j))
           DO k = 1, ncomp
              IF (myfac(k) /= 0.0.OR.q_fac(k) /= 0.0) THEN
                 xijkfa=x(i)*x(j)*x(k)/(sig_ij(i,j)*sig_ij(i,k)*sig_ij(j,k))**2  &
                      *( myfac(i)*q_fac(j)*myfac(k) + myfac(i)*q_fac(j)*q_fac(k)*1.193735 )
                 fdq3  =fdq3   + factor3 * xijkfa*Idq3(i,j,k)
                 fdq3z =fdq3z  + factor3 * xijkfa*Idq3z(i,j,k)
                 fdq3z2=fdq3z2 + factor3 * xijkfa*Idq3z2(i,j,k)
                 fdq3z3=fdq3z3 + factor3 * xijkfa*Idq3z3(i,j,k)
                 fdq3z4=fdq3z4 + factor3 * xijkfa*Idq3z4(i,j,k)
              END IF
           END DO
        END IF
     END DO
  END DO
  IF (fdq2 < -1.E-50 .AND. fdq3 /= 0.0 .AND. fdq2z /= 0.0 .AND. fdq3z /= 0.0) THEN
     fdqdr = fdq2* (fdq2*fdq2z - 2.0*fdq3*fdq2z+fdq2*fdq3z) /(fdq2-fdq3)**2
     fdqd2= (2.0*fdq2*fdq2z*fdq2z +fdq2*fdq2*fdq2z2  &
          -2.0*fdq2z**2 *fdq3-2.0*fdq2*fdq2z2*fdq3+fdq2*fdq2*fdq3z2)  &
          /(fdq2-fdq3)**2 + fdqdr * 2.0*(fdq3z-fdq2z)/(fdq2-fdq3)
     fdqd3=(2.0*fdq2z**3  +6.0*fdq2*fdq2z*fdq2z2+fdq2*fdq2*fdq2z3  &
          -6.0*fdq2z*fdq2z2*fdq3-2.0*fdq2z**2 *fdq3z  &
          -2.0*fdq2*fdq2z3*fdq3 -2.0*fdq2*fdq2z2*fdq3z  &
          +2.0*fdq2*fdq2z*fdq3z2+fdq2*fdq2*fdq3z3) /(fdq2-fdq3)**2  &
          + 2.0/(fdq2-fdq3)* ( 2.0*fdqd2*(fdq3z-fdq2z)  &
          +     fdqdr*(fdq3z2-fdq2z2) -     fdqdr/(fdq2-fdq3)*(fdq3z-fdq2z)**2 )
     fdqd4=( 12.0*fdq2z**2 *fdq2z2+6.0*fdq2*fdq2z2**2  &
          +8.0*fdq2*fdq2z*fdq2z3+fdq2*fdq2*fdq2z4-6.0*fdq2z2**2 *fdq3  &
          -12.0*fdq2z*fdq2z2*fdq3z -8.0*fdq2z*fdq2z3*fdq3  &
          -2.0*fdq2*fdq2z4*fdq3-4.0*fdq2*fdq2z3*fdq3z  &
          +4.0*fdq2*fdq2z*fdq3z3+fdq2**2 *fdq3z4 ) /(fdq2-fdq3)**2  &
          + 6.0/(fdq2-fdq3)* ( fdqd3*(fdq3z-fdq2z)  &
          -fdqd2/(fdq2-fdq3)*(fdq3z-fdq2z)**2  &
          - fdqdr/(fdq2-fdq3)*(fdq3z-fdq2z)*(fdq3z2-fdq2z2)  &
          + fdqd2*(fdq3z2-fdq2z2) +1.0/3.0*fdqdr*(fdq3z3-fdq2z3) )
     zdq   = fdqdr*eta
     zdqz  = fdqd2*eta + fdqdr
     zdqz2 = fdqd3*eta + 2.0* fdqd2
     zdqz3 = fdqd4*eta + 3.0* fdqd3
  END IF


END SUBROUTINE P_DQ_VRABEC_GROSS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_POLAR ( k, z3_rk, fdd_rk, fqq_rk, fdq_rk )

  USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: k
  REAL, INTENT(IN)                       :: z3_rk
  REAL, INTENT(OUT)                      :: fdd_rk, fqq_rk, fdq_rk

  !-----------------------------------------------------------------------------
  INTEGER                                :: dipole
  INTEGER                                :: quadrupole
  INTEGER                                :: dipole_quad
  !-----------------------------------------------------------------------------

  fdd_rk = 0.0
  fqq_rk = 0.0
  fdq_rk = 0.0

  dipole      = 0
  quadrupole  = 0
  dipole_quad = 0
  IF ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
  IF ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
  IF ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

  !-----------------------------------------------------------------------------
  ! dipole-dipole term
  !-----------------------------------------------------------------------------
  IF (dipole == 1) THEN

     IF (dd_term == 'GV') CALL PHI_DD_GROSS_VRABEC( k, z3_rk, fdd_rk )
     ! IF (dd_term == 'SF') CALL PHI_DD_SAAGER_FISCHER( k )

     IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! quadrupole-quadrupole term
  !-----------------------------------------------------------------------------
  IF (quadrupole == 1) THEN

     !IF (qq_term == 'SF') CALL PHI_QQ_SAAGER_FISCHER( k )
     IF (qq_term == 'JG') CALL PHI_QQ_GROSS( k, z3_rk, fqq_rk )

     IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

  ENDIF

  !-----------------------------------------------------------------------------
  ! dipole-quadrupole cross term
  !-----------------------------------------------------------------------------
  IF (dipole_quad == 1) THEN

     IF (dq_term == 'VG') CALL PHI_DQ_VRABEC_GROSS( k, z3_rk, fdq_rk )

     IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

  ENDIF

END SUBROUTINE PHI_POLAR


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_DD_GROSS_VRABEC( k, z3_rk, fdd_rk )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: k
  REAL, INTENT(IN)                       :: z3_rk
  REAL, INTENT(IN OUT)                   :: fdd_rk

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, l, m

  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
  REAL                                   :: fdd2, fdd3, fdd2x, fdd3x
  REAL, DIMENSION(nc)                    :: my2dd
  REAL, DIMENSION(nc,nc)                 :: Idd2, Idd4, Idd2x, Idd4x
  REAL, DIMENSION(nc,nc,nc)              :: Idd3, Idd3x
  !-----------------------------------------------------------------------------


  fdd_rk = 0.0
  z3 = eta
  DO i = 1, ncomp
     IF ( uij(i,i) == 0.0 ) write (*,*) 'PHI_DD_GROSS_VRABEC: do not use dimensionless units'
     IF ( uij(i,i) == 0.0 ) stop
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Idd2(i,j)  = 0.0
        Idd4(i,j)  = 0.0
        Idd2x(i,j) = 0.0
        Idd4x(i,j) = 0.0
        IF (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) THEN
           DO m=0,4
              Idd2(i,j)  =Idd2(i,j) + ddp2(i,j,m)*z3**m
              Idd4(i,j)  =Idd4(i,j) + ddp4(i,j,m)*z3**m
              Idd2x(i,j) =Idd2x(i,j)+ ddp2(i,j,m)*REAL(m)*z3**(m-1)*z3_rk
              Idd4x(i,j) =Idd4x(i,j)+ ddp4(i,j,m)*REAL(m)*z3**(m-1)*z3_rk
           END DO
           DO l = 1, ncomp
              Idd3(i,j,l)  = 0.0
              Idd3x(i,j,l) = 0.0
              IF (parame(l,6) /= 0.0) THEN
                 DO m=0,4
                    Idd3(i,j,l) =Idd3(i,j,l) +ddp3(i,j,l,m)*z3**m
                    Idd3x(i,j,l)=Idd3x(i,j,l)+ddp3(i,j,l,m)*REAL(m)*z3**(m-1)*z3_rk
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO

  factor2= -PI
  factor3= -4.0/3.0*PI**2

  fdd2  = 0.0
  fdd3  = 0.0
  fdd2x = 0.0
  fdd3x = 0.0
  DO i = 1, ncomp
     xijfa_x = 2.0*x(i)*rho*uij(i,i)*my2dd(i)*sig_ij(i,i)**3 /t  &
          *uij(k,k)*my2dd(k)*sig_ij(k,k)**3 /t/sig_ij(i,k)**3
     eij = (parame(i,3)*parame(k,3))**0.5
     fdd2x = fdd2x + factor2*xijfa_x*( Idd2(i,k) + eij/t*Idd4(i,k) )
     DO j = 1, ncomp
        IF (parame(i,6) /= 0.0 .AND. parame(j,6) /= 0.0) THEN
           xijfa =x(i)*rho*uij(i,i)*my2dd(i)*sig_ij(i,i)**3 /t  &
                *x(j)*rho*uij(j,j)*my2dd(j)*sig_ij(j,j)**3 /t/sig_ij(i,j)**3
           eij = (parame(i,3)*parame(j,3))**0.5
           fdd2=  fdd2  +factor2*xijfa*(Idd2(i,j) +eij/t*Idd4(i,j) )
           fdd2x =fdd2x +factor2*xijfa*(Idd2x(i,j)+eij/t*Idd4x(i,j))
           !---------------------
           xijkf_x=x(i)*rho*uij(i,i)*my2dd(i)*sig_ij(i,i)**3 /t/sig_ij(i,j)  &
                *x(j)*rho*uij(j,j)*my2dd(j)*sig_ij(j,j)**3 /t/sig_ij(i,k)  &
                *3.0*   uij(k,k)*my2dd(k)*sig_ij(k,k)**3 /t/sig_ij(j,k)
           fdd3x=fdd3x+factor3*xijkf_x*Idd3(i,j,k)
           DO l=1,ncomp
              IF (parame(l,6) /= 0.0) THEN
                 xijkfa= x(i)*rho*uij(i,i)/t*my2dd(i)*sig_ij(i,i)**3   &
                      *x(j)*rho*uij(j,j)/t*my2dd(j)*sig_ij(j,j)**3   &
                      *x(l)*rho*uij(l,l)/t*my2dd(l)*sig_ij(l,l)**3   &
                      /sig_ij(i,j)/sig_ij(i,l)/sig_ij(j,l)
                 fdd3  =fdd3  + factor3 * xijkfa *Idd3(i,j,l)
                 fdd3x =fdd3x + factor3 * xijkfa *Idd3x(i,j,l)
              END IF
           END DO
        END IF
     END DO
  END DO

  IF (fdd2 < -1.E-50 .AND. fdd3 /= 0.0 .AND. fdd2x /= 0.0 .AND. fdd3x /= 0.0)THEN

     fdd_rk = fdd2* (fdd2*fdd2x - 2.0*fdd3*fdd2x+fdd2*fdd3x) / (fdd2-fdd3)**2

  END IF

END SUBROUTINE PHI_DD_GROSS_VRABEC



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_QQ_GROSS( k, z3_rk, fqq_rk )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: k
  REAL, INTENT(IN)                       :: z3_rk
  REAL, INTENT(IN OUT)                   :: fqq_rk

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, l, m

  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
  REAL                                   :: fqq2, fqq3, fqq2x, fqq3x
  REAL, DIMENSION(nc)                    :: qq2
  REAL, DIMENSION(nc,nc)                 :: Iqq2, Iqq4, Iqq2x, Iqq4x
  REAL, DIMENSION(nc,nc,nc)              :: Iqq3, Iqq3x
  !-----------------------------------------------------------------------------


  fqq_rk = 0.0
  z3 = eta
  DO i = 1, ncomp
     IF ( uij(i,i) == 0.0 ) write (*,*) 'PHI_QQ_GROSS: do not use dimensionless units'
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Iqq2(i,j)  = 0.0
        Iqq4(i,j)  = 0.0
        Iqq2x(i,j) = 0.0
        Iqq4x(i,j) = 0.0
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           DO m = 0, 4
              Iqq2(i,j)  = Iqq2(i,j)  + qqp2(i,j,m) * z3**m
              Iqq4(i,j)  = Iqq4(i,j)  + qqp4(i,j,m) * z3**m
              Iqq2x(i,j) = Iqq2x(i,j) + qqp2(i,j,m) * REAL(m)*z3**(m-1)*z3_rk
              Iqq4x(i,j) = Iqq4x(i,j) + qqp4(i,j,m) * REAL(m)*z3**(m-1)*z3_rk
           END DO
           DO l = 1, ncomp
              Iqq3(i,j,l)  = 0.0
              Iqq3x(i,j,l) = 0.0
              IF (parame(l,7) /= 0.0) THEN
                 DO m = 0, 4
                    Iqq3(i,j,l)  = Iqq3(i,j,l)  + qqp3(i,j,l,m)*z3**m
                    Iqq3x(i,j,l) = Iqq3x(i,j,l) + qqp3(i,j,l,m)*REAL(m) *z3**(m-1)*z3_rk
                 END DO
              END IF
           END DO
        END IF
     END DO
  END DO

  factor2= -9.0/16.0*PI
  factor3=  9.0/16.0*PI**2

  fqq2  = 0.0
  fqq3  = 0.0
  fqq2x = 0.0
  fqq3x = 0.0
  DO i = 1, ncomp
     xijfa_x = 2.0*x(i)*rho*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
          *uij(k,k)*qq2(k)*sig_ij(k,k)**5 /t/sig_ij(i,k)**7
     eij = (parame(i,3)*parame(k,3))**0.5
     fqq2x =fqq2x +factor2*xijfa_x*(Iqq2(i,k)+eij/t*Iqq4(i,k))
     DO j=1,ncomp
        IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
           xijfa =x(i)*rho*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
                *x(j)*rho*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7
           eij = (parame(i,3)*parame(j,3))**0.5
           fqq2=  fqq2  +factor2*xijfa*(Iqq2(i,j) +eij/t*Iqq4(i,j) )
           fqq2x =fqq2x +factor2*xijfa*(Iqq2x(i,j)+eij/t*Iqq4x(i,j))
           ! ------------------
           xijkf_x=x(i)*rho*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t/sig_ij(i,j)**3   &
                *x(j)*rho*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,k)**3   &
                *3.0*   uij(k,k)*qq2(k)*sig_ij(k,k)**5 /t/sig_ij(j,k)**3
           fqq3x = fqq3x + factor3*xijkf_x*Iqq3(i,j,k)
           DO l = 1, ncomp
              IF (parame(l,7) /= 0.0) THEN
                 xijkfa=x(i)*rho*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t/sig_ij(i,j)**3   &
                      *x(j)*rho*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,l)**3   &
                      *x(l)*rho*uij(l,l)*qq2(l)*sig_ij(l,l)**5 /t/sig_ij(j,l)**3
                 fqq3  =fqq3  + factor3 * xijkfa *Iqq3(i,j,l)
                 fqq3x =fqq3x + factor3 * xijkfa *Iqq3x(i,j,l)
              END IF
           END DO
        END IF
     END DO
  END DO

  IF (fqq2 < -1.E-50 .AND. fqq3 /= 0.0 .AND. fqq2x /= 0.0 .AND. fqq3x /= 0.0) THEN
     fqq_rk = fqq2* (fqq2*fqq2x - 2.0*fqq3*fqq2x+fqq2*fqq3x) / (fqq2-fqq3)**2
  END IF

END SUBROUTINE PHI_QQ_GROSS

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_DQ_VRABEC_GROSS( k, z3_rk, fdq_rk )

  USE PARAMETERS, ONLY: PI, KBOL
  USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
  USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: k
  REAL, INTENT(IN)                       :: z3_rk
  REAL, INTENT(IN OUT)                   :: fdq_rk

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, l, m

  REAL                                   :: factor2, factor3, z3
  REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
  REAL                                   :: fdq2, fdq3, fdq2x, fdq3x
  REAL, DIMENSION(nc)                    :: my2dd, myfac, qq2, q_fac
  REAL, DIMENSION(nc,nc)                 :: Idq2, Idq4, Idq2x, Idq4x
  REAL, DIMENSION(nc,nc,nc)              :: Idq3, Idq3x
  !-----------------------------------------------------------------------------

  fdq_rk = 0.0
  z3 = eta
  DO i=1,ncomp
     my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*KBOL*mseg(i)*sig_ij(i,i)**3 *1.E-30)
     myfac(i) = parame(i,3)/t*parame(i,2)**4 *my2dd(i)
     qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
     q_fac(i) = parame(i,3)/t*parame(i,2)**4 *qq2(i)
  END DO

  DO i = 1, ncomp
     DO j = 1, ncomp
        Idq2(i,j)  = 0.0
        Idq4(i,j)  = 0.0
        Idq2x(i,j) = 0.0
        Idq4x(i,j) = 0.0
        DO m = 0, 4
           Idq2(i,j)  = Idq2(i,j)  + dqp2(i,j,m)*z3**m
           Idq4(i,j)  = Idq4(i,j)  + dqp4(i,j,m)*z3**m
           Idq2x(i,j) = Idq2x(i,j) + dqp2(i,j,m)*REAL(m)*z3**(m-1) *z3_rk
           Idq4x(i,j) = Idq4x(i,j) + dqp4(i,j,m)*REAL(m)*z3**(m-1) *z3_rk
        END DO
        DO l = 1, ncomp
           Idq3(i,j,l)  = 0.0
           Idq3x(i,j,l) = 0.0
           DO m = 0, 4
              Idq3(i,j,l) =Idq3(i,j,l) +dqp3(i,j,l,m)*z3**m
              Idq3x(i,j,l)=Idq3x(i,j,l)+dqp3(i,j,l,m)*REAL(m)*z3**(m-1)*z3_rk
           END DO
        END DO
     END DO
  END DO

  factor2= -9.0/4.0*PI
  factor3=  PI**2

  fdq2  = 0.0
  fdq3  = 0.0
  fdq2x = 0.0
  fdq3x = 0.0
  DO i = 1, ncomp
     xijfa_x = x(i)*rho*( myfac(i)*q_fac(k) + myfac(k)*q_fac(i) ) / sig_ij(i,k)**5
     eij = (parame(i,3)*parame(k,3))**0.5
     fdq2x =fdq2x +factor2*xijfa_x*(Idq2(i,k)+eij/t*Idq4(i,k))
     DO j=1,ncomp
        xijfa =x(i)*rho*myfac(i) * x(j)*rho*q_fac(j) /sig_ij(i,j)**5
        eij = (parame(i,3)*parame(j,3))**0.5
        fdq2=  fdq2  +factor2*xijfa*(Idq2(i,j)  +eij/t*Idq4(i,j) )
        fdq2x =fdq2x +factor2*xijfa*(Idq2x(i,j) +eij/t*Idq4x(i,j))
        !---------------------
        xijkf_x=x(i)*rho*x(j)*rho/(sig_ij(i,j)*sig_ij(i,k)*sig_ij(j,k))**2  &
             *( myfac(i)*q_fac(j)*myfac(k) + myfac(i)*q_fac(k)*myfac(j)  &
             + myfac(k)*q_fac(i)*myfac(j) +myfac(i)*q_fac(j)*q_fac(k)*1.1937350  &
             +myfac(i)*q_fac(k)*q_fac(j)*1.193735  &
             +myfac(k)*q_fac(i)*q_fac(j)*1.193735      )
        fdq3x = fdq3x + factor3*xijkf_x*Idq3(i,j,k)
        DO l = 1, ncomp
           xijkfa=x(i)*rho*x(j)*rho*x(l)*rho/(sig_ij(i,j)*sig_ij(i,l)*sig_ij(j,l))**2  &
                *( myfac(i)*q_fac(j)*myfac(l)  &
                +myfac(i)*q_fac(j)*q_fac(l)*1.193735 )
           fdq3  =fdq3  + factor3 * xijkfa *Idq3(i,j,l)
           fdq3x =fdq3x + factor3 * xijkfa *Idq3x(i,j,l)
        END DO
     END DO
  END DO

  IF (fdq2 < -1.E-50 .AND. fdq3 /= 0.0 .AND. fdq2x /= 0.0 .AND. fdq3x /= 0.0)THEN

     fdq_rk = fdq2* (fdq2*fdq2x - 2.0*fdq3*fdq2x+fdq2*fdq3x) / (fdq2-fdq3)**2

  END IF

END SUBROUTINE PHI_DQ_VRABEC_GROSS

END MODULE EOS_POLAR
