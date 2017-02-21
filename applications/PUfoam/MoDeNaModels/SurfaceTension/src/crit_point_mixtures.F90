!>This file contains the subroutines which determine the critical point of 
!!a mixture.


SUBROUTINE CRIT_POINT_MIX(tc,user)

!PETSc module
USE PetscManagement

!VLE and DFT modules
USE BASIC_VARIABLES

IMPLICIT NONE

#include <finclude/petscsys.h>

type (userctx)         user 
REAL            :: tc
!local
REAL            :: tc_L, tc_V, xi_save(np,nc),p_save
PetscErrorCode ierr


!!J:save value of pressure because it is changed during the calculation
p_save = p

! ---------------------------------------------------------------------
 ! determine the critical temp. to xi of liquid and to xi of vapor
 ! ---------------------------------------------------------------------


 xi_save(:,:) = xi(:,:)
 dense(1) = 0.15
 
 !call MPI_Barrier(PETSC_COMM_WORLD,ierr)

 !only proc 0 reads in the estimate and then sends it to all other procs using MPI_Bcast
!  IF(user%rank == 0) THEN
!   WRITE (*,*) 'provide an estimate of the crit. Temp. of the mixture'
!   READ (*,*) t
!  END IF

 t = 1.2 * t  !initial estimate of critical temperature

 CALL MPI_Bcast(t,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr)

 xiF(1:ncomp) = xi_save(1,1:ncomp)
 CALL Heidemann_Khalil
 tc_L = t
!  IF(user%rank == 0) THEN
!    WRITE (*,*) 'critical temperature to xi_liquid',tc_L
!  END IF

 xiF(1:ncomp) = xi_save(2,1:ncomp)
 ! dense(1) = 0.15
 ! WRITE (*,*) 'provide an estimate of the crit. Temp. of the mixture'
 ! READ (*,*) t
 CALL Heidemann_Khalil
 tc_V = t
!  IF(user%rank == 0) THEN
!    WRITE (*,*) 'critical temperature to xi_vapor ',tc_V
!  END IF

 ! tc_L = 600.0
 ! WRITE (*,*) 'I have tentitativly set tc=700 '
 ! pause


 !tc = ( tc_L + tc_V ) / 2.0
 tc = tc_L
 xi(:,:) = xi_save(:,:)
 densta(1:nphas) = val_conv(1:nphas)
 dense(1:nphas) = val_conv(1:nphas)
 t = val_conv(3)
 IF(user%rank == 0) THEN
  WRITE (*,*) 'estimate of critical temperature:',tc,'K'
  WRITE (*,*) ' '
 END IF

!!J: set pressure to its regular value
p = p_save


END SUBROUTINE CRIT_POINT_MIX









!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Heidemann_Khalil
!
! This subroutine ....
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE Heidemann_Khalil
!
 USE PARAMETERS, ONLY: PI
 USE BASIC_VARIABLES
 USE Module_Heidemann_Khalil
 USE Solve_NonLin
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTERFACE
   SUBROUTINE Heidemann_Khalil_obj ( iter_no, y, residu, dummy )
     INTEGER, INTENT(IN)                :: iter_no
     REAL, INTENT(IN)                   :: y(iter_no)
     REAL, INTENT(OUT)                  :: residu(iter_no)
     INTEGER, INTENT(IN OUT)            :: dummy
   END SUBROUTINE Heidemann_Khalil_obj
 END INTERFACE
!
 INTEGER                                :: info
 REAL, ALLOCATABLE                      :: y(:),diag(:),residu(:)

 INTEGER                                :: i, count, nphas_save
 REAL                                   :: error0, error1, dense0, dfdx, delta_rho
 CHARACTER  (LEN=2)                     :: ensemble_save
! ----------------------------------------------------------------------

 ensemble_save =  ensemble_flag
 ensemble_flag = 'tv'

 nphas_save = nphas
 nphas = 1

 ! xiF(2) = 0.5 * ( 0.71928411 + 0.72025642 )
 ! xiF(1) = 1.0 - xiF(2)
 ! dense(1) = 0.5 * ( 0.159315 + 0.158817 ) + 0.02
 ! t = 500.0

 dense0 = dense(1)

 info = 1
 n_unkw = ncomp + 1
 acc_a  = 5.E-8
 step_a = 1.E-8


 ALLOCATE( y(n_unkw), diag(n_unkw), residu(n_unkw) )
 DO i = 1,ncomp
    y(i) = 1.0 / SQRT( REAL( ncomp ) ) 
 END DO
 y(ncomp+1) = t
 count = 0
 error0 = 1.0

 DO WHILE ( ABS( error0) > 0.001 .AND. count < 20 )
    count = count + 1

    dense(1) = dense0 + 0.0001
    CALL hbrd (Heidemann_Khalil_obj, n_unkw, y, residu, step_a, acc_a, info, diag)
    error1 = error_condition2
    IF (SUM( ABS( residu(1:n_unkw) ) ) > 1.E-5) write (*,*) 'caution: error 1st inner loop', SUM( ABS( residu(1:n_unkw) ) )

    dense(1) = dense0
    CALL hbrd (Heidemann_Khalil_obj, n_unkw, y, residu, step_a, acc_a, info, diag)
    IF (SUM( ABS( residu(1:n_unkw) ) ) > 1.E-5) write (*,*) 'caution: error 2nd inner loop', SUM( ABS( residu(1:n_unkw) ) )
    error0 = error_condition2

    ! write (*,'(a,4G18.10)') ' t, p, eta error', t, p, dense(1), error0
    ! pause
    dfdx = ( error1 - error0 ) / 0.0001
    delta_rho = MIN( error0 / dfdx, 0.02)
    delta_rho = MAX( delta_rho, -0.02)
    dense0 = dense0 - delta_rho
    dense(1) = dense0

 END DO

 !tc = t
 !pc = p

 DEALLOCATE( y, diag, residu )

 ensemble_flag = ensemble_save
 nphas = nphas_save
 IF ( ABS( error_condition2 ) > 1.E-1 ) write (*,*) 'caution: error outer loop', ABS( error_condition2 )

END SUBROUTINE Heidemann_Khalil


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Heidemann_Khalil
!
! This subroutine ....
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE Heidemann_Khalil_obj ( iter_no, y, residu, dummy )
!
 USE PARAMETERS, ONLY: PI
 USE BASIC_VARIABLES
 USE Module_Heidemann_Khalil
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: iter_no
 REAL, INTENT(IN)                       :: y(iter_no)
 REAL, INTENT(OUT)                      :: residu(iter_no)
 INTEGER, INTENT(IN OUT)                :: dummy
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k
 REAL, DIMENSION(nc)                    :: dn
 REAL, DIMENSION(nc)                    :: dhs, rhoi00, rhoi0
 REAL                                   :: rho, d_rho
 ! REAL                                 :: lnphi(np,nc)
 REAL                                   :: qij(nc,nc), qij0(nc,nc), qijk(nc,nc,nc)
 ! CHARACTER (LEN=3)                      :: char_len
! ----------------------------------------------------------------------

 dn(1:ncomp) = y(1:ncomp)
 t           = y(ncomp+1)
 !dense(1)    = y(ncomp+2)

 dhs(1:ncomp) = parame(1:ncomp,2) * ( 1.0 - 0.12 *EXP( -3.0*parame(1:ncomp,3)/t ) )
 rho = dense(1) / SUM( PI/6.0*xiF(1:ncomp)*parame(1:ncomp,1)* dhs(1:ncomp)**3 )
 rhoi00(1:ncomp) = xiF(1:ncomp)*rho

 d_rho = 0.000001

 DO k = 1, ncomp

    rhoi0(1:ncomp) = rhoi00(1:ncomp)
    rhoi0(k) = rhoi00(k) + d_rho
    CALL qij_matrix ( rhoi0, dhs, d_rho, qij )
    qij0(:,:) = qij(:,:)

    rhoi0(1:ncomp) = rhoi00(1:ncomp)
    CALL qij_matrix ( rhoi0, dhs, d_rho, qij )

    DO i = 1, ncomp
       DO j = 1, ncomp
          qijk(i,j,k) = ( qij0(i,j) - qij(i,j) ) / d_rho
          ! write(*,*) i,j,k,qijk(i,j,k)
       END DO
    END DO

 END DO

 ! write (*,'(a,4G18.10)') 'det',qij(2,2)*qij(1,1) - qij(2,1)*qij(1,2)
 ! write (*,'(a,3G18.10)') ' t,p,eta ', t, p, dense(1)
 DO j = 1, ncomp
    residu(j) = SUM( qij(1:ncomp,j)*dn(1:ncomp) )
 END DO
 residu(ncomp+1) = 1.0 - SUM( dn(1:ncomp)*dn(1:ncomp) )

 error_condition2 = 0.0
 DO k = 1, ncomp
    DO i = 1, ncomp
       DO j = 1, ncomp
          error_condition2 = error_condition2 + qijk(i,j,k) * dn(i) * dn(j) * dn(k)
       END DO
    END DO
 END DO
 error_condition2 = error_condition2 * 1.E-4     ! the values are scaled down to prevent numerical dominance of this error
 !residu(ncomp+2) = error_condition2

 !write (char_len,'(I3)') ncomp+2
 !write (*,'(a,'//char_len//'G18.10)') ' error',residu(1:ncomp+1)
 !write (*,*) ' '
 !pause

END SUBROUTINE Heidemann_Khalil_obj


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Heidemann_Khalil
!
! This subroutine calculates the spinodal for a binary mixture to given
! T, p and for a given starting value of the density vector rhoi
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE binary_tp_spinodal ( rhoi, p_sp )
!
 USE PARAMETERS, ONLY: PI
 USE BASIC_VARIABLES
 USE Module_Heidemann_Khalil
 USE Solve_NonLin
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTERFACE
   SUBROUTINE binary_spinodal_obj ( iter_no, y, residu, dummy )
     INTEGER, INTENT(IN)                :: iter_no
     REAL, INTENT(IN)                   :: y(iter_no)
     REAL, INTENT(OUT)                  :: residu(iter_no)
     INTEGER, INTENT(IN OUT)            :: dummy
   END SUBROUTINE binary_spinodal_obj
 END INTERFACE
 
 REAL, dimension(nc)                    :: rhoi
 REAL                                   :: p_sp
! ----------------------------------------------------------------------

 INTEGER                                :: info
 REAL, ALLOCATABLE                      :: y(:), diag(:), residu(:)

 INTEGER                                :: i, nphas_save
 CHARACTER  (LEN=2)                     :: ensemble_save
! ----------------------------------------------------------------------

 ensemble_save =  ensemble_flag
 ensemble_flag = 'tv'

 nphas_save = nphas
 nphas = 1

 info = 1
 n_unkw = ncomp + 2
 acc_a  = 5.E-8
 step_a = 1.E-8


 ALLOCATE( y(n_unkw), diag(n_unkw), residu(n_unkw) )
 DO i = 1, ncomp
    y(i) = 1.0 / SQRT( REAL( ncomp ) ) 
 END DO
 DO i = 1, ncomp
    y( ncomp + i ) = rhoi( i ) 
 END DO

 CALL hbrd (binary_spinodal_obj, n_unkw, y, residu, step_a, acc_a, info, diag)

 DO i = 1, ncomp
    rhoi( i ) = y( ncomp + i )
 END DO
 write (*,*) 'info',info
 write (*,*) 'rhoi',rhoi(1:ncomp)
 !write (*,*) 'eta',PI / 6.0 * sum( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
 write (*,*) 'x',rhoi(1:ncomp) / sum( rhoi(1:ncomp) )

 DEALLOCATE( y, diag, residu )

 ensemble_flag = ensemble_save
 nphas = nphas_save

END SUBROUTINE binary_tp_spinodal


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Heidemann_Khalil
!
! This subroutine ....
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE binary_spinodal_obj ( iter_no, y, residu, dummy )
!
 USE PARAMETERS, ONLY: PI
 USE BASIC_VARIABLES
 USE Module_Heidemann_Khalil
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: iter_no
 REAL, INTENT(IN)                       :: y(iter_no)
 REAL, INTENT(OUT)                      :: residu(iter_no)
 INTEGER, INTENT(IN OUT)                :: dummy
!
! ----------------------------------------------------------------------
 INTEGER                                :: j, k
 REAL                                   :: p_calculated, zges, p_sp
 REAL                                   :: d_rho
 REAL, DIMENSION(nc)                    :: dn
 REAL, DIMENSION(nc)                    :: dhs, rhoi
 REAL, DIMENSION(nc,nc)                 :: qij
! ----------------------------------------------------------------------

 dn(1:ncomp)   = y(1:ncomp)
 rhoi(1:ncomp) = y( (ncomp+1) : (ncomp+ncomp) )
 dhs(1:ncomp) = parame(1:ncomp,2) * ( 1.0 - 0.12 *EXP( -3.0*parame(1:ncomp,3)/t ) )  ! is this calc. needed?

 call p_calc ( p_calculated, zges )

 d_rho = 0.000001

 DO k = 1, ncomp

    CALL qij_matrix ( rhoi, dhs, d_rho, qij )

 END DO

 !write (*,'(a,4G18.10)') 'det',qij(2,2)*qij(1,1) - qij(2,1)*qij(1,2)
 write (*,'(a,3G18.10)') ' t,p,eta ', t, p, dense(1)
 DO j = 1, ncomp
    residu(j) = SUM( qij(1:ncomp,j)*dn(1:ncomp) )
 END DO
 residu(ncomp+1) = 1.0 - SUM( dn(1:ncomp)*dn(1:ncomp) )
 write (*,*) 'p_sp has to be handed over to the obj fct properly!'
 write (*,*) 'Can I simply set p = p_sp? In other words is p altered during the calculation?'
 stop 5
 residu(ncomp+2) = p_sp - p_calculated

 write (*,'(a,4G18.10)') ' error',residu(1:4)
 write (*,*) ' '
 pause

END SUBROUTINE binary_spinodal_obj


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Heidemann_Khalil
!
! This subroutine ....
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE qij_matrix ( rhoi0, dhs, d_rho, qij )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE BASIC_VARIABLES
 USE EOS_VARIABLES, ONLY: pges
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, DIMENSION(nc), INTENT(IN)        :: rhoi0
 REAL, DIMENSION(nc), INTENT(IN)        :: dhs
 REAL, INTENT(IN)                       :: d_rho
 REAL, DIMENSION(nc,nc), INTENT(OUT)    :: qij
!
! ----------------------------------------------------------------------
 INTEGER                                 :: i
 REAL, DIMENSION(nc)                     :: rhoi, lnf0, lnf1, lnf2
 REAL                                    :: lnphi(np,nc), zges
! ----------------------------------------------------------------------


 DO i = 1, ncomp
    rhoi(1:ncomp) = rhoi0(1:ncomp)
    rhoi(i) = rhoi0(i) + d_rho
    dense(1) = SUM( PI/6.0*rhoi(1:ncomp)*parame(1:ncomp,1)* dhs(1:ncomp)**3 )
    densta(1) = dense(1)
    xi(1,1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
    CALL FUGACITY ( lnphi )
    p = pges
    zges = (pges * 1.d-30) / ( KBOL*t*SUM(rhoi(1:ncomp)) )
    lnf1(1:ncomp) = lnphi(1,1:ncomp) + LOG( rhoi(1:ncomp) )

    IF (rhoi0(i) - d_rho > 0.0 ) THEN
       rhoi(1:ncomp) = rhoi0(1:ncomp)
       rhoi(i) = rhoi0(i) - d_rho
       dense(1) = SUM( PI/6.0*rhoi(1:ncomp)*parame(1:ncomp,1)* dhs(1:ncomp)**3 )
       densta(1) = dense(1)
       xi(1,1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
       CALL FUGACITY ( lnphi )
       p = pges
       zges = (pges * 1.d-30) / ( KBOL*t*SUM(rhoi(1:ncomp)) )
       lnf2(1:ncomp) = lnphi(1,1:ncomp) + LOG( rhoi(1:ncomp) )
    END IF

    rhoi(1:ncomp) = rhoi0(1:ncomp)
    dense(1) = SUM( PI/6.0*rhoi(1:ncomp)*parame(1:ncomp,1)* dhs(1:ncomp)**3 )
    densta(1) = dense(1)
    xi(1,1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
    CALL FUGACITY ( lnphi )
    p = pges
    zges = (pges * 1.d-30) / ( KBOL*t*SUM(rhoi(1:ncomp)) )
    lnf0(1:ncomp) = lnphi(1,1:ncomp) + LOG( rhoi(1:ncomp) )

    IF (rhoi0(i) - d_rho > 0.0 ) THEN
       qij(i,1:ncomp) = ( lnf1(1:ncomp) - lnf2(1:ncomp) ) / (2.0*d_rho)  ! qij = d(F/V) / (d_rho_i*d_rho_j)
    ELSE
       qij(i,1:ncomp) = ( lnf1(1:ncomp) - lnf0(1:ncomp) ) / d_rho        ! qij = d(F/V) / (d_rho_i*d_rho_j)
    END IF
    ! write (*,*) i,1,qij(i,1)
    ! write (*,*) i,2,qij(i,2)
 END DO

END SUBROUTINE qij_matrix

