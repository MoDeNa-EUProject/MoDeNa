!> This file contains the subroutines which perform and control the 
!! phase equilibrium calculation.





!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE start_var
!!
!! This subroutine generates a converged solution for binary systems
!! or performes a flash calculation for mixtues. This routine is a
!! fairly weak point of the program.
!!
!! IF a polymer is considered, starting values for mole fractions
!! are determined from the SUBROUTINGE POLY_STA_VAR (see below). The
!! polymer needs to be placed as component 1 (first line) in INPUT
!! file.
!!
!! A phase equilib. iteration is started at the end of this routine.
!! If no solution is found (converg=0), the program will stop within
!! this routine.
!!
!! Currently, this routine assumes two-phase equilibrium and derives
!! starting values (xi,density) only for two phases.
!!
!! Prerequisites are:
!! SUBROUTINE INPUT needs to be called prior to this routine, because
!! all pure comp. parameters as well as (T,P,kij) need to be in place.
!! Also, the variable to be iterated "it(i)" and the variables to be
!! calculated through the summation relation "sum_rel(i)" have to be
!! defined.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 SUBROUTINE start_var(converg,user)

!Petsc modules
 USE PetscManagement

!VLE modules
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
type (userctx)         user 
 INTEGER, INTENT(IN OUT)                 :: converg
!
! ----------------------------------------------------------------------
 INTEGER                                 :: ph, i, k
 INTEGER                                 :: ncompsav, n_unkwsav, ph_split
 LOGICAL                                 :: lle_check, flashcase, renormalize
 REAL                                    :: den1, den2, x_1, x_2
 CHARACTER (LEN=50)                      :: filename
! ----------------------------------------------------------------------

 converg = 0
  
  renormalize = .false.       ! for renormalization group theory (RGT)
  IF (num == 2) renormalize = .true.
  IF (num == 2) num = 0       ! if RGT: initial phase equilibr. is for non-renormalized model

  flashcase = .false.         ! .true. when a specific feed conc. xif is given
  IF (xif(1) /= 0.0) flashcase = .true.

  lle_check = .true.
    
    DO i=1,ncomp              ! setting mole-fractions for the case that
                              ! anything goes wrong in the coming routines
      xi(1,i) = 1.0 / REAL(ncomp)
      xi(2,i) = 1.0 / REAL(ncomp)
    END DO
    

    ! ------------------------------------------------------------------
    ! determine an initial conc. (phase 1) that will phase split
    ! ------------------------------------------------------------------
    IF( ncomp == 2 .AND. .NOT.flashcase ) THEN
      CALL vle_min( lle_check )
      !IF(user%rank == 0) THEN
      !WRITE(*,*)' INITIAL FEED-COMPOSITION',(xi(1,i), i=1,ncomp),converg
      !END IF   
   END IF
    
    ! ------------------------------------------------------------------
    ! perform a phase stability test
    ! ------------------------------------------------------------------
    ph_split = 0
    CALL phase_stability ( .false., flashcase, ph_split,user )
    !IF(user%rank == 0) THEN
    ! write (*,*) 'stability analysis I indicates phase-split is:',ph_split
    !END IF

    ! ------------------------------------------------------------------
    ! determine species i, for which x(i) is calc from summation relation
    ! ------------------------------------------------------------------
    CALL select_sum_rel (1,0,1)         ! synthax (m,n,o): phase m
                                        ! exclude comp. n
                                        ! assign it(o) and higher
    CALL select_sum_rel (2,0,2)         ! for ncomp>=3, the quantities
                                        ! to be iterated will be overwritten

    ! ------------------------------------------------------------------
    ! if 2 phases (VLE)
    ! ------------------------------------------------------------------
    IF (ph_split == 1) THEN

      ! ---  perform tangent plane minimization ------------------------
      CALL tangent_plane
      ph_split = 0

      ! --- determine, for which substance summation relation is used --
      IF (flashcase) THEN
        CALL determine_flash_it2
      ELSE
        CALL select_sum_rel (1,0,1)
        CALL select_sum_rel (2,0,2)
      END IF

      ! --- do full phase equilibr. calculation ------------------------
      n_unkw = ncomp       ! number of quantities to be iterated
      CALL objective_ctrl (converg)
      !IF(user%rank == 0) THEN
      ! IF (converg == 1 ) write (*,*) ' converged (maybe a VLE)',dense(1),dense(2)
      !END IF
    END IF

    ! ------------------------------------------------------------------
    ! test for LLE
    ! ------------------------------------------------------------------
    ph_split = 0

    IF (lle_check) CALL phase_stability (lle_check,flashcase,ph_split)
    !IF(user%rank == 0) THEN
    ! IF (lle_check) write (*,*) 'stability analysis II, phase-split is:',ph_split
    !END IF

    ! ------------------------------------------------------------------
    ! if two phases (LLE)
    ! ------------------------------------------------------------------
    IF (ph_split == 1) THEN
     ! IF(user%rank == 0) THEN
     !  write (*,*) ' LLE-stability test indicates 2 phases (VLE or LLE)'
     ! END IF

      ! ---  perform tangent plane minimization ------------------------
      IF (flashcase) CALL select_sum_rel (1,0,1)
      IF (flashcase) CALL select_sum_rel (2,0,2)

      CALL tangent_plane

      ! --- determine, for which substance summation relation ----------
      IF (flashcase) THEN
        CALL determine_flash_it2
      ELSE
        CALL select_sum_rel (1,0,1)
        CALL select_sum_rel (2,0,2)
      END IF

      ! --- do full phase equilibr. calculation ------------------------
      n_unkw = ncomp       ! number of quantities to be iterated
      val_conv(2) = 0.0
      CALL objective_ctrl (converg)
      !IF(user%rank == 0) THEN
      ! IF (converg == 1 ) write (*,*) ' converged (maybe an LLE)',dense(1),dense(2)
      !END IF
    END IF

    ! ------------------------------------------------------------------
    ! equilibr. calc. converged: set initial var. for further calc.
    ! ------------------------------------------------------------------
    IF (converg == 1) THEN
        val_init = val_conv
        DO ph = 1,nphas
          DO i = 1,ncomp
            xi(ph,i) = EXP( val_conv(4+i+(ph-1)*ncomp) )
          END DO
        END DO
        dense(1:2) = val_conv(1:2)
    ELSE
      !IF(user%rank == 0) THEN
      !WRITE (*,*) 'Surface Tension Code, Phase equilibrium calculation: NO SOLUTION FOUND FOR THE STARTING VALUES'
      !END IF
      !STOP 5
    END IF

END SUBROUTINE start_var




!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE objective_ctrl
!!
!! This subroutine controls the iso-fugacity iteration. It uses
!! the variables defined in the array "val_init". If successfull,
!! the converged values are written to "val_conv", and the flag
!! converg is set to 1.
!! See also above desciption for subroutine PHASE_EQUILIB
!! This routine calls SUBROUTINE HYBRID, which is a solver (modified
!! POWELL HYBRID METHOD). HYBRID is freely available for non-commercial
!! applications. HYBRID requires three definitions:
!! 1.the number of equations to be solved (=No. of variables to be
!!   iterated). The appropriate parameter is: "n_unkw"
!! 2.the equations to be iterated, they are here gathered in the SUB-
!!   ROUTINE OBJEC_FCT (see below). Since HYBRID is a root finder,
!!   these objective functions are iterated to be zero (essentially,
!!   OBJEC_FCT contains the iso-fugacity relation.
!! 3.an array of variables is required, containing the quatities to be
!!   iterated. This array is termed "y(i)"
!!
!!   INPUT VARIABLES:
!! val_init(i)   array containing (densities,T,P,lnx's) serving as
!!               starting values for the phase equilibrium calculation
!! it(i)         contains the information, which variable is  deter-
!!               mined iteratively. For syntax refer e.g.to SUB BINMIX.
!! sum_rel(i)    indicates, which mole fraction is determined from the
!!               summation relation sum(xi)=1
!!
!!   OUTPUT VARIABLES:
!! val_conv(i)   array containing the converged system variables
!!               analogous to "val_init"
!! converg     0 if no convergence achieved, 1 if converged solution
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 SUBROUTINE objective_ctrl (converg)
!
 USE BASIC_VARIABLES
 USE Solve_NonLin
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(OUT)                   :: converg
!
! ----------------------------------------------------------------------
 INTERFACE
   SUBROUTINE objec_fct ( iter_no, y, residu, dummy )
     INTEGER, INTENT(IN)                :: iter_no
     REAL, INTENT(IN)                   :: y(iter_no)
     REAL, INTENT(OUT)                  :: residu(iter_no)
     INTEGER, INTENT(IN OUT)            :: dummy
   END SUBROUTINE objec_fct
 END INTERFACE
!
 INTEGER                                :: info,k,posn,i
 INTEGER, PARAMETER                     :: mxr = nc*(nc+1)/2
 REAL, ALLOCATABLE                      :: y(:),diag(:),residu(:)
 REAL                                   :: x_init, x_solut, r_diff1, r_diff2, totres
 REAL                                   :: r_thrash, x_thrash
 CHARACTER (LEN=2)                      :: compon
 LOGICAL                                :: convergence
! ----------------------------------------------------------------------

info=1

ALLOCATE( y(n_unkw), diag(n_unkw), residu(n_unkw) )

IF (num == 0) acc_a  = 1.E-7
IF (num == 0) step_a = 2.E-8
IF (num == 1) acc_a  = 1.E-7
IF (num == 1) step_a = 2.E-8
IF (num == 2) acc_a  = 5.E-7
IF (num == 2) step_a = 1.E-7

posn = 0
DO i = 1,n_unkw
  posn = posn + 1
  IF (it(i) == 't') y(posn) = val_init(3)
  IF (it(i) == 'p') y(posn) = val_init(4)
  IF (it(i) == 'lnp') y(posn) = LOG( val_init(4) )
  IF (it(i) == 'fls') y(posn) = alpha
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') y(posn) = val_init(4+k)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') y(posn) = val_init(4+ncomp+k)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') y(posn) = val_init(4+ncomp+ncomp+k)
END DO

CALL init_vars

x_init  = 0.0
DO i = 1,ncomp
  IF (lnx(1,i) /= 0.0 .AND. lnx(2,i) /= 0.0) THEN
    x_init  = x_init  + ABS( 1.0 - lnx(2,i)/lnx(1,i) )
  ELSE
    x_init  = x_init  + ABS( 1.0 - EXP(lnx(2,i))/EXP(lnx(1,i)) )
  END IF
END DO

CALL hbrd (objec_fct, n_unkw, y, residu, step_a, acc_a, info, diag)

x_solut = 0.0
DO i = 1,ncomp
    IF ( lnx(1,i) /= 0.0 .AND. lnx(2,i) /= 0.0 ) THEN
      x_solut = x_solut + ABS( 1.0 - lnx(2,i)/lnx(1,i) )
    ELSE
      IF (lnx(1,i) < 1E300 .AND. lnx(1,i) > -1.E300 )  &
          x_solut = x_solut + ABS( 1.0 - EXP(lnx(2,i))/EXP(lnx(1,i)) )
    END IF
END DO
r_diff1 = ABS( 1.0 - dense(1)/dense(2) )
IF ( val_conv(2) > 0.0 ) THEN
    r_diff2 = ABS( 1.0 - val_conv(1)/val_conv(2) )
ELSE
    r_diff2 = 0.0
END IF

totres = SUM( ABS( residu(1:n_unkw) ) )

r_thrash = 0.0005
x_thrash = 0.0005
if (num > 0 ) r_thrash = r_thrash * 10.0
if (num > 0 ) x_thrash = x_thrash * 100.0

convergence = .true.

IF ( info >= 2 ) convergence = .false.
IF ( ABS( 1.0- dense(1)/dense(2) ) < r_thrash .AND. x_solut < x_thrash ) THEN
  IF ( x_init > 0.050 ) convergence = .false.
  IF ( ( ABS( 1.0- dense(1)/dense(2) ) + x_solut ) < 1.E-7 ) convergence = .false.
ENDIF
IF ( r_diff2 /= 0.0 .AND. r_diff2 > (4.0*r_diff1) .AND. bindiag == 1 ) convergence = .false.
IF ( ncomp == 1 .AND. totres > 100.0*acc_a ) convergence = .false.
IF ( totres > 1000.0*acc_a ) convergence = .false.
IF ( ncomp == 1 .AND. r_diff1 < 1.d-5      ) convergence = .false.

IF ( convergence ) THEN
    converg = 1
    ! write (*,*) residu(1),residu(2)
    CALL converged
    IF (num <= 1) CALL enthalpy_etc
ELSE
    converg = 0
END IF

DEALLOCATE( y, diag, residu )

END SUBROUTINE objective_ctrl



!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE objec_fct
!!
!! This subroutine contains the equations to be solved numerically
!! (iso-fugacity: fi'-fi''=0) as well as other dependent equations,
!! which can be solved analytically, namely the summation relation
!! xi=1-sum(xj) or the condition of equal charge for electrolyte
!! solutions.
!! This subroutine is required and controlled by the solver HBRD !
!! HBRD varies the variables "y(i)" and eveluates the result of
!! these changes from this routine.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 SUBROUTINE objec_fct ( iter_no, y, residu, dummy )
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES, ONLY: density_error
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: iter_no
 REAL, INTENT(IN)                       :: y(iter_no)
 REAL, INTENT(OUT)                      :: residu(iter_no)
 INTEGER, INTENT(IN OUT)                :: dummy
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, ph,k,posn, skip,phase
 REAL                                   :: lnphi(np,nc),isofugacity
 CHARACTER (LEN=2)                      :: compon
! ----------------------------------------------------------------------


posn = 0
DO i = 1,n_unkw
  posn = posn + 1
  IF (it(i) == 't') t = y(posn)
  IF (it(i) == 'p') p = y(posn)
  IF (it(i) == 'lnp') p = EXP( y(posn) )
  IF (it(i) == 'fls') alpha = y(posn)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') lnx(1,k) = y(posn)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') lnx(2,k) = y(posn)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') compon = it(i)(3:3)
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') READ(compon,*) k
  IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') lnx(3,k) = y(posn)
END DO

DO k = 1,ncomp
  IF (lnx(1,k) > 0.0) lnx(1,k) = 0.0
  IF (lnx(2,k) > 0.0) lnx(2,k) = 0.0
END DO

IF (p < 1.E-100) p = 1.E-12
!IF ( IsNaN( p ) ) p = 1000.0           ! rebounce for the case of NaN-solver output
!IF ( IsNaN( t ) ) t = 300.0            ! rebounce for the case of NaN-solver output
!IF ( IsNaN( alpha ) ) alpha = 0.5      ! rebounce for the case of NaN-solver output
IF ( p /= p ) p = 1000.0                ! rebounce for the case of NaN-solver output
IF ( t /= t ) t = 300.0                 ! rebounce for the case of NaN-solver output
IF ( alpha /= alpha ) alpha = 0.5       ! rebounce for the case of NaN-solver output

! --- setting of mole fractions ----------------------------------------
DO ph = 1, nphas
  DO i = 1, ncomp
    IF ( lnx(ph,i) < -300.0 ) THEN
      xi(ph,i) = 0.0
    ELSE
      xi(ph,i) = EXP( lnx(ph,i) )
    END IF
  END DO
END DO

IF (ncomp > 1) CALL x_summation

CALL fugacity (lnphi)

phase = 2
DO i = 1,n_unkw
  skip = 0  !for ions/polymers, the isofug-eq. is not always solved
  IF (n_unkw < (ncomp*(nphas-1))) skip = ncomp*(nphas-1) - n_unkw
  IF ((i+skip-ncomp*(phase-2)) > ncomp) phase = phase + 1
  residu(i) = isofugacity((i+skip-ncomp*(phase-2)),phase,lnphi)
  if ( density_error(phase) /= 0.0 ) residu(i) = residu(i) + SIGN( density_error(phase),  residu(i) ) * 0.001
END DO

END SUBROUTINE objec_fct



!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! REAL FUNCTION isofugacity
!!
!! calculates the deviation from the condition of equal fugacities in
!! logarithmic form.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 REAL FUNCTION isofugacity (i,phase,lnphi)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: i
 INTEGER, INTENT(IN)                    :: phase
 REAL, INTENT(IN)                       :: lnphi(np,nc)
!
! ----------------------------------------------------------------------
 INTEGER                                :: p1, p2
! ----------------------------------------------------------------------


! p1=1
p1 = phase-1
p2 = phase

isofugacity = scaling(i) *( lnphi(p2,i)+lnx(p2,i)-lnx(p1,i)-lnphi(p1,i) )
! write (*,'(a,      4G18.8)') ' t, p ',t,p,dense(1),dense(2)
! write (*,'(a,i3,i3,3G18.8)') ' phi_V',i,p2,lnx(p2,i),lnphi(p2,i),dense(p2)
! write (*,'(a,i3,i3,3G18.8)') ' phi_L',i,p1,lnx(p1,i),lnphi(p1,i),dense(p1)
! write (*,*) ' ISOFUGACITY',i,ISOFUGACITY, scaling(i)
! write (*,'(a,i3,4G18.8)') ' ISOFUGACITY',i,ISOFUGACITY, lnphi(p2,i)+lnx(p2,i), -lnx(p1,i)-lnphi(p1,i)
! pause

END FUNCTION isofugacity

















!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE vle_min(lle_check)
!
 USE PARAMETERS, ONLY: RGAS
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 LOGICAL, INTENT(OUT)                   :: lle_check

 INTEGER                                :: i,j,k,phasen(0:40),steps
 REAL                                   :: lnphi(np,nc)
 REAL                                   :: vlemin(0:40),llemin(0:40),xval(0:40)
 REAL                                   :: start_xv(0:40),start_xl(0:40),x_sav,dg_dx2
! ----------------------------------------------------------------------



start_xl = 0.
start_xv = 0.


j = 0
k = 0
nphas = 2

steps = 40

x_sav = xi(1,1)
sum_rel(1) = 'x12' ! summation relation
sum_rel(2) = 'x22' ! summation relation

DO i = 0, steps
  densta(1) = 0.45
  densta(2) = 1.d-6
  xi(1,1) = 1.0 - REAL(i) / REAL(steps)
  IF ( xi(1,1) <= 1.E-50 ) xi(1,1) = 1.E-50
  xi(2,1)  = xi(1,1)
  lnx(1,1) = LOG(xi(1,1))
  lnx(2,1) = LOG(xi(2,1))
  
  CALL x_summation
  CALL fugacity (lnphi)
  CALL enthalpy_etc         !!KANN DAS RAUS????
  



  xval(i) = xi(1,1)
  llemin(i)= gibbs(1) +(xi(1,1)*lnx(1,1)+xi(1,2)*lnx(1,2))*RGAS*t

  IF ( ABS(1.0-dense(1)/dense(2)) > 0.0001 ) THEN
    vlemin(i)= gibbs(1) +(xi(1,1)*lnx(1,1)+xi(1,2)*lnx(1,2))*RGAS*t  &
        - ( gibbs(2) +(xi(2,1)*lnx(2,1)+xi(2,2)*lnx(2,2))*RGAS*t )
    phasen(i) = 2
  ELSE
    phasen(i) = 1
  END IF
  
  IF (i > 0 .AND. phasen(i) == 2) THEN
    IF (phasen(i-1) == 2 .AND. ABS(vlemin(i)+vlemin(i-1)) <  &
          ABS(vlemin(i))+ABS(vlemin(i-1))) THEN
      j = j + 1
      start_xv(j)=xval(i-1) + (xval(i)-xval(i-1))  &
          * ABS(vlemin(i-1))/ABS(vlemin(i)-vlemin(i-1))
    END IF
  END IF
  
END DO


DO i=2,steps-2
  dg_dx2 = (-llemin(i-2)+16.0*llemin(i-1)-30.0*llemin(i)  &
      +16.0*llemin(i+1)-llemin(i+2)) / (12.0*((xval(i)-xval(i-1))**2))
  IF (dg_dx2 < 0.0) THEN
    k = k + 1
    start_xl(k)=xval(i)
  END IF
END DO


IF (start_xl(1) == 0.0 .AND. start_xv(1) /= 0.0) THEN
  xi(1,1) = start_xv(1)
  xi(1,2) = 1.0-xi(1,1)
  lle_check=.false.
  ! write (*,*) 'VLE is likely', xi(1,1),xi(1,2)
ELSE IF (start_xl(1) /= 0.0 .AND. start_xv(1) == 0.0) THEN
  xi(1,1) = start_xl(1)
  xi(1,2) = 1.0-xi(1,1)
  ! write (*,*) 'LLE is likely', xi(1,1),xi(1,2)
  lle_check=.true.
ELSE IF (start_xl(1) /= 0.0 .AND. start_xv(1) /= 0.0) THEN
  xi(1,1) = start_xv(1)
  xi(1,2) = 1.0-xi(1,1)
  ! write(*,*) 'starting with VLE and check for LLE'
  lle_check=.true.
ELSE
  xi(1,1) = x_sav
  xi(1,2) = 1.0 - xi(1,1)
END IF


CALL x_summation

END SUBROUTINE vle_min


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE phase_stability
!!
!! the index 'LLE_check' is for the starting density (which determines
!! whether a liquid or vapor phase is found) of the trial phase. The
!! feed-point exits either as a vapor or as a liquid. If it can exist as
!! both (feedphases=2), then both states are tested.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 SUBROUTINE phase_stability ( lle_check, flashcase, ph_split,user )

!Petsc modules
 USE PetscManagement

!VLE modules
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 USE EOS_VARIABLES, ONLY: dhs, PI, x, eta, eta_start, z3t, fres
 USE EOS_NUMERICAL_DERIVATIVES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 type (userctx)         user 
 LOGICAL                                 :: lle_check
 LOGICAL, INTENT(IN OUT)                 :: flashcase
 INTEGER, INTENT(OUT)                    :: ph_split
! ----------------------------------------------------------------------

 INTERFACE
   REAL FUNCTION F_STABILITY ( optpara, n )
     INTEGER, INTENT(IN)                     :: n
     REAL, INTENT(IN OUT)                    :: optpara(n)
   END FUNCTION
 END INTERFACE

!INTERFACE
!  SUBROUTINE F_STABILITY (fmin, optpara, n)
!    REAL, INTENT(IN OUT)                :: fmin
!    REAL, INTENT(IN)                    :: optpara(:)
!    INTEGER, INTENT(IN)                 :: n
!  END SUBROUTINE F_STABILITY
!
!  SUBROUTINE stability_grad (g, optpara, n)
!    REAL, INTENT(IN OUT)                :: g(:)
!    REAL, INTENT(IN)                    :: optpara(:)
!    INTEGER, INTENT(IN)                 :: n
!  END SUBROUTINE stability_grad
!
!  SUBROUTINE stability_hessian (hessian, g, fmin, optpara, n)
!    REAL, INTENT(IN OUT)                :: hessian(:,:)
!    REAL, INTENT(IN OUT)                :: g(:)
!    REAL, INTENT(IN OUT)                :: fmin
!    REAL, INTENT(IN)                    :: optpara(:)
!    INTEGER, INTENT(IN)                 :: n
!  END SUBROUTINE stability_hessian
!END INTERFACE

 INTEGER                                 :: n, PRIN
 REAL                                    :: fmin, t0, h0, MACHEP, PRAXIS
 REAL, ALLOCATABLE                       :: optpara(:)

 INTEGER                                 :: i, feedphases, trial
 REAL                                    :: rhoi(nc),rho_start
 REAL                                    :: feeddens, rho_phas(np)
 REAL                                    :: fden
 REAL                                    :: dens
 REAL                                    :: rhot
 REAL                                    :: lnphi(np,nc)
 REAL                                    :: w(np,nc), mean_mass
! ----------------------------------------------------------------------

n = ncomp
ALLOCATE( optpara(n) )


IF(user%rank == 0) THEN
 IF (lle_check) WRITE (*,*) ' stability test starting with dense phase'
END IF

DO i = 1, ncomp     ! setting feed-phase x's
  IF (.NOT.flashcase) xif(i) = xi(1,i)
  IF (flashcase) xi(1,i) = xif(i)
  xi(2,i) = xif(i)  ! feed is tested for both: V and L density
END DO

densta(1) = 0.45
densta(2) = 1.d-6

CALL dens_calc(rho_phas)
IF ( ABS(1.0-dense(1)/dense(2)) > 0.0005 ) THEN
  feedphases=2      ! feed-composition can exist both, in V and L
ELSE
  feedphases=1      ! feed-composition can exist either in V or L
END IF
densta(1) = dense(1)
feeddens  = dense(2)
!write (*,*) 'feedphases',dense(1), dense(2),feedphases

10  CONTINUE        ! IF FeedPhases=2 THEN there is a second cycle

  trial = 1
  
  ! --------------------------------------------------------------------
  ! setting trial-phase mole-fractions
  ! if there is no phase-split then further trial-phases are
  ! considered (loop: 20 CONTINUE)
  ! --------------------------------------------------------------------
  DO i = 1, ncomp
    w(2,i) = 1.0 / REAL(ncomp)
  END DO
  mean_mass = 1.0 /  SUM( w(2,1:ncomp)/mm(1:ncomp) )
  xi(2,1:ncomp) = w(2,1:ncomp)/mm(1:ncomp) * mean_mass

  20  CONTINUE
  
  DO i = 1, ncomp
    rhoif(i) = rho_phas(1) * xif(i)
    rhoi(i)  = rhoif(i)
  END DO
  
  !write (*,'(a,6G16.8)') 'startval',rho_phas(2),xi(2,1:ncomp)

  ! --------------------------------------------------------------------
  ! calc Helmholtz energy density and derivative (numerical) to rhoif(i).
  ! The derivative is taken around the "feed-point" not the trial phase
  ! --------------------------------------------------------------------

  rhot = SUM( rhoi(1:ncomp) )
  x(1:ncomp) = rhoi(1:ncomp) / rhot
  CALL PERTURBATION_PARAMETER
  xi(1,1:ncomp) = x(1:ncomp)
  eta = rhot * z3t
  eta_start = eta
  densta(1) = eta_start
  ensemble_flag = 'tv'
  CALL FUGACITY (lnphi)
  ensemble_flag = 'tp'

  call fden_calc ( fden, rhoi )
  fdenf = fden

  grad_fd(1:ncomp) = lnphi(1,1:ncomp) + LOG( rhoi(1:ncomp) )


  ! --------------------------------------------------------------------
  ! starting values for iteration (optpara)
  ! --------------------------------------------------------------------
  rho_start = 1.E-5
  IF (lle_check) THEN
    densta(2) = 0.45
    CALL dens_calc(rho_phas)
    rho_start = rho_phas(2)*0.45/dense(2)
  END IF
  DO i = 1,ncomp
    rhoi(i) = xi(2,i)*rho_start
    optpara(i) = LOG( rhoi(i) )
  END DO
  
  ! --------------------------------------------------------------------
  ! minimizing the objective fct. Phase split for values of fmin < 0.0
  ! --------------------------------------------------------------------
  t0 = 5.E-5
  h0 = 0.5
  PRIN = 0
  MACHEP = 1.E-15

  fmin = PRAXIS( t0, MACHEP, h0, n, PRIN, optpara, F_STABILITY, fmin )


  ! --------------------------------------------------------------------
  ! updating the ln(x) valus from optpara. The optimal optpara-vector is
  ! not necessarily the one that was last evaluated. At the very end, 
  ! cg_decent writes the best values to optpara
  ! --------------------------------------------------------------------
  fmin = F_STABILITY( optpara, n )



  ! IF ( n == 2 ) THEN
  !   CALL Newton_Opt_2D ( stability_hessian, F_stability, optpara, n, 1.E-8, 1.E-8, g, fmin)
  ! ELSE
  !   CALL cg_descent (1.d-5, optpara, n, F_STABILITY, stability_grad, STATUS, &
  !                    gnorm, fmin,  iter, nfunc, ngrad, d, g, xtemp, gtemp)
  ! ENDIF
  ! CALL F_STABILITY (fmin, optpara, n)

  
  ! --------------------------------------------------------------------
  ! determine instability & non-trivial solution
  ! --------------------------------------------------------------------
  ph_split = 0
  IF (fmin < -1.E-7 .AND.  &
      ABS( 1.0 - maxval(EXP(optpara),mask=optpara /= 0.0) /maxval(rhoif) ) > 0.0005) THEN
    ph_split = 1
  END IF

  IF (ph_split == 1) THEN

    ! ------------------------------------------------------------------
    ! here, there should be IF FeedPhases=2 THEN GOTO 10
    ! and test for another phase (while saving optpara)
    ! ------------------------------------------------------------------

    rhoi2(1:ncomp) = EXP( optpara(1:ncomp) )
    dens = PI/6.0 * SUM( rhoi2(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
    rhot = SUM( rhoi2(1:ncomp) )
    xi(2,1:ncomp) = rhoi2(1:ncomp) / rhot

  ELSE

    IF (trial <= ncomp + ncomp) THEN
      ! ----------------------------------------------------------------
      ! setting trial-phase x's
      ! ----------------------------------------------------------------
      IF (trial <= ncomp) THEN
         DO i=1,ncomp
            w(2,i) = 1.0 / REAL(ncomp-1) * 0.05
         END DO
         w(2,trial) = 0.95
      ELSE
         DO i=1,ncomp
            w(2,i) = 1.0 / REAL(ncomp-1) * 0.00001
         END DO
         w(2,trial-ncomp) = 0.99999
      END IF
      mean_mass = 1.0 / SUM( w(2,1:ncomp)/mm(1:ncomp) )
      xi(2,1:ncomp) = w(2,1:ncomp)/mm(1:ncomp) * mean_mass
      trial = trial + 1
      GO TO 20
    END IF
    ! IF (.NOT.LLE_check) write (*,*) 'no phase split detected'
    ! IF (.NOT.LLE_check) pause
    IF (feedphases > 1 .AND. .NOT.lle_check .AND. densta(1) > 0.2) THEN
      densta(1) = feeddens         ! this will be the lower-valued density (vapor)
      CALL dens_calc(rho_phas)
      ! WRITE (*,*) 'try feed as vapor-phase'
      GO TO 10
    END IF

  END IF

DEALLOCATE( optpara )

END SUBROUTINE phase_stability


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE select_sum_rel
!!
!! This subroutine determines which component of a phase "ph" is calculated
!! from the summation relation x_i = 1 - sum(x_j). The other components are,
!! by default, said to be iterated during the phase equilibrium calculation.
!!
!! Note that for flash calculations not all of these mole fractions are in
!! fact iterated - this is raken care of in "determine_flash_it".
!!
!! ph           phase
!! excl         exclude comp. n
!! startindex   assign it(startindex) for quantities to be iterated
!!              (further it(startindex+1) is assigned, for a ternary
!!              mixture, etc.)
!!
!! sum_index    indicates the component, with the largest mole
!!              fraction. If ph=1 and sum_index=2, we define
!!              sum_rel(ph=1)='x12', so that this component is
!!              calculated from the summation relation.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 SUBROUTINE select_sum_rel (ph,excl,startindex)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: ph
 INTEGER, INTENT(IN)                    :: excl
 INTEGER, INTENT(IN)                    :: startindex
! ----------------------------------------------------------------------
 INTEGER                                :: i,j, sum_index
 REAL                                   :: xmax(np)
 ! CHARACTER                              :: compNo*2,phasNo*2
! ----------------------------------------------------------------------

xmax(ph) = 0.0
DO i = 1, ncomp

  IF ( xi(ph,i) > xmax(ph) ) THEN
    xmax(ph) = xi(ph,i)
    sum_index = i
    
    IF (ph == 1 .AND. i == 1) sum_rel(1) = 'x11'
    IF (ph == 1 .AND. i == 2) sum_rel(1) = 'x12'
    IF (ph == 1 .AND. i == 3) sum_rel(1) = 'x13'
    IF (ph == 1 .AND. i == 4) sum_rel(1) = 'x14'
    IF (ph == 1 .AND. i == 5) sum_rel(1) = 'x15'
    
    IF (ph == 2 .AND. i == 1) sum_rel(2) = 'x21'
    IF (ph == 2 .AND. i == 2) sum_rel(2) = 'x22'
    IF (ph == 2 .AND. i == 3) sum_rel(2) = 'x23'
    IF (ph == 2 .AND. i == 4) sum_rel(2) = 'x24'
    IF (ph == 2 .AND. i == 5) sum_rel(2) = 'x25'
    
    IF (ph == 3 .AND. i == 1) sum_rel(3) = 'x31'
    IF (ph == 3 .AND. i == 2) sum_rel(3) = 'x32'
    IF (ph == 3 .AND. i == 3) sum_rel(3) = 'x33'
    IF (ph == 3 .AND. i == 4) sum_rel(3) = 'x34'
    IF (ph == 3 .AND. i == 5) sum_rel(3) = 'x35'
!    write (*,*) ph,i,xi(ph,i),sum_rel(ph)
  END IF

END DO

j = 0
DO i = 1, ncomp

  IF ( i /= sum_index .AND. i /= excl ) THEN
    IF (ph == 1 .AND. i == 1) it(startindex+j) = 'x11'
    IF (ph == 1 .AND. i == 2) it(startindex+j) = 'x12'
    IF (ph == 1 .AND. i == 3) it(startindex+j) = 'x13'
    IF (ph == 1 .AND. i == 4) it(startindex+j) = 'x14'
    IF (ph == 1 .AND. i == 5) it(startindex+j) = 'x15'
    
    IF (ph == 2 .AND. i == 1) it(startindex+j) = 'x21'
    IF (ph == 2 .AND. i == 2) it(startindex+j) = 'x22'
    IF (ph == 2 .AND. i == 3) it(startindex+j) = 'x23'
    IF (ph == 2 .AND. i == 4) it(startindex+j) = 'x24'
    IF (ph == 2 .AND. i == 5) it(startindex+j) = 'x25'
    
    IF (ph == 3 .AND. i == 1) it(startindex+j) = 'x31'
    IF (ph == 3 .AND. i == 2) it(startindex+j) = 'x32'
    IF (ph == 3 .AND. i == 3) it(startindex+j) = 'x33'
    IF (ph == 3 .AND. i == 4) it(startindex+j) = 'x34'
    IF (ph == 3 .AND. i == 5) it(startindex+j) = 'x35'
!    write (*,*) 'iter  ',it(startindex+j)
    j = j + 1
  END IF

END DO

END SUBROUTINE select_sum_rel




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE tangent_plane
!
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
!!$ INTERFACE
!!$   SUBROUTINE tangent_value (fmin, optpara, n)
!!$     INTEGER, INTENT(IN)                 :: n
!!$     REAL, INTENT(IN OUT)                :: fmin
!!$     REAL, INTENT(IN)                    :: optpara(:)
!!$   END SUBROUTINE tangent_value
!!$
!!$   SUBROUTINE tangent_grad (g, optpara, n)
!!$     INTEGER, INTENT(IN)                 :: n
!!$     REAL, INTENT(IN OUT)                :: g(:)
!!$     REAL, INTENT(IN)                    :: optpara(:)
!!$   END SUBROUTINE tangent_grad
!!$ END INTERFACE

!
! ----------------------------------------------------------------------
 INTERFACE
   REAL FUNCTION PRAXIS( t0, MACHEP, h0, n, PRIN, optpara, TANGENT_VALUE, fmin )
     REAL, INTENT(IN OUT)               :: t0
     REAL, INTENT(IN)                   :: machep
     REAL, INTENT(IN)                   :: h0
     INTEGER                            :: n
     INTEGER, INTENT(IN OUT)            :: prin
     REAL, INTENT(IN OUT)               :: optpara(n)
     REAL, EXTERNAL                     :: TANGENT_VALUE
     REAL, INTENT(IN OUT)               :: fmin
   END FUNCTION

   REAL FUNCTION TANGENT_VALUE2 ( optpara, n )
     INTEGER, INTENT(IN)                :: n
     REAL, INTENT(IN)                   :: optpara(n)
   END FUNCTION
 END INTERFACE
!
! ----------------------------------------------------------------------
 INTEGER                                :: n
 INTEGER                                :: i, k, ph
 INTEGER                                :: small_i, min_ph, other_ph
 INTEGER                                :: PRIN
 REAL                                   :: fmin , t0, h0, MACHEP
 REAL                                   :: lnphi(np,nc)
 REAL, ALLOCATABLE                      :: optpara(:)

! INTEGER                                :: STATUS,  iter, nfunc, ngrad
! REAL                                   :: gnorm
! REAL, ALLOCATABLE                      :: d(:), g(:), xtemp(:), gtemp(:)
! ----------------------------------------------------------------------

n = ncomp
t0 = 1.E-4
h0 = 0.1
PRIN = 0
MACHEP = 1.E-15

ALLOCATE( optpara(n) )
!ALLOCATE( d(n) )
!ALLOCATE( g(n) )
!ALLOCATE( xtemp(n) )
!ALLOCATE( gtemp(n) )

DO i = 1,ncomp
  rhoi1(i) = rhoif(i)
  lnx(1,i) = LOG(xi(1,i))
  lnx(2,i) = LOG(xi(2,i))
END DO

DO i = 1,ncomp
  optpara(i) = LOG( xi(2,i) * 0.001 )
END DO

! CALL cg_descent (1.d-4, optpara, n, tangent_value, tangent_grad, STATUS, &
!                  gnorm, fmin,  iter, nfunc, ngrad, d, g, xtemp, gtemp)
!
! updating the ln(x) valus from optpara. The optimal optpara-vector is not necessarily
! the one that was last evaluated. At the very end, cg_decent writes the best values to optpara
! CALL tangent_value (fmin, optpara, n)



fmin = PRAXIS( t0, MACHEP, h0, n, PRIN, optpara, TANGENT_VALUE2, fmin )

! The optimal optpara-vector is not necessarily the one that was last evaluated.
! TANGENT_VALUE is reexecuted with the optimal vector optpara, in order to update the ln(x) values 
fmin = TANGENT_VALUE2( optpara, n )


! ----------------------------------------------------------------------
! If one component is a polymer (indicated by a low component-density)
! then get an estimate of the polymer-lean composition, by solving for
! xi_p1 = ( xi_p2 * phii_p2) / phii_p1     (phase equilibrium condition,
! with p1 for phase 1)
! ----------------------------------------------------------------------
IF ( MINVAL( lnx(1,1:ncomp) ) < MINVAL( lnx(2,1:ncomp) ) ) THEN
  min_ph   = 1
  other_ph = 2
ELSE
  min_ph   = 2
  other_ph = 1
ENDIF
small_i = MINLOC( lnx(min_ph,1:ncomp), 1 )
! --- if one component is a polymer ------------------------------------
IF ( MINVAL( lnx(min_ph,1:ncomp) ) < -20.0 ) THEN
  CALL FUGACITY ( lnphi )
  lnx(min_ph,small_i) = lnx(other_ph,small_i)+lnphi(other_ph,small_i) - lnphi(min_ph,small_i)
  optpara(small_i) = lnx(2,small_i) + LOG( SUM( EXP( optpara(1:ncomp) ) ) )
END IF

! ----------------------------------------------------------------------
! caution: these initial values are for a flashcase overwritten in
! SUBROUTINE determine_flash_it2, because in that case, the lnx-values
! treated as ln(mole_number).
! ----------------------------------------------------------------------
val_init(1) = dense(1)
val_init(2) = dense(2)
val_init(3) = t
val_init(4) = p
DO ph = 1,nphas
  DO k = 1,ncomp
    val_init(4+k+(ph-1)*ncomp) = lnx(ph,k)
  END DO
END DO
!alpha = optpara(1)


!DEALLOCATE( optpara, d, g, xtemp, gtemp )
DEALLOCATE( optpara )

END SUBROUTINE tangent_plane


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE determine_flash_it2
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, k, ph
 REAL                                   :: n_phase1, n_phase2, max_x_diff
! ----------------------------------------------------------------------

 IF ( MINVAL( lnx(1,1:ncomp) ) < MINVAL( lnx(2,1:ncomp) ) ) THEN
    it(1) = 'x11'
    it(2) = 'x12'
    IF (ncomp >= 3) it(3) = 'x13'
    IF (ncomp >= 4) it(4) = 'x14'
    IF (ncomp >= 5) it(5) = 'x15'
    sum_rel(1) = 'nfl'
 ELSE
    it(1) = 'x21'
    it(2) = 'x22'
    IF (ncomp >= 3) it(3) = 'x23'
    IF (ncomp >= 4) it(4) = 'x24'
    IF (ncomp >= 5) it(5) = 'x25'
    sum_rel(2) = 'nfl'
 ENDIF
 max_x_diff = 0.0
 DO i = 1,ncomp
    IF ( ABS( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) ) > max_x_diff ) THEN
       max_x_diff = ABS( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) )
       n_phase1 = ( xif(i) - EXP( lnx(2,i) ) ) / ( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) )
       n_phase2 = 1.0 - n_phase1
    END IF
 END DO
 lnx(1,1:ncomp) = lnx(1,1:ncomp) + LOG( n_phase1 )   ! these x's are treated as mole numbers
 lnx(2,1:ncomp) = lnx(2,1:ncomp) + LOG( n_phase2 )   ! these x's are treated as mole numbers


 val_init(1) = dense(1)
 val_init(2) = dense(2)
 val_init(3) = t
 val_init(4) = p
 DO ph = 1,nphas
    DO k = 1,ncomp
       val_init(4+k+(ph-1)*ncomp) = lnx(ph,k) ! - LOG( SUM( EXP( lnx(ph,1:ncomp) ) ) )
       !            write (*,*) ph,k, lnx(ph,k)
    END DO
 END DO

END SUBROUTINE determine_flash_it2

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE poly_sta_var
!
! This subroutine generates starting values for mole fractons of
! polymer-solvent systems.
! The determination of these starting values follows a two-step
! procedure. Fist, the equilibrium concentration of the polymer-rich
! phase is estimated with the assumption of zero concentration
! of polymer in the polymer-lean-phase. This is achieved in the
! SUBROUTINE POLYMER_FREE. (Only one equation has to be iterated
! for this case). Once this is achieved, the rigorous calculation
! is triggered. If it converges, fine! If no solution is obtained,
! the pressure is somewhat reduced, the procedure is repeated and
! a calculation is started to approach the originally specified
! pressure.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE poly_sta_var (converg)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN OUT)                 :: converg
!
! ----------------------------------------------------------------------
 INTEGER                                 :: k,ph,sol
 REAL                                    :: p_spec,solution(10,4+nc*np)
! ----------------------------------------------------------------------

 p_spec = p

 find_equilibrium: DO

    CALL polymer_free(p_spec,sol,solution)

    WRITE (*,*) ' '
    WRITE (*,*) ' GENERATING STARTING VALUES'

    val_init(1) = solution(1,1)   ! approx.solutions for next iteration
    val_init(2) = solution(1,2)   ! approx.solutions for next iteration
    val_init(3) = solution(1,3)   ! approx.solutions for next iteration
    val_init(4) = solution(1,4)   ! approx.solutions for next iteration
    DO ph = 1,nphas
       DO k = 1,ncomp
          val_init(4+k+(ph-1)*ncomp) = solution(1,4+k+(ph-1)*ncomp)
       END DO
    END DO
    val_init(7) = -10000.0       ! start.val. for lnx(2,1) for iterat.

    IF (p /= p_spec)  &
         WRITE (*,*) ' INITIAL EQUILIBRIUM CALC. FAILD. NEXT STEP STARTS'

    IF (p == p_spec) THEN
       n_unkw = ncomp       ! number of quantities to be iterated
       it(1)='x11'          ! iteration of mol fraction of comp.1 phase 1
       it(2)='x21'          ! iteration of mol fraction of comp.1 phase 2
       CALL objective_ctrl (converg)
    ELSE
       outp = 0            ! output to terminal
       running ='p'        ! Pressure is running var. in PHASE_EQUILIB
       CALL phase_equilib(p_spec,5.0,converg)
    END IF

    IF (converg == 1) EXIT find_equilibrium
    p = p * 0.9
    IF ( p < (0.7*p_spec) ) WRITE (*,*) 'Surface Tension Code, Phase equilibrium calculation: NO SOLUTION FOUND'
    IF ( p < (0.7*p_spec) ) STOP 5

 END DO find_equilibrium

 WRITE (*,*) ' FINISHED: POLY_STA_VAR'

END SUBROUTINE poly_sta_var


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE x_summation
!!
!! This subroutine solves the summation relation: xi=1-sum(xj)
!! The variable "sum_rel(i)" contains the information, which mole
!! fraction is the one to be calculated here. Consider the example
!! sum_rel(1)='x12'. The fist letter 'x' of this variable indicates,
!! that this subroutine needs to be executed and that the mole
!! fraction of a component has to be calculated. The second letter
!! of the string points to phase 1, the third letter to component 2.
!!     If the fist letter is 'e', not 'x', then the subroutine
!! NEUTR_CHARGE is called. This is the case of electrolyte solutions,
!! neutral charges have to be enforced in all phases (see below).
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!!
 SUBROUTINE x_summation
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, comp_i, ph_i
 REAL                                   :: sum_x
 CHARACTER (LEN=2)                      :: phasno
 CHARACTER (LEN=2)                      :: compno
 LOGICAL                                :: flashcase2
! ----------------------------------------------------------------------

DO j = 1, nphas
  IF (sum_rel(j)(1:3) == 'nfl') THEN
     CALL new_flash (j)
     RETURN
  END IF
END DO



flashcase2 = .false.

DO j = 1, nphas

  IF (sum_rel(j)(1:1) == 'x') THEN

    phasno = sum_rel(j)(2:2)
    READ(phasno,*) ph_i
    compno = sum_rel(j)(3:3)
    READ(compno,*) comp_i
    IF ( sum_rel(nphas+j)(1:1) == 'e' ) CALL neutr_charge(nphas+j)

    sum_x = 0.0
    DO i = 1, ncomp
      IF ( i /= comp_i ) sum_x = sum_x + xi(ph_i,i)
    END DO
    xi(ph_i,comp_i) = 1.0 - sum_x
    IF ( xi(ph_i,comp_i ) <  0.0 ) xi(ph_i,comp_i) = 0.0
    IF ( xi(ph_i,comp_i ) /= 0.0 ) THEN
      lnx(ph_i,comp_i) = LOG( xi(ph_i,comp_i) )
    ELSE
      lnx(ph_i,comp_i) = -100000.0
    END IF
    ! write (*,*) 'sum_x',ph_i,comp_i,lnx(ph_i,comp_i),xi(ph_i,comp_i)

  ELSE IF ( sum_rel(j)(1:2) == 'fl' ) THEN

    flashcase2 = .true.
    ! ------------------------------------------------------------------
    ! This case is true when all molefractions of one phase are
    ! determined from a component balance. What is needed to
    ! calculate all molefractions of that phase are all mole-
    ! fractions of the other phase (nphas=2, so far) and the
    ! phase fraction alpha.
    ! Alpha is calculated (in FLASH_ALPHA) from the mole fraction
    ! of component {sum_rel(j)(3:3)}. IF sum_rel(2)='fl3', then
    ! the alpha is determined from the molefraction of comp. 3 and
    ! the molefraction of phase 2 is then completely determined        ELSE
    ! ------------------------------------------------------------------

  ELSE
    WRITE (*,*) 'Surface Tension Code, Phase equilibrium calculation: summation relation not defined'
    STOP 5
  END IF

END DO

IF ( it(1) == 'fls' ) CALL flash_sum
IF ( flashcase2 ) CALL flash_alpha

END SUBROUTINE x_summation



!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE FUGACITY
!!
!! This subroutine serves as an interface to the eos-subroutines.
!! (1) case 1, when ensemble_flag = 'tp'
!!      mu_i^res(T,p,x)/kT = ln( phi_i )      
!!     and in addition, the densities that satisfy the specified p
!! (2) case 2, when ensemble_flag = 'tv'
!!     The subroutine gives the residual chemical potential:
!!     -->   mu_i^res(T,rho,x)/kT
!!     and in addition the resulting pressure for the given density.
!! The term "residual" means difference of the property and the same
!! property for an ideal gas mixture.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! 
 SUBROUTINE FUGACITY (ln_phi)
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES, ONLY: phas, x, eta, eta_start, lnphi, fres, rho, pges, KBOL
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: ln_phi(np,nc)
!
! --- local variables --------------------------------------------------
 INTEGER                                :: ph
! ----------------------------------------------------------------------
!
IF (eos < 2) THEN
  
  DO ph = 1,nphas
 
    phas = ph
    eta_start = densta(ph)
    x(1:ncomp)   = xi(ph,1:ncomp)

    IF (p < 1.E-100) THEN
      !WRITE(*,*) ' FUGACITY: PRESSURE TOO LOW', p
      p = 1.E-6
    END IF

    IF (num == 0) CALL PHI_EOS
    IF (num == 1) CALL PHI_NUMERICAL
    !!IF (num == 2) CALL PHI_CRITICAL_RENORM
    IF (num == 2) write(*,*) 'CRITICAL RENORM NOT INCLUDED YET'

    dense(ph) = eta
    ln_phi(ph,1:ncomp) = lnphi(1:ncomp)
    ! gibbs(ph) = fres + sum( xi(ph,1:ncomp)*( log( xi(ph,1:ncomp)*rho) - 1.0 ) )  &
    !            + (pges * 1.d-30) / (KBOL*t*rho)     ! includes ideal gas contribution

    ! f_res(ph) = fres
    ! write (*,'(i3,4G20.11)') ph,eta,lnphi(1),lnphi(2)
    ! DO i = 1,ncomp
    !   DO j=1,NINT(parame(i,12))
    !     mxx(ph,i,j) = mx(i,j)
    !   END DO
    ! END DO

  END DO
  
ELSE

!  IF (eos == 2) CALL srk_eos (ln_phi)
!  IF (eos == 3) CALL  pr_eos (ln_phi)
!  dense(1) = 0.01
!  dense(2) = 0.3
!  IF (eos == 4.OR.eos == 5.OR.eos == 6.OR.eos == 8) CALL lj_fugacity(ln_phi)
!  IF (eos == 7) CALL sw_fugacity(ln_phi)
!  IF (eos == 9) CALL lj_bh_fug(ln_phi)

END IF

END SUBROUTINE FUGACITY







!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE FUGACITY
!!
!! This subroutine serves as an interface to the eos-subroutines.
!! (1) case 1, when ensemble_flag = 'tp'
!!     The subroutine gives the residual chemical potential:
!!      mu_i^res(T,p,x)/kT = ln( phi_i )      
!!     and in addition, the densities that satisfy the specified p
!! (2) case 2, when ensemble_flag = 'tv'
!!     The subroutine gives the residual chemical potential:
!!     -->   mu_i^res(T,rho,x)/kT
!!     and in addition the resulting pressure for the given density.
!! The term "residual" means difference of the property and the same
!! property for an ideal gas mixture.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! 
 SUBROUTINE FUGACITY2 (ln_phi)
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES, ONLY: phas, x, eta, eta_start, lnphi, fres, rho, pges, KBOL
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: ln_phi(np,nc)
!
! --- local variables --------------------------------------------------
 INTEGER                                :: ph
! ----------------------------------------------------------------------
!
IF (eos < 2) THEN
  
 ! DO ph = 1,nphas
 
    phas = 2! ph
    eta_start = densta(2)
    x(1:ncomp)   = xi(2,1:ncomp)

    IF (p < 1.E-100) THEN
      WRITE(*,*) ' FUGACITY: PRESSURE TOO LOW', p
      p = 1.E-6
    END IF

    IF (num == 0) CALL PHI_EOS2
    IF (num == 1) CALL PHI_NUMERICAL
    !!IF (num == 2) CALL PHI_CRITICAL_RENORM
    IF (num == 2) write(*,*) 'CRITICAL RENORM NOT INCLUDED YET'

    dense(2) = eta
    ln_phi(2,1:ncomp) = lnphi(1:ncomp)
    ! gibbs(ph) = fres + sum( xi(ph,1:ncomp)*( log( xi(ph,1:ncomp)*rho) - 1.0 ) )  &
    !            + (pges * 1.d-30) / (KBOL*t*rho)     ! includes ideal gas contribution

    ! f_res(ph) = fres
    ! write (*,'(i3,4G20.11)') ph,eta,lnphi(1),lnphi(2)
    ! DO i = 1,ncomp
    !   DO j=1,NINT(parame(i,12))
    !     mxx(ph,i,j) = mx(i,j)
    !   END DO
    ! END DO

 ! END DO
  
ELSE

!  IF (eos == 2) CALL srk_eos (ln_phi)
!  IF (eos == 3) CALL  pr_eos (ln_phi)
!  dense(1) = 0.01
!  dense(2) = 0.3
!  IF (eos == 4.OR.eos == 5.OR.eos == 6.OR.eos == 8) CALL lj_fugacity(ln_phi)
!  IF (eos == 7) CALL sw_fugacity(ln_phi)
!  IF (eos == 9) CALL lj_bh_fug(ln_phi)

END IF

END SUBROUTINE FUGACITY2










!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE enthalpy_etc
!
! This subroutine serves as an interface to the EOS-routines. The
! residual enthalpy h_res, residual entropy s_res, residual Gibbs
! enthalpy g_res, and residual heat capacity at constant pressure
! (cp_res) corresponding to converged conditions are calculated.
! The conditions in (T,P,xi,rho) need to be converged equilibrium
! conditions  !!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE enthalpy_etc
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES
 IMPLICIT NONE
!
 INTEGER                                :: ph
! ------------------------------------------------------------------

IF (eos <= 1) THEN

  DO ph=1,nphas
  
    phas = ph
    eta       = dense(ph)
!    eta_start = dense(ph)
    x(1:ncomp)   = xi(ph,1:ncomp)
  
    IF(num == 0) THEN
      CALL H_EOS
    ELSE
      IF(num == 1) CALL H_numerical
      IF(num == 2) write (*,*) 'Surface Tension Code, Phase equilibrium calculation: enthalpy_etc: incorporate H_EOS_RN'
      IF(num == 2) stop 5
!      IF(num == 2) CALL H_EOS_rn
    END IF
    enthal(ph) = h_res
    entrop(ph) = s_res
    ! gibbs(ph)  = h_res - t * s_res  ! already defined in eos.f90 (including ideal gas)
    cpres(ph)  = cp_res
  
  END DO
  IF (nphas == 2) h_lv = enthal(2)-enthal(1)

ENDIF

END SUBROUTINE enthalpy_etc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE dens_calc
!
! This subroutine serves as an interface to the EOS-routines. The
! densities corresponding to given (P,T,xi) are calculated.
! (Note: the more common interface is SUBROUTINE FUGACITY.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE dens_calc(rho_phas)
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES
 IMPLICIT NONE
!
!
!------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: rho_phas(np)
!
 INTEGER                                :: ph
!------------------------------------------------------------------


DO ph = 1, nphas
  
  IF (eos < 2) THEN
    
    phas = ph
    eta       = densta(ph)
    eta_start = densta(ph)
    x(1:ncomp)   = xi(ph,1:ncomp)
    
    CALL PERTURBATION_PARAMETER
    CALL DENSITY_ITERATION

    dense(ph)= eta
    rho_phas(ph) = eta/z3t
    
  ELSE
    write (*,*) 'Surface Tension Code, Phase equilibrium calculation: SUBROUTINE DENS_CALC not available for cubic EOS'
    stop 5
  END IF
  
END DO

END SUBROUTINE dens_calc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE fden_calc (fden, rhoi)
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: fden
 REAL, INTENT(IN OUT)                   :: rhoi(nc)
! ----------------------------------------------------------------------
 REAL                                   :: rhot, fden_id
! ----------------------------------------------------------------------


IF (eos < 2) THEN

  rhot = SUM( rhoi(1:ncomp) )
  x(1:ncomp) = rhoi(1:ncomp) / rhot
  
  CALL PERTURBATION_PARAMETER
  eta = rhot * z3t
  eta_start = eta

  IF (num == 0) THEN
    CALL F_EOS
  ELSE IF(num == 1) THEN
    CALL F_NUMERICAL
  ELSE
    write (*,*) 'Surface Tension Code, Phase equilibrium calculation: deactivated this line when making a transition to f90'
    stop 5
    ! CALL F_EOS_rn
  END IF

  fden_id = SUM( rhoi(1:ncomp) * ( LOG( rhoi(1:ncomp) ) - 1.0 ) )

  fden = fres * rhot  +  fden_id
  
ELSE
  write (*,*) 'Surface Tension Code, Phase equilibrium calculation: SUBROUTINE FDEN_CALC not available for cubic EOS'
  stop 5
END IF

END SUBROUTINE fden_calc




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE polymer_free
!
! This subroutine performes a phase equilibrium calculation assuming
! the polymer-lean hase to be polymer-free (x_poly=0). Only the
! equality of the solvent-fugacities has to be ensured (only one
! equation to be iterated). This procedure delivers very good
! appoximations for the polymer-rich phase up-to fairly close to the
! mixture critical point. Both, liquid-liquid and vapor-liquid
! equilibria can be calculated.
! See also comments to SUBROUTINE POLY_STA_VAR.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE polymer_free (p_spec,sol,solution)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                    :: p_spec
 INTEGER, INTENT(OUT)                    :: sol
 REAL, INTENT(OUT)                       :: solution(10,4+nc*np)
!
! ----------------------------------------------------------------------
 INTEGER                                 :: k,j,ph, converg
 REAL                                    :: grid(10)
! ----------------------------------------------------------------------

 sol = 0

 grid(1)=0.98
 grid(2)=0.9
 grid(3)=0.7
 grid(4)=0.5
 grid(5)=0.3
 grid(6)=0.2
 grid(7)=0.1
 grid(8)=0.05

 DO WHILE ( sol == 0 )

    DO j = 1,8
       ! Phase 2 is solvent-phase
       ! starting value for xi(1,1) of polymer-phase 1: w_polymer=0.95 to 0.05
       ! from simple approximate equation
       xi(1,1) = grid(j) / ( (1.0-grid(j)) * mm(1) / mm(2) )            !xi(1,1) Phase 1 Komponente 1
       IF ( mm(1) < 5000.0 ) xi(1,1) = xi(1,1) * 0.8
       xi(1,2) = 1.0 - xi(1,1)                                          !xi(1,2) Phase 1 Komponente 2
       lnx(1,1) = LOG(xi(1,1))
       lnx(1,2) = LOG(xi(1,2))
       lnx(2,1) = -1.E10                                                !ln(xi) Phase 2 Komponente 1
       lnx(2,2) = 0.0                                                   !ln(xi) Phase 2 Komponente 2



       val_init(1) = 0.45               ! starting density targeting at a liquid phase
       val_init(2) = 0.0001             ! starting density targeting at a vapor phase
       ! val_init(2) = 0.45
       val_init(3) = t
       val_init(4) = p
       DO ph = 1,nphas
          DO k = 1,ncomp
             val_init(4+k+(ph-1)*ncomp) = lnx(ph,k)
          END DO
       END DO




       n_unkw = ncomp-1                 ! number of quantities to be iterated
       it(1) = 'x11'                    ! iteration of mol fraction of comp.1 phase 1
       it(2) = ' '
       sum_rel(1) = 'x12'               ! summation relation: x12 = 1 - sum(x1j)
       sum_rel(2) = 'x22'

       CALL objective_ctrl (converg)

       IF (converg == 1 .AND. ABS(dense(1)/dense(2)-1.0) > 1.d-3 .AND. dense(1) > 0.1) THEN
          IF (sol == 0) THEN
             sol = sol + 1
             DO k = 1,4+ncomp*nphas
                solution(sol,k) = val_conv(k)
             END DO
          ELSE IF (ABS(solution(sol,5)/lnx(1,1)-1.0) > 1.d-2) THEN
             sol = sol + 1
             DO k = 1,4+ncomp*nphas
                solution(sol,k) = val_conv(k)
             END DO
          END IF
       END IF

    END DO





    IF (sol == 0) THEN
       WRITE (*,*) ' no initial solution found'
       p = p*0.9
       IF (p < (0.7*p_spec)) WRITE (*,*) 'Surface Tension Code, Phase equilibrium calculation: NO SOLUTION FOUND'
       IF (p < (0.7*p_spec)) STOP 5
    ELSE IF (sol > 1) THEN
       ! write (*,*) ' '
       ! write (*,*) ' ',sol,' solutions found:'
       ! write (*,*) ' lnx(1,1),   dichte_1,   dichte_2'
       ! DO k = 1,sol
       !    write (*,*) solution(k,5),solution(k,1),solution(k,2)
       ! END DO
    END IF
 END DO

 n_unkw = ncomp               ! number of quantities to be iterated
 it(1) = 'x11'                ! iteration of mol fraction of comp.1 phase 1
 it(2) = 'x21'                ! iteration of mol fraction of comp.1 phase 2
 sum_rel(1) = 'x12'           ! summation relation: x12 = 1 - sum(x1j)
 sum_rel(2) = 'x22'           ! summation relation: x22 = 1 - sum(x2j)


 END SUBROUTINE polymer_free





!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE phase_equilib
!!
!! This subroutine varies a predefined "running variable" and
!! organizes phase equilibrium calculations. For an isotherm
!! calculation e.g., the running variable is often the pressure. The
!! code is designed to deliver only converged solutions. In order to
!! enforce convergence, a step-width adjustment (reduction) is
!! implemented.
!!
!!   VARIABLE LIST:
!! running     defines the running variable. For example: if you want
!!             to calculate the vapor pressure curve of a component
!!             starting from 100C to 200C, then running is 't'. The
!!             temperature is step-wise increased until the end-
!!             -temperature of 200C is reached.
!!             (in this example end_x=200+273.15)
!! end_x       end point for running variable
!! steps       No. of calculation steps towards the end point of calc.
!! converg     0 if no convergence achieved, 1 if converged solution
!!
!!   PREREQUISITES:
!! prior to execution of this routine, the follwing variables have to
!! be defined: "val_init" an array containing the starting values for
!! this iteration, "it(i)" provides the information, which variable is
!! determined iteratively, "sum_rel(i)" indicates, which mole fraction
!! is determined from the summation relation sum(xi)=1. Furthermore,
!! the number of phases and the variables provided by the subroutine
!! INPUT are required.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! 
 SUBROUTINE phase_equilib (end_x,steps,converg)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN)                       :: end_x
 REAL, INTENT(IN)                       :: steps
 INTEGER, INTENT(OUT)                   :: converg
!
! ----------------------------------------------------------------------
 INTEGER                                :: k, count1,count2,runindex,maxiter
 REAL                                   :: delta_x,delta_org,val_org,runvar
 CHARACTER (LEN=2)                      :: compon
 LOGICAL                                :: continue_cycle
! ----------------------------------------------------------------------

IF (running(1:2) == 'd1') runindex = 1
IF (running(1:2) == 'd2') runindex = 2
IF (running(1:1) == 't')  runindex = 3
IF (running(1:1) == 'p')  runindex = 4
IF (running(1:2) == 'x1') compon = running(3:3)
IF (running(1:2) == 'x1') READ(compon,*) k
IF (running(1:2) == 'x1') runindex = 4+k
IF (running(1:2) == 'x2') compon = running(3:3)
IF (running(1:2) == 'x2') READ(compon,*) k
IF (running(1:2) == 'x2') runindex = 4+ncomp+k
IF (running(1:2) == 'l1') compon = running(3:3)
IF (running(1:2) == 'l1') READ(compon,*) k
IF (running(1:2) == 'l1') runindex = 4+k
IF (running(1:2) == 'l2') compon = running(3:3)
IF (running(1:2) == 'l2') READ(compon,*) k
IF (running(1:2) == 'l2') runindex = 4+ncomp+k

maxiter = 200
IF ( ncomp >= 3 ) maxiter = 1000
count1 = 0
count2 = 0
delta_x   = ( end_x - val_init(runindex) ) / steps        !J: calc increment in running var = (phi_end - phi_init)/steps
delta_org = ( end_x - val_init(runindex) ) / steps
val_org = val_init(runindex)
IF ( running(1:1) == 'x' ) THEN
  delta_x   = ( end_x - EXP(val_init(runindex)) ) / steps
  delta_org = ( end_x - EXP(val_init(runindex)) ) / steps
  val_org = EXP(val_init(runindex))
END IF

continue_cycle = .true.

DO WHILE ( continue_cycle )

   count1 = count1 + 1
   count2 = count2 + 1
   ! val_org = val_init(runindex)


   CALL objective_ctrl (converg)

   IF (converg == 1) THEN
      val_init( 1:(4+ncomp*nphas) ) = val_conv( 1:(4+ncomp*nphas) )
      IF (outp == 1 .AND. (ABS(delta_x) > 0.1*ABS(delta_org) .OR. count2 == 2)) CALL output
   ELSE
      delta_x = delta_x / 2.0
      IF (num == 2) delta_x = delta_x / 2.0
      val_init(runindex) = val_org
      IF (running(1:1) == 'x') val_init(runindex) = LOG(val_org)
      continue_cycle = .true.
      count2 = 0
   END IF
   runvar = val_init(runindex)
   IF (running(1:1) == 'x') runvar = EXP(val_init(runindex))

   IF ( end_x == 0.0 .AND. running(1:1) /= 'x' ) THEN
      IF ( ABS(runvar-end_x) < 1.E-8 ) continue_cycle = .false.
   ELSE IF ( ABS((runvar-end_x)/end_x) < 1.E-8 ) THEN
      ! IF(delta_org.NE.0.0) WRITE (*,*)' FINISHED ITERATION',count1
      continue_cycle = .false.
   ELSE IF ( count1 == maxiter ) THEN
      WRITE (*,*) ' MAX. NO OF ITERATIONS',count1
      converg = 0
      continue_cycle = .false.
   ELSE IF ( ABS(delta_x) < 1.E-5*ABS(delta_org) ) THEN
      ! WRITE (*,*) ' CLOSEST APPROACH REACHED',count1
      converg = 0
      continue_cycle = .false.
   ELSE
      continue_cycle = .true.
      val_org = runvar
      IF (ABS(runvar+delta_x-end_x) > ABS(runvar-end_x)) delta_x = end_x - runvar    ! if end-point passed
      val_init(runindex) = runvar + delta_x
      IF (running(1:1) == 'x') val_init(runindex) = LOG(runvar+delta_x)
   END IF

   IF (ABS(delta_x) < ABS(delta_org) .AND. count2 >= 5) THEN
      delta_x = delta_x * 2.0
      count2 = 0
   END IF

END DO            ! continue_cycle

END SUBROUTINE phase_equilib

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE new_flash (ph_it)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: ph_it
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, ph_cal
 REAL, DIMENSION(nc)                    :: ni_1, ni_2
! ----------------------------------------------------------------------

  ph_cal = 3 - ph_it          ! for two phases only

  DO i = 1, ncomp
    IF ( lnx(ph_it,i) < -300.0 ) THEN
      ni_2(i) = 0.0
    ELSE
      ni_2(i) = EXP( lnx(ph_it,i) )
    END IF
  END DO

  DO i = 1, ncomp
    ni_1(i) = xif(i)-ni_2(i)
    IF ( ni_2(i) > xif(i) ) THEN
       ni_2(i) = xif(i)
       ni_1(i) = xif(i) * 1.E-20
    ENDIF
  END DO

  xi(ph_it,1:ncomp) = ni_2(1:ncomp) / SUM( ni_2(1:ncomp) )
  DO i = 1, ncomp
    IF ( xi(ph_it,i) >= 1.E-300 ) lnx(ph_it,i) = LOG( xi(ph_it,i) )
  END DO
  xi(ph_cal,1:ncomp) = ni_1(1:ncomp) / SUM( ni_1(1:ncomp) )
  lnx(ph_cal,1:ncomp) = LOG( xi(ph_cal,1:ncomp) )

END SUBROUTINE new_flash


!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! SUBROUTINE PHI_EOS
!!
!! This subroutine gives the residual chemical potential:
!! -->   mu_i^res(T,p,x)/kT = ln( phi_i )       when ensemble_flag = 'tp'
!! The required input for this case (T, p, x(nc)) and as a starting value
!! eta_start
!!
!! or it gives
!!
!! -->   mu_i^res(T,rho,x)/kT                   when ensemble_flag = 'tv'
!! The required input for this case (T, eta_start, x(nc)). Note that
!! eta_start is the specified density (packing fraction) in this case.
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! 
 SUBROUTINE PHI_EOS
!
 USE PARAMETERS
 USE EOS_VARIABLES
 USE EOS_CONSTANTS


 IMPLICIT NONE
!
! --- local variables---------------------------------------------------
 INTEGER                                :: i, j, k, ki, l, m
 REAL                                   :: z0, z1, z2, z3, z0_rk, z1_rk, z2_rk, z3_rk
 REAL                                   :: zms, m_mean
 REAL, DIMENSION(nc)                    :: mhs, mdsp, mhc, myres
 REAL, DIMENSION(nc)                    :: m_rk
 REAL                                   :: gij_rk(nc,nc)
 REAL                                   :: zres, zges
 REAL                                   :: dpdz, dpdz2

 REAL                                   :: I1, I2, I1_rk, I2_rk
 REAL                                   :: ord1_rk, ord2_rk
 REAL                                   :: c1_con, c2_con, c1_rk
 REAL                                   :: zmr, nmr, zmr_rk, nmr_rk, um_rk
 REAL, DIMENSION(nc,0:6)                :: ap_rk, bp_rk

 LOGICAL                                :: assoc
 REAL                                   :: ass_s2, m_hbon(nc)

 REAL                                   :: fdd_rk, fqq_rk, fdq_rk
 REAL, DIMENSION(nc)                    :: my_dd, my_qq, my_dq
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! obtain parameters and density independent expressions
! ----------------------------------------------------------------------
CALL PERTURBATION_PARAMETER


! ----------------------------------------------------------------------
! density iteration: (pTx)-ensemble   OR   p calc.: (pvx)-ensemble
! ----------------------------------------------------------------------
IF (ensemble_flag == 'tp') THEN
   CALL DENSITY_ITERATION
ELSEIF (ensemble_flag == 'tv') THEN
  eta = eta_start
  CALL P_EOS
ELSE
  write (*,*) 'Surface Tension Code, Phase equilibrium calculation: PHI_EOS: define ensemble, ensemble_flag == (pv) or (pt)'
  stop 5
END IF


! --- Eq.(A.8) ---------------------------------------------------------
rho = eta / z3t
z0  = z0t * rho
z1  = z1t * rho
z2  = z2t * rho
z3  = z3t * rho

m_mean = z0t / (PI/6.0)
zms    = 1.0 - eta

! ----------------------------------------------------------------------
! compressibility factor z = p/(kT*rho)
! ----------------------------------------------------------------------
zges = (p * 1.d-30) / (KBOL*t*rho)
IF ( ensemble_flag == 'tv' ) zges = (pges * 1.d-30) / (KBOL*t*rho)
zres = zges - 1.0



! ======================================================================
! calculate the derivatives of f to mole fraction x ( d(f)/d(x) )
! ======================================================================

DO  k = 1, ncomp
  
  z0_rk = PI/6.0 * mseg(k)
  z1_rk = PI/6.0 * mseg(k) * dhs(k)
  z2_rk = PI/6.0 * mseg(k) * dhs(k)*dhs(k)
  z3_rk = PI/6.0 * mseg(k) * dhs(k)**3
  
! --- derivative d(m_mean)/d(x) ----------------------------------------
  m_rk(k) = ( mseg(k) - m_mean ) / rho
  ! lij(1,2)= -0.050
  ! lij(2,1)=lij(1,2)
  ! r_m2dx(k)=0.0
  ! m_mean2=0.0
  ! DO i =1,ncomp
  !    r_m2dx(k)=r_m2dx(k)+2.0*x(i)*(mseg(i)+mseg(k))/2.0*(1.0-lij(i,k))
  !    DO j =1,ncomp
  !       m_mean2=m_mean2+x(i)*x(j)*(mseg(i)+mseg(j))/2.0*(1.0-lij(i,j))
  !    ENDDO
  ! ENDDO

  ! --------------------------------------------------------------------
  ! d(f)/d(x) : hard sphere contribution
  ! --------------------------------------------------------------------
  mhs(k) =  6.0/PI* (  3.0*(z1_rk*z2+z1*z2_rk)/zms + 3.0*z1*z2*z3_rk/zms/zms  &
                + 3.0*z2*z2*z2_rk/z3/zms/zms + z2**3 *z3_rk*(3.0*z3-1.0)/z3/z3/zms**3   &
                + ((3.0*z2*z2*z2_rk*z3-2.0*z2**3 *z3_rk)/z3**3 -z0_rk) *LOG(zms)  &
                + (z0-z2**3 /z3/z3)*z3_rk/zms  )

  ! --------------------------------------------------------------------
  ! d(f)/d(x) : chain term
  ! --------------------------------------------------------------------
  DO i = 1, ncomp
    DO j = 1, ncomp
      gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
      gij_rk(i,j) = z3_rk/zms/zms  &
                      + 3.0*dij_ab(i,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
                      + dij_ab(i,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
    END DO
  END DO
  
  mhc(k) = 0.0
  DO i = 1, ncomp
    mhc(k) = mhc(k) + x(i)*rho * (1.0-mseg(i)) / gij(i,i) * gij_rk(i,i)
  END DO
  mhc(k) = mhc(k) + ( 1.0-mseg(k)) * LOG( gij(k,k) )

  
  ! --------------------------------------------------------------------
  ! PC-SAFT:  d(f)/d(x) : dispersion contribution
  ! --------------------------------------------------------------------
  IF (eos == 1) THEN
    
    ! --- derivatives of apar, bpar to rho_k ---------------------------
    DO m = 0, 6
      ap_rk(k,m) = m_rk(k)/m_mean**2 * ( ap(m,2) + (3.0 -4.0/m_mean) *ap(m,3) )
      bp_rk(k,m) = m_rk(k)/m_mean**2 * ( bp(m,2) + (3.0 -4.0/m_mean) *bp(m,3) )
    END DO

    I1    = 0.0
    I2    = 0.0
    I1_rk = 0.0
    I2_rk = 0.0
    DO m = 0, 6
      I1  = I1 + apar(m)*eta**REAL(m)
      I2  = I2 + bpar(m)*eta**REAL(m)
      I1_rk = I1_rk + apar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + ap_rk(k,m)*eta**REAL(m)
      I2_rk = I2_rk + bpar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + bp_rk(k,m)*eta**REAL(m)
    END DO
    
    ord1_rk  = 0.0
    ord2_rk  = 0.0
    DO i = 1,ncomp
      ord1_rk = ord1_rk + 2.0*mseg(k)*rho*x(i)*mseg(i)*sig_ij(i,k)**3  *uij(i,k)/t
      ord2_rk = ord2_rk + 2.0*mseg(k)*rho*x(i)*mseg(i)*sig_ij(i,k)**3 *(uij(i,k)/t)**2 
    END DO
    
    c1_con= 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3*z3)/zms**4   &
                    + (1.0 - m_mean)*(20.0*z3-27.0*z3*z3 +12.0*z3**3 -2.0*z3**4 )  &
                    /(zms*(2.0-z3))**2  )
    c2_con= - c1_con*c1_con *(  m_mean*(-4.0*z3*z3+20.0*z3+8.0)/zms**5   &
                                + (1.0 - m_mean) *(2.0*z3**3 +12.0*z3*z3-48.0*z3+40.0)  &
                                  /(zms*(2.0-z3))**3  )
    c1_rk= c2_con*z3_rk - c1_con*c1_con*m_rk(k)   *  ( (8.0*z3-2.0*z3*z3)/zms**4   &
           - (-2.0*z3**4 +12.0*z3**3 -27.0*z3*z3+20.0*z3) / (zms*(2.0-z3))**2  )
    
    mdsp(k) = -2.0*PI* ( order1*rho*rho*I1_rk + ord1_rk*I1 )  &
              -    PI* c1_con*m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
              -    PI* ( c1_con*m_rk(k) + c1_rk*m_mean ) * order2*rho*rho*I2
    
  ! --------------------------------------------------------------------
  ! SAFT:  d(f)/d(x) : dispersion contribution
  ! --------------------------------------------------------------------
  ELSE
    
    zmr    = 0.0
    nmr    = 0.0
    zmr_rk = 0.0
    nmr_rk = 0.0
    DO i = 1, ncomp
      DO j = 1, ncomp
        zmr = zmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)*uij(i,j)
        nmr = nmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)
      END DO
      zmr_rk = zmr_rk + 2.0*mseg(k) * x(i)*mseg(i)*vij(k,i)*uij(k,i)
      nmr_rk = nmr_rk + 2.0*mseg(k) * x(i)*mseg(i)*vij(k,i)
    END DO
    
    um_rk = 1.0/nmr**2 * ( nmr*zmr_rk - zmr*nmr_rk )
    
    mdsp(k) = 0.0
    DO i = 1,4
      DO j = 1,9
        mdsp(k) = mdsp(k) + dnm(i,j)*(um/t)**REAL(i)*(eta/tau)**REAL(j)  &
                 * ( 1.0 + z3_rk*rho/eta*REAL(j) + um_rk*rho/um*REAL(i) )
      END DO
    END DO

  END IF
  ! --- end of dispersion contribution------------------------------------
  
  
  ! --------------------------------------------------------------------
  ! TPT-1-association according to Chapman et al.
  ! --------------------------------------------------------------------
  m_hbon(k) = 0.0
  assoc = .false.
  DO i = 1,ncomp
    IF (nhb_typ(i) /= 0) assoc = .true.
  END DO
  IF (assoc) THEN

    ass_s2  = 0.0
    DO l = 1, nhb_typ(k)
      ass_s2  = ass_s2  + nhb_no(k,l) * LOG(mx(k,l))
    END DO

    m_hbon(k)=ass_s2
    DO i = 1, ncomp
      DO ki = 1, nhb_typ(i)
        DO j = 1, ncomp
          DO l = 1, nhb_typ(j)
            m_hbon(k)= m_hbon(k) - rho*rho/2.0*x(i)*x(j) *mx(i,ki)*mx(j,l) *nhb_no(i,ki)*nhb_no(j,l)  &
                                   * gij_rk(i,j) * ass_d(i,j,ki,l)
          END DO
        END DO
      END DO
    END DO

  END IF
  ! --- end of TPT-1-association accord. to Chapman --------------------
  
  
  ! --------------------------------------------------------------------
  ! polar terms
  ! --------------------------------------------------------------------
  CALL PHI_POLAR ( k, z3_rk, fdd_rk, fqq_rk, fdq_rk )
  my_dd(k) = fdd_rk
  my_qq(k) = fqq_rk
  my_dq(k) = fdq_rk

  
  ! --------------------------------------------------------------------
  ! d(f)/d(x) : summation of all contributions
  ! --------------------------------------------------------------------
  myres(k) = mhs(k) +mhc(k) +mdsp(k) +m_hbon(k) +my_dd(k) +my_qq(k) +my_dq(k)




END DO


! ----------------------------------------------------------------------
! finally calculate
! mu_i^res(T,p,x)/kT = ln( phi_i )       when ensemble_flag = 'tp'
! mu_i^res(T,rho,x)/kT                   when ensemble_flag = 'tv'
! ----------------------------------------------------------------------

DO k = 1, ncomp
  ! write (*,*) k,myres(k) +LOG(rho*x(k)),rho
  IF (ensemble_flag == 'tp' ) lnphi(k) = myres(k) - LOG(zges)
  IF (ensemble_flag == 'tv' ) lnphi(k) = myres(k)
  ! write (*,*) 'in',k,EXP(lnphi(k)),LOG(zges),eta
END DO
!write (*,'(5G18.10)') lnphi(1),rho

dpdz  = pgesdz
dpdz2 = pgesd2

END SUBROUTINE PHI_EOS








!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PHI_EOS
!
! This subroutine gives the residual chemical potential:
! -->   mu_i^res(T,p,x)/kT = ln( phi_i )       when ensemble_flag = 'tp'
! The required input for this case (T, p, x(nc)) and as a starting value
! eta_start
!
! or it gives
!
! -->   mu_i^res(T,rho,x)/kT                   when ensemble_flag = 'tv'
! The required input for this case (T, eta_start, x(nc)). Note that
! eta_start is the specified density (packing fraction) in this case.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE PHI_EOS2
!
 USE PARAMETERS
 USE EOS_VARIABLES
 USE EOS_CONSTANTS


 IMPLICIT NONE
!
! --- local variables---------------------------------------------------
 INTEGER                                :: i, j, k, ki, l, m
 REAL                                   :: z0, z1, z2, z3, z0_rk, z1_rk, z2_rk, z3_rk
 REAL                                   :: zms, m_mean
 REAL, DIMENSION(nc)                    :: mhs, mdsp, mhc, myres
 REAL, DIMENSION(nc)                    :: m_rk
 REAL                                   :: gij_rk(nc,nc)
 REAL                                   :: zres, zges
 REAL                                   :: dpdz, dpdz2

 REAL                                   :: I1, I2, I1_rk, I2_rk
 REAL                                   :: ord1_rk, ord2_rk
 REAL                                   :: c1_con, c2_con, c1_rk
 REAL                                   :: zmr, nmr, zmr_rk, nmr_rk, um_rk
 REAL, DIMENSION(nc,0:6)                :: ap_rk, bp_rk

 LOGICAL                                :: assoc
 REAL                                   :: ass_s2, m_hbon(nc)

 REAL                                   :: fdd_rk, fqq_rk, fdq_rk
 REAL, DIMENSION(nc)                    :: my_dd, my_qq, my_dq
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! obtain parameters and density independent expressions
! ----------------------------------------------------------------------
CALL PERTURBATION_PARAMETER


! ----------------------------------------------------------------------
! density iteration: (pTx)-ensemble   OR   p calc.: (pvx)-ensemble
! ----------------------------------------------------------------------
IF (ensemble_flag == 'tp') THEN
   CALL DENSITY_ITERATION
ELSEIF (ensemble_flag == 'tv') THEN
  eta = eta_start
  CALL P_EOS
ELSE
  write (*,*) 'Surface Tension Code, Phase equilibrium calculation: PHI_EOS: define ensemble, ensemble_flag == (pv) or (pt)'
  stop 5
END IF


! --- Eq.(A.8) ---------------------------------------------------------
rho = eta / z3t
z0  = z0t * rho
z1  = z1t * rho
z2  = z2t * rho
z3  = z3t * rho

m_mean = z0t / (PI/6.0)
zms    = 1.0 - eta

! ----------------------------------------------------------------------
! compressibility factor z = p/(kT*rho)
! ----------------------------------------------------------------------
zges = (p * 1.d-30) / (KBOL*t*rho)
IF ( ensemble_flag == 'tv' ) zges = (pges * 1.d-30) / (KBOL*t*rho)
zres = zges - 1.0



! ======================================================================
! calculate the derivatives of f to mole fraction x ( d(f)/d(x) )
! ======================================================================

DO  k = 1, ncomp
  
  z0_rk = PI/6.0 * mseg(k)
  z1_rk = PI/6.0 * mseg(k) * dhs(k)
  z2_rk = PI/6.0 * mseg(k) * dhs(k)*dhs(k)
  z3_rk = PI/6.0 * mseg(k) * dhs(k)**3
  
! --- derivative d(m_mean)/d(x) ----------------------------------------
  m_rk(k) = ( mseg(k) - m_mean ) / rho
  ! lij(1,2)= -0.050
  ! lij(2,1)=lij(1,2)
  ! r_m2dx(k)=0.0
  ! m_mean2=0.0
  ! DO i =1,ncomp
  !    r_m2dx(k)=r_m2dx(k)+2.0*x(i)*(mseg(i)+mseg(k))/2.0*(1.0-lij(i,k))
  !    DO j =1,ncomp
  !       m_mean2=m_mean2+x(i)*x(j)*(mseg(i)+mseg(j))/2.0*(1.0-lij(i,j))
  !    ENDDO
  ! ENDDO

  ! --------------------------------------------------------------------
  ! d(f)/d(x) : hard sphere contribution
  ! --------------------------------------------------------------------
  mhs(k) =  6.0/PI* (  3.0*(z1_rk*z2+z1*z2_rk)/zms + 3.0*z1*z2*z3_rk/zms/zms  &
                + 3.0*z2*z2*z2_rk/z3/zms/zms + z2**3 *z3_rk*(3.0*z3-1.0)/z3/z3/zms**3   &
                + ((3.0*z2*z2*z2_rk*z3-2.0*z2**3 *z3_rk)/z3**3 -z0_rk) *LOG(zms)  &
                + (z0-z2**3 /z3/z3)*z3_rk/zms  )

  ! --------------------------------------------------------------------
  ! d(f)/d(x) : chain term
  ! --------------------------------------------------------------------
  DO i = 1, ncomp
    DO j = 1, ncomp
      gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
      gij_rk(i,j) = z3_rk/zms/zms  &
                      + 3.0*dij_ab(i,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
                      + dij_ab(i,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
    END DO
  END DO
  
  mhc(k) = 0.0
  DO i = 1, ncomp
    mhc(k) = mhc(k) + x(i)*rho * (1.0-mseg(i)) / gij(i,i) * gij_rk(i,i)
  END DO
  mhc(k) = mhc(k) + ( 1.0-mseg(k)) * LOG( gij(k,k) )

  
  ! --------------------------------------------------------------------
  ! PC-SAFT:  d(f)/d(x) : dispersion contribution
  ! --------------------------------------------------------------------
  IF (eos == 1) THEN
    
    ! --- derivatives of apar, bpar to rho_k ---------------------------
    DO m = 0, 6
      ap_rk(k,m) = m_rk(k)/m_mean**2 * ( ap(m,2) + (3.0 -4.0/m_mean) *ap(m,3) )
      bp_rk(k,m) = m_rk(k)/m_mean**2 * ( bp(m,2) + (3.0 -4.0/m_mean) *bp(m,3) )
    END DO

    I1    = 0.0
    I2    = 0.0
    I1_rk = 0.0
    I2_rk = 0.0
    DO m = 0, 6
      I1  = I1 + apar(m)*eta**REAL(m)
      I2  = I2 + bpar(m)*eta**REAL(m)
      I1_rk = I1_rk + apar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + ap_rk(k,m)*eta**REAL(m)
      I2_rk = I2_rk + bpar(m)*REAL(m)*eta**REAL(m-1)*z3_rk + bp_rk(k,m)*eta**REAL(m)
    END DO
    
    ord1_rk  = 0.0
    ord2_rk  = 0.0
    DO i = 1,ncomp
      ord1_rk = ord1_rk + 2.0*mseg(k)*rho*x(i)*mseg(i)*sig_ij(i,k)**3  *uij(i,k)/t
      ord2_rk = ord2_rk + 2.0*mseg(k)*rho*x(i)*mseg(i)*sig_ij(i,k)**3 *(uij(i,k)/t)**2 
    END DO
    
    c1_con= 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3*z3)/zms**4   &
                    + (1.0 - m_mean)*(20.0*z3-27.0*z3*z3 +12.0*z3**3 -2.0*z3**4 )  &
                    /(zms*(2.0-z3))**2  )
    c2_con= - c1_con*c1_con *(  m_mean*(-4.0*z3*z3+20.0*z3+8.0)/zms**5   &
                                + (1.0 - m_mean) *(2.0*z3**3 +12.0*z3*z3-48.0*z3+40.0)  &
                                  /(zms*(2.0-z3))**3  )
    c1_rk= c2_con*z3_rk - c1_con*c1_con*m_rk(k)   *  ( (8.0*z3-2.0*z3*z3)/zms**4   &
           - (-2.0*z3**4 +12.0*z3**3 -27.0*z3*z3+20.0*z3) / (zms*(2.0-z3))**2  )
    
    mdsp(k) = -2.0*PI* ( order1*rho*rho*I1_rk + ord1_rk*I1 )  &
              -    PI* c1_con*m_mean * ( order2*rho*rho*I2_rk + ord2_rk*I2 )  &
              -    PI* ( c1_con*m_rk(k) + c1_rk*m_mean ) * order2*rho*rho*I2
    
  ! --------------------------------------------------------------------
  ! SAFT:  d(f)/d(x) : dispersion contribution
  ! --------------------------------------------------------------------
  ELSE
    
    zmr    = 0.0
    nmr    = 0.0
    zmr_rk = 0.0
    nmr_rk = 0.0
    DO i = 1, ncomp
      DO j = 1, ncomp
        zmr = zmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)*uij(i,j)
        nmr = nmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)
      END DO
      zmr_rk = zmr_rk + 2.0*mseg(k) * x(i)*mseg(i)*vij(k,i)*uij(k,i)
      nmr_rk = nmr_rk + 2.0*mseg(k) * x(i)*mseg(i)*vij(k,i)
    END DO
    
    um_rk = 1.0/nmr**2 * ( nmr*zmr_rk - zmr*nmr_rk )
    
    mdsp(k) = 0.0
    DO i = 1,4
      DO j = 1,9
        mdsp(k) = mdsp(k) + dnm(i,j)*(um/t)**REAL(i)*(eta/tau)**REAL(j)  &
                 * ( 1.0 + z3_rk*rho/eta*REAL(j) + um_rk*rho/um*REAL(i) )
      END DO
    END DO

  END IF
  ! --- end of dispersion contribution------------------------------------
  
  
  ! --------------------------------------------------------------------
  ! TPT-1-association according to Chapman et al.
  ! --------------------------------------------------------------------
  m_hbon(k) = 0.0
  assoc = .false.
  DO i = 1,ncomp
    IF (nhb_typ(i) /= 0) assoc = .true.
  END DO
  IF (assoc) THEN

    ass_s2  = 0.0
    DO l = 1, nhb_typ(k)
      ass_s2  = ass_s2  + nhb_no(k,l) * LOG(mx(k,l))
    END DO

    m_hbon(k)=ass_s2
    DO i = 1, ncomp
      DO ki = 1, nhb_typ(i)
        DO j = 1, ncomp
          DO l = 1, nhb_typ(j)
            m_hbon(k)= m_hbon(k) - rho*rho/2.0*x(i)*x(j) *mx(i,ki)*mx(j,l) *nhb_no(i,ki)*nhb_no(j,l)  &
                                   * gij_rk(i,j) * ass_d(i,j,ki,l)
          END DO
        END DO
      END DO
    END DO

  END IF
  ! --- end of TPT-1-association accord. to Chapman --------------------
  
  
  ! --------------------------------------------------------------------
  ! polar terms
  ! --------------------------------------------------------------------
  CALL PHI_POLAR ( k, z3_rk, fdd_rk, fqq_rk, fdq_rk )
  my_dd(k) = fdd_rk
  my_qq(k) = fqq_rk
  my_dq(k) = fdq_rk

  
  ! --------------------------------------------------------------------
  ! d(f)/d(x) : summation of all contributions
  ! --------------------------------------------------------------------
  myres(k) = mhs(k) +mhc(k) +mdsp(k) +m_hbon(k) +my_dd(k) +my_qq(k) +my_dq(k)


END DO




muhs(1:ncomp) = mhs(1:ncomp)
muhc(1:ncomp) = mhc(1:ncomp)
mudisp(1:ncomp) = mdsp(1:ncomp)




! ----------------------------------------------------------------------
! finally calculate
! mu_i^res(T,p,x)/kT = ln( phi_i )       when ensemble_flag = 'tp'
! mu_i^res(T,rho,x)/kT                   when ensemble_flag = 'tv'
! ----------------------------------------------------------------------

DO k = 1, ncomp
  ! write (*,*) k,myres(k) +LOG(rho*x(k)),rho
  IF (ensemble_flag == 'tp' ) lnphi(k) = myres(k) - LOG(zges)
  IF (ensemble_flag == 'tv' ) lnphi(k) = myres(k)
  ! write (*,*) 'in',k,EXP(lnphi(k)),LOG(zges),eta
END DO
!write (*,'(5G18.10)') lnphi(1),rho

dpdz  = pgesdz
dpdz2 = pgesd2


write(*,*)'tp?',myres(1), lnphi(1)

END SUBROUTINE PHI_EOS2

































!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE PHI_NUMERICAL
!
 USE EOS_VARIABLES
 USE EOS_CONSTANTS
 USE DFT_MODULE, ONLY: z_ges, fres_temp
 IMPLICIT NONE
!
!-----local variables-------------------------------------------------
 INTEGER                                :: k
 REAL                                   :: zres, zges
 REAL                                   :: fres1, fres2, fres3, fres4, fres5
 REAL                                   :: delta_rho
 REAL, DIMENSION(nc)                    :: myres
 REAL, DIMENSION(nc)                    :: rhoi, rhoi_0
 REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4, tfr_5
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! density iteration or pressure calculation
!-----------------------------------------------------------------------

IF (ensemble_flag == 'tp') THEN
  CALL DENSITY_ITERATION
ELSEIF (ensemble_flag == 'tv') THEN
  eta = eta_start
  CALL P_NUMERICAL
ELSE
  write (*,*) 'Surface Tension Code, Phase equilibrium calculation: PHI_EOS: define ensemble, ensemble_flag == (tv) or (tp)'
  stop 5
END IF

!-----------------------------------------------------------------------
! compressibility factor z = p/(kT*rho)
!-----------------------------------------------------------------------

zges = (p * 1.E-30) / (kbol*t*eta/z3t)
IF ( ensemble_flag == 'tv' ) zges = (pges * 1.E-30) / (kbol*t*eta/z3t)
zres = zges - 1.0
z_ges = zges

rhoi_0(1:ncomp) = x(1:ncomp) * eta/z3t
rhoi(1:ncomp) = rhoi_0(1:ncomp)


!-----------------------------------------------------------------------
! derivative to rho_k (keeping other rho_i's constant
!-----------------------------------------------------------------------

DO  k = 1, ncomp

   IF ( rhoi_0(k) > 1.d-9 ) THEN
      delta_rho = 1.E-13 * 10.0**(0.5*(15.0+LOG10(rhoi_0(k))))
   ELSE
      delta_rho = 1.E-10
   END IF

   rhoi(k) = rhoi_0(k) + delta_rho
   eta = PI/6.0 * SUM( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
   x(1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
   rho = SUM( rhoi(1:ncomp) )
   CALL F_NUMERICAL
   fres1 = fres*rho
   tfr_1 = tfr*rho

   rhoi(k) = rhoi_0(k) + 0.5 * delta_rho
   eta = PI/6.0 * SUM( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
   x(1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
   rho = SUM( rhoi(1:ncomp) )
   CALL F_NUMERICAL
   fres2 = fres*rho
   tfr_2 = tfr*rho

   IF ( rhoi_0(k) > 1.E-9 ) THEN
      rhoi(k) = rhoi_0(k) - 0.5 * delta_rho
      eta = PI/6.0 * SUM( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
      x(1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
      rho = SUM( rhoi(1:ncomp) )
      CALL F_NUMERICAL
      fres4 = fres*rho
      tfr_4 = tfr*rho

      rhoi(k) = rhoi_0(k) - delta_rho
      eta = PI/6.0 * SUM( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
      x(1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
      rho = SUM( rhoi(1:ncomp) )
      CALL F_NUMERICAL
      fres5 = fres*rho
      tfr_5 = tfr*rho
   END IF

   rhoi(k) = rhoi_0(k)
   eta = PI/6.0 * SUM( rhoi(1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
   x(1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )
   rho = SUM( rhoi(1:ncomp) )
   CALL F_NUMERICAL
   fres3 = fres*rho
   tfr_3 = tfr*rho

   IF ( rhoi_0(k) > 1.E-9 ) THEN
      myres(k) = ( fres5 - 8.0*fres4 + 8.0*fres2 - fres1 ) / ( 6.0*delta_rho )
   ELSE
      myres(k) = ( -3.0*fres3 + 4.0*fres2 - fres1 ) / delta_rho
   END IF

END DO


!-----------------------------------------------------------------------
! residual Helmholtz energy
!-----------------------------------------------------------------------

fres_temp = fres

!-----------------------------------------------------------------------
! residual chemical potential
!-----------------------------------------------------------------------

DO  k = 1, ncomp
  IF (ensemble_flag == 'tp') lnphi(k) = myres(k) - LOG(zges)
  IF (ensemble_flag == 'tv' .AND. eta >= 0.0) lnphi(k) = myres(k) !+LOG(rho)
  ! write (*,*) 'in',k,EXP(lnphi(k)),LOG(zges),eta
  ! IF (DFT.GE.98) write (*,*) dft
  ! write (*,*) 'lnphi',k,LNPHI(k),x(k),MYRES(k), -LOG(ZGES)
  ! pause
  ! write (*,*) k, myres(k), fres, ZRES
END DO

END SUBROUTINE PHI_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!      SUBROUTINE H_EOS (phas,h_res,s_res,cp_res,X,T,P,PARAME,
!     1                  KIJ,lij,NCOMP,ETA_START,ETA,eos,tfr)
!      IMPLICIT NONE
!      INTEGER nc
!      PARAMETER (nc=20)
!      INTEGER phas,ncomp,eos,i
!      REAL kij(nc,nc),lij(nc,nc),x(nc),t,p,parame(nc,25)
!      REAL eta_start,eta,tfr,h_res,cp_res,s_res


!      i=1

!      IF (i.EQ.1) THEN
!      CALL H_EOS_1(phas,h_res,s_res,cp_res,X,T,P,PARAME,
!     1                  KIJ,lij,NCOMP,ETA_START,ETA,eos,tfr)
!      ELSE
!      CALL H_EOS_NUM (phas,h_res,s_res,cp_res,X,T,P,PARAME,
!     1                  KIJ,lij,NCOMP,ETA_START,ETA,eos,tfr)
!      ENDIF

!      RETURN
!      END


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE H_EOS
!
 USE PARAMETERS, ONLY: RGAS
 USE EOS_CONSTANTS
 USE EOS_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL                                   :: zges, df_dt, dfdr, ddfdrdr
 REAL                                   :: cv_res, df_dt2, df_drdt
 REAL                                   :: fact, dist, t_tmp, rho_0
 REAL                                   :: fdr1, fdr2, fdr3, fdr4

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
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Initializing
! ----------------------------------------------------------------------
CALL PERTURBATION_PARAMETER

rho = eta/z3t
z0 = z0t*rho
z1 = z1t*rho
z2 = z2t*rho
z3 = z3t*rho

sumseg = z0t / (PI/6.0)
zms    = 1.0 - z3


! ----------------------------------------------------------------------
! first and second derivative of f to density (dfdr,ddfdrdr)
! ----------------------------------------------------------------------
CALL P_EOS

zges = (pges * 1.E-30)/(kbol*t*rho)

dfdr    = pges/(eta*rho*(kbol*t)/1.E-30)
ddfdrdr = pgesdz/(eta*rho*(kbol*t)/1.E-30) - dfdr*2.0/eta - 1.0/eta**2 


! ----------------------------------------------------------------------
! Helmholtz Energy f/kT = fres
! ----------------------------------------------------------------------
CALL F_EOS


! ----------------------------------------------------------------------
! derivative of some auxilliary properties to temperature
! ----------------------------------------------------------------------
DO i = 1,ncomp
  dhsdt(i)=parame(i,2) *(-3.0*parame(i,3)/t/t)*0.12*EXP(-3.0*parame(i,3)/t)
  dhsdt2(i) = dhsdt(i)*3.0*parame(i,3)/t/t  &
              + 6.0*parame(i,2)*parame(i,3)/t**3  *0.12*EXP(-3.0*parame(i,3)/t)
END DO

z1tdt = 0.0
z2tdt = 0.0
z3tdt = 0.0
DO i = 1,ncomp
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
DO i = 1,ncomp
  z1tdt2 = z1tdt2 + x(i)*mseg(i)*dhsdt2(i)
  z2tdt2 = z2tdt2 + x(i)*mseg(i)*2.0 *( dhsdt(i)*dhsdt(i) +dhs(i)*dhsdt2(i) )
  z3tdt2 = z3tdt2 + x(i)*mseg(i)*3.0 *( 2.0*dhs(i)*dhsdt(i)*  &
           dhsdt(i) +dhs(i)*dhs(i)*dhsdt2(i) )
END DO
z1dt2  = PI / 6.0*z1tdt2 *rho
z2dt2  = PI / 6.0*z2tdt2 *rho
z3dt2  = PI / 6.0*z3tdt2 *rho


! ----------------------------------------------------------------------
! 1st & 2nd derivative of f/kT hard spheres to temp. (fhsdt)
! ----------------------------------------------------------------------
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


! ----------------------------------------------------------------------
! 1st & 2nd derivative of f/kT of chain term to T (fchdt)
! ----------------------------------------------------------------------
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


! ----------------------------------------------------------------------
! 1st & 2nd derivative of f/kT dispersion term to T (fdspdt)
! ----------------------------------------------------------------------
I1 = 0.0
I2 = 0.0
I1dt = 0.0
I2dt = 0.0
I1dt2= 0.0
I2dt2= 0.0
DO m = 0, 6
  I1   = I1   + apar(m)*z3**REAL(m)
  I2   = I2   + bpar(m)*z3**REAL(m)
  I1dt = I1dt + apar(m)*z3dt*REAL(m)*z3**REAL(m-1)
  I2dt = I2dt + bpar(m)*z3dt*REAL(m)*z3**REAL(m-1)
  I1dt2= I1dt2+ apar(m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
              + apar(m)*z3dt*z3dt *REAL(m)*REAL(m-1)*z3**REAL(m-2)
  I2dt2= I2dt2+ bpar(m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
              + bpar(m)*z3dt*z3dt *REAL(m)*REAL(m-1)*z3**REAL(m-2)
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


! ----------------------------------------------------------------------
! 1st & 2nd derivative of f/kT association term to T (fhbdt)
! ----------------------------------------------------------------------
fhbdt  = 0.0
fhbdt2 = 0.0
assoc = .false.
DO i = 1,ncomp
  IF ( nhb_typ(i) /= 0 ) assoc = .true.
END DO
IF (assoc) THEN
  
  DO i = 1,ncomp
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
  
  DO i = 1,ncomp
    DO j = 1,ncomp
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
          ass_d_dt(i,j,k,l)= - kap_hb(i,j) *sig_ij(i,j)**3 * eps_hb(i,j,k,l)/t/t*EXP(eps_hb(i,j,k,l)/t)
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
  
  
! ------ constants for iteration ---------------------------------------
  attenu = 0.7
  tol = 1.E-10
  IF (eta < 0.2) tol = 1.E-11
  IF (eta < 0.01) tol = 1.E-12
  max_eval = 200
  
! ------ initialize mxdt(i,j) ------------------------------------------
  DO i = 1, ncomp
    DO k = 1, nhb_typ(i)
      mxdt(i,k) = 0.0
      mxdt2(i,k) = 0.0
    END DO
  END DO
  
  
! ------ iterate over all components and all sites ---------------------
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
  WRITE(6,*) 'Surface Tension Code, Phase equilibrium calculation: CAL_PCSAFT: max_eval violated err_sum = ',err_sum,tol
  STOP 5
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


! ----------------------------------------------------------------------
! derivatives of f/kT of dipole-dipole term to temp. (fdddt)
! ----------------------------------------------------------------------
fdddt  = 0.0
fdddt2 = 0.0
dipole = 0
DO i = 1,ncomp
  my2dd(i) = 0.0
  IF ( parame(i,6) /= 0.0 .AND. uij(i,i) /= 0.0 ) THEN
    dipole = 1
    my2dd(i) = (parame(i,6))**2 *1.E-49 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**3 *1.E-30)
  END IF
  my0(i) = my2dd(i)        ! needed for dipole-quadrupole-term
END DO

IF (dipole == 1) THEN
  DO i = 1,ncomp
    DO j = 1,ncomp
      idd2(i,j)   =0.0
      idd4(i,j)   =0.0
      idd2dt(i,j) =0.0
      idd4dt(i,j) =0.0
      idd2dt2(i,j)=0.0
      idd4dt2(i,j)=0.0
      DO m=0,4
        idd2(i,j)  = idd2(i,j)   +ddp2(i,j,m)*z3**REAL(m)
        idd4(i,j)  = idd4(i,j)   +ddp4(i,j,m)*z3**REAL(m)
        idd2dt(i,j)= idd2dt(i,j) +ddp2(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
        idd4dt(i,j)= idd4dt(i,j) +ddp4(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
        idd2dt2(i,j)=idd2dt2(i,j)+ddp2(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                     + ddp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
        idd4dt2(i,j)=idd4dt2(i,j)+ddp4(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                     + ddp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
      END DO
      DO k = 1,ncomp
        idd3(i,j,k)   =0.0
        idd3dt(i,j,k) =0.0
        idd3dt2(i,j,k)=0.0
        DO m = 0, 4
          idd3(i,j,k)   = idd3(i,j,k)  +ddp3(i,j,k,m)*z3**REAL(m)
          idd3dt(i,j,k) = idd3dt(i,j,k)+ddp3(i,j,k,m)*z3dt*REAL(m) *z3**REAL(m-1)
          idd3dt2(i,j,k)= idd3dt2(i,j,k)+ddp3(i,j,k,m)*z3dt2*REAL(m)  &
              *z3**REAL(m-1) +ddp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1) *z3**REAL(m-2)
        END DO
      END DO
    END DO
  END DO
  
  
  factor2= -PI *rho
  factor3= -4.0/3.0*PI**2 * rho**2 
  
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


! ----------------------------------------------------------------------
! derivatives f/kT of quadrupole-quadrup. term to T  (fqqdt)
! ----------------------------------------------------------------------
fqqdt  = 0.0
fqqdt2 = 0.0
qudpole = 0
DO i = 1, ncomp
  qq2(i) = (parame(i,7))**2 *1.E-69 / (uij(i,i)*kbol*mseg(i)*sig_ij(i,i)**5 *1.E-50)
  IF (qq2(i) /= 0.0) qudpole = 1
END DO

IF (qudpole == 1) THEN
  
  DO i = 1,ncomp
    DO j = 1,ncomp
      iqq2(i,j)   = 0.0
      iqq4(i,j)   = 0.0
      iqq2dt(i,j) = 0.0
      iqq4dt(i,j) = 0.0
      iqq2dt2(i,j)= 0.0
      iqq4dt2(i,j)= 0.0
      DO m = 0, 4
        iqq2(i,j)   = iqq2(i,j)  + qqp2(i,j,m)*z3**REAL(m)
        iqq4(i,j)   = iqq4(i,j)  + qqp4(i,j,m)*z3**REAL(m)
        iqq2dt(i,j) = iqq2dt(i,j)+ qqp2(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
        iqq4dt(i,j) = iqq4dt(i,j)+ qqp4(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
        iqq2dt2(i,j)= iqq2dt2(i,j)+qqp2(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                      + qqp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
        iqq4dt2(i,j)= iqq4dt2(i,j)+qqp4(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                      + qqp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
      END DO
      DO k = 1,ncomp
        iqq3(i,j,k)   =0.0
        iqq3dt(i,j,k) =0.0
        iqq3dt2(i,j,k)=0.0
        DO m = 0, 4
          iqq3(i,j,k)   = iqq3(i,j,k)  + qqp3(i,j,k,m)*z3**REAL(m)
          iqq3dt(i,j,k) = iqq3dt(i,j,k)+ qqp3(i,j,k,m)*z3dt*REAL(m) * z3**REAL(m-1)
          iqq3dt2(i,j,k)= iqq3dt2(i,j,k)+qqp3(i,j,k,m)*z3dt2*REAL(m) * z3**REAL(m-1)  &
                      + qqp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
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
  DO i = 1,ncomp
    DO j = 1,ncomp
      xijmt = x(i)*uij(i,i)*qq2(i)*sig_ij(i,i)**5   &
              * x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /sig_ij(i,j)**7.0
      eij = (parame(i,3)*parame(j,3))**0.5
      fqq2  = fqq2   +factor2* xijmt/t/t*(iqq2(i,j)+eij/t*iqq4(i,j))
      fqq2dt= fqq2dt +factor2* xijmt/t/t*(iqq2dt(i,j)-2.0*iqq2(i,j)/t  &
              + eij/t*iqq4dt(i,j)-3.0*eij/t/t*iqq4(i,j))
      fqq2dt2=fqq2dt2+factor2*xijmt/t/t*(iqq2dt2(i,j)-4.0*iqq2dt(i,j)/t  &
              + 6.0*iqq2(i,j)/t/t+eij/t*iqq4dt2(i,j)  &
              - 6.0*eij/t/t*iqq4dt(i,j)+12.0*eij/t**3 *iqq4(i,j))
      DO k = 1,ncomp
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


! ----------------------------------------------------------------------
! derivatives f/kT of dipole-quadruppole term to T  (fdqdt)
! ----------------------------------------------------------------------
fdqdt = 0.0
fdqdt2= 0.0
dip_quad = 0
DO i = 1,ncomp
  DO j = 1,ncomp
    IF (parame(i,6) /= 0.0 .AND. parame(j,7) /= 0.0) dip_quad = 1
  END DO
  myfac(i) = parame(i,3)*parame(i,2)**4 *my0(i)
  q_fac(i) = parame(i,3)*parame(i,2)**4 *qq2(i)
END DO

IF (dip_quad == 1) THEN
  
  DO i = 1,ncomp
    DO j = 1,ncomp
      idq2(i,j)   = 0.0
      idq4(i,j)   = 0.0
      idq2dt(i,j) = 0.0
      idq4dt(i,j) = 0.0
      idq2dt2(i,j)= 0.0
      idq4dt2(i,j)= 0.0
      IF ( myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0 ) THEN
        DO m = 0, 4
          idq2(i,j)   = idq2(i,j)  + dqp2(i,j,m)*z3**REAL(m)
          idq4(i,j)   = idq4(i,j)  + dqp4(i,j,m)*z3**REAL(m)
          idq2dt(i,j) = idq2dt(i,j)+ dqp2(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
          idq4dt(i,j) = idq4dt(i,j)+ dqp4(i,j,m)*z3dt*REAL(m)*z3**REAL(m-1)
          idq2dt2(i,j)= idq2dt2(i,j)+dqp2(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                        + dqp2(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
          idq4dt2(i,j)= idq4dt2(i,j)+dqp4(i,j,m)*z3dt2*REAL(m)*z3**REAL(m-1)  &
                        + dqp4(i,j,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
        END DO
        
        DO k = 1,ncomp
          idq3(i,j,k)   = 0.0
          idq3dt(i,j,k) = 0.0
          idq3dt2(i,j,k)= 0.0
          IF ( myfac(k) /= 0.0 .OR. q_fac(k) /= 0.0 ) THEN
            DO m = 0, 4
              idq3(i,j,k)  = idq3(i,j,k)  + dqp3(i,j,k,m)*z3**REAL(m)
              idq3dt(i,j,k)= idq3dt(i,j,k)+ dqp3(i,j,k,m)*z3dt*REAL(m) *z3**REAL(m-1)
              idq3dt2(i,j,k)= idq3dt2(i,j,k)+dqp3(i,j,k,m)*z3dt2*REAL(m) *z3**REAL(m-1)  &
                              + dqp3(i,j,k,m)*z3dt**2 *REAL(m)*REAL(m-1)*z3**REAL(m-2)
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
  DO i = 1,ncomp
    DO j = 1,ncomp
      IF ( myfac(i) /= 0.0 .AND. q_fac(j) /= 0.0 ) THEN
        xijmt = x(i)*myfac(i) * x(j)*q_fac(j) /sig_ij(i,j)**5 
        eij = (parame(i,3)*parame(j,3))**0.5
        fdq2  = fdq2  + factor2* xijmt/t/t*(idq2(i,j)+eij/t*idq4(i,j))
        fdq2dt= fdq2dt+ factor2* xijmt/t/t*(idq2dt(i,j)-2.0*idq2(i,j)/t  &
                + eij/t*idq4dt(i,j)-3.0*eij/t/t*idq4(i,j))
        fdq2dt2 = fdq2dt2+factor2*xijmt/t/t*(idq2dt2(i,j)-4.0*idq2dt(i,j)/t  &
                + 6.0*idq2(i,j)/t/t+eij/t*idq4dt2(i,j)  &
                - 6.0*eij/t/t*idq4dt(i,j)+12.0*eij/t**3 *idq4(i,j))
        DO k = 1,ncomp
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
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
! total derivative of fres/kT to temperature
! ----------------------------------------------------------------------

df_dt = fhsdt + fchdt + fdspdt + fhbdt + fdddt + fqqdt + fdqdt



! ----------------------------------------------------------------------
! second derivative of fres/kT to T
! ----------------------------------------------------------------------

df_dt2 = fhsdt2 +fchdt2 +fdspdt2 +fhbdt2 +fdddt2 +fqqdt2 +fdqdt2



! ----------------------------------------------------------------------
! ------ derivatives of fres/kt to density and to T --------------------
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! the analytic derivative of fres/kT to (density and T) (df_drdt)
! is still missing. A numerical differentiation is implemented.
! ----------------------------------------------------------------------
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





! ----------------------------------------------------------------------
! thermodynamic properties
! ----------------------------------------------------------------------

s_res = ( - df_dt *t - fres )*RGAS + RGAS * LOG(zges)
h_res = ( - t*df_dt  +  zges-1.0  ) * RGAS *t
cv_res = - (   t*df_dt2 + 2.0*df_dt  ) * RGAS*t
cp_res = cv_res - RGAS  + RGAS*(zges +eta*t*df_drdt)**2   &
                          / (1.0 + 2.0*eta*dfdr +eta**2  *ddfdrdr)

! write (*,*) 'df_... ', df_dt,df_dt2
! write (*,*) 'kreuz ', zges,eta*t*df_drdt,eta*dfdr, eta**2 *ddfdrdr
! write (*,*) 'h,cv,cp', h_res,cv_res,cp_res


END SUBROUTINE H_EOS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE H_EOS_num
!
! This subroutine calculates enthalpies and heat capacities (cp) by
! taking numerical derivatieves.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE H_EOS_num
!
 USE PARAMETERS, ONLY: RGAS
 USE EOS_VARIABLES
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL                                   :: dist, fact, rho_0
 REAL                                   :: fres1, fres2, fres3, fres4, fres5
 REAL                                   :: f_1, f_2, f_3, f_4
 REAL                                   :: cv_res, t_tmp, zges
 REAL                                   :: df_dt, df_dtdt, df_drdt, dfdr, ddfdrdr

!-----------------------------------------------------------------------


CALL PERTURBATION_PARAMETER
rho_0 = eta/z3t


fact = 1.0
dist = t * 100.E-5 * fact

t_tmp  = t

t  = t_tmp - 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_EOS
fres1  = fres
t  = t_tmp - dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_EOS
fres2  = fres
t  = t_tmp + dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_EOS
fres3  = fres
t  = t_tmp + 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_EOS
fres4  = fres
t  = t_tmp
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
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
eta = z3t*rho_0
CALL P_EOS
f_1  = pges/(eta*rho_0*(kbol*t)/1.E-30)

t  = t_tmp - dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_EOS
f_2  = pges/(eta*rho_0*(kbol*t)/1.E-30)

t  = t_tmp + dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_EOS
f_3  = pges/(eta*rho_0*(kbol*t)/1.E-30)

t  = t_tmp + 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_EOS
f_4  = pges/(eta*rho_0*(kbol*t)/1.E-30)

t  = t_tmp
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_EOS

dfdr    = pges / (eta*rho_0*(kbol*t)/1.E-30)
ddfdrdr = pgesdz/(eta*rho_0*(kbol*t)/1.E-30) - dfdr*2.0/eta - 1.0/eta**2 

df_drdt   = ( -f_4 +8.0*f_3 -8.0*f_2 +f_1) / (12.0*dist)

cp_res = cv_res - RGAS +RGAS*(zges+eta*t*df_drdt)**2   &
                  * 1.0/(1.0 + 2.0*eta*dfdr + eta**2  *ddfdrdr)

! write (*,*) 'n',df_dt,df_dtdt
! write (*,*) 'kreuz ', zges,eta*t*df_drdt,eta*dfdr, eta**2 *ddfdrdr
! write (*,*) 'h, cv', h_res, cv_res
! write (*,*) h_res - t*s_res
! write (*,*) cv_res,cp_res,eta
! pause

END SUBROUTINE H_EOS_num


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DENSITY_ITERATION
!
! iterates the density until the calculated pressure 'pges' is equal to
! the specified pressure 'p'. A Newton-scheme is used for determining
! the root to the objective function  f(eta) = (pges / p ) - 1.0. 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE DENSITY_ITERATION
!
 USE BASIC_VARIABLES, ONLY: num
 USE EOS_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, start, max_i
 REAL                                   :: eta_iteration
 REAL                                   :: error, dydx, acc_i, delta_eta
! ----------------------------------------------------------------------


IF ( densav(phas) /= 0.0 .AND. eta_start == denold(phas) ) THEN
  denold(phas) = eta_start
  eta_start = densav(phas)
ELSE
  denold(phas) = eta_start
  densav(phas) = eta_start
END IF


acc_i = 1.d-9
max_i = 30
density_error(:) = 0.0

i = 0
eta_iteration = eta_start

! ----------------------------------------------------------------------
! iterate density until p_calc = p
! ----------------------------------------------------------------------
iterate_density: DO

  i = i + 1
  eta = eta_iteration

  IF ( num == 0 ) THEN
     CALL P_EOS 
  ELSE IF ( num == 1 ) THEN
    CALL P_NUMERICAL
  ELSE IF ( num == 2 ) THEN
   WRITE(*,*) 'CRITICAL RENORM NOT INCLUDED YET' 
   !!CALL F_EOS_RN
  ELSE
    write (*,*) 'define calculation option (num)'
  END IF

  error = (pges / p ) - 1.0

  ! --- instable region correction -------------------------------------
  IF ( pgesdz < 0.0 .AND. i < max_i ) THEN
    IF ( error > 0.0 .AND. pgesd2 > 0.0 ) THEN                           ! no liquid density
      CALL PRESSURE_SPINODAL
      eta_iteration = eta
      error  = (pges / p ) - 1.0
      IF ( ((pges/p ) -1.0) > 0.0 ) eta_iteration = 0.001                ! no solution possible
      IF ( ((pges/p ) -1.0) <=0.0 ) eta_iteration = eta_iteration * 1.1  ! no solution found so far
    ELSE IF ( error < 0.0 .AND. pgesd2 < 0.0 ) THEN                      ! no vapor density
      CALL PRESSURE_SPINODAL
      eta_iteration = eta
      error  = (pges / p ) - 1.0
      IF ( ((pges/p ) -1.0) < 0.0 ) eta_iteration = 0.5                  ! no solution possible
      IF ( ((pges/p ) -1.0) >=0.0 ) eta_iteration = eta_iteration * 0.9  ! no solution found so far
    ELSE
      eta_iteration = (eta_iteration + eta_start) / 2.0
      IF (eta_iteration == eta_start) eta_iteration = eta_iteration + 0.2
    END IF
    CYCLE iterate_density
  END IF


  dydx = pgesdz/p
  delta_eta = error/ dydx
  IF ( delta_eta >  0.05 ) delta_eta = 0.05
  IF ( delta_eta < -0.05 ) delta_eta = -0.05

  eta_iteration   = eta_iteration - delta_eta

  IF (eta_iteration > 0.9)  eta_iteration = 0.6
  IF (eta_iteration <= 0.0) eta_iteration = 1.E-16
  start = 1

  IF ( ABS(error*p/pgesdz) < 1.d-12 ) start = 0
  IF ( ABS(error) < acc_i ) start = 0
  IF ( i > max_i ) THEN
    start = 0
    density_error(phas) = ABS( error )
    ! write (*,*) 'density iteration failed'
  END IF

  IF (start /= 1) EXIT iterate_density

END DO iterate_density

eta = eta_iteration

IF ((eta > 0.3 .AND. densav(phas) > 0.3) .OR.  &
    (eta < 0.1 .AND. densav(phas) < 0.1)) densav(phas) = eta

END SUBROUTINE DENSITY_ITERATION




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE F_EOS
!
! calculates the Helmholtz energy f/kT. The input to the subroutine is
! (T,eta,x), where eta is the packing fraction.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_EOS
!
 USE PARAMETERS, ONLY: nc, nsite
 USE EOS_VARIABLES
 USE EOS_CONSTANTS
 IMPLICIT NONE
!
! --- local variables --------------------------------------------------
 INTEGER                                :: i, j, k, l, m, n
 REAL                                   :: z0, z1, z2, z3
 REAL                                   :: zms, m_mean   ! ,lij(nc,nc)
 REAL                                   :: I1,I2, c1_con
 REAL                                   :: fhs, fdsp, fhc

 LOGICAL                                :: assoc
 INTEGER                                :: ass_cnt,max_eval
 REAL                                   :: delta(nc,nc,nsite,nsite)
 REAL                                   :: mx_itr(nc,nsite), err_sum, sum, attenu, tol, fhb
 REAL                                   :: ass_s1, ass_s2

 REAL                                   :: fdd, fqq, fdq
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! abbreviations
! ----------------------------------------------------------------------
rho = eta/z3t
z0 = z0t*rho
z1 = z1t*rho
z2 = z2t*rho
z3 = z3t*rho

m_mean = z0t / ( PI / 6.0 )
zms    = 1.0 - eta

! m_mean2  = 0.0
! lij(1,2) = -0.05
! lij(2,1) = lij(1,2)
! DO i = 1, ncomp
!    DO j = 1, ncomp
!       m_mean2 = m_mean2 + x(i)*x(j)*(mseg(i)+mseg(j))/2.0*(1.0-lij(i,j))
!    ENDDO
! ENDDO


! ----------------------------------------------------------------------
! radial distr. function at contact,  gij
! ----------------------------------------------------------------------
DO  i = 1, ncomp
  DO  j=1,ncomp
    gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
  END DO
END DO


! ----------------------------------------------------------------------
! Helmholtz energy : hard sphere contribution
! ----------------------------------------------------------------------
fhs= m_mean*(  3.0*z1*z2/zms + z2**3 /z3/zms/zms + (z2**3 /z3/z3-z0)*LOG(zms)  )/z0


! ----------------------------------------------------------------------
! Helmholtz energy : chain term
! ----------------------------------------------------------------------
fhc = 0.0
DO i = 1, ncomp
  fhc = fhc + x(i) *(1.0- mseg(i)) *LOG(gij(i,i))
END DO


! ----------------------------------------------------------------------
! Helmholtz energy : PC-SAFT dispersion contribution
! ----------------------------------------------------------------------
IF (eos == 1) THEN
  
  I1 = 0.0
  I2 = 0.0
  DO m = 0, 6
    I1 = I1 + apar(m)*eta**REAL(m)
    I2 = I2 + bpar(m)*eta**REAL(m)
  END DO
  
  c1_con= 1.0/ (  1.0 + m_mean*(8.0*eta-2.0*eta**2 )/zms**4   &
          + (1.0 - m_mean)*(20.0*eta-27.0*eta**2   &
           + 12.0*eta**3 -2.0*eta**4 ) /(zms*(2.0-eta))**2   )
  
  fdsp  = -2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2
  
! ----------------------------------------------------------------------
! Helmholtz energy : SAFT (Chen, Kreglewski) dispersion contribution
! ----------------------------------------------------------------------
ELSE
  
  fdsp  = 0.0
  DO  n = 1,4
    DO  m = 1,9
      fdsp =  fdsp + dnm(n,m) * (um/t)**REAL(n) *(eta/tau)**REAL(m)
    END DO
  END DO
  fdsp = m_mean * fdsp

END IF


! ----------------------------------------------------------------------
! TPT-1-association according to Chapman et al.
! ----------------------------------------------------------------------
fhb = 0.0
assoc = .false.
DO i = 1, ncomp
  IF (nhb_typ(i) /= 0) assoc = .true.
END DO
IF (assoc) THEN
  
  DO i = 1, ncomp
    DO k = 1, nhb_typ(i)
      IF (mx(i,k) == 0.0) mx(i,k) = 1.0        !  Initialize mx(i,j)
      DO j = 1, ncomp
        DO l = 1, nhb_typ(j)
          delta(i,j,k,l) = gij(i,j) * ass_d(i,j,k,l)
        END DO
      END DO
    END DO
  END DO
  
  
! --- constants for iteration ------------------------------------------
  attenu = 0.70
  tol = 1.d-10
  IF (eta < 0.2)  tol = 1.d-12
  IF (eta < 0.01) tol = 1.d-13
  max_eval = 200
  
! --- iterate over all components and all sites ------------------------
  ass_cnt = 0
  iterate_TPT1: DO

    ass_cnt = ass_cnt + 1
    
    DO i = 1, ncomp
      DO k = 1, nhb_typ(i)
        sum = 0.0
        DO j = 1, ncomp
          DO l = 1, nhb_typ(j)
            sum = sum +  x(j)* mx(j,l)*nhb_no(j,l) *delta(i,j,k,l)
!            if (ass_cnt == 1) write (*,*) j,l,x(j), mx(j,l)
          END DO
        END DO
        mx_itr(i,k) = 1.0 / (1.0 + sum * rho)
!        if (ass_cnt == 1) write (*,*) 'B',ass_cnt,sum, rho
      END DO
    END DO
    
    err_sum = 0.0
    DO i = 1, ncomp
      DO k = 1, nhb_typ(i)
        err_sum = err_sum + ABS(mx_itr(i,k) - mx(i,k))    ! / ABS(mx_itr(i,k))
        mx(i,k) = mx_itr(i,k) * attenu + mx(i,k) * (1.0 - attenu)
      END DO
    END DO

    IF ( err_sum <= tol .OR. ass_cnt >= max_eval ) THEN
      IF (ass_cnt >= max_eval) WRITE(*,'(a,2G15.7)') 'F_EOS: Max_eval violated (mx) Err_Sum = ',err_sum,tol
      EXIT iterate_TPT1
    END IF
    
  END DO iterate_TPT1
  
  DO i = 1, ncomp
    ass_s1  = 0.0
    ass_s2  = 0.0
    DO k = 1, nhb_typ(i)
      ass_s1  = ass_s1  + nhb_no(i,k) * ( 1.0 - mx(i,k) )
      ass_s2  = ass_s2  + nhb_no(i,k) * LOG( mx(i,k) )
    END DO
    fhb = fhb + x(i) * ( ass_s2 + ass_s1 / 2.0 )
  END DO
  
END IF
! --- TPT-1-association accord. to Chapman -----------------------------


! ----------------------------------------------------------------------
! polar terms
! ----------------------------------------------------------------------
 CALL F_POLAR ( fdd, fqq, fdq )


! ----------------------------------------------------------------------
! resid. Helmholtz energy f/kT
! ----------------------------------------------------------------------
fres = fhs + fhc + fdsp + fhb + fdd + fqq + fdq

tfr= fres

END SUBROUTINE F_EOS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_NUMERICAL
!
 USE EOS_VARIABLES
 USE EOS_CONSTANTS
 USE EOS_NUMERICAL_DERIVATIVES, ONLY: ideal_gas, hard_sphere, chain_term, &
                                      disp_term, hb_term, LC_term, branch_term,  &
                                      II_term, ID_term, subtract1, subtract2
 IMPLICIT NONE
!
!---------------------------------------------------------------------

!-----local variables-------------------------------------------------
 INTEGER                                :: i, j
 REAL                                   :: m_mean2
 REAL                                   :: fid, fhs, fdsp, fhc
 REAL                                   :: fhb, fdd, fqq, fdq
 REAL                                   :: fhend, fcc
 REAL                                   :: fbr, flc
 REAL                                   :: fref

 REAL                                   :: eps_kij, k_kij
!---------------------------------------------------------------------

eps_kij = 0.0
k_kij   = 0.0

fid = 0.0
fhs = 0.0
fhc = 0.0
fdsp= 0.0
fhb = 0.0
fdd = 0.0
fqq = 0.0
fdq = 0.0
fcc = 0.0
fbr = 0.0
flc = 0.0


CALL PERTURBATION_PARAMETER

! ----------------------------------------------------------------------
! overwrite the standard mixing rules by those published by Tang & Gross
! using an additional lij parameter
! WARNING : the lij parameter is set to lij = - lji in 'para_input'
! ----------------------------------------------------------------------
order1 = 0.0
order2 = 0.0
DO i = 1, ncomp
  DO j = 1, ncomp
    order1 = order1 + x(i)*x(j)* mseg(i)*mseg(j) * sig_ij(i,j)**3 * uij(i,j)/t
    order2 = order2 + x(i)*x(j)* mseg(i)*mseg(j) * sig_ij(i,j)**3 * (uij(i,j)/t)**2
  END DO
END DO
DO i = 1, ncomp
  DO j = 1, ncomp
    order1 = order1 + x(i)*mseg(i)/t*( x(j)*mseg(j)  &
                      *sig_ij(i,j)*(uij(i,i)*uij(j,j))**(1.0/6.0) )**3 *lij(i,j)
  END DO
END DO


! ----------------------------------------------------------------------
! a non-standard mixing rule scaling the hard-sphere term
! WARNING : the lij parameter is set to lij = - lji in 'para_input'
! (uses an additional lij parameter)
! ----------------------------------------------------------------------
m_mean2 = 0.0
DO i = 1, ncomp
  DO j = 1, ncomp
    m_mean2 = m_mean2 + x(i)*x(j)*(mseg(i)+mseg(j))/2.0
  END DO
END DO
DO i = 1, ncomp
  DO j = 1, ncomp
    ! m_mean2=m_mean2+x(i)*(x(j)*((mseg(i)+mseg(j))*0.5)**(1.0/3.0) *lij(i,j) )**3
  END DO
END DO

! --- ideal gas contribution -------------------------------------------
IF ( ideal_gas == 'yes' )       CALL F_IDEAL_GAS ( fid )
! ----------------------------------------------------------------------

! --- hard-sphere contribution -----------------------------------------
IF ( hard_sphere == 'CSBM' )    CALL F_HARD_SPHERE ( m_mean2, fhs )
! ----------------------------------------------------------------------

! -- chain term --------------------------------------------------------
IF ( chain_term == 'TPT1' )     CALL F_CHAIN_TPT1 ( fhc )
IF ( chain_term == 'TPT2' )     CALL F_CHAIN_TPT_D ( fhc )
IF ( chain_term == 'HuLiu' )    CALL F_CHAIN_HU_LIU ( fhc )
IF ( chain_term == 'HuLiu_rc' ) CALL F_CHAIN_HU_LIU_RC ( fhs, fhc )
!!IF ( chain_term == 'SPT' )      CALL F_SPT ( fhs, fhc )
IF ( chain_term == 'SPT' )    WRITE(*,*) 'SPT NOT INCLUDED YET'
! ----------------------------------------------------------------------

! --- dispersive attraction --------------------------------------------
IF ( disp_term == 'PC-SAFT')    CALL F_DISP_PCSAFT ( fdsp )
IF ( disp_term == 'CK')         CALL F_DISP_CKSAFT ( fdsp )
IF ( disp_term(1:2) == 'PT')    CALL F_pert_theory ( fdsp )
! ----------------------------------------------------------------------

! --- H-bonding contribution / Association -----------------------------
IF ( hb_term == 'TPT1_Chap')    CALL F_ASSOCIATION( eps_kij, k_kij, fhb )
! ----------------------------------------------------------------------

! --- polar terms ------------------------------------------------------
 CALL F_POLAR ( fdd, fqq, fdq )
! ----------------------------------------------------------------------

! --- ion-dipole term --------------------------------------------------
IF ( ID_term == 'TBH')          CALL F_ION_DIPOLE_TBH ( fhend )
! ----------------------------------------------------------------------

! --- ion-ion term -----------------------------------------------------
IF ( II_term == 'primMSA')      CALL F_ION_ION_PrimMSA ( fcc )
IF ( II_term == 'nprMSA')       CALL F_ION_ION_nonPrimMSA ( fdd, fqq, fdq, fcc )
! ----------------------------------------------------------------------

! --- liquid-crystal term ----------------------------------------------
IF ( LC_term == 'MSaupe')       CALL F_LC_MayerSaupe ( flc )

!!IF ( LC_term == 'OVL')          fref = fhs + fhc
IF ( LC_term == 'OVL') WRITE(*,*) 'OVL NOT INCLUDED YET'
!IF ( LC_term == 'OVL')          CALL F_LC_OVL ( fref, flc )
!! IF ( LC_term == 'SPT')          fref = fhs + fhc
IF ( LC_term == 'SPT') WRITE(*,*) 'SPT NOT INCLUDED YET'
!!IF ( LC_term == 'SPT')          CALL F_LC_SPT( fref, flc )
! ----------------------------------------------------------------------

! ======================================================================
! SUBTRACT TERMS (local density approximation) FOR DFT
! ======================================================================

!IF ( subtract1 == '1PT')     CALL F_subtract_local_pert_theory ( subtract1, fdsp )
!IF ( subtract1 == '2PT')     CALL F_subtract_local_pert_theory ( subtract1, fdsp )
!IF ( subtract2 =='ITTpolar') CALL F_local_ITT_polar ( fdd )
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! residual Helmholtz energy F/(NkT)
! ----------------------------------------------------------------------
fres = fid + fhs + fhc + fdsp + fhb + fdd + fqq + fdq + fcc + flc

tfr = 0.0

END SUBROUTINE F_NUMERICAL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE P_EOS
!
! calculates the pressure in units (Pa).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE P_EOS
!
! ----------------------------------------------------------------------
 USE PARAMETERS, ONLY: nc, nsite
 USE EOS_VARIABLES
 USE EOS_CONSTANTS
 IMPLICIT NONE
!
! --- local variables --------------------------------------------------
 INTEGER                                :: i, j, k, l, m, n
 INTEGER                                :: ass_cnt,max_eval
 LOGICAL                                :: assoc
 REAL                                   :: z0, z1, z2, z3
 REAL                                   :: zms, m_mean
 REAL                                   :: zges, zgesdz, zgesd2, zgesd3
 REAL                                   :: zhs, zhsdz, zhsd2, zhsd3
 REAL                                   :: zhc, zhcdz, zhcd2, zhcd3
 REAL, DIMENSION(nc,nc)                 :: dgijdz, dgijd2, dgijd3, dgijd4
 REAL                                   :: zdsp, zdspdz, zdspd2, zdspd3
 REAL                                   :: c1_con, c2_con, c3_con, c4_con, c5_con
 REAL                                   :: I2, edI1dz, edI2dz, edI1d2, edI2d2
 REAL                                   :: edI1d3, edI2d3, edI1d4, edI2d4
 REAL                                   :: fdspdz,fdspd2
 REAL                                   :: zhb, zhbdz, zhbd2, zhbd3
 REAL, DIMENSION(nc,nc,nsite,nsite)     :: delta, dq_dz, dq_d2, dq_d3, dq_d4
 REAL, DIMENSION(nc,nsite)              :: mx_itr, dmx_dz, ndmxdz, dmx_d2, ndmxd2
 REAL, DIMENSION(nc,nsite)              :: dmx_d3, ndmxd3, dmx_d4, ndmxd4
 REAL                                   :: err_sum, sum0, sum1, sum2, sum3, sum4, attenu, tol
 REAL                                   :: sum_d1, sum_d2, sum_d3, sum_d4
 REAL                                   :: zdd, zddz, zddz2, zddz3
 REAL                                   :: zqq, zqqz, zqqz2, zqqz3
 REAL                                   :: zdq, zdqz, zdqz2, zdqz3
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! abbreviations
! ----------------------------------------------------------------------
rho = eta/z3t
z0 = z0t*rho
z1 = z1t*rho
z2 = z2t*rho
z3 = z3t*rho

m_mean = z0t/(PI/6.0)
zms    = 1.0 -eta

! m_mean2=0.0
! lij(1,2)= -0.050
! lij(2,1)=lij(1,2)
! DO i =1,ncomp
!   DO j =1,ncomp
!     m_mean2=m_mean2+x(i)*x(j) * (mseg(i)+mseg(j))/2.0*(1.0-lij(i,j))
!   ENDDO
! ENDDO


! ----------------------------------------------------------------------
! radial distr. function at contact,  gij , and derivatives
! dgijdz=d(gij)/d(eta)   and   dgijd2 = dd(gij)/d(eta)**2
! ----------------------------------------------------------------------
DO  i = 1, ncomp
  DO  j=1,ncomp
    ! j=i
    gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms + 2.0*(dij_ab(i,j)*z2)**2 /zms**3 
    dgijdz(i,j)= 1.0/zms/zms + 3.0*dij_ab(i,j)*z2*(1.0+z3)/z3/zms**3   &
        + (dij_ab(i,j)*z2/zms/zms)**2 *(4.0+2.0*z3)/z3
    dgijd2(i,j) = 2.0/zms**3   &
        + 6.0*dij_ab(i,j)*z2/z3/zms**4 *(2.0+z3)  &
        + (2.0*dij_ab(i,j)*z2/z3)**2 /zms**5  *(1.0+4.0*z3+z3*z3)
    dgijd3(i,j) = 6.0/zms**4   &
        + 18.0*dij_ab(i,j)*z2/z3/zms**5 *(3.0+z3)  &
        + 12.0*(dij_ab(i,j)*z2/z3/zms**3 )**2  *(3.0+6.0*z3+z3*z3)
    dgijd4(i,j) = 24.0/zms**5   &
        + 72.0*dij_ab(i,j)*z2/z3/zms**6 *(4.0+z3)  &
        + 48.0*(dij_ab(i,j)*z2/z3)**2 /zms**7  *(6.0+8.0*z3+z3*z3)
  END DO
END DO


! ----------------------------------------------------------------------
! p : hard sphere contribution
! ----------------------------------------------------------------------
zhs   = m_mean* ( z3/zms + 3.0*z1*z2/z0/zms/zms + z2**3 /z0*(3.0-z3)/zms**3 )
zhsdz = m_mean*(  1.0/zms/zms + 3.0*z1*z2/z0/z3*(1.0+z3)/zms**3   &
    + 6.0*z2**3 /z0/z3/zms**4  )
zhsd2 = m_mean*(  2.0/zms**3  + 6.0*z1*z2/z0/z3*(2.0+z3)/zms**4   &
    + 6.0*z2**3 /z0/z3/z3*(1.0+3.0*z3)/zms**5  )
zhsd3 = m_mean*(  6.0/zms**4  + 18.0*z1*z2/z0/z3*(3.0+z3)/zms**5   &
    + 24.0*z2**3 /z0/z3/z3*(2.0+3.0*z3)/zms**6  )


! ----------------------------------------------------------------------
! p : chain term
! ----------------------------------------------------------------------
zhc   = 0.0
zhcdz = 0.0
zhcd2 = 0.0
zhcd3 = 0.0
DO i= 1, ncomp
  zhc = zhc + x(i)*(1.0-mseg(i))*eta/gij(i,i)* dgijdz(i,i)
  zhcdz = zhcdz + x(i)*(1.0-mseg(i)) *(-eta*(dgijdz(i,i)/gij(i,i))**2   &
      + dgijdz(i,i)/gij(i,i) + eta/gij(i,i)*dgijd2(i,i))
  zhcd2 = zhcd2 + x(i)*(1.0-mseg(i))  &
      *( 2.0*eta*(dgijdz(i,i)/gij(i,i))**3   &
      -2.0*(dgijdz(i,i)/gij(i,i))**2   &
      -3.0*eta/gij(i,i)**2 *dgijdz(i,i)*dgijd2(i,i)  &
      +2.0/gij(i,i)*dgijd2(i,i) +eta/gij(i,i)*dgijd3(i,i) )
  zhcd3 = zhcd3 + x(i)*(1.0-mseg(i)) *( 6.0*(dgijdz(i,i)/gij(i,i))**3   &
      -6.0*eta*(dgijdz(i,i)/gij(i,i))**4   &
      +12.0*eta/gij(i,i)**3 *dgijdz(i,i)**2 *dgijd2(i,i)  &
      -9.0/gij(i,i)**2 *dgijdz(i,i)*dgijd2(i,i) +3.0/gij(i,i)*dgijd3(i,i)  &
      -3.0*eta*(dgijd2(i,i)/gij(i,i))**2   &
      -4.0*eta/gij(i,i)**2 *dgijdz(i,i)*dgijd3(i,i)  &
      +eta/gij(i,i)*dgijd4(i,i) )
END DO

! ----------------------------------------------------------------------
! p : PC-SAFT dispersion contribution
!     note: edI1dz is equal to d(eta*I1)/d(eta), analogous for edI2dz
! ----------------------------------------------------------------------
IF (eos == 1) THEN
  
  I2     = 0.0
  edI1dz = 0.0
  edI2dz = 0.0
  edI1d2 = 0.0
  edI2d2 = 0.0
  edI1d3 = 0.0
  edI2d3 = 0.0
  edI1d4 = 0.0
  edI2d4 = 0.0
  DO  m=0,6
    I2    = I2 + bpar(m)*z3**REAL(m)
    edI1dz= edI1dz+apar(m)*REAL(m+1)*z3**REAL(m)
    edI2dz= edI2dz+bpar(m)*REAL(m+1)*z3**REAL(m)
    edI1d2= edI1d2+apar(m)*REAL((m+1)*m)*z3**REAL(m-1)
    edI2d2= edI2d2+bpar(m)*REAL((m+1)*m)*z3**REAL(m-1)
    edI1d3= edI1d3+apar(m)*REAL((m+1)*m*(m-1))*z3**REAL(m-2)
    edI2d3= edI2d3+bpar(m)*REAL((m+1)*m*(m-1))*z3**REAL(m-2)
    edI1d4= edI1d4+apar(m)*REAL((m+1)*m*(m-1)*(m-2))*z3**REAL(m-3)
    edI2d4= edI2d4+bpar(m)*REAL((m+1)*m*(m-1)*(m-2))*z3**REAL(m-3)
  END DO
  
  c1_con= 1.0/ (  1.0 + m_mean*(8.0*eta-2.0*eta**2 )/zms**4   &
      + (1.0 - m_mean)*(20.0*eta-27.0*eta**2   &
      + 12.0*eta**3 -2.0*eta**4 ) /(zms*(2.0-eta))**2   )
  c2_con= - c1_con*c1_con  &
      *(m_mean*(-4.0*eta**2 +20.0*eta+8.0)/zms**5  + (1.0 - m_mean)  &
      *(2.0*eta**3 +12.0*eta**2 -48.0*eta+40.0)  &
      /(zms*(2.0-eta))**3  )
  c3_con= 2.0 * c2_con*c2_con/c1_con - c1_con*c1_con  &
      *( m_mean*(-12.0*eta**2 +72.0*eta+60.0)/zms**6   &
      + (1.0 - m_mean)  &
      *(-6.0*eta**4 -48.0*eta**3 +288.0*eta**2   &
      -480.0*eta+264.0) /(zms*(2.0-eta))**4  )
  c4_con= 6.0*c2_con*c3_con/c1_con -6.0*c2_con**3 /c1_con**2   &
      - c1_con*c1_con  &
      *( m_mean*(-48.0*eta**2 +336.0*eta+432.0)/zms**7   &
      + (1.0 - m_mean)  &
      *(24.0*eta**5 +240.0*eta**4 -1920.0*eta**3   &
      +4800.0*eta**2 -5280.0*eta+2208.0) /(zms*(2.0-eta))**5  )
  c5_con= 6.0*c3_con**2 /c1_con - 36.0*c2_con**2 /c1_con**2 *c3_con  &
      + 8.0*c2_con/c1_con*c4_con+24.0*c2_con**4 /c1_con**3   &
      - c1_con*c1_con  &
      *( m_mean*(-240.0*eta**2 +1920.0*eta+3360.0)/zms**8   &
      + (1.0 - m_mean)  &
      *(-120.0*eta**6 -1440.0*eta**5 +14400.0*eta**4   &
      -48000.0*eta**3 +79200.0*eta**2  -66240.0*eta+22560.0)  &
      /(zms*(2.0-eta))**6  )
  
  zdsp  = - 2.0*PI*rho*edI1dz*order1  &
      - PI*rho*order2*m_mean*(c2_con*I2*eta + c1_con*edI2dz)
  zdspdz= zdsp/eta - 2.0*PI*rho*edI1d2*order1  &
      - PI*rho*order2*m_mean*(c3_con*I2*eta  &
      + 2.0*c2_con*edI2dz + c1_con*edI2d2)
  zdspd2= -2.0*zdsp/eta/eta +2.0*zdspdz/eta  &
      - 2.0*PI*rho*edI1d3*order1 - PI*rho*order2*m_mean*(c4_con*I2*eta  &
      + 3.0*c3_con*edI2dz +3.0*c2_con*edI2d2 +c1_con*edI2d3)
  zdspd3= 6.0*zdsp/eta**3  -6.0*zdspdz/eta/eta  &
      + 3.0*zdspd2/eta - 2.0*PI*rho*edI1d4*order1  &
      - PI*rho*order2*m_mean*(c5_con*I2*eta  &
      + 4.0*c4_con*edI2dz +6.0*c3_con*edI2d2  &
      + 4.0*c2_con*edI2d3 + c1_con*edI2d4)
  
  
! ----------------------------------------------------------------------
! p : SAFT (Chen & Kreglewski) dispersion contribution
! ----------------------------------------------------------------------
ELSE
  
  fdspdz = 0.0
  fdspd2 = 0.0
  DO  n = 1,4
    DO m = 1,9
      fdspdz = fdspdz + m_mean/tau * dnm(n,m) * (um/t)**REAL(n) *REAL(m)*(eta/tau)**REAL(m-1)
    END DO
    DO m= 2,9
      fdspd2= fdspd2 + m_mean/tau * dnm(n,m)*(um/t)**REAL(n) *REAL(m)*REAL(m-1)  &
          * (eta/tau)**REAL(m-2) * 1.0/tau
    END DO
  END DO
  zdsp = eta * fdspdz
  zdspdz = (2.0*fdspdz + eta*fdspd2) - zdsp/z3
  
END IF
! --- end of dispersion contribution -----------------------------------


! ----------------------------------------------------------------------
! p: TPT-1-association accord. to Chapman et al.
! ----------------------------------------------------------------------
zhb   = 0.0
zhbdz = 0.0
zhbd2 = 0.0
zhbd3 = 0.0
assoc = .false.
DO i = 1,ncomp
  IF (nhb_typ(i) /= 0) assoc = .true.
END DO
IF (assoc) THEN
  
  DO j = 1, ncomp
    DO i = 1, nhb_typ(j)
      DO k = 1, ncomp
        DO l = 1, nhb_typ(k)
          delta(j,k,i,l) = gij(j,k)    * ass_d(j,k,i,l)
          dq_dz(j,k,i,l) = dgijdz(j,k) * ass_d(j,k,i,l)
          dq_d2(j,k,i,l) = dgijd2(j,k) * ass_d(j,k,i,l)
          dq_d3(j,k,i,l) = dgijd3(j,k) * ass_d(j,k,i,l)
          dq_d4(j,k,i,l) = dgijd4(j,k) * ass_d(j,k,i,l)
        END DO
      END DO
    END DO
  END DO
  
! --- constants for iteration ------------------------------------------
  attenu = 0.7
  tol = 1.d-10
  IF ( eta < 0.2  ) tol = 1.d-12
  IF ( eta < 0.01 ) tol = 1.d-13
  IF ( eta < 1.E-6) tol = 1.d-15
  max_eval = 1000

! --- initialize mx(i,j) -----------------------------------------------
  DO i = 1, ncomp
    DO j = 1, nhb_typ(i)
      mx(i,j) = 1.0
      dmx_dz(i,j) = 0.0
      dmx_d2(i,j) = 0.0
      dmx_d3(i,j) = 0.0
      dmx_d4(i,j) = 0.0
    END DO
  END DO
  
! --- iterate over all components and all sites ------------------------
  ass_cnt = 0
  err_sum = tol + 1.0
  DO WHILE ( err_sum > tol .AND. ass_cnt <= max_eval)
  ass_cnt = ass_cnt + 1
  DO i = 1, ncomp
    DO j = 1, nhb_typ(i)
      sum0 = 0.0
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
      sum4 = 0.0
      DO k = 1, ncomp
        DO l = 1, nhb_typ(k)
          sum0 =sum0 +x(k)*nhb_no(k,l)*     mx(k,l)* delta(i,k,j,l)
          sum1 =sum1 +x(k)*nhb_no(k,l)*(    mx(k,l)* dq_dz(i,k,j,l)  &
              +      dmx_dz(k,l)* delta(i,k,j,l))
          sum2 =sum2 +x(k)*nhb_no(k,l)*(    mx(k,l)* dq_d2(i,k,j,l)  &
              + 2.0*dmx_dz(k,l)* dq_dz(i,k,j,l)  &
              +      dmx_d2(k,l)* delta(i,k,j,l))
          sum3 =sum3 +x(k)*nhb_no(k,l)*(    mx(k,l)* dq_d3(i,k,j,l)  &
              + 3.0*dmx_dz(k,l)* dq_d2(i,k,j,l)  &
              + 3.0*dmx_d2(k,l)* dq_dz(i,k,j,l)  &
              +      dmx_d3(k,l)* delta(i,k,j,l))
          sum4 =sum4 + x(k)*nhb_no(k,l)*(   mx(k,l)* dq_d4(i,k,j,l)  &
              + 4.0*dmx_dz(k,l)* dq_d3(i,k,j,l)  &
              + 6.0*dmx_d2(k,l)* dq_d2(i,k,j,l)  &
              + 4.0*dmx_d3(k,l)* dq_dz(i,k,j,l)  &
              +      dmx_d4(k,l)* delta(i,k,j,l))
        END DO
      END DO
      mx_itr(i,j)= 1.0 / (1.0 + sum0 * rho)
      ndmxdz(i,j)= -(mx_itr(i,j)*mx_itr(i,j))* (sum0/z3t +sum1*rho)
      ndmxd2(i,j)= + 2.0/mx_itr(i,j)*ndmxdz(i,j)*ndmxdz(i,j)  &
          - (mx_itr(i,j)*mx_itr(i,j)) * (2.0/z3t*sum1 + rho*sum2)
      ndmxd3(i,j)= - 6.0/mx_itr(i,j)**2 *ndmxdz(i,j)**3   &
          + 6.0/mx_itr(i,j)*ndmxdz(i,j)*ndmxd2(i,j) - mx_itr(i,j)*mx_itr(i,j)  &
          * (3.0/z3t*sum2 + rho*sum3)
      ndmxd4(i,j)= 24.0/mx_itr(i,j)**3 *ndmxdz(i,j)**4   &
          -36.0/mx_itr(i,j)**2 *ndmxdz(i,j)**2 *ndmxd2(i,j)  &
          +6.0/mx_itr(i,j)*ndmxd2(i,j)**2   &
          +8.0/mx_itr(i,j)*ndmxdz(i,j)*ndmxd3(i,j) - mx_itr(i,j)**2   &
          *(4.0/z3t*sum3 + rho*sum4)
    END DO
  END DO

  err_sum = 0.0
  DO i = 1, ncomp
    DO j = 1, nhb_typ(i)
      err_sum = err_sum + ABS(mx_itr(i,j) - mx(i,j))  &
          + ABS(ndmxdz(i,j) - dmx_dz(i,j)) + ABS(ndmxd2(i,j) - dmx_d2(i,j))
      mx(i,j)     = mx_itr(i,j)*attenu +     mx(i,j) * (1.0-attenu)
      dmx_dz(i,j) = ndmxdz(i,j)*attenu + dmx_dz(i,j) * (1.0-attenu)
      dmx_d2(i,j) = ndmxd2(i,j)*attenu + dmx_d2(i,j) * (1.0-attenu)
      dmx_d3(i,j) = ndmxd3(i,j)*attenu + dmx_d3(i,j) * (1.0-attenu)
      dmx_d4(i,j) = ndmxd4(i,j)*attenu + dmx_d4(i,j) * (1.0-attenu)
    END DO
  END DO
  END DO
  
  IF ( ass_cnt >= max_eval .AND. err_sum > SQRT(tol) ) THEN
   ! WRITE (*,'(a,2G15.7)') 'P_EOS: Max_eval violated (mx) Err_Sum= ',err_sum,tol
    ! stop
  END IF

  
  ! --- calculate the hydrogen-bonding contribution --------------------
  DO i = 1, ncomp
    sum_d1 = 0.0
    sum_d2 = 0.0
    sum_d3 = 0.0
    sum_d4 = 0.0
    DO j = 1, nhb_typ(i)
      sum_d1= sum_d1 +nhb_no(i,j)* dmx_dz(i,j)*(1.0/mx(i,j)-0.5)
      sum_d2= sum_d2 +nhb_no(i,j)*(dmx_d2(i,j)*(1.0/mx(i,j)-0.5)  &
          -(dmx_dz(i,j)/mx(i,j))**2 )
      sum_d3= sum_d3 +nhb_no(i,j)*(dmx_d3(i,j)*(1.0/mx(i,j)-0.5)  &
          -3.0/mx(i,j)**2 *dmx_dz(i,j)*dmx_d2(i,j) + 2.0*(dmx_dz(i,j)/mx(i,j))**3 )
      sum_d4= sum_d4 +nhb_no(i,j)*(dmx_d4(i,j)*(1.0/mx(i,j)-0.5)  &
          -4.0/mx(i,j)**2 *dmx_dz(i,j)*dmx_d3(i,j)  &
          + 12.0/mx(i,j)**3 *dmx_dz(i,j)**2 *dmx_d2(i,j)  &
          - 3.0/mx(i,j)**2 *dmx_d2(i,j)**2  - 6.0*(dmx_dz(i,j)/mx(i,j))**4 )
    END DO
    zhb   = zhb   + x(i) * eta * sum_d1
    zhbdz = zhbdz + x(i) * eta * sum_d2
    zhbd2 = zhbd2 + x(i) * eta * sum_d3
    zhbd3 = zhbd3 + x(i) * eta * sum_d4
  END DO
  zhbdz = zhbdz + zhb/eta
  zhbd2 = zhbd2 + 2.0/eta*zhbdz-2.0/eta**2 *zhb
  zhbd3 = zhbd3 - 6.0/eta**2 *zhbdz+3.0/eta*zhbd2 + 6.0/eta**3 *zhb
END IF
! --- TPT-1-association accord. to Chapman -----------------------------


! ----------------------------------------------------------------------
! p: polar terms
! ----------------------------------------------------------------------
CALL P_POLAR ( zdd, zddz, zddz2, zddz3, zqq, zqqz, zqqz2, zqqz3, zdq, zdqz, zdqz2, zdqz3 )


! ----------------------------------------------------------------------
! compressibility factor z and total p
! as well as derivatives d(z)/d(eta) and d(p)/d(eta) with unit [Pa]
! ----------------------------------------------------------------------
zges   = 1.0 + zhs + zhc + zdsp + zhb + zdd + zqq + zdq
zgesdz = zhsdz + zhcdz + zdspdz + zhbdz + zddz + zqqz + zdqz
zgesd2 = zhsd2 + zhcd2 + zdspd2 + zhbd2 + zddz2 +zqqz2+zdqz2
zgesd3 = zhsd3 + zhcd3 + zdspd3 + zhbd3 + zddz3 +zqqz3+zdqz3

pges   =   zges  *rho *(kbol*t)/1.d-30
pgesdz = ( zgesdz*rho + zges*rho/z3 )*(kbol*t)/1.d-30
pgesd2 = ( zgesd2*rho + 2.0*zgesdz*rho/z3 )*(kbol*t)/1.d-30
pgesd3 = ( zgesd3*rho + 3.0*zgesd2*rho/z3 )*(kbol*t)/1.d-30

END SUBROUTINE P_EOS




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE PHI_POLAR ( k, z3_rk, fdd_rk, fqq_rk, fdq_rk )
!
 USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: k
 REAL, INTENT(IN)                       :: z3_rk
 REAL, INTENT(OUT)                      :: fdd_rk, fqq_rk, fdq_rk
!
! --- local variables---------------------------------------------------
 INTEGER                                :: dipole
 INTEGER                                :: quadrupole
 INTEGER                                :: dipole_quad
! ----------------------------------------------------------------------

 fdd_rk = 0.0
 fqq_rk = 0.0
 fdq_rk = 0.0

 dipole      = 0
 quadrupole  = 0
 dipole_quad = 0
 IF ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
 IF ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
 IF ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

 ! --------------------------------------------------------------------
 ! dipole-dipole term
 ! --------------------------------------------------------------------
 IF (dipole == 1) THEN

    IF (dd_term == 'GV') CALL PHI_DD_GROSS_VRABEC( k, z3_rk, fdd_rk )
    ! IF (dd_term == 'SF') CALL PHI_DD_SAAGER_FISCHER( k )

    IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

 ENDIF

 ! --------------------------------------------------------------------
 ! quadrupole-quadrupole term
 ! --------------------------------------------------------------------
 IF (quadrupole == 1) THEN

    !IF (qq_term == 'SF') CALL PHI_QQ_SAAGER_FISCHER( k )
    IF (qq_term == 'JG') CALL PHI_QQ_GROSS( k, z3_rk, fqq_rk )

    IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

 ENDIF

 ! --------------------------------------------------------------------
 ! dipole-quadrupole cross term
 ! --------------------------------------------------------------------
 IF (dipole_quad == 1) THEN

    IF (dq_term == 'VG') CALL PHI_DQ_VRABEC_GROSS( k, z3_rk, fdq_rk )

    IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

 ENDIF

END SUBROUTINE PHI_POLAR


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE PHI_DD_GROSS_VRABEC( k, z3_rk, fdd_rk )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: k
 REAL, INTENT(IN)                       :: z3_rk
 REAL, INTENT(IN OUT)                   :: fdd_rk
!
! --- local variables---------------------------------------------------
 INTEGER                                :: i, j, l, m

 REAL                                   :: factor2, factor3, z3
 REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
 REAL                                   :: fdd2, fdd3, fdd2x, fdd3x
 REAL, DIMENSION(nc)                    :: my2dd
 REAL, DIMENSION(nc,nc)                 :: Idd2, Idd4, Idd2x, Idd4x
 REAL, DIMENSION(nc,nc,nc)              :: Idd3, Idd3x
! ----------------------------------------------------------------------


 fdd_rk = 0.0
 z3 = eta
 DO i = 1, ncomp
    IF ( uij(i,i) == 0.0 ) write (*,*) 'Surface Tension Code: PHI_DD_GROSS_VRABEC: do not use dimensionless units'
    IF ( uij(i,i) == 0.0 ) stop 5
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



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE PHI_QQ_GROSS( k, z3_rk, fqq_rk )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: k
 REAL, INTENT(IN)                       :: z3_rk
 REAL, INTENT(IN OUT)                   :: fqq_rk
!
! --- local variables---------------------------------------------------
 INTEGER                                :: i, j, l, m

 REAL                                   :: factor2, factor3, z3
 REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
 REAL                                   :: fqq2, fqq3, fqq2x, fqq3x
 REAL, DIMENSION(nc)                    :: qq2
 REAL, DIMENSION(nc,nc)                 :: Iqq2, Iqq4, Iqq2x, Iqq4x
 REAL, DIMENSION(nc,nc,nc)              :: Iqq3, Iqq3x
! ----------------------------------------------------------------------


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
         *uij(k,k)*qq2(k)*sig_ij(k,k)**5 /t/sig_ij(i,k)**7.0
    eij = (parame(i,3)*parame(k,3))**0.5
    fqq2x =fqq2x +factor2*xijfa_x*(Iqq2(i,k)+eij/t*Iqq4(i,k))
    DO j=1,ncomp
       IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
          xijfa =x(i)*rho*uij(i,i)*qq2(i)*sig_ij(i,i)**5 /t  &
               *x(j)*rho*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7.0
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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE PHI_DQ_VRABEC_GROSS( k, z3_rk, fdq_rk )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: k
 REAL, INTENT(IN)                       :: z3_rk
 REAL, INTENT(IN OUT)                   :: fdq_rk
!
! --- local variables---------------------------------------------------
 INTEGER                                :: i, j, l, m

 REAL                                   :: factor2, factor3, z3
 REAL                                   :: xijfa, xijkfa, xijfa_x, xijkf_x, eij
 REAL                                   :: fdq2, fdq3, fdq2x, fdq3x
 REAL, DIMENSION(nc)                    :: my2dd, myfac, qq2, q_fac
 REAL, DIMENSION(nc,nc)                 :: Idq2, Idq4, Idq2x, Idq4x
 REAL, DIMENSION(nc,nc,nc)              :: Idq3, Idq3x
! ----------------------------------------------------------------------

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



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE P_NUMERICAL
!
 USE EOS_VARIABLES
 IMPLICIT NONE
!
!-----local variables-------------------------------------------------
 REAL                                   :: dzetdv, eta_0, dist, fact
 REAL                                   :: fres1, fres2, fres3, fres4, fres5
 REAL                                   :: df_dr, df_drdr, pideal, dpiddz
 REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4, tfr_5
!-----------------------------------------------------------------------


IF (eta > 0.1) THEN
  fact = 1.0
ELSE IF (eta <= 0.1 .AND. eta > 0.01) THEN
  fact = 10.0
ELSE
  fact = 10.0
END IF
dist = eta*3.d-3 *fact
!      dist = eta*4.d-3 *fact
!*****************************
! fuer MC simulation: neues dist:
!      dist = eta*5.d-3*fact

eta_0  = eta
eta  = eta_0 - 2.0*dist
CALL F_NUMERICAL
fres1  = fres
tfr_1  = tfr
eta  = eta_0 - dist
CALL F_NUMERICAL
fres2  = fres
tfr_2  = tfr
eta  = eta_0 + dist
CALL F_NUMERICAL
fres3  = fres
tfr_3  = tfr
eta  = eta_0 + 2.0*dist
CALL F_NUMERICAL
fres4  = fres
tfr_4  = tfr
eta  = eta_0
CALL F_NUMERICAL
fres5  = fres
tfr_5  = tfr

!---------------------------------------------------------
!      ptfr   = (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)/(12.0*dist)
!     &           *dzetdv*(KBOL*T)/1.E-30
!      ztfr =ptfr /( rho * (KBOL*t) / 1.E-30)
!      ptfrdz = (-tfr_4+16.0*tfr_3-3.d1*tfr_5+16.0*tfr_2-tfr_1)
!     &             /(12.0*(dist**2 ))* dzetdv*(KBOL*T)/1.E-30
!     &         + (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)
!     &            /(12.0*dist) * 2.0 *eta*6.0/PI/D
!     &                               *(KBOL*T)/1.E-30
!      ztfrdz=ptfrdz/( rho*(kbol*T)/1.E-30 ) -  ztfr/eta
!      write (*,*) eta,ztfr,ztfrdz

!      dtfr_dr   = (-tfr_4+8.0*tfr_3-8.0*tfr_2+tfr_1)/(12.0*dist)
!      write (*,*) eta,dtfr_dr
!      stop
!---------------------------------------------------------

df_dr   = (-fres4+8.0*fres3-8.0*fres2+fres1) / (12.0*dist)
df_drdr = (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1)  &
    /(12.0*(dist**2 ))


dzetdv = eta*rho

pges   = (-fres4+8.0*fres3-8.0*fres2+fres1)  &
    /(12.0*dist) *dzetdv*(kbol*t)/1.E-30

dpiddz  = 1.0/z3t*(kbol*t)/1.E-30
pgesdz = (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1)  &
    /(12.0*(dist**2 ))* dzetdv*(kbol*t)/1.E-30  &
    + (-fres4+8.0*fres3-8.0*fres2+fres1) /(12.0*dist) * 2.0 *rho  &
    *(kbol*t)/1.E-30 + dpiddz

pgesd2 = (fres4-2.0*fres3+2.0*fres2-fres1) /(2.0*dist**3 )  &
    * dzetdv*(kbol*t)/1.E-30  &
    + (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1) /(12.0*(dist**2 ))  &
    * 4.0 *rho *(kbol*t)/1.E-30 + (-fres4+8.0*fres3-8.0*fres2+fres1)  &
    /(12.0*dist) * 2.0 /z3t *(kbol*t)/1.E-30
pgesd3 = (fres4-4.0*fres3+6.0*fres5-4.0*fres2+fres1) /(dist**4 )  &
    * dzetdv*(kbol*t)/1.E-30 + (fres4-2.0*fres3+2.0*fres2-fres1)  &
    /(2.0*dist**3 ) * 6.0 *rho *(kbol*t)/1.E-30  &
    + (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1)  &
    /(12.0*dist**2 )* 6.0 /z3t *(kbol*t)/1.E-30

!------------------p ideal------------------------------------
pideal = rho * (kbol*t) / 1.E-30

!------------------p summation, p comes out in Pa ------------
pges   = pideal + pges

END SUBROUTINE P_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE H_numerical
!
 USE PARAMETERS, ONLY: RGAS
 USE EOS_VARIABLES
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL                                   :: dist, fact, rho_0
 REAL                                   :: fres1,fres2,fres3,fres4,fres5
 REAL                                   :: f_1, f_2, f_3, f_4
 REAL                                   :: cv_res, t_tmp, zges
 REAL                                   :: f_dt, f_dtdt, f_dr, f_drdr, f_drdt
!-----------------------------------------------------------------------


CALL PERTURBATION_PARAMETER
rho_0 = eta/z3t


fact = 1.0
dist = t * 100.E-5 * fact

t_tmp  = t

t  = t_tmp - 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_NUMERICAL
fres1  = fres
t  = t_tmp - dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_NUMERICAL
fres2  = fres
t  = t_tmp + dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_NUMERICAL
fres3  = fres
t  = t_tmp + 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_NUMERICAL
fres4  = fres
t  = t_tmp
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL F_NUMERICAL
fres5  = fres
! *(KBOL*T)/1.E-30

zges = (p * 1.E-30)/(kbol*t*rho_0)


f_dt   = (-fres4+8.0*fres3-8.0*fres2+fres1)/(12.0*dist)
f_dtdt = (-fres4+16.0*fres3-3.d1*fres5+16.0*fres2-fres1) /(12.0*(dist**2 ))

s_res  =  (- f_dt -fres/t)*RGAS*t  + RGAS * LOG(zges)
h_res  =  ( - t*f_dt  +  zges-1.0  ) * RGAS*t
cv_res = -(   t*f_dtdt + 2.0*f_dt  ) * RGAS*t



t  = t_tmp - 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_NUMERICAL
f_1  = pges/(eta*rho_0*(KBOL*T)/1.E-30)

t  = t_tmp - dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_NUMERICAL
f_2  = pges/(eta*rho_0*(KBOL*T)/1.E-30)

t  = t_tmp + dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_NUMERICAL
f_3  = pges/(eta*rho_0*(KBOL*T)/1.E-30)

t  = t_tmp + 2.0*dist
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_NUMERICAL
f_4  = pges/(eta*rho_0*(KBOL*T)/1.E-30)

t  = t_tmp
CALL PERTURBATION_PARAMETER
eta = z3t*rho_0
CALL P_NUMERICAL

f_dr   = pges  / (eta*rho_0*(KBOL*T)/1.E-30)
f_drdr = pgesdz/ (eta*rho_0*(KBOL*T)/1.E-30) - f_dr*2.0/eta - 1.0/eta**2 

f_drdt   = ( - f_4 + 8.0*f_3 - 8.0*f_2 + f_1 ) / ( 12.0*dist )

cp_res = cv_res - RGAS + RGAS*( zges + eta*t*f_drdt)**2  /  (1.0 + 2.0*eta*f_dr + eta**2 *f_drdr)
! write (*,*) cv_res,cp_res,eta


END SUBROUTINE H_numerical











!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE F_POLAR ( fdd, fqq, fdq )
!
 USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: fdd, fqq, fdq
!
! --- local variables---------------------------------------------------
 INTEGER                                :: dipole
 INTEGER                                :: quadrupole
 INTEGER                                :: dipole_quad
! ----------------------------------------------------------------------

  fdd = 0.0
  fqq = 0.0
  fdq = 0.0

  dipole      = 0
  quadrupole  = 0
  dipole_quad = 0
  IF ( SUM( parame(1:ncomp,6) ) /= 0.0 ) dipole = 1
  IF ( SUM( parame(1:ncomp,7) ) /= 0.0 ) quadrupole = 1
  IF ( dipole == 1 .AND. quadrupole == 1 ) dipole_quad = 1

  ! --------------------------------------------------------------------
  ! dipole-dipole term
  ! --------------------------------------------------------------------
  IF (dipole == 1) THEN

     IF (dd_term == 'GV') CALL F_DD_GROSS_VRABEC( fdd )
     ! IF (dd_term == 'SF') CALL F_DD_SAAGER_FISCHER( k )
     IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

  ENDIF

  ! --------------------------------------------------------------------
  ! quadrupole-quadrupole term
  ! --------------------------------------------------------------------
  IF (quadrupole == 1) THEN

     !IF (qq_term == 'SF') CALL F_QQ_SAAGER_FISCHER( k )
     IF (qq_term == 'JG') CALL F_QQ_GROSS( fqq )
     IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

  ENDIF

  ! --------------------------------------------------------------------
  ! dipole-quadrupole cross term
  ! --------------------------------------------------------------------
  IF (dipole_quad == 1) THEN

     IF (dq_term == 'VG') CALL F_DQ_VRABEC_GROSS( fdq )
     IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

  ENDIF

END SUBROUTINE F_POLAR



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PRESSURE_SPINODAL
!
! iterates the density until the derivative of pressure 'pges' to
! density is equal to zero. A Newton-scheme is used for determining
! the root to the objective function. 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE PRESSURE_SPINODAL
!
 USE BASIC_VARIABLES, ONLY: num
 USE EOS_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, max_i
 REAL                                   :: eta_iteration
 REAL                                   :: error, acc_i, delta_eta
! ----------------------------------------------------------------------

acc_i = 1.d-6
max_i = 30

i = 0
eta_iteration = eta_start

! ----------------------------------------------------------------------
! iterate density until p_calc = p
! ----------------------------------------------------------------------

error = acc_i + 1.0
DO WHILE ( ABS(error) > acc_i .AND. i < max_i )

  i = i + 1
  eta = eta_iteration

  IF ( num == 0 ) THEN
     CALL P_EOS 
  ELSE IF ( num == 1 ) THEN
    CALL P_NUMERICAL
  ELSE IF ( num == 2 ) THEN
  WRITE(*,*) 'CRITICAL RENORM NOT INCLUDED YET'
    !!CALL F_EOS_RN
  ELSE
    write (*,*) 'define calculation option (num)'
  END IF

  error = pgesdz

  delta_eta = error/ pgesd2
  IF ( delta_eta >  0.02 ) delta_eta = 0.02
  IF ( delta_eta < -0.02 ) delta_eta = -0.02

  eta_iteration   = eta_iteration - delta_eta
  ! write (*,'(a,i3,3G18.10)') 'iter',i, error, eta_iteration, pgesdz

  IF (eta_iteration > 0.9)  eta_iteration = 0.5
  IF (eta_iteration <= 0.0) eta_iteration = 1.E-16

END DO

eta = eta_iteration

END SUBROUTINE PRESSURE_SPINODAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_IDEAL_GAS ( fid )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, t, x, rho, PI, KBOL, NAv
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   ::  fid
!---------------------------------------------------------------------
 INTEGER                                :: i
 REAL, DIMENSION(nc)                    :: rhoi
!----------------------------------------------------------------------

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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_HARD_SPHERE ( m_mean2, fhs )
!
 USE EOS_VARIABLES, ONLY: z0t, z1t, z2t, z3t, eta, rho
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN)                       ::  m_mean2
 REAL, INTENT(IN OUT)                   ::  fhs
!---------------------------------------------------------------------
 REAL                                   :: z0, z1, z2, z3, zms
!----------------------------------------------------------------------

 rho = eta / z3t
 z0 = z0t * rho
 z1 = z1t * rho
 z2 = z2t * rho
 z3 = z3t * rho
 zms = 1.0 - z3

 fhs= m_mean2*( 3.0*z1*z2/zms + z2**3 /z3/zms/zms + (z2**3 /z3/z3-z0)*LOG(zms) )/z0


 END SUBROUTINE F_HARD_SPHERE

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_CHAIN_TPT1 ( fhc )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, z0t, z1t, z2t, z3t,  &
                          rho, eta, dij_ab, gij
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   ::  fhc
!---------------------------------------------------------------------
 INTEGER                                :: i, j
 REAL                                   :: z0, z1, z2, z3, zms
!---------------------------------------------------------------------

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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_CHAIN_TPT_D ( fhc )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, z0t, z1t, z2t, z3t, rho, eta,  &
                          dhs, mseg, dij_ab, gij
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(OUT)                      ::  fhc
!---------------------------------------------------------------------
 INTEGER                                :: i, j
 REAL, DIMENSION(nc)                    :: gij_hd
 REAL                                   :: z0, z1, z2, z3, zms
!---------------------------------------------------------------------

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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_CHAIN_HU_LIU ( fhc )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, mseg, x, rho, eta
 IMPLICIT NONE
!
! This subroutine calculates the hard chain contribution of the TPT-Liu-Hu Eos.
!---------------------------------------------------------------------
 REAL, INTENT(OUT)                      ::  fhc
!---------------------------------------------------------------------
 REAL                                   :: a2, b2, c2, a3, b3, c3
 REAL                                   :: a20, b20, c20, a30, b30, c30
 REAL                                   :: sum1, sum2, am, bm, cm
 REAL                                   :: zms
!---------------------------------------------------------------------

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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_CHAIN_HU_LIU_RC ( fhs, fhc )
!
 USE EOS_VARIABLES, ONLY: mseg, chiR, eta
 IMPLICIT NONE
!
! This subroutine calculates the hard chain contribution of the TPT-Liu-Hu Eos.
!---------------------------------------------------------------------
 REAL, INTENT(IN)                       ::  fhs
 REAL, INTENT(OUT)                      ::  fhc
!---------------------------------------------------------------------
 REAL                                   :: a2, b2, c2, a3, b3, c3
 REAL                                   :: para1,para2,para3,para4
 REAL                                   :: aLH,bLH,cLH
!---------------------------------------------------------------------

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




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_DISP_PCSAFT ( fdsp )
!
 USE EOS_VARIABLES, ONLY: PI, rho, eta, z0t, apar, bpar, order1, order2
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdsp
!---------------------------------------------------------------------
 INTEGER                                :: m
 REAL                                   :: I1, I2, c1_con, z3, zms, m_mean
!---------------------------------------------------------------------

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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_DISP_CKSAFT ( fdsp )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, PI, TAU, t, rho, eta, x, z0t, mseg, vij, uij, parame, um
 USE EOS_CONSTANTS, ONLY: DNM
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdsp
!---------------------------------------------------------------------
 INTEGER                                :: i, j, n, m
 REAL                                   :: zmr, nmr, m_mean
!---------------------------------------------------------------------

  m_mean = z0t / ( PI / 6.0 )

  DO i = 1, ncomp
    DO j = 1, ncomp
      vij(i,j)=(0.5*((parame(i,2)*(1.0-0.12 *EXP(-3.0*parame(i,3)/t))**3 )**(1.0/3.0)  &
                    +(parame(j,2)*(1.0-0.12 *EXP(-3.0*parame(j,3)/t))**3 )**(1.0/3.0)))**3
    END DO
  END DO
  zmr = 0.0
  nmr = 0.0
  DO i = 1, ncomp
    DO j = 1, ncomp
      zmr = zmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)*uij(i,j)
      nmr = nmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)
    END DO
  END DO
  um = zmr / nmr
  fdsp  = 0.0
  DO n = 1, 4
    DO m = 1, 9
      fdsp =  fdsp + DNM(n,m) * (um/t)**n *(eta/TAU)**m
    END DO
  END DO
  fdsp = m_mean * fdsp
  

 END SUBROUTINE F_DISP_CKSAFT

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_ASSOCIATION ( eps_kij, k_kij, fhb )
!
 USE EOS_VARIABLES, ONLY: nc, nsite, ncomp, t, z0t, z1t, z2t, z3t, rho, eta, x,  &
                          parame, sig_ij, dij_ab, gij, nhb_typ, mx, nhb_no
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN)                       :: eps_kij, k_kij
 REAL, INTENT(IN OUT)                   :: fhb
!---------------------------------------------------------------------
 LOGICAL                                :: assoc
 INTEGER                                :: i, j, k, l, no, ass_cnt, max_eval
 REAL, DIMENSION(nc,nc)                 :: kap_hb
 REAL, DIMENSION(nc,nc,nsite,nsite)     :: eps_hb
 REAL, DIMENSION(nc,nsite,nc,nsite)     :: delta
 REAL, DIMENSION(nc,nsite)              :: mx_itr
 REAL                                   :: err_sum, sum0, amix, tol, ass_s1, ass_s2
 REAL                                   :: z0, z1, z2, z3, zms
!---------------------------------------------------------------------

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
  
!-----setting the self-association to zero for ionic compounds------
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
          delta(i,k,j,l)=gij(i,j)*kap_hb(i,j)*(EXP(eps_hb(i,j,k,l)/t)-1.0) *sig_ij(i,j)**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              IF ((i+j).EQ.3) delta(i,k,j,l)=94.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
      END DO
      IF ( mx(i,k) == 0.0 ) mx(i,k) = 1.0
    END DO
  END DO
  
!------constants for Iteration ---------------------------------------
  amix = 0.7
  tol = 1.E-10
  IF (eta < 0.2) tol = 1.E-11
  IF (eta < 0.01) tol = 1.E-14
  max_eval = 200
  
! --- Iterate over all components and all sites ------------------------
  ass_cnt = 0
  iterate_TPT1: DO

    ass_cnt = ass_cnt + 1

    DO i = 1, ncomp
      DO k = 1, nhb_typ(i)
        sum0 = 0.0
        DO j = 1, ncomp
          DO l = 1, nhb_typ(j)
            sum0 = sum0 +  x(j)* mx(j,l)*nhb_no(j,l) *delta(i,k,j,l)
          END DO
        END DO
        mx_itr(i,k) = 1.0 / (1.0 + sum0*rho)
      END DO
    END DO

    err_sum = 0.0
    DO i = 1, ncomp
      DO k = 1, nhb_typ(i)
        err_sum = err_sum + ABS( mx_itr(i,k) - mx(i,k) )    ! / ABS(mx_itr(i,k))
        mx(i,k) = mx_itr(i,k) * amix + mx(i,k) * (1.0 - amix)
        IF ( mx(i,k) <= 0.0 ) mx(i,k)=1.E-50
        IF ( mx(i,k) > 1.0 )  mx(i,k)=1.0
      END DO
    END DO

    IF ( err_sum <= tol .OR. ass_cnt >= max_eval ) THEN
      IF ( ass_cnt >= max_eval ) WRITE(*,*) 'F_NUMERICAL: Max_eval violated = ',err_sum,tol
      EXIT iterate_TPT1
    END IF
    
  END DO iterate_TPT1

  DO i = 1, ncomp
    ass_s1 = 0.0
    ass_s2 = 0.0
    DO k = 1, nhb_typ(i)
      ass_s1 = ass_s1 + nhb_no(i,k) * ( 1.0 - mx(i,k) )
      ass_s2 = ass_s2 + nhb_no(i,k) * LOG(mx(i,k))
    END DO
    fhb = fhb + x(i) * ( ass_s2 + ass_s1/2.0 )
  END DO

 END IF

 END SUBROUTINE F_ASSOCIATION


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_ION_DIPOLE_TBH ( fhend )
!
 USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, t, rho, eta, x, z0t, parame, uij, sig_ij
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fhend
!---------------------------------------------------------------------
 INTEGER                                :: i, dipole, ions
 REAL                                   :: m_mean
 REAL                                   :: fh32, fh2, fh52, fh3
 REAL                                   :: e_elem, eps_cc0, rho_sol, dielec
 REAL                                   :: polabil, ydd, kappa, x_dipol, x_ions
 REAL, DIMENSION(nc)                    :: my2dd, z_ii, e_cd, x_dd, x_ii
 REAL                                   :: sig_c, sig_d, sig_cd, r_s
 REAL                                   :: I0cc, I1cc, I2cc, Icd, Idd
 REAL                                   :: Iccc, Iccd, Icdd, Iddd
!---------------------------------------------------------------------

m_mean = z0t / ( PI / 6.0 )

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
         + 2.78E1*(t/293.15))*rho_sol**2   &
         + (-9.63E1/(t/293.15)+4.18E1*(t/293.15)  &
         - 1.02E1*(t/293.15)**2 )*rho_sol**3  +(-4.52E1/(t/293.15)**2   &
         + 8.46E1/(t/293.15)-3.59E1)*rho_sol**4 


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dielec = 1.0

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!----------------Dipole-Ion Term-----------------------------------
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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_ION_ION_PrimMSA ( fcc )
!
 USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, t, rho, x, parame, mx
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fcc
!---------------------------------------------------------------------
 INTEGER                                :: i, j, cc_it, ions
 REAL                                   :: e_elem, eps_cc0, rho_sol, dielec
 REAL                                   :: x_ions
 REAL                                   :: cc_sig1, cc_sig2, cc_sig3
 REAL, DIMENSION(nc)                    :: z_ii, x_ii, sigm_i, my2dd
 REAL                                   :: alpha_2, kappa, ii_par
 REAL                                   :: cc_omeg, p_n, q2_i, cc_q2, cc_gam
 REAL                                   :: cc_error(2), cc_delt
 REAL                                   :: rhs, lambda, lam_s
!---------------------------------------------------------------------

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
  ! ii_par = DSQRT( e_elem**2 /eps_cc0/dielec/kbol/T  &
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
    WRITE (*,*) 'Surface Tension Code, Phase equilibrium calculation: cc error'
    STOP 5
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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_ION_ION_nonPrimMSA ( fdd, fqq, fdq, fcc )
!
 USE EOS_VARIABLES, ONLY: nc, ncomp, t, eta, x, parame, mseg
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdd, fqq, fdq, fcc
!---------------------------------------------------------------------
 INTEGER                                :: dipole
 !REAL                                   :: A_MSA !, A_CC, A_CD, A_DD, U_MSA, chempot
 REAL, DIMENSION(nc)                    :: x_export, msegm
!---------------------------------------------------------------------

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
    write (*,*) 'Surface Tension Code, Phase equilibrium calculation: why are individual contrib. A_CC,A_CD,A_DD not used'
    stop 5
 END IF

 END SUBROUTINE F_ION_ION_nonPrimMSA


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_LC_MayerSaupe ( flc )
!
 USE EOS_VARIABLES, ONLY: nc, PI, KBOL, NAv, ncomp, phas, t, rho, eta,  &
                          x, mseg, parame, E_lc, S_lc, dhs
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: flc
!---------------------------------------------------------------------
 INTEGER                                :: i, j, k
 INTEGER                                :: liq_crystal, count_lc, steps_lc
 REAL                                   :: alpha_lc, tolerance, deltay
 REAL                                   :: integrand1, integrand2, accel_lc
 REAL                                   :: error_lc, u_term, sphase
 REAL, DIMENSION(nc)                    :: z_lc, S_lc1, S_lc2, sumu
 REAL, DIMENSION(nc,nc)                 :: u_lc, klc
!---------------------------------------------------------------------
INTEGER :: stabil
COMMON /stabil / stabil
!---------------------------------------------------------------------


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

 1 CONTINUE

   DO i = 1, ncomp
      IF (S_lc2(i) <= 0.3) S_lc1(i) = S_lc2(i)
      IF (S_lc2(i) > 0.3) S_lc1(i) = S_lc1(i) + (S_lc2(i)-S_lc1(i))*accel_lc
   END DO

   count_lc = count_lc + 1

   ! --- single-particle orientation partition function Z_LC in liquid crystals

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

   ! --- order parameter S_lc in liquid crystals -----------------------

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

   ! --- the anisotropic contribution to the Helmholtz energy ----------

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
   ! pause

 END IF
 ! write (*,'(i2,i2,4(f15.8))') phas,stabil,flc,eta,S_lc2(1),x(1)


 END SUBROUTINE F_LC_MayerSaupe



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE P_POLAR ( zdd, zddz, zddz2, zddz3, zqq, zqqz, zqqz2, zqqz3, zdq, zdqz, zdqz2, zdqz3 )
!
 USE EOS_VARIABLES, ONLY: ncomp, parame, dd_term, qq_term, dq_term
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: zdd, zddz, zddz2, zddz3
 REAL, INTENT(OUT)                      :: zqq, zqqz, zqqz2, zqqz3
 REAL, INTENT(OUT)                      :: zdq, zdqz, zdqz2, zdqz3
!
! --- local variables---------------------------------------------------
 INTEGER                                :: dipole
 INTEGER                                :: quadrupole
 INTEGER                                :: dipole_quad
! ----------------------------------------------------------------------

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

 ! --------------------------------------------------------------------
 ! dipole-dipole term
 ! --------------------------------------------------------------------
 IF (dipole == 1) THEN

    IF (dd_term == 'GV') CALL P_DD_GROSS_VRABEC( zdd, zddz, zddz2, zddz3 )
    ! IF (dd_term == 'SF') CALL F_DD_SAAGER_FISCHER( k )
    IF (dd_term /= 'GV' .AND. dd_term /= 'SF') write (*,*) 'specify dipole term !'

 ENDIF

 ! --------------------------------------------------------------------
 ! quadrupole-quadrupole term
 ! --------------------------------------------------------------------
 IF (quadrupole == 1) THEN

    !IF (qq_term == 'SF') CALL F_QQ_SAAGER_FISCHER( k )
    IF (qq_term == 'JG') CALL P_QQ_GROSS( zqq, zqqz, zqqz2, zqqz3 )
    IF (qq_term /= 'JG' .AND. qq_term /= 'SF') write (*,*) 'specify quadrupole term !'

 ENDIF

 ! --------------------------------------------------------------------
 ! dipole-quadrupole cross term
 ! --------------------------------------------------------------------
 IF (dipole_quad == 1) THEN

    IF (dq_term == 'VG') CALL P_DQ_VRABEC_GROSS( zdq, zdqz, zdqz2, zdqz3 )
    IF (dq_term /= 'VG' ) write (*,*) 'specify DQ-cross term !'

 ENDIF

END SUBROUTINE P_POLAR


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE P_DD_GROSS_VRABEC( zdd, zddz, zddz2, zddz3 )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: zdd, zddz, zddz2, zddz3
! ----------------------------------------------------------------------
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
! ----------------------------------------------------------------------


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



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE P_QQ_GROSS( zqq, zqqz, zqqz2, zqqz3 )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: zqq, zqqz, zqqz2, zqqz3
! ----------------------------------------------------------------------
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
! ----------------------------------------------------------------------

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
               *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7.0
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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE P_DQ_VRABEC_GROSS( zdq, zdqz, zdqz2, zdqz3 )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: zdq, zdqz, zdqz2, zdqz3
! ----------------------------------------------------------------------
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
! ----------------------------------------------------------------------

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



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE F_pert_theory ( fdsp )
!
 USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, p, rho, eta,  &
                          x, z0t, mseg, parame, order1, order2
 USE EOS_NUMERICAL_DERIVATIVES, ONLY: disp_term
 USE DFT_MODULE
 IMPLICIT NONE
!
!---------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdsp
!---------------------------------------------------------------------
 REAL                                   :: I1, I2
 REAL                                   :: z3, zms, c1_con, m_mean
!---------------------------------------------------------------------
 
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
    I1 = - ( - 8.0/9.0 - 4.0/9.0*(rc**(-9) -3.0*rc**(-3) ) - tau_cut/3.0*(rc**3 -1.0) )
    fdsp  = + ( - 2.0*PI*rho*I1*order1 )
    write (*,*) 'caution: not thoroughly checked and tested'

 ELSE
    write (*,*) 'Surface Tension Code, Phase equilibrium calculation: define the type of perturbation theory'
    stop 5
 END IF

 ! I1 = I1 + 4.0/9.0*(2.5**-9 -3.0*2.5**-3 )
 ! fdsp  = + ( - 2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2 )

 END SUBROUTINE F_pert_theory




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE f_pert_theory_mix ( fdsp )
!
 USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, rho, eta, x, parame, mseg, dhs, sig_ij, uij
 USE DFT_MODULE
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdsp
!
! ----------------------------------------------------------------------
 INTEGER                                :: k, ih
 INTEGER                                :: l, m
 REAL                                   :: z3
 REAL                                   :: ua, ua_c, rm
 REAL, DIMENSION(nc,nc)                 :: I1
 REAL                                   :: int10, int11
 REAL                                   :: d_ij, dzr_local
 REAL                                   :: rad, xg, rdf
 REAL                                   :: dg_dz3, dg_dr
 ! REAL                                 :: intgrid(0:5000),intgri2(0:5000), utri(5000),I1_spline
! ----------------------------------------------------------------------

! -----constants--------------------------------------------------------
 ua_c = 4.0 * ( rc**(-12) - rc**(-6) )
 rm = 2.0**(1.0/6.0)

 I1(:,:) = 0.0

 DO l = 1, ncomp
    DO m = 1, ncomp

       rad = rc

       int10 = rc * rc * ua_c
       ! intgrid(0)= int10

       k = 0
       ih = 85

       DO WHILE ( rad /= 1.0 )

          dzr_local = dzr
          IF ( rad - dzr_local <= 1.0 ) dzr_local = rad - 1.0

          rad = rad - dzr_local

          k = k + 1

          d_ij = 0.5*(dhs(l)+dhs(m)) / sig_ij(l,m)   ! dimensionless effective hs-diameter d(T)/sig
          xg = rad / d_ij
          z3 = eta
          rdf  = 1.0
          dg_dz3 = 0.0
          IF ( rad <= rg ) THEN
             IF ( l == 1 .AND. m == 1 ) CALL BI_CUB_SPLINE (z3,xg,ya_11,x1a_11,x2a_11,y1a_11,y2a_11,y12a_11,  &
                  c_bicub_11,rdf,dg_dz3,dg_dr,den_step,ih,k)
             IF ( l /= m ) CALL BI_CUB_SPLINE (z3,xg,ya_12,x1a_12,x2a_12,y1a_12,y2a_12,y12a_12,  &
                  c_bicub_12,rdf,dg_dz3,dg_dr,den_step,ih,k)
             IF ( l == 2 .AND. m == 2 ) CALL BI_CUB_SPLINE (z3,xg,ya_22,x1a_22,x2a_22,y1a_22,y2a_22,y12a_22,  &
                  c_bicub_22,rdf,dg_dz3,dg_dr,den_step,ih,k)
          END IF

          ua = 4.0 * ( rad**(-12) - rad**(-6) )

          int11 = rdf * rad * rad * ua
          I1(l,m) = I1(l,m) + dzr_local * ( int11 + int10 ) / 2.0

          int10 = int11
          ! intgrid(k)= int11

       END DO

       ! stepno = k
       ! CALL SPLINE_PARA (dzr,intgrid,utri,stepno)
       ! CALL SPLINE_INT  (I1_spline,dzr,intgrid,utri,stepno)


       ! caution: 1st order integral is in F_EOS.f defined with negative sign
       ! ---------------------------------------------------------------
       ! cut-off corrections
       ! ---------------------------------------------------------------
       ! I1(l,m) = I1(l,m) + ( 4.0/9.0 * rc**-9 - 4.0/3.0 * rc**-3 )
       ! I2(l,m) = I2(l,m) + 16.0/21.0 * rc**-21 - 32.0/15.0 * rc**-15 + 16.0/9.0 * rc**-9

    END DO
 END DO


 fdsp = 0.0
 DO l = 1, ncomp
    DO m = 1, ncomp
       fdsp = fdsp + 2.0*PI*rho*x(l)*x(m)* mseg(l)*mseg(m)*sig_ij(l,m)**3 * uij(l,m)/t *I1(l,m)
                 ! ( 2.0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2 )
    END DO
 END DO


!!$ IF (disp_term == 'PT1') THEN
!!$    c1_con = 0.0
!!$    I2 = 0.0
!!$ ELSEIF (disp_term == 'PT2') THEN
!!$    zms = 1.0 - z3
!!$    c1_con = 1.0/ (  1.0 + m_mean*(8.0*z3-2.0*z3**2 )/zms**4  &
!!$             + (1.0 - m_mean)*( 20.0*z3 -27.0*z3**2 +12.0*z3**3 -2.0*z3**4 )  &
!!$               /(zms*(2.0-z3))**2 )
!!$ ELSE
!!$    write (*,*) 'define the type of perturbation theory'
!!$    stop
!!$ END IF


END SUBROUTINE f_pert_theory_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE mu_pert_theory_mix ( mu_dsp )
!
 USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, rho, eta, x, parame, mseg, dhs, sig_ij, uij
 USE DFT_MODULE
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: mu_dsp(nc)
!
! ----------------------------------------------------------------------
 INTEGER                                :: k, ih
 INTEGER                                :: l, m
 REAL                                   :: z3
 REAL                                   :: ua, ua_c, rm
 REAL, DIMENSION(nc,nc)                 :: I1, I2
 REAL                                   :: int1_0, int1_1, int2_0, int2_1
 REAL                                   :: d_ij, dzr_local
 REAL                                   :: rad, xg, rdf
 REAL                                   :: dg_dz3, dg_dr
 REAL                                   :: term1(nc), term2
 ! REAL                                 :: intgrid(0:5000),intgri2(0:5000), utri(5000),I1_spline
! ----------------------------------------------------------------------

! -----constants--------------------------------------------------------
 ua_c = 4.0 * ( rc**(-12) - rc**(-6) )
 rm = 2.0**(1.0/6.0)

 I1(:,:) = 0.0
 I2(:,:) = 0.0

 DO l = 1, ncomp

    term1(l) = 0.0

    DO m = 1, ncomp

       rad = rc

       int1_0 = rc * rc * ua_c
       int2_0 = 0.0

       k = 0
       ih = 85

       DO WHILE ( rad /= 1.0 )

          dzr_local = dzr
          IF ( rad - dzr_local <= 1.0 ) dzr_local = rad - 1.0

          rad = rad - dzr_local
          k = k + 1

          d_ij = 0.5*(dhs(l)+dhs(m)) / sig_ij(l,m)   ! dimensionless effective hs-diameter d(T)/sig
          xg = rad / d_ij
          z3 = eta
          rdf  = 1.0
          dg_dz3 = 0.0
          IF ( rad <= rg ) THEN
             IF ( l == 1 .AND. m == 1 ) CALL BI_CUB_SPLINE (z3,xg,ya_11,x1a_11,x2a_11,y1a_11,y2a_11,y12a_11,  &
                  c_bicub_11,rdf,dg_dz3,dg_dr,den_step,ih,k)
             IF ( l /= m ) CALL BI_CUB_SPLINE (z3,xg,ya_12,x1a_12,x2a_12,y1a_12,y2a_12,y12a_12,  &
                  c_bicub_12,rdf,dg_dz3,dg_dr,den_step,ih,k)
             IF ( l == 2 .AND. m == 2 ) CALL BI_CUB_SPLINE (z3,xg,ya_22,x1a_22,x2a_22,y1a_22,y2a_22,y12a_22,  &
                  c_bicub_22,rdf,dg_dz3,dg_dr,den_step,ih,k)
          END IF

          ua = 4.0 * ( rad**(-12) - rad**(-6) )

          int1_1 = rdf * rad * rad * ua
          int2_1 = dg_dz3 * rad * rad * ua
          I1(l,m) = I1(l,m) + dzr_local * ( int1_1 + int1_0 ) / 2.0
          I2(l,m) = I2(l,m) + dzr_local * ( int2_1 + int2_0 ) / 2.0

          int1_0 = int1_1
          int2_0 = int2_1

          term1(l) = term1(l) +4.0*PI*rho*x(m)* mseg(l)*mseg(m) *sig_ij(l,m)**3 *uij(l,m)/t* dzr_local*(int1_1+int1_0)/2.0

       END DO

    END DO
 END DO


 ! DO l = 1, ncomp
 !    term1(l) = 0.0
 !    DO m = 1, ncomp
 !       term1(l) = term1(l) + 4.0*PI*rho*x(m)* mseg(l)*mseg(m) * sig_ij(l,m)**3 * uij(l,m)/t *I1(l,m)
 !    END DO
 ! END DO

 term2 = 0.0
 DO l = 1, ncomp
    DO m = 1, ncomp
       term2 = term2 + 2.0*PI*rho*x(l) * rho*x(m)* mseg(l)*mseg(m) * sig_ij(l,m)**3 * uij(l,m)/t *I2(l,m)
    END DO
 END DO

 DO l = 1, ncomp
    mu_dsp(l) = term1(l) + term2 * PI/ 6.0 * mseg(l)*dhs(l)**3
 END DO

END SUBROUTINE mu_pert_theory_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE F_DD_GROSS_VRABEC( fdd )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: ddp2, ddp3, ddp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdd
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k, m
 INTEGER                                :: ddit, ddmax
 REAL                                   :: factor2, factor3
 REAL                                   :: xijfa, xijkfa, xijf_j, xijkf_j, eij
 REAL                                   :: fdd2, fdd3
 REAL, DIMENSION(nc)                    :: my2dd, my0, alph_tst, z1dd, z2dd, dderror
 REAL, DIMENSION(nc)                    :: fdd2m, fdd3m, fdd2m2, fdd3m2, fddm, fddm2
 REAL, DIMENSION(nc,nc)                 :: Idd2, Idd4
 REAL, DIMENSION(nc,nc,nc)              :: Idd3
! ----------------------------------------------------------------------

 fdd    = 0.0
 ddit   = 0
 ddmax  = 0   ! value assigned, if polarizable compound is present
 fddm(:) = 0.0
 DO i = 1, ncomp
    IF ( uij(i,i) == 0.0 ) write (*,*) 'Surface Tension Code: F_DD_GROSS_VRABEC: do not use dimensionless units'
    IF ( uij(i,i) == 0.0 ) stop 5
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



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE F_QQ_GROSS( fqq )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: qqp2, qqp3, qqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fqq
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k, m
 REAL                                   :: factor2, factor3
 REAL                                   :: xijfa, xijkfa, eij
 REAL                                   :: fqq2, fqq3
 REAL, DIMENSION(nc)                    :: qq2
 REAL, DIMENSION(nc,nc)                 :: Iqq2, Iqq4
 REAL, DIMENSION(nc,nc,nc)              :: Iqq3
! ----------------------------------------------------------------------

 
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
               *x(j)*uij(j,j)*qq2(j)*sig_ij(j,j)**5 /t/sig_ij(i,j)**7.0
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

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE F_DQ_VRABEC_GROSS( fdq )
!
 USE PARAMETERS, ONLY: PI, KBOL
 USE EOS_VARIABLES, ONLY: nc, ncomp, uij, parame, mseg, sig_ij, rho, eta, x, t
 USE EOS_CONSTANTS, ONLY: dqp2, dqp3, dqp4
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: fdq
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k, m
 REAL                                   :: factor2, factor3
 REAL                                   :: xijfa, xijkfa, eij
 REAL                                   :: fdq2, fdq3
 REAL, DIMENSION(nc)                    :: my2dd, myfac, qq2, q_fac
 REAL, DIMENSION(nc,nc)                 :: Idq2, Idq4
 REAL, DIMENSION(nc,nc,nc)              :: Idq3
! ----------------------------------------------------------------------

 
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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE f_dft ( I1_dft, I2_dft )
!
 USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, rho, eta, x, mseg, parame
 USE DFT_MODULE
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: I1_dft
 REAL, INTENT(OUT)                      :: I2_dft
!
! ----------------------------------------------------------------------
 INTEGER                                :: k,ih
 ! REAL                                   :: z3
 REAL                                   :: ua, ua_c, ua_2, ua_c_2, rm
 REAL                                   :: int10, int11, int20, int21
 REAL                                   :: dg_drho
 REAL                                   :: rad, xg, rdf, rho_st, msegm
 REAL                                   :: sig_ij
 REAL                                   :: dg_dr, dzr_org !,rdf_d
 ! REAL                                 :: intgrid(0:NDFT),intgri2(0:NDFT)
! ----------------------------------------------------------------------

! -----constants--------------------------------------------------------
msegm = parame(1,1)
rho_st = rho * parame(1,2)**3

ua_c = 4.0 * ( rc**(-12) - rc**(-6) )
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
   ua = 4.0 * ( rad**(-12) - rad**(-6) )
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
I1_dft= - I1_dft - ( 4.0/9.0 * rc**(-9) - 4.0/3.0 * rc**(-3) )

! CALL SPLINE_PARA (dzr,intgri2,utri,stepno)
! CALL SPLINE_INT (I2,dzr,intgri2,utri,stepno)

I2_dft = I2_dft + 16.0/21.0 * rc**(-21) - 32.0/15.0 * rc**(-15) + 16.0/9.0 * rc**(-9)


END SUBROUTINE f_dft




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 REAL FUNCTION TANGENT_VALUE2 ( optpara, n )
! SUBROUTINE TANGENT_VALUE ( fmin, optpara, n )
!
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN)                       :: optpara(n)
 !REAL, INTENT(IN)                       :: optpara(:)
 !REAL, INTENT(IN OUT)                   :: fmin
!
! ----------------------------------------------------------------------
 INTEGER                                :: i
 REAL                                   :: lnphi(np,nc),ph_frac, gibbs_full(np),xlnx1,xlnx2
 REAL, DIMENSION(nc)                    :: ni_1, ni_2
! ----------------------------------------------------------------------


 ! --- setting of mole fractions ---------------------------------------
 DO i = 1, ncomp
   IF ( optpara(i) < -300.0 ) THEN
     ni_2(i) = 0.0
   ELSE
     ni_2(i) = EXP( optpara(i) )
   END IF
 END DO

 DO i = 1, ncomp
   ni_1(i) = xif(i) - ni_2(i)
   IF ( ni_2(i) > xif(i) ) THEN
      ni_2(i) = xif(i)
      ni_1(i) = xif(i) * 1.E-20
   ENDIF
 END DO

 xi(2,1:ncomp) = ni_2(1:ncomp) / SUM( ni_2(1:ncomp) )
 lnx(2,1:ncomp) = optpara(1:ncomp) - LOG( SUM( ni_2(1:ncomp) ) )

 ph_frac = SUM( ni_1(1:ncomp) )
 xi(1,1:ncomp) = ni_1(1:ncomp) / ph_frac
 lnx(1,1:ncomp) = LOG( ni_1(1:ncomp) ) - LOG( ph_frac )
 ! write (*,'(a,4G18.8)') 'FF',(xif(i),i=1,ncomp)
 ! write (*,'(a,4G18.8)') 'AA',(xi(1,i),i=1,ncomp)
 ! write (*,'(a,3G18.8)') 'BB',(xi(2,i),i=1,ncomp)

 CALL fugacity (lnphi)
 !CALL enthalpy_etc

 gibbs(1) = SUM( xi(1,1:ncomp) * lnphi(1,1:ncomp) )   ! dimensionless g/RT
 gibbs(2) = SUM( xi(2,1:ncomp) * lnphi(2,1:ncomp) )

 xlnx1 = SUM( xi(1,1:ncomp)*lnx(1,1:ncomp) )          ! dimensionless s/RT
 xlnx2 = SUM( xi(2,1:ncomp)*lnx(2,1:ncomp) )

 gibbs_full(1) = gibbs(1) + xlnx1
 gibbs_full(2) = gibbs(2) + xlnx2

 TANGENT_VALUE2 = gibbs_full(1)*ph_frac + gibbs_full(2)*(1.0-ph_frac)
 !fmin = gibbs_full(1)*ph_frac + gibbs_full(2)*(1.0-ph_frac)
  !write (*,'(a,4G18.8)') 'TP',TANGENT_VALUE2,(lnx(1,i),i=1,ncomp)
  !write (*,'(a,4G18.8)') 'al',ph_frac,(lnx(2,i), i=1,ncomp)
  !write (*,*) ' '
  !pause

END FUNCTION TANGENT_VALUE2








