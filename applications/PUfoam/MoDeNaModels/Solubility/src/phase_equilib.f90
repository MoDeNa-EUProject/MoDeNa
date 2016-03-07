
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE phase_equilib
!
! This subroutine varies a predefined "running variable" and
! organizes phase equilibrium calculations. For an isotherm
! calculation e.g., the running variable is often the pressure. The
! code is designed to deliver only converged solutions. In order to
! enforce convergence, a step-width adjustment (reduction) is
! implemented.

!   VARIABLE LIST:
! running     defines the running variable. For example: if you want
!             to calculate the vapor pressure curve of a component
!             starting from 100C to 200C, then running is 't'. The
!             temperature is step-wise increased until the end-
!             -temperature of 200C is reached.
!             (in this example end_x=200+273.15)
! end_x       end point for running variable
! steps       No. of calculation steps towards the end point of calc.
! converg     0 if no convergence achieved, 1 if converged solution
!
!   PREREQUISITES:
! prior to execution of this routine, the follwing variables have to
! be defined: "val_init" an array containing the starting values for
! this iteration, "it(i)" provides the information, which variable is
! determined iteratively, "sum_rel(i)" indicates, which mole fraction
! is determined from the summation relation sum(xi)=1. Furthermore,
! the number of phases and the variables provided by the subroutine
! INPUT are required.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE phase_equilib (end_x,steps,converg)

  USE BASIC_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: end_x
  REAL, INTENT(IN)                       :: steps
  INTEGER, INTENT(OUT)                   :: converg

  !-----------------------------------------------------------------------------
  INTEGER                                :: k, count1,count2,runindex,maxiter
  REAL                                   :: delta_x,delta_org,val_org,runvar
  CHARACTER (LEN=2)                      :: compon
  LOGICAL                                :: continue_cycle
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! determine variable that is varied ("running variable")
  !-----------------------------------------------------------------------------

  IF (running(1:2) == 'd1') runindex = 1
  IF (running(1:2) == 'd2') runindex = 2
  IF (running(1:1) == 't')  runindex = 3
  IF (running(1:1) == 'p')  runindex = 4
  IF (running(1:2) == 'x1') compon = running(3:3)
  IF (running(1:2) == 'x1') READ (compon,*) k
  IF (running(1:2) == 'x1') runindex = 4+k
  IF (running(1:2) == 'x2') compon = running(3:3)
  IF (running(1:2) == 'x2') READ (compon,*) k
  IF (running(1:2) == 'x2') runindex = 4+ncomp+k
  IF (running(1:2) == 'l1') compon = running(3:3)
  IF (running(1:2) == 'l1') READ (compon,*) k
  IF (running(1:2) == 'l1') runindex = 4+k
  IF (running(1:2) == 'l2') compon = running(3:3)
  IF (running(1:2) == 'l2') READ (compon,*) k
  IF (running(1:2) == 'l2') runindex = 4+ncomp+k

  !-----------------------------------------------------------------------------
  ! define control variables of the running variable
  !-----------------------------------------------------------------------------

  maxiter = 200
  IF ( ncomp >= 3 ) maxiter = 1000
  count1 = 0
  count2 = 0
  delta_x   = ( end_x - val_init(runindex) ) / steps
  delta_org = ( end_x - val_init(runindex) ) / steps
  val_org = val_init(runindex)
  IF ( running(1:1) == 'x' ) THEN
     delta_x   = ( end_x - EXP(val_init(runindex)) ) / steps
     delta_org = ( end_x - EXP(val_init(runindex)) ) / steps
     val_org = EXP(val_init(runindex))
  END IF

  continue_cycle = .true.

  !-----------------------------------------------------------------------------
  ! stepwise change of the running variables
  !-----------------------------------------------------------------------------

  DO WHILE ( continue_cycle )

     count1 = count1 + 1
     count2 = count2 + 1
     ! val_org = val_init(runindex)


     CALL objective_ctrl (converg)

     IF (converg == 1) THEN
        val_init( 1:(4+ncomp*nphas) ) = val_conv( 1:(4+ncomp*nphas) )
        IF (outp >= 1 .AND. (ABS(delta_x) > 0.1*ABS(delta_org) .OR. count2 == 2)) CALL output
        ! write (17,'(4G18.8)') rhoi_cal(1,1:2)*32000.,rhoi_cal(2,1:2)*32000.
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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE objective_ctrl
!
! This subroutine controls the iso-fugacity iteration. It uses
! the variables defined in the array "val_init". If successfull,
! the converged values are written to "val_conv", and the flag
! converg is set to 1.
! See also above desciption for subroutine PHASE_EQUILIB
! This routine calls SUBROUTINE HYBRID, which is a solver (modified
! POWELL HYBRID METHOD). HYBRID is freely available for non-commercial
! applications. HYBRID requires three definitions:
! 1.the number of equations to be solved (=No. of variables to be
!   iterated). The appropriate parameter is: "n_unkw"
! 2.the equations to be iterated, they are here gathered in the SUB-
!   ROUTINE OBJEC_FCT (see below). Since HYBRID is a root finder,
!   these objective functions are iterated to be zero (essentially,
!   OBJEC_FCT contains the iso-fugacity relation.
! 3.an array of variables is required, containing the quatities to be
!   iterated. This array is termed "y(i)"
!
!   INPUT VARIABLES:
! val_init(i)   array containing (densities,T,P,lnx's) serving as
!               starting values for the phase equilibrium calculation
! it(i)         contains the information, which variable is  deter-
!               mined iteratively. For syntax refer e.g.to SUB BINMIX.
! sum_rel(i)    indicates, which mole fraction is determined from the
!               summation relation sum(xi)=1
!
!   OUTPUT VARIABLES:
! val_conv(i)   array containing the converged system variables
!               analogous to "val_init"
! converg     0 if no convergence achieved, 1 if converged solution
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE objective_ctrl (converg)

  USE BASIC_VARIABLES
  USE STARTING_VALUES
  USE EOS_VARIABLES, only: PI, KBOL, mseg, dhs
  USE Solve_NonLin
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(OUT)                   :: converg

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE objec_fct ( iter_no, y, residu, dummy )
       INTEGER, INTENT(IN)                :: iter_no
       REAL, INTENT(IN)                   :: y(iter_no)
       REAL, INTENT(OUT)                  :: residu(iter_no)
       INTEGER, INTENT(IN OUT)            :: dummy
     END SUBROUTINE objec_fct
  END INTERFACE

  INTEGER                                :: info, k, posn, i
  REAL, ALLOCATABLE                      :: y(:), diag(:), residu(:)
  REAL                                   :: x_init, x_solut, r_diff1, r_diff2, totres
  REAL                                   :: r_thrash, x_thrash
  CHARACTER (LEN=2)                      :: compon

  INTEGER                                :: ph_split, promising_min
  REAL                                   :: eta_trial, phi_2_start
  REAL, dimension(nc)                    :: xif_save, rhoi_feed, rhoi1, rhoi2
  LOGICAL                                :: convergence
  !-----------------------------------------------------------------------------

  info=1

  ALLOCATE( y(n_unkw), diag(n_unkw), residu(n_unkw) )

  IF (num == 0) acc_a  = 1.E-7
  IF (num == 0) step_a = 2.E-8
  IF (num == 1) acc_a  = 1.E-7
  IF (num == 1) step_a = 2.E-8
  IF (num == 2) acc_a  = 5.E-7
  IF (num == 2) step_a = 1.E-7

  !-----------------------------------------------------------------------------
  ! assign starting values for the iterated variables (to vector y)
  !-----------------------------------------------------------------------------

  posn = 0
  DO i = 1,n_unkw
     posn = posn + 1
     IF (it(i) == 't') y(posn) = val_init(3)
     IF (it(i) == 'p') y(posn) = val_init(4)
     IF (it(i) == 'lnp') y(posn) = LOG( val_init(4) )
     IF (it(i) == 'fls') y(posn) = alpha
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') y(posn) = val_init(4+k)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') y(posn) = val_init(4+ncomp+k)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') y(posn) = val_init(4+ncomp+ncomp+k)
  END DO

  !-----------------------------------------------------------------------------
  ! initialize some variables
  !-----------------------------------------------------------------------------

  CALL init_vars

  x_init  = 0.0
  DO i = 1,ncomp
     IF (lnx(1,i) /= 0.0 .AND. lnx(2,i) /= 0.0) THEN
        x_init  = x_init  + ABS( 1.0 - lnx(2,i)/lnx(1,i) )
     ELSE
        x_init  = x_init  + ABS( 1.0 - EXP(lnx(2,i))/EXP(lnx(1,i)) )
     END IF
  END DO

  !-----------------------------------------------------------------------------
  ! solve the phase equilibrium conditions
  !-----------------------------------------------------------------------------

  CALL hbrd (objec_fct, n_unkw, y, residu, step_a, acc_a, info, diag)
  IF ( info <= 1 .AND. SUM( ABS( residu(1:n_unkw) ) ) > 100.0*acc_a )  &
       CALL hbrd (objec_fct, n_unkw, y, residu, step_a, acc_a, info, diag)

  !-----------------------------------------------------------------------------
  ! determine, whether the solution is acceptable (if yes: convergence = .true.)
  !-----------------------------------------------------------------------------

  x_solut = 0.0
  DO i = 1,ncomp
     IF ( lnx(1,i) /= 0.0 .AND. lnx(2,i) /= 0.0 ) THEN
        x_solut = x_solut + ABS( 1.0 - lnx(2,i)/lnx(1,i) )
     ELSE
        IF (lnx(1,i) < 300.0 .AND. lnx(1,i) > -300.0 )  &
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
  IF ( ncomp == 1 .AND. r_diff1 < 1.E-5      ) convergence = .false.
  IF ( totres > 1000.0*acc_a .AND. .NOT. generate_starting_val ) convergence = .false.
  IF ( totres > 10000.0*acc_a .AND. generate_starting_val ) convergence = .false.

  !-----------------------------------------------------------------------------
  ! if convergence: test stability, save result to vector (val_conv) and output result
  !-----------------------------------------------------------------------------

  IF ( convergence ) THEN

     !-----------------------------------------------------------------------------
     ! store result and output result to terminal
     !-----------------------------------------------------------------------------

     converg = 1
     CALL converged
     CALL enthalpy_etc
     ! IF (outp == 2 ) CALL output

     !-----------------------------------------------------------------------------
     ! test phase stability
     !-----------------------------------------------------------------------------

     if ( check_stability_of_phases ) then
        my_f( 1:ncomp ) = my_cal(1,1:ncomp)
        rhoi_feed( 1:ncomp ) = rhoi_cal( MAXLOC(dense(1:2), 1), 1:ncomp )
        xiF_save(1:ncomp) = xiF(1:ncomp)
        xiF(1:ncomp) = 0.95 * rhoi_cal( MAXLOC(dense(1:2), 1), 1:ncomp )  &
                              / sum( rhoi_cal( MAXLOC(dense(1:2), 1), 1:ncomp ) )  &
                     + 0.05 * rhoi_cal( MINLOC(dense(1:2), 1), 1:ncomp )  &
                              / sum( rhoi_cal( MINLOC(dense(1:2), 1), 1:ncomp ) )
        if ( minval( dense(1:2) ) < 0.18 ) eta_trial = 0.40
        if ( minval( dense(1:2) ) >=0.18 ) eta_trial = p / ( KBOL*1.E30 * t )  &
             * PI/6.0 * SUM( rhoi_feed( 1:ncomp ) / sum(rhoi_feed( 1:ncomp )) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
        call phase_stability ( rhoi_feed, eta_trial, ph_split, rhoi_trial )

        if ( ph_split == 1 ) then
           rhoi1( 1:ncomp ) = rhoi_trial( 1:ncomp )
           rhoi2( 1:ncomp ) = rhoi_feed( 1:ncomp )
           call tangent_plane_line_search ( rhoi1, rhoi2, phi_2_start, promising_min )

           dense( MINLOC(dense(1:2), 1) ) = PI/6.0 *sum( rhoi_trial(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
           call rachford_rice ( converg, rhoi1, rhoi2 )
           if ( converg == 1 ) then
              call occupy_val_init ( rhoi1, rhoi2 )
              if ( ncomp == 2 ) CALL select_sum_rel (1,0,1)
              if ( ncomp == 2 ) CALL select_sum_rel (2,0,2)
              val_conv = val_init
              if ( outp > 0 ) write (*,*) ' '
              if ( outp > 0 ) write (*,*) '===== the type (LLE, VLE) of phase equilibrium changed ====='
              if ( outp > 0 ) write (*,*) ' '
           else
              call restore_converged
              converg = 1
              if ( outp > 0 .AND. ph_split == 1 ) write (*,*) 'calc. phase equilibr. not stable (press enter)'
              if ( outp > 0 .AND. ph_split == 1 ) read (*,*)
           end if
        else
           call restore_converged
           converg = 1
        end if
        xiF(1:ncomp) = xiF_save(1:ncomp)
     end if

  ELSE

     converg = 0

  END IF

  DEALLOCATE( y, diag, residu )

END SUBROUTINE objective_ctrl



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE objec_fct
!
! This subroutine contains the equations to be solved numerically
! (iso-fugacity: fi'-fi''=0) as well as other dependent equations,
! which can be solved analytically, namely the summation relation
! xi=1-sum(xj) or the condition of equal charge for electrolyte
! solutions.
! This subroutine is required and controlled by the solver HBRD !
! HBRD varies the variables "y(i)" and eveluates the result of
! these changes from this routine.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE objec_fct ( iter_no, y, residu, dummy )

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: density_error
  use utilities, only: x_summation
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: iter_no
  REAL, INTENT(IN)                       :: y(iter_no)
  REAL, INTENT(OUT)                      :: residu(iter_no)
  INTEGER, INTENT(IN OUT)                :: dummy

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, ph,k,posn, skip,phase
  REAL                                   :: lnphi(np,nc),isofugacity
  CHARACTER (LEN=2)                      :: compon
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! assign values for the iterated quantities
  !-----------------------------------------------------------------------------

  posn = 0
  DO i = 1,n_unkw
     posn = posn + 1
     IF (it(i) == 't') t = y(posn)
     IF (it(i) == 'p') p = y(posn)
     IF (it(i) == 'lnp') p = EXP( y(posn) )
     IF (it(i) == 'fls') alpha = y(posn)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '1') lnx(1,k) = y(posn)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '2') lnx(2,k) = y(posn)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') compon = it(i)(3:3)
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') READ (compon,*) k
     IF (it(i)(1:1) == 'x' .AND. it(i)(2:2) == '3') lnx(3,k) = y(posn)
  END DO

  DO k = 1,ncomp
     IF (lnx(1,k) > 0.0) lnx(1,k) = 0.0
     IF (lnx(2,k) > 0.0) lnx(2,k) = 0.0
  END DO

  !-----------------------------------------------------------------------------
  ! saveguard real values
  !-----------------------------------------------------------------------------

  IF (p < 1.E-100) p = 1.E-12
  IF ( p /= p ) p = 1000.0                ! rebounce for the case of NaN-solver output
  IF ( t /= t ) t = 300.0                 ! rebounce for the case of NaN-solver output
  IF ( alpha /= alpha ) alpha = 0.5       ! rebounce for the case of NaN-solver output

  !-----------------------------------------------------------------------------
  ! setting of mole fractions
  !-----------------------------------------------------------------------------

  DO ph = 1, nphas
     DO i = 1, ncomp
        IF ( lnx(ph,i) < -300.0 ) THEN
           xi(ph,i) = 0.0
        ELSE
           xi(ph,i) = EXP( lnx(ph,i) )
        END IF
     END DO
  END DO

  !-----------------------------------------------------------------------------
  ! use summation relation to determine missing mole fractions
  !-----------------------------------------------------------------------------

  IF (ncomp > 1) CALL x_summation

  !-----------------------------------------------------------------------------
  ! calculate chemical potentials
  !-----------------------------------------------------------------------------

  CALL fugacity (lnphi)

  !-----------------------------------------------------------------------------
  ! determine the residuum to isofugacity relation
  !-----------------------------------------------------------------------------

  phase = 2
  DO i = 1,n_unkw
     skip = 0  !for ions/polymers, the isofug-eq. is not always solved
     IF (n_unkw < (ncomp*(nphas-1))) skip = ncomp*(nphas-1) - n_unkw
     IF ((i+skip-ncomp*(phase-2)) > ncomp) phase = phase + 1
     residu(i) = isofugacity((i+skip-ncomp*(phase-2)),phase,lnphi)
     if ( density_error(phase) /= 0.0 ) residu(i) = residu(i) + SIGN( density_error(phase),  residu(i) ) * 0.001
  END DO

END SUBROUTINE objec_fct



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! REAL FUNCTION isofugacity
!
! calculates the deviation from the condition of equal fugacities in
! logarithmic form.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION isofugacity (i,phase,lnphi)

  USE BASIC_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: phase
  REAL, INTENT(IN)                       :: lnphi(np,nc)

  !-----------------------------------------------------------------------------
  INTEGER                                :: p1, p2
  !-----------------------------------------------------------------------------


  ! p1=1
  p1 = phase-1
  p2 = phase

  isofugacity = scaling(i) *( lnphi(p2,i)+lnx(p2,i)-lnx(p1,i)-lnphi(p1,i) )
  ! write (*,'(a,      4G18.8)') ' t, p ',t,p,dense(1),dense(2)
  ! write (*,'(a,i3,i3,3G24.14)') ' phi_V',i,p2,lnx(p2,i),lnphi(p2,i),dense(p2)
  ! write (*,'(a,i3,i3,3G24.14)') ' phi_L',i,p1,lnx(p1,i),lnphi(p1,i),dense(p1)
  ! write (*,*) ' isofugacity',i,isofugacity, scaling(i)
  ! write (*,'(a,i3,4G18.8)') ' isofugacity',i,isofugacity, lnphi(p2,i)+lnx(p2,i), -lnx(p1,i)-lnphi(p1,i)
  ! read (*,*)

END FUNCTION isofugacity
