module utilities

  implicit none

  PRIVATE
  PUBLIC :: converged, restore_converged, switch, init_vars, SI_DENS, select_sum_rel,  &
            determine_flash_it, determine_flash_it2, new_flash, x_summation, flash_alpha,  &
            flash_sum, neutr_charge, molefrac, file_open, p_calc, dens_calc, enthalpy_etc,  &
            CRITICAL, PURE_CRIT_POINT_ITERATION, fden_calc, ONLY_ONE_TERM_EOS_NUMERICAL,  &
            RESTORE_PREVIOUS_EOS_NUMERICAL, paus, error_message

contains

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine converged
!
! Once a converged solution for an equilibrium calculation is found,
! this subroutine writes variables to an array "val_conv".
! (= short for converged values)
! The array comprises the following elements
! element 0              density of third phase
! element 1              density of first phase
! element 2              density of second phase
! element 3              temperature [K]
! element 4              pressure [Pa]
! element 5,..(4+NCOMP)  mole-fraction of comp. i in log. from (phase 1)
! element ...            mole-fraction of comp. i in log. from (phase 2)
! element ...            mole-fraction of comp. i in log. from (phase 3)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine converged

  use BASIC_VARIABLES

  integer                                :: i, ph
  !-----------------------------------------------------------------------------

  val_conv(0)  = dense(3)
  val_conv(1)  = dense(1)
  val_conv(2)  = dense(2)
  val_conv(3)  = t
  val_conv(4)  = p
  DO ph = 1, nphas
     DO i = 1, ncomp
        val_conv(4+i+(ph-1)*ncomp) = lnx(ph,i)
     END DO
  END DO

end subroutine converged



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine init_vars
!
! This subroutine writes variables from an array "val_init" to the
! system-variables. Those variables
! include some specifications but also some starting values for a
! phase equilibrium calculation. (val_init stands for 'initial value')

! The array comprises the following elements
! element 0              density of third phase
! element 1              density of first phase
! element 2              density of second phase
! element 3              temperature [K]
! element 4              pressure [Pa]
! element 5,..(5+NCOMP)  mole-fraction of comp. i in log. from (phase 1)
! element ...            mole-fraction of comp. i in log. from (phase 2)
! element ...            mole-fraction of comp. i in log. from (phase 3)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine init_vars

  use BASIC_VARIABLES

  integer                                :: i, ph
  !-----------------------------------------------------------------------------

  densta(3) = val_init(0)
  densta(1) = val_init(1)
  densta(2) = val_init(2)
  t = val_init(3)
  p = val_init(4)
  DO ph = 1, nphas
     DO i = 1, ncomp
        lnx(ph,i) = val_init(4+i+(ph-1)*ncomp)
     END DO
  END DO

end subroutine init_vars


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine switch
!
! This subroutine switches the phase-index of phase 1 and phase 2.
! Sometimes this is handy.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine switch

  use BASIC_VARIABLES

  integer                                :: i
  real                                   :: savval(nc*np+6)
  !-----------------------------------------------------------------------------

  savval(1) = val_init(1)
  savval(2) = val_init(2)
  val_init(1) = savval(2)
  val_init(2) = savval(1)
  DO i = 1, ncomp
     savval(i) = val_init(4+i)
     val_init(4+i) = val_init(4+i+ncomp)
  END DO
  DO i = 1, ncomp
     val_init(4+i+ncomp) = savval(i)
  END DO

end subroutine switch



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine restore_converged

  use BASIC_VARIABLES

  integer                                :: i, ph
  !-----------------------------------------------------------------------------

  dense(3) = val_conv(0)
  dense(1) = val_conv(1)
  dense(2) = val_conv(2)
  t = val_conv(3)
  p = val_conv(4)
  DO ph = 1, nphas
     DO i = 1, ncomp
        lnx(ph,i) = val_conv(4+i+(ph-1)*ncomp)
        xi(ph,i) = 0.0
        if ( lnx(ph,i) > -300.0 ) xi(ph,i) = EXP( lnx(ph,i) )
     END DO
  END DO

end subroutine restore_converged



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine SI_DENS (density,w)
!
! This subroutine calculates the (macroskopic) fluid-density in
! units [kg/m3] from the dimensionless density (eta=zeta3).
! Further, mass fractions are calculated from mole fractions.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine SI_DENS (density,w)

  use PARAMETERS, ONLY: pi, nav, tau
  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(OUT)                      :: density(np)
  real, INTENT(OUT)                      :: w(np,nc)

  !-----------------------------------------------------------------------------
  integer                                :: i, ph
  real                                   :: mm_mean, rho, z3t
  real                                   :: dhs(nc), d00(nc), t_p, pcon, l_st
  !-----------------------------------------------------------------------------


  DO i = 1, ncomp
     IF (eos == 1) THEN
        dhs(i) = parame(i,2) * ( 1.0 - 0.12 *EXP( -3.0*parame(i,3)/t ) )
     ELSE IF (eos == 0) THEN
        d00(i) = ( 1.E30/1.E6*tau*parame(i,2)*6.0/pi/nav )**(1.0/3.0)
        dhs(i) = d00(i) * ( 1.0 - 0.12 *EXP( -3.0*parame(i,3)/t ) )
     ELSE IF (eos == 4) THEN
        dhs(i) = (  0.1617/0.3107 / ( 0.689+0.311*(t/parame(i,3)/1.328)**0.3674 )  &
             / ( pi/6.0 )  )**(1.0/3.0) * parame(i,2)
     ELSE IF (eos == 5.OR.eos == 6) THEN
        l_st = parame(1,25)
        IF (ncomp /= 1) write (*,*) ' ERROR for EOS = 5'
        t_p  =((34.037352+17.733741*l_st) /(1.0+0.53237307*l_st+12.860239*l_st**2 ))**0.5
        IF (l_st == 0.0) t_p = t_p/4.0
        IF (eos == 5 .AND. l_st /= 0.0) t_p = t_p/4.0*parame(1,1)**2
        t_p = t/parame(i,3)/t_p
        pcon =0.5256+3.2088804*l_st**2 -3.1499114*l_st**2.5 +0.43049357*l_st**4
        dhs(i) = (  pcon/(0.67793+0.32207*(t_p)**0.3674) /(pi/6.0)  )**(1.0/3.0) *parame(i,2)
     ELSE IF (eos == 8) THEN
        dhs(i) = parame(i,2)*(1.0+0.2977*t/parame(i,3))  &
             /(1.0+0.33163*t/parame(i,3) +1.0477E-3*(t/parame(i,3))**2 )
     ELSE
        write (*,*) 'define EOS in subroutine: SI_DENS'
        stop 5
     END IF
  END DO

  DO ph = 1, nphas
     mm_mean = 0.0
     z3t = 0.0
     DO i = 1, ncomp
        mm_mean = mm_mean + xi(ph,i)*mm(i)
        z3t = z3t + xi(ph,i) * parame(i,1) * dhs(i)**3
     END DO
     z3t = pi/6.0 * z3t
     rho = dense(ph)/z3t
     density(ph) = rho * mm_mean * 1.E27 /nav
     DO i = 1, ncomp
        w(ph,i) = xi(ph,i) * mm(i)/mm_mean
     END DO
     !  write (*,*) density(ph),rho,mm_mean,1.d27 /NAV
  END DO

end subroutine SI_DENS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine select_sum_rel
!
! This subroutine determines which component of a phase "ph" is calculated
! from the summation relation x_i = 1 - sum(x_j). The other components are,
! by default, said to be iterated during the phase equilibrium calculation.
!
! Note that for flash calculations not all of these mole fractions are in
! fact iterated - this is raken care of in "determine_flash_it".
!
! ph           phase
! excl         exclude comp. n
! startindex   assign it(startindex) for quantities to be iterated
!              (further it(startindex+1) is assigned, for a ternary
!              mixture, etc.)
!
! sum_index    indicates the component, with the largest mole
!              fraction. If ph=1 and sum_index=2, we define
!              sum_rel(ph=1)='x12', so that this component is
!              calculated from the summation relation.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine select_sum_rel ( ph, excl, startindex )

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer, INTENT(IN)                    :: ph
  integer, INTENT(IN)                    :: excl
  integer, INTENT(IN)                    :: startindex
  !-----------------------------------------------------------------------------
  integer                                :: i,j
  integer                                :: sum_index = 0
  real                                   :: xmax(np)
  ! character                              :: compNo*2,phasNo*2
  !-----------------------------------------------------------------------------

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
        ! write (*,*) ph,i,xi(ph,i),sum_rel(ph)
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
        ! write (*,*) 'iter  ',it(startindex+j)
        j = j + 1
     END IF

  END DO

end subroutine select_sum_rel


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine determine_flash_it

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer                                :: i, j, ph
  integer                                :: min_index = 0
  integer                                :: max_index = 0
  integer                                :: phase = 0
  integer                                :: phase_from_flash
  real                                   :: d_xmax, xmin, xmax
  !-----------------------------------------------------------------------------

  write (*,*) 'entering determine_flash_it, press return'
  read (*,*)

  IF (nphas /= 2) WRITE (*,*) 'DETERMINE_FLASH_IT: only 2 phases!'
  IF (nphas /= 2) STOP 2

  !-----------------------------------------------------------------------------
  ! determine component (max_index) with largest difference in xi of both phases
  !-----------------------------------------------------------------------------

  d_xmax = 0.0
  DO i = 1, ncomp
     IF ( ABS( xi(1,i)-xi(2,i) ) > d_xmax ) THEN
        d_xmax = ABS( xi(1,i)-xi(2,i) )
        max_index = i
     END IF
  END DO

  !-----------------------------------------------------------------------------
  ! determine the phase that is determined from component balance
  !-----------------------------------------------------------------------------

  IF ( MINVAL( lnx(1,1:ncomp) ) < MINVAL( lnx(2,1:ncomp) ) ) THEN
     !IF ( ABS(xi(1,max_index)-xif(max_index)) > ABS(xif(max_index)-xi(2,max_index)) ) THEN
     phase_from_flash = 2
     IF (max_index == 1) sum_rel(2) = 'fl1'
     IF (max_index == 2) sum_rel(2) = 'fl2'
     IF (max_index == 3) sum_rel(2) = 'fl3'
     IF (max_index == 4) sum_rel(2) = 'fl4'
     IF (max_index == 5) sum_rel(2) = 'fl5'
     IF (max_index == 6) sum_rel(2) = 'fl6'
     IF (max_index > 6) WRITE (*,*) 'DETERMINE_FLASH_IT:extend list!'
  ELSE
     phase_from_flash = 1
     IF (max_index == 1) sum_rel(1) = 'fl1'
     IF (max_index == 2) sum_rel(1) = 'fl2'
     IF (max_index == 3) sum_rel(1) = 'fl3'
     IF (max_index == 4) sum_rel(1) = 'fl4'
     IF (max_index == 5) sum_rel(1) = 'fl5'
     IF (max_index == 6) sum_rel(1) = 'fl6'
     IF (max_index > 6) WRITE (*,*) 'DETERMINE_FLASH_IT:extend list!'
  END IF

  DO ph = 1, nphas

     IF (ph == phase_from_flash) THEN

        !-----------------------------------------------------------------------------
        ! find substance of lowest conc. This substance is iterated (in log)
        !-----------------------------------------------------------------------------
        xmin = 10.0
        DO i = 1, ncomp
           IF ( xi(ph,i) < xmin ) THEN
              xmin = xi(ph,i)
              min_index = i
           END IF
        END DO
        IF (phase_from_flash == 1) it(1) = 'x1'       ! define phase of it(1)
        IF (phase_from_flash == 2) it(1) = 'x2'
        IF (min_index == 1) it(1) = trim(it(1))//'1'  ! define substance of it(1)
        IF (min_index == 2) it(1) = trim(it(1))//'2'
        IF (min_index == 3) it(1) = trim(it(1))//'3'
        IF (min_index == 4) it(1) = trim(it(1))//'4'
        IF (min_index == 5) it(1) = trim(it(1))//'5'
        IF (min_index == 6) it(1) = trim(it(1))//'6'

     ELSE             ! for the other phase ...

        !-----------------------------------------------------------------------------
        ! find substance of highest conc. ==> calc. from summation relation
        !-----------------------------------------------------------------------------
        xmax=0.0
        DO i = 1, ncomp
           IF (xi(ph,i) > xmax) THEN
              xmax = xi(ph,i)
              max_index = i
           END IF
        END DO
        IF (phase_from_flash == 1) phase = 2
        IF (phase_from_flash == 1) sum_rel(2) = 'x2'  ! define phase of sum_rel(...)
        IF (phase_from_flash == 2) phase = 1
        IF (phase_from_flash == 2) sum_rel(1) = 'x1'

        IF (max_index == 1) sum_rel(phase) = trim(sum_rel(phase))//'1' ! define substance of sum_rel(...)
        IF (max_index == 2) sum_rel(phase) = trim(sum_rel(phase))//'2'
        IF (max_index == 3) sum_rel(phase) = trim(sum_rel(phase))//'3'
        IF (max_index == 4) sum_rel(phase) = trim(sum_rel(phase))//'4'
        IF (max_index == 5) sum_rel(phase) = trim(sum_rel(phase))//'5'
        IF (max_index == 6) sum_rel(phase) = trim(sum_rel(phase))//'6'

        j = 0
        DO i = 1, ncomp
           IF ( i /= max_index ) THEN
              IF (phase == 1) it(2+j) = 'x1'
              IF (phase == 2) it(2+j) = 'x2'
              IF (i == 1) it(2+j) = trim(it(2+j))//'1'
              IF (i == 2) it(2+j) = trim(it(2+j))//'2'
              IF (i == 3) it(2+j) = trim(it(2+j))//'3'
              IF (i == 4) it(2+j) = trim(it(2+j))//'4'
              IF (i == 5) it(2+j) = trim(it(2+j))//'5'
              IF (i == 6) it(2+j) = trim(it(2+j))//'6'
              j = j + 1
           END IF
        END DO

     END IF
  END DO

end subroutine determine_flash_it



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine determine_flash_it2

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer                                :: i, k, ph
  real                                   :: n_phase1, n_phase2, max_x_diff
  !-----------------------------------------------------------------------------

  n_phase1 = 0.0

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
  DO i = 1, ncomp
     IF ( ABS( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) ) > max_x_diff ) THEN
        max_x_diff = ABS( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) )
        n_phase1 = ( xif(i) - EXP( lnx(2,i) ) ) / ( EXP( lnx(1,i) ) - EXP( lnx(2,i) ) )
        n_phase2 = 1.0 - n_phase1
     END IF
  END DO
  if ( n_phase1 < 1.E-12 ) n_phase1 = 1.E-12
  if ( n_phase1 > ( 1.0 - 1.E-12 ) ) n_phase1 = 1.0 - 1.E-12
  n_phase2 = 1.0 - n_phase1
  lnx(1,1:ncomp) = lnx(1,1:ncomp) + LOG( n_phase1 )   ! these x's are treated as mole numbers
  lnx(2,1:ncomp) = lnx(2,1:ncomp) + LOG( n_phase2 )   ! these x's are treated as mole numbers

  val_init(1) = dense(1)
  val_init(2) = dense(2)
  val_init(3) = t
  val_init(4) = p
  DO ph = 1, nphas
     DO k = 1, ncomp
        val_init(4+k+(ph-1)*ncomp) = lnx(ph,k) ! - LOG( SUM( EXP( lnx(ph,1:ncomp) ) ) )
        !            write (*,*) ph,k, lnx(ph,k)
     END DO
  END DO

end subroutine determine_flash_it2



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine new_flash (ph_it)

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer, INTENT(IN)                    :: ph_it

  !-----------------------------------------------------------------------------
  integer                                :: i, ph_cal
  real, DIMENSION(nc)                    :: ni_1, ni_2
  !-----------------------------------------------------------------------------

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
     IF ( ni_2(i) > xif(i) .OR. SUM( ni_1(1:ncomp)) <= 0.0 ) THEN
        ni_2(i) = xif(i)
        ni_1(i) = xif(i) * 1.E-20
     ENDIF
  END DO

  !  if ( SUM( ni_2(1:ncomp) ) > (1.0 - 1.E-12) ) then
  !     ni_2(1:ncomp) = ni_2(1:ncomp) * ( 1.0 - 1.E-12 ) / SUM( ni_2(1:ncomp) )
  !  end if
  !  if ( SUM( ni_1(1:ncomp) ) > (1.0 - 1.E-12) ) then
  !     ni_1(1:ncomp) = ni_1(1:ncomp) * ( 1.0 - 1.E-12 ) / SUM( ni_1(1:ncomp) )
  !  end if

  xi(ph_it,1:ncomp) = ni_2(1:ncomp) / SUM( ni_2(1:ncomp) )
  DO i = 1, ncomp
     IF ( xi(ph_it,i) >= 1.E-300 ) lnx(ph_it,i) = LOG( xi(ph_it,i) )
     IF ( xi(ph_it,i) <  1.E-300 ) lnx(ph_it,i) = 1.E-300
  END DO
  xi(ph_cal,1:ncomp) = ni_1(1:ncomp) / SUM( ni_1(1:ncomp) )
  DO i = 1, ncomp
     IF ( xi(ph_cal,i) >= 1.E-300 ) lnx(ph_cal,i) = LOG( xi(ph_cal,i) )
     IF ( xi(ph_cal,i) <  1.E-300 ) lnx(ph_cal,i) = 1.E-300
  END DO

end subroutine new_flash



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine x_summation
!
! This subroutine solves the summation relation: xi=1-sum(xj)
! The variable "sum_rel(i)" contains the information, which mole
! fraction is the one to be calculated here. Consider the example
! sum_rel(1)='x12'. The fist letter 'x' of this variable indicates,
! that this subroutine needs to be executed and that the mole
! fraction of a component has to be calculated. The second letter
! of the string points to phase 1, the third letter to component 2.
!     If the fist letter is 'e', not 'x', then the subroutine
! NEUTR_CHARGE is called. This is the case of electrolyte solutions,
! neutral charges have to be enforced in all phases (see below).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine x_summation

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer                                :: i, j, comp_i, ph_i
  real                                   :: sum_x
  character (LEN=2)                      :: phasno
  character (LEN=2)                      :: compno
  logical                                :: flashcase2
  !-----------------------------------------------------------------------------

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
        READ (phasno,*) ph_i
        compno = sum_rel(j)(3:3)
        READ (compno,*) comp_i
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
        WRITE (*,*) 'summation relation not defined'
        STOP 3
     END IF

  END DO

  IF ( it(1) == 'fls' ) CALL flash_sum
  IF ( flashcase2 ) CALL flash_alpha

end subroutine x_summation


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine flash_alpha
!
! This subroutine calculates all molefractions of one phase
! from a component balance. What is needed for this calculation
! are all molefractions of the other phase (nphas=2, so far)
! and the phase fraction alpha.
! Alpha is calculated from the mole fraction
! of component {sum_rel(j)(3:3)}. If for example sum_rel(2)='fl3',
! then the alpha is determined from the molefraction of comp. 3 and
! all molefractions of one phase are calculated using that alpha-value.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine flash_alpha

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer                                :: i, j, comp_i, phase1, phase2
  character (LEN=2)                      :: compno
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! first calculate the phase fraction alpha from a known composition
  ! of component sum_rel(j)(3:3).
  !-----------------------------------------------------------------------------

  DO j = 1, nphas
     IF ( sum_rel(j)(1:2) == 'fl' ) THEN
        compno = sum_rel(j)(3:3)
        READ (compno,*) comp_i
        IF ( (xi(1,comp_i)-xi(2,comp_i)) /= 0.0 ) THEN
           alpha = (xif(comp_i)-xi(2,comp_i)) / (xi(1,comp_i)-xi(2,comp_i))
           write (*,*) 'flsh',(xif(comp_i)-xi(2,comp_i)),(xi(1,comp_i)-xi(2,comp_i))
        ELSE
           alpha = 0.5
           WRITE (*,*) 'FLASH_ALPHA:error in calc. of phase fraction',comp_i
        END IF
        ! IF (alpha <= 0.0 .OR. alpha >= 1.0) WRITE(*,*) 'FLASH_ALPHA: error',alpha
        IF (alpha > 1.0) alpha = 1.0
        IF (alpha < 0.0) alpha = 0.0
     END IF
  END DO

  !-----------------------------------------------------------------------------
  ! determine which phase is fully determined by iterated molefractions (+ summation relation)
  !-----------------------------------------------------------------------------
  phase1 = 0
  phase2 = 0
  DO i = 1, ncomp
     IF ( it(i)(2:2) == '1' ) phase1 = phase1 + 1
     IF ( it(i)(2:2) == '2' ) phase2 = phase2 + 1
  END DO
  DO i = 1, ncomp
     IF ( sum_rel(i)(2:2) == '1' ) phase1 = phase1 + 1
     IF ( sum_rel(i)(2:2) == '2' ) phase2 = phase2 + 1
  END DO


  IF ( phase1 == ncomp ) THEN
     ! --------------------------------------------------------------------
     ! phase 1 is defined by iterated molefractions + summation relation
     ! phase 2 is determined from componennt balance (using alpha)
     ! --------------------------------------------------------------------
     IF ( alpha == 1.0 ) alpha = 1.0 - 1.0E-10
     DO i = 1, ncomp
        xi(2,i) = ( xif(i) - alpha*xi(1,i) ) / (1.0-alpha)
        IF ( xi(2,i) < 0.0 ) xi(2,i) = 0.0
        IF ( xi(2,i) > 1.0 ) xi(2,i) = 1.0
        IF ( xi(2,i) /= 0.0 ) THEN
           lnx(2,i) = LOG( xi(2,i) )
        ELSE
           lnx(2,i) = -100000.0
        END IF
        write (*,*) 'fl_cal ph=2',i,lnx(2,i),xi(2,i)
     END DO
  ELSE IF ( phase2 == ncomp ) THEN
     ! --------------------------------------------------------------------
     ! phase 2 is defined by iterated molefractions + summation relation
     ! phase 1 is determined from componennt balance (using alpha)
     ! --------------------------------------------------------------------
     DO i = 1, ncomp
        xi(1,i) = ( xif(i) - (1.0-alpha)*xi(2,i) ) /alpha
        IF ( xi(1,i) < 0.0 ) xi(1,i) = 0.0
        IF ( xi(1,i) > 1.0 ) xi(1,i) = 1.0
        IF ( xi(1,i) /= 0.0 ) THEN
           lnx(1,i) = LOG( xi(1,i) )
        ELSE
           lnx(1,i) = -100000.0
        END IF
        write (*,*) 'fl_cal ph=1',i,lnx(1,i),xi(1,i),alpha
     END DO
  ELSE
     WRITE (*,*) ' FLASH_ALPHA: undefined flash-case'
     STOP 2
  END IF

end subroutine flash_alpha



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine flash_sum

  use BASIC_VARIABLES

  integer                                :: i, j, ph_i, phase1, phase2
  !-----------------------------------------------------------------------------

  phase1=0
  phase2=0
  DO j = 1, ncomp
     IF (it(j)(2:2) == '1') phase1=phase1+1
     IF (it(j)(2:2) == '2') phase2=phase2+1
  END DO

  IF (phase1 == ncomp-1) THEN
     ph_i = 1
  ELSE IF (phase2 == ncomp-1) THEN
     ph_i = 2
  ELSE
     WRITE (*,*) ' FLASH_SUM: undefined flash-case'
     STOP
  END IF



  IF (ph_i == 1) THEN
     DO i = 1, ncomp
        IF (alpha > DMIN1(1.0,xif(i)/xi(1,i),  &
             (xif(i)-1.0)/(xi(1,i)-1.0),alpha)) THEN
           WRITE (*,*) ' FLASH_SUM: exeeded 1st alpha-bound'
           alpha=DMIN1(1.0,xif(i)/xi(1,i),(xif(i)-1.0)/(xi(1,i)-1.0))
        END IF
     END DO
     DO i = 1, ncomp
        xi(2,i) = ( xif(i) - alpha*xi(1,i) ) / (1.0-alpha)
        IF (xi(2,i) > 0.0) THEN
           lnx(2,i) = LOG(xi(2,i))
        ELSE
           xi(2,i) = 0.0
           lnx(2,i) = -100000.0
        END IF
     END DO
  ELSE IF (ph_i == 2) THEN
     DO i = 1, ncomp
        IF (alpha > DMAX1(0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)),  &
             1.0-xif(i)/xi(2,i),alpha)) THEN
           WRITE (*,*) ' FLASH_SUM: exeeded 2nd alpha-bound'
           WRITE (*,*) 0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)), 1.0-xif(i)/xi(2,i)
           alpha=DMAX1(0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)), 1.0-xif(i)/xi(2,i))
        END IF
     END DO
     DO i = 1, ncomp
        xi(1,i) = ( xif(i) - (1.0-alpha)*xi(2,i) ) / alpha
        !    write (*,*) 'x1,i',xi(1,i),xi(2,i),alpha
        IF (xi(1,i) > 0.0) THEN
           lnx(1,i) = LOG(xi(1,i))
        ELSE
           xi(1,i) = 0.0
           lnx(1,i) = -100000.0
        END IF
     END DO
  END IF

end subroutine flash_sum



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine neutr_charge
!
! This subroutine is called for electrolye solutions, where a
! neutral overall-charge has to be enforced in all phases. The basic
! philosophy is similar to the above described routine X_SUMMATION.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine neutr_charge(i)

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  integer, INTENT(IN)                    :: i

  !-----------------------------------------------------------------------------
  integer                                :: comp_e, ph_e
  real                                   :: sum_c
  character (LEN=2)                      :: phasno
  character (LEN=2)                      :: compno
  !-----------------------------------------------------------------------------

  phasno = sum_rel(i)(2:2)
  READ (phasno,*) ph_e
  compno = sum_rel(i)(3:3)
  READ (compno,*) comp_e

  sum_c = 0.0
  write (*,*) 'there must be an error in neutr_charge'
  stop
  ! there is an error in the following passage. The index i is an
  ! argument to this subroutine - I guess it is INTENT(IN), so the
  ! index in the following loop can not be i.
  !
  ! I have commented the loop until I check the code.
  !DO i = 1, ncomp
  !  IF ( comp_e /= i .AND. parame(i,10) /= 0.0)  &
  !      sum_c = sum_c + xi(ph_e,i)*parame(i,10)
  !END DO

  xi(ph_e,comp_e) = - sum_c
  IF (xi(ph_e,comp_e) < 0.0) xi(ph_e,comp_e)=0.0
  IF (xi(ph_e,comp_e) /= 0.0) THEN
     lnx(ph_e,comp_e) = LOG(xi(ph_e,comp_e))
  ELSE
     lnx(ph_e,comp_e) = -100000.0
  END IF

  ! xi(2,1) = xi(2,2)
  ! IF (xi(2,1).NE.0.0) lnx(2,1) = LOG(xi(2,1))

end subroutine neutr_charge



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine molefrac (w,phas1,phas2)
!
! This subroutine calculates mole fractions from mass fractions.
! The calculation is performed for phase 1 if phas1=1 and for
! phase 2 if phas2=1.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine molefrac ( w, phas1, phas2 )

  use BASIC_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(IN)                     :: w(np,nc)
  integer, INTENT(IN)                  :: phas1
  integer, INTENT(IN)                  :: phas2

  !-----------------------------------------------------------------------------
  integer :: i,ph
  real :: mm_denom
  !-----------------------------------------------------------------------------


  DO ph = 1, nphas
     IF ((ph == 1.AND.phas1 == 1).OR.(ph == 2.AND.phas2 == 1)) THEN
        mm_denom = 0.0
        DO i = 1, ncomp
           mm_denom = mm_denom + w(ph,i)/mm(i)
        END DO
        IF (mm_denom > 0.0) THEN
           DO i = 1, ncomp
              xi(ph,i)=w(ph,i)/mm(i)/mm_denom
           END DO
        END IF
     END IF
  END DO

end subroutine molefrac



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine file_open
!
! This subroutine opens files for reading. Beforehand, it checks
! whether this file is available.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine file_open(filename,file_number)

  !-----------------------------------------------------------------------------
  character (LEN=*)                        :: filename
  integer                                  :: file_number
  logical                                  :: filefound
  !-----------------------------------------------------------------------------

  INQUIRE (FILE=filename, EXIST = filefound)
  IF (filefound) THEN
     OPEN (file_number, FILE = filename)
  ELSE
     write (*,*) ' Solubility Code: FOLLOWING FILE CAN NOT BE OPENED', filename
     stop 6
  END IF

end subroutine file_open


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine p_calc (pges_transfer, zges)
!
! This subroutine serves as an iterface to the EOS-routines. The
! system pressure corresponding to given (desity,T,xi) is calculated.
! (Note: the more common interface is subroutine FUGACITY. This
! routine is only used for one-phase systems, e.g. calculation of
! virial coefficients)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine p_calc (pges_transfer, zges)

  use BASIC_VARIABLES
  use EOS_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(IN OUT)                   :: pges_transfer
  real, INTENT(OUT)                      :: zges
  !-----------------------------------------------------------------------------

  IF (nphas /= 1 ) THEN
     write (*,*) 'P_CALC: can only be called for single phases'
     stop 5
  ENDIF

  IF (eos < 2) THEN

     phas = 1
     eta = dense(1)
     x(1:ncomp) = xi(1,1:ncomp)

     CALL PERTURBATION_PARAMETER
     CALL p_eos_interface

     pges_transfer = pges
     rho = eta / z3t
     zges = (pges * 1.E-30) / ( kbol * t * rho )

  ELSE
     write (*,*) ' subroutine P_CALC not available for cubic EOS'
     stop 6
  END IF

end subroutine p_calc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine dens_calc
!
! This subroutine serves as an interface to the EOS-routines. The
! densities corresponding to given (P,T,xi) are calculated.
! (Note: the more common interface is subroutine FUGACITY.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine dens_calc ( rho_phas )

  use BASIC_VARIABLES
  use EOS_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(OUT)                      :: rho_phas(np)

  integer                                :: ph
  !-----------------------------------------------------------------------------


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
        write (*,*) ' subroutine DENS_CALC not available for cubic EOS'
        stop 6
     END IF

  END DO

end subroutine dens_calc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine enthalpy_etc
!
! This subroutine serves as an interface to the EOS-routines. The
! residual enthalpy h_res, residual entropy s_res, residual Gibbs
! enthalpy g_res, and residual heat capacity at constant pressure
! (cp_res) corresponding to converged conditions are calculated.
! The conditions in (T,P,xi,rho) need to be converged equilibrium
! conditions  !!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine enthalpy_etc

  use BASIC_VARIABLES
  use EOS_VARIABLES

  integer                                :: ph
  !-----------------------------------------------------------------------------

  IF ( eos <= 1 ) THEN

     DO ph = 1, nphas

        phas = ph
        eta = dense(ph)
        ! eta_start = dense(ph)
        x(1:ncomp) = xi(ph,1:ncomp)

       CALL h_eos_interface
        enthal(ph) = h_res
        entrop(ph) = s_res
        ! gibbs(ph)  = h_res - t * s_res  ! already defined in eos.f90 (including ideal gas)
        cpres(ph)  = cp_res
        speed_of_sound(ph) = speed_sound

     END DO
     IF (nphas == 2) h_lv = enthal(2)-enthal(1)

  ENDIF

end subroutine enthalpy_etc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine CRITICAL
!
! This subroutine serves as an interface for the calculation of
! a pure component's critical point.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine CRITICAL ( tc, pc, rhoc )

  use BASIC_VARIABLES
  use EOS_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(IN OUT)                   :: tc
  real, INTENT(IN OUT)                   :: pc
  real, INTENT(OUT)                      :: rhoc

  !-----------------------------------------------------------------------------
  integer                                :: phassav
  real                                   :: density(np), w(np,nc)
  !-----------------------------------------------------------------------------

  phassav= nphas
  nphas  = 1
  xi(1,1) = 1.0         ! only for pure components


  IF (eos < 2.AND.ncomp == 1) THEN

     t = tc
     p = pc
     phas = 1
     eta       = 0.2
     eta_start = 0.2
     x(1:ncomp)   = xi(nphas,1:ncomp)

     CALL PERTURBATION_PARAMETER
     CALL PURE_CRIT_POINT_ITERATION
     dense(1)= eta
     CALL SI_DENS (density,w)
     rhoc = density(1)
     tc   = t
     pc   = p

  ELSE
     write (*,*) ' PURE_CRIT_POINT_ITERATION not available for cubic EOS'
     write (*,*) ' and not for mixtures'
     stop 5
  END IF

  nphas = phassav

end subroutine CRITICAL

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine PURE_CRIT_POINT_ITERATION
!
! This subroutine determines the critical point of a pure substance.
! The density and temperature of the pure component is iterated
! until the derivative of (dP/d_rho) and the second derivative
! (dP2/d2_rho) are zero.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine PURE_CRIT_POINT_ITERATION

  use BASIC_VARIABLES, ONLY: num
  use EOS_VARIABLES

  !-----------------------------------------------------------------------------
  integer :: count, count_t, max_eta_count
  real    :: tdelt, pdz1, pdz2, t_tol, eta_tol, zc, r_damp, t_damp, eta_l1
  !-----------------------------------------------------------------------------

  r_damp = 1.0
  IF (num == 2) r_damp = 0.15
  tdelt   = 0.1
  t_tol   = 0.005

  eta_tol= 0.2
  IF (num == 1) eta_tol = 2.0
  IF (num == 2) eta_tol = 0.5

  count_t = 0
  max_eta_count = 200
  t = 2.0 * t            ! prefer a larger starting value for t over a too low value.
  pgesdz = t_tol + 1.0

  DO WHILE (ABS(pgesdz) > t_tol .AND. count_t < 30)
     t = t - tdelt
     count = 0
     count_t =count_t + 1
     t_damp = 1.0
     pgesd2 = eta_tol + 1.0

     DO WHILE ( ABS(pgesd2) > eta_tol .AND. count < max_eta_count )
        count =count + 1
        CALL PERTURBATION_PARAMETER
        CALL p_eos_interface
        eta = eta - r_damp*pgesd2/pgesd3
        IF (eta < 0.0) eta=0.001
        IF (eta > 0.7) eta=0.5
        ! write(*,'(a,i4,4(E16.7))') 'c_1',count,eta,pges,pgesd2,pgesd3
        eta_l1 = eta
     END DO
     pdz1 = pgesdz



     t = t + tdelt
     count = 0
     pgesd2 = eta_tol + 1.0

     DO WHILE ( ABS(pgesd2) > eta_tol .AND. count < max_eta_count )
        count =count +1
        CALL PERTURBATION_PARAMETER
        CALL p_eos_interface
        eta = eta - r_damp*pgesd2/pgesd3
        IF (eta < 0.0) eta = 0.001
        IF (eta > 0.7) eta = 0.5
     END DO
     pdz2 = pgesdz

     IF ( ABS(pdz2 / ( (pdz2-pdz1)/tdelt ))/t >= 0.05 ) t_damp = 0.5
     t = t - t_damp * pdz2  /  ( (pdz2-pdz1)/tdelt )

     IF (.NOT.( t > 0.0 )) t = 300.0

  END DO

  IF (ABS(pgesdz) > t_tol) write (*,*) 'PURE_CRIT_POINT_ITERATION: T-iteration not converged'
  IF (ABS(pgesd2) > eta_tol) write (*,*) 'PURE_CRIT_POINT: density-iteration not converged',pgesd2

  p  = pges
  rho= eta/z3t
  zc = (p * 1.E-30)/(kbol*t*rho)

end subroutine PURE_CRIT_POINT_ITERATION


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine fden_calc (fden, rhoi)

  use BASIC_VARIABLES
  use EOS_VARIABLES

  !-----------------------------------------------------------------------------
  real, INTENT(OUT)                      :: fden
  real, INTENT(IN OUT)                   :: rhoi(nc)
  !-----------------------------------------------------------------------------
  real                                   :: rhot, fden_id
  !-----------------------------------------------------------------------------


  IF (eos < 2) THEN

     rhot = SUM( rhoi(1:ncomp) )
     x(1:ncomp) = rhoi(1:ncomp) / rhot

     CALL PERTURBATION_PARAMETER
     eta = rhot * z3t
     eta_start = eta

     call f_eos_interface

     fden_id = SUM( rhoi(1:ncomp) * ( LOG( rhoi(1:ncomp) ) - 1.0 ) )

     fden = fres * rhot  +  fden_id

  ELSE
     write (*,*) ' subroutine FDEN_CALC not available for cubic EOS'
     stop 5
  END IF

end subroutine fden_calc


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine ONLY_ONE_TERM_EOS_NUMERICAL ( only_term, type_of_term )

  use EOS_NUMERICAL_DERIVATIVES

  character (LEN=9)                        :: only_term, type_of_term
  !-----------------------------------------------------------------------------

  save_eos_terms(1) = ideal_gas
  save_eos_terms(2) = hard_sphere
  save_eos_terms(3) = chain_term
  save_eos_terms(4) = disp_term
  save_eos_terms(5) = hb_term
  save_eos_terms(6) = LC_term
  save_eos_terms(7) = branch_term
  save_eos_terms(8) = II_term
  save_eos_terms(9) = ID_term
  save_eos_terms(10)= dd_term
  save_eos_terms(11)= qq_term
  save_eos_terms(12)= dq_term

  ideal_gas   = 'no'
  hard_sphere = 'no'
  chain_term  = 'no'
  disp_term   = 'no'
  hb_term     = 'no'
  LC_term     = 'no'
  branch_term = 'no'
  II_term     = 'no'
  ID_term     = 'no'
  dd_term     = 'no'
  qq_term     = 'no'
  dq_term     = 'no'

  IF ( only_term == 'ideal_gas' )   ideal_gas   = trim( adjustl( type_of_term ) )
  IF ( only_term == 'hard_sphere' ) hard_sphere = trim( adjustl( type_of_term ) )
  IF ( only_term == 'chain_term' )  chain_term  = trim( adjustl( type_of_term ) )
  IF ( only_term == 'disp_term' )   disp_term   = trim( adjustl( type_of_term ) )
  IF ( only_term == 'hb_term' )     hb_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'LC_term' )     LC_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'branch_term' ) branch_term = trim( adjustl( type_of_term ) )
  IF ( only_term == 'II_term' )     II_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'ID_term' )     ID_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'dd_term' )     dd_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'qq_term' )     qq_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'dq_term' )     dq_term     = trim( adjustl( type_of_term ) )

end subroutine ONLY_ONE_TERM_EOS_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine RESTORE_PREVIOUS_EOS_NUMERICAL

  use EOS_NUMERICAL_DERIVATIVES

  !-----------------------------------------------------------------------------
  ideal_gas   = trim( adjustl( save_eos_terms(1) ) )
  hard_sphere = trim( adjustl( save_eos_terms(2) ) )
  chain_term  = trim( adjustl( save_eos_terms(3) ) )
  disp_term   = trim( adjustl( save_eos_terms(4) ) )
  hb_term     = trim( adjustl( save_eos_terms(5) ) )
  LC_term     = trim( adjustl( save_eos_terms(6) ) )
  branch_term = trim( adjustl( save_eos_terms(7) ) )
  II_term     = trim( adjustl( save_eos_terms(8) ) )
  ID_term     = trim( adjustl( save_eos_terms(9) ) )
  dd_term     = trim( adjustl( save_eos_terms(10) ) )
  qq_term     = trim( adjustl( save_eos_terms(11) ) )
  dq_term     = trim( adjustl( save_eos_terms(12) ) )

end subroutine RESTORE_PREVIOUS_EOS_NUMERICAL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine paus ( message_to_screen )

  character (LEN=*), optional, intent(in)         :: message_to_screen
  !-----------------------------------------------------------------------------

  if ( present( message_to_screen ) ) then
     write (*,*) 'MESSAGE: ',message_to_screen
  else
     write (*,*) 'PAUSE: press any key to continue'
  end if

  read (*,*)

end subroutine paus

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine error_message ( message_to_screen )

  implicit none
  character (LEN=*)                               :: message_to_screen

  write (*,*) message_to_screen
  stop 6

end subroutine error_message

end module utilities
