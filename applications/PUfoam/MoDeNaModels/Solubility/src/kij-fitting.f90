!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE KIJ_FITTING
!
! vectors of experimental data needed for fitting kij-parameters. A vextor
! giving deviations between experiment and calculations is also included
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module KIJ_FITTING

  implicit none
  save

  PRIVATE
  INTEGER                              :: n_exp
  REAL, DIMENSION(400)                 :: x_exp, y_exp, t_exp, p_exp
  REAL, DIMENSION(400,2)               :: deviation


  PUBLIC :: KIJ_FIT_LIJ, KIJ_FIT

contains


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE KIJ_FIT_LIJ
!
! routine for fitting kij (and lij, if needed) to binary systems
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE KIJ_FIT_LIJ

  USE BASIC_VARIABLES
  use cg_minimization
  use optimizer_2D
  use utilities, only: file_open, paus

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, nadjust, PRIN
  REAL, ALLOCATABLE                      :: optpars(:)
  REAL                                   :: fmin, t0, h0, MACHEP
  REAL                                   :: dummy2, dummy3, dummy4, dummy5
  REAL                                   :: rms1, rms2, aad_m1, aad_m2
  CHARACTER                              :: fitin*50, fitin2*71, fitout*80, dummy1*30
  LOGICAL                                :: go_up

  ! INTEGER                                :: STATUS,  iter, nfunc, ngrad
  ! REAL                                   :: gnorm
  ! REAL, ALLOCATABLE                      :: d(:), g(:), xtemp(:), gtemp(:)
  REAL, ALLOCATABLE                      :: g(:)
  !-----------------------------------------------------------------------------

  ncomp = 2
  nphas = 2
  u_in_P = 1.0
  u_out_P = 1.E6


  write (*,*) ' CHOOSE THE INPUT FILE (./input_file/kij-dat/...) :'
  write (*,*) ' SPECIFY FULL NAME:'
  READ (*,*) fitin
  ! fitin = 'co2-methanol.dat'
  ! fitin = 'co2-heptane.dat'
  ! fitin = 'methane-tetracontane.dat'
  ! fitin = 'n2-acetone.dat'
  ! fitin = 'meoh-octane-lle.dat'
  ! fitin = 'ethygly-hexane-lle.dat'
  ! fitin = 'co2-benzene_all.dat'
  ! fitin = '2-propanol-benzene_318K.dat'
  fitin2 = './input_file/kij-dat/'//fitin

  CALL FILE_OPEN(fitin2,87)

  i = 0
  READ (87,*) compna(1)
  READ (87,*) compna(2)
  go_up  = .true.

  DO WHILE ( go_up )
     READ (87,*) dummy1,dummy2,dummy3,dummy4,dummy5
     go_up  = .false.
     IF (dummy1 /= 'end') THEN
        if ( .NOT.( dummy2 == 1.0 .OR. dummy3 == 1.0 ) .AND.  &
             .NOT.(dummy2 == 0.0 .AND. dummy3 == 0.0 ) ) then
           i = i + 1
           x_exp(i)  = dummy2
           y_exp(i)  = dummy3
           p_exp(i)  = dummy4 * 1.E5
           t_exp(i)  = dummy5  + 273.15
        end if
        go_up     = .true.
     END IF
  END DO
  CLOSE (87)
  n_exp = i

  write (*,*) ' CHOOSE THE EOS :'
  write (*,*) '      0:    SAFT EOS '
  write (*,*) '      1: PC-SAFT EOS '
  write (*,*) '      2:     SRK EOS '
  write (*,*) '      3:      PR EOS '
  read (*,*) eos

  IF (eos == 1) THEN
     write (*,*) ' FURTHER SPECIFY THE EOS :'
     write (*,*) '      0: PC-SAFT EOS '
     write (*,*) '      1: PCP-SAFT EOS '
     write (*,*) '      2: PCIP-SAFT EOS '
     read (*,*) pol
  END IF

  write (*,*) ' CHOOSE THE BINARY PARAMETERS TO BE ADJUSTED :'
  write (*,*) '      1:    kij '
  write (*,*) '      2:    kij AND lij'
  read (*,*) nadjust
  IF ( nadjust == 2 ) THEN
     num = 1
     CALL SET_DEFAULT_EOS_NUMERICAL
     write (*,*) 'calculation conducted with numerical derivatives (press enter)'
     call paus
  END IF
  ALLOCATE ( optpars(nadjust) )

  fitout = fitin
  IF ( eos == 1 ) fitout = './output_file/kij-pc-saft/'//fitout
  IF ( eos == 0 ) fitout = './output_file/kij-saft/'//fitout
  IF ( eos == 2 ) fitout = './output_file/kij-srk/'//fitout
  IF ( eos == 3 ) fitout = './output_file/kij-pr/'//fitout
  OPEN ( 88, FILE=fitout )


  CALL PARA_INPUT

  write (*,'(a,F12.5)') ' current value of kij=',kij(1,2)
  write (*,'(a)') ' do you want to give another starting value? ( 1 : yes, 0 : no)'
  read(*,*) i
  if ( i == 1 ) read (*,*) kij(1,2)

  !-----------------------------------------------------------------------------
  t0 = 1.E-4
  h0 = 0.01
  PRIN = 0
  MACHEP = 1.E-15

  optpars(1) = kij(1,2) + 1.0
  IF (nadjust == 2) optpars(2) = lij(1,2) + 1.0
  !optpars(1) = kij(1,2)
  !IF (nadjust == 2) optpars(2) = lij(1,2)


  IF (nadjust > 1) THEN
     !call PRAXIS2( t0, MACHEP, h0, nadjust, PRIN, optpars, BINDEV_KIJ_p, fmin )
     ALLOCATE ( g(nadjust) )
     CALL Newton_Opt_2D ( BINDEV_KIJ_p, optpars, nadjust, 1.E-8, 1.E-8, g, fmin)
     DEALLOCATE ( g )
     !ALLOCATE( d(nadjust) )
     !ALLOCATE( g(nadjust) )
     !ALLOCATE( xtemp(nadjust) )
     !ALLOCATE( gtemp(nadjust) )
     !CALL cg_descent (1.d-5, optpars, nadjust, BINDEV_KIJ, DEV_KIJ_grad, STATUS, &
     !                  gnorm, fmin,  iter, nfunc, ngrad, d, g, xtemp, gtemp)
  ELSE
     call DIM1MIN( t0, h0, nadjust, optpars, BINDEV_KIJ_p, fmin )
  END IF


  write (*,*)' Ergebnis: kij=',kij(1,2),'  AAD=',fmin / real( n_exp )
  write(88,*)' Ergebnis: kij=',kij(1,2),'  AAD=',fmin / real( n_exp )
  write (*,*)' Ergebnis: lij=',lij(1,2)
  write(88,*)' Ergebnis: lij=',lij(1,2)


  rms1   = 0.0
  rms2   = 0.0
  aad_m1 = 0.0
  aad_m2 = 0.0
  write (88,*) 'experimental data'
  DO i = 1, n_exp
     write (88,'(i4,4(2x,G13.6))') i,x_exp(i),y_exp(i), p_exp(i)/1.E5,t_exp(i)-273.15
  END DO
  write (88,*) '  '
  write (88,*) 'calculated values'
  DO i = 1, n_exp
     write (88,'(i4,4(2x,G13.6))') i,x_exp(i)+deviation(i,1), &
          y_exp(i)+deviation(i,2), &
          p_exp(i)/1.E5,t_exp(i)-273.15
  END DO
  write (88,*) '  '
  DO i = 1, n_exp
     write (*,*)  ' DEV% ph_1=',100.0*deviation(i,1),' DEV% ph_2=',100.0*deviation(i,2)
     write (88,*) ' DEVIATION% phase_1=',100.0*deviation(i,1), &
          '  DEVIATION% phase_2=',100.0*deviation(i,2)
     rms1   = rms1   + deviation(i,1)**2
     rms2   = rms2   + deviation(i,2)**2
     aad_m1 = aad_m1 + ABS( deviation(i,1) )
     aad_m2 = aad_m2 + ABS( deviation(i,2) )
  END DO
  rms1 = 100.0 * ( rms1 / real( n_exp ) )**0.5
  rms2 = 100.0 * ( rms2 / real( n_exp ) )**0.5
  aad_m1 = 100.0 * aad_m1 / real( n_exp )
  aad_m2 = 100.0 * aad_m2 / real( n_exp )

  ! write(*,*)' '
  ! write(*,*)' MAX.DEV% ph_1=',maxerr(1,1),' AT',(maxerr(1,i),i=2,3)
  ! write(*,*)' MAX.DEV% ph_2=',maxerr(2,1),' AT',(maxerr(2,i),i=2,3)
  write(*,*)' '
  write(*,*)' AAD% phase_1=',aad_m1
  write(*,*)' AAD% phase_2=',aad_m2
  write(*,*)' '
  write(*,*)' RMS% phase_1=',rms1
  write(*,*)' RMS% phase_2=',rms2

  ! write(88,*)' '
  ! write(88,*)' MAX. DEVIATION% phase_1=',maxerr(1,1),'  AT',(maxerr(1,i),i=2,3)
  ! write(88,*)' MAX. DEVIATION% phase_2=',maxerr(2,1),'  AT',(maxerr(2,i),i=2,3)
  write(88,*)' '
  write(88,*)' AAD% phase_1=',aad_m1
  write(88,*)' AAD% phase_2=',aad_m2
  write(88,*)' '
  write(88,*)' RMS% phase_1=',rms1
  write(88,*)' RMS% phase_2=',rms2

  write (*,*) ' '
  write(*,*)' Result: kij=',kij(1,2)
  IF ( nadjust == 2 ) write(*,*)'         lij=',lij(1,2)
  write(88,*)' Result: kij=',kij(1,2)
  IF ( nadjust == 2 ) write(88,*)'         lij=',lij(1,2)
  write (*,*) ' '
  write (88,*) ' '

  CLOSE (88)
  DEALLOCATE ( optpars )
  ! DEALLOCATE( d, g, xtemp, gtemp )

END SUBROUTINE KIJ_FIT_LIJ




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BINDEV_KIJ_p ( fmin, optpars, nadjust )

  USE BASIC_VARIABLES
  USE STARTING_VALUES
  USE utilities

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: nadjust
  REAL, INTENT(IN)                       :: optpars(:)
  REAL, INTENT(IN OUT)                   :: fmin

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, converg, iterate_t
  REAL                                   :: error
  REAL                                   :: x_sta, y_sta
  REAL                                   :: dev_x_i, dev_y_i, weigh_x_p, weigh_x_y
  REAL, DIMENSION(400)                   :: dev
  REAL                                   :: xi_save
  real, dimension(nc)                    :: rhoi1, rhoi2
  CHARACTER (LEN=1)                      :: report_success_step
  CHARACTER (LEN=30)                     :: report_character
  !-----------------------------------------------------------------------------


  kij(1,2) = optpars(1) - 1.0
  kij(2,1) = kij(1,2)

  IF (nadjust == 2) lij(1,2) = optpars(2)-1.0
  IF (nadjust == 2) lij(2,1) = - lij(1,2)

  weigh_x_p = 0.25   ! weight of errors in x to errors in p
  weigh_x_y = 0.8    ! weight of errors in x to errors in y
  ! weigh_x_p = 0.0   ! weight of errors in x to errors in p
  ! weigh_x_y = 1.0    ! weight of errors in x to errors in y

  !-----------------------------------------------------------------------------
  ! loop for experimental mixture data
  !-----------------------------------------------------------------------------

  DO i = 1, n_exp

     report_character = ' '
     report_success_step = ' '

     t = t_exp(i)
     p = p_exp(i)

     n_unkw = ncomp            ! number of quantities to be iterated 
     it(1) = 'x11'             ! iteration of mol fraction of comp.1 phase 1
     it(2) = 'x21'             ! iteration of mol fraction of comp.1 phase 2
     sum_rel(1) = 'x12'        ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2) = 'x22'        ! summation relation: x22 = 1 - sum(x2j)

     x_sta = x_exp(i)
     y_sta = y_exp(i)
     if ( x_sta == 0.0 ) x_sta = 0.001   ! these are only starting values
     if ( x_sta == 1.0 ) x_sta = 0.999   ! these are only starting values
     if ( y_sta == 0.0 ) y_sta = 0.001   ! these are only starting values
     if ( y_sta == 1.0 ) y_sta = 0.999   ! these are only starting values

     xi(1,2) = x_sta
     xi(2,2) = y_sta
     xi(1,1) = 1.0 - x_sta
     xi(2,1) = 1.0 - y_sta

     !--------------------------------------------------------------------------
     ! step A for solving the pT-flash problem: Rachford-Rice
     !--------------------------------------------------------------------------

     xiF(1) = 0.5 * ( xi(1,1) + xi(2,1) )
     xiF(2) = 1.0 - xiF(1)

     dense(1) = 0.4
     dense(2) = 0.001

     !outp = 2
     call rachford_rice ( converg, rhoi1, rhoi2 )
     if ( converg == 1 ) then
        CALL occupy_val_init ( rhoi1, rhoi2 )
        val_conv = val_init
        report_success_step = 'A'
     end if


     !--------------------------------------------------------------------------
     ! step B for solving the pT-flash problem: VLE-scan
     !--------------------------------------------------------------------------

     IF ( converg /= 1 ) THEN

        call vle_min
        xiF( 1:ncomp ) = xi( 1, 1:ncomp )
        call start_var_fixed_composition ( converg )
        if ( converg == 1 ) then
           val_conv = val_init
           call restore_converged
           report_success_step = 'B'
        end if

     END IF

     !--------------------------------------------------------------------------
     ! step C for solving the pT-flash problem: x-scan and flash
     !--------------------------------------------------------------------------

     IF ( converg /= 1 ) THEN

        outp = 0                    ! output to terminal
        call scan_compositions ( converg )
        if ( converg == 1 ) then
           val_conv = val_init
           call restore_converged
           report_success_step = 'C'
        end if

     END IF


     !==========================================================================
     ! calculate entry to the objective function
     !==========================================================================

     IF (converg == 1 .AND. ABS(val_conv(1)/val_conv(2)-1.0) > 0.25 ) THEN

        !-----------------------------------------------------------------------
        ! converged to VLE, first calculate the x-deviations and y-deviations
        !-----------------------------------------------------------------------

        if ( val_init(1)/val_init(2) < 0.6 ) then
           write (*,*) 'swap A',i
           xi_save = xi(1,2)
           xi(1,2) = xi(2,2)
           xi(2,2) = xi_save
           xi_save = dense(1)
           dense(1) = dense(2)
           dense(2) = xi_save
        end if

        dev(i) = 0.0
        dev_x_i = 0.0
        dev_y_i = 0.0
        IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) dev_x_i = (xi(1,2)-x_exp(i))**2
        IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) dev_y_i = (xi(2,2)-y_exp(i))**2

        IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) deviation(i,1) = xi(1,2)-x_exp(i)
        IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) deviation(i,2) = xi(2,2)-y_exp(i)

        !-----------------------------------------------------------------------
        ! secondly, calculate the deviation in bubble point pressure
        !-----------------------------------------------------------------------

        xi(1,2) = x_exp(i)
        xi(1,1) = 1.0 - xi(1,2)
        iterate_t = 0

        call bubble_point_rachford_rice ( iterate_t, converg, rhoi1, rhoi2 )
        IF ( converg == 1 ) CALL occupy_val_init ( rhoi1, rhoi2 )
        IF ( converg == 1 ) val_conv = val_init

        !-----------------------------------------------------------------------
        ! now calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------

        if ( converg == 1 ) then
           dev(i) = ( weigh_x_p * 1.0/dev_x_i   &
             + ( 1.0-weigh_x_p )*( val_conv(4)/p_exp(i)-1.0 )**(-2) )**(-1)
           dev(i) = weigh_x_y * dev(i) + ( 1.0 - weigh_x_y ) * dev_y_i
        else
           dev(i) = weigh_x_y * dev_x_i + ( 1.0 - weigh_x_y ) * dev_y_i
        end if
        if ( converg /= 1 ) write (*,*) 'second step not converged, i=',i

     ELSE IF (converg == 1 .AND. ABS(val_conv(1)/val_conv(2)-1.0) <= 0.25 ) THEN

        !-----------------------------------------------------------------------
        ! converged to LLE, first check, whether the phase index has to be exchanged
        !-----------------------------------------------------------------------

        if ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 .AND. abs(xi(1,2)-x_exp(i)) > abs(xi(2,2)-x_exp(i)) ) then
           write (*,*) 'swap1',i
           xi_save = xi(1,2)
           xi(1,2) = xi(2,2)
           xi(2,2) = xi_save
        end if
        if ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 .AND. abs(xi(2,2)-y_exp(i)) > abs(xi(1,2)-y_exp(i)) ) then
           write (*,*) 'swap2',i
           xi_save = xi(1,2)
           xi(1,2) = xi(2,2)
           xi(2,2) = xi_save
        end if

        !-----------------------------------------------------------------------
        ! secondly, calculate the x- and y-deviations as the objectiv function
        !-----------------------------------------------------------------------

        dev(i) = 0.0
        IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) dev(i) = (xi(1,2)-x_exp(i))**2
        IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) dev(i) = dev(i) + (xi(2,2)-y_exp(i))**2

        IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) deviation(i,1) = xi(1,2)-x_exp(i)
        IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) deviation(i,2) = xi(2,2)-y_exp(i)

     ELSE

        !-----------------------------------------------------------------------
        ! no solution found for (p,T) of point i. Now calculate delta(p)
        !-----------------------------------------------------------------------

        dense(1) = 0.4
        dense(2) = 1.E-5
        xi(1,2) = x_exp(i)
        xi(1,1) = 1.0 - xi(1,2)
        xi(2,1:2) = xi(1,1:2)
        t = t_exp(i)
        p = p_exp(i)
        iterate_t = 0
        call bubble_point_rachford_rice ( iterate_t, converg, rhoi1, rhoi2 )
        if ( converg == 1 ) then
           CALL occupy_val_init ( rhoi1, rhoi2 )
           val_conv = val_init
        end if

        !-----------------------------------------------------------------------
        ! calculate the contribution of point i to the objective function
        !-----------------------------------------------------------------------

        if ( converg == 1 ) dev(i) = 5.0*( val_conv(4)/p_exp(i)-1.0 )**2
        if ( converg == 1 ) report_character = 'only delta(p) converged'
        if ( converg /= 1 ) write (*,'(a,i4,a,2G18.11)') ' point number ',i, &
                                      ' did not converge, x,y=',x_exp(i),y_exp(i)
        IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) deviation(i,1) = 0.0  ! temporary
        IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) deviation(i,2) = 0.0  ! temporary


     END IF
     !write (*,'(a,i3,G20.12,2x,a,x,a)') 'deviation',i,dev(i), report_success_step, report_character
     !read (*,*)

  END DO

  error = SUM( dev(1:n_exp) )
  fmin = error

  write(*,'(a,f14.8,a,3(f14.8))') 'deviation=',error,'     kij=',optpars(1:nadjust) - 1.0

END SUBROUTINE BINDEV_KIJ_p



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE KIJ_FIT
!
! routine for fitting kij to binary systems
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE KIJ_FIT

  USE BASIC_VARIABLES
  use utilities, only: file_open, paus

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, opt
  REAL                                   :: dummy2, dummy3, dummy4
  REAL                                   :: kijorg, devges(5)
  REAL                                   :: dist, kijsav(5), devopt, dummy5, rms1, rms2
  REAL                                   :: aad1(400,5), aad2(400,5), aad_m1, aad_m2, maxerr(2,3)
  REAL                                   :: max_t, min_t, max_p, min_p, aver_y
  REAL                                   :: fmin, optpars(1)
  CHARACTER                              :: fitin*50, fitin2*71, fitout*80, dummy1*30
  LOGICAL                                :: go_up
  !-----------------------------------------------------------------------------

  ncomp  = 2
  nphas  = 2
  n_unkw = ncomp
  it(1) = 'x11'                ! iteration of mol fraction of comp.1 phase 1
  it(2) = 'x21'                ! iteration of mol fraction of comp.1 phase 2
  sum_rel(1) = 'x12'           ! summation relation: x12 = 1 - sum(x1j)
  sum_rel(2) = 'x22'           ! summation relation: x22 = 1 - sum(x2j)

  dist   = 0.01

  ! write (*,*) ' CHOOSE YOUR OPTIMIZATION STRATEGY :'
  ! write (*,*) ' OPTIMIZE THE   x AND   y VALUES ...... CHOOSE: 1'
  ! write (*,*) ' OPTIMIZE THE  K1 AND  K2 VALUES ...... CHOOSE: 2'
  ! READ (*,*) c_option

  write (*,*) ' CHOOSE THE INPUT FILE (./input_file/kij-dat/...) :'
  write (*,*) ' SPECIFY FULL NAME:'
  READ (*,*) fitin
  ! fitin='h2s-n-decane.dat'
  ! fitin='methane-tetracontane.dat'
  fitin2 = './input_file/kij-dat/'//fitin

  CALL FILE_OPEN(fitin2,87)

  n_exp = -1
  i = 0
  READ (87,*) compna(1)
  READ (87,*) compna(2)

  go_up  = .true.
  DO WHILE ( go_up )
     READ (87,*) dummy1,dummy2,dummy3,dummy4,dummy5
     n_exp = n_exp + 1
     IF (dummy1 == 'end') THEN
        go_up = .false.
     ELSE
        i = i + 1
        x_exp(i)  = dummy2
        y_exp(i)  = dummy3
        p_exp(i)  = dummy4 * 1.E5
        t_exp(i)  = dummy5 + 273.15
     ENDIF
  END DO
  CLOSE (87)

  write (*,*) ' CHOOSE THE EOS :'
  write (*,*) '      0:    SAFT EOS '
  write (*,*) '      1: PC-SAFT EOS '
  write (*,*) '      2:     SRK EOS '
  write (*,*) '      3:      PR EOS '
  read (*,*) eos

  IF (eos == 1) THEN
     write (*,*) ' FURTHER SPECIFY THE EOS :'
     write (*,*) '      0: PC-SAFT EOS '
     write (*,*) '      1: PCP-SAFT EOS '
     write (*,*) '      2: PCIP-SAFT EOS '
     read (*,*) pol
     IF(pol == 2 .AND. num == 0) write(*,*)'switch to num. derivatives'
     IF(pol == 2 .AND. num == 0) num=1
  ENDIF

  fitout = fitin
  IF (eos == 1) fitout = './output_file/kij-pc-saft/'//fitout
  IF (eos == 0) fitout = './output_file/kij-saft/'//fitout
  IF (eos == 2) fitout = './output_file/kij-srk/'//fitout
  IF (eos == 3) fitout = './output_file/kij-pr/'//fitout
  OPEN (88,FILE=fitout)

  CALL PARA_INPUT


  kij_adjust: DO

     kijorg = kij(1,2)

     DO j = 1,5

        IF (j == 1) THEN
           kij(1,2) = kijorg - 1.0 * dist
        ELSEIF (j == 2) THEN
           kij(1,2) = kijorg - 0.5 * dist
        ELSEIF (j == 3) THEN
           kij(1,2) = kijorg 
        ELSEIF (j == 4) THEN
           kij(1,2) = kijorg + 0.5 * dist
        ELSE
           kij(1,2) = kijorg + 1.0 * dist
        ENDIF
        kijsav(j) = kij(1,2)
        kij(2,1) = kij(1,2)
        optpars(1) = kij(1,2) + 1.0
        ! optpars(1) = kij(1,2)

        call BINDEV_KIJ ( fmin, optpars, 1 )
        devges(j) = fmin
        write (*,'(a,i3,2(f14.8))') '   deviation',j,devges(j),kij(1,2)

        DO i = 1, n_exp
           aad1(i,j) = deviation(i,1)
           aad2(i,j) = deviation(i,2)
        ENDDO

     END DO

     devopt   = devges(1)
     kij(1,2) = kijsav(1)
     kij(2,1) = kijsav(1)
     opt      = 1
     DO k = 2,5
        IF (devges(k) < devopt) THEN
           kij(1,2) = kijsav(k)
           kij(2,1) = kijsav(k)
           devopt   = devges(k)
           opt      = k
        ENDIF
     END DO

     IF ( opt /= 1 .AND. opt /= 5 ) dist = 0.25*dist

     write(*,'(a,f14.8,a,2(f14.8))')' new kij',kij(1,2),' +/-',dist,devopt



     IF (dist <= 0.00005) EXIT kij_adjust

  END DO kij_adjust

  rms1   = 0.0
  rms2   = 0.0
  aad_m1 = 0.0
  aad_m2 = 0.0
  maxerr(1,1) = 100.0 * aad1(1,opt)
  maxerr(1,2) = t_exp(1)
  maxerr(1,3) = p_exp(1)
  maxerr(2,1) = 100.0 * aad2(1,opt)
  maxerr(2,2) = t_exp(1)
  maxerr(2,3) = p_exp(1)
  DO i = 1, n_exp
     write (88,'(i4,4(2x,G13.6))') i,x_exp(i),y_exp(i),p_exp(i)/1.E5,t_exp(i)-273.15
  ENDDO
  write (88,*) '  '
  DO i = 1, n_exp
     IF( (100.0 * ABS(aad1(i,opt))) > maxerr(1,1) ) THEN
        maxerr(1,1) = 100.0 * ABS(aad1(i,opt))
        maxerr(1,2) = t_exp(i)
        maxerr(1,3) = p_exp(i)
     ENDIF
     IF( (100.0 * ABS(aad2(i,opt))) > maxerr(2,1) ) THEN
        maxerr(2,1) = 100.0 * ABS(aad2(i,opt))
        maxerr(2,2) = t_exp(i)
        maxerr(2,3) = p_exp(i)
     ENDIF
     write (*,*)  ' DEV% ph_1=', 100.0 * aad1(i,opt),' DEV% ph_2=',100.0 * aad2(i,opt)
     write (88,*) ' DEVIATION% phase_1=',100.0 * aad1(i,opt), &
          '  DEVIATION% phase_2=',100.0 * aad2(i,opt)
     rms1   = rms1   + aad1(i,opt)**2
     rms2   = rms2   + aad2(i,opt)**2
     aad_m1 = aad_m1 + ABS(aad1(i,opt))
     aad_m2 = aad_m2 + ABS(aad2(i,opt))
  END DO
  rms1 = 100.0 * ( rms1 / real( n_exp ) )**0.5
  rms2 = 100.0 * ( rms2 / real( n_exp ) )**0.5
  aad_m1 = 100.0 * aad_m1 / real( n_exp )
  aad_m2 = 100.0 * aad_m2 / real( n_exp )

  write(*,*)' '
  write(*,*)' MAX.DEV% ph_1=',maxerr(1,1),' AT',(maxerr(1,i),i=2,3)
  write(*,*)' MAX.DEV% ph_2=',maxerr(2,1),' AT',(maxerr(2,i),i=2,3)
  write(*,*)' '
  write(*,*)' AAD% phase_1=',aad_m1
  write(*,*)' AAD% phase_2=',aad_m2
  write(*,*)' '
  write(*,*)' RMS% phase_1=',rms1
  write(*,*)' RMS% phase_2=',rms2

  write(88,*)' '
  write(88,*)' MAX. DEVIATION% phase_1=',maxerr(1,1),'  AT',(maxerr(1,i),i=2,3)
  write(88,*)' MAX. DEVIATION% phase_2=',maxerr(2,1),'  AT',(maxerr(2,i),i=2,3)
  write(88,*)' '
  write(88,*)' AAD% phase_1=',aad_m1
  write(88,*)' AAD% phase_2=',aad_m2
  write(88,*)' '
  write(88,*)' RMS% phase_1=',rms1
  write(88,*)' RMS% phase_2=',rms2

  write (*,*) ' '
  write(*,*)' Result: kij=',kij(1,2)
  write(88,*)' Result: kij=',kij(1,2)
  write (*,*) ' '
  write (88,*) ' '

  max_t=t_exp(1)
  min_t=t_exp(1)
  max_p=p_exp(1)
  min_p=p_exp(1)
  aver_y = y_exp(1) / real( n_exp )
  DO i = 2, n_exp
     IF (t_exp(i) > max_t) max_t = t_exp(i)
     IF (t_exp(i) < min_t) min_t = t_exp(i)
     IF (p_exp(i) > max_p) max_p = p_exp(i)
     IF (p_exp(i) < min_p) min_p = p_exp(i)
     aver_y = aver_y + y_exp(i) / real( n_exp )
  END DO
  write (88,'(a,1x,12(1x,e12.5),1x,a)') fitin,kij(1,2),aad_m1 &
       ,aad_m2,rms1,rms2,maxerr(1,1),maxerr(2,1),min_t-273.15 &
       ,max_t-273.15,min_p/1.E5,max_p/1.E5,aver_y
  write (*,*) ' '

  CLOSE (88)

END SUBROUTINE KIJ_FIT


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BINDEV_KIJ ( fmin, optpars, nadjust )

  USE BASIC_VARIABLES
  USE STARTING_VALUES

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: nadjust
  REAL, INTENT(IN)                       :: optpars(:)
  REAL, INTENT(IN OUT)                   :: fmin

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, k, ph, ph_split, converg
  REAL                                   :: error
  REAL                                   :: x_sta, y_sta
  REAL, DIMENSION(400)                   :: dev
  REAL                                   :: end_x, steps
  !-----------------------------------------------------------------------------


  kij(1,2) = optpars(1) - 1.0
  ! kij(1,2) = optpars(1)
  kij(2,1) = kij(1,2)
  IF (nadjust == 2) lij(1,2) = optpars(2)-1.0
  ! IF (nadjust == 2) lij(1,2) = optpars(2)
  IF (nadjust == 2) lij(2,1) = - lij(1,2)

  DO i = 1, n_exp

     t = t_exp(i)
     p = p_exp(i)

     n_unkw = ncomp            ! number of quantities to be iterated 
     it(1) = 'x11'             ! iteration of mol fraction of comp.1 phase 1
     it(2) = 'x21'             ! iteration of mol fraction of comp.1 phase 2
     sum_rel(1) = 'x12'        ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2) = 'x22'        ! summation relation: x22 = 1 - sum(x2j)

     x_sta = x_exp(i)
     y_sta = y_exp(i)
     if ( x_sta == 0.0 ) x_sta = 0.01   ! these are only starting values
     if ( x_sta == 1.0 ) x_sta = 0.99   ! these are only starting values
     if ( y_sta == 0.0 ) y_sta = 0.01   ! these are only starting values
     if ( y_sta == 1.0 ) y_sta = 0.99   ! these are only starting values

     lnx(1,2) = LOG( x_sta )
     lnx(2,2) = LOG( y_sta )
     lnx(1,1) = LOG( 1.0 - x_sta )
     lnx(2,1) = LOG( 1.0 - y_sta )
     val_init(1) = 0.4
     val_init(2) = 0.001
     val_init(3) = t
     val_init(4) = p
     DO ph = 1,nphas
        DO k = 1,ncomp
           val_init(4+k+(ph-1)*ncomp) = lnx(ph,k)
        ENDDO
     ENDDO

     CALL OBJECTIVE_CTRL (converg)

     IF (converg /= 1) THEN

        outp = 0                    ! output to terminal
        call vle_min
        xiF( 1:ncomp ) = xi( 1, 1:ncomp )
        call start_var_fixed_composition ( converg )

        IF (converg /= 1) THEN
           write (*,*) ' NO SOLUTION FOUND FOR REG. PRESSURE',i
           flashcase = .true.
           xiF(2) = ( x_exp(i) + y_exp(i) ) / 2.0
           xiF(1) = 1.0 - xiF(2)
           p = p_exp(i) - 0.1*p_exp(i)
           val_init(4) = p
           ph_split = 0
           call start_var_fixed_composition ( converg )
           IF (converg == 1) THEN
              steps = 10.0
              end_x = p_exp(i)
              running = 'p'
              call phase_equilib (end_x,steps,converg)
           END IF
           IF (converg == 1) write (*,*) ' NOW A SOLUTION WAS FOUND ?, CHECK !',i,p_exp(i)
        END IF
        ! IF (converg /= 1) write (*,*) ' NO SOLUTION FOUND',i,p_exp(i)
     ENDIF


     dev(i) = 0.0
     IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) dev(i) =          (xi(1,2)-x_exp(i))**2
     IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) dev(i) = dev(i) + (xi(2,2)-y_exp(i))**2

     IF ( x_exp(i) /= 0.0 .AND. x_exp(i) /= 1.0 ) deviation(i,1) = xi(1,2)-x_exp(i)
     IF ( y_exp(i) /= 0.0 .AND. y_exp(i) /= 1.0 ) deviation(i,2) = xi(2,2)-y_exp(i)

     ! write (*,*) 'converg A',converg, xi(1,2)-x_exp(i), deviation(i,1)

  END DO
  !call paus

  error = SUM( dev(1:n_exp) )
  fmin = error

  write(*,'(a,f14.8,a,3(f14.8))') 'deviation=',error,'     kij=',optpars(1:nadjust) - 1.0

END SUBROUTINE BINDEV_KIJ


end module KIJ_FITTING

