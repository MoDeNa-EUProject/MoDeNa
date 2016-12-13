!> \file starting_value.f90
!! \brief This subroutine performs a phase equilibrium calculation.




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! module STARTING_VALUES
!
! This module contains parameters and variables for a phase stability
! analyis as part of a flash calculation.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

module STARTING_VALUES

  use PARAMETERS, only: nc
  implicit none
  save

  integer                               :: scan_index
  real, dimension(nc), public           :: rhoi_trial
  real                                  :: fdenf
  logical                               :: first_call = .true.
  logical, public                       :: generate_starting_val
  logical, public                       :: flashcase

  real, dimension(nc)                   :: rhoi_best
  real                                  :: fmin_best

  PRIVATE
  PUBLIC :: start_var, start_var_fixed_composition, phase_stability, vle_min,  &
            tangent_plane_line_search, rachford_rice, occupy_val_init,  &
            scan_compositions, bubble_point_rachford_rice, stability_hessian

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!> subroutine start_var
!>
!> This subroutine generates a converged solution for binary systems
!> or performes a flash calculation for mixtues.
!>
!> IF a polymer is considered, starting values for mole fractions
!> are determined from the SUBROUTINGE POLY_STA_VAR (see below). The
!> polymer needs to be placed as component 1 (first line) in INPUT
!> file.
!>
!> A phase equilib. iteration is started at the end of this routine.
!> If no solution is found (converg=0), the program will stop within
!> this routine.
!>
!> Currently, this routine assumes two-phase equilibrium and derives
!> starting values (xi,density) only for two phases.
!>
!> Prerequisites are:
!> subroutine INPUT needs to be called prior to this routine, because
!> all pure comp. parameters as well as (T,P,kij) need to be in place.
!> Also, the variable to be iterated "it(i)" and the variables to be
!> calculated through the summation relation "sum_rel(i)" have to be
!> defined.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine start_var (converg)

  use BASIC_VARIABLES
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                 :: converg

  !-----------------------------------------------------------------------------
  logical                                 :: renormalize
  !-----------------------------------------------------------------------------

  converg = 0
  nphas = 2
  n_unkw = ncomp
  outp = 0                    ! output to terminal
  generate_starting_val = .true.

  ensemble_flag = 'tp'

  if ( first_call .AND. xif(1) == 0.0 ) call read_mode_staring_value ( scan_index, outp )

  renormalize = .false.       ! for renormalization group theory (RGT)
  if (num == 2) renormalize = .true.
  if (num == 2) num = 0       ! if RGT: initial phase equilibr. is for non-renormalized model

  if ( xif(1) /= 0.0 ) then

     flashcase = .true.

     call start_var_fixed_composition ( converg )

  else

     flashcase = .true.

     if ( ncomp > 2 ) then
        write (*,*) 'Solubility Code: SR starting_value requires a defined feed composition for all but binary mixtures'
        stop 5
     end if

     if ( scan_index == 2 ) then
        call vle_min
        xiF( 1:ncomp ) = xi( 1, 1:ncomp )
        call start_var_fixed_composition ( converg )
     else
        call scan_compositions ( converg )
     end if

     xif(:) = 0.0

  end if

  IF (renormalize) num = 2
  generate_starting_val = .false.

end subroutine start_var


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine scan_compositions
!
! this SR initiates initiates flash calculations along a grid of compositions.
! Only for binary mixtures and for given t, p.
!
! input:    t, p
! output:   converg
!           val_init ( !!! )    if a converged solution is found
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine scan_compositions ( converg )

  use BASIC_VARIABLES
  use utilities
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                 :: converg

  !-----------------------------------------------------------------------------
  integer                                 :: i, i_conv, steps
  integer                                 :: converg_overall
  real                                    :: x1_right
  real                                    :: diff_para
  real                                    :: xiF_in_between_x
  real, dimension( 4, 0:nc*np+6 )         :: val_save
  logical                                 :: different_equilibrium_as_previous
  real                                    :: density(np),w(np,nc)
  !-----------------------------------------------------------------------------

  steps = 100
  x1_right = 0.0
  converg_overall = 0
  if ( ncomp > 2 ) write (*,*) 'SR scan_compositions only suited for binary systems'
  if ( ncomp > 2 ) stop 5

  do i = 0, steps

     xiF(1) = REAL(i) / REAL(steps)
     if ( xiF(1) <= 1.E-30 ) xiF(1) = 1.E-30
     if ( xiF(1) >= (1.0 - 1.E-12) ) xiF(1) = 1.0 - 1.E-12
     xiF(2)  = 1.0 - xiF(1)

     if ( xiF(1) > x1_right ) then

        if ( outp >= 1 ) write (*,'(a,G19.11)') ' scan x, with x(1)=',xiF(1)

        call start_var_fixed_composition ( converg )

        xiF_in_between_x = ( EXP(val_init(5)) - xiF(1) ) / (xiF(1) - EXP(val_init(7)) )

        if ( converg == 1 .AND. xiF_in_between_x > 0.0 ) then

           different_equilibrium_as_previous = .true.
           if ( converg_overall > 0 ) then
              diff_para =  ABS( exp(val_save(converg_overall,5))-exp(val_init(5)) )  &
                         + ABS( exp(val_save(converg_overall,7))-exp(val_init(7)) )
              if ( diff_para < 1.E-5 ) different_equilibrium_as_previous = .false.
           end if

           if ( different_equilibrium_as_previous ) then
             converg_overall = converg_overall + 1
             val_save( converg_overall, : ) = val_init(:)
             x1_right = max( xi(1,1), xi(2,1) )
             if ( outp >= 1 ) write (*,'(a,2i4,5G19.11)') ' found equil.',i,converg_overall,  &
                                      EXP(val_init(5) ),xiF(1),EXP(val_init(7) ),val_init(1:2)
             if ( outp >= 1 ) write (*,*) ' '
             ! call paus ('end of scan_compositions cycle')

           end if
        end if

     end if

  end do

  if ( converg_overall > 0) converg = 1

  if ( outp >= 1 .OR. converg_overall > 1 ) then
  if ( converg_overall >= 1) write (*,*) ' '
  if ( converg_overall >= 1) write (*,*) '================================='
  if ( converg_overall == 1) write (*,*) 'phase equilibrium:'
  if ( converg_overall == 2) write (*,*) 'two phase equilibria found:'
  if ( converg_overall >= 1) write (*,*) '================================='
  do i_conv = 1, converg_overall
     dense(1:2) = val_save( i_conv, 1:2 )
     lnx(1,1:2) = val_save( i_conv, 5:6 )
     lnx(2,1:2) = val_save( i_conv, 7:8 )
     xi(1,1:2) = exp( lnx(1,1:2) )
     xi(2,1:2) = exp( lnx(2,1:2) )
     call si_dens ( density, w )
     write (*,*) ' '
     if ( converg_overall > 1) write (*,*) 'EQUILIBRIUM', i_conv
     if ( converg_overall > 1) write (*,*) '---------------------------------------------'
     write(*,'(t20,a,f7.2,a,f10.5,a)') 'T =',t-u_out_t,'   P =',p/u_out_p
     write (*,'(t26,a)') 'PHASE I           PHASE II   '
     write (*,'(a,t20,F13.4,t39,F13.4)') ' Density     ',density(1),density(2)
     write (*,'(x,a,t20,a3,E14.6,t42,E14.6)') compna(1),'x =',EXP( lnx(1:2,1) )
     write (*,'(x,a,t20,a3,E14.6,t42,E14.6)') compna(2),'x =',EXP( lnx(1:2,2) )
     write (*,'(x,a,t20,a3,E14.6,t42,E14.6)') compna(1),'w =',w(1:2,1)
     write (*,'(x,a,t20,a3,E14.6,t42,E14.6)') compna(2),'w =',w(1:2,2)
  end do
  end if

  val_init( : ) = val_save( 1, : )
  if ( converg_overall > 1) then
     write (*,*) ' '
     write (*,*) 'choose a phase equilibrium'
     write (*,*) '---------------------------------------------'
     read (*,*) i_conv
     val_init( : ) = val_save( i_conv, : )
  end if

end subroutine scan_compositions


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine calculate_equilibrium_feed
!
! input:    t, p
!           rhoi1, rhoi2 ( intent(in), starting values )
!
! output:   converg, val_init (meaningful, if converg=1)
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine calculate_equilibrium_feed ( rhoi1, rhoi2, converg )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: PI, mseg, dhs
  use utilities
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                 :: converg
  real, dimension(nc), intent(in out)     :: rhoi1
  real, dimension(nc), intent(in out)     :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                 :: promising_min
  real                                    :: phi_2_start
  real, dimension(np,nc)                  :: rho_enter
  !-----------------------------------------------------------------------------

  rho_enter(1,1:ncomp) = rhoi1(1:ncomp)
  rho_enter(2,1:ncomp) = rhoi2(1:ncomp)

  xi(1,1:ncomp) = rhoi1(1:ncomp) / sum(rhoi1(1:ncomp))
  xi(2,1:ncomp) = rhoi2(1:ncomp) / sum(rhoi2(1:ncomp))
  dense(1) = PI/6.0 * SUM( rhoi1(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
  dense(2) = PI/6.0 * SUM( rhoi2(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  CALL rachford_rice ( converg, rhoi1, rhoi2 )

  if ( converg == 0 ) then

     rhoi1(1:ncomp) = rho_enter(1,1:ncomp)
     rhoi2(1:ncomp) = rho_enter(2,1:ncomp)
     call tangent_plane_line_search ( rhoi1, rhoi2, phi_2_start, promising_min )

     if ( promising_min == 1 .AND. maxval( ABS( xi(1,1:ncomp) - xi(2,1:ncomp) ) ) > 0.0001 ) then

        call tangent_plane_2 ( phi_2_start, rhoi1, rhoi2 )
        xi( 1, 1:ncomp ) = rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
        xi( 2, 1:ncomp ) = rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
        if ( maxval( ABS( xi( 1, 1:ncomp ) - xi( 2, 1:ncomp ) ) ) > 0.0001 ) then

           CALL occupy_val_init ( rhoi1, rhoi2 )
           if ( .NOT. flashcase ) then
              CALL select_sum_rel (1,0,1)
              CALL select_sum_rel (2,0,2)
           else
              CALL determine_flash_it2
           end if

           CALL objective_ctrl ( converg )

           if ( converg == 1 ) val_init = val_conv

        end if

     end if

  else

     CALL occupy_val_init ( rhoi1, rhoi2 )
     val_conv = val_init

  end if

  if ( outp >= 2 ) write (*,*) ' '
  if ( outp >= 2 ) write (*,*) 'leaving calculate_equilibrium_feed: convergence=',converg

end subroutine calculate_equilibrium_feed



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine start_var_fixed_composition
!
! This subroutine performs a flash calculation for a mixtue with defined compo-
! sition. The density is not a priori defined.
!
! input:    t, p, xiF(:)
! output:   converg
!           val_init ( !!! )    if a converged solution is found
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine start_var_fixed_composition ( converg )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: PI, kbol
  use utilities
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                 :: converg

  !-----------------------------------------------------------------------------
  integer                                 :: ph_split
  real, dimension(np,nc)                  :: lnphi
  real, dimension(nc)                     :: rhoi_feed
  real, dimension(nc)                     :: rhoi1, rhoi2
  real, dimension(nc)                     :: rhoi_trial
  real                                    :: eta_trial
  real                                    :: phi_1, fmin
  !-----------------------------------------------------------------------------

  if ( sum( xiF(1:ncomp) ) /= 1.0 ) write (*,*) 'Solubility Code: feed composition /= 1.0', sum( xiF(1:ncomp) )
  if ( sum( xiF(1:ncomp) ) /= 1.0 ) stop 5

  xi(1,1:ncomp) = xiF(1:ncomp)
  xi(2,1:ncomp) = xiF(1:ncomp)

  !-----------------------------------------------------------------------------
  ! start with a Rachford-Rice iteration that searches for VLE. The attempt is
  ! naive, because standard starting values for dense(1), dense(2), and
  ! xi(2,:)=xiF(:) are used.
  !-----------------------------------------------------------------------------

  dense(1) = 0.4
  dense(2) = 1.E-6
  CALL rachford_rice ( converg, rhoi1, rhoi2 )

  !-----------------------------------------------------------------------------
  ! distinguish: case (1), RR did not converge. Continue with a conventional
  !              stability analysis
  !              case (2), RR converged. Then we have to verify that the liquid
  !              phase is stable.
  !-----------------------------------------------------------------------------

  if ( converg /= 1 ) then

     !--------------------------------------------------------------------------
     ! No solution was found with the Rachford-Rice procedure. Now do a
     ! traditional stability test first. If two densities are possible at
     ! xiF(1:ncomp), then do the phase stability for both "feed" conditions.
     !--------------------------------------------------------------------------
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
     if ( outp >= 1 ) write (*,*) 'do phase stability analysis'
     if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
     if ( outp >= 2 ) write (*,*) ' '

     densta(1) = 0.4
     densta(2) = 1.E-6
     xi( 1, 1:ncomp ) = xiF( 1:ncomp )
     xi( 2, 1:ncomp ) = xiF( 1:ncomp )
     CALL fugacity (lnphi)

     !--------------------------------------------------------------------------
     ! phase stability based on liquid "feed"
     !--------------------------------------------------------------------------
     my_f( 1:ncomp ) = my_cal(1,1:ncomp)
     rhoi_feed( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )

     eta_trial = 0.4
     ! to do: check wheather eta_trial = dense(1), or (2) is better choice
     CALL phase_stability ( rhoi_feed, eta_trial, ph_split, rhoi_trial )
     rhoi1( 1:ncomp ) = rhoi_feed( 1:ncomp )
     rhoi2( 1:ncomp ) = rhoi_trial( 1:ncomp )
     if ( outp > 0 ) write (*,*) 'phase split =', ph_split
     if ( ph_split == 0 .AND. outp > 0 ) write (*,*) 'no phase split detected, stable phase is likely'
     if ( ph_split == 2 .AND. outp > 0 ) write (*,*) 'condition could be close to critical point'

     if ( ph_split >= 1 ) then

        !-----------------------------------------------------------------------
        ! If phase_stability gave ph_split = 1 or 2: do tangent plane iteration
        !-----------------------------------------------------------------------

        call calculate_equilibrium_feed ( rhoi1, rhoi2, converg )

        if ( converg == 1 .AND. outp > 0 ) write (*,*) 'end of starting value - solution found'

     end if

     if ( converg == 0 ) then

        !-----------------------------------------------------------------------
        ! test phase stability based on vapor "feed"
        !-----------------------------------------------------------------------
!!$        densta(1) = 0.4
!!$        densta(2) = 1.E-6
!!$        xi( 1, 1:ncomp ) = xiF( 1:ncomp )
!!$        xi( 2, 1:ncomp ) = xiF( 1:ncomp )
!!$        CALL fugacity (lnphi)
!!$
!!$        my_f( 1:ncomp ) = my_cal( 2, 1:ncomp )
!!$        eta_trial = dense(2)
!!$        rhoi_feed( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )
!!$
!!$        CALL phase_stability ( rhoi_feed, eta_trial, ph_split, rhoi_trial )
!!$        rhoi1( 1:ncomp ) = rhoi_feed( 1:ncomp )
!!$        rhoi2( 1:ncomp ) = rhoi_trial( 1:ncomp )
!!$        if ( outp > 0 ) write (*,*) 'phase split =', ph_split
!!$        if ( ph_split == 0 .AND. outp > 0 ) write (*,*) 'no phase split',  &
!!$                                                        ' detected, most probably a stable phase'
!!$        if ( ph_split == 2 .AND. outp > 0 ) write (*,*) 'condition could be close to critical point'
!!$
!!$        if ( ph_split >= 1 ) then
!!$
!!$           call calculate_equilibrium_feed ( rhoi1, rhoi2, converg )
!!$
!!$        end if
!!$
!!$        if ( converg == 0 .AND. ph_split == 1 .AND. outp > 0 ) &
!!$                       write (*,*) 'vapor phase is unstable, but phase split did not converge'
     end if

  else

     !--------------------------------------------------------------------------
     ! RR was successful. VLE was calculated. Now, test stability of liquid phase
     !--------------------------------------------------------------------------
     if ( outp >= 2 ) write (*,*) ' now test phase stability'
     if ( outp >= 2 ) write (*,*) ' '

     my_f( 1:ncomp ) = my_cal(1,1:ncomp)
     eta_trial = 0.45
     rhoi_feed( 1:ncomp ) = rhoi1( 1:ncomp )
     if ( dense(2) > dense(1) ) rhoi_feed( 1:ncomp ) = rhoi2( 1:ncomp )

     CALL phase_stability ( rhoi_feed, eta_trial, ph_split, rhoi_trial )
     if ( outp > 0 ) write (*,*) 'result of stability-test for L-phase, ph_split =', ph_split

     if ( ph_split /= 1 ) then

        !-----------------------------------------------------------------------
        ! for the calculated VLE: liquid phase is stable
        !-----------------------------------------------------------------------
        if ( outp >= 2 ) write (*,*) 'the calculated VLE is stable'
        CALL occupy_val_init ( rhoi1, rhoi2 )

     else

        !-----------------------------------------------------------------------
        ! for the calculated VLE: liquid phase is not stable.
        ! Determine stable equilibrium. Use the density vector rhoi_trial (de-
        ! termined from phase_stability) as a starting value. The starting value
        ! of second phase is determined from the component balance (using the
        ! L-phase and afterwards the V-phase of the previously determined VLE)
        !-----------------------------------------------------------------------
        if ( outp >= 2 ) write (*,*) '========================================================'
        if ( outp >  0 ) write (*,*) ' converged VLE (or LLE) not stable: searching (other) LLE'
        if ( outp >= 2 ) write (*,*) '========================================================'
        if ( outp >= 2 ) write (*,*) ' '

        !-----------------------------------------------------------------------
        ! temporarily save the converged VLE solution into val_init
        !-----------------------------------------------------------------------
        call  occupy_val_init ( rhoi1, rhoi2 )

        rhoi1( 1:ncomp ) = rhoi_feed( 1:ncomp )
        rhoi2( 1:ncomp ) = rhoi_trial( 1:ncomp )
        ! I am not sure the following line (and the condition 0<phi_1<1) is helpful. Why did I introduce this?
        phi_1 = ( xif(1) - rhoi2(1)/sum(rhoi2(1:ncomp)) )  &
              / ( rhoi1(1)/sum(rhoi1(1:ncomp)) - rhoi2(1)/sum(rhoi2(1:ncomp)) )

        if ( outp > 0 ) write (*,*) 'initial phase fraction', phi_1, kij(1,2)

        !if ( phi_1 > 0.0 .AND. phi_1 < 1.0 ) then

           call calculate_equilibrium_feed ( rhoi1, rhoi2, converg )

           if ( converg == 1 ) then
              phi_1 = ( xif(1) - xi(2,1) ) / ( xi(1,1) - xi(2,1) )
              fmin = phi_1 * gibbs(1) + ( 1.0 - phi_1 ) * gibbs(2)
              if ( outp >= 2 ) write (*,*) ' '
              if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
              if ( outp > 0  ) write (*,*) 'first phase equilibrium trial converged ',fmin
              if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
              if ( outp >= 2 ) write (*,*) ' '
           else
              if ( outp >= 2 ) write (*,*) 'search for additional phase split gave no solution'
           end if
        !else
        !   if ( outp >= 2 ) write (*,*) 'first phase equilibrium trial did not qualify'
        !end if

!!$        !rhoi_test1( 1:ncomp ) = rhoi_trial( 1:ncomp )
!!$        !rhoi_test2( 1:ncomp ) = rhoi_save( 1:ncomp )
!!$        !phi_1 = ( xif(1) - rhoi_test2(1)/sum(rhoi_test2(1:ncomp)) )  &
!!$        !     / ( rhoi_test1(1)/sum(rhoi_test1(1:ncomp)) - rhoi_test2(1)/sum(rhoi_test2(1:ncomp)) )
!!$        if( outp > 0 ) write (*,*) 'initial phase fraction', phi_1
!!$
!!$        if ( phi_1 > 0.0 .AND. phi_1 < 1.0 ) then
!!$
!!$           !rhoi1( 1:ncomp ) = rhoi_test1( 1:ncomp )
!!$           !rhoi2( 1:ncomp ) = rhoi_test2( 1:ncomp )
!!$           call calculate_equilibrium_feed ( rhoi1, rhoi2, converg ) ! rhoi_test1, rhoi_test2
!!$
!!$           converg_all = converg_all + converg
!!$
!!$           if ( converg == 1 ) then
!!$              phi_1 = ( xif(1) - xi(2,1) ) / ( xi(1,1) - xi(2,1) )
!!$              fmin = phi_1 * gibbs(1) + ( 1.0 - phi_1 ) * gibbs(2)
!!$              if ( outp >= 2 ) write (*,*) ' '
!!$              if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
!!$              if ( outp > 0 ) write (*,*) 'second phase equilibrium trial converged ',fmin
!!$              if ( outp >= 2 ) write (*,*) '--------------------------------------------------------'
!!$              if ( outp >= 2 ) write (*,*) ' '
!!$           else
!!$              if ( outp >= 2 ) write (*,*) 'second phase equilibrium trial gave no solution'
!!$           end if
!!$        else
!!$           if ( outp >= 2 ) write (*,*) 'second phase equilibrium trial did not qualify'
!!$        end if

        if ( converg == 0 ) then  ! this is done to recover the VLE solution of Rachford-Rice
           !--------------------------------------------------------------------
           ! the search for additional phase split did not deliver stable phases.
           ! Restore the initially converged results -- no particular action
           ! needed because the results are already stored in "val_init".
           !--------------------------------------------------------------------
           converg = 1
        end if

     end if

  end if
  if ( outp > 0 ) write (*,*) ' '

end subroutine start_var_fixed_composition



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!  input:   rhoi1, rhoi2
!
!  output:  rhoi1, rhoi2, xi(1,:), xi(2,:)
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_plane_line_search ( rhoi1, rhoi2, phi_2_opt, promising_min )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: PI, mseg, dhs
  implicit none

  !-----------------------------------------------------------------------------
  real, dimension(nc), intent(in out)     :: rhoi1
  real, dimension(nc), intent(in out)     :: rhoi2
  real, intent(out)                       :: phi_2_opt
  integer, intent(out)                    :: promising_min

  !-----------------------------------------------------------------------------
  integer                                 :: n, i, steps
  real                                    :: fmin
  real                                    :: phi_2, phi_2_min, phi_2_max, fmin_scan
  real, allocatable                       :: optpara(:)
  !-----------------------------------------------------------------------------

  phi_2_opt = 0.5
  promising_min = 0

  n = ncomp
  ALLOCATE( optpara(n) )

  xi( 1, 1:ncomp ) = rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
  xi( 2, 1:ncomp ) = rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
  lnx( 1, 1:ncomp ) = LOG( xi( 1, 1:ncomp ) )
  lnx( 2, 1:ncomp ) = LOG( xi( 2, 1:ncomp ) )

  densta(1) = PI / 6.0 * sum( rhoi1( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )
  densta(2) = PI / 6.0 * sum( rhoi2( 1:ncomp ) * mseg( 1:ncomp ) * dhs( 1:ncomp )**3 )
  !pause

  phi_2_min = 0.001
  phi_2_max = 1.0
  do i = 1, ncomp
    phi_2_max = min( phi_2_max, xif(i) / xi(2,i) )
  end do

  steps = 50
  fmin_scan = 1.E10

  do i = 0, steps

     phi_2 = phi_2_min + ( phi_2_max - phi_2_min ) * real(i) / real(steps)
     if ( phi_2 < 1.E-6 ) phi_2 = 1.E-6
     if ( phi_2 > (1.0 - 1.E-6) ) phi_2 = 1.0 - 1.E-6

     optpara( 1:ncomp ) = LOG( xi( 2, 1:ncomp ) * phi_2 )

     call tangent_value ( fmin, optpara, n )

     if ( fmin < fmin_scan ) then
        fmin_scan = fmin
        phi_2_opt = phi_2
        rhoi1( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )
        rhoi2( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )
        if ( i > 0 ) promising_min = 1
     end if

  end do

  DEALLOCATE( optpara )

  xi( 1, 1:ncomp ) = rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
  xi( 2, 1:ncomp ) = rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )

  if ( outp >= 1 ) write (*,*) 'finished line search'
  if ( outp >= 2 ) write (*,*) ' '

end subroutine tangent_plane_line_search



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_value ( fmin, optpara, n )

  use BASIC_VARIABLES, only: ncomp, np, outp, gibbs, xif, xi, lnx, my_f
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(IN)                    :: n
  real, intent(IN)                       :: optpara(:)
  real, intent(IN OUT)                   :: fmin

  !-----------------------------------------------------------------------------
  integer                                :: i
  real                                   :: lnphi(np,nc)
  real                                   :: ph_1_frac, punish
  real, dimension(nc)                    :: ni_1, ni_2
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! saveguard that the entering parameters are reasonable
  !-----------------------------------------------------------------------------
  if ( maxval( optpara(1:n) ) > 5.0 ) then
     fmin = 1.E8 * maxval( optpara(1:n) )
     return
  end if
  if ( maxval( optpara(1:n) ) < -100.0 ) then
     fmin = 1.E8
     return
  end if
  do i = 1, n
     if ( optpara(i) /= optpara(i) ) then
        fmin = 1.E8
        return
     end if
  end do

  !-----------------------------------------------------------------------------
  ! setting of mole numbers of all species in second phase
  !-----------------------------------------------------------------------------

  do i = 1, ncomp
     if ( optpara(i) < -100.0 ) then
        ni_2(i) = 0.0
     else
        ni_2(i) = EXP( optpara(i) )
     end if
  end do

  punish = 0.0
  do i = 1, ncomp
     ni_1(i) = xif(i) - ni_2(i)
     if ( ni_2(i) > xif(i) ) then
        punish = punish - 1.E8*ni_1(i) + 0.1
        !ni_2(i) = xif(i)
        ni_1(i) = xif(i) * 1.E-50
     end if
     if ( ni_1(i) <= 0.0 ) ni_1(i) = xif(i) * 1.E-50
  end do


  !-----------------------------------------------------------------------------
  ! calculate mole fractions of second and first phase (from component balance)
  !-----------------------------------------------------------------------------

  xi(2,1:ncomp) = ni_2(1:ncomp) / SUM( ni_2(1:ncomp) )
  lnx(2,1:ncomp) = optpara(1:ncomp) - LOG( SUM( ni_2(1:ncomp) ) )

  ph_1_frac = SUM( ni_1(1:ncomp) )
  xi( 1, 1:ncomp ) = ni_1( 1:ncomp ) / ph_1_frac
  lnx( 1, 1:ncomp ) = LOG( ni_1( 1:ncomp ) ) - LOG( ph_1_frac )

  !-----------------------------------------------------------------------------
  ! calculate total gibbs energy (for given T.p, and x(2,:) of trial phase)
  !-----------------------------------------------------------------------------

  ! write (*,*) 'f_TP x', xi( 1, 1:ncomp )
  ! write (*,*) 'f_TP x', xi( 2, 1:ncomp )
  CALL fugacity (lnphi)

  fmin = gibbs(1) * ph_1_frac + gibbs(2) * ( 1.0 - ph_1_frac )  + punish
  ! write (*,*) 'f_TP f=',gibbs(1) , gibbs(2)
  fmin = fmin - sum( xif(1:ncomp)*my_f( 1:ncomp ))
  if ( outp >= 2 ) write (*,'(a,4G18.8)') 'TP',fmin, ni_1(1:ncomp) !,xi(1,1:ncomp)
  ! write (*,'(a,4G18.8)') 'al', ph_1_frac, (xi(2,i), i=1,ncomp)
  ! write (*,*) ' '
  ! pause

end subroutine tangent_value



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_grad ( g, optpara, n )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: KBOL
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(IN)                    :: n
  real, intent(IN)                       :: optpara(:)
  real, intent(in out)                   :: g(:)

  !-----------------------------------------------------------------------------
  integer                                :: i
  real                                   :: lnphi(np,nc), ph_1_frac
  real, dimension(nc)                    :: ni_1, ni_2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! saveguard that the entering parameters are reasonable
  !-----------------------------------------------------------------------------
  if ( maxval( optpara(1:n) ) < -100.0 ) then
     g(:) = 0.0
     return
  end if
  do i = 1, n
     if ( optpara(i) /= optpara(i) .OR. optpara(i) > 5.0 ) then
        g(:) = 0.0
        return
     end if
  end do

  !-----------------------------------------------------------------------------
  ! setting of mole numbers of all species in second phase
  !-----------------------------------------------------------------------------

  do i = 1, n
     if ( optpara(i) < -300.0 ) then
        ni_2(i) = 0.0
     else
        ni_2(i) = EXP( optpara(i) )
     end if
  end do

  do i = 1, ncomp
     ni_1(i) = xif(i) - ni_2(i)
     if ( ni_2(i) > xif(i) ) then
        !write (*,*) 'now active bound ================================'
        !ni_2(i) = xif(i)
        ni_1(i) = xif(i) * 1.E-100
     end if
  end do

  !-----------------------------------------------------------------------------
  ! calculate mole fractions of second and first phase (from component balance)
  !-----------------------------------------------------------------------------

  xi(2,1:ncomp) = ni_2(1:ncomp) / SUM( ni_2(1:ncomp) )
  lnx(2,1:ncomp) = optpara(1:ncomp) - LOG( SUM( ni_2(1:ncomp) ) )

  ph_1_frac = SUM( ni_1(1:ncomp) )
  xi( 1, 1:ncomp ) = ni_1( 1:ncomp ) / ph_1_frac
  lnx( 1, 1:ncomp ) = LOG( ni_1( 1:ncomp ) ) - LOG( ph_1_frac )

  !-----------------------------------------------------------------------------
  ! calculate chemical potentials (for given T.p, and x(2,:) of trial phase)
  !-----------------------------------------------------------------------------

  CALL fugacity (lnphi)

  !g( 1:ncomp ) =  - ( lnphi(1,1:ncomp) + log( rhoi_cal(1,1:ncomp) )   &
  !                   + LOG(p/(kbol*t*1.E30*sum(rhoi_cal(1,1:ncomp))))  )  &
  !                + ( lnphi(2,1:ncomp) + log( rhoi_cal(2,1:ncomp) )   &
  !                + LOG(p/(kbol*t*1.E30*sum(rhoi_cal(2,1:ncomp))))   )
  !g( 1:ncomp ) = g( 1:ncomp ) * EXP( optpara(1:ncomp) )

  g( 1:ncomp ) =  - ( lnphi(1,1:ncomp) + log( xi(1,1:ncomp) )  )  &
                  + ( lnphi(2,1:ncomp) + log( xi(2,1:ncomp) )  )
  g( 1:ncomp ) = g( 1:ncomp ) * EXP( optpara(1:ncomp) )

  ! write (*,'(a,4G18.8)') 'g', g(1:ncomp)

end subroutine tangent_grad

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_hessian (hessian, gtrans, optpara, n)

  use BASIC_VARIABLES
  implicit none

  integer, intent(IN)        :: n
  real, intent(IN)           :: optpara(:)
  real, intent(IN OUT)       :: gtrans(:)
  real, intent(IN OUT)       :: hessian(:,:)

  !-----------------------------------------------------------------------------
  integer                                 :: i, j
  real                                    :: delta
  real                                    :: optpara_mod(n), g(n), gi(n,n), gi_left(n,n)
  !-----------------------------------------------------------------------------


  delta = 1.0E-4

  optpara_mod = optpara

  do i = 1, n

     optpara_mod = optpara
     optpara_mod(i) = optpara(i)*(1.0+delta)

     call tangent_grad (g, optpara_mod, n)
     gi(i,1:n) = g(1:n)

     optpara_mod(i) = optpara(i)*(1.0-delta)

     call tangent_grad (g, optpara_mod, n)
     gi_left(i,1:n) = g(1:n)

  end do

  call tangent_grad (g, optpara, n)


  do i = 1, n
     do j = 1, n

        hessian(i,j) = ( gi(i,j) - gi_left(i,j) ) / ( 2.0*optpara(i)*delta )
        !  hessian(j,i) = hessian(i,j)
        !  write (*,*) i,j,hessian(i,j)

     end do
  end do

  gtrans = g
  ! call tangent_value (fmin, optpara, n)
  ! write (*,'(4G20.12)') g(1:n),xi(1,1),xi(2,1)

end subroutine tangent_hessian




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!  input:   xi(1,:), xi(2,:), densta(1:2)
!
!  output:  rhoi1, rhoi2, xi(1,:), xi(2,:), dense(:)
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_plane_2 ( phi_2_start, rhoi1, rhoi2 )

  use BASIC_VARIABLES
  use cg_minimization
  implicit none

  !-----------------------------------------------------------------------------
  real, intent(in)                       :: phi_2_start
  real, dimension(nc), intent(out)       :: rhoi1
  real, dimension(nc), intent(out)       :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                :: n
  integer                                :: small_i, min_ph, other_ph
  real                                   :: lnphi(np,nc)
  real, allocatable                      :: optpara(:)
  integer                                :: STATUS,  iter, nfunc, ngrad
  real                                   :: gnorm, fmin
  real, allocatable                      :: d(:), g(:), xtemp(:), gtemp(:)
  !-----------------------------------------------------------------------------

  n = ncomp

  ALLOCATE( optpara(n) )
  ALLOCATE( d(n) )
  ALLOCATE( g(n) )
  ALLOCATE( xtemp(n) )
  ALLOCATE( gtemp(n) )

  lnx( 1, 1:ncomp ) = LOG( xi( 1, 1:ncomp ) )
  lnx( 2, 1:ncomp ) = LOG( xi( 2, 1:ncomp ) )

  optpara( 1:ncomp ) = LOG( xi( 2, 1:ncomp ) * phi_2_start )

  !if ( n == 2 ) then
  !   CALL Newton_Opt_2D ( tangent_hessian, tangent_value, optpara, n, 1.E-8, 1.E-8, gtemp, fmin )
  !else
     CALL cg_descent ( 1.E-6, optpara, n, tangent_value, tangent_grad, STATUS, &
                       gnorm, fmin,  iter, nfunc, ngrad, d, g, xtemp, gtemp )
  !end if

  !-----------------------------------------------------------------------------
  ! If one component is a polymer (indicated by a low component-density)
  ! then get an estimate of the polymer-lean composition, by solving for
  ! xi_p1 = ( xi_p2 * phii_p2) / phii_p1     (phase equilibrium condition,
  ! with p1 for phase 1)
  !-----------------------------------------------------------------------------
  if ( MINVAL( lnx( 1, 1:ncomp ) ) < MINVAL( lnx( 2, 1:ncomp ) ) ) then
     min_ph   = 1
     other_ph = 2
  else
     min_ph   = 2
     other_ph = 1
  end if
  small_i = MINLOC( lnx(min_ph,1:ncomp), 1 )
  !---- if one component is a polymer ------------------------------------------
  if ( MINVAL( lnx(min_ph,1:ncomp) ) < -20.0 ) then
     CALL FUGACITY ( lnphi )
     lnx(min_ph,small_i) = lnx(other_ph,small_i)+lnphi(other_ph,small_i) - lnphi(min_ph,small_i)
     xi(min_ph,1:ncomp) = EXP( lnx(min_ph,1:ncomp) ) / SUM( EXP( lnx(min_ph,1:ncomp) ) )
     call FUGACITY ( lnphi )
  end if

  rhoi1( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )
  rhoi2( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )

  DEALLOCATE( optpara, d, g, xtemp, gtemp )


end subroutine tangent_plane_2




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine tangent_plane ( phi_2_start, rhoi1, rhoi2 )

  use BASIC_VARIABLES
  implicit none

  !-----------------------------------------------------------------------------
  real, intent(in)                       :: phi_2_start
  real, dimension(nc), intent(out)       :: rhoi1
  real, dimension(nc), intent(out)       :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                :: n
  integer                                :: small_i, min_ph, other_ph
  integer                                :: PRIN
  real                                   :: fmin , t0, h0, MACHEP
  real                                   :: lnphi(np,nc)
  real, allocatable                      :: optpara(:)
  !-----------------------------------------------------------------------------

  n = ncomp
  t0 = 1.E-4
  h0 = 0.1
  PRIN = 0
  MACHEP = 1.E-15

  ALLOCATE( optpara(n) )

  lnx( 1, 1:ncomp ) = LOG( xi( 1, 1:ncomp ) )
  lnx( 2, 1:ncomp ) = LOG( xi( 2, 1:ncomp ) )

  optpara( 1:ncomp ) = LOG( xi( 2, 1:ncomp ) * phi_2_start )


  call PRAXIS( t0, MACHEP, h0, n, PRIN, optpara, tangent_value, fmin )

  ! The optimal optpara-vector is not necessarily the one that was last evaluated.
  ! tangent_value is reexecuted with the optimal vector optpara, in order to update the ln(x) values
  call tangent_value( fmin, optpara, n )


  !-----------------------------------------------------------------------------
  ! If one component is a polymer (indicated by a low component-density)
  ! then get an estimate of the polymer-lean composition, by solving for
  ! xi_p1 = ( xi_p2 * phii_p2) / phii_p1     (phase equilibrium condition,
  ! with p1 for phase 1)
  !-----------------------------------------------------------------------------
  if ( MINVAL( lnx( 1, 1:ncomp ) ) < MINVAL( lnx( 2, 1:ncomp ) ) ) then
     min_ph   = 1
     other_ph = 2
  else
     min_ph   = 2
     other_ph = 1
  end if
  small_i = MINLOC( lnx(min_ph,1:ncomp), 1 )
  !---- if one component is a polymer ------------------------------------------
  if ( MINVAL( lnx(min_ph,1:ncomp) ) < -20.0 ) then
     CALL FUGACITY ( lnphi )
     lnx(min_ph,small_i) = lnx(other_ph,small_i)+lnphi(other_ph,small_i) - lnphi(min_ph,small_i)
     optpara(small_i) = lnx(2,small_i) + LOG( SUM( EXP( optpara(1:ncomp) ) ) )
  end if

  !-----------------------------------------------------------------------------
  ! caution: these initial values are for a flashcase overwritten in
  ! SUBROUTINE determine_flash_it2, because in that case, the lnx-values
  ! treated as ln(mole_number).
  !-----------------------------------------------------------------------------
  rhoi1( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )
  rhoi2( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )

  DEALLOCATE( optpara )

end subroutine tangent_plane



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine Helmholtz_flash_grads ( n, x_in, f_tpd, grad, hessian, diagonal )

  use PARAMETERS, only: kbol, pi
  use BASIC_VARIABLES, only: np, nc, ncomp, t, p, xi, xiF, densta, my_cal, parame
  use EOS_VARIABLES, only: dhs
  use utilities, only: fden_calc
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real, dimension(n), intent(in)         :: x_in
  real, intent(out)                      :: f_tpd
  real, dimension(n), intent(out)        :: grad
  real, dimension(n,n), intent(out)      :: hessian
  real, dimension(n,n), intent(out)      :: diagonal

  !-----------------------------------------------------------------------------
  integer                                :: i, j
  real                                   :: phi
  real                                   :: V
  real, dimension(ncomp)                 :: w
  real, dimension(ncomp)                 :: rhoi1, rhoi2
  real                                   :: f1, f2
  real, dimension(ncomp)                 :: my1, my2
  real                                   :: lnphi(np,nc)
  real                                   :: df_phi, df_V
  real, dimension(ncomp)                 :: df_wi
  real, allocatable, dimension(:,:)      :: A1_rr, Aig1_rr, A2_rr, Aig2_rr
  real, dimension(ncomp,ncomp)           :: df_wi_wj
  real, dimension(ncomp)                 :: df_wi_V, df_wi_phi
  real                                   :: df_phi_phi, df_phi_V, df_V_V
  real, dimension(2)                     :: sum_f_rhoi_rhoj_w_w, sum_f_rhoi_rhoj_w_r, sum_f_rhoi_rhoj_r_r
  !-----------------------------------------------------------------------------

  w(1:ncomp) = x_in(1:ncomp)
  phi    = x_in(ncomp+1)
  V      = x_in(ncomp+2)

  ! V1 = V * ( 0.5 + phi )
  ! V2 = V * ( 0.5 - phi )

  do i = 1, ncomp
     rhoi1(i) = xiF(i) / V + ( phi - 0.5 ) * w(i)
     rhoi2(i) = xiF(i) / V + ( phi + 0.5 ) * w(i)
     if ( rhoi1(i) < 0.0 ) rhoi1(i) = 0.0
     if ( rhoi2(i) < 0.0 ) rhoi2(i) = 0.0
  end do

  call fden_calc (f1, rhoi1)
  call fden_calc (f2, rhoi2)

  densta(1) = PI / 6.0 * sum( rhoi1(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
  densta(2) = PI / 6.0 * sum( rhoi2(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
  xi(1,1:ncomp) = rhoi1(1:ncomp) / sum( rhoi1(1:ncomp) )
  xi(2,1:ncomp) = rhoi2(1:ncomp) / sum( rhoi2(1:ncomp) )
  call FUGACITY ( lnphi )

  my1( 1:ncomp ) = my_cal(1,1:ncomp) ! lnphi(1,1:ncomp) + log( xi(1,1:ncomp) )
  my2( 1:ncomp ) = my_cal(2,1:ncomp) ! lnphi(2,1:ncomp) + log( xi(2,1:ncomp) )

  ! f_tpd = V1 * f1  +  V2 * f2  +  p / ( KBOL*1.E30 * t ) * V
  f_tpd = V * ( ( 0.5 + phi ) * f1  +  ( 0.5 - phi ) * f2  +  p / ( KBOL*1.E30 * t ) )
  ! write (*,*) 'f_tpd',f_tpd

  df_phi = V * ( f1 -f2 + ( 0.5 + phi )*sum( my1(1:ncomp)*w(1:ncomp) )  &
                        + ( 0.5 - phi )*sum( my2(1:ncomp)*w(1:ncomp) ) )
  df_V =  ( f_tpd - ( 0.5 + phi )*sum( my1(1:ncomp)*xiF(1:ncomp) )  &
                  - ( 0.5 - phi )*sum( my2(1:ncomp)*xiF(1:ncomp) ) ) / V
  do i = 1, ncomp
     df_wi(i) =  V * ( phi*phi - 0.25 ) * ( my1(i) - my2(i) )
  end do

  grad(1:ncomp) = df_wi(1:ncomp)
  grad(ncomp+1) = df_phi
  grad(ncomp+2) = df_V

  allocate( A1_rr( ncomp, ncomp ), Aig1_rr( ncomp, ncomp ) )
  allocate( A2_rr( ncomp, ncomp ), Aig2_rr( ncomp, ncomp ) )

  call ddA_drhoi_drhoj_EOS ( ncomp, rhoi1, A1_rr, Aig1_rr )
  A1_rr(:,:) = A1_rr(:,:) + Aig1_rr(:,:)

  sum_f_rhoi_rhoj_w_w(1) = 0.0
  sum_f_rhoi_rhoj_w_r(1) = 0.0
  sum_f_rhoi_rhoj_r_r(1) = 0.0
  do i = 1, ncomp
     do j = 1, ncomp
        sum_f_rhoi_rhoj_w_w(1) = sum_f_rhoi_rhoj_w_w(1) + A1_rr(i,j)*w(i)*w(j)
        sum_f_rhoi_rhoj_w_r(1) = sum_f_rhoi_rhoj_w_r(1) + A1_rr(i,j)*w(i)*xiF(j)
        sum_f_rhoi_rhoj_r_r(1) = sum_f_rhoi_rhoj_r_r(1) + A1_rr(i,j)*xiF(i)*xiF(j)
     end do
  end do

  call ddA_drhoi_drhoj_EOS ( ncomp, rhoi2, A2_rr, Aig2_rr )
  A2_rr(:,:) = A2_rr(:,:) + Aig2_rr(:,:)


  sum_f_rhoi_rhoj_w_w(2) = 0.0
  sum_f_rhoi_rhoj_w_r(2) = 0.0
  sum_f_rhoi_rhoj_r_r(2) = 0.0
  do i = 1, ncomp
     do j = 1, ncomp
        sum_f_rhoi_rhoj_w_w(2) = sum_f_rhoi_rhoj_w_w(2) + A2_rr(i,j)*w(i)*w(j)
        sum_f_rhoi_rhoj_w_r(2) = sum_f_rhoi_rhoj_w_r(2) + A2_rr(i,j)*w(i)*xiF(j)
        sum_f_rhoi_rhoj_r_r(2) = sum_f_rhoi_rhoj_r_r(2) + A2_rr(i,j)*xiF(i)*xiF(j)
     end do
  end do


  do i = 1, ncomp
     do j = 1, ncomp
        df_wi_wj(i,j) = V * (phi*phi - 0.25 ) * ( (phi-0.5) * A1_rr(i,j) - (phi+0.5) * A2_rr(i,j) )
     end do
  end do

  do i = 1, ncomp
     df_wi_phi(i) = V * ( 2.0*phi*( my1(i)-my2(i) ) &
                           + (phi*phi-0.25)*sum( w(1:ncomp)*(A1_rr(i,1:ncomp)-A2_rr(i,1:ncomp)) ) )
     df_wi_V(i) = ( df_wi(i) - (phi*phi-0.25)*sum( xiF(1:ncomp)*(A1_rr(i,1:ncomp)-A2_rr(i,1:ncomp)) ) ) / V
  end do

  df_phi_phi = V * ( 2.0*sum( ( my1(1:ncomp)-my2(1:ncomp))*w(1:ncomp) )  &
                        + ( 0.5 + phi ) * sum_f_rhoi_rhoj_w_w(1)  &
                        + ( 0.5 - phi ) * sum_f_rhoi_rhoj_w_w(2) )
  df_phi_V = ( df_phi - sum( ( my1(1:ncomp)-my2(1:ncomp))*xiF(1:ncomp) )  &
                        - ( 0.5 + phi ) * sum_f_rhoi_rhoj_w_r(1)  &
                        - ( 0.5 - phi ) * sum_f_rhoi_rhoj_w_r(2) ) / V
  df_V_V = ( (0.5+phi) * sum_f_rhoi_rhoj_r_r(1) + (0.5-phi) * sum_f_rhoi_rhoj_r_r(2) ) / V**3

  do i = 1, ncomp
     do j = 1, ncomp
        hessian(i,j) = df_wi_wj(i,j)
     end do
     hessian(i,ncomp+1) = df_wi_phi(i)
     hessian(i,ncomp+2) = df_wi_V(i)
     hessian(ncomp+1,i) = df_wi_phi(i)
     hessian(ncomp+2,i) = df_wi_V(i)
  end do
  hessian(ncomp+1,ncomp+1) = df_phi_phi
  hessian(ncomp+1,ncomp+2) = df_phi_V
  hessian(ncomp+2,ncomp+1) = df_phi_V
  hessian(ncomp+2,ncomp+2) = df_V_V

  do i = 1, ncomp
        diagonal(i,i) = V * (phi*phi - 0.25 ) * ( (phi-0.5)/rhoi1(i) - (phi+0.5)/rhoi2(i) )
  end do
  diagonal(ncomp+1,ncomp+1) = V * ( 2.0*sum( ( log( rhoi1(1:ncomp) )- log(rhoi2(1:ncomp)) )*w(1:ncomp) )  &
                              + (0.5+phi) * sum(w(1:ncomp)*w(1:ncomp)/rhoi1(1:ncomp))  &
                              + (0.5-phi) * sum(w(1:ncomp)*w(1:ncomp)/rhoi2(1:ncomp)) )
  diagonal(ncomp+2,ncomp+2) = ( (0.5+phi) * sum(xiF(1:ncomp)*xiF(1:ncomp)/rhoi1(1:ncomp))  &
                              + (0.5-phi) * sum(xiF(1:ncomp)*xiF(1:ncomp)/rhoi2(1:ncomp)) ) / V**3

  !write (*,*) ' '
  !write (*,*) 'G1_i:',df_wi(1:ncomp)
  !write (*,*) 'G2  :',df_phi
  !write (*,*) 'G3  :',df_V

  deallocate( A1_rr, Aig1_rr, A2_rr, Aig2_rr )

end subroutine Helmholtz_flash_grads



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine phase_stability
!
! input: my_f (not in argument list) as the chem. potential of the feed-phase.
!        rhoi_feed (argument list) density vector of feed phase
!        eta_trial (argument list) as the starting density for trial phase.
!
! output is: ph_split=0 for stable phase or ph_split=1 for a phase split. The
! array rhoi_trial( 1:ncomp ) is then an output, defining the minimum of the
! tangent plane distance.
! The value ph_split=2 is assigned for cases, where the TPD is not negative, but
! very close to zero and where the compositions are different to the feed-phase.
! This is encountered close to the critical point.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine phase_stability ( rhoi_feed, eta_trial, ph_split, rhoi_trial )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: dhs, PI, mseg
  use EOS_NUMERICAL_DERIVATIVES
  use cg_minimization
  use optimizer_2D
  implicit none

  !-----------------------------------------------------------------------------
  real, dimension(nc), intent(IN)            :: rhoi_feed
  real, intent(IN)                           :: eta_trial
  integer, intent(OUT)                       :: ph_split
  real, dimension(nc), intent(OUT)           :: rhoi_trial
  !-----------------------------------------------------------------------------

  integer                                    :: n
  real                                       :: fmin
  real, allocatable                          :: optpara(:)

  integer                                    :: i
  integer                                    :: i_trial_composition
  real                                       :: rhoi(nc)
  real                                       :: delta_rhoi(nc)
  real                                       :: w(np,nc), mean_mass
  integer                                    :: ph_split_i( 2*ncomp+1 )
  real                                       :: fmin_vec( 2*ncomp+1 )
  real                                       :: rhoi2_vec(2*ncomp+1, 1:ncomp)
  character (LEN=2)                          :: ensemble_flag_save

  integer                                    :: STATUS,  iter, nfunc, ngrad
  real                                       :: gnorm
  real, allocatable                          :: d(:), g(:), xtemp(:), gtemp(:)
  !-----------------------------------------------------------------------------


  ph_split = 0
  ensemble_flag_save = ensemble_flag
  ensemble_flag = 'tv'

  n = ncomp
  ALLOCATE( optpara(n) )
  ALLOCATE( d(n) )
  ALLOCATE( g(n) )
  ALLOCATE( xtemp(n) )
  ALLOCATE( gtemp(n) )


  if ( outp >= 2 ) then
     write (*,'(a)') ' reference x_i, tested for stability:'
     write (*,'(a,4G21.12)') ' x_i feed:', rhoi_feed(1:ncomp) / sum( rhoi_feed(1:ncomp) )
     write (*,*) ' '
  end if


  do i_trial_composition = 0, (ncomp + ncomp)

     !--------------------------------------------------------------------------
     ! setting trial-phase mole-fractions
     !--------------------------------------------------------------------------
     if ( i_trial_composition == 0 ) w(2,1:ncomp) = 1.0 / REAL(ncomp)
     if ( i_trial_composition >= 1 .AND. i_trial_composition <= ncomp ) then
        w( 2, 1:ncomp ) = 1.0 / REAL( ncomp-1 ) * 0.05
        w( 2, i_trial_composition ) = 0.95
     else if ( i_trial_composition >= ncomp+1 ) then
        w( 2, 1:ncomp ) = 1.0 / REAL( ncomp-1 ) * 0.00001
        w( 2, i_trial_composition - ncomp ) = 0.99999
     end if

     mean_mass = 1.0 /  SUM( w(2,1:ncomp)/mm(1:ncomp) )
     xi(2,1:ncomp) = w(2,1:ncomp)/mm(1:ncomp) * mean_mass
     if ( outp >= 2 ) write (*,'(a,4G20.12)') ' trial phase x  ---- ', xi(2,1:ncomp)

     !--------------------------------------------------------------------------
     ! starting values for iteration (optpara)
     !--------------------------------------------------------------------------
     do i = 1,ncomp
        rhoi(i) = xi(2,i)*eta_trial/sum( PI/6.0*xi(2,1:ncomp)*mseg(1:ncomp)*dhs(1:ncomp)**3 )
        optpara(i) = LOG( rhoi(i) )
        !optpara(i) = rhoi(i)
     end do

     !--------------------------------------------------------------------------
     ! minimizing the objective fct. Phase split for values of fmin < 0.0
     !--------------------------------------------------------------------------

     if ( n == 2 ) then
        CALL Newton_Opt_2D ( f_stability, optpara, n, 1.E-8, 1.E-8, gtemp, fmin )
        fmin_best = fmin
        do i = 1, ncomp
        if ( optpara(i) > -50.0 ) then
           rhoi_best(i) = exp( optpara(i) )
        else
           rhoi_best(i) = 1.E-50
        end if
        end do
     else
        fmin_best = 1.E30
        CALL cg_descent ( 1.E-5, optpara, n, f_stability, stability_grad, STATUS, &
                          gnorm, fmin,  iter, nfunc, ngrad, d, g, xtemp, gtemp)
        fmin = fmin_best
     end if

     !--------------------------------------------------------------------------
     ! save result to vector and determine stability of resulting phase
     !--------------------------------------------------------------------------

     delta_rhoi( 1:n ) = ABS( 1.0 - rhoi_best(1:n) / rhoi_feed(1:n) )

     rhoi2_vec( i_trial_composition + 1, 1:ncomp ) = rhoi_best(1:ncomp)
     fmin_vec( i_trial_composition + 1 ) = 1.E8

     ph_split_i( i_trial_composition + 1 ) = 0
     if ( fmin < 1.E-3 .AND. maxval( delta_rhoi(1:n) ) > 0.1 ) then
        ph_split_i( i_trial_composition + 1 ) = 2
        fmin_vec( i_trial_composition + 1 ) = fmin
     end if

     if (fmin < -1.E-8 .AND. maxval( delta_rhoi(1:n) ) > 0.001 ) then
        ph_split_i( i_trial_composition + 1 ) = 1
        fmin_vec( i_trial_composition + 1 ) = fmin
     end if


     !--------------------------------------------------------------------------
     ! output to terminal
     !--------------------------------------------------------------------------
     if ( outp >= 2 ) then
        write (*,'(a, 5G20.12)') ' eta, xi_trialphase  ',PI/6.0  &
               * SUM( rhoi_best(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 ),  &
                 rhoi_best(1:ncomp) / sum( rhoi_best(1:ncomp) )
        write (*,'(a,G20.12,i4,G20.12,i3)') ' trial-result: fmin, N, d(x), ph-split:',  &
                 fmin, i_trial_composition, maxval( delta_rhoi( 1:n ) ),  &
                 ph_split_i( i_trial_composition+1 )
        write (*,*) ' '
     end if


  end do

  !-----------------------------------------------------------------------------
  ! write vector (rhoi_trial). Vector rhoi_trial gives the (overall) minimum of the TPD
  !-----------------------------------------------------------------------------
  i_trial_composition = minloc( fmin_vec( 1: (ncomp+ncomp)+1 ), 1 )
  rhoi_trial( 1:ncomp ) = rhoi2_vec( i_trial_composition, 1:ncomp )
  ph_split = ph_split_i( i_trial_composition )

  if ( outp >= 2 ) write (*,*) 'selected trial-result',i_trial_composition,fmin_vec( i_trial_composition )

  DEALLOCATE( optpara )
  DEALLOCATE( d, g, xtemp, gtemp )
  ensemble_flag = ensemble_flag_save

end subroutine phase_stability



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
subroutine f_stability ( f_tpd, optpara, n )

  use PARAMETERS, only: PI, KBOL
  use BASIC_VARIABLES
  use EOS_VARIABLES, only: dhs
  use utilities, only: fden_calc
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real, intent(in)                       :: optpara(:)
  real, intent(in out)                   :: f_tpd

  !-----------------------------------------------------------------------------
  integer                                :: i
  real                                   :: rhoi(nc),gradterm
  real                                   :: fden,punish
  real                                   :: dens
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! saveguard that the entering parameters are reasonable
  !-----------------------------------------------------------------------------
  if ( maxval( optpara(1:n) ) > 2.0 ) then
     f_tpd = 100.0 * maxval( optpara(1:n) )
     return
  end if
  do i = 1, n
     if ( optpara(i) /= optpara(i) ) then
        f_tpd = 1.E8
        return
     end if
  end do
  if ( maxval( optpara(1:n) ) < -50.0 ) then
     f_tpd = - 100.0 * sum( optpara(1:n) )
     return
  end if

  !-----------------------------------------------------------------------------
  ! assign density of each component
  !-----------------------------------------------------------------------------

  punish = 0.0

  do i = 1, n
     if ( optpara(i) > -50.0 ) then
        rhoi(i) = exp( optpara(i) )
     else
        rhoi(i) = 1.E-50
     end if
     !if ( optpara(i) < 0.0) then
     !   rhoi(i) = 1.E-50 / (1.0-optpara(i)*1.E10)
     !   punish = punish - optpara(i)*1.E6
     !end if
     !if ( optpara(i) >= 0.6) rhoi(i) = 0.6
  end do

  dens = PI/6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )

  !-----------------------------------------------------------------------------
  ! correcting too high densities
  !-----------------------------------------------------------------------------

  if (dens > 0.65) then
     punish = punish + (dens-0.65)*1000.0
     rhoi(1:n) = rhoi(1:n)*0.65/dens
  end if

  !-----------------------------------------------------------------------------
  ! calculate the Helmholtz energ and the tangent plane distance (TPD)
  !-----------------------------------------------------------------------------

  call fden_calc (fden, rhoi)

  gradterm = sum( my_f(1:n) * rhoi(1:n) )

  f_tpd = fden + (p * 1.E-30) / (KBOL*t) - gradterm + punish
  f_tpd = f_tpd * 1000.0
  ! write (*,'(a,5G16.8)') 'f_tpd', f_tpd, optpara(1:n), punish

  if ( f_tpd < fmin_best ) then
     fmin_best = f_tpd
     rhoi_best(1:n) = rhoi(1:n)
  end if

end subroutine f_stability



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_grad (g, optpara, n)

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: dhs, PI
  implicit none

  integer, intent(in)        :: n
  real, intent(in)           :: optpara(:)
  real, intent(in out)       :: g(:)

  !-----------------------------------------------------------------------------
  integer                                 :: i
  integer                                 :: nphas_save
  real                                    :: rhoi(nc)
  real                                    :: lnphi(np,nc)
  !-----------------------------------------------------------------------------

  nphas_save = nphas
  nphas = 1

  !-----------------------------------------------------------------------------
  ! saveguard that the entering parameters are reasonable
  !-----------------------------------------------------------------------------

  do i = 1, n
     if ( optpara(i) /= optpara(i) .OR. optpara(i) > 5.0 ) then
        g(:) = 0.0
        nphas = nphas_save
        return
     end if
  end do
  if ( maxval( optpara(1:n) ) < -50.0 ) then
     g(:) = 0.0
     nphas = nphas_save
     return
  end if

  !-----------------------------------------------------------------------------
  ! assign density of all species
  !-----------------------------------------------------------------------------

  do i = 1, n
     if ( optpara(i) > -50.0 ) then
        rhoi(i) = exp( optpara(i) )
     else
        rhoi(i) = 1.E-50
     end if
     !if ( optpara(i) < 0.0) rhoi(i) = 1.E-50 / (1.0-optpara(i)*1.E10)
     !if ( optpara(i) >= 0.6) rhoi(i) = 0.6
  end do

  densta(1) = PI/6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )

  !-----------------------------------------------------------------------------
  ! correcting too high densities
  !-----------------------------------------------------------------------------

  if (densta(1) > 0.65) then
     rhoi(1:n) = rhoi(1:n) * 0.65 / densta(1)
     densta(1) = PI/6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
  end if

  xi(1,1:ncomp) = rhoi(1:ncomp) / SUM( rhoi(1:ncomp) )

  !-----------------------------------------------------------------------------
  ! calculate the chemical potential of all species
  !-----------------------------------------------------------------------------

  call fugacity (lnphi)

  g( 1:n ) = lnphi( 1, 1:n ) + LOG( rhoi( 1:n ) )   -  my_f( 1:n )
  g( 1:n ) = g( 1:n ) * rhoi( 1:n ) * 1000.0
  ! write (*,'(a,5G16.8)') 'f_tpd_grad',g(1:n)

  nphas = nphas_save

end subroutine stability_grad




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine stability_hessian (hessian, gtrans, fmin, optpara, n)

  use BASIC_VARIABLES
  implicit none

  integer, intent(IN)        :: n
  real, intent(IN)           :: optpara(:)
  real, intent(IN OUT)       :: fmin
  real, intent(IN OUT)       :: gtrans(:)
  real, intent(IN OUT)       :: hessian(:,:)

  !-----------------------------------------------------------------------------
  integer                                 :: i, j
  real                                    :: delta
  real                                    :: optpara_mod(n), g(n), gi_right(n,n), gi_left(n,n)
  real, allocatable, dimension(:,:)       :: A_rr, Aig_rr
  !-----------------------------------------------------------------------------


  delta = 1.0E-4

  optpara_mod = optpara

  do i = 1, n

     optpara_mod = optpara
     optpara_mod(i) = optpara(i)*(1.0+delta)

     call stability_grad (g, optpara_mod, n)
     gi_right(i,1:n) = g(1:n)

     optpara_mod(i) = optpara(i)*(1.0-delta)

     call stability_grad (g, optpara_mod, n)
     gi_left(i,1:n) = g(1:n)

  end do

  call stability_grad (g, optpara, n)

  do i = 1, n
     do j = 1, n

        hessian(i,j) = ( gi_right(i,j) - gi_left(i,j) ) / ( 2.0*optpara(i)*delta )
        !  hessian(j,i) = hessian(i,j)
          write (*,*) i,j,hessian(i,j)

     end do
  end do
  write (*,*) ' '

  gtrans = g
  call f_stability (fmin, optpara, n)

  allocate( A_rr(ncomp,ncomp), Aig_rr(ncomp,ncomp) )
  call ddA_drhoi_drhoj_EOS ( ncomp, exp( optpara(:) ), A_rr, Aig_rr )
  do i = 1, n
     do j = 1, n
          hessian(i,j) = A_rr(i,j)
          if (i==j) hessian(i,j) = hessian(i,j) + 1.0 / exp( optpara(i) )
          !--- convert to second derivative ddA / dln(rho_k) dln(rho_l) ( from ... d(rho_l) )
          if (i/=j) hessian(i,j) = hessian(i,j)*exp( optpara(i) )*exp( optpara(j) )*1000.0
          if (i==j) hessian(i,j) = hessian(i,j)*exp( optpara(i) )*exp( optpara(j) )*1000.0 + g(i)
          if (i/=j) write (*,*) i,j,hessian(i,j)
          if (i==j) write (*,*) i,j,hessian(i,j)
     end do
  end do
  deallocate( A_rr, Aig_rr )
  read (*,*)

end subroutine stability_hessian



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine vle_min

  use PARAMETERS, only: RGAS
  use BASIC_VARIABLES
  use utilities
  implicit none

  !-----------------------------------------------------------------------------
  integer, parameter                     :: steps = 40
  logical                                :: lle_check  ! currently not any more used
  integer                                :: i, j, k, phasen(0:steps)
  real                                   :: lnphi(np,nc)
  real, dimension(0:steps)               :: vlemin, llemin, xval, start_xv, start_xl
  real                                   :: x_sav,dg_dx2
  !-----------------------------------------------------------------------------

  j = 0
  k = 0
  nphas = 2

  start_xv(:) = 0.0
  start_xl(:) = 0.0

  x_sav = xi(1,1)
  sum_rel(1) = 'x12' ! summation relation
  sum_rel(2) = 'x22' ! summation relation

  do i = 0, steps

     densta(1) = 0.45
     densta(2) = 1.E-6
     xi(1,1) = 1.0 - REAL(i) / REAL(steps)
     if ( xi(1,1) <= 1.E-50 ) xi(1,1) = 1.E-50
     xi(2,1)  = xi(1,1)
     lnx(1,1) = LOG(xi(1,1))
     lnx(2,1) = LOG(xi(2,1))

     CALL x_summation
     CALL fugacity ( lnphi )

     xval(i) = xi(1,1)
     llemin(i)= gibbs(1) * RGAS * t

     if ( ABS(1.0-dense(1)/dense(2)) > 0.0001 ) then
        vlemin(i) = ( gibbs(1) - gibbs(2) ) * RGAS * t
        phasen(i) = 2
     else
        phasen(i) = 1
     end if

     if (i > 0 .AND. phasen(i) == 2) then

        if ( phasen(i-1) == 2 .AND. ABS( vlemin(i) + vlemin(i-1) ) <  &
                                    ABS( vlemin(i) ) + ABS( vlemin(i-1) ) ) then
           j = j + 1
           start_xv(j) = xval(i-1) + (xval(i)-xval(i-1))  &
                * ABS( vlemin(i-1) ) / ABS( vlemin(i) - vlemin(i-1) )
        end if

     end if

  end do

  do i = 2, steps - 2

     dg_dx2 = (-llemin(i-2) + 16.0*llemin(i-1) - 30.0*llemin(i)  &
          + 16.0*llemin(i+1) - llemin(i+2)) / ( 12.0 * ((xval(i)-xval(i-1))**2) )
     if ( dg_dx2 < 0.0 ) then
        k = k + 1
        start_xl(k) = xval(i)
     end if

  end do

  if ( start_xl(1) == 0.0 .AND. start_xv(1) /= 0.0 ) then
     xi(1,1) = start_xv(1)
     xi(1,2) = 1.0 - xi(1,1)
     lle_check = .false.
     if ( outp >= 1) write (*,*) 'VLE is likely', xi(1,1),xi(1,2)
  else if ( start_xl(1) /= 0.0 .AND. start_xv(1) == 0.0 ) then
     xi(1,1) = start_xl(1)
     xi(1,2) = 1.0 - xi(1,1)
     if ( outp >= 1) write (*,*) 'LLE is likely', xi(1,1),xi(1,2)
     lle_check = .true.
  else if ( start_xl(1) /= 0.0 .AND. start_xv(1) /= 0.0 ) then
     xi(1,1) = start_xv(1)
     xi(1,2) = 1.0 - xi(1,1)
     if ( outp >= 1) write(*,*) 'starting with VLE and check for LLE'
     lle_check = .true.
  else
     xi(1,1) = x_sav
     xi(1,2) = 1.0 - xi(1,1)
  end if

  CALL x_summation

end subroutine vle_min


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine rachford_rice
!
! This subroutine performs a p-T flash for a defined composition (xiF) in a
! feed-trial phase by iterating the Rachford-Rice equation.
!
! input:   xiF(:), xi(1,:), xi(2,:), dense(1), dense(2), t, p
!
! output:  rhoi1(:), rhoi2(:)
!          xi(1,:), xi(2,:)
!          values of: dense(:), lnphi(:,:) are also converged values, if converg=1
!          ph_frac (phase fraction) could be defined as output, but is not yet.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine rachford_rice ( converg, rhoi1, rhoi2 )

  use BASIC_VARIABLES
  use optimizer_2D, only: Newton_Opt_analytic
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in out)                 :: converg
  real, dimension(nc), intent(out)        :: rhoi1
  real, dimension(nc), intent(out)        :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                 :: k1, k2
  integer                                 :: nr_inside_steps, nr_outside_steps
  real, parameter                         :: tol_out = 1.E-11
  real, parameter                         :: tol_in = 1.E-13
  real                                    :: ph_frac, ph_frac0
  real                                    :: lnphi(np,nc), Ki(nc)
  real                                    :: f1, d_ph, dfdph
  real                                    :: error1, error2 = 0.0
  real                                    :: xi_compare(2,nc)
  logical                                 :: comp_balance_violation = .false.

  integer                                 :: i
  real, allocatable, dimension(:)         :: w
  real, allocatable, dimension(:)         :: x_in, grad
  real, allocatable, dimension(:,:)       :: hessian
  real                                    :: phi, V1, V2, V, phi_00, f_tpd
  !-----------------------------------------------------------------------------
  outp = 3
  converg = 0
  nphas = 2                   ! number of phases

  nr_outside_steps = 80
  nr_inside_steps = 10

  xi_compare(1:2,:) = xi(1:2,:)

  ensemble_flag = 'tp'
  densta(1) = dense(1)        ! Index 1 is for liquid density (here: packing fraction eta)
  densta(2) = dense(2)        ! Index 2 is for vapour density (here: packing fraction eta)

  ph_frac = 0.5
  d_ph = 0.00001
  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'RR, N_out, N_inner, xi(1,1), xiF(1), xi(2,1), ph_frac, error1'

  !-----------------------------------------------------------------------------
  ! start iteration
  !-----------------------------------------------------------------------------

  k1 = 0
  error1 = tol_out + 1.0
  do while ( error1 > tol_out .AND. k1 < nr_outside_steps )

     !---------------------------------------------------------------------------
     ! outer loop (converging the mole fractions)
     !---------------------------------------------------------------------------

     k1 = k1 + 1

     CALL FUGACITY ( lnphi )

     Ki(1:ncomp) = EXP( lnphi(1,1:ncomp) - lnphi(2,1:ncomp) )
     write (*,*) 'dense',dense(1),dense(2)

     if ( ABS(1.0-dense(1)/dense(2)) < 1.E-6 .AND. SUM(ABS(Ki(1:ncomp)-1.0)) < 1.E-8 ) exit

     k2 = 0
     error2 = 1.0
     do while ( error2 > tol_in .AND. k2 < nr_inside_steps )

        !------------------------------------------------------------------------
        ! inner loop (converging the fraction of one phase)
        !------------------------------------------------------------------------

        k2 = k2 + 1

        ph_frac0 = ph_frac
        f1 = SUM( xiF(1:ncomp)*(Ki(1:ncomp)-1.0) / ( 1.0 + (Ki(1:ncomp)-1.0)*ph_frac ) )

        dfdph = - SUM( xiF(1:ncomp) / ( 1.0/(Ki(1:ncomp)-1.0) + ph_frac )**2 )
        if ( dfdph == 0.0 ) dfdph = 0.0000001

        ph_frac = ph_frac0 - f1 / dfdph

        error2 = ABS( ph_frac - ph_frac0 )
        if ( outp >= 3 ) write (*,'(a,3G20.8)') 'phase_fraction, error', ph_frac, ph_frac0, error2

     end do

     !---------------------------------------------------------------------------
     ! constrain the phase fraction to ( 0 < ph_frac < 1 )
     !---------------------------------------------------------------------------
     comp_balance_violation = .false.
     if ( ph_frac > 1.0 ) comp_balance_violation = .true.
     if ( ph_frac < 0.0 ) comp_balance_violation = .true.
     if ( ph_frac > 1.0 ) ph_frac = 0.9999999
     if ( ph_frac < 0.0 ) ph_frac = 0.0000001
     !if ( k2 >= nr_inside_steps ) write (*,*) 'Rachford-Rice: Newton-loop not converged'

     !---------------------------------------------------------------------------
     ! determine x and error of outer loop
     !---------------------------------------------------------------------------
     xi(1,1:ncomp) =              xiF(1:ncomp) / ( 1.0 + ph_frac *( Ki(1:ncomp)-1.0 ) )
     xi(2,1:ncomp) = Ki(1:ncomp)* xiF(1:ncomp) / ( 1.0 + ph_frac *( Ki(1:ncomp)-1.0 ) )

     if ( sum( xi(2,1:ncomp) ) < 0.2 ) xi(2,1:ncomp) = xi(2,1:ncomp) / sum( xi(2,1:ncomp) )
     if ( sum( xi(1,1:ncomp) ) < 0.2 ) xi(1,1:ncomp) = xi(1,1:ncomp) / sum( xi(1,1:ncomp) )

     error1 =  SUM( ABS( xi_compare(1,1:ncomp) - xi(1,1:ncomp) ) )  &
          + SUM( ABS( xi_compare(2,1:ncomp) - xi(2,1:ncomp) ) )
     xi_compare(1:2,:) = xi(1:2,:)

     if ( outp >= 2 ) write (*,'(a,2i5,5G20.12)') 'RR',k1,k2,xi(1,1), xiF(1),xi(2,1), ph_frac, error1

     if ( k1 == 100 .AND. outp > 0 ) write (*,*) 'Rachford-Rice: outside loop slow convergence'

  end do


  !-----------------------------------------------------------------------------
  ! if convergence: accept solution
  !-----------------------------------------------------------------------------

  if ( error1 < tol_out .AND. error2 < tol_in .AND. .NOT. comp_balance_violation ) then

     rhoi1( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )
     rhoi2( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )
     converg = 1
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'Rachford-Rice converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) 'eta 1,2 = ',dense(1), dense(2)
     if ( outp >= 2 ) write (*,*) ' '
     alpha = ph_frac

  end if

  !-----------------------------------------------------------------------------
  ! for slow convergence, take second order approach.
  ! algorithm proposed by Nagarajan, Cullick, Griewank (FPE, 1991, 191-210).
  ! So far, a rough criterion for the approach to convergence is used.
  !-----------------------------------------------------------------------------

  if ( converg /= 1 .AND. error1 < 1.E-4 ) then

     if( outp >= 1) write (*,*) 'Rachf.-R.: not converged'

     phi_00 = 0.5
     ensemble_flag = 'tv'
     allocate ( w(ncomp), x_in(ncomp+2), grad(ncomp+2), hessian(ncomp+2,ncomp+2) )
     rhoi1(1:ncomp) = rhoi_cal(1,1:ncomp)
     rhoi2(1:ncomp) = rhoi_cal(2,1:ncomp)
     w(1:ncomp) = rhoi2(1:ncomp) - rhoi1(1:ncomp)
     V1 = (1.0-phi_00) * sum( xiF(1:ncomp) ) / sum( rhoi1(1:ncomp) )
     V2 = phi_00       * sum( xiF(1:ncomp) ) / sum( rhoi2(1:ncomp) )
     write (*,*) 'phi_00, xiF(1:ncomp)', phi_00, xiF(1:ncomp)
     write (*,*) 'V1, V2', V1, V2
     write (*,*) 'check:',1, (1.0-ph_frac)*xi(1,1) + ph_frac*xi(2,1), xif(1)
     write (*,*) 'check:',2, (1.0-ph_frac)*xi(1,2) + ph_frac*xi(2,2), xif(2)
     V = V1 + V2
     phi = 0.5 * ( V1 - V2 ) / V
     x_in(1:ncomp) = w(1:ncomp)
     x_in(ncomp+1) = phi
     x_in(ncomp+2) = V
     call Newton_Opt_analytic ( Helmholtz_flash_grads, (ncomp+2), x_in, 1.E-7, 1.E-12, grad, f_tpd, converg )

     w(1:ncomp) = x_in(1:ncomp)
     phi    = x_in(ncomp+1)
     V      = x_in(ncomp+2)
     do i = 1, ncomp
        rhoi1(i) = xiF(i) / V + ( phi - 0.5 ) * w(i)
        rhoi2(i) = xiF(i) / V + ( phi + 0.5 ) * w(i)
        if ( rhoi1(i) < 0.0 ) rhoi1(i) = 0.0
        if ( rhoi2(i) < 0.0 ) rhoi2(i) = 0.0
     end do

     if ( maxval(abs(grad(:))) < 1.E-7 ) then

        converg = 1
        if ( outp >= 2 ) write (*,*) ' '
        if ( outp >= 2 ) write (*,*) '========================================================'
        if ( outp >= 1 ) write (*,*) 'Newton iteration (subsequent to Rachford-Rice) converged'
        if ( outp >= 2 ) write (*,*) '========================================================'
        if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
        if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
        if ( outp >= 2 ) write (*,*) 'eta 1,2 = ',dense(1), dense(2)
        if ( outp >= 2 ) write (*,*) ' '
        alpha = V * ( 0.5 - phi )*sum( rhoi2(1:ncomp) )

     end if

     deallocate( w, x_in, grad, hessian )

  end if


end subroutine rachford_rice



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine bubble_point_rachford_rice
!
! This subroutine performs a bubble- or dew-point calculation for a defined
! composition (xiF). A bubble point calculation for a given liquid composition
! is done, when dense(1) is a liquid density (starting value, needs to be
! possible to converge). A dew point calculation is done, when dense(1) (low
! starting value) can converge to a vapor density.
!
! Either, the temperature or the pressure is iterated. The algorith is based on
! the Rachford-Rice equation.
!
! input:   xi(1,:),  t or p
!          And starting values: dense(1), dense(2), p or t
!
! output:  p or t, xi(2,:),  rhoi1(:), rhoi1(:)
!          values of: dense(:), lnphi(:,:) are also converged values, if converg=1
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine bubble_point_rachford_rice ( iterate_t, converg, rhoi1, rhoi2 )

  use BASIC_VARIABLES
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(in)                     :: iterate_t
  integer, intent(in out)                 :: converg
  real, dimension(nc), intent(out)        :: rhoi1
  real, dimension(nc), intent(out)        :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                 :: k1, k2
  integer                                 :: nr_inside_steps, nr_outside_steps
  real, parameter                         :: tol_out = 1.E-11
  real, parameter                         :: tol_in_p = 1.E-10
  real, parameter                         :: tol_in_t = 1.E-8
  real                                    :: tol_in
  real                                    :: p0, t0, dp, dt, deltap
  real                                    :: p_sav, t_sav
  real                                    :: damping, acceleration
  real                                    :: lnphi(np,nc), Ki(nc)
  real                                    :: f0, f1, dfdp
  real                                    :: error1, error2 = 0.0
  real                                    :: xi_compare(nc)
  logical                                 :: SR_exit
  !-----------------------------------------------------------------------------

  t_sav = t
  p_sav = p

  converg = 0
  nphas = 2                   ! number of phases

  nr_outside_steps = 1000
  nr_inside_steps = 10
  if ( iterate_t == 1 ) tol_in = tol_in_t
  if ( iterate_t == 0 ) tol_in = tol_in_p
  SR_exit = .false.

  xi_compare(:) = xi(2,:)

  ensemble_flag = 'tp'
  densta(1) = dense(1)        ! Index 1 is for liquid density (here: packing fraction eta)
  densta(2) = dense(2)        ! Index 2 is for vapour density (here: packing fraction eta)

  dp = 1.0                    ! pressure step in unit [Pa]
  dt = 1.E-3                  ! temperature step in unit [K]
  if ( outp >= 2 ) write (*,'(a,4G20.12)') 'bbp-RR,   N_out, N_inner, xi(1,1),   xi(2,1),   p,   t,   error1'

  !-----------------------------------------------------------------------------
  ! start iteration
  !-----------------------------------------------------------------------------

  k1 = 0
  error1 = tol_out + 1.0
  do while ( error1 > tol_out .AND. k1 < nr_outside_steps )

     !--------------------------------------------------------------------------
     ! outer loop (converging the mole fractions)
     !--------------------------------------------------------------------------

     k1 = k1 + 1

     k2 = 0
     error2 = tol_in + 1.0
     do while ( error2 > tol_in .AND. k2 < nr_inside_steps )

        !-----------------------------------------------------------------------
        ! inner loop (converging the pressure)
        !-----------------------------------------------------------------------

        k2 = k2 + 1

        p0 = p
        t0 = t
        if ( iterate_t == 1 ) then
           t = t0 - dt
        else
           p = p0 - dp
        end if

        CALL FUGACITY ( lnphi )

        if ( ABS(1.0-dense(1)/dense(2)) < 1.E-6 .AND. SUM(ABS(Ki(1:ncomp)-1.0)) < 1.E-8 ) exit
        Ki(1:ncomp) = EXP( lnphi(1,1:ncomp) - lnphi(2,1:ncomp) )
        f0 = SUM( xi(1,1:ncomp) * Ki(1:ncomp) )   - 1.0

        p = p0
        t = t0
        CALL FUGACITY ( lnphi )

        Ki(1:ncomp) = EXP( lnphi(1,1:ncomp) - lnphi(2,1:ncomp) )
        f1 = SUM( xi(1,1:ncomp) * Ki(1:ncomp) )   - 1.0

        if ( iterate_t == 1 ) then

           dfdp = ( f1 - f0 ) / dt
           if ( dfdp == 0.0 ) dfdp = 0.0000001

           !--- Newton ---------------------------------------------------------
           t = t0 - f1 / dfdp

           !--- error ----------------------------------------------------------
           error2 = ABS( t - t0 )
           if ( outp >= 3 ) write (*,'(a,3G20.12)') 'inner loop error', t, t0, error2
           if ( error2 > 2.0 .OR. t < 0.0 ) SR_exit = .true.
           if ( error2 > 2.0 .OR. t < 0.0 ) exit

        else

           dfdp = ( f1 - f0 ) / ( LOG(p0) - LOG(p0 - dp) )
           if ( dfdp == 0.0 ) dfdp = 0.0000001

           !--- Newton ---------------------------------------------------------
           ! p = EXP( LOG(p0) - f1 / dfdp  )
           deltap = f1 / dfdp
           if ( deltap >  2.0 ) deltap =  2.0
           if ( deltap < -2.0 ) deltap = -2.0
           p = LOG( p0) - deltap
           p = EXP( p )
           p = min( p, 1.E9 )

           !--- damping --------------------------------------------------------
           damping = 0.6
           if ( ABS( p / p0 - 1.0 ) > 0.3 ) p = damping * p + ( 1.0 - damping ) * p0

           !--- error ----------------------------------------------------------
           error2 = ABS( p / p0 - 1.0 )

           if ( outp >= 3 ) write (*,'(a,3G20.12)') 'inner loop error', p, p0, error2
           if ( error2 > 10.0 ) SR_exit = .true.
           if ( error2 > 10.0 ) exit

        end if

     end do

     if ( SR_exit ) exit

     !--------------------------------------------------------------------------
     ! determine x and error of outer loop
     !--------------------------------------------------------------------------
     xi(2,1:ncomp) = Ki(1:ncomp)* xi(1,1:ncomp) / SUM( xi(1,1:ncomp) * Ki(1:ncomp) )
     xi(2,1:ncomp) = max( xi(2,1:ncomp), 1.E-50 )

     if ( sum( xi(2,1:ncomp) ) < 0.4 ) xi(2,1:ncomp) = xi(2,1:ncomp) / sum( xi(2,1:ncomp) )

     error1 =  SUM( ABS( xi_compare(1:ncomp) - xi(2,1:ncomp) ) )

     if ( k1 > 100 ) then
        acceleration = real(k1) / 100.0 + 1.0
        if ( k1 > 900 ) acceleration = acceleration * 2.0
        xi(2,1:ncomp) = xi(2,1:ncomp) + acceleration * ( xi(2,1:ncomp) - xi_compare(1:ncomp) )  ! acceleration
        xi(2,1:ncomp) = max( xi(2,1:ncomp), 1.E-50 )
        ! write (*,'(a,4g20.12)') 'accel',acceleration,xi(2,1), xi(2,2), sum( xi(2,1:ncomp) )
        xi(2,1:ncomp) = xi(2,1:ncomp) / sum( xi(2,1:ncomp) )
     end if

     xi_compare(:) = xi(2,:)

     if ( error2 > tol_in ) error1 = error1 + tol_out

     if ( k1 >= 500 ) error1 = error1 / 10.0
     if ( k1 >= 1000 ) error1 = error1 / 100.0

     if ( outp >= 2 ) write (*,'(a,2i5,5G18.10)') 'bbp-RR',k1,k2,xi(1,1), xi(2,1), p, t, error1

     if ( k1 == 100 .AND. outp > 0 ) write (*,*) 'bbp-Rachford-Rice: outside loop slow convergence'

  end do


  !-----------------------------------------------------------------------------
  ! if convergence: accept solution
  !-----------------------------------------------------------------------------

  if ( error1 < tol_out .AND. error2 < tol_in ) then

     rhoi1( 1:ncomp ) = rhoi_cal( 1, 1:ncomp )
     rhoi2( 1:ncomp ) = rhoi_cal( 2, 1:ncomp )
     converg = 1
     if ( outp >= 2 ) write (*,*) ' '
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 1 ) write (*,*) 'bubble point Rachford-Rice converged'
     if ( outp >= 2 ) write (*,*) '========================================================'
     if ( outp >= 2 ) write (*,*) '  x p1   = ',rhoi1( 1:ncomp ) / sum( rhoi1( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) '  x p2   = ',rhoi2( 1:ncomp ) / sum( rhoi2( 1:ncomp ) )
     if ( outp >= 2 ) write (*,*) 'eta 1,2 = ',dense(1), dense(2)
     if ( outp >= 2 ) write (*,*) ' '

  else

     t = t_sav
     p = p_sav
     if ( outp >= 1) write (*,*) 'bbp-Rachf.-R.: not converged'

  end if

end subroutine bubble_point_rachford_rice



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine T-x rachford_rice
!
! This subroutine solves the a T-x flash problem, i.e. for specified T and x of
! species index 2 in the phase 1 (usually liquid).
!
! input:   xiF(1:ncomp), xi(1,1:ncomp), xi(2,1:ncomp), dense(1), dense(2), t, p
!
! output:  rhoi1( 1:ncomp ), rhoi2( 1:ncomp )
!          xi( 1, 1:ncomp ), xi( 2, 1:ncomp )
!          values of: dense(:), lnphi(:,:) are also converged values, if converg=1
!          ph_frac (phase fraction) could be defined as output, but is not yet.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine T_x_rachford_rice ( x12, converg, rhoi1, rhoi2 )

  use BASIC_VARIABLES
  implicit none

  !-----------------------------------------------------------------------------
  real, intent(in)                        :: x12
  integer, intent(in out)                 :: converg
  real, dimension(nc), intent(out)        :: rhoi1
  real, dimension(nc), intent(out)        :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                 :: count
  real                                    :: f, fr, dfdx
  real                                    :: delta
  real                                    :: p0
  !-----------------------------------------------------------------------------

  delta = 0.00001


  count = 0
  f = 1.0
  do while ( ABS( f ) > 1.E-6 .AND. count < 10 )

     count = count + 1
     p0 = p

     p = p0 * ( 1.0 + delta )
     call rachford_rice ( converg, rhoi1, rhoi2 )
     fr = xi(1,2) - x12

     p = p0
     call rachford_rice ( converg, rhoi1, rhoi2 )
     f = xi(1,2) - x12

     dfdx = ( fr - f ) / ( p0 * delta )

     p = p0 - f / dfdx
     !write (*,*) f, p0

  end do

end subroutine T_x_rachford_rice



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine occupy_val_init ( rhoi1, rhoi2 )

  use BASIC_VARIABLES
  use EOS_VARIABLES, only: PI, dhs, mseg
  implicit none

  !-----------------------------------------------------------------------------
  real, dimension(nc), intent(in)         :: rhoi1
  real, dimension(nc), intent(in)         :: rhoi2

  !-----------------------------------------------------------------------------
  integer                                 :: i, ph
  !-----------------------------------------------------------------------------

  dense(1) = PI/6.0 * SUM( rhoi1(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
  dense(2) = PI/6.0 * SUM( rhoi2(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  xi( 1, 1:ncomp ) = rhoi1( 1:ncomp ) / sum( rhoi1(1:ncomp) )
  xi( 2, 1:ncomp ) = rhoi2( 1:ncomp ) / sum( rhoi2(1:ncomp) )
  lnx( 1, 1:ncomp ) = LOG( xi( 1, 1:ncomp ) )
  lnx( 2, 1:ncomp ) = LOG( xi( 2, 1:ncomp ) )

  val_init(1)  = dense(1)
  val_init(2)  = dense(2)
  val_init(3)  = t
  val_init(4)  = p
  do ph = 1, 2
     do i = 1, ncomp
        val_init(4+i+(ph-1)*ncomp) = lnx( ph, i )
     end do
  end do

end subroutine occupy_val_init



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine read_mode_staring_value ( scan_index, outp )

  use utilities, only: file_open
  implicit none

  !-----------------------------------------------------------------------------
  integer, intent(out)                    :: scan_index
  integer, intent(in out)                 :: outp

  !-----------------------------------------------------------------------------
  character (LEN=50)                      :: textstring
  character (LEN=50)                      :: start_value_file
  !-----------------------------------------------------------------------------

  first_call = .false.
  start_value_file = './input_file/starting_value_default.inp'
  call file_open ( start_value_file, 84 )
  read (84,*) textstring
  read (84,*) scan_index, outp
  close (84)
  if ( outp == 3 ) then
     if ( scan_index == 1 ) outp = 0
     if ( scan_index == 2 ) outp = 1
  end if
  if ( outp >= 1 ) write (*,*) 'level of info given to terminal:',outp
  write (*,*) 'calculation option composition-scan:',scan_index

end subroutine read_mode_staring_value



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine output_rho_info ( rhoi )

  use EOS_VARIABLES, only: ncomp, nc, PI, dhs, mseg
  implicit none
  real, dimension(nc), intent(IN)        :: rhoi
  real                                   :: dens
  !-----------------------------------------------------------------------------

  dens = PI/6.0 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  write (*,*) ' eta = ', dens
  write (*,*) '   x = ', rhoi( 1:ncomp ) / sum( rhoi(1:ncomp) )

end subroutine output_rho_info



end module STARTING_VALUES
