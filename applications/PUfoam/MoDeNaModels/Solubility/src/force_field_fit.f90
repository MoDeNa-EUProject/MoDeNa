!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE pure_fit_parameters
!
! This module contains parameters and variables needed for the pure
! component parameter regression.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

MODULE FF_FIT_PARAMETERS

  implicit none
  save

  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER          :: MMX = 3500

  INTEGER                     :: nlv
  INTEGER                     :: polar
  INTEGER                     :: comp_index_t

  REAL, DIMENSION(MMX)        :: plvcal, rl_cal,rv_cal, dh_calc
  REAL, DIMENSION(MMX)        :: plv, tlv
  REAL, DIMENSION(MMX)        :: rliq, rvap, dh_vap
  REAL, DIMENSION(MMX)        :: std_plv, std_rliq, std_rvap, std_caloric
  REAL, DIMENSION(MMX)        :: charge_val

  REAL, DIMENSION(10)         :: phi_ff_s, phi_ff_e, phi_ff_m
  REAL, DIMENSION(10)         :: dipole_scale

  CHARACTER (LEN=30), DIMENSION(MMX)  :: comp_name
  INTEGER, DIMENSION(MMX)     :: comp_index
  CHARACTER (LEN=10)          :: fit_what

  REAL, DIMENSION(MMX,4)      :: eos_misfit, misfit

  INTEGER                     :: number_substances
  REAL, DIMENSION(20)         :: T_low_substances, T_high_substances
  CHARACTER (LEN=30), DIMENSION(20)  :: name_substances

  REAL, DIMENSION(MMX,4)      :: ff_weights

  REAL                        :: size_scale_mu

End Module FF_FIT_PARAMETERS




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE OBJECTIVE_FF
!
! Routine for fitting force field parameters. The residual between
! calculated values and simulated data is calculated and stored in
! vector 'fvec'. The vector 'fvec' serves as the objective function for
! the parameter.
! OBJECTIVE_FF is called by FITTING via the Levenberg-Marquardt routine.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE OBJECTIVE_FF (lm_m_dat, n_fit, para, fvec, iflag)

  USE BASIC_VARIABLES
  USE FF_FIT_PARAMETERS
  USE MOD_PARA_TRANSFER
  use utilities, only: SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! INTEGER, PARAMETER                     :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, INTENT(IN)                    :: lm_m_dat
  INTEGER, INTENT(IN)                    :: n_fit
  REAL, INTENT(IN)                       :: para(:)
  REAL, INTENT(IN OUT)                   :: fvec(:)
  INTEGER, INTENT(IN OUT)                :: iflag

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE eos_from_sim (comp_name, comp_index, fit_what, n_fit, para, phi_ff_s, phi_ff_e)
       IMPLICIT NONE
       CHARACTER (LEN=30), INTENT(IN)     :: comp_name
       INTEGER, INTENT(IN)                :: comp_index
       CHARACTER (LEN=10), INTENT(IN)     :: fit_what
       INTEGER, INTENT(IN)                :: n_fit
       REAL, INTENT(IN)                   :: para(n_fit)
       REAL, INTENT(OUT)                  :: phi_ff_s(10)
       REAL, INTENT(OUT)                  :: phi_ff_e(10)
     END SUBROUTINE eos_from_sim
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  INTEGER                                :: j_dat
  INTEGER                                :: converg, crit_dat

  REAL                                   :: density(np), w(np,nc)
  REAL                                   :: tc, v_dum !, pc, rhoc
  ! REAL                                   :: t_org
  REAL                                   :: error

  CHARACTER (LEN=30)                     :: comp_name_t
  CHARACTER (LEN=4)                      :: char_len
  !-----------------------------------------------------------------------------


  nphas  = 2
  n_unkw = ncomp            ! number of quantities to be iterated
  it(1)  = 'lnp'            ! iteration of pressure

  j_dat = 0                 ! counter for the entries to the objective function

  tc = 0.0

  fvec(:) = 0.0             ! initialize objective-function vector

  ! IF ( polar == 1 ) dipole_scale=para(N_fit) !*4.926006705428530

  WRITE (char_len,'(I3)') n_fit
  IF (iflag == 1) WRITE(*,'(a,'//char_len//'(G15.8))') ' parameter ',(para(i), i = 1, n_fit)



  !-----------------------------------------------------------------------------
  ! start calculating the phase equilibrium. For all rows (i = 1,nlv) in
  ! GEMC.dat, VLE calculations are carried out. After this the complete
  ! objective function (fvec) can be calculated.
  !-----------------------------------------------------------------------------

  DO i = 1, nlv

     crit_dat  = 0
     plvcal(i) = 0.0
     rl_cal(i) = 0.0             ! liquid density in SI-units in [kg/m**3]
     rv_cal(i) = 0.0             ! vapor  density in SI-units in [kg/m**3]
     dh_calc(i)= 0.0             ! enthalpy of vaporization in [J/kg]

     parame(1,6) = 0.0

!!$     IF ( fit_what == 'phi_to_sim' ) para_transfer(1) = para(N_fit) * 3500.0

     !--------------------------------------------------------------------------
     ! calculate EOS parameter from force-field parameters
     !--------------------------------------------------------------------------
     comp_name_t = comp_name(i)
     comp_index_t = comp_index(i)

     CALL eos_from_sim ( comp_name_t, comp_index_t, fit_what, n_fit, para, phi_ff_s, phi_ff_e )
     ! write (*,'(a,a,5G18.10)') 'pure_pars',comp_name_t,parame(1,1:5)
     ! call paus (' ')

     IF ( polar == 1 ) THEN
        IF ( fit_what == 'phi_to_sim' ) THEN
           dipole_scale( comp_index_t ) = 4.5 * para( N_fit - number_substances + comp_index_t )
           parame(1,6) = charge_val(i) * dipole_scale( comp_index_t ) ! overwrite dipole moment
        ELSE
           ! parame(1,6) = dipole_scale( comp_index_t ) * para(N_fit)
           parame(1,6) = dipole_scale( comp_index_t ) * charge_val(i)
           ! if (i==1) write (*,*) ' did not allow mu to change',para(N_fit), charge_val(i)
           ! IF (iflag == 1 .AND. i==1) write (*,*) ' dipole_scale=',dipole_scale( comp_index_t ),parame(1,6)
           parame(1,6) = SQRT( size_scale_mu ) * parame(1,6)
           if (i==1) write (*,*) ' size_scale_mu', size_scale_mu, num, parame(1,3)
        END IF
     END IF
     ! if ( i == 1) write (*,*) 'eps_assoc',parame(1,15)
     IF (iflag == 1 .AND. i == 1 ) WRITE(*,'(a,G15.8)') ' parame(1,3) ',parame(1,3)

     val_init(1) = 0.5           ! initial guess for a liquid density
     val_init(2) = 1.E-8         ! initial guess for a vapor density
     val_init(3) = tlv(i)
     val_init(4) = plv(i)
     val_init(5) = 0.0           ! mole fraction lnx(1,1)
     val_init(6) = 0.0           ! mole fraction lnx(2,1)

     ! IF (iflag == 1 ) write (*,'(i3,x,a,2G18.10)') i,comp_name_t,parame(1,6),charge_val(i)


     !--------------------------------------------------------------------------
     ! calculate phase equilibrium
     !--------------------------------------------------------------------------
     CALL pure_equilibrium_fit ( crit_dat, tc, tlv(i), plv(i), converg, v_dum, plvcal(i) )


     IF ( converg == 1 ) THEN
        CALL SI_dens (density,w)
        rl_cal(i) = density(1)  ! liquid density in SI-units in [kg/m**3]
        rv_cal(i) = density(2)  ! vapor  density in SI-units in [kg/m**3]
        dh_calc(i) = h_lv / mm(1) * 1.d3    ! enthalpy of vaporization in [J/kg]
     ELSE
        IF ( crit_dat == 1 .AND. tc > tlv(i) ) THEN
           plvcal(i) = plv(i) * 0.75
           WRITE (*,'(a,i3,2G14.5,G15.8)') ' loose ! ',i,tlv(i),t,plv(i)
        ELSE
           plvcal(i) = plv(i) * 0.75
        END IF
     END IF


     !--------------------------------------------------------------------------
     ! liquid density data
     !--------------------------------------------------------------------------
     j_dat = j_dat + 1
     fvec(j_dat) = ff_weights(i,2) * ABS((rl_cal(i)-rliq(i)+misfit(i,1))/rliq(i)) ! / std_rliq(i)
     if (converg == 1) fvec(j_dat) = ff_weights(i,2) * ABS((rl_cal(i)-rliq(i)+misfit(i,1))/SQRT(rliq(i)*rl_cal(i)))
     eos_misfit(i,1) = rliq(i)-rl_cal(i)
     ! write (*,'(a,4G18.10)') 'rho',rl_cal(i),rl_cal(i)+misfit(i,1),rliq(i) !,ff_weights(i,2)

     ! IF (iflag == 1) write (*,*) i,comp_index_t,comp_name_t,parame(1,6)


     !--------------------------------------------------------------------------
     ! vapor pressure data
     !--------------------------------------------------------------------------
     j_dat = j_dat + 1
     fvec(j_dat) = ff_weights(i,1) * ABS((plvcal(i)-plv(i)+misfit(i,2))/plv(i))
     eos_misfit(i,2) = plv(i) - plvcal(i)
     ! write (*,*) 'p',plvcal(i),plvcal(i)+misfit(i,2),plv(i)

     ! penalty function for calc. crit. point < exp. crit. point
     !IF ( tc < tlv(i) .AND. crit_dat == 1 ) THEN
     !  fvec(j_dat) = fvec(j_dat) + 5.0 * ( tlv(i) + 2.0 - tc )
     !  WRITE (*,'(a,i3,3(f15.5))') ' skipped ! ', i, tlv(i), tc
     !END IF

     ! IF (iflag == 1 .AND. fit_what == 'ffp_to_exp' ) write (*,'(i3,3G18.10)') i, &
     !        ((plvcal(i)+misfit(i,2))/plv(i) - 1.0)*1.E2,  plvcal(i)+misfit(i,2),  plv(i)


     !--------------------------------------------------------------------------
     ! enthalpy of vaporization
     !--------------------------------------------------------------------------
     j_dat = j_dat + 1
     fvec(j_dat) = ff_weights(i,4) * ABS((dh_calc(i)-dh_vap(i)+misfit(i,3))/dh_vap(i))
     eos_misfit(i,3) = dh_vap(i) - dh_calc(i)
     ! write (*,*) 'h', dh_calc(i), dh_calc(i)+misfit(i,3), dh_vap(i)


     !--------------------------------------------------------------------------
     ! vapor density data
     !--------------------------------------------------------------------------
     j_dat = j_dat + 1
     fvec(j_dat) = ff_weights(i,3) * ABS((rv_cal(i)-rvap(i)+misfit(i,4))/rvap(i))
     if (converg == 1) fvec(j_dat) = ff_weights(i,3) * ABS((rv_cal(i)-rvap(i)+misfit(i,4))/SQRT(rvap(i)*rv_cal(i)))
     eos_misfit(i,4) = rvap(i) - rv_cal(i)
     ! write (*,*) 'r_v',rv_cal(i),rv_cal(i)+misfit(i,4),rvap(i)
     ! if ( fit_what == 'ffp_to_exp' ) write (*,*) (plvcal(i)+misfit(i,2)),  tlv(i),  rl_cal(i)+misfit(i,1)
     ! write (*,'(i5,6G18.9)') i, charge_val(i), plv(i),  tlv(i),  rliq(i),  plvcal(i), rl_cal(i)

  END DO


  !-----------------------------------------------------------------------------
  ! calculating the Absolute Average Deviation in % (AAD)
  !-----------------------------------------------------------------------------
  OPEN (60, FILE = './output_file/AAD')

  error = SUM( fvec(1:lm_m_dat) / ( REAL(nlv)*SUM(ff_weights(1,1:4)) ) * 100.0 )        ! AAD-%
  WRITE (*,*) ' AAD (%)', error,' RMS (%)',SQRT( SUM( fvec(1:lm_m_dat)*fvec(1:lm_m_dat)  &
       / (REAL(nlv)*SUM(ff_weights(1,1:4))) ) )* 100.0, nlv*NINT(SUM(ff_weights(1,1:4)))

  WRITE(60,*) fit_what, error
  ! IF (iflag == 1 .AND. fit_what /= 'phi_to_sim' ) call paus (' ')

END SUBROUTINE objective_ff



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE eos_from_sim
!
! Subroutine calculates PC_SAFT parameters from molecular
! simulation-force fields (ff)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE eos_from_sim (comp_name, comp_index, fit_what, n_fit, para, phi_ff_s, phi_ff_e)

  USE BASIC_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  CHARACTER (LEN=30), INTENT(IN)         :: comp_name
  INTEGER, INTENT(IN)                    :: comp_index
  CHARACTER (LEN=10), INTENT(IN)         :: fit_what
  INTEGER, INTENT(IN)                    :: n_fit
  REAL, INTENT(IN)                       :: para(n_fit)
  REAL, INTENT(OUT)                      :: phi_ff_s(10)
  REAL, INTENT(OUT)                      :: phi_ff_e(10)

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE gc_par (comp_name, comp_index, fit_what, x_types, x_n, x_m, x_sig, x_eps,  &
          n_fit, para, phi_ff_s, phi_ff_e)
       IMPLICIT NONE
       CHARACTER (LEN=30), INTENT(IN)     :: comp_name
       INTEGER, INTENT(IN)                :: comp_index
       CHARACTER (LEN=10), INTENT(IN)     :: fit_what
       INTEGER, INTENT(OUT)               :: x_types
       REAL, DIMENSION(10), INTENT(OUT)   :: x_n, x_m, x_sig, x_eps
       INTEGER, INTENT(IN)                :: n_fit
       REAL, INTENT(IN)                   :: para(n_fit)
       REAL, DIMENSION(10), INTENT(OUT)   :: phi_ff_s, phi_ff_e
     END SUBROUTINE gc_par
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i,j
  INTEGER                                :: x_types
  REAL, DIMENSION(10)                    :: x_n, x_m, x_sig, x_eps
  !-----------------------------------------------------------------------------

  CALL gc_par (comp_name, comp_index, fit_what, x_types, x_n, x_m, x_sig, x_eps,  &
       n_fit, para, phi_ff_s, phi_ff_e)

  !-----------------------------------------------------------------------------
  ! calculate m
  !-----------------------------------------------------------------------------
  parame(1,1) = 0.0
  DO i = 1, 10
     parame(1,1) = parame(1,1) + x_m(i) * x_n(i)
  END DO

  !-----------------------------------------------------------------------------
  ! calculate sigma
  !-----------------------------------------------------------------------------
  parame(1,2) = 0.0
  DO i=1,10
     parame(1,2) = parame(1,2) +x_m(i)*x_n(i)*(x_sig(i))**3
  END DO
  parame(1,2) = parame(1,2) / parame(1,1)
  parame(1,2) = (parame(1,2))**(1.0/3.0)

  !-----------------------------------------------------------------------------
  ! calculate epsilon/k
  !-----------------------------------------------------------------------------
  parame(1,3) = 0.0
  DO i = 1, 10
     DO j = 1, 10
        parame(1,3) = parame(1,3) + x_n(i)*x_n(j) *((x_sig(i)+x_sig(j))/2.0)**3  &
             * (x_eps(i)*x_eps(j))**0.5
     END DO
  END DO
  parame(1,3) = parame(1,3)/parame(1,1)**2 / (parame(1,2))**3

!!$  parame(1,3) = 0.0
!!$  DO i = 1, 10
!!$     parame(1,3) = parame(1,3) + x_n(i) * x_eps(i)
!!$  END DO
!!$  parame(1,3) = parame(1,3)/parame(1,1)
!!$  ! write (*,*) 'BBB',parame(1,3)
!!$  ! if ( fit_what /= 'phi_to_sim' ) call paus (' ')


END SUBROUTINE eos_from_sim



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE gc_par
!
! Obtain group-contribution parameters
! comp_name : name
! x_types   : number of different types of groups making up a molecule
! x_N       : number of groups of one type (double precision !!)
! x_m       : number of tangent-sphere segments a group is equivalent with
! x_sig     : LJ-sigma parameter from MS-force field
! x_eps     : LJ-epsilon parameter from MS-force field
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE gc_par (comp_name, comp_index, fit_what, x_types, x_n_out, x_m, x_sig, x_eps,  &
     n_fit, para, phi_ff_s, phi_ff_e)

  USE BASIC_VARIABLES
  USE GC_GROUP
  USE FF_FIT_PARAMETERS, only: size_scale_mu
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  CHARACTER (LEN=30), INTENT(IN)         :: comp_name
  INTEGER, INTENT(IN)                    :: comp_index
  CHARACTER (LEN=10), INTENT(IN)         :: fit_what
  INTEGER, INTENT(OUT)                   :: x_types
  REAL, DIMENSION(10), INTENT(OUT)       :: x_n_out, x_m, x_sig, x_eps
  INTEGER, INTENT(IN)                    :: n_fit
  REAL, INTENT(IN)                       :: para(n_fit)
  REAL, DIMENSION(10), INTENT(OUT)       :: phi_ff_s, phi_ff_e

  !-----------------------------------------------------------------------------
  INTEGER                                :: i,j
  INTEGER                                :: m_fit_index, sig_fit_index, eps_fit_index

  REAL                                   :: x_n(10)
  REAL                                   :: temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling
  CHARACTER (LEN=20)                     :: x_typ_name(10)
  CHARACTER (LEN=20)                     :: temp_name
  CHARACTER (LEN=50)                     :: filename
  CHARACTER (LEN=30)                     :: comp_name_t(N_comp)
  !-----------------------------------------------------------------------------

  x_types = 0
  DO i = 1, 10
     x_n(i) = 0.0
     x_n_out(i) = 0.0
     x_m(i) = 0.0
     x_sig(i) = 0.0
     x_eps(i) = 0.0
     x_typ_name(i) = ' '
  END DO

  gc_tot_nr_groups(:) = 0.0
  gc_groups(:,:) = 0.0
  gc_grouptype(:,:) = ' '

  comp_name_t(1) = comp_name
  CALL GC_PARAMETER ( 1, comp_name_t )
  x_types = gc_tot_nr_groups(1)
  x_typ_name(1:x_types) = gc_grouptype(1,1:x_types)
  x_n(1:x_types) = REAL( gc_groups(1,1:x_types) )
  mm(1) = SUM( mm_gc(1,1:x_types) * REAL( gc_groups(1,1:x_types) ) )


  !-----------------------------------------------------------------------------
  ! reading (initial guess) forcefield parameters from input file
  !-----------------------------------------------------------------------------
  filename = 'gc.inp'
  filename = './input_file/'//filename
  CALL file_open(filename,77)

  i = 0   ! index for all groups present in the problem
  ! note, that there are possibly more groups in the problem
  ! than in the considered molecule
  read_loop: DO
     i = i + 1
     READ (77,*) temp_name, temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling, &
          m_fit_index, sig_fit_index, eps_fit_index
     IF ( temp_name == 'end' ) EXIT read_loop

     DO j = 1, x_types     ! index for groups present in the considered molecule

        IF (temp_name == x_typ_name(j)) THEN

           !--------------------------------------------------------------------
           ! set m-parameter
           !--------------------------------------------------------------------
           x_m(i) = temp_x_m

           !--------------------------------------------------------------------
           ! set sigma-parameter
           !--------------------------------------------------------------------
           size_scale_mu = 1.0
           IF ( fit_what == 'phi_to_sim' ) THEN
              phi_ff_s(comp_index) = para(2*comp_index-1)
              x_sig(i) = temp_x_sig * phi_ff_s(comp_index)
           ELSE IF ( sig_fit_index /= 0 .AND. fit_what == 'ffp_to_exp' ) THEN
              x_sig(i) = para(sig_fit_index) * phi_ff_s(comp_index)
              size_scale_mu = para(sig_fit_index) / temp_x_sig
           ELSE
              x_sig(i)= temp_x_sig * phi_ff_s(comp_index)
           END IF
           !--------------------------------------------------------------------
           ! set epsilon-parameter
           !--------------------------------------------------------------------
           ! eps_scale_mu = 1.0
           IF ( fit_what == 'phi_to_sim' ) THEN
              phi_ff_e(comp_index) = para(2*comp_index)
              x_eps(i)= temp_x_eps * LJ_scaling * phi_ff_e(comp_index)  ! note LJ_scaling carries an index i (just as temp_x_eps)
           ELSE IF ( eps_fit_index /= 0 .AND. fit_what == 'ffp_to_exp' ) THEN
              x_eps(i) = para(eps_fit_index) * LJ_scaling * phi_ff_e(comp_index)
              ! eps_scale_mu = para(eps_fit_index) / temp_x_eps
           ELSE
              x_eps(i)= temp_x_eps * LJ_scaling * phi_ff_e(comp_index)
           END IF
           !--------------------------------------------------------------------
           ! set x_N_out: number of times the group is present in molec.
           !--------------------------------------------------------------------
           x_n_out(i) = x_n(j)

        END IF

     END DO   ! loop over j=1,x_types (groups present in the molec.)

  END DO read_loop  ! loop over i (overall number of groups present in problem)
  CLOSE(77)

END SUBROUTINE gc_par



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE fitting_forcefield

  USE BASIC_VARIABLES
  USE FF_FIT_PARAMETERS
  USE Levenberg_Marquardt
  USE utilities, only: file_open, paus
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE OBJECTIVE_FF(lm_m_dat, n_fit, para, fvec, iflag)
       IMPLICIT NONE
       INTEGER, INTENT(IN)                :: lm_m_dat, n_fit
       REAL, INTENT(IN)                   :: para(:)
       REAL, INTENT(IN OUT)               :: fvec(:)
       INTEGER, INTENT(IN OUT)            :: iflag
     END SUBROUTINE OBJECTIVE_FF
  END INTERFACE

  INTERFACE
     SUBROUTINE data_in_forcefield( lm_m_dat, n_fit, para )
       IMPLICIT NONE
       INTEGER, INTENT(IN OUT)            :: n_fit
       INTEGER, INTENT(OUT)               :: lm_m_dat
       REAL, ALLOCATABLE, INTENT(OUT)     :: para(:)
     END SUBROUTINE data_in_forcefield
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                :: info
  INTEGER                                :: n_fit            ! number of adjustable parameters
  INTEGER                                :: lm_m_dat         ! number of data points
  INTEGER, ALLOCATABLE                   :: ipvt(:)
  REAL, ALLOCATABLE                      :: para(:)
  REAL, ALLOCATABLE                      :: fvec(:)
  REAL                                   :: tol    = 1.0E-8
  REAL                                   :: epsfcn
  REAL                                   :: lm_factor

  INTEGER                                :: i
  INTEGER                                :: m_fit_index, sig_fit_index, eps_fit_index
  REAL                                   :: temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling
  CHARACTER (LEN=50)                     :: filename
  CHARACTER (LEN=20)                     :: temp_name
  !-----------------------------------------------------------------------------

  parame(1,4) = 0.0
  parame(1,5) = 0.12

  misfit(:,:) = 0.0

  ncomp = 1
  eos   = 1
  pol   = 1

  !-----------------------------------------------------------------------------
  ! numerical constants for Levenberg-Marquardt algorithm
  !-----------------------------------------------------------------------------
  scaling(1) = 1.E2
  IF ( num == 0 ) epsfcn = 1.0E-6**2  ! sqrt of relat. step size (finite differences)
  IF ( num == 1 ) epsfcn = 1.0E-5**2  ! sqrt of relat. step size (finite differences)
  lm_factor = 0.2                     ! maximum initial step size (parameter iteration)


  !-----------------------------------------------------------------------------
  ! read weights for p_sat, rho_L, rho_V, h_LV
  !-----------------------------------------------------------------------------
  ! weights are defined as first line of exp. data-file
  ! in ./input_file/ff_db/... Values are used for the sim. & the exp. data
  CALL read_ff_weights




  !=============================================================================
  ! adjust correction parameters to GEMC simulation data
  !=============================================================================

  fit_what = 'phi_to_sim'

  !-----------------------------------------------------------------------------
  ! retrieving simulation values, and set N_fit (number of adjustable parameters)
  !-----------------------------------------------------------------------------
  CALL data_in_forcefield ( lm_m_dat, n_fit, para )

  ALLOCATE( ipvt(n_fit) )

  para( 1:n_fit ) = 1.0

  !para(1) = 0.9
  !para(2) = 1.6
  !para(3) = 0.9
  !para(4) = 1.6
  !para(5) = 0.9
  !para(6) = 1.6

  ALLOCATE( fvec(lm_m_dat) )


  !-----------------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  !-----------------------------------------------------------------------------

  CALL lmdif1(OBJECTIVE_FF,lm_m_dat,n_fit,para,fvec,tol,epsfcn,lm_factor,info,ipvt)

  WRITE (*,*) 'finished with phi_to_sim'

  IF ( polar == 1 ) WRITE (*,*) 'dipole-scale', para(n_fit)

  DEALLOCATE( para, ipvt )

  call paus (' ')




  !=============================================================================
  ! adjust force-field parameters to experimental data
  !=============================================================================

  fit_what = 'ffp_to_exp'

  DO i = 1, MMX
     misfit(i,1) = eos_misfit(i,1)
     misfit(i,2) = eos_misfit(i,2)
     misfit(i,3) = eos_misfit(i,3)
     misfit(i,4) = eos_misfit(i,4)
  END DO

  !-----------------------------------------------------------------------------
  ! retrieving experimental values, and set N_fit (number of adjustable parameters)
  ! and get the initial values for the adjustable ff parameters
  !-----------------------------------------------------------------------------
  CALL data_in_forcefield ( lm_m_dat, n_fit, para )

  ALLOCATE( ipvt(n_fit) )

  !-----------------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  !-----------------------------------------------------------------------------

  CALL lmdif1(OBJECTIVE_FF,lm_m_dat,n_fit,para,fvec,tol,epsfcn,lm_factor,info,ipvt)

  !-----------------------------------------------------------------------------
  ! output
  !-----------------------------------------------------------------------------
  WRITE(*,*) ' '
  WRITE(*,*) ' ---------------------------------------------------'
  WRITE(*,*) ' Parameters of first step: fitting of simulated data'
  WRITE(*,*) ' ---------------------------------------------------'

  DO i = 1, number_substances
     WRITE(*,*) ' '
     WRITE (*,*) ' Parameters of component ',trim(name_substances(i))
     WRITE (*,'(2(a,G18.10))') '    sigma-corr.: ', phi_ff_s(i),' epsilon-corr.:',phi_ff_e(i)
  END DO

  filename = 'gc.inp'
  filename = './input_file/'//filename
  CALL file_open(filename,77)

  WRITE(*,*) ' '
  WRITE(*,*) ' ---------------------------------------------------'
  WRITE(*,*) ' Second step: fitting force field to exp. data'
  WRITE(*,*) ' ---------------------------------------------------'

  read_loop2: DO
     READ (77,*) temp_name, temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling,  &
          m_fit_index, sig_fit_index, eps_fit_index
     IF ( temp_name == 'end' ) EXIT read_loop2

     WRITE(*,*) ' '
     WRITE (*,*) ' force field of the ',trim(temp_name),' group:'
     WRITE (*,*) ' --------------------------------'
     IF ( sig_fit_index /= 0 ) WRITE (*,*) '    sigma-parameter:  ',para(sig_fit_index)
     IF ( eps_fit_index /= 0 ) WRITE (*,*) '    epsilon-parameter:',para(eps_fit_index)
     IF ( polar == 1 ) WRITE (*,*) '    charge-parameter:  ',para(n_fit)
  END DO read_loop2
  CLOSE(77)

  !WRITE (*,*) para(1)
  !WRITE (*,*) para(2)
  !WRITE (*,*) 'CH2'
  !WRITE (*,*) para(3)
  !WRITE (*,*) para(4)

  DEALLOCATE( para, fvec, ipvt )

END SUBROUTINE fitting_forcefield


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE read_ff_weights

  USE FF_FIT_PARAMETERS
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, i0
  INTEGER                                :: i_left
  REAL                                   :: dummy(5)
  CHARACTER (LEN=50)                     :: filename

  CHARACTER (LEN=30)                     :: substance_name_1, substance_name_2
  !-----------------------------------------------------------------------------

  filename = './input_file/GEMC.dat'
  CALL file_open( filename, 22 )
  READ(22,*) substance_name_1,dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)
  IF ( index(substance_name_1,'_q') > 0 ) THEN    ! the charge is entered as part of the name
     i_left  = index( substance_name_1, '_q' )
     substance_name_1 = trim( adjustl( substance_name_1( 1: i_left-1 ) ) )
  END IF
  i = 1
  i0 = 1

  reading_data: DO

     READ(22,*) substance_name_2,dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)
     IF ( index(substance_name_2,'_q') > 0 ) THEN    ! the charge is entered as part of the name
        i_left  = index( substance_name_2, '_q' )
        substance_name_2 = trim( adjustl( substance_name_2( 1: i_left-1 ) ) )
     END IF

     IF ( substance_name_2 /= substance_name_1 ) THEN
        ! at this point all data-points from i0 to i of one substance are read.
        ! Either: substance_name_2 =='end' OR substance_name_2 points at a new substance
        filename = ' '
        filename = './input_file/ff_db/'//trim( adjustl( substance_name_1 ) )
        CALL file_open( filename, 23 )
        READ(23,*) dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)
        CLOSE (23)
        ff_weights(i0:i,1) = dummy(1)
        ff_weights(i0:i,2) = dummy(3)
        ff_weights(i0:i,3) = dummy(4)
        ff_weights(i0:i,4) = dummy(5)
        i0 = i + 1
     END IF

     IF ( substance_name_2 == 'end' ) EXIT reading_data

     i = i + 1                    ! i is the number of rows in data-file per substance
     substance_name_1 = substance_name_2

  END DO reading_data
  CLOSE (22)

END SUBROUTINE read_ff_weights




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE data_in_forcefield ( lm_m_dat, n_fit, para )

  USE BASIC_VARIABLES
  USE FF_FIT_PARAMETERS
  USE MOD_PARA_TRANSFER
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN OUT)                :: n_fit
  INTEGER, INTENT(OUT)                   :: lm_m_dat
  REAL, ALLOCATABLE, INTENT(OUT)         :: para(:)

  !-----------------------------------------------------------------------------
  REAL                                   :: dummy(10)
  CHARACTER (LEN=50)                     :: filename

  INTEGER                                :: i, k, kk
  INTEGER                                :: i_left, i_right
  CHARACTER (LEN=30)                     :: dum_name

  INTEGER                                :: m_fit_index, sig_fit_index, eps_fit_index
  REAL                                   :: temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling
  CHARACTER (LEN=20)                     :: temp_name
  REAL                                   :: para_temp(10)
  LOGICAL                                :: other_charge = .false.
  !-----------------------------------------------------------------------------


  DO i = 1, MMX
     plv(i) = 0.0
     rliq(i)= 0.0
     rvap(i)= 0.0
     dh_vap(i)  = 0.0
     IF ( fit_what == 'phi_to_sim' ) THEN
        tlv(i) = 0.0
        charge_val(i) = 0.0
        comp_name(i)  = ' '
        std_plv(i) = 1.0
        std_rliq(i)= 1.0
        std_rvap(i)= 1.0
        std_caloric(i)= 1.0
     END IF
  END DO

  dummy(10) = 0.0

  nlv = 0
  !-----------------------------------------------------------------------------
  ! simulation data together with all standard deviations is read
  ! from GEMC.dat ;experimental data from exp.dat
  !-----------------------------------------------------------------------------
  IF ( fit_what == 'phi_to_sim' ) THEN

     polar = 0

     number_substances = 1
     filename = 'GEMC.dat'
     filename = './input_file/'//filename
     CALL file_open(filename,22)
     reading_data1: DO
        READ(22,*) dum_name,dummy(1),dummy(2),dummy(3),dummy(4),dummy(5) !,dummy(6),dummy(7),dummy(8),dummy(9)
        IF ( index(dum_name,'_q') > 0 ) THEN
           polar = 1
           i_left  = index( dum_name, '_q' )
           i_right = index( dum_name, '_', .true. )
           read( dum_name(i_left+2:i_right-1), * )  dummy(10)       ! charge value
           dum_name = trim( adjustl( dum_name( 1: i_left-1 ) ) )
        END IF
        IF (dum_name == 'end') EXIT reading_data1
        nlv = nlv + 1                    ! nlv is the number of rows in GEMC.dat
        comp_name(nlv)   = dum_name
        plv(nlv)         = dummy(1)      ! vapour pressure (Pa)
        tlv(nlv)         = dummy(2)      ! temperature (K)
        rliq(nlv)        = dummy(3)      ! liquid density (kg/m**3)
        rvap(nlv)        = dummy(4)      ! vapour density (kg/m**3)
        dh_vap(nlv)      = dummy(5)      ! enthalpy of vaporization (J/kg?)
        ! std_plv(nlv)     = dummy(6)      ! standard deviation of the vapour pressure (fraction)
        ! std_rliq(nlv)    = dummy(7)      ! standard deviation of the liquid density (fraction)
        ! std_rvap(nlv)    = dummy(8)      ! standard deviation of the vapour density (fraction)
        ! std_caloric(nlv) = dummy(9)      ! standard deviation of h_LV (exp.dat) or u_LV (GEMC.dat) (fraction)
        charge_val(nlv)  = dummy(10)     ! characteristic (e.g. highest) charge of a polar substance

        IF (nlv == 1 ) THEN
           name_substances(1) = comp_name(nlv)
           T_low_substances(1) = tlv(nlv)
           other_charge = .false.
        ELSE              ! nlv is > 1
           IF ( comp_name(nlv) /= comp_name(nlv-1) ) THEN
              number_substances = number_substances + 1
              name_substances(number_substances) = comp_name(nlv)
              T_low_substances(number_substances) = tlv(nlv)
              other_charge = .false.
           ELSEIF ( charge_val(nlv) /= charge_val(nlv-1) ) THEN
              other_charge = .true.
              ! number_substances = number_substances + 1
              ! name_substances(number_substances) = comp_name(nlv)
              ! T_low_substances(number_substances) = tlv(nlv)
           END IF
        END IF
        IF ( .NOT. other_charge ) T_high_substances(number_substances) = tlv(nlv)
        comp_index(nlv) = number_substances
        write (*,*) '# of substances',number_substances,nlv, charge_val(nlv)

     END DO reading_data1

     CLOSE (22)

  ELSE IF ( fit_what == 'ffp_to_exp' ) THEN


     DO k = 1, number_substances
        write (*,*) k, T_low_substances(k),T_high_substances(k)
        ! read (*,*)
        filename = ' '
        filename = './input_file/ff_db/'//trim( adjustl( name_substances(k) ) )
        CALL file_open(filename,22)

        READ(22,*) dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)  ! first line specifies weights (here ignored)

        reading_data2: DO

           READ(22,*) dummy(1),dummy(2),dummy(3),dummy(4),dummy(5)
           ! write (*,*) k,nlv,dummy(2),T_low_substances(k)

           IF ( dummy(2) >= T_low_substances(k) ) THEN
              nlv = nlv + 1                    ! nlv is the number of rows in exp. data-file
              ! write (*,*) nlv,comp_index(nlv), dummy(2),comp_name(nlv)
              IF ( dummy(2) /= tlv(nlv) ) THEN
                 IF ( charge_val(nlv) /= charge_val(nlv-1) ) THEN
                    nlv = nlv - 1
                    kk = nlv
                    DO WHILE ( comp_index(kk) == comp_index(nlv) )
                       kk = kk + 1
                    END DO
                    tlv(nlv+1:MMX+nlv-kk+1)  = tlv(kk:MMX)
                    plv(nlv+1:MMX+nlv-kk+1)  = plv(kk:MMX)       ! vapour pressure (Pa)
                    rliq(nlv+1:MMX+nlv-kk+1) = rliq(kk:MMX)      ! liquid density (kg/m**3)
                    rvap(nlv+1:MMX+nlv-kk+1) = rvap(kk:MMX)      ! vapour density (kg/m**3)
                    dh_vap(nlv+1:MMX+nlv-kk+1) = dh_vap(kk:MMX)  ! enthalpy of vaporization (J/kg?)
                    misfit(nlv+1:MMX+nlv-kk+1,1) = misfit(kk:MMX,1)
                    misfit(nlv+1:MMX+nlv-kk+1,2) = misfit(kk:MMX,2)
                    misfit(nlv+1:MMX+nlv-kk+1,3) = misfit(kk:MMX,3)
                    misfit(nlv+1:MMX+nlv-kk+1,4) = misfit(kk:MMX,4)
                    comp_index(nlv+1:MMX+nlv-kk+1) = comp_index(kk:MMX)
                    comp_name(nlv+1:MMX+nlv-kk+1)  = comp_name(kk:MMX)
                    charge_val(nlv+1:MMX+nlv-kk+1) = charge_val(kk:MMX)
                    write (*,*) 'new array',nlv+1,tlv(nlv+1), comp_name(nlv+1)
                    EXIT reading_data2
                 ELSE
                    write (*,*) 'The data in file exp.dat does not match the ',  &
                         'simulated conditions',comp_name(nlv)
                    write (*,*) dummy(2), tlv(nlv),nlv
                    stop
                 END IF
              END IF
              plv(nlv)  = dummy(1)      ! vapour pressure (Pa)
              rliq(nlv) = dummy(3)      ! liquid density (kg/m**3)
              rvap(nlv) = dummy(4)      ! vapour density (kg/m**3)
              dh_vap(nlv) = dummy(5)      ! enthalpy of vaporization (J/kg?)
              IF ( dummy(2) == T_high_substances(k) )  THEN
                 IF ( charge_val(nlv) /= charge_val(nlv+1) .AND. comp_name(nlv) == comp_name(nlv+1) ) THEN
                    kk = nlv
                    DO WHILE ( comp_index(kk) == comp_index(nlv) )
                       kk = kk + 1
                    END DO
                    tlv(nlv+1:MMX+nlv-kk+1)  = tlv(kk:MMX)
                    plv(nlv+1:MMX+nlv-kk+1)  = plv(kk:MMX)       ! vapour pressure (Pa)
                    rliq(nlv+1:MMX+nlv-kk+1) = rliq(kk:MMX)      ! liquid density (kg/m**3)
                    rvap(nlv+1:MMX+nlv-kk+1) = rvap(kk:MMX)      ! vapour density (kg/m**3)
                    dh_vap(nlv+1:MMX+nlv-kk+1) = dh_vap(kk:MMX)  ! enthalpy of vaporization (J/kg?)
                    comp_index(nlv+1:MMX+nlv-kk+1) = comp_index(kk:MMX)
                    misfit(nlv+1:MMX+nlv-kk+1,1) = misfit(kk:MMX,1)
                    misfit(nlv+1:MMX+nlv-kk+1,2) = misfit(kk:MMX,2)
                    misfit(nlv+1:MMX+nlv-kk+1,3) = misfit(kk:MMX,3)
                    misfit(nlv+1:MMX+nlv-kk+1,4) = misfit(kk:MMX,4)
                    comp_name(nlv+1:MMX+nlv-kk+1)  = comp_name(kk:MMX)
                    charge_val(nlv+1:MMX+nlv-kk+1) = charge_val(kk:MMX)
                    write (*,*) 'new array',nlv+1,tlv(nlv+1), comp_name(nlv+1)
                 END IF
                 EXIT reading_data2
              END IF
           END IF

        END DO reading_data2

        CLOSE (22)
     END DO

  END IF
  ! write (*,*) 'total number of points',nlv
  ! read (*,*)

  !-----------------------------------------------------------------------------
  ! number of "data-points" considered
  ! ( =nlv * vapor pressure, nlv * liquid density,
  !    nlv * vapor density,  nlv * enthalpy of vaporiz. )
  !-----------------------------------------------------------------------------
  lm_m_dat = 4 * nlv



  !-----------------------------------------------------------------------------
  ! determine n_fit (the number of parameters that are adjusted)
  ! for the 2nd fitting step: set the initial values of parameters para(i)
  !-----------------------------------------------------------------------------

  IF ( fit_what == 'phi_to_sim' ) THEN

     n_fit = number_substances * 2
     IF ( polar == 1 ) n_fit = n_fit + number_substances

!!$     n_fit = n_fit + 1
!!$     ALLOCATE ( para_transfer( 1 ) )
!!$     call paus ('these 2 lines are temporary, delete soon')


  ELSE IF ( fit_what == 'ffp_to_exp' ) THEN

     n_fit = 0
     filename = 'gc.inp'
     filename = './input_file/'//filename
     CALL file_open(filename,77)

     reading_data3: DO
        READ (77,*) temp_name, temp_x_m, temp_x_sig, temp_x_eps, LJ_scaling, &
             m_fit_index, sig_fit_index, eps_fit_index
        IF ( temp_name == 'end' ) EXIT reading_data3
        IF (sig_fit_index /= 0) THEN
           n_fit = n_fit + 1
           para_temp(sig_fit_index) = temp_x_sig
        END IF
        IF (eps_fit_index /= 0) THEN
           n_fit = n_fit + 1
           para_temp(eps_fit_index) = temp_x_eps
        END IF
     END DO reading_data3
     CLOSE(77)
     IF ( polar == 1 ) n_fit = n_fit + 1
     IF ( polar == 1 ) para_temp(n_fit) = charge_val(1)

  END IF

  !-----------------------------------------------------------------------------
  ! allocate the vector para(:)
  !-----------------------------------------------------------------------------
  ALLOCATE( para(n_fit) )
  IF ( fit_what == 'ffp_to_exp' ) para(1:n_fit) = para_temp(1:n_fit)

END SUBROUTINE data_in_forcefield




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE feed_gcmc_simulation
!
! This subroutine calculates properties needed for a multiple histogram GCMC
! simulation. It calculates the chemical potential and temperatures at several
! conditions tracing the vapor liquid envelope of a pure substance.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE feed_gcmc_simulation

  USE PARAMETERS, ONLY: PI, KBOL
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  USE DFT_MODULE, ONLY: fres_temp
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, crit_dat, converg, calc_option
  INTEGER                                :: loop, t_steps, dividing_n(10), n_max(10), n_l(10)
  REAL                                   :: rhob(2), my0, zges(np), end_x, steps, pbulk
  REAL                                   :: tc, pc, rhoc, temp1(10,0:8000), temp2(0:8000)
  REAL                                   :: lnphi(np,nc)
  REAL                                   :: rhopt, box_l, fres_l, fres_v, t_stage(10)
  REAL                                   :: t_ink_stage(10), mu_stage(10), min_distance, volume
  REAL                                   :: v_dum, plv_dum
  REAL                                   :: teilung
  CHARACTER (LEN=4)                      :: char_len
  !-----------------------------------------------------------------------------
  ! note: variable    ensemble_flag = 'tp'    is set for calculating: my_res=my_res(T,p)
  ! note: variable    ensemble_flag = 'tv'    is set for calculating: my_res=my_res(T,rho)


  write (*,*) 'Choose your calculation option:'
  write (*,*) '(0) fit bias potential to simulation results'
  write (*,*) '(1) Calculate chemical potential & bias potential to a single T: pure comp.'
  write (*,*) '(2) Calculate chemical potential & bias potential to a single T: mixtures'
  write (*,*) '(3) Calculate input for multiple histogram VLE-method'
  write (*,*) '(4) Calculate input for multiple histogram isotherms'
  read (*,*) calc_option
  IF (calc_option == 0) CALL match_bias_potential
  IF (calc_option == 1) CALL VLE_bias_potential
  IF (calc_option == 2) CALL VLE_bias_potential_mix
  IF (calc_option == 4) CALL isotherm_bias_potential
  IF (calc_option /= 3) RETURN


  !-----------------------------------------------------------------------------
  ! define number of GCMC steps (T_steps) and reduced Temp. T/Tc
  !-----------------------------------------------------------------------------

  t_steps = 10                  ! number of GCMC conditions

  teilung = 5.0                 ! T-spacing in unit [K] for defining the temperature conditions

  !-----------------------------------------------------------------------------
  ! the first two values are for a vapor phase !
  !-----------------------------------------------------------------------------
  t_stage(1) = 0.80             ! reduced temperature (final Kelvin-value is set divisible by 10 k)
  t_stage(2) = 0.90
  t_stage(3) = 0.97
  t_stage(4) = 0.90
  t_stage(5) = 0.84
  t_stage(6) = 0.78
  t_stage(7) = 0.72
  t_stage(8) = 0.67
  t_stage(9) = 0.62
  t_stage(10)= 0.57


  !-----------------------------------------------------------------------------
  ! read substance name and input specifications
  !-----------------------------------------------------------------------------

  WRITE(*,*) 'which substance?'
  WRITE(*,*) '(for available substances, see para_input.f)'
  READ(*,*) compna(1)

  WRITE (*,*) 'enter an estimate for crit. temp.'
  READ (*,*) tc

  WRITE(*,*) 'specify volume of GCMC simulation box, [Angstrom**3]'
  READ(*,*) volume
  box_l = volume**(1.0/3.0)

  !-----------------------------------------------------------------------------
  ! draw pure component parameters from "database"
  !-----------------------------------------------------------------------------

  ncomp = 1
  num = 1
  eos = 1
  pol = 1
  CALL SET_DEFAULT_EOS_NUMERICAL

  CALL para_input            ! retriev pure comp. parameters

  !z3t = PI/6.0* parame(1,1) * ( 1.0  &
  !    *(1.0-0.12*EXP(-3.0*parame(1,3)/t)) )**3    ! divided by parame(i,2)**3


  !-----------------------------------------------------------------------------
  ! prepare phase equilibrium calculation
  !-----------------------------------------------------------------------------

  nphas = 2
  n_unkw = ncomp                ! number of quantities to be iterated
  it(1) = 'lnp'                 ! iteration of pressure

  running='t'                   ! T is running variable in PHASE_EQUILIB - here T=const.
  end_x  = t                    ! end temperature - here T=const.
  steps = 1.0                   ! number of steps of running var. - here 1, since T=const.

  outp = 0                      ! output to terminal

  OPEN (68,FILE = './output_file/GCMC_conditions.dat')

  !-----------------------------------------------------------------------------
  ! calculate the critical temperature (according to PC-SAFT)
  !-----------------------------------------------------------------------------

  ensemble_flag = 'tp'
  CALL critical (tc,pc,rhoc)
  WRITE (*,'(a,3(f16.4))') 'critical point',tc,pc/1.d5,rhoc
  WRITE (*,*) ' '

  !=============================================================================
  ! LOOP AROUND ALL TEMPERATURES
  !=============================================================================

  DO loop = 1, t_steps

     t = REAL( NINT( t_stage(loop)/teilung * tc ) )*teilung   ! set temperature divisible by 10K
     val_init(3) = t
     t_ink_stage(loop) = t
     end_x = t

     ensemble_flag = 'tp'    ! is set for 'regular calculations'
     !z3t = PI/6.0* parame(1,1) * ( 1.0  &
     !    *(1.0-0.12*EXP(-3.0*parame(1,3)/t)) )**3   ! divided by parame(i,2)**3


     !--------------------------------------------------------------------------
     ! calculate phase equilibrium for given T
     !--------------------------------------------------------------------------

     !  find_equilibrium: DO
     nphas = 2
     val_init(1) = 0.45                  ! starting value for liquid density
     val_init(2) = 1.E-5                 ! starting value for vapor density
     val_init(3) = t                     ! value of temperature       NOTE: value assigned below!
     val_init(4) = 1.E-3                  ! starting value for p in [Pa]
     val_init(5) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 1
     val_init(6) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 2

     crit_dat = 1
     CALL pure_equilibrium_fit ( crit_dat, tc, t, 1.E-1, converg, v_dum, plv_dum )
     ! write (*,*) converg, t,p, val_conv(3:4)
     ! call paus (' ')

     IF (converg /= 1) THEN
        write (*,*) 'no phase equilibrium found for T=',t
     END IF


     !--------------------------------------------------------------------------
     ! coexisting molecular density in unit (1 / Angstrom**3)
     !--------------------------------------------------------------------------
     rhob(1) = dense(1) / z3t           ! coexisting bulk density L
     rhob(2) = dense(2) / z3t           ! coexisting bulk density V
     WRITE (*,*) 'temperature     ',t, p
     WRITE (*,*) 'densities       ',rhob(1), rhob(2)

     !--------------------------------------------------------------------------
     ! (re-)calculate resid. chem. potential (but for given T,rho)
     !--------------------------------------------------------------------------
     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     densta(1) = dense(1)
     densta(2) = dense(2)
     CALL fugacity (lnphi)
     my0 = lnphi(1,1) + LOG(rhob(1))   ! my0 = my_res(T,rho_bulk_L) + ln(rho_bulk_l)
     mu_stage(loop) = my0
     zges(1) = p/KBOL/t/rhob(1)/1.E30
     zges(2) = p/KBOL/t/rhob(2)/1.E30

     pbulk = zges(1)*rhob(1)     ! pressure  p/kT (= Z*rho)   in unit (1 / Angstrom**3)
     WRITE (*,*) 'chem. potentials', my0
     WRITE (*,*) ' '



     !--------------------------------------------------------------------------
     ! (re-)calculate the Helmholtz energy of vap. and liq. phase
     !--------------------------------------------------------------------------

     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     dense(1)  = rhob(1)*z3t
     densta(1) = rhob(1)*z3t
     nphas = 1
     CALL fugacity (lnphi)
     fres_l = fres_temp*rhob(1) + rhob(1)*(LOG(rhob(1))-1.0)
     fres_l = fres_l * parame(1,2)**3
     dense(1)  = rhob(2)*z3t
     densta(1) = rhob(2)*z3t
     CALL fugacity (lnphi)
     fres_v = fres_temp*rhob(2) + rhob(2)*(LOG(rhob(2))-1.0)
     fres_v = fres_v * parame(1,2)**3

     !--------------------------------------------------------------------------
     ! calculate an estimate of bias potential between V and L
     !--------------------------------------------------------------------------

     n_max(loop) = NINT( rhob(1)*1.4*box_l**3 )
     n_l(loop)   = NINT( rhob(1)*box_l**3 )

     ! --- loop = 1 is for the vapor phase --------------------------------
     IF (loop <= 2) THEN
        n_max(loop) = NINT( rhob(1)*0.3*box_l**3 )
        n_l(loop)   = NINT( rhob(2)*box_l**3 )
     ELSE IF (loop == 3) THEN    ! this condition is also for the first (vapor) condition.
        n_l(loop)   = NINT( rhob(2)*box_l**3 )
     END IF

     DO i = 0, n_max(loop)
        rhopt = REAL(i) / box_l**3
        IF (rhopt == 0.0) rhopt = 1.d-5
        dense(1)  = rhopt * z3t
        densta(1) = rhopt * z3t
        ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
        CALL fugacity (lnphi)
        temp2(i) = fres_temp*rhopt + rhopt*( LOG(rhopt) - 1.0 )
        temp2(i) = temp2(i) * parame(1,2)**3

        temp1(loop,i)= -(  fres_l + (fres_v-fres_l)/(rhob(2)-rhob(1))*(rhopt-rhob(1))  &
             - temp2(i) )/2.0 *box_l**2
        WRITE(68,'(i4,G18.10)') i, temp1(loop,i)
     END DO
     WRITE (68,*) ' '

  END DO

  !-----------------------------------------------------------------------------
  ! determine a suitable division of N-space for each condition
  !-----------------------------------------------------------------------------

  DO loop = 2, t_steps
     min_distance = 1.E10
     ! write (*,*) loop,N_L(loop-1),N_L(loop)
     DO i = n_l(loop-1), n_l(loop)
        !IF (loop == 3) WRITE (*,*) ABS( temp1(loop-1,i)-temp1(loop,i) ), min_distance
        IF ( ABS( temp1(loop-1,i)-temp1(loop,i) ) < min_distance ) THEN
           min_distance = ABS( temp1(loop-1,i)-temp1(loop,i) )
           dividing_n(loop-1) = i
        END IF
     END DO
     n_max(loop-1) = dividing_n(loop-1)
     ! write (*,*) 'fff',loop-1,N_max(loop-1)
     IF ( loop == t_steps ) THEN
        i = n_l(loop)
        DO WHILE ( temp1(loop,i) < 10.0 )
           i = i + 1
           n_max(loop) = i
        END DO
        ! write (*,*) 'fff',loop,N_max(loop)
     END IF
  END DO

  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' ------ bias potential for temperatures given below -'
  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' '


  DO i = 0, n_max(1)
     WRITE (*,*) i, temp1(1,i)
  END DO
  DO loop = 2, t_steps
     DO i = n_max(loop-1), n_max(loop)
        WRITE (*,*) i, temp1(loop,i)
     END DO
  END DO
  WRITE (*,*) ' '


  WRITE (*,*) ' '
  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' ------ proposed conditions for GCMC simulation -----'
  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' '
  WRITE (*,'(a,I6,a)') 'foreach volume (',NINT(box_l**3),' )'
  WRITE (*,'(a,G10.1,a)') '@  counter = $counter + 1'
  WRITE (*,*) ' '
  WRITE (*,'(a)') 'if ( $counter == 1 ) then'
  WRITE (*,'(a,I4)') '@ maxmol = ', n_max(t_steps)
  WRITE (*,'(a)') 'endif'
  WRITE (*,'(a)') 'if ( $counter == 2 ) then'
  WRITE (*,'(a)') '@ maxmol = ...'
  WRITE (*,'(a)') 'endif'
  WRITE (*,*) ' '
  WRITE (char_len,'(I3)') T_steps
  WRITE(*,'(a,'//char_len//'G10.5,a)') 'set  temp      = ( ',(t_ink_stage(loop),loop=1,t_steps),' )'
  WRITE(*,'(a,'//char_len//'G10.4,a)') 'set  mu_sim    = ( ',(mu_stage(loop),loop=1,t_steps),' )'
  WRITE(*,'(a,'//char_len//'(4x,I6),a)') 'set  N_step1   = (    '  &
       ,(n_max(loop-1),loop=2,t_steps),100000, ' )   # last value is dummy'
  WRITE (*,*) ' '

END SUBROUTINE feed_gcmc_simulation


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE VLE_bias_potential
!
! This subroutine calculates the chemical potential as well as a biasing
! function (allowing a simple transition between liquid and vapor phases)
! for parameterizing GCMC simulations.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE VLE_bias_potential

  USE PARAMETERS, ONLY: PI, KBOL
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  use utilities, only: critical
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, crit_dat, converg
  INTEGER                                :: n_max
  REAL                                   :: temp
  REAL                                   :: rhob(2), end_x, steps, pbulk
  REAL                                   :: tc, pc, rhoc
  REAL                                   :: f_mean, F_bias(0:8000), f_dft(0:8000)
  REAL                                   :: lnphi(np,nc)
  REAL                                   :: rhopt, box_l, fres_l, fres_v
  REAL                                   :: mu_coex, volume
  REAL                                   :: v_dum, plv_dum
  !-----------------------------------------------------------------------------
  ! note: variable    ensemble_flag = 'tp'    is set for calculating: my_res=my_res(T,p)
  ! note: variable    ensemble_flag = 'tv'    is set for calculating: my_res=my_res(T,rho)


  !-----------------------------------------------------------------------------
  ! read substance name and input specifications
  !-----------------------------------------------------------------------------

  WRITE(*,*) 'which substance?'
  WRITE(*,*) '(for available substances, see para_input.f)'
  READ(*,*) compna(1)

  WRITE(*,*) 'specify Temperature, [K]'
  READ(*,*) temp

  WRITE (*,*) 'enter an estimate for crit. temp.'
  READ (*,*) tc

  WRITE(*,*) 'specify volume of GCMC simulation box, [Angstrom**3]'
  READ(*,*) volume
  box_l = volume**(1.0/3.0)

  !-----------------------------------------------------------------------------
  ! draw pure component parameters from "database"
  !-----------------------------------------------------------------------------

  ncomp = 1
  num = 1
  eos = 1
  pol = 1
  CALL SET_DEFAULT_EOS_NUMERICAL

  CALL para_input            ! retriev pure comp. parameters


  !-----------------------------------------------------------------------------
  ! prepare phase equilibrium calculation
  !-----------------------------------------------------------------------------

  nphas = 2
  n_unkw = ncomp                ! number of quantities to be iterated
  it(1) = 'lnp'                 ! iteration of pressure

  running='t'                   ! T is running variable in PHASE_EQUILIB - here T=const.
  end_x  = t                    ! end temperature - here T=const.
  steps = 1.0                   ! number of steps of running var. - here 1, since T=const.

  outp = 0                      ! output to terminal

  OPEN (68,FILE = './output_file/GCMC_bias.dat')

  !-----------------------------------------------------------------------------
  ! calculate the critical temperature (according to PC-SAFT)
  !-----------------------------------------------------------------------------

  ensemble_flag = 'tp'
  CALL critical (tc,pc,rhoc)
  WRITE (*,'(a,3(f16.4))') 'critical point',tc,pc/1.d5,rhoc
  WRITE (*,*) ' '

  ensemble_flag = 'tp'    ! is set for 'regular calculations'


  !-----------------------------------------------------------------------------
  ! calculate phase equilibrium for given T
  !-----------------------------------------------------------------------------

  nphas = 2
  val_init(1) = 0.45                  ! starting value for liquid density
  val_init(2) = 1.E-5                 ! starting value for vapor density
  val_init(3) = temp                  ! value of temperature       NOTE: value assigned below!
  val_init(4) = 1.E4                  ! starting value for p in [Pa]
  val_init(5) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 1
  val_init(6) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 2

  end_x = temp

  crit_dat = 1
  CALL pure_equilibrium_fit ( crit_dat, tc, t, 1.E4, converg, v_dum, plv_dum )

  IF (converg /= 1) THEN
     write (*,*) 'no phase equilibrium found for T=',t
  END IF


  !-----------------------------------------------------------------------------
  ! coexisting molecular density in unit (1 / Angstrom**3)
  !-----------------------------------------------------------------------------
  rhob(1) = dense(1) / z3t           ! coexisting bulk density L
  rhob(2) = dense(2) / z3t           ! coexisting bulk density V
  WRITE (*,*) 'temperature     ',t, p
  WRITE (*,*) 'densities       ',rhob(1), rhob(2)

  mu_coex = my_cal( 1, 1 )   ! mu_coex = my_res(T,rho_bulk_L) + ln(rho_bulk_l)
  pbulk = z_cal(1)*rhob(1)     ! pressure  p/kT (= Z*rho)   in unit (1 / Angstrom**3)
  fres_l = f_res(1) !* parame(1,2)**3
  fres_v = f_res(2) !* parame(1,2)**3

  WRITE (*,*) 'chem. potentials', mu_coex
  WRITE (*,*) ' '
  WRITE (68,'(a)') ' T      and       mu_coex(chem.potential.at.coexistence)'
  WRITE (68,'(2G18.10)') t, mu_coex
  WRITE (68,'(a)') '-------------------------------------------------'
  WRITE (68,*) ' '



  !-----------------------------------------------------------------------------
  ! calculate an estimate of bias potential between V and L
  !-----------------------------------------------------------------------------

  n_max = NINT( rhob(1)*1.4*volume )

  nphas = 1
  write (*,*) n_max, rhob(1),rhob(2),box_l

  DO i = 0, n_max
     rhopt = REAL(i) / volume
     IF (rhopt == 0.0) rhopt = 1.d-5
     dense(1)  = rhopt * z3t
     densta(1) = rhopt * z3t
     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     CALL fugacity (lnphi)
     f_dft(i) = f_res(1)
     !temp2(i) = temp2(i) * parame(1,2)**3

     f_mean = fres_l + (fres_v-fres_l)/(rhob(2)-rhob(1))*(rhopt-rhob(1))
     F_bias(i)= ( f_dft(i) - f_mean ) * box_l**3   !/2.0 *box_l**2
     WRITE(68,'(i4,G18.10)') i, F_bias(i)
  END DO
  WRITE (68,*) ' '


  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' bias potential written to file: output_file/GCMC_bias.dat'
  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' '

END SUBROUTINE VLE_bias_potential


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE VLE_bias_potential
!
! This subroutine calculates the chemical potential as well as a biasing
! function (allowing a simple transition between liquid and vapor phases)
! for parameterizing GCMC simulations.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE VLE_bias_potential_mix

  USE PARAMETERS, ONLY: PI, KBOL
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  USE STARTING_VALUES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, converg
  INTEGER                                :: n_max
  REAL                                   :: rhob(2)
  REAL                                   :: f_mean, F_bias(0:8000), f_dft(0:8000)
  REAL                                   :: lnphi(np,nc)
  REAL                                   :: rhopt, box_l, fres_l, fres_v
  REAL                                   :: volume
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! call input specifications
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2

  WRITE(*,*) 'specify volume of GCMC simulation box, [Angstrom**3]'
  READ(*,*) volume
  box_l = volume**(1.0/3.0)

  !-----------------------------------------------------------------------------
  ! phase equilibrium calculation
  !-----------------------------------------------------------------------------

  CALL START_VAR(converg)     ! gets starting values, sets "val_init"
  CALL objective_ctrl ( converg )

  OPEN (68,FILE = './output_file/GCMC_bias.dat')

  WRITE (*,*) 'chem. potentials comp. 1', my_cal( 1:2, 1 )
  WRITE (*,*) 'chem. potentials comp. 2', my_cal( 1:2, 2 )
  WRITE (*,*) ' '
  WRITE (68,'(a)') ' T,      p      and       mu_coex(chem.potential.at.coexistence)'
  WRITE (68,'(4G18.10)') t, p,  my_cal( 1, 1:2 )
  WRITE (68,'(a, 2G18.10)') 'N1,N2 phase_1', rhoi_cal(1,1:ncomp)*volume
  WRITE (68,'(a, 2G18.10)') 'N1,N2 phase_2', rhoi_cal(2,1:ncomp)*volume
  WRITE (68,'(a)') '-------------------------------------------------'
  WRITE (68,*) ' '
  fres_l = f_res(1)
  fres_v = f_res(2)

  write (*,*) 'finish the part below !!!!!'
  stop

  !-----------------------------------------------------------------------------
  ! calculate an estimate of bias potential between V and L
  !-----------------------------------------------------------------------------

  nphas = 1

  n_max = NINT( rhob(1) * 1.4 * volume )
  write (*,*) n_max, rhob(1),rhob(2),box_l

  DO i = 0, n_max
     rhopt = REAL(i) / volume
     IF (rhopt == 0.0) rhopt = 1.d-5
     dense(1)  = rhopt * z3t
     densta(1) = rhopt * z3t
     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     CALL fugacity (lnphi)
     f_dft(i) = f_res(1)
     !temp2(i) = temp2(i) * parame(1,2)**3

     f_mean = fres_l + (fres_v-fres_l)/(rhob(2)-rhob(1))*(rhopt-rhob(1))
     F_bias(i)= ( f_dft(i) - f_mean ) * box_l**3   !/2.0 *box_l**2
     WRITE(68,'(i4,G18.10)') i, F_bias(i)
  END DO
  WRITE (68,*) ' '

  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' bias potential written to file: output_file/GCMC_bias.dat'
  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' '

END SUBROUTINE VLE_bias_potential_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE isotherm_bias_potential
!
! This subroutine calculates the chemical potential as well as a biasing
! function (allowing a simple transition between liquid and vapor phases)
! for parameterizing GCMC simulations.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine isotherm_bias_potential

  USE PARAMETERS, ONLY: Nav
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t, x, rho, pges, mseg
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  integer                                :: i
  integer                                :: dens_option
  integer                                :: N_max_isotherm
  real                                   :: lnphi(np,nc)
  real                                   :: volume
  real                                   :: rho_max, eta_max
  real, dimension(5000)                  :: mu_isotherm

  integer                                :: N_index, N_division, loop
  integer, dimension(0:12)               :: N_step
  real, dimension(12)                    :: N_step_fraction, t_step, mu_step
  real                                   :: delta_window
  character (LEN=4)                      :: char_len
  !-----------------------------------------------------------------------------
  ! note: variable    ensemble_flag = 'tp'    is set for calculating: my_res=my_res(T,p)
  ! note: variable    ensemble_flag = 'tv'    is set for calculating: my_res=my_res(T,rho)


  !-----------------------------------------------------------------------------
  ! read substance name and input specifications
  !-----------------------------------------------------------------------------

  write (*,*) 'which substance?'
  write (*,*) '(for available substances, see para_input.f)'
  read (*,*) compna(1)

  write (*,*) 'specify Temperature, [K]'
  read (*,*) t

  write (*,*) 'specify volume of GCMC simulation box, [Angstrom**3]'
  read (*,*) volume

  !-----------------------------------------------------------------------------
  ! draw pure component parameters from "database"
  !-----------------------------------------------------------------------------

  ncomp = 1
  num = 1
  eos = 1
  pol = 1
  call SET_DEFAULT_EOS_NUMERICAL

  call para_input            ! retriev pure comp. parameters

  nphas = 1
  ensemble_flag = 'tv'
  outp = 0                      ! output to terminal
  u_out_t = 0.0
  u_out_p = 1.0

  xi(1,1) = 1.0

  x(1) = xi(1,1)
  call PERTURBATION_PARAMETER


  write (*,*) ' In what unit will the maximum density be given:'
  write (*,*) ' (1) in [kg/m^3]'
  write (*,*) ' (2) in [mol/m^3]'
  write (*,*) ' (3) dimensionless'
  read (*,*) dens_option

  if ( dens_option == 1 ) then
     write (*,*) ' maximum density in [kg/m^3]'
     read (*,*) rho_max
     rho_max = rho_max / (  mm(1) * 1.E27 / Nav   )
  else if ( dens_option == 2 ) then
     write (*,*) ' maximum density in [mol/m^3]'
     read (*,*) rho_max
     rho_max = rho_max / (  1.E30 / Nav   )
  else
     write (*,*) ' maximum dimensionless segment density rho*, in [-]'
     read (*,*) rho_max
     rho_max = rho_max / mseg(1)
  end if

  eta_max = rho_max * z3t

  N_max_isotherm = INT( rho_max * volume ) + 1

  do i = 1, N_max_isotherm

     densta(1) = REAL(i) / real(N_max_isotherm) * eta_max
     dense(1) = densta(1)

     CALL FUGACITY (lnphi)
     mu_isotherm(i) = lnphi(1,1) + log( rho )
     write (*,'(i5,4G18.8)') i, mu_isotherm(i), pges/1.E5, rho, rho * (  mm(1) * 1.E27 / Nav   )

  end do

  write (*,*) ' maximum number of molecules:', N_max_isotherm
  write (*,*) ' define window width:'
  read (*,*) delta_window


  !N_division = 8
  !N_step_fraction(1) = 0.30
  !N_step_fraction(2) = 0.47
  !N_step_fraction(3) = 0.60
  !N_step_fraction(4) = 0.72
  !N_step_fraction(5) = 0.82
  !N_step_fraction(6) = 0.90
  !N_step_fraction(7) = 0.95
  !N_step_fraction(8) = 1.00
  N_division = 11
  N_step_fraction(1) = 0.26
  N_step_fraction(2) = 0.45
  N_step_fraction(3) = 0.58
  N_step_fraction(4) = 0.68
  N_step_fraction(5) = 0.76
  N_step_fraction(6) = 0.82
  N_step_fraction(7) = 0.87
  N_step_fraction(8) = 0.91
  N_step_fraction(9) = 0.94
  N_step_fraction(10) = 0.97
  N_step_fraction(11) = 1.00

  do i = 1, N_division
     t_step(i) = t
     N_step(i) = delta_window * NINT( N_step_fraction(i) * real(N_max_isotherm) / real(delta_window) )
  end do

  N_step(0) = 0
  do i = 0, N_division - 1
     N_index = NINT( 1.0/3.0 * real( N_step(i) ) + 2.0/3.0 * real( N_step(i+1) ) )
     mu_step( i+1 ) = mu_isotherm( N_index )
     write (*,*) N_step(i), N_index, mu_isotherm( N_index )
  end do


  WRITE (*,*) ' '
  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' ------ proposed conditions for GCMC simulation -----'
  WRITE(*,*)' ----------------------------------------------------'
  WRITE(*,*)' '
  WRITE (*,'(a,I6,a)') 'foreach volume (',NINT(volume),' )'
  WRITE (*,'(a,G10.1,a)') '@  counter = $counter + 1'
  WRITE (*,*) ' '
  WRITE (*,'(a)') 'if ( $counter == 1 ) then'
  WRITE (*,'(a,I4)') '@ maxmol = ', N_step( N_division )
  WRITE (*,'(a)') 'endif'
  WRITE (*,'(a)') 'if ( $counter == 2 ) then'
  WRITE (*,'(a)') '@ maxmol = ...'
  WRITE (*,'(a)') 'endif'
  WRITE (*,*) ' '
  WRITE (char_len,'(I3)') N_division
  WRITE(*,'(a,'//char_len//'G10.5,a)') 'set  temp      = ( ',(t_step(loop),loop=1,N_division),' )'
  WRITE(*,'(a,'//char_len//'G10.4,a)') 'set  mu_sim    = ( ',(mu_step(loop),loop=1,N_division),' )'
  WRITE(*,'(a,'//char_len//'(4x,I6),a)') 'set  N_step1   = (    '  &
       ,(N_step(loop-1),loop=2,N_division),100000, ' )   # last value is dummy'
  WRITE (*,*) ' '


end subroutine isotherm_bias_potential


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE match_fit_parameters
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module match_fit_parameters

  implicit none
  save
  REAL                                   :: rhob(2)
  REAL                                   :: fres_l, fres_v
  REAL                                   :: box_l, normalize_sigma

End Module match_fit_parameters



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE match_bias_potential
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE match_bias_potential

  USE PARAMETERS, ONLY: PI, KBOL
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  USE DFT_MODULE, ONLY: fres_temp
  USE Levenberg_Marquardt
  USE match_fit_parameters
  use utilities, only: file_open, CRITICAL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE match_bias_obj (lm_m_dat, lm_n_par, para, fvec, iflag)
       IMPLICIT NONE
       INTEGER, INTENT(IN)                 :: lm_m_dat, lm_n_par
       REAL, INTENT(IN)                    :: para(:)
       REAL, INTENT(IN OUT)                :: fvec(:)
       INTEGER, INTENT(IN OUT)             :: iflag
     END SUBROUTINE match_bias_obj
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                :: crit_dat, converg
  INTEGER                                :: n_max
  REAL                                   :: temp
  REAL                                   :: zges(np), end_x, steps, pbulk
  REAL                                   :: tc, pc, rhoc
  REAL                                   :: lnphi(np,nc)
  REAL                                   :: mu_coex, volume
  REAL                                   :: v_dum, plv_dum
  CHARACTER (LEN=50)                     :: filename
  !-----------------------------------------------------------------------------
  INTEGER                                :: info
  INTEGER                                :: lm_n_par         ! number of adjustable parameters
  INTEGER                                :: lm_m_dat         ! number of data points
  INTEGER, ALLOCATABLE                   :: ipvt(:)
  REAL, ALLOCATABLE                      :: para(:)
  REAL, ALLOCATABLE                      :: fvec(:)
  REAL                                   :: tol    = 1.0E-8
  REAL                                   :: epsfcn
  REAL                                   :: lm_factor
  !-----------------------------------------------------------------------------
  ! note: variable    ensemble_flag = 'tp'    is set for calculating: my_res=my_res(T,p)
  ! note: variable    ensemble_flag = 'tv'    is set for calculating: my_res=my_res(T,rho)


  !-----------------------------------------------------------------------------
  ! read substance name and input specifications
  !-----------------------------------------------------------------------------

  WRITE(*,*) 'which substance?'
  WRITE(*,*) '(for available substances, see para_input.f)'
  !READ(*,*) compna(1)
  compna(1) = 'ethane'

  WRITE(*,*) 'specify Temperature, [K]'
  !READ(*,*) temp
  temp = 183.0

  WRITE (*,*) 'enter an estimate for crit. temp.'
  !READ (*,*) tc
  tc = 300.0

  WRITE(*,*) 'specify volume of GCMC simulation box, [Angstrom**3]'
  !READ(*,*) volume
  volume = 100000.0
  box_l = volume**(1.0/3.0)

  !-----------------------------------------------------------------------------
  ! draw pure component parameters from "database"
  !-----------------------------------------------------------------------------

  ncomp = 1
  num = 1
  eos = 1
  pol = 1
  CALL SET_DEFAULT_EOS_NUMERICAL

  CALL para_input            ! retriev pure comp. parameters


  !-----------------------------------------------------------------------------
  ! prepare phase equilibrium calculation
  !-----------------------------------------------------------------------------

  nphas = 2
  n_unkw = ncomp                ! number of quantities to be iterated
  it(1) = 'lnp'                 ! iteration of pressure

  running='t'                   ! T is running variable in PHASE_EQUILIB - here T=const.
  end_x  = t                    ! end temperature - here T=const.
  steps = 1.0                   ! number of steps of running var. - here 1, since T=const.

  outp = 0                      ! output to terminal

  OPEN (68,FILE = './output_file/GCMC_bias.dat')

  !-----------------------------------------------------------------------------
  ! calculate the critical temperature (according to PC-SAFT)
  !-----------------------------------------------------------------------------

  ensemble_flag = 'tp'
  CALL critical (tc,pc,rhoc)
  WRITE (*,'(a,3(f16.4))') 'critical point',tc,pc/1.d5,rhoc
  WRITE (*,*) ' '

  ensemble_flag = 'tp'    ! is set for 'regular calculations'


  !-----------------------------------------------------------------------------
  ! calculate phase equilibrium for given T
  !-----------------------------------------------------------------------------

  nphas = 2
  val_init(1) = 0.45                  ! starting value for liquid density
  val_init(2) = 1.E-5                 ! starting value for vapor density
  val_init(3) = temp                  ! value of temperature       NOTE: value assigned below!
  val_init(4) = 1.E4                  ! starting value for p in [Pa]
  val_init(5) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 1
  val_init(6) = 0.0                   ! logar. mole fraction: lnx=0, x=1, phase 2

  end_x = temp

  crit_dat = 1
  CALL pure_equilibrium_fit ( crit_dat, tc, t, 1.E4, converg, v_dum, plv_dum )

  IF (converg /= 1) THEN
     write (*,*) 'no phase equilibrium found for T=',t
  END IF


  !-----------------------------------------------------------------------------
  ! coexisting molecular density in unit (1 / Angstrom**3)
  !-----------------------------------------------------------------------------
  rhob(1) = dense(1) / z3t           ! coexisting bulk density L
  rhob(2) = dense(2) / z3t           ! coexisting bulk density V
  WRITE (*,*) 'temperature     ',t, p
  WRITE (*,*) 'densities       ',rhob(1), rhob(2)

  !-----------------------------------------------------------------------------
  ! (re-)calculate resid. chem. potential (but for given T,rho)
  !-----------------------------------------------------------------------------
  ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
  densta(1) = dense(1)
  densta(2) = dense(2)
  CALL fugacity (lnphi)
  mu_coex = lnphi(1,1) + LOG(rhob(1))   ! mu_coex = my_res(T,rho_bulk_L) + ln(rho_bulk_l)
  zges(1) = p/KBOL/t/rhob(1)/1.E30
  zges(2) = p/KBOL/t/rhob(2)/1.E30

  pbulk = zges(1)*rhob(1)     ! pressure  p/kT (= Z*rho)   in unit (1 / Angstrom**3)
  WRITE (*,*) 'chem. potentials', mu_coex
  WRITE (*,*) ' '
  WRITE (68,'(a)') ' T      and       mu_coex(chem.potential.at.coexistence)'
  WRITE (68,'(2G18.10)') t, mu_coex
  WRITE (68,'(a)') '-------------------------------------------------'
  WRITE (68,*) ' '



  !-----------------------------------------------------------------------------
  ! (re-)calculate the Helmholtz energy of vap. and liq. phase
  !-----------------------------------------------------------------------------

  ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
  dense(1)  = rhob(1)*z3t
  densta(1) = rhob(1)*z3t
  nphas = 1
  CALL fugacity (lnphi)
  fres_l = fres_temp*rhob(1) + rhob(1)*(LOG(rhob(1))-1.0)
  fres_l = fres_l * parame(1,2)**3
  dense(1)  = rhob(2)*z3t
  densta(1) = rhob(2)*z3t
  CALL fugacity (lnphi)
  fres_v = fres_temp*rhob(2) + rhob(2)*(LOG(rhob(2))-1.0)
  fres_v = fres_v * parame(1,2)**3

  !-----------------------------------------------------------------------------
  ! calculate an estimate of bias potential between V and L
  !-----------------------------------------------------------------------------

  n_max = NINT( rhob(1)*1.4*box_l**3 )
  write (*,*) n_max, rhob(1),rhob(2),box_l

  filename = './input_file/bias_pot_sim.dat'
  CALL file_open( filename, 67 )
  read (67,*) lm_m_dat
  rewind( 67 )

  lm_n_par = 3
  ! lm_m_dat = n_max_read - n_min_read + 1

  ALLOCATE( fvec(lm_m_dat) )
  ALLOCATE( para(lm_n_par) )
  ALLOCATE( ipvt(lm_n_par) )

  normalize_sigma = parame(1,2)
  para(1) = 1.0
  para(2) = 0.0
  para(3) = 0.0
  !para(3) = parame(1,2)
  !para(1) = parame(1,2)
  !para(2) = parame(1,3)

  epsfcn = 1.0E-6**2     ! sqrt of relat. step size (finite differences)
  lm_factor = 1.0        ! maximum initial step size (parameter iteration)

  CALL lmdif1(match_bias_obj,lm_m_dat,lm_n_par,para,fvec,tol,epsfcn,lm_factor,info,ipvt)

  DEALLOCATE( fvec, para, ipvt )
  CLOSE( 67 )

  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' bias potential written to file: output_file/GCMC_bias.dat'
  WRITE(*,*)' ----------------------------------------------------------'
  WRITE(*,*)' '

END SUBROUTINE match_bias_potential



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE match_bias_potential
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!SUBROUTINE match_bias_obj ( n_max, box_l, fres_v, fres_l, rhob )
SUBROUTINE match_bias_obj (lm_m_dat, lm_n_par, para, fvec, iflag)

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  USE DFT_MODULE, ONLY: fres_temp
  USE match_fit_parameters
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! INTEGER, PARAMETER                     :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, INTENT(IN)                    :: lm_m_dat
  INTEGER, INTENT(IN)                    :: lm_n_par
  REAL, INTENT(IN)                       :: para(:)
  REAL, INTENT(IN OUT)                   :: fvec(:)
  INTEGER, INTENT(IN OUT)                :: iflag

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: rhopt

  INTEGER                                :: n_min, n_tmp
  REAL                                   :: temp1(0:8000), temp2(0:8000)
  REAL                                   :: lnphi(np,nc)
  REAL                                   :: p_n_sim(4000)
  !-----------------------------------------------------------------------------

  IF (iflag == 1) WRITE (*,'(a,2(G15.8))') ' parameter ',( para(i), i = 1, lm_n_par )

  !parame(1,2) = para(1)
  !parame(1,3) = para(2)

  fvec = 0.0

  read (67,*) n_min

  DO i = 1, lm_m_dat

     read (67,*) n_tmp, p_n_sim(i)
     rhopt = REAL(n_tmp) / box_l**3
     IF (rhopt == 0.0) rhopt = 1.d-6
     dense(1)  = rhopt * z3t
     densta(1) = rhopt * z3t
     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     CALL fugacity (lnphi)
     temp2(i) = fres_temp*rhopt + rhopt*( LOG(rhopt) - 1.0 )
     temp2(i) = temp2(i) * normalize_sigma**3

     temp1(i)= -(  fres_l + (fres_v-fres_l)/(rhob(2)-rhob(1))*(rhopt-rhob(1))  &
          - temp2(i) )/2.0 *box_l**2


     temp1(i) = para(1) * (temp1(i) - (REAL(n_tmp)*para(2)) ) + para(3)
     !temp1(i) = (temp1(i) - (REAL(n_tmp)*para(3)) )
     ! WRITE (68,'(i4,G18.10)') i, temp1(i)

     fvec(i) = 1.0 * ( temp1(i) - p_n_sim(i) ) !/ p_n_sim(i)
     !if (iflag==1) write (*,*) n_tmp,p_n_sim(i),temp1(i)
     !IF (iflag == 1 .AND. i==2) write (*,*) i, n_tmp,dense(1), temp1(i), p_n_sim(i)

  END DO

  IF (iflag == 1) WRITE (*,*) 'AAD (%)', SUM( ABS(fvec(1:lm_m_dat)) / REAL(lm_m_dat) * 100.0 ), &
       ' RMS (%)',SQRT( SUM( fvec(1:lm_m_dat)*fvec(1:lm_m_dat)  &
       / (REAL(lm_m_dat)) ) )* 100.0

  rewind( 67 )

  rewind( 68 )
  DO i = 0, 1180
     rhopt = REAL(i) / box_l**3
     IF (rhopt == 0.0) rhopt = 1.d-6
     dense(1)  = rhopt * z3t
     densta(1) = rhopt * z3t
     ensemble_flag = 'tv'                  ! is set for calculating: my_res=my_res(T,rho)
     CALL fugacity (lnphi)
     temp2(i) = fres_temp*rhopt + rhopt*( LOG(rhopt) - 1.0 )
     temp2(i) = temp2(i) * normalize_sigma**3

     temp1(i)= -(  fres_l + (fres_v-fres_l)/(rhob(2)-rhob(1))*(rhopt-rhob(1))  &
          - temp2(i) )/2.0 *box_l**2

     temp1(i) = para(1) * (temp1(i) - (REAL(i)*para(2)) ) + para(3)
     !temp1(i) = (temp1(i) - (REAL(i)*para(3)) )
     !IF (iflag == 1 .AND. i==950) write (*,*) i, i, dense(1), temp1(i)
     !IF (iflag == 1 .AND. i==950) call paus (' ')
     WRITE(68,'(i5,G18.10)') i, temp1(i)
  END DO

END SUBROUTINE match_bias_obj
