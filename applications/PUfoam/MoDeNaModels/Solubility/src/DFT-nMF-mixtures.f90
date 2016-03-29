!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DFT_nMFT2
!
! This subroutine performs a DFT calculation of the PCP-SAFT model using
! the fundamental measure theory (FMT) for the hard-sphere fluid, the
! TPT1-theory for chain formation and a first (or second) order perturba-
! tion theory for the dispersive interactions.
!
! Background to the code:
! The functional derivative of d(F)/d(rho) is done in two parts, as
! 'myrho(i) + mu_att(i)' , where 'mu_att' comprises the perturbation
! theory (dispersive attraction) as well as the non-local part of the hard-
! sphere term and chain term. For the hard-sphere term mu_att contains
! contributions of the form   mu_FMT = d(F^hs)/d(rho*) - (dF^hs/drho*)_LDA,
! where the index 'hs' is for hard-sphere (analogous for the chain term)
! and 'LDA' is for 'local density approximation', i.e. the values for the
! local density.
!
! The term 'myrho' contains the hard-sphere term and the chain term in
! the LDA. Because the PC-SAFT dispersion term is not identical to the
! Barker-Henderson (BH) perturtation theory (first or second order) for
! LJ fluids, the term 'myrho' also contains the difference of the PC-SAFT
! dispersion contribution and the BH perturbation theory used here the
! LDA. This difference is obtained by using the flag {subtract1 = '1PT'}
! the effect of the flag is that the BH-perturbation theory is subtracted
! from the PC-SAFT model for calculating 'myrho'.
!
! Caution: the indication for using the second order term in the pertur-
! bation expansion is currently soley given by the flag {subtract1 = '2PT'}
! The second order case needs to be re-tested before usage !
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE DFT_nMFT_mix

  USE parameters, ONLY: PI, RGAS, KBOL
  USE basic_variables
  USE EOS_VARIABLES, ONLY: fres, eta, dhs, mseg, uij, sig_ij, rho, x
  USE STARTING_VALUES
  USE DFT_MODULE
  use utilities
  use EOS_polar, only: f_polar, phi_polar
  USE EOS_NUMERICAL_DERIVATIVES

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, ik
  INTEGER                                :: l, m
  INTEGER                                :: kk, converg, read_file
  INTEGER                                :: discret, fa(nc), outpt, irc(nc), irc_j, maxiter, ih
  LOGICAL                                :: diagram
  REAL, DIMENSION(-NDFT:NDFT)            :: zp
  REAL, DIMENSION(-NDFT:NDFT,2)          :: rhop, rhop_o, r_corr
  REAL, DIMENSION(-NDFT:NDFT)            :: f_tot, f_corr
  REAL                                   :: f_L, f_R
  REAL                                   :: f_att
  REAL, DIMENSION(-NDFT:NDFT,2)          :: dF_drho_tot
  REAL, DIMENSION(2)                     :: dF_drho_att
  REAL, DIMENSION(2,0:nc)                :: rhob
  REAL, DIMENSION(nc)                    :: dhs_star
  REAL                                   :: sigij, dij
  REAL                                   :: f_disp_1PT, f_disp_pcsaft
  REAL, DIMENSION(2)                     :: mu_disp_1PT, mu_disp_pcsaft
  REAL                                   :: zges(np), pbulk
  REAL                                   :: delta_st, tanhfac
  REAL                                   :: my0(nc)
  REAL                                   :: f_int_r, mu_int1_r, mu_int2_r
  REAL                                   :: zz1
  REAL                                   :: dz_local
  REAL                                   :: dev, maxdev
  REAL                                   :: delta_f, free
  REAL                                   :: surftens(0:240), st_macro(240)
  REAL                                   :: f01, f02, f03
  ! REAL                                   :: z3,z2,z1,z0,zms,z0_rk,z1_rk,z2_rk,mhc(2)       ! temporary
  ! REAL                                   :: gij(2,2),gij_rk(2,2),dij_ab(2,2),tgij(2),lngij_rk(2) ! temporary
  ! REAL                                   :: mu_dsp(nc)                                           ! temporary
  INTEGER                                :: dum_i
  REAL                                   :: dum_r
  CHARACTER*100                          :: dum_text
  REAL                                   :: tc !, pc, rhoc, tsav, psav
  REAL                                   :: density(np), w(np,nc), lnphi(np,nc)
  REAL                                   :: damppic(2)
  REAL                                   :: box_l_no_unit
  REAL, DIMENSION(-NDFT:NDFT, 2)         :: lambda, rhobar
  REAL, DIMENSION(-NDFT:NDFT)            :: phi_dn0, phi_dn1, phi_dn2
  REAL, DIMENSION(-NDFT:NDFT)            :: phi_dn3, phi_dn4, phi_dn5
  REAL, DIMENSION(0:5,-NDFT:NDFT)        :: ni
  REAL                                   :: zs
  REAL                                   :: term1(nc), term2, I2_up(nc,nc), last_z = 0.0
  REAL                                   :: f_ch, dF_drho_ch(nc)
  REAL                                   :: f_fmt, dF_drho_fmt(nc)
  REAL                                   :: tc_L, tc_V, xi_save(np,nc), m_average
  REAL                                   :: f_assoc, mu_assoc(nc)
  CHARACTER (LEN=4)                      :: char_len
  ! REAL                                   :: steps, end_x
  REAL                                   :: fres_polar, fdd, fqq, fdq
  REAL                                   :: mu_polar(nc), fdd_rk, fqq_rk, fdq_rk, z3_rk
  !-----------------------------------------------------------------------------
  ! note: variable    ensemble_flag = 'tp'    is set for calculating: mu_res=mu_res(T,p)
  ! note: variable    ensemble_flag = 'tv'    is set for calculating: mu_res=mu_res(T,rho)
  ! note: variable    subtract1 = 'no'        is set for regular calculations
  ! note: variable    subtract1 = '1PT'       is set for calculating properties subtracting a
  !                                           first order perturbation term (dispersion)

  num = 1
  CALL SET_DEFAULT_EOS_NUMERICAL


  WCA = .false.
  shift = .false.
  diagram = .true.

  CALL READ_INPUT


  dhs(1:ncomp) = parame(1:ncomp,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1:ncomp,3)/t ) )  ! needed for rdf_matrix
  dhs_star(1:ncomp) = dhs(1:ncomp)/parame(1:ncomp,2)
  ! z3t_st = PI/6.0* parame(1,1) * d_hs**3

  IF ( ncomp > 2 ) THEN
     write (*,*) 'SPECIFY ONLY ONE OR TWO COMPONENTS IN THE INPUT-FILE:'
     write (*,*) '    ./input_file/INPUT.INP'
     stop
  END IF
  OPEN (68,FILE = './output_file/DFT_profiles.xlo')
  OPEN (69,FILE = './output_file/DFT_sigma.xlo')
  OPEN (71,FILE = './output_file/DFT_iteration.xlo')
  OPEN (72,FILE = './output_file/DFT_h.xlo')


  !-----------------------------------------------------------------------------
  ! The the cut-off distance rc is the max. distance for integrating inter-
  ! actions. Beyond rc, tail-corrections are added.
  !-----------------------------------------------------------------------------
  rc = 9.0                      ! dimensionless cut-off distance rc = r_c/sigma


  !-----------------------------------------------------------------------------
  ! basic definitions for calculating density profile,
  !-----------------------------------------------------------------------------
  ! grid-size  dzp = zp(1)-zp(0)

  box_l_no_unit = 250.0         ! lenth of simulation box (Angstrom)
  discret= 1000                  ! number of spacial discretizations for profile
  dzp    = box_l_no_unit / REAL(discret)    ! grid-distance (Angstrom)
  fa(1:ncomp) = NINT( parame(1:ncomp,2) / dzp + 1 )    ! number of steps per sigma (currently of component 1)
  outpt  = 100                  ! number output-points = discret/outpt


  !-----------------------------------------------------------------------------
  ! definitions for the numerical algorithm
  !-----------------------------------------------------------------------------

  maxiter = 200                 ! maximum number of iterations per temp.step
  maxdev  = 1.E-6               ! maximum deviation
  damppic(1) = 0.00001                ! damping of Picard-iteration
  damppic(2) = 0.00001               ! damping of Picard-iteration


  !-----------------------------------------------------------------------------
  ! get a matrix of values for the pair correlation function g(r,rho)
  !-----------------------------------------------------------------------------

  rg = 4.0
  den_step = 40
  !dzr = dzp / 2.0              ! dimensionless grid-distance for integrating
  !                             ! the Barker-Henderson attraction term
  dzr = 0.02                    ! dimensionless grid-distance for integrating
  ! the Barker-Henderson attraction term
  CALL rdf_matrix_mix


  !-----------------------------------------------------------------------------
  ! prepare for phase equilibrium calculation for given T
  !-----------------------------------------------------------------------------

  ! Diagram_for_various_T_Loop: DO

  nphas  = 2
  outp = 0                      ! output to terminal


  !-----------------------------------------------------------------------------
  ! calculate phase equilibrium for given T
  !-----------------------------------------------------------------------------

  CALL START_VAR (converg)      ! gets starting values, sets "val_init"

  IF ( converg /= 1 ) THEN
     WRITE (*,*) 'no VLE found'
     RETURN
  END IF
  ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!$    n_unkw = ncomp  ! number of quantities to be iterated
!!$    it(1) = 'p'   ! iteration of temperature
!!$    it(2) = 'x21'               ! iteration of mol fraction of comp.1 phase 2
!!$    sum_rel(1) = 'x12'          ! summation relation: x12 = 1 - sum(x1j)
!!$    sum_rel(2) = 'x22'          ! summation relation: x22 = 1 - sum(x2j)
!!$    end_x = 0.4
!!$    running = 'x11'               ! xi(1,1) is running var. in PHASE_EQUILIB
!!$    steps = 5.0
!!$    ensemble_flag = 'tp'                 ! this flag is for: 'regular calculation'
!!$    CALL PHASE_EQUILIB( end_x, steps, converg )
  ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  ! rhob(phase,0): molecular density
  rhob(1,0) = dense(1) / (  PI/6.0* SUM( xi(1,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
  rhob(2,0) = dense(2) / (  PI/6.0* SUM( xi(2,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
  ! rhob(phase,i): molecular component density (with i=(1,...ncomp) ) in units (1/A^3)
  rhob(1,1:ncomp) = rhob(1,0)*xi(1,1:ncomp)
  rhob(2,1:ncomp) = rhob(2,0)*xi(2,1:ncomp)
  WRITE (*,*) ' '
  WRITE (*,*) 'temperature  ',t,      'K,  and  p=', p/1.E5,' bar'
  WRITE (*,*) 'x1_liquid    ',xi(1,1),'   x1_vapor', xi(2,1)
  WRITE (*,*) 'densities    ',rhob(1,0), rhob(2,0)
  WRITE (*,*) 'dense        ',dense(1), dense(2)


  !-----------------------------------------------------------------------------
  ! get density in SI-units (kg/m**3)
  !-----------------------------------------------------------------------------
  CALL SI_DENS ( density, w )


  !-----------------------------------------------------------------------------
  ! (re-)calculate residual chemical potential of both phases
  !-----------------------------------------------------------------------------
  ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)
  densta(1) = dense(1)                    ! Index 1 is for liquid density (here: packing fraction eta)
  densta(2) = dense(2)                    ! Index 2 is for vapour density (here: packing fraction eta)
  CALL fugacity (lnphi)
  my0(1:ncomp) = lnphi(1,1:ncomp)! + LOG( rhob(1,1:ncomp) )     ! my0 = mu_res(T,rho_bulk_L) + ln(rho_bulk_l)
  !zges(1) = p / ( RGAS*t*density(1) ) * mm_average/1000.0
  !zges(2) = p / ( RGAS*t*density(2) ) * mm_average/1000.0
  zges(1) = ( p * 1.E-30 ) / ( KBOL*t* rhob(1,0) )
  zges(2) = ( p * 1.E-30 ) / ( KBOL*t* rhob(2,0) )

  pbulk = zges(1) * rhob(1,0)              ! pressure  p/kT (= Z*rho)
  WRITE (*,*) ' '
  WRITE (*,*) 'chem. potentials comp. 1', lnphi(1,1) + LOG( rhob(1,1) ),  &
       lnphi(2,1) + LOG( rhob(2,1) )
  WRITE (*,*) 'chem. potentials comp. 2', lnphi(1,2) + LOG( rhob(1,2) ),  &
       lnphi(2,2) + LOG( rhob(2,2) )
  WRITE (*,*) ' '


  !-----------------------------------------------------------------------------
  ! determine the critical temp. to xi of liquid and to xi of vapor
  !-----------------------------------------------------------------------------

  num = 0
  xi_save(:,:) = xi(:,:)
  dense(1) = 0.15
  WRITE (*,*) 'provide an estimate of the crit. Temp. of the mixture'
  READ (*,*) t

  xiF(1:ncomp) = xi_save(1,1:ncomp)
  CALL Heidemann_Khalil
  tc_L = t
  WRITE (*,*) 'critical temperature to xi_liquid',tc_L

  xiF(1:ncomp) = xi_save(2,1:ncomp)
  ! dense(1) = 0.15
  ! WRITE (*,*) 'provide an estimate of the crit. Temp. of the mixture'
  ! READ (*,*) t
  CALL Heidemann_Khalil
  tc_V = t
  WRITE (*,*) 'critical temperature to xi_vapor ',tc_V
  num = 1

  ! tc_L = 600.0
  ! WRITE (*,*) 'I have tentitativly set tc=700 '
  ! pause


  !tc = ( tc_L + tc_V ) / 2.0
  tc = tc_L
  xi(:,:) = xi_save(:,:)
  densta(1:nphas) = val_conv(1:nphas)
  dense(1:nphas) = val_conv(1:nphas)
  t = val_conv(3)
  WRITE (*,*) 'estimate of critical temperature:',tc,'K'
  WRITE (*,*) ' '

  !tc = (  xi(1,1)*630.39 + xi(1,2)*551.92  &
  !      + xi(2,1)*630.39 + xi(2,2)*551.92 ) / 2.0     - 10.0
  !tc = (  xi(1,1)*191.406 + xi(1,2)*133.68  &
  !      + xi(2,1)*191.406 + xi(2,2)*133.68 ) / 2.0


  !-----------------------------------------------------------------------------
  ! update z3t, the T-dependent quantity that relates eta and rho, as eta = z3t*rho
  !-----------------------------------------------------------------------------
  CALL PERTURBATION_PARAMETER


  !-----------------------------------------------------------------------------
  ! define initial density profile rhop(i)
  ! and dimensionless space coordinates zp(i)
  !
  ! discret  : number of grid-points within "the box"
  ! irc      : number of grid-points extending the the box to left and right
  !            the box is extended in order to allow for numerical integration
  !            'irc' is determined by the cut-off distance 'rc'
  !-----------------------------------------------------------------------------

  irc(1:ncomp) = NINT(rc*parame(1:ncomp,2)/dzp) + 1
  tanhfac = -2.3625*t/tc + 2.4728
  ! tanhfac = 2.7*(1.0 - t/tc)   ! this parameterization was published (Gross, JCP, 2009)

  irc_j = MAXVAL( irc(1:ncomp) )
  DO j = 1, ncomp
     DO i = -irc_j, (discret+irc_j)
        zp(i) = REAL(i) * dzp
     END DO
     DO i = -irc_j, (discret+irc_j)
        !rhop(i,j) = TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) * (rhob(1,j)-rhob(2,j))/2.0  &
        !        + (rhob(1,j)+rhob(2,j))/2.0
        rhop(i,j) =  ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) + 1.0 ) * rhob(1,j)/2.0 &
             - ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) - 1.0 ) * rhob(2,j)/2.0
        ! rhop(i,j) = rhob(1,j)
     END DO
  END DO


  !-----------------------------------------------------------------------------
  ! Optional: read starting profile of densities from file: DFT_profile.xlo
  !-----------------------------------------------------------------------------
  write (*,*) ' read starting profile of rho(z) from file: DFT_profile.xlo? (0: no, 1: yes)'
  read (*,*) read_file
  IF (read_file == 1 ) THEN
     READ (68,*) dum_text
     READ (68,*) dum_text
     READ (68,*) dum_text
     READ (68,*) dum_text
     DO i = 0, discret-1
        READ(68,*) dum_i, dum_r, rhop(i,1), rhop(i,2)
     END DO
     REWIND (68)
  END IF


  !-----------------------------------------------------------------------------
  ! Initialize the DENS_INV subroutine
  !-----------------------------------------------------------------------------

  nphas = 1
  ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)



  !=============================================================================
  ! Start iterating the density profile
  !=============================================================================

  kk = 1
  ih = 85

  DFT_Convergence_Loop: DO

     !--------------------------------------------------------------------------
     ! Getting auxilliary quantities along the profile
     !--------------------------------------------------------------------------

     CALL aux_chain_mix ( discret, fa, dzp, zp, rhop, rhobar, lambda )

     CALL fmt_dens_mix ( discret, fa, dzp, zp, rhop, ni,  &
          phi_dn0, phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5 )


     !--------------------------------------------------------------------------
     ! Start loop for density profile
     !--------------------------------------------------------------------------

     DO i = 0, discret

        !-----------------------------------------------------------------------
        ! hard sphere contribution
        !    f_fmt  is the Helmholtz energy density F_FMT/(VkT) = F_FMT/(NkT)*rho
        !    dF_drho_fmt is the functional derivative of F_FMT/(kT) to rho
        !-----------------------------------------------------------------------
        zs = 1.0 - ni(3,i)
        f_fmt  = - ni(0,i)*LOG(zs) + ni(1,i)*ni(2,i)/zs - ni(4,i)*ni(5,i)/zs  &
             + (ni(2,i)**3 -3.0*ni(2,i)*ni(5,i)*ni(5,i)) *(ni(3,i)+zs*zs*LOG(zs))  &
             /36.0/PI/zs/zs/ni(3,i)**2

        CALL dF_FMT_drho_mix ( i, fa, dzp, zp, phi_dn0,  &
             phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5, dF_drho_fmt )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$       z0 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) )
!!$       z1 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp) )
!!$       z2 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp)**2 )
!!$       z3 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
!!$       z0_rk = PI/6.0 * mseg(1)
!!$       z1_rk = PI/6.0 * mseg(1) * dhs(1)
!!$       z2_rk = PI/6.0 * mseg(1) * dhs(1)*dhs(1)
!!$       z3_rk = PI/6.0 * mseg(1) * dhs(1)**3
!!$       zms = 1.0 - z3
!!$       write(*,*)  f_fmt,SUM( rhop(i,1:ncomp)*parame(1:ncomp,1) )*(3.0*z1*z2/zms  &
!!$            + z2**3 /z3/zms/zms + (z2**3 /z3/z3-z0)*LOG(zms) )/z0
!!$       write(*,*)  'error A',dF_drho_fmt(1)- 6.0/PI* (  3.0*(z1_rk*z2+z1*z2_rk)/zms + 3.0*z1*z2*z3_rk/zms/zms  &
!!$                + 3.0*z2*z2*z2_rk/z3/zms/zms + z2**3 *z3_rk*(3.0*z3-1.0)/z3/z3/zms**3   &
!!$                + ((3.0*z2*z2*z2_rk*z3-2.0*z2**3 *z3_rk)/z3**3 -z0_rk) *LOG(zms)  &
!!$                + (z0-z2**3 /z3/z3)*z3_rk/zms  )
!!$       z0_rk = PI/6.0 * mseg(2)
!!$       z1_rk = PI/6.0 * mseg(2) * dhs(2)
!!$       z2_rk = PI/6.0 * mseg(2) * dhs(2)*dhs(2)
!!$       z3_rk = PI/6.0 * mseg(2) * dhs(2)**3
!!$       write(*,*)  'error A',dF_drho_fmt(2)- 6.0/PI* (  3.0*(z1_rk*z2+z1*z2_rk)/zms + 3.0*z1*z2*z3_rk/zms/zms  &
!!$                + 3.0*z2*z2*z2_rk/z3/zms/zms + z2**3 *z3_rk*(3.0*z3-1.0)/z3/z3/zms**3   &
!!$                + ((3.0*z2*z2*z2_rk*z3-2.0*z2**3 *z3_rk)/z3**3 -z0_rk) *LOG(zms)  &
!!$                + (z0-z2**3 /z3/z3)*z3_rk/zms  )
!!$       pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !-----------------------------------------------------------------------
        ! chain term
        !-----------------------------------------------------------------------
        CALL dF_chain_drho_mix ( i, fa, dzp, zp, rhop, lambda, rhobar, f_ch, dF_drho_ch )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$       z0 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) )
!!$       z1 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp) )
!!$       z2 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp)**2 )
!!$       z3 = SUM( rhop(i,1:ncomp) * PI/6.0 * parame(1:ncomp,1) * dhs(1:ncomp)**3 )
!!$       zms = 1.0 - z3
!!$       z1_rk = PI/6.0 * mseg(1) * dhs(1)
!!$       z2_rk = PI/6.0 * mseg(1) * dhs(1)*dhs(1)
!!$       z3_rk = PI/6.0 * mseg(1) * dhs(1)**3
!!$       DO l = 1, ncomp
!!$          DO j = 1, ncomp
!!$             dij_ab(l,j)=dhs(l)*dhs(j)/(dhs(l)+dhs(j))
!!$             gij(l,j) = 1.0/zms + 3.0*dij_ab(l,j)*z2/zms/zms + 2.0*(dij_ab(l,j)*z2)**2 /zms**3
!!$             gij_rk(l,j) = z3_rk/zms/zms  &
!!$                  + 3.0*dij_ab(l,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
!!$                  + dij_ab(l,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
!!$             if (l == j ) tgij(l) = gij(l,j)
!!$             if (l == j ) lngij_rk(l) = gij_rk(l,j) / gij(l,j)
!!$          END DO
!!$       END DO
!!$       write(*,*)  f_ch,SUM( rhop(i,1:ncomp) *(1.0- mseg(1:ncomp)) *LOG(tgij(1:ncomp)) )
!!$
!!$       mhc(1) = ( 1.0-mseg(1)) * LOG( gij(1,1) ) + SUM( rhop(i,1:ncomp) * (1.0-mseg(1:ncomp)) * lngij_rk(1:ncomp) )
!!$       write(*,*)  'error 1', dF_drho_ch(1) - mhc(1)
!!$       z1_rk = PI/6.0 * mseg(2) * dhs(2)
!!$       z2_rk = PI/6.0 * mseg(2) * dhs(2)*dhs(2)
!!$       z3_rk = PI/6.0 * mseg(2) * dhs(2)**3
!!$       DO l = 1, ncomp
!!$          DO j = 1, ncomp
!!$             gij_rk(l,j) = z3_rk/zms/zms  &
!!$                  + 3.0*dij_ab(l,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
!!$                  + dij_ab(l,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
!!$             if (l == j ) lngij_rk(l) = gij_rk(l,j) / gij(l,j)
!!$          END DO
!!$       END DO
!!$       mhc(2) = ( 1.0-mseg(2)) * LOG( gij(2,2) ) + SUM( rhop(i,1:ncomp) * (1.0-mseg(1:ncomp)) * lngij_rk(1:ncomp) )
!!$       write(*,*)  'error 2',dF_drho_ch(2) - mhc(2)
!!$       pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !-----------------------------------------------------------------------
        ! dispersive attraction
        !-----------------------------------------------------------------------
        f_att   = 0.0
        dF_drho_att(:) = 0.0
        term2 = 0.0
        I2_up(:,:) = 0.0
        DO l = 1, ncomp
           term1(l) = 0.0
           DO m = 1, ncomp
              f01 = 0.0
              f02 = 0.0
              f03 = 0.0
              term2 = 0.0
              sigij = sig_ij(l,m)
              dij = ( dhs(l) + dhs(m) ) / 2.0
              irc_j = ( irc(l) + irc(m) + 1 ) / 2
              DO j = (i-irc_j), (i+irc_j)             ! integration in z-coordinate
                 zz1 = ABS( zp(j) - zp(i) )           ! distance z12 between 1 and 2
                 dz_local = dzp
                 IF ( zz1 < rc*sigij ) THEN
                    IF ( -( zp(j) - zp(i) ) >= rc*sigij - dzp ) dz_local = rc*sigij - zz1
                    CALL dft_rad_int_mix ( i, j, l, m, ih, zz1, rhop, sigij, dij, f_int_r, mu_int1_r, mu_int2_r )

                    f_att   = f_att    + PI*mseg(l)*mseg(m) *uij(l,m)/t * rhop(i,l) *  &
                         dz_local * ( rhop(j,m)* f_int_r + f01 ) / 2.0
                    term1(l)= term1(l) + 2.0*PI*mseg(l)*mseg(m) *uij(l,m)/t *  &
                         dz_local * ( rhop(j,m)*mu_int1_r + f02 ) / 2.0
                    !                   term2   = term2    + PI*mseg(l)*mseg(m) *uij(l,m)/t * rhop(i,l) *  &
                    !                                        dz_local * ( rhop(j,m)* mu_int2_r + f03 ) / 2.0
                    I2_up(l,m) = I2_up(l,m) + dz_local * ( rhop(j,m)* mu_int2_r + f03 ) / 2.0
                    f01 = rhop(j,m)* f_int_r
                    f02 = rhop(j,m)* mu_int1_r
                    f03 = rhop(j,m)* mu_int2_r
                    last_z = zz1
                 ELSE
                    f01 = rc * sig_ij(l,m)  *  4.0 * ( rc**-12 - rc**-6 ) * rhop(j,m)
                    f02 = rc * sig_ij(l,m)  *  4.0 * ( rc**-12 - rc**-6 ) * rhop(j,m)
                    f03 = 0.0
                 END IF

              END DO
              dz_local = rc*sigij - last_z
              f_att  = f_att + PI*mseg(l)*mseg(m) *uij(l,m)/t * rhop(i,l) * dz_local * f01
              term1(l) = term1(l) + 2.0*PI*mseg(l)*mseg(m) *uij(l,m)/t *  dz_local * f02

           END DO
        END DO

        term2 = 0.0
        DO l = 1, ncomp
           DO m = 1, ncomp
              term2 = term2 + PI * rhop(i,l) * mseg(l)*mseg(m) * uij(l,m)/t *I2_up(l,m)
           END DO
        END DO


        DO l = 1, ncomp
           dF_drho_att(l) = term1(l) + term2 * PI / 6.0 * mseg(l) * dhs(l)**3
        END DO



        !-----------------------------------------------------------------------
        ! cut-off corrections
        !-----------------------------------------------------------------------

        ! CALL cutoff_corrections_mix ( i, irc, rhop, rhob, f_att, dF_drho_att )



        !-----------------------------------------------------------------------
        ! The integration of the DFT-integral can be done with cubic splines
        !-----------------------------------------------------------------------
        ! stepno = discret + irc*2
        ! CALL SPLINE_PARA (dzp, intgrid, utri, stepno)
        ! CALL SPLINE_INT (f_int_z, dzp, intgrid, utri, stepno)
        ! f_int_z = f_int_z + rhob(1)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )
        ! f_int_z = f_int_z + rhob(2)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )


        !-----------------------------------------------------------------------
        ! Calculation of 'dF_drho_att', which denotes non-local contributions
        ! of d(F)/d(rho*). Apart from the dispersive attraction, it contains
        ! non-local contributions of the hard-sphere term, and the chain term.
        ! 'dF_drho_att' is defined as:      dF_drho_att = d(F)/d(rho*) - (dF/drho*)_LDA
        ! where LDA is for 'local density approximation', i.e. the values for
        ! the local density.
        !
        ! F is the intrinsic Helmholtz energy
        ! rho* = rhop(i) denotes the dimensionless density.
        ! parame(1,3) is epsilon/k  (LJ-energy parameter)
        ! parame(1,2) is sigma      (LJ-size parameter)
        !
        ! For a non-local second order term uncomment lines starting with '!2'
        !-----------------------------------------------------------------------
        ! remember: factor 2*PI because of the cylindrical coordinates

        !       dF_drho_att(1:ncomp) = 2.0*PI*parame(1,1)**2 * parame(1,3)/t * mu_int_z(1:ncomp)



        !-----------------------------------------------------------------------
        ! make the dispersive attraction term consistent with PC-SAFT, by
        ! adding the difference ( PC-SAFT - 1PT )_dispersion locally (LDA)
        !-----------------------------------------------------------------------

        !*****************************************************************************
        ! under construction
        !*****************************************************************************
        dense(1) = PI / 6.0 * SUM( rhop(i,1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
        densta(1) = dense(1)
        xi(1,1:ncomp) = rhop(i,1:ncomp) / SUM( rhop(i,1:ncomp) )

        ensemble_flag = 'tv'

        call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PT_MIX   ' )
        call FUGACITY ( lnphi )
        call RESTORE_PREVIOUS_EOS_NUMERICAL
        f_disp_1PT  = fres * SUM( rhop(i,1:ncomp) )
        mu_disp_1PT(1:ncomp) = lnphi(1,1:ncomp)

        call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PC-SAFT  ' )
        call FUGACITY ( lnphi )
        call RESTORE_PREVIOUS_EOS_NUMERICAL
        f_disp_pcsaft  = fres * SUM( rhop(i,1:ncomp) )
        mu_disp_pcsaft(1:ncomp) = lnphi(1,1:ncomp)

        fres_polar  = 0.0
        mu_polar(:) = 0.0
        IF ( SUM( parame(1:ncomp,6) ) > 1.E-10 .OR. SUM( parame(1:ncomp,7) ) > 1.E-10 ) THEN
           eta = densta(1)
           rho = SUM(rhop(i,1:ncomp))
           x(1:ncomp) = xi(1,1:ncomp)
           call F_POLAR ( fdd, fqq, fdq )
           fres_polar  = ( fdd + fqq + fdq ) * SUM( rhop(i,1:ncomp) )
           DO ik = 1, ncomp
              z3_rk = PI/6.0 * mseg(ik) * dhs(ik)**3
              call PHI_POLAR ( ik, z3_rk, fdd_rk, fqq_rk, fdq_rk )
              mu_polar(ik) = fdd_rk + fqq_rk + fdq_rk
           END DO
        END IF

        f_assoc  = 0.0
        mu_assoc(:) = 0.0
        IF ( SUM( NINT(parame(1:ncomp,12)) ) > 0) THEN
           call ONLY_ONE_TERM_EOS_NUMERICAL ( 'hb_term  ', 'TPT1_Chap' )
           call FUGACITY ( lnphi )
           call RESTORE_PREVIOUS_EOS_NUMERICAL
           f_assoc  = fres * SUM( rhop(i,1:ncomp) )
           mu_assoc(1:ncomp) = lnphi(1,1:ncomp)
        END IF

        ensemble_flag = 'tp'

        !       CALL mu_pert_theory_mix ( mu_dsp )
        !       write (*,*) dF_drho_att(1) - mu_disp_1PT(1), dF_drho_att(1), mu_disp_1PT(1)
        !       write (*,*) dF_drho_att(2) - mu_disp_1PT(2), dF_drho_att(2),  mu_disp_1PT(2)
        !       write (*,*) f_att - f_disp_1PT, f_att, f_disp_1PT
        !       write (*,*) 'err dis a',dF_drho_att(1) - mu_dsp(1)
        !       write (*,*) 'err dis b',dF_drho_att(2) - mu_dsp(2)
        !       write (*,*) 'err dis 1',dF_drho_att(1) - mu_disp_1PT(1)
        !       write (*,*) 'err dis 2',dF_drho_att(2) - mu_disp_1PT(2)
        !       pause

        dF_drho_att(1:ncomp)= dF_drho_att(1:ncomp) + ( mu_disp_pcsaft(1:ncomp) - mu_disp_1PT(1:ncomp) )  &
             + mu_assoc(1:ncomp) + mu_polar(1:ncomp)
        f_att = f_att + ( f_disp_pcsaft - f_disp_1PT )  + f_assoc + fres_polar

        !*****************************************************************************
        !*****************************************************************************


        !-----------------------------------------------------------------------
        ! collect the total Helmholtz energy density 'f_tot' and the functional
        ! derivative to rhop(z) 'dF_drho_tot' (both are residual quantities)
        !-----------------------------------------------------------------------

        dF_drho_tot(i,1:ncomp) = dF_drho_fmt(1:ncomp) + dF_drho_ch(1:ncomp) + dF_drho_att(1:ncomp)
        f_tot(i) = f_fmt + f_ch + f_att


        !-----------------------------------------------------------------------
        ! on-the-fly report on local errors in the profile
        !-----------------------------------------------------------------------
        IF ( MOD(i, outpt) == 0 )write(*,'(i5,4(G15.6))') i, rhop(i,1), rhop(i,2),  &
             my0(1)-dF_drho_tot(i,1)+LOG( rhob(1,1)/rhop(i,1) ),  &
             my0(2)-dF_drho_tot(i,2)+LOG( rhob(1,2)/rhop(i,2) )


     END DO          ! end of loop (i = 0, discret) along the profile


     !--------------------------------------------------------------------------
     ! correction for numerical inaccuracies in calculating the density profile
     !--------------------------------------------------------------------------

     DO j = 1, ncomp
        f_L = EXP( my0(j) - dF_drho_tot(10,j) +LOG( rhob(1,j)/rhop(10,j) ) )
        f_R = EXP( my0(j) - dF_drho_tot(discret-10,j) +LOG( rhob(1,j)/rhop(discret-10,j) ))
        DO i = 0, discret
           r_corr(i,j) =  ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) + 1.0 ) * f_L/2.0 &
                - ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) - 1.0 ) * f_R/2.0
        END DO
     END DO

     !--------------------------------------------------------------------------
     ! update the density profile using either a Picard-iteration scheme
     ! (i.e. a direct substitution scheme, or a LDA - Inversion Procedure
     ! similar to R.Evans (Bristol), which is done in function DENS_INVERS2.
     !--------------------------------------------------------------------------
     dev = 0.0
     DO j = 1, ncomp
        DO i = 0, discret
           rhop_o(i,j) = rhop(i,j)
           rhop(i,j)   = rhob(1,j) * EXP( my0(j) - dF_drho_tot(i,j) ) /r_corr(i,j)
           rhop(i,j)   = rhop_o(i,j) + ( rhop(i,j) - rhop_o(i,j) )*damppic(j)
           dev = dev + ( ( rhop(i,j) - rhop_o(i,j) )*parame(j,2)**3 / damppic(j) )**2
        END DO
     END DO
     ! if ( kk == 4 )  damppic(1) = damppic(1) * 2.0
     if ( kk == 8 )  damppic(:) = damppic(:) * 10.0
     if ( kk == 16 )  damppic(:) = damppic(:) * 5.0
     dev = dev / real(discret)


     !--------------------------------------------------------------------------
     ! correction for numerical inaccuracies in calculating the surface tension
     !--------------------------------------------------------------------------

     f_L = f_tot(10) + SUM( rhop(10,1:ncomp)*( LOG(rhop(10,1:ncomp)/rhob(1,1:ncomp))-1.0 ) )  &
          - ( SUM(rhop(10,1:ncomp)*my0(1:ncomp)) - pbulk)
     f_R = f_tot(discret-10) + SUM( rhop(discret-10,1:ncomp)*( LOG(rhop(discret-10,1:ncomp)/rhob(1,1:ncomp))-1.0 ) )  &
          - ( SUM(rhop(discret-10,1:ncomp)*my0(1:ncomp)) - pbulk)
     j=1
     DO i = 1, discret
        f_corr(i) =  ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) + 1.0 ) * f_L/2.0 &
             - ( TANH(-(zp(i)-zp(INT(discret/2))) / parame(j,2) *tanhfac) - 1.0 ) * f_R/2.0
     END DO

     !--------------------------------------------------------------------------
     ! calculate surface tension
     !--------------------------------------------------------------------------
     free = 0.0
     DO i = 1, discret

        !-----------------------------------------------------------------------
        ! now add the ideal gas term
        !-----------------------------------------------------------------------
        f_tot(i) = f_tot(i) + SUM( rhop(i,1:ncomp)*( LOG(rhop(i,1:ncomp)/rhob(1,1:ncomp))-1.0 ) )

        delta_f = f_tot(i)  -   ( SUM(rhop(i,1:ncomp)*my0(1:ncomp)) - pbulk)  - f_corr(i) ! all quantities .../(kT)
        free  = free  + delta_f*dzp

     END DO
     surftens(kk) = KBOL * t *1.E20*1000.0 *free


     !--------------------------------------------------------------------------
     ! add an approximate capillary wave contribution to the surface tension
     !--------------------------------------------------------------------------
     m_average = SUM( rhob(1,1:ncomp)*parame(1:ncomp,1) ) / rhob(1,0)
     !m_average = 0.5 * (  SUM( rhob(1,1:ncomp)*parame(1:ncomp,1) ) / rhob(1,0)  &
     !                   + SUM( rhob(2,1:ncomp)*parame(1:ncomp,1) ) / rhob(2,0) )
     st_macro(kk) = surftens(kk) / ( 1.0 + 3.0/8.0/PI *t/tc  &
          * (1.0/2.55)**2  / (0.0674*m_average+0.0045) )


     delta_st = 1.0
     !IF ( kk > 1 ) delta_st = ABS( surftens(kk)-surftens(kk-1) ) / surftens(kk)

     WRITE (*,*)  '-----------------------------------------------------------'
     WRITE (*,*)  ' #      error-sum      intrinsic ST    total ST'
     WRITE (*,'(i3,3F15.6)')  kk, dev, surftens(kk), st_macro(kk)
     WRITE (*,*)  '-----------------------------------------------------------'
     WRITE (71,'(i3,3G18.8)') kk, dev, surftens(kk), st_macro(kk)
     kk = kk + 1

     !--------------------------------------------------------------------------
     ! convergence criterion
     !--------------------------------------------------------------------------
     ! write (*,*) 'error measure',dev, delta_st
     IF ( dev < maxdev .OR. delta_st < 1.E-8 ) THEN
        EXIT DFT_Convergence_Loop
     ELSE
        IF ( kk > maxiter ) THEN
           WRITE(*,*)' no convergence in ',maxiter,' steps'
           EXIT DFT_Convergence_Loop
        ENDIF
     ENDIF

  ENDDO DFT_Convergence_Loop

  !-----------------------------------------------------------------------------
  ! write resulting density profile
  !-----------------------------------------------------------------------------
  WRITE(68,'(a,4(G18.8))') 't,p,x11,x21',t, val_conv(4)/1.E5, rhob(1,1)/SUM(rhob(1,1:ncomp)),  &
       rhob(2,1)/SUM(rhob(2,1:ncomp))
  WRITE(68,'(a,2G18.8)') 'rho1,rho2',density(1), density(2)
  WRITE(68,'(a,2G18.8)') 'gamma',surftens(kk-1), st_macro(kk-1)
  WRITE(68,*) ' '
  DO i = 0, discret
     WRITE (char_len,'(I3)') 2*ncomp+1
     WRITE (68,'(i6,'//char_len//'(G18.10))') i,zp(i)-zp(INT(discret/2)),  &
          rhop(i,1:ncomp), -( dF_drho_tot(i,1:ncomp)-my0(1:ncomp) )
     ! IF ( MOD(i, outpt) == 0 ) WRITE(*,'(i6,3(f18.12))') i,zp(i)-zp(INT(discret/2)),rhop(i)
     ! write (69,*) zp(i)*parame(1,2), pN_pT(i)/parame(1,2)**3 /1E-30*parame(1,3)*KBOL/1.E6 ! in MPa
  END DO
  WRITE (68,*) ' '


!!$ ! ---------------------------------------------------------------------
!!$ ! summarize results
!!$ ! ---------------------------------------------------------------------
!!$ WRITE (*,*) ' '
!!$ WRITE (*,*) 'SUMMARY FOR A SINGLE TEMPERATURE'
!!$ WRITE (*,*) ' '
!!$ WRITE (*,*) 'Temp. [K], Pressure [bar]    ',t,val_conv(4)/1.E5
!!$ WRITE (*,*) 'Critical point Temp., Press. ',tc,pc/1.E5
!!$ WRITE (*,*) 'Density [kg/m**3]            ',density(1),density(2)
!!$ WRITE (*,*) 'Dimensionless Density (rho*) ',rhob(1),rhob(2)
!!$ WRITE (*,*) 'Excess Grand Potential       ',free
!!$ WRITE (*,*) 'Intrinsic Interf. Tension [mN/m] ',surftens(kk-1)
!!$ WRITE (*,*) 'Macroscop.Interf. Tension [mN/m] ',st_macro(kk-1)
!!$ WRITE (*,*) '============================================================'
!!$ WRITE (*,*) ' '
!!$ WRITE (69,'(9(f18.10))') t, val_conv(4)/1.E5,  &
!!$                          rhob(1),rhob(2),surftens(kk-1),st_macro(kk-1),free,dev
!!$
!!$ ! ---------------------------------------------------------------------
!!$ ! when calc. a phase diagram & diagram of surface tens., loop for higher T
!!$ ! ---------------------------------------------------------------------
!!$ ensemble_flag = 'tp'          ! this flag is for 'regular calculations'
!!$ subtract1     = 'no'          ! this flag is for 'regular calculations'
!!$ IF ( diagram ) THEN
!!$    IF ( (t+8.0) <= tc ) THEN
!!$       t = t + 5.0
!!$       IF ( (t+15.0) <= tc ) t = t + 5.0
!!$       IF ( (t+25.0) <= tc ) t = t + 10.0
!!$       IF ( (t+45.0) <= tc ) t = t + 20.0
!!$       nphas = 2
!!$       n_unkw = ncomp            ! number of quantities to be iterated
!!$       it(1) = 'p'               ! iteration of pressure
!!$       val_init(3) = t           ! value of temperature
!!$       running = 't'             ! T is running variable in PHASE_EQUILIB - here T=const.
!!$       end_x  = t                ! end temperature - here T=const.
!!$       steps = 1.0               ! number of steps of running var. - here 1, since T=const.
!!$       d_hs = parame(1,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1,3)/t ) )
!!$       z3t_st = PI/6.0* parame(1,1) * d_hs**3
!!$       dhs_st = d_hs/parame(1,2)
!!$       CALL rdf_matrix_units
!!$    ELSE
!!$       EXIT Diagram_for_various_T_Loop
!!$    END IF
!!$    WRITE (69,'(7(f18.10))') tc, pc/1.E5, rhoc, rhoc, 0., 0., 0.
!!$ ELSE
!!$    EXIT Diagram_for_various_T_Loop
!!$ END IF
!!$
!!$ ENDDO Diagram_for_various_T_Loop

END SUBROUTINE DFT_nMFT_mix




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE aux_chain_mix ( discret, fa, dzp, zp, rhop, rhobar, lambda )

  USE basic_variables
  USE EOS_VARIABLES, ONLY: dhs
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: discret
  INTEGER, INTENT(IN)                    :: fa(nc)
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT,2)
  ! REAL, INTENT(IN)                       :: rhob(2,0:nc)
  REAL, INTENT(OUT)                      :: rhobar(-NDFT:NDFT,2)
  REAL, INTENT(OUT)                      :: lambda(-NDFT:NDFT,2)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, nn
  REAL                                   :: int1, int2, zz1, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, rx, ry1, ry2
  !-----------------------------------------------------------------------------

  ! DO i=(-fa),(discret+fa)
  ! lambda(i) = 0.0
  ! rhobar(i) = 0.0
  ! f01 = 0.0  !  = 0.75*rhop(i-fa)*( 1.**2 -1.**2 )
  ! f02 = 0.5 *rhop(i-fa)
  ! DO j=(i-fa+1),(i+fa)      ! fa=NINT(1./dzp)
  !   zz1 = ABS(zp(j)-zp(i))  ! distance z12 between 1 and 2
  !   int_r= 0.75*rhop(j)*(1.0-zz1*zz1)
  !   int_l= 0.5 *rhop(j)
  !   rhobar(i) = rhobar(i) + dzp * ( int_r+f01 )/2.0
  !   lambda(i) = lambda(i) + dzp * ( int_l+f02 )/2.0
  !   f01=int_r
  !   f02=int_l
  ! ENDDO
  ! write(69,*) i,rhobar(i)
  ! ENDDO


  DO k = 1, ncomp
     ry1(:) = 0.0
     ry2(:) = 0.0
     rx(:)  = 0.0
     y2(:)  = 0.0
     DO i = (-fa(k)), (discret+fa(k))
        nn = 1
        DO j = (i-fa(k)), (i+fa(k))      ! fa(k) = NINT(sigma(k)/dzp)
           IF ( zp(j+1) > (zp(i)-dhs(k)) .AND. (zp(i)-dhs(k)) >= zp(j) ) THEN
              dz =  zp(j+1) - (zp(i)-dhs(k))
              rx(1)  = 0.0
              ry1(1) = 0.0  !  = 0.75*rhop(i-fa,k)*( d.**2 -d.**2 )
              ry2(1) = 0.5/dhs(k)*rhop(j,k) + ((zp(i)-dhs(k))-zp(j))/dzp  &
                   *( 0.5/dhs(k)*rhop(j+1,k) - 0.5/dhs(k)*rhop(j,k) )
           ELSE IF ( zp(j) > (zp(i)-dhs(k)) .AND. zp(j) <= (zp(i)+dhs(k)) ) THEN
              zz1 = ABS( zp(j)-zp(i) )          ! distance z12 between 1 and 2
              nn = nn + 1
              ry1(nn) = 0.75/dhs(k)**3 * rhop(j,k) * (dhs(k)**2 -zz1*zz1)
              ry2(nn) = 0.5/dhs(k) * rhop(j,k)
              rx(nn)  = rx(nn-1) + dz
              dz = dzp
              IF ( (zp(i)+dhs(k)) < zp(j+1) ) THEN
                 dz = (zp(i)+dhs(k)) - zp(j)
                 nn = nn + 1
                 rx(nn) = rx(nn-1) + dz
                 ry1(nn) = 0.0
                 ry2(nn) = 0.5/dhs(k)*rhop(j,k) + ((zp(i)+dhs(k))-zp(j))/dzp  &
                      *( 0.5/dhs(k)*rhop(j+1,k) - 0.5/dhs(k)*rhop(j,k) )
              END IF
           END IF
        END DO
        xl = rx(1)
        xh = rx(nn)
        CALL spline          ( rx(1:nn), ry1(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
        CALL splint_integral ( rx(1:nn), ry1(1:nn), y2(1:nn), nn, xl, xh, int1 )
        rhobar(i,k) = int1
        if ( rhobar(i,k) < 0.0 ) then
           rhobar(i,k) = 0.0
           do j = 2, nn
              rhobar(i,k) = rhobar(i,k) + (ry1(j)+ry1(j-1))/2.0 *(rx(j)-rx(j-1))
           end do
        end if
        if ( rhobar(i,k) < 0.0 ) rhobar(i,k) = rhop(i,k)
        CALL spline          ( rx(1:nn), ry2(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
        CALL splint_integral ( rx(1:nn), ry2(1:nn), y2(1:nn), nn, xl, xh, int2 )
        lambda(i,k) = int2
        if ( lambda(i,k) < 0.0 ) then
           lambda(i,k) = 0.0
           do j = 2, nn
              lambda(i,k) = lambda(i,k) + (ry2(j)+ry2(j-1))/2.0 *(rx(j)-rx(j-1))
           end do
        end if
        if ( lambda(i,k) < 0.0 ) lambda(i,k) = rhop(i,k)

        ! write (*,'(i5,i5,3G19.10)') k,i,rhob(1,k),rhob(2,k),lambda(i,k)

        ! if ( k == 1 .AND. lambda(i,k) < 0.5*rhob(2,k) ) write (*,*) 'warning: lambda too low',i,lambda(i,k),rhob(2,k)
        ! if ( k == 1 .AND. lambda(i,k) < 0.5*rhob(2,k) ) lambda(i,k) = 0.5*rhob(2,k)
     END DO
  END DO

END SUBROUTINE aux_chain_mix



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE fmt_dens_mix ( discret, fa, dzp, zp, rhop, ni,  &
     phi_dn0, phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5 )

  USE parameters, ONLY: PI
  USE basic_variables
  USE eos_variables, ONLY: dhs
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: discret
  INTEGER, INTENT(IN)                    :: fa(nc)
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT,2)
  REAL, INTENT(OUT)                      :: ni(0:5,-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn0(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn1(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn2(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn3(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn4(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn5(-NDFT:NDFT)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k, nn, fa2
  REAL                                   :: int2, int3, int5
  REAL                                   :: zz1, zs, d2, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy2, hy3, hy5
  !-----------------------------------------------------------------------------

  ni(2,:) = 0.0     ! corresponds to z2*6
  ni(3,:) = 0.0     ! corresponds to z3
  ni(5,:) = 0.0
  ni(0,:) = 0.0     ! corresponds to z0/PI*6.0
  ni(1,:) = 0.0     ! corresponds to z1/PI*3.0
  ni(4,:) = 0.0

  DO k = 1, ncomp
     hx(:)  = 0.0
     hy2(:) = 0.0
     hy3(:) = 0.0            ! = rhop(i-fa/2)*parame(1,1)*(0.25-(zp(i-fa/2)-zp(i))**2 )
     hy5(:) = 0.0
     DO i = (-fa(k)), (discret+fa(k))
        nn = 1
        d2 = dhs(k) / 2.0
        fa2= ( fa(k) + 1 ) / 2
        ! if (i==-fa(k)) write (*,*) 'even?',k,fa(k),fa2
        ! if (i==-fa(k)) pause
        DO j = (i-fa2), (i+fa2)
           IF ( zp(j+1) > (zp(i)-d2) .AND. (zp(i)-d2) >= zp(j) ) THEN
              dz =  zp(j+1) - (zp(i)-d2)
              hx(1) = 0.0
              hy2(1) = rhop(j,k)  *parame(k,1) + ((zp(i)-d2)-zp(j))/dzp  &
                   * ( rhop(j+1,k)*parame(k,1) - rhop(j,k)*parame(k,1) )
              hy3(1) = 0.0            ! = rhop(i-fa2)*parame(k,1)*(0.25-(zp(i-fa2)-zp(i))**2 )
              hy5(1) = rhop(j,k)*parame(k,1)*(zp(j)-zp(i)) +((zp(i)-d2)-zp(j))/dzp  &
                   * ( rhop(j+1,k)*parame(k,1)*(zp(j+1)-zp(i)) - rhop(j,k)  *parame(k,1)*(zp(j)-zp(i)) )
           ELSE IF ( zp(j) > (zp(i)-d2) .AND. zp(j) <= (zp(i)+d2) ) THEN
              zz1 = zp(j)-zp(i)       ! distance z12 between 1 and 2
              nn = nn + 1
              hy2(nn) = rhop(j,k)*parame(k,1)
              hy3(nn) = rhop(j,k)*parame(k,1) * ( d2**2 - zz1**2 )
              hy5(nn) = rhop(j,k)*parame(k,1) * zz1
              hx(nn)  = hx(nn-1) + dz
              dz = dzp
              IF ( (zp(i)+d2) < zp(j+1) ) THEN
                 dz = (zp(i)+d2) - zp(j)
                 nn = nn + 1
                 hy2(nn) = rhop(j,k)*parame(k,1) + ((zp(i)+d2)-zp(j))/dzp  &
                      * ( rhop(j+1,k)*parame(k,1) - rhop(j,k)*parame(k,1) )
                 hy3(nn) = 0.0
                 hy5(nn) = rhop(j,k)*parame(k,1)*(zp(j)-zp(i)) +((zp(i)+d2)-zp(j))/dzp  &
                      * ( rhop(j+1,k)*parame(k,1)*(zp(j+1)-zp(i)) - rhop(j,k)*parame(k,1)*(zp(j) -zp(i)) )
                 hx(nn) = hx(nn-1) + dz
              END IF
           END IF
        END DO
        xl = hx(1)
        xh = hx(nn)
        CALL spline          ( hx(1:nn), hy2(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
        CALL splint_integral ( hx(1:nn), hy2(1:nn), y2(1:nn), nn, xl, xh, int2 )
        CALL spline          ( hx(1:nn), hy3(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
        CALL splint_integral ( hx(1:nn), hy3(1:nn), y2(1:nn), nn, xl, xh, int3 )
        CALL spline          ( hx(1:nn), hy5(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
        CALL splint_integral ( hx(1:nn), hy5(1:nn), y2(1:nn), nn, xl, xh, int5 )
        ni(2,i) = ni(2,i) + PI * dhs(k) *int2                         ! corresponds to z2*6
        ni(3,i) = ni(3,i) + PI * int3                                 ! corresponds to z3
        ni(5,i) = ni(5,i) + 2.0* PI *int5
        ni(0,i) = ni(0,i) + (PI * dhs(k) *int2) / (PI*dhs(k)**2 )     ! corresponds to z0/PI*6.0
        ni(1,i) = ni(1,i) + (PI * dhs(k) *int2) / (2.0*PI*dhs(k))     ! corresponds to z1/PI*3.0
        ni(4,i) = ni(4,i) + (2.0* PI *int5) / (2.0*PI*dhs(k))
     END DO
  END DO

  !     derivatives of phi to ni
  DO i = (-MAXVAL(fa(1:ncomp))), (discret+MAXVAL(fa(1:ncomp)))
     zs = 1.0 - ni(3,i)
     phi_dn0(i) = - LOG(zs)   !!!!
     phi_dn1(i) = ni(2,i)/zs
     phi_dn2(i) = ni(1,i)/zs + 3.0*(ni(2,i)*ni(2,i) - ni(5,i)*ni(5,i))  &
          * (ni(3,i)+zs*zs*LOG(zs))/(36.0*PI*ni(3,i)*ni(3,i)*zs*zs)
     phi_dn3(i) = ni(0,i)/zs + (ni(1,i)*ni(2,i)-ni(4,i)*ni(5,i))/zs/zs  &
          - (ni(2,i)**3 - 3.0*ni(2,i)*ni(5,i)*ni(5,i))  &
          *( ni(3,i)*(ni(3,i)*ni(3,i)-5.0*ni(3,i)+2.0)  &
          + 2.0*zs**3 *LOG(zs) ) / ( 36.0*PI*(ni(3,i)*zs)**3 )
     phi_dn4(i) = - ni(5,i)/zs
     phi_dn5(i) = - ni(4,i)/zs - 6.0*ni(2,i)*ni(5,i)  &
          * (ni(3,i)+zs*zs*LOG(zs))/(36.0*PI*ni(3,i)*ni(3,i)*zs*zs)
  END DO

END SUBROUTINE fmt_dens_mix





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DFT_RAD_INT_mix
!
! This subroutine integrates the kernel of the perturbation theory in
! radial direction (more precisely in distance coordinate r^ hat). A
! non-mean field approach is taken, using a radial distribution function
! (currently at a Percus-Yevick level).
!
! The first and second order contributions to the perturbation theory
! are calculated - although the second order contribution is rarely used.
!
! ToDo: comment the subroutine and remove the GOTO constructs
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE DFT_RAD_INT_mix ( i, j, l, m, ih, zz_1, rhop, sigij, dij, f_int_r, my_int1_r, my_int2_r )

  USE DFT_MODULE
  USE EOS_VARIABLES, ONLY: PI, ncomp, mseg, dhs
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: j
  INTEGER, INTENT(IN)                    :: l
  INTEGER, INTENT(IN)                    :: m
  INTEGER, INTENT(IN OUT)                :: ih
  REAL, INTENT(IN)                       :: zz_1
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT,2)
  REAL, INTENT(IN)                       :: sigij
  REAL, INTENT(IN)                       :: dij
  REAL, INTENT(OUT)                      :: f_int_r
  REAL, INTENT(OUT)                      :: my_int1_r
  REAL, INTENT(OUT)                      :: my_int2_r

  !-----------------------------------------------------------------------------
  INTEGER                                :: k

  REAL                                   :: zz1
  REAL                                   :: dzr_local
  REAL                                   :: fint0, fint1
  REAL                                   :: myint1_0, myint1_1
  REAL                                   :: myint2_0, myint2_1
  REAL                                   :: dg_dz3, dg_dr
  REAL                                   :: rad, xg, rdf, z3_bar, ua, rs
  REAL                                   :: analytic1, tau_rs
  LOGICAL                                :: shortcut
  !-----------------------------------------------------------------------------

  shortcut = .true.
  fint0    = rc * tau_cut
  myint1_0 = rc * tau_cut
  myint2_0 = 0.0

  f_int_r   = 0.0
  my_int1_r = 0.0
  my_int2_r = 0.0

  !-----------------------------------------------------------------------------
  ! for mixtures it is advantageous to write all distances in dimensionless
  ! form, e.g. r^hat = r^hat / sigma.
  ! For mixtures, sigij = ( sig_i + sig_j ) / 2.
  !-----------------------------------------------------------------------------
  zz1 = zz_1 / sigij         ! here, dimensionless

  !-----------------------------------------------------------------------------
  ! this block only speeds up the integration
  !-----------------------------------------------------------------------------
  IF ( shortcut ) THEN
     rs = MAX( rg, zz1 ) ! +dzr
     IF ( rs > rc ) WRITE (*,*) 'error !!!!'
     analytic1 = 0.4*rs**-10 - 0.4*rc**-10 - rs**-4  + rc**-4
     f_int_r   = f_int_r   + analytic1
     my_int1_r = my_int1_r + analytic1
     IF ( rs == zz1 ) GO TO 10
     tau_rs   = 4.0 * ( rs**-12 - rs**-6 )
     fint0    = rs * tau_rs
     myint1_0 = rs * tau_rs
     rad = rs                      ! the simple integration scheme: set to rc
     k = 0 + NINT( (rc-rs)/dzr )   ! in simple scheme: set to 0
  ELSE
     rad = rc
     k = 0
  END IF
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !rad = rc
  !k = 0


  DO WHILE ( rad /= MAX( 1.0, zz1 ) )

     dzr_local = dzr                     ! dzr is set in f_dft (dimensionless)

     IF ( rad - dzr_local <= MAX( 1.0, zz1 ) ) THEN
        dzr_local = rad - MAX( 1.0, zz1 )
        rad = MAX( 1.0, zz1 )
     ELSE
        rad = rad - dzr_local
     END IF

     k = k + 1
     ua = 4.0 * ( rad**-12 -rad**-6 )
     xg = rad / dij * sigij
     !rho_bar = 0.5 * ( rhop(i) + rhop(j) )       * sigij**3
     z3_bar = PI / 6.0 * SUM(  0.5*( rhop(i,1:ncomp) + rhop(j,1:ncomp) )  * mseg(1:ncomp) * dhs(1:ncomp)**3 )
     rdf = 1.0
     dg_dz3 = 0.0
     IF ( rad <= rg ) THEN
        IF ( l == 1 .AND. m == 1 ) CALL BI_CUB_SPLINE (z3_bar,xg,ya_11,x1a_11,x2a_11,y1a_11,y2a_11,y12a_11,  &
             c_bicub_11,rdf,dg_dz3,dg_dr,den_step,ih,k)
        IF ( l /= m ) CALL BI_CUB_SPLINE (z3_bar,xg,ya_12,x1a_12,x2a_12,y1a_12,y2a_12,y12a_12,  &
             c_bicub_12,rdf,dg_dz3,dg_dr,den_step,ih,k)
        IF ( l == 2 .AND. m == 2 ) CALL BI_CUB_SPLINE (z3_bar,xg,ya_22,x1a_22,x2a_22,y1a_22,y2a_22,y12a_22,  &
             c_bicub_22,rdf,dg_dz3,dg_dr,den_step,ih,k)
     END IF

     fint1    = rdf * rad * ua
     myint1_1 = rdf * rad * ua
     myint2_1 = dg_dz3 * rad * ua

     f_int_r   = f_int_r   + dzr_local * ( fint1    + fint0   ) / 2.0
     my_int1_r = my_int1_r + dzr_local * ( myint1_1 + myint1_0) / 2.0
     my_int2_r = my_int2_r + dzr_local * ( myint2_1 + myint2_0) / 2.0

     fint0   = fint1
     myint1_0  = myint1_1
     myint2_0  = myint2_1

  ENDDO


10 CONTINUE

  ! analytic1 = 4.0/10.0*rc**-10  - rc**-4
  ! f_int_r   = f_int_r   + analytic1
  ! my_int1_r = my_int1_r + analytic1

  f_int_r   = f_int_r   * sigij*sigij
  my_int1_r = my_int1_r * sigij*sigij
  my_int2_r = my_int2_r * sigij*sigij

END SUBROUTINE DFT_RAD_INT_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE f_pert_theory_mix ( fdsp )

  USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, rho, eta, x, mseg, dhs, sig_ij, uij
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: fdsp

  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! constants
  !-----------------------------------------------------------------------------
  ua_c = 4.0 * ( rc**-12 - rc**-6 )
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

           ua = 4.0 * ( rad**-12 - rad**-6 )

           int11 = rdf * rad * rad * ua
           I1(l,m) = I1(l,m) + dzr_local * ( int11 + int10 ) / 2.0

           int10 = int11
           ! intgrid(k)= int11

        END DO

        ! stepno = k
        ! CALL SPLINE_PARA (dzr,intgrid,utri,stepno)
        ! CALL SPLINE_INT  (I1_spline,dzr,intgrid,utri,stepno)


        ! caution: 1st order integral is in F_EOS.f defined with negative sign
        !-----------------------------------------------------------------------
        ! cut-off corrections
        !-----------------------------------------------------------------------
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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE mu_pert_theory_mix ( mu_dsp )

  USE EOS_VARIABLES, ONLY: nc, PI, ncomp, t, rho, eta, x, mseg, dhs, sig_ij, uij
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN OUT)                   :: mu_dsp(nc)

  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! constants
  !-----------------------------------------------------------------------------
  ua_c = 4.0 * ( rc**-12 - rc**-6 )
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

           ua = 4.0 * ( rad**-12 - rad**-6 )

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


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE rdf_matrix_mix

  USE PARAMETERS, ONLY: PI
  USE EOS_VARIABLES, ONLY: parame, dhs
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                 :: i, j, k = 0
  REAL                                    :: rdf, rad, xg, rho_rdf, z3, mseg_ij, d_ij = 0.0
  !-----------------------------------------------------------------------------


  DO j = 1,3

     IF( j == 1) mseg_ij = 0.5 * ( parame(1,1) + parame(1,1) )
     IF( j == 2) mseg_ij = 0.5 * ( parame(1,1) + parame(2,1) )
     IF( j == 3) mseg_ij = 0.5 * ( parame(2,1) + parame(2,1) )

     DO i = 1,den_step
        ! rho_rdf= rhob(1) + (rhob(2)-rhob(1))*REAL(i)/den_step
        rho_rdf = 1.E-5 + (1.0)*REAL(i-1)/REAL(den_step-1)   ! segment density, rho_s*sigma**3
        rad = rc
        k = 0
        !write (*,*) 'eta=',rho_rdf*mseg_ij * PI/6.0  * dhs_st**3
        DO WHILE ( rad - 1.E-8 > 0.95 )
           rad = rad - dzr
           k = k + 1
           IF( j == 1) d_ij = dhs(1) / parame(1,2)
           IF( j == 2) d_ij = 0.5*(dhs(1)+dhs(2)) / ( 0.5*(parame(1,2)+parame(2,2)) )
           IF( j == 3) d_ij = dhs(2) / parame(2,2)

           xg = rad / d_ij
           z3 = rho_rdf * PI/6.0 * d_ij**3    ! d_ij is the dimensionless effective diameter d*(T)=d(T)/sigma

           IF ( j == 1) ya_11(i,k) = 1.0
           IF ( j == 2) ya_12(i,k) = 1.0
           IF ( j == 3) ya_22(i,k) = 1.0

           IF ( xg <= rg .AND. z3 > 0.0 ) CALL rdf_int ( z3, mseg_ij, xg, rdf )

           IF ( xg <= rg .AND. z3 > 0.0 .AND. j == 1 ) ya_11(i,k) = rdf
           IF ( xg <= rg .AND. z3 > 0.0 .AND. j == 2 ) ya_12(i,k) = rdf
           IF ( xg <= rg .AND. z3 > 0.0 .AND. j == 3 ) ya_22(i,k) = rdf
           ! ya(i,k) = y(x1a(i), x2a(k))  with x1a: density-vector, x2a: r-vector
           IF ( j == 1) x1a_11(i) = z3
           IF ( j == 1) x2a_11(k) = xg
           IF ( j == 2) x1a_12(i) = z3
           IF ( j == 2) x2a_12(k) = xg
           IF ( j == 3) x1a_22(i) = z3
           IF ( j == 3) x2a_22(k) = xg
        END DO
     END DO

     if ( xg > 1.0 ) stop 'rdf_matrix_mix: 0.95*sigma is too high for lower bound'

     WRITE (*,*) ' done with calculating g(r)',dhs_st

     kmax = k
     IF( j == 1) CALL bicub_derivative ( ya_11, x1a_11, x2a_11, y1a_11, y2a_11, y12a_11, den_step, kmax )
     IF( j == 1) CALL bicub_c ( ya_11, x1a_11, x2a_11, y1a_11, y2a_11, y12a_11, c_bicub_11, den_step, kmax )
     IF( j == 2) CALL bicub_derivative ( ya_12, x1a_12, x2a_12, y1a_12, y2a_12, y12a_12, den_step, kmax )
     IF( j == 2) CALL bicub_c ( ya_12, x1a_12, x2a_12, y1a_12, y2a_12, y12a_12, c_bicub_12, den_step, kmax )
     IF( j == 3) CALL bicub_derivative ( ya_22, x1a_22, x2a_22, y1a_22, y2a_22, y12a_22, den_step, kmax )
     IF( j == 3) CALL bicub_c ( ya_22, x1a_22, x2a_22, y1a_22, y2a_22, y12a_22, c_bicub_22, den_step, kmax )

  END DO

END SUBROUTINE rdf_matrix_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE dF_FMT_drho_mix ( i, fa, dzp, zp, phi_dn0,  &
     phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5, dF_drho_fmt )

  USE parameters, ONLY: PI
  USE basic_variables
  USE eos_variables, ONLY: dhs
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: fa(nc)
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn0(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn1(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn2(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn3(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn4(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn5(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: dF_drho_fmt(nc)

  !-----------------------------------------------------------------------------
  INTEGER                                :: j, k
  INTEGER                                :: nn, fa2
  REAL                                   :: int0, int1, int2, int3, int4, int5
  REAL                                   :: zz1, d2, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy0, hy1, hy2, hy3, hy4, hy5
  !-----------------------------------------------------------------------------


  DO k = 1, ncomp

     hx(:)  = 0.0
     hy0(:) = 0.0
     hy1(:) = 0.0
     hy2(:) = 0.0
     hy3(:) = 0.0
     hy4(:) = 0.0
     hy5(:) = 0.0
     y2(:)  = 0.0

     nn = 1
     d2 = dhs(k) / 2.0
     fa2 = ( fa(k) + 1 ) / 2
     DO j = (i-fa2), (i+fa2)
        IF ( zp(j+1) > (zp(i)-d2) .AND. (zp(i)-d2) >= zp(j) ) THEN
           hx(1)  = 0.0
           hy0(1) = phi_dn0(j) + ((zp(i)-d2)-zp(j))/dzp *( phi_dn0(j+1) - phi_dn0(j) )
           hy1(1) = phi_dn1(j) + ((zp(i)-d2)-zp(j))/dzp *( phi_dn1(j+1) - phi_dn1(j) )
           hy2(1) = phi_dn2(j) + ((zp(i)-d2)-zp(j))/dzp *( phi_dn2(j+1) - phi_dn2(j) )
           hy3(1) = 0.0 ! = rhop(i-fa/2)*parame(1,1)*(0.25-(zp(i-fa/2)-zp(i))**2)
           hy4(1) = phi_dn4(j)*(zp(j)-zp(i)) +((zp(i)-d2)-zp(j))/dzp  &
                * ( phi_dn4(j+1)*(zp(j+1)-zp(i)) - phi_dn4(j)*(zp(j)-zp(i)) )
           hy5(1) = phi_dn5(j)*(zp(j)-zp(i)) +((zp(i)-d2)-zp(j))/dzp  &
                * ( phi_dn5(j+1)*(zp(j+1)-zp(i)) - phi_dn5(j)*(zp(j)-zp(i)) )
           dz =  zp(j+1) - (zp(i)-d2)
        ELSE IF ( zp(j) > (zp(i)-d2).AND.zp(j) <= (zp(i)+d2) ) THEN
           zz1 = zp(j) - zp(i)                 ! distance z12 between 1 and 2
           nn = nn + 1
           hy0(nn) = phi_dn0(j)
           hy1(nn) = phi_dn1(j)
           hy2(nn) = phi_dn2(j)
           hy3(nn) = phi_dn3(j) * ( d2**2 - zz1**2 )
           hy4(nn) = phi_dn4(j) * zz1
           hy5(nn) = phi_dn5(j) * zz1
           hx(nn)  = hx(nn-1) + dz
           dz = dzp
           IF ( (zp(i)+d2) < zp(j+1) ) THEN
              dz = (zp(i)+d2) - zp(j)
              nn = nn + 1
              hx(nn)  = hx(nn-1) + dz
              hy0(nn) = phi_dn0(j) + ( (zp(i)+d2) - zp(j) ) / dzp * ( phi_dn0(j+1) - phi_dn0(j) )
              hy1(nn) = phi_dn1(j) + ( (zp(i)+d2) - zp(j) ) / dzp * ( phi_dn1(j+1) - phi_dn1(j) )
              hy2(nn) = phi_dn2(j) + ( (zp(i)+d2) - zp(j) ) / dzp * ( phi_dn2(j+1) - phi_dn2(j) )
              hy3(nn) = 0.0
              hy4(nn) = phi_dn4(j) * (zp(j)-zp(i)) + ((zp(i)+d2)-zp(j)) / dzp  &
                   * ( phi_dn4(j+1)*(zp(j+1)-zp(i)) -phi_dn4(j)*(zp(j)-zp(i)) )
              hy5(nn) = phi_dn5(j) * (zp(j)-zp(i)) + ((zp(i)+d2)-zp(j))/dzp  &
                   * ( phi_dn5(j+1)*(zp(j+1)-zp(i)) -phi_dn5(j)*(zp(j)-zp(i)) )
           END IF
        END IF
     END DO
     xl = hx(1)
     xh = hx(nn)
     CALL spline         ( hx(1:nn), hy0(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy0(1:nn), y2(1:nn), nn, xl, xh, int0 )
     CALL spline         ( hx(1:nn), hy1(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy1(1:nn), y2(1:nn), nn, xl, xh, int1 )
     CALL spline         ( hx(1:nn), hy2(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy2(1:nn), y2(1:nn), nn ,xl, xh, int2 )
     CALL spline         ( hx(1:nn), hy3(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy3(1:nn), y2(1:nn), nn, xl, xh, int3 )
     CALL spline         ( hx(1:nn), hy4(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy4(1:nn), y2(1:nn), nn, xl, xh, int4 )
     CALL spline         ( hx(1:nn), hy5(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy5(1:nn), y2(1:nn), nn, xl, xh, int5 )
     ! write (*,*) int3,int4,int5
     ! pause
     ! int0 = int0/d_hs
     ! int1 = 0.5*parame(1,2)*int1
     ! int2 = PI*d_hs*parame(1,2)**2 *int2
     ! int3 = PI*parame(1,2)**3 *int3
     ! int4 = parame(1,2)*int4
     ! int5 = 2.0*PI*parame(1,2)**2 *int5
     int0 = int0 / dhs(k)
     int1 = 0.5 * int1
     int2 = PI * dhs(k) * int2
     int3 = PI * int3
     int4 = int4 / dhs(k)
     int5 = 2.0 * PI * int5
     ! dF_drho_fmt=dF_hs/drho = dF_hs/drho_segment *m_mean
     dF_drho_fmt(k) = (int0+int1+int2+int3+int4+int5)*parame(k,1)
  END DO

END SUBROUTINE dF_FMT_drho_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE dF_chain_drho_mix
!
! dF/kT[rho] / d(rho_k(r)) of the chain term is here calculated as
! =     + (m_k - 1) * ln(rho_k)  -  (m_k - 1)* [ ln( y_kk(rho^bar)*lambda_k ) - 1 ]
!       - SUM_j[ (m_j - 1) * INTEGRAL[ d(ln(y_jj))/r(rho^bar_k) *rho_j * 0.75/d_k^3*(d_k^2-zp^2) dz ] ]
!       - (m_k - 1) * INTEGRAL[ 1/lambda_k * rho_k * 0.5/d_k  dz ]
! The second last line is abbreviated as I_ch1. The last line is
! abbreviated as - (m_k - 1)*I_ch1.
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE dF_chain_drho_mix ( i, fa, dzp, zp, rhop, lambda, rhobar, f_ch, dF_drho_ch )

  USE parameters, ONLY: PI
  USE basic_variables
  USE eos_variables, ONLY: dhs, mseg
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER                                :: fa(nc)
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT,2)
  REAL, INTENT(IN)                       :: lambda(-NDFT:NDFT,2)
  REAL, INTENT(IN)                       :: rhobar(-NDFT:NDFT,2)
  REAL, INTENT(OUT)                      :: f_ch
  REAL, INTENT(OUT)                      :: dF_drho_ch(nc)

  !-----------------------------------------------------------------------------
  INTEGER                                :: j, k, nn
  REAL                                   :: zz1, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy0, hy1
  REAL, DIMENSION(2)                     :: ycorr
  REAL, DIMENSION(2,2)                   :: dlnydr
  REAL                                   :: I_ch1, I_ch2
  !-----------------------------------------------------------------------------

  CALL cavity_mix ( rhobar(i,:), ycorr, dlnydr )     ! this is: d( ln( yij ) ) / d( rho(k) ) used with i=j

  f_ch = 0.0

  DO k = 1, ncomp
     hx(:)  = 0.0
     hy0(:) = 0.0
     hy1(:) = 0.0
     nn = 1
     DO j = (i-fa(k)), (i+fa(k))         ! fa=NINT(1./dzp)
        IF ( zp(j+1) > (zp(i)-dhs(k)) .AND. zp(j)-1E-11 <= (zp(i)-dhs(k)) ) THEN
           dz =  zp(j+1) - (zp(i)-dhs(k))
           hx(1)  = 0.0
           hy0(1) = 0.0  !  = 0.75*rhop(i-fa)*DLNYDR(d_hs,rhobar(i-fa))*(d.**2 -d.**2 )
           hy1(1) = 0.5/dhs(k)*rhop(j,k)/lambda(j,k) + ((zp(i)-dhs(k))-zp(j))/dzp  &
                * ( 0.5/dhs(k)*rhop(j+1,k)/lambda(j+1,k) - 0.5/dhs(k)*rhop(j,k)/lambda(j,k) )
        ELSE IF ( zp(j) > (zp(i)-dhs(k)) .AND. zp(j)+1.E-11 < (zp(i)+dhs(k)) ) THEN
           zz1 = zp(j) - zp(i)       ! distance z12 between 1 and 2
           nn = nn + 1

           hy0(nn) = - SUM( ( mseg(1:ncomp)-1.0 ) * rhop(j,1:ncomp)*dlnydr(1:ncomp,k) ) &
                * 0.75/dhs(k)**3 * (dhs(k)*dhs(k)-zz1*zz1)
           hy1(nn) = 0.5/dhs(k)*rhop(j,k)/lambda(j,k)
           hx(nn)  = hx(nn-1) + dz
           dz = dzp
           IF ( zp(j+1)+1.E-11 >= (zp(i)+dhs(k)) ) THEN
              dz = (zp(i)+dhs(k)) - zp(j)
              nn = nn + 1
              hy0(nn) = 0.0
              hy1(nn) = 0.5/dhs(k)*rhop(j,k)/lambda(j,k) +((zp(i)+dhs(k))-zp(j))/dzp  &
                   * ( 0.5/dhs(k)*rhop(j+1,k)/lambda(j+1,k)-0.5/dhs(k)*rhop(j,k)/lambda(j,k) )
              hx(nn) = hx(nn-1) + dz
           END IF
        END IF
     END DO
     xl = hx(1)
     xh = hx(nn)
     CALL spline         ( hx(1:nn), hy0(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy0(1:nn), y2(1:nn), nn, xl, xh, I_ch1 )
     CALL spline         ( hx(1:nn), hy1(1:nn), nn, 1.E30, 1.E30, y2(1:nn) )
     CALL splint_integral( hx(1:nn), hy1(1:nn), y2(1:nn), nn, xl, xh, I_ch2 )

     f_ch = f_ch + ( parame(k,1) - 1.0 ) * rhop(i,k) * ( LOG(rhop(i,k)) - 1.0 )  &
          - ( parame(k,1) - 1.0 ) * rhop(i,k) * ( LOG(ycorr(k)*lambda(i,k)) - 1.0 )
     dF_drho_ch(k) = + ( parame(k,1) - 1.0 ) * LOG(rhop(i,k))  &
          - ( parame(k,1) - 1.0 ) * ( LOG( ycorr(k)*lambda(i,k) ) - 1.0 + I_ch2 )  &
          + I_ch1
  END DO

END SUBROUTINE dF_chain_drho_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE cavity_mix ( rhoi, ycorr, dlnydr )

  USE parameters, ONLY: PI
  USE basic_variables
  USE eos_variables, ONLY: dhs, mseg
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rhoi(2)
  REAL, INTENT(OUT)                      :: ycorr(2)
  REAL, INTENT(OUT)                      :: dlnydr(2,2)     ! this is: d( ln( yij ) ) / d( rho(k) ) used with i=j

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, k
  REAL                                   :: z0, z1, z2, z3, zms, z1_rk, z2_rk, z3_rk
  REAL, DIMENSION(nc,nc)                 :: dij_ab, gij, gij_rk
  !-----------------------------------------------------------------------------

  z0 = PI / 6.0 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) )
  z1 = PI / 6.0 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp) )
  z2 = PI / 6.0 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**2 )
  z3 = PI / 6.0 * SUM( rhoi(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  zms    = 1.0 - z3

  DO i = 1,ncomp
     DO j=1,ncomp
        dij_ab(i,j)=dhs(i)*dhs(j)/(dhs(i)+dhs(j))
     ENDDO
  END DO

  DO k = 1, ncomp
     DO i = 1, ncomp
        z1_rk = PI/6.0 * mseg(k) * dhs(k)
        z2_rk = PI/6.0 * mseg(k) * dhs(k)*dhs(k)
        z3_rk = PI/6.0 * mseg(k) * dhs(k)**3
        !DO j = 1, ncomp
        j = i
        gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms  &
             + 2.0*(dij_ab(i,j)*z2)**2 /zms**3
        !dgijdz(i,j)= 1.0/zms/zms + 3.0*dij_ab(i,j)*z2*(1.0+z3)/z3/zms**3   &
        !           + (dij_ab(i,j)*z2/zms/zms)**2 *(4.0+2.0*z3)/z3
        gij_rk(i,j) = z3_rk/zms/zms  &
             + 3.0*dij_ab(i,j)*(z2_rk+2.0*z2*z3_rk/zms)/zms/zms  &
             + dij_ab(i,j)**2 *z2/zms**3  *(4.0*z2_rk+6.0*z2*z3_rk/zms)
        !END DO

        ycorr(i)  = gij(i,i)
        dlnydr(i,k) = gij_rk(i,i) / gij(i,i)

     END DO
  END DO

END SUBROUTINE cavity_mix


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE cutoff_corrections_mix ( i, irc, rhop, rhob, f_att, dF_drho_att )

  USE EOS_VARIABLES, ONLY: nc, ncomp, PI, mseg, sig_ij, uij, t
  USE DFT_Module, ONLY: rc, NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: irc(nc)
  REAL, DIMENSION(-NDFT:NDFT,2), INTENT(IN)     :: rhop
  REAL, DIMENSION(2,0:nc), INTENT(IN)    :: rhob
  REAL, INTENT(IN OUT)                   :: f_att
  REAL, DIMENSION(2), INTENT(IN OUT)     :: dF_drho_att

  !-----------------------------------------------------------------------------
  INTEGER                                :: l, m
  INTEGER                                :: irc_j
  REAL                                   :: cutoffz, rho_l, rho_m
  !-----------------------------------------------------------------------------

  cutoffz = 4.0/90.0 * rc**-9 - 1.0/3.0 * rc**-3

  !-----------------------------------------------------------------------------
  ! cutoff-corrections to the left
  !-----------------------------------------------------------------------------
  DO l = 1, ncomp
     DO m = 1, ncomp
        irc_j = ( irc(l) + irc(m) + 1 ) / 2
        rho_l = rhob(1,l)
        rho_m = rhob(1,m)
        IF ( ( ABS(rhop(i-irc_j,l)-rhob(1,l))/rhob(1,l) ) > 1.E-3 ) rho_l = rhop(i-irc_j,l)  ! a crude approximation !
        IF ( ( ABS(rhop(i-irc_j,m)-rhob(1,m))/rhob(1,m) ) > 1.E-3 ) rho_m = rhop(i-irc_j,m)  ! a crude approximation !
        write (*,*) 'le',i, rho_l,rhop(i,l)
        write (*,*) 'le',i, rho_m,rhop(i,m)

        f_att  = f_att + PI*mseg(l)*mseg(m) *rho_l* rho_m *sig_ij(l,m) *uij(l,m)/t *cutoffz
        write(*,*) 'AAA', l,m,f_att,PI*mseg(l)*mseg(m) *rho_l* rho_m *sig_ij(l,m) *uij(l,m)/t *cutoffz
        dF_drho_att(l) = dF_drho_att(l) + 2.0*PI*mseg(l)*mseg(m) *rho_m *sig_ij(l,m) *uij(l,m)/t *cutoffz
     END DO
  END DO

  !-----------------------------------------------------------------------------
  ! cutoff-corrections to the right
  !-----------------------------------------------------------------------------
  DO l = 1, ncomp
     DO m = 1, ncomp
        irc_j = ( irc(l) + irc(m) + 1 ) / 2
        rho_l = rhob(2,l)
        rho_m = rhob(2,m)
        IF ( ( ABS(rhop(i+irc_j,l)-rhob(2,l))/rhob(2,l) ) > 1.E-3 ) rho_l = rhop(i+irc_j,l)  ! a crude approximation !
        IF ( ( ABS(rhop(i+irc_j,m)-rhob(2,m))/rhob(2,m) ) > 1.E-3 ) rho_m = rhop(i+irc_j,m)  ! a crude approximation !

        f_att  = f_att + PI*mseg(l)*mseg(m) *rho_l* rho_m *sig_ij(l,m)**3 *uij(l,m)/t *cutoffz
        dF_drho_att(l) = dF_drho_att(l) + 2.0*PI*mseg(l)*mseg(m) *rho_m *sig_ij(l,m)**3 *uij(l,m)/t *cutoffz
        write(*,*) 'BBB', l,m,f_att,PI*mseg(l)*mseg(m) *rho_l* rho_m *sig_ij(l,m) *uij(l,m)/t *cutoffz
     END DO
  END DO

END SUBROUTINE cutoff_corrections_mix
