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

SUBROUTINE DFT_nMFT_units

  USE parameters, ONLY: PI, RGAS, KBOL
  USE basic_variables
  USE EOS_VARIABLES, ONLY: fres, eta, rho, x, z3t
  USE DFT_MODULE
  use utilities
  use EOS_polar, only: f_polar, phi_polar
  USE EOS_NUMERICAL_DERIVATIVES

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, kk, converg
  INTEGER                                :: discret, fa, outpt, irc, maxiter, ih
  LOGICAL                                :: diagram
  REAL, DIMENSION(-NDFT:NDFT)            :: zp, rhop, rhop_o
  REAL, DIMENSION(-NDFT:NDFT)            :: dF_drho_tot, f_tot
  REAL, DIMENSION(-NDFT:NDFT)            :: dF_drho_att, f_att
  REAL, DIMENSION(2)                     :: rhob
  REAL                                   :: f_disp_1PT, f_disp_pcsaft
  REAL                                   :: mu_disp_1PT, mu_disp_pcsaft
  REAL                                   :: zges(np), pbulk
  REAL                                   :: msegm, delta_st, tanhfac
  REAL                                   :: end_x, steps, my0
  REAL                                   :: f_int_z, mu_int_z
  REAL                                   :: f_int_r, mu_int_r
  REAL                                   :: f_int2_z, mu_int2_z, mu_int3_z
  REAL                                   :: f_int2_r, mu_int2_r, mu_int3_r
  REAL                                   :: zz1
  REAL                                   :: dz_local
  REAL                                   :: dev, maxdev
  REAL                                   :: delta_f, free
  REAL                                   :: surftens(0:200), st_macro(200)
  REAL                                   :: f01, f02, f03, f04, f05
  REAL                                   :: c1_con, c2_con
  REAL                                   :: zms, z3
  REAL                                   :: tsav, psav, tc, pc, rhoc
  REAL                                   :: density(np), w(np,nc), lnphi(np,nc)
  REAL                                   :: damppic
  REAL                                   :: box_l_no_unit
  REAL, DIMENSION(-NDFT:NDFT)            :: lambda, rhobar, phi_dn0, phi_dn1
  REAL, DIMENSION(-NDFT:NDFT)            :: phi_dn2, phi_dn3, phi_dn4, phi_dn5
  REAL, DIMENSION(0:5,-NDFT:NDFT)        :: ni
  REAL                                   :: zs
  REAL                                   :: f_ch, dF_drho_ch
  REAL                                   :: f_fmt, dF_drho_fmt
  REAL                                   :: mu_assoc, f_assoc
  REAL                                   :: fres_polar, fdd, fqq, fdq
  REAL                                   :: mu_polar, fdd_rk, fqq_rk, fdq_rk, z3_rk
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
  msegm = parame(1,1)


  d_hs = parame(1,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1,3)/t ) )
  dhs_st = d_hs/parame(1,2)
  z3t_st = PI/6.0* parame(1,1) * d_hs**3

  IF ( ncomp /= 1 ) THEN
     write (*,*)'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE:'
     write (*,*)'    ./input_file/INPUT.INP'
     stop
  END IF
  OPEN (68,FILE = './output_file/DFT_profiles.xlo')
  OPEN (69,FILE = './output_file/DFT_sigma.xlo')
  OPEN (71,FILE = './output_file/DFT_iteration.xlo')
  OPEN (72,FILE = './output_file/DFT_h.xlo')


  !-----------------------------------------------------------------------------
  ! The cut-off distance rc is the max. distance for integrating inter-
  ! actions. Beyond rc, tail-corrections are added.
  !-----------------------------------------------------------------------------
  rc = 9.0                      ! dimensionless cut-off distance rc = r_c/sigma


  !-----------------------------------------------------------------------------
  ! basic definitions for calculating density profile,
  !-----------------------------------------------------------------------------
  ! grid-size  dzp = zp(1)-zp(0)

  box_l_no_unit = 160.0         ! lenth of simulation box (Angstrom)
  discret= 1000                 ! number of spacial discretizations for profile
  dzp    = box_l_no_unit / REAL(discret)    ! grid-distance (Angstrom)
  fa     = NINT( parame(1,2) / dzp + 1 )    ! number of steps per sigma (currently of component 1)
  outpt  = 100                  ! number output-points = discret/outpt


  !-----------------------------------------------------------------------------
  ! definitions for the numerical algorithm
  !-----------------------------------------------------------------------------

  maxiter = 100                  ! maximum number of iterations per temp.step
  maxdev  = 2.E-7               ! maximum deviation
  damppic = 0.005               ! damping of Picard-iteration


  !-----------------------------------------------------------------------------
  ! get a matrix of values for the pair correlation function g(r,rho)
  !-----------------------------------------------------------------------------

  rg = 4.0
  den_step = 40
  !dzr = dzp / 2.0               ! dimensionless grid-distance for integrating
  !                              ! the Barker-Henderson attraction term
  dzr = 0.025                   ! dimensionless grid-distance for integrating
  ! the Barker-Henderson attraction term
  CALL rdf_matrix_units (msegm)


  !-----------------------------------------------------------------------------
  ! prepare for phase equilibrium calculation for given T
  !-----------------------------------------------------------------------------

  Diagram_for_various_T_Loop: DO

     nphas  = 2
     n_unkw = ncomp                ! number of quantities to be iterated
     it(1)  = 'p'                  ! iteration of pressure

     val_init(1) = 0.45            ! starting value for liquid density
     val_init(2) = 1.E-5           ! starting value for vapor density
     val_init(3) = t               ! value of temperature
     val_init(4) = 10.0            ! default starting value for p in (Pa)
     IF ( eos >= 4 ) val_init(4) = 1.E-3       ! default starting value for LJ-model
     val_init(5) = 0.0             ! logar. mole fraction: lnx=0, x=1, phase 1
     val_init(6) = 0.0             ! logar. mole fraction: lnx=0, x=1, phase 2

     running = 't'                 ! T is running variable in PHASE_EQUILIB - here T=const.
     end_x   = t                   ! end temperature - here T=const.
     steps   = 1.0                 ! number of steps of running var. - here 1, since T=const.

     outp = 0                      ! output to terminal


     !--------------------------------------------------------------------------
     ! calculate phase equilibrium for given T
     !--------------------------------------------------------------------------

     converg = 0
     DO WHILE ( converg == 0 )
        ensemble_flag = 'tp'                 ! this flag is for: 'regular calculation'
        CALL phase_equilib( end_x, steps, converg )
        IF ( converg /= 1 ) THEN
           val_conv(2) = 0.0
           val_init(3) = t - 10.0            ! value of temperature
           IF ( val_conv(4) /= 0.0 ) val_init(4) = val_conv(4) / 4.0
           ! caution: the intermediate points calculated here are not correct, because
           ! the temperature-dependent props like d_BH, g(r), or z3t are not evaluated
           ! for every temp-step. The final converged result, however, is correct.
           steps = 2.0  ! number of steps of running var.
           WRITE (*,*) 'no VLE found'
           IF ( val_init(3) <= 10.0 ) RETURN
        END IF
     END DO

     ! rhob: molecular density times sigma**3 (rho_molec*sigma**3)
     rhob(1) = dense(1)/z3t_st               ! coexisting bulk density L
     rhob(2) = dense(2)/z3t_st               ! coexisting bulk density V
     WRITE (*,*) 'temperature  ',t,p/1.E5
     WRITE (*,*) 'densities    ',rhob(1), rhob(2)
     WRITE (*,*) 'dense        ',dense(1),dense(2)


     ! --- get density in SI-units (kg/m**3) -----------------------------------
     CALL SI_DENS ( density, w )


     ! --- (re-)calculate residual chemical potential of both phases -----------
     ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)
     densta(1) = dense(1)                    ! Index 1 is for liquid density (here: packing fraction eta)
     densta(2) = dense(2)                    ! Index 2 is for vapour density (here: packing fraction eta)
     CALL fugacity (lnphi)
     my0 = lnphi(1,1) + LOG( rhob(1) )     ! my0 = mu_res(T,rho_bulk_L) + ln(rho_bulk_l)
     zges(1) = p / ( RGAS*t*density(1) ) * mm(1)/1000.0
     zges(2) = p / ( RGAS*t*density(2) ) * mm(1)/1000.0

     pbulk = zges(1) * rhob(1)               ! pressure  p/kT (= Z*rho)
     WRITE (*,*) 'chem. potentials', lnphi(1,1) + LOG( rhob(1) ),  &
          lnphi(2,1) + LOG( rhob(2) )
     WRITE (*,*) ' '
     tc = t
     IF (num == 2) WRITE (*,*) 'enter an estimate for crit. temp.'
     IF (num == 2) READ (*,*) tc
     ensemble_flag = 'tp'                    ! this flag is for: 'regular calculation'
     tsav = t
     psav = p
     IF ( eos < 4 ) CALL critical (tc,pc,rhoc)
     !IF ( eos >= 4 .AND. eos /= 7 ) CALL lj_critical (tc,pc,rhoc)
     t = tsav
     p = psav
     WRITE (*,'(a,3(f16.4))') 'critical point',tc, pc/1.E5, rhoc
     WRITE (*,*) ' '

     !--------------------------------------------------------------------------
     ! update z3t, the T-dependent quantity that relates eta and rho, as eta = z3t*rho
     !--------------------------------------------------------------------------
     CALL PERTURBATION_PARAMETER


     !--------------------------------------------------------------------------
     ! define initial density profile rhop(i)
     ! and dimensionless space coordinates zp(i)
     !
     ! discret  : number of grid-points within "the box"
     ! irc      : number of grid-points extending the the box to left and right
     !            the box is extended in order to allow for numerical integration
     !            'irc' is determined by the cut-off distance 'rc'
     !--------------------------------------------------------------------------

     irc = NINT(rc*parame(1,2)/dzp) + 1
     tanhfac = -2.3625*t/tc + 2.4728
     ! tanhfac = 2.7*(1.0 - t/tc)   ! this parameterization was published (Gross, JCP, 2009)
     DO i = -irc, (discret+irc)
        zp(i) = REAL(i) * dzp
     END DO
     DO i = -irc, (discret+irc)
        rhop(i) = TANH(-(zp(i)-zp(INT(discret/2))) / parame(1,2) *tanhfac) * (rhob(1)-rhob(2))/2.0  &
             + (rhob(1)+rhob(2))/2.0
        ! rhop(i) = rhob(1)
     END DO


     !--------------------------------------------------------------------------
     ! Initialize the DENS_INV subroutine
     !--------------------------------------------------------------------------

     nphas = 1

     ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)



     !==========================================================================
     ! Start iterating the density profile
     !==========================================================================

     kk = 1
     ih = 85

     DFT_Convergence_Loop: DO

        !-----------------------------------------------------------------------
        ! Getting auxilliary quantities along the profile
        !-----------------------------------------------------------------------

        CALL aux_chain ( discret, fa, dzp, d_hs, zp, rhop, rhobar, lambda )

        CALL fmt_dens ( discret, fa, dzp, d_hs, zp, rhop, ni,  &
             phi_dn0, phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5 )


        !-----------------------------------------------------------------------
        ! Start loop for density profile
        !-----------------------------------------------------------------------

        DO i = 0, discret

           !--------------------------------------------------------------------
           ! hard sphere contribution
           !    f_fmt  is the Helmholtz energy density F_FMT/(VkT) = F_FMT/(NkT)*rho
           !    dF_drho_fmt is the functional derivative of F_FMT/(kT) to rho
           !--------------------------------------------------------------------
           zs = 1.0 - ni(3,i)
           f_fmt  = - ni(0,i)*LOG(zs) + ni(1,i)*ni(2,i)/zs - ni(4,i)*ni(5,i)/zs  &
                + (ni(2,i)**3 -3.0*ni(2,i)*ni(5,i)*ni(5,i)) *(ni(3,i)+zs*zs*LOG(zs))  &
                /36.0/PI/zs/zs/ni(3,i)**2

           CALL dF_FMT_dr ( i, fa, dzp, d_hs, zp, phi_dn0,  &
                phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5, dF_drho_fmt )


           !--------------------------------------------------------------------
           ! chain term
           !--------------------------------------------------------------------
           CALL dF_chain_dr ( i, fa, dzp, d_hs, zp, rhop, lambda,  &
                rhobar, z3t_st, f_ch, dF_drho_ch )


           !--------------------------------------------------------------------
           ! dispersive attraction
           !--------------------------------------------------------------------
           f_int_z   = 0.0
           f_int2_z  = 0.0
           mu_int_z  = 0.0
           mu_int2_z = 0.0
           mu_int3_z = 0.0
           f01 = 0.0
           f02 = 0.0
           f03 = 0.0
           f04 = 0.0
           f05 = 0.0
           DO j = (i-irc), (i+irc)             ! integration in z-coordinate
              zz1 = ABS( zp(j) - zp(i) )       ! distance z12 between 1 and 2
              dz_local = dzp
              IF ( zz1 < rc*parame(1,2) ) THEN
                 IF ( zz1 >= rc*parame(1,2) - dzp ) dz_local = rc*parame(1,2) - zz1
                 CALL dft_rad_int_units ( i, j, ih, zz1, rhop, f_int_r, mu_int_r,  &
                      f_int2_r, mu_int2_r, mu_int3_r )
              ELSE
                 f_int_r  = 0.0
                 f_int2_r = 0.0
                 mu_int_r = 0.0
                 mu_int2_r= 0.0
                 mu_int3_r= 0.0
              END IF

              f_int_z  = f_int_z  + dz_local * (rhop(j)* f_int_r + f01)/2.0
              mu_int_z = mu_int_z + dz_local * (rhop(j)*mu_int_r + f02)/2.0
              f_int2_z = f_int2_z + dz_local * ( f_int2_r + f03)/2.0
              mu_int2_z= mu_int2_z+ dz_local * (mu_int2_r + f04)/2.0
              mu_int3_z= mu_int3_z+ dz_local * (mu_int3_r + f05)/2.0
              f01 = rhop(j)* f_int_r
              f02 = rhop(j)* mu_int_r
              f03 = f_int2_r
              f04 = mu_int2_r
              f05 = mu_int3_r
           END DO


           !--------------------------------------------------------------------
           ! cut-off corrections
           !--------------------------------------------------------------------

           CALL cutoff_corrections_units ( i, irc, rhop, rhob, f_int_z, mu_int_z, f_int2_z, mu_int2_z )


           !--------------------------------------------------------------------
           ! for second order dispersive attraction ( usually not used !!! )
           !--------------------------------------------------------------------

           z3 = rhop(i) * z3t_st
           zms = 1.0 - z3
           c1_con= 1.0/ (  1.0 + parame(1,1)*(8.0*z3-2.0*z3**2 )/zms**4   &
                + (1.0 - parame(1,1))*(20.0*z3-27.0*z3**2   &
                + 12.0*z3**3 -2.0*z3**4 ) /(zms*(2.0-z3))**2   )
           c2_con= - c1_con*c1_con  &
                * ( parame(1,1)*(-4.0*z3**2 +20.0*z3+8.0)/zms**5   &
                + (1.0 - parame(1,1)) *(2.0*z3**3 +12.0*z3**2 -48.0*z3+40.0)  &
                / (zms*(2.0-z3))**3  )
           mu_int2_z = mu_int2_z /4.0 * ( 2.0*rhop(i)*c1_con + rhop(i)*z3*c2_con )
           mu_int3_z = mu_int3_z /4.0  *rhop(i)*rhop(i)*c1_con
           f_int2_z  = f_int2_z  /4.0  *rhop(i)*rhop(i)*c1_con

           !--------------------------------------------------------------------
           ! The integration of the DFT-integral can be done with cubic splines
           !--------------------------------------------------------------------
           ! stepno = discret + irc*2
           ! CALL SPLINE_PARA (dzp, intgrid, utri, stepno)
           ! CALL SPLINE_INT (f_int_z, dzp, intgrid, utri, stepno)
           ! f_int_z = f_int_z + rhob(1)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )
           ! f_int_z = f_int_z + rhob(2)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )


           !--------------------------------------------------------------------
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
           !--------------------------------------------------------------------
           ! remember: factor 2*PI because of the cylindrical coordinates

           dF_drho_att(i) = 2.0*PI*parame(1,1)**2 * parame(1,3)/t * mu_int_z
           f_att(i) = PI*parame(1,1)**2 *parame(1,3)/t*rhop(i)*f_int_z


           ! --- if a second order perturbation theory is considered -----------
           IF (subtract1 =='2PT') THEN
              dF_drho_att(i)= dF_drho_att(i) -2.0*PI*parame(1,1)**3 *(parame(1,3)/t)**2 *(mu_int2_z+mu_int3_z)
              f_att(i) = f_att(i)  -2.0*PI*parame(1,1)**3 *(parame(1,3)/t)**2 *f_int2_z
           END IF

           !--------------------------------------------------------------------
           ! make the dispersive attraction term consistent with PC-SAFT, by
           ! adding the difference ( PC-SAFT - 1PT )_dispersion locally (LDA)
           !--------------------------------------------------------------------

           !*****************************************************************************
           ! under construction
           !*****************************************************************************
           dense(1) = rhop(i) * z3t_st
           densta(1) = dense(1)

           ensemble_flag = 'tv'

           call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PT1      ' )
           call FUGACITY ( lnphi )
           call RESTORE_PREVIOUS_EOS_NUMERICAL
           f_disp_1PT  = fres*rhop(i)
           mu_disp_1PT = lnphi(1,1)

           call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PC-SAFT  ' )
           call FUGACITY ( lnphi )
           call RESTORE_PREVIOUS_EOS_NUMERICAL
           f_disp_pcsaft  = fres*rhop(i)
           mu_disp_pcsaft = lnphi(1,1)

           eta = densta(1)
           rho = eta/z3t
           z3_rk = z3t         ! only for pure substances
           x(1) = 1.0
           call F_POLAR ( fdd, fqq, fdq )
           call PHI_POLAR ( 1, z3_rk, fdd_rk, fqq_rk, fdq_rk )
           fres_polar  = ( fdd + fqq + fdq ) * rhop(i)
           mu_polar = fdd_rk + fqq_rk + fdq_rk

           f_assoc  = 0.0
           mu_assoc = 0.0
           IF ( SUM( NINT(parame(1:ncomp,12)) ) > 0) THEN
              call ONLY_ONE_TERM_EOS_NUMERICAL ( 'hb_term  ', 'TPT1_Chap' )
              call FUGACITY ( lnphi )
              call RESTORE_PREVIOUS_EOS_NUMERICAL
              f_assoc  = fres*rhop(i)
              mu_assoc = lnphi(1,1)
           END IF

           ensemble_flag = 'tp'

           dF_drho_att(i)= dF_drho_att(i) + ( mu_disp_pcsaft - mu_disp_1PT )  + mu_assoc + mu_polar
           f_att(i)      = f_att(i)       + ( f_disp_pcsaft  - f_disp_1PT  )  + f_assoc + fres_polar
           !*****************************************************************************
           !*****************************************************************************


           !--------------------------------------------------------------------
           ! collect the total Helmholtz energy density 'f_tot' and the
           ! functional derivative to rhop(z) 'dF_drho_tot' (including ideal gas)
           ! For f_tot, it is numerically advantageous to add the ideal gas term
           ! after updating rhop(i)
           !--------------------------------------------------------------------

           dF_drho_tot(i) = LOG( rhop(i) ) + dF_drho_fmt + dF_drho_ch + dF_drho_att(i)
           f_tot(i) = f_fmt + f_ch + f_att(i)


           ! --- on-the-fly report on local errors in the profile --------------
           IF ( MOD(i, outpt) == 0 )write(*,'(i7,2(G15.6))') i,rhop(i), my0 - dF_drho_tot(i)

        END DO          ! end of loop (i = 0, discret) along the profile


        !-----------------------------------------------------------------------
        ! update the density profile using either a Picard-iteration scheme
        ! (i.e. a direct substitution scheme, or a LDA - Inversion Procedure
        ! similar to R.Evans (Bristol), which is done in function DENS_INVERS2.
        !-----------------------------------------------------------------------
        dev = 0.0
        DO i = 0, discret
           rhop_o(i) = rhop(i)
           rhop(i)   = rhop(i) * EXP( my0 - dF_drho_tot(i) )
           rhop(i)   = rhop_o(i) + (rhop(i)-rhop_o(i))*damppic
           dev = dev + ( (rhop(i)-rhop_o(i))*parame(1,2)**3 )**2
        END DO


        !-----------------------------------------------------------------------
        ! calculate surface tension
        !-----------------------------------------------------------------------
        free = 0.0
        DO i = 1, discret
           f_tot(i) = f_tot(i) + rhop(i)*( LOG(rhop(i))-1.0 )  ! now add the ideal gas term
           delta_f = f_tot(i)  -   (rhop(i)*my0 - pbulk)  ! all quantities .../(kT)
           free  = free  + delta_f*dzp
        END DO
        surftens(kk) = KBOL * t *1.E20*1000.0 *free


        !-----------------------------------------------------------------------
        ! add an approximate capillary wave contribution to the surface tension
        !-----------------------------------------------------------------------
        st_macro(kk) = surftens(kk) / ( 1.0 + 3.0/8.0/PI *t/tc  &
             * (1.0/2.55)**2  / (0.0674*parame(1,1)+0.0045) )


        delta_st = 1.0
        IF ( kk > 1 ) delta_st = ABS( surftens(kk)-surftens(kk-1) ) / surftens(kk)

        WRITE (*,*)  '-----------------------------------------------------------'
        WRITE (*,*)  ' #      error-sum      intrinsic ST    total ST'
        WRITE (*,'(i3,3(F16.8))')  kk, dev, surftens(kk), st_macro(kk)
        WRITE (*,*)  '-----------------------------------------------------------'
        WRITE (71,'(i3,4(E18.8))') kk, dev, surftens(kk), st_macro(kk)
        kk = kk + 1

        !-----------------------------------------------------------------------
        ! add an approximate capillary wave contribution to the surface tension
        !-----------------------------------------------------------------------
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

     !--------------------------------------------------------------------------
     ! write resulting density profile
     !--------------------------------------------------------------------------
     WRITE(68,'(i3,6(f18.8))') 0, zp(0)-zp(INT(discret/2)), rhop(0), t,  &
          val_conv(4)/1.E5, density(1), density(2)
     DO i = 0, discret
        IF ( MOD(i, outpt) == 0 ) WRITE (68,'(i6,3(f18.12))') i,zp(i)-zp(INT(discret/2)),  &
             rhop(i), -( dF_drho_tot(i)-my0 )
        ! IF ( MOD(i, outpt) == 0 ) WRITE(*,'(i6,3(f18.12))') i,zp(i)-zp(INT(discret/2)),rhop(i)
        ! write (69,*) zp(i)*parame(1,2), pN_pT(i)/parame(1,2)**3 /1E-30*parame(1,3)*KBOL/1.E6 ! in MPa
     END DO
     WRITE (68,*) ' '


     !--------------------------------------------------------------------------
     ! summarize results
     !--------------------------------------------------------------------------
     WRITE (*,*) ' '
     WRITE (*,*) 'SUMMARY FOR A SINGLE TEMPERATURE'
     WRITE (*,*) ' '
     WRITE (*,*) 'Temp. [K], Pressure [bar]    ',t,val_conv(4)/1.E5
     WRITE (*,*) 'Critical point Temp., Press. ',tc,pc/1.E5
     WRITE (*,*) 'Density [kg/m**3]            ',density(1),density(2)
     WRITE (*,*) 'Dimensionless Density (rho*) ',rhob(1),rhob(2)
     WRITE (*,*) 'Excess Grand Potential       ',free
     WRITE (*,*) 'Intrinsic Interf. Tension [mN/m] ',surftens(kk-1)
     WRITE (*,*) 'Macroscop.Interf. Tension [mN/m] ',st_macro(kk-1)
     WRITE (*,*) '============================================================'
     WRITE (*,*) ' '
     WRITE (69,'(9(f18.10))') t, val_conv(4)/1.E5,  &
          rhob(1),rhob(2),surftens(kk-1),st_macro(kk-1),free,dev

     !--------------------------------------------------------------------------
     ! when calc. a phase diagram & diagram of surface tens., loop for higher T
     !--------------------------------------------------------------------------
     ensemble_flag = 'tp'          ! this flag is for 'regular calculations'
     subtract1     = 'no'          ! this flag is for 'regular calculations'
     IF ( diagram ) THEN
        IF ( (t+8.0) <= tc ) THEN
           t = t + 5.0
           IF ( (t+15.0) <= tc ) t = t + 5.0
           IF ( (t+25.0) <= tc ) t = t + 10.0
           IF ( (t+45.0) <= tc ) t = t + 20.0
           nphas = 2
           n_unkw = ncomp            ! number of quantities to be iterated
           it(1) = 'p'               ! iteration of pressure
           val_init(3) = t           ! value of temperature
           running = 't'             ! T is running variable in PHASE_EQUILIB - here T=const.
           end_x  = t                ! end temperature - here T=const.
           steps = 1.0               ! number of steps of running var. - here 1, since T=const.
           d_hs = parame(1,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1,3)/t ) )
           z3t_st = PI/6.0* parame(1,1) * d_hs**3
           dhs_st = d_hs/parame(1,2)
           CALL rdf_matrix_units (msegm)
        ELSE
           EXIT Diagram_for_various_T_Loop
        END IF
     ELSE
        EXIT Diagram_for_various_T_Loop
     END IF

  ENDDO Diagram_for_various_T_Loop
  WRITE (69,'(7(f18.10))') tc, pc/1.E5, rhoc, rhoc, 0., 0., 0.

END SUBROUTINE DFT_nMFT_units





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE cutoff_corrections_units ( i, irc, rhop, rhob, f_int_z, mu_int_z, f_int2_z, mu_int2_z )

  USE DFT_Module, ONLY: NDFT, rc
  USE EOS_VARIABLES, ONLY: parame
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: irc
  REAL, DIMENSION(-NDFT:NDFT), INTENT(IN)     :: rhop
  REAL, DIMENSION(2), INTENT(IN)         :: rhob
  REAL, INTENT(IN OUT)                   :: f_int_z, f_int2_z
  REAL, INTENT(IN OUT)                   :: mu_int_z, mu_int2_z

  !-----------------------------------------------------------------------------
  REAL                                   :: cutoffz, cutoffz2, rhotemp
  !-----------------------------------------------------------------------------

  cutoffz  = ( 4.0/90.0 * rc**-9 - 1.0/3.0 * rc**-3 ) * parame(1,2)**3
  cutoffz2 = ( 16.0/22.0/21.0 * rc**-21 - 2.0/15.0 * rc**-15 + 16.0/90.0 * rc**-9 ) * parame(1,2)**3
  IF ( ABS( rhop(i+irc)-rhob(2) ) > 1.E-5 ) THEN
     rhotemp = rhop(i+irc)               ! this is a crude approximation !
     ! rhotemp = rhob(2)
     f_int_z  = f_int_z  + rhotemp*cutoffz
     mu_int_z = mu_int_z + rhotemp*cutoffz
     f_int2_z = f_int2_z + cutoffz2
     mu_int2_z= mu_int2_z+ cutoffz2
  ELSE
     f_int_z  = f_int_z  + rhob(2)*cutoffz
     mu_int_z = mu_int_z + rhob(2)*cutoffz
     f_int2_z = f_int2_z + cutoffz2
     mu_int2_z= mu_int2_z+ cutoffz2
  END IF
  IF ( ABS( rhop(i-irc)-rhob(1) ) > 1.E-3 ) THEN
     rhotemp = rhop(i-irc)
     ! rhotemp=rhob(1)
     f_int_z  = f_int_z  + rhotemp*cutoffz
     mu_int_z = mu_int_z + rhotemp*cutoffz
     f_int2_z = f_int2_z + cutoffz2
     mu_int2_z= mu_int2_z+ cutoffz2
  ELSE
     f_int_z  = f_int_z  + rhob(1)*cutoffz
     mu_int_z = mu_int_z + rhob(1)*cutoffz
     f_int2_z = f_int2_z + cutoffz2
     mu_int2_z= mu_int2_z+ cutoffz2
  END IF

END SUBROUTINE cutoff_corrections_units





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DFT_RAD_INT
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

SUBROUTINE DFT_RAD_INT_units (i,j,ih,zz_1,rhop,f_int_r,my_int_r,  &
     f_int2_r,my_int2_r,my_int3_r)

  USE DFT_MODULE
  USE EOS_VARIABLES, ONLY: parame
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: j
  INTEGER, INTENT(IN OUT)                :: ih
  REAL, INTENT(IN)                       :: zz_1
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: f_int_r
  REAL, INTENT(OUT)                      :: my_int_r
  REAL, INTENT(OUT)                      :: f_int2_r
  REAL, INTENT(OUT)                      :: my_int2_r
  REAL, INTENT(OUT)                      :: my_int3_r

  !-----------------------------------------------------------------------------
  INTEGER                                :: k

  REAL                                   :: zz1
  REAL                                   :: sig_ij
  REAL                                   :: dzr_local
  REAL                                   :: fint0, fint0_2, fint1, fint1_2
  REAL                                   :: myint0, myint0_2, myint0_3, myint1
  REAL                                   :: myint1_2, myint1_3
  REAL                                   :: dg_drho, dg_dr
  REAL                                   :: rad, xg, rdf, rho_bar, ua, rs
  REAL                                   :: analytic1, analytic2, tau_rs
  LOGICAL                                :: shortcut
  !-----------------------------------------------------------------------------

  shortcut = .true.
  fint0    = rc * tau_cut                 ! first order term
  fint0_2  = rc * tau_cut*tau_cut         ! 2nd order
  myint0   = rc * tau_cut                 ! first order term
  myint0_2 = rc * tau_cut*tau_cut         ! 2nd order
  myint0_3 = 0.0                          ! 2nd order

  f_int_r  = 0.0
  f_int2_r = 0.0
  my_int_r = 0.0
  my_int2_r= 0.0
  my_int3_r= 0.0

  !-----------------------------------------------------------------------------
  ! for mixtures it is advantageous to write all distances in dimensionless
  ! form, e.g. r^hat = r^hat / sigma.
  ! For mixtures, sig_ij = ( sig_i + sig_j ) / 2.
  !-----------------------------------------------------------------------------
  sig_ij = parame(1,2)                   ! for pure substances
  zz1 = zz_1 / sig_ij

  ! --- this block only speeds up the integration -------------------------------
  IF ( shortcut ) THEN
     rs = MAX( rg, zz1 ) ! +dzr
     IF ( rs > rc ) WRITE (*,*) 'error !!!!'
     analytic1 = 0.4*rs**-10 - 0.4*rc**-10 - rs**-4  + rc**-4
     analytic2 = 16.0/22.0 * (rs**-22 - rc**-22 ) - 2.0*rs**-16 +2.0*rc**-16  +1.6*rs**-10 - 1.6*rc**-10
     f_int_r  = f_int_r  + analytic1
     f_int2_r = f_int2_r + analytic2
     my_int_r = my_int_r + analytic1
     my_int2_r= my_int2_r+ analytic2
     IF ( rs == zz1 ) GO TO 10
     tau_rs  = 4.0 * ( rs**-12 - rs**-6 )
     fint0   = rs * tau_rs
     fint0_2 = rs * tau_rs*tau_rs
     myint0  = rs * tau_rs
     myint0_2= rs * tau_rs*tau_rs
     rad = rs                      ! the simple integration scheme: set to rc
     k = 0 + NINT( (rc-rs)/dzr )   ! in simple scheme: set to 0
  ELSE
     rad = rc
     k = 0
  END IF
  !-------------
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
     xg = rad / dhs_st
     rho_bar = ( rhop(i) + rhop(j) )/2.0       * sig_ij**3
     rdf = 1.0
     dg_drho = 0.0
     IF ( rad <= rg ) THEN
        CALL BI_CUB_SPLINE (rho_bar,xg,ya,x1a,x2a,y1a,y2a,y12a,  &
             c_bicub,rdf,dg_drho,dg_dr,den_step,ih,k)
        dg_drho = dg_drho * parame(1,2)**3        ! caution: introduced with length dimensions in Angstrom
     END IF

     fint1    = rdf * rad * ua
     fint1_2  = rdf * rad * ua * ua
     myint1   = rad * (rdf + 0.5*rhop(i)*dg_drho) * ua
     myint1_2 = rdf * rad * ua * ua
     myint1_3 = dg_drho * rad * ua * ua
     ! intf(k) = fint1
     f_int_r  = f_int_r  + dzr_local * (fint1 + fint0)/2.0
     f_int2_r = f_int2_r + dzr_local * (fint1_2 + fint0_2)/2.0
     my_int_r = my_int_r + dzr_local * (myint1 + myint0)/2.0
     my_int2_r= my_int2_r+ dzr_local * (myint1_2 + myint0_2)/2.0
     my_int3_r= my_int3_r+ dzr_local * (myint1_3 + myint0_3)/2.0

     fint0   = fint1
     fint0_2 = fint1_2
     myint0  = myint1
     myint0_2= myint1_2
     myint0_3= myint1_3

     !IF ( zz1 >= 1.0 .AND. rad-dzr_local+1.E-8 >= zz1 ) integration_cycle = .true.  ! integration down to ABS(zz1)
     !IF ( zz1 < 1.0  .AND. rad-dzr_local+1.E-8 >= 1.0 ) integration_cycle = .true.  ! integration down to ABS(zz1) but for r^hat<1, g(r)=0 (so stop at r^hat=1)
     !integration_cycle = .false.

  ENDDO

  ! IF (k.GT.30) THEN
  !   stepno = k
  !   CALL SPLINE_PARA (dzr,intf,utri,stepno)
  !   CALL SPLINE_INT (f_int_r,dzr,intf,utri,stepno)
  ! ENDIF

10 CONTINUE

  analytic1 = 4.0/10.0*rc**-10  - rc**-4
  analytic2 = 16.0/22.0*rc**-22 - 2.0*rc**-16 + 1.6*rc**-10
  f_int_r  = f_int_r  + analytic1
  f_int2_r = f_int2_r + analytic2
  my_int_r = my_int_r + analytic1
  my_int2_r= my_int2_r+ analytic2

  f_int_r  = f_int_r  * sig_ij*sig_ij
  f_int2_r = f_int2_r * sig_ij*sig_ij
  my_int_r = my_int_r * sig_ij*sig_ij
  my_int2_r= my_int2_r* sig_ij*sig_ij

END SUBROUTINE DFT_RAD_INT_units


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE rdf_matrix_units (msegm)

  USE PARAMETERS, ONLY: PI
  USE DFT_MODULE
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                        :: msegm

  !-----------------------------------------------------------------------------
  INTEGER                                 :: i, k = 0
  REAL                                    :: rdf, rad, xg, rho_rdf, z3
  !-----------------------------------------------------------------------------


  do i = 1, den_step

     ! rho_rdf= rhob(1) + (rhob(2)-rhob(1))*REAL(i)/den_step
     rho_rdf = 1.E-5 + (1.0) * REAL(i-1) / REAL(den_step-1)  ! segment density, rho_s*sigma**3
     rho_rdf = rho_rdf / msegm                               ! molar density, rho_m*sigma**3
     rad = rc
     k = 0
     !write (*,*) 'eta=',rho_rdf*msegm * PI/6.0  * dhs_st**3

     do while ( rad - 1.E-8 > 0.95 )

        rad = rad - dzr
        k = k + 1
        xg = rad / dhs_st
        z3 = rho_rdf * msegm * PI/6.0  * dhs_st**3      ! dhs_st is dim.less effective diam. d*(T)=d(T)/sigma
        ya(i,k) = 1.0
        IF ( xg <= rg .AND. z3 > 0.0 ) CALL rdf_int ( z3, msegm, xg, rdf )
        IF ( xg <= rg .AND. z3 > 0.0 ) ya(i,k) = rdf
        ! ya(i,k) = y(x1a(i), x2a(k))  with x1a: density-vector, x2a: r-vector
        x1a(i) = rho_rdf
        x2a(k) = xg

     end do

  end do

  if ( xg > 1.0 ) stop 'rdf_matrix_units: 0.95*sigma is too high for lower bound'

  WRITE (*,*) ' done with calculating g(r)',dhs_st

  kmax = k
  CALL bicub_derivative ( ya, x1a, x2a, y1a, y2a, y12a, den_step, kmax )
  CALL bicub_c ( ya, x1a, x2a, y1a, y2a, y12a, c_bicub, den_step, kmax )

END SUBROUTINE rdf_matrix_units
