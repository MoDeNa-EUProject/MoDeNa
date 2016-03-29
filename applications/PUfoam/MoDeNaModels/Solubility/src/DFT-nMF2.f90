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

SUBROUTINE DFT_nMFT2

  USE parameters, ONLY: PI, RGAS, KBOL
  USE basic_variables
  USE EOS_VARIABLES, ONLY: fres, eta, rho, x, z3t
  USE DFT_MODULE
  use utilities
  use EOS_POLAR, only: f_polar, phi_polar
  USE EOS_NUMERICAL_DERIVATIVES

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     REAL FUNCTION DENS_INVERS2 ( rhob, mu_rho )
       IMPLICIT NONE
       REAL, INTENT(IN)                   :: rhob(2)
       REAL, INTENT(IN OUT)               :: mu_rho
     END FUNCTION DENS_INVERS2

     SUBROUTINE DENS_INV_COEFF2 ( rhob )
       IMPLICIT NONE
       REAL, INTENT(IN)                   :: rhob(2)
     END SUBROUTINE DENS_INV_COEFF2
  END INTERFACE

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
  REAL                                   :: dev, maxdev
  REAL                                   :: delta_f, free
  REAL                                   :: surftens(0:200), st_macro(200)
  REAL                                   :: f01, f02, f03, f04, f05
  REAL                                   :: c1_con, c2_con
  REAL                                   :: zms
  REAL                                   :: tsav, psav, tc, pc, rhoc
  REAL                                   :: density(np), w(np,nc), lnphi(np,nc)
  REAL                                   :: damppic
  REAL                                   :: box_l_no_unit
  REAL, DIMENSION(-NDFT:NDFT)            :: lambda, rhobar, phi_dn0, phi_dn1
  REAL, DIMENSION(-NDFT:NDFT)            :: phi_dn2, phi_dn3, phi_dn4, phi_dn5
  REAL, DIMENSION(0:5,-NDFT:NDFT)        :: ni
  REAL                                   :: f_ch, dF_drho_ch
  REAL                                   :: dlngijdr, gij
  REAL                                   :: z0, z1, z2, z3
  REAL                                   :: f_fmt, dF_drho_fmt
  REAL                                   :: pN_pT(0:NDFT), zs
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


  z3t_st = PI/6.0* parame(1,1) * ( 1.0  &          ! divided by parame(i,2)**3
       * ( 1.0 - 0.12*EXP( -3.0*parame(1,3)/t ) ) )**3
  d_hs = ( 1.0 - 0.12*EXP( -3.0*parame(1,3)/t ) )  ! divided by parame(i,2)
  dhs_st = d_hs

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
  ! The the cut-off distance rc is the max. distance for integrating inter-
  ! actions. Beyond rc, tail-corrections are added.
  !-----------------------------------------------------------------------------
  rc = 9.0


  !-----------------------------------------------------------------------------
  ! basic definitions for calculating density profile,
  !-----------------------------------------------------------------------------
  ! grid-size  dzp = zp(1)-zp(0)

  box_l_no_unit = 50.0          ! lenth of simulation box (dimensionless, L/sigma)
  discret= 1000                 ! number of spacial discretizations for profile
  dzp    = box_l_no_unit / REAL(discret)    ! dimensionless grid-distance dzp*=dzp/sigma
  fa     = NINT( 1.0 / dzp )    ! size of the box is discret/fa * sigma   (sigma: LJ diameter)
  WRITE (*,*) fa, 1.0/dzp
  outpt  = 100                  ! number output-points = discret/outpt


  !-----------------------------------------------------------------------------
  ! definitions for the numerical algorithm
  !-----------------------------------------------------------------------------

  maxiter = 180                  ! maximum number of iterations per temp.step
  maxdev  = 2.E-7               ! maximum deviation
  damppic = 0.002               ! damping of Picard-iteration


  !-----------------------------------------------------------------------------
  ! get a matrix of values for the pair correlation function g(r,rho)
  !-----------------------------------------------------------------------------

  rg = 4.0
  den_step = 40
  dzr = dzp / 2.0               ! dimensionless grid-distance for integrating
  ! the Barker-Henderson attraction term
  CALL rdf_matrix (msegm)


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


     !--------------------------------------------------------------------------
     ! get density in SI-units (kg/m**3)
     !--------------------------------------------------------------------------
     CALL SI_DENS ( density, w )


     !--------------------------------------------------------------------------
     ! (re-)calculate residual chemical potential of both phases
     !--------------------------------------------------------------------------
     ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)
     densta(1) = dense(1)                    ! Index 1 is for liquid density (here: packing fraction eta)
     densta(2) = dense(2)                    ! Index 2 is for vapour density (here: packing fraction eta)
     CALL fugacity (lnphi)
     my0 = lnphi(1,1) + LOG( rhob(1)/parame(1,2)**3 )     ! my0 = mu_res(T,rho_bulk_L) + ln(rho_bulk_l)
     zges(1) = p / ( RGAS*t*density(1) ) * mm(1)/1000.0
     zges(2) = p / ( RGAS*t*density(2) ) * mm(1)/1000.0

     pbulk = zges(1) * rhob(1)               ! pressure  p/kT (= Z*rho)
     WRITE (*,*) 'chem. potentials', lnphi(1,1) + LOG( rhob(1)/parame(1,2)**3 ),  &
          lnphi(2,1) + LOG( rhob(2)/parame(1,2)**3 )
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
     ! define initial density profile rhop(i)
     ! and dimensionless space coordinates zp(i)
     !
     ! discret  : number of grid-points within "the box"
     ! irc      : number of grid-points extending the the box to left and right
     !            the box is extended in order to allow for numerical integration
     !            'irc' is determined by the cut-off distance 'rc'
     !--------------------------------------------------------------------------

     irc = NINT(rc/dzp) + 1
     tanhfac = -2.3625*t/tc + 2.4728
     ! tanhfac = 2.7*(1.0 - t/tc)   ! this parameterization was published (Gross, JCP, 2009)
     DO i = -irc, (discret+irc)
        zp(i) = REAL(i) * dzp
     END DO
     DO i = -irc, (discret+irc)
        rhop(i) = TANH(-(zp(i)-zp(INT(discret/2)))*tanhfac) * (rhob(1)-rhob(2))/2.0  &
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

           ! ---------------------------------------------------------------
           ! some auxiliary quantities used for compensating LDA contributions
           ! ---------------------------------------------------------------

           z0 = rhop(i) * PI/6.0 * parame(1,1)
           z1 = rhop(i) * PI/6.0 * parame(1,1) * d_hs
           z2 = rhop(i) * PI/6.0 * parame(1,1) * d_hs**2
           z3 = rhop(i) * PI/6.0 * parame(1,1) * d_hs**3
           zms = 1.0 - z3
           gij = 1.0/zms + 3.0*(d_hs/2.0)*z2/zms/zms + 2.0*((d_hs/2.0)*z2)**2 /zms**3
           dlngijdr = (1.0/zms/zms + 3.0*(d_hs/2.0)*z2*(1.0+z3)/z3/zms**3   &
                + ((d_hs/2.0)*z2/zms/zms)**2 *(4.0+2.0*z3)/z3 ) /gij*z3t_st

           ! ---------------------------------------------------------------
           ! hard sphere contribution
           !    f_fmt  is the Helmholtz energy density F_FMT/(VkT) = F_FMT/(NkT)*rho
           !    dF_drho_fmt is the functional derivative of F_FMT/(kT) to rho
           ! ---------------------------------------------------------------
           zs = 1.0 - ni(3,i)
           f_fmt  = - ni(0,i)*LOG(zs) + ni(1,i)*ni(2,i)/zs - ni(4,i)*ni(5,i)/zs  &
                + (ni(2,i)**3 -3.0*ni(2,i)*ni(5,i)*ni(5,i)) *(ni(3,i)+zs*zs*LOG(zs))  &
                /36.0/PI/zs/zs/ni(3,i)**2

           CALL dF_FMT_dr ( i, fa, dzp, d_hs, zp, phi_dn0,  &
                phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5, dF_drho_fmt )


           ! ---------------------------------------------------------------
           ! chain term
           ! ---------------------------------------------------------------
           CALL dF_chain_dr ( i, fa, dzp, d_hs, zp, rhop, lambda,  &
                rhobar, z3t_st, f_ch, dF_drho_ch )


           ! ---------------------------------------------------------------
           ! dispersive attraction
           ! ---------------------------------------------------------------
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
              IF ( zz1 <= rc-dzp+1.E-9 ) THEN
                 CALL dft_rad_int ( i, j, ih, zz1, rhop, f_int_r, mu_int_r,  &
                      f_int2_r, mu_int2_r, mu_int3_r )
              ELSE
                 f_int_r  = 0.0
                 f_int2_r = 0.0
                 mu_int_r = 0.0
                 mu_int2_r= 0.0
                 mu_int3_r= 0.0
              END IF

              f_int_z  = f_int_z  + dzp * (rhop(j)* f_int_r + f01)/2.0
              mu_int_z = mu_int_z + dzp * (rhop(j)*mu_int_r + f02)/2.0
              f_int2_z = f_int2_z + dzp * ( f_int2_r + f03)/2.0
              mu_int2_z= mu_int2_z+ dzp * (mu_int2_r + f04)/2.0
              mu_int3_z= mu_int3_z+ dzp * (mu_int3_r + f05)/2.0
              f01 = rhop(j)* f_int_r
              f02 = rhop(j)* mu_int_r
              f03 = f_int2_r
              f04 = mu_int2_r
              f05 = mu_int3_r
           END DO

           ! ---------------------------------------------------------------
           ! cut-off corrections
           ! ---------------------------------------------------------------

           CALL cutoff_corrections ( i, irc, rhop, rhob, f_int_z, mu_int_z, f_int2_z, mu_int2_z )


           ! ---------------------------------------------------------------
           ! for second order dispersive attraction ( usually not used !!! )
           ! ---------------------------------------------------------------

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

           ! ---------------------------------------------------------------
           ! The integration of the DFT-integral can be done with cubic splines
           ! ---------------------------------------------------------------
           ! stepno = discret + irc*2
           ! CALL SPLINE_PARA (dzp, intgrid, utri, stepno)
           ! CALL SPLINE_INT (f_int_z, dzp, intgrid, utri, stepno)
           ! f_int_z = f_int_z + rhob(1)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )
           ! f_int_z = f_int_z + rhob(2)*( 4.0/90.0*rc**-9 -1.0/3.0*rc**-3 )


           ! ---------------------------------------------------------------
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
           ! ---------------------------------------------------------------
           ! remember: factor 2*PI because of the cylindrical coordinates

           dF_drho_att(i) = 2.0*PI*parame(1,1)**2 * parame(1,3)/t * mu_int_z
           f_att(i) = PI*parame(1,1)**2 *parame(1,3)/t*rhop(i)*f_int_z


           ! --- if a second order perturbation theory is considered -------
           IF (subtract1 =='2PT') THEN
              dF_drho_att(i)= dF_drho_att(i) -2.0*PI*parame(1,1)**3 *(parame(1,3)/t)**2 *(mu_int2_z+mu_int3_z)
              f_att(i) = f_att(i)  -2.0*PI*parame(1,1)**3 *(parame(1,3)/t)**2 *f_int2_z
           END IF

           ! ---------------------------------------------------------------
           ! make the dispersive attraction term consistent with PC-SAFT, by
           ! adding the difference ( PC-SAFT - 1PT )_dispersion locally (LDA)
           ! ---------------------------------------------------------------

           !*****************************************************************************
           ! under construction
           !*****************************************************************************
           dense(1) = rhop(i) * z3t_st
           densta(1) = dense(1)

           ensemble_flag = 'tv'

           call ONLY_ONE_TERM_EOS_NUMERICAL ( 'disp_term', 'PT1      ' )
           call FUGACITY ( lnphi )
           call RESTORE_PREVIOUS_EOS_NUMERICAL
           !call f_PT (..... f_disp_1PT)
           !call mu_PT ( .... mu_disp_1PT)
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


           ! ---------------------------------------------------------------
           ! apart from the perturbation theory, dF_drho_att and f_att also
           ! collect the non-local parts of the hs and chain fluid
           ! ---------------------------------------------------------------

           dF_drho_tot(i) = LOG( rhop(i)/parame(1,2)**3 ) + dF_drho_fmt + dF_drho_ch + dF_drho_att(i)
           f_tot(i) = f_fmt + f_ch + f_att(i)


           ! --- on-the-fly report on local errors in the profile ----------
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
           dev = dev + (rhop(i)-rhop_o(i))**2
        END DO


        !-----------------------------------------------------------------------
        ! calculate surface tension
        !-----------------------------------------------------------------------
        free = 0.0
        DO i = 1, discret
           delta_f = f_tot(i) + rhop(i)*(LOG(rhop(i)/parame(1,2)**3 )-1.0)   -   (rhop(i)*my0 - pbulk)  ! all quantities .../(kT)
           pN_pT(i) = delta_f
           free  = free  + delta_f*dzp
        END DO
        surftens(kk) = KBOL * t *1.E20*1000.0 *free / parame(1,2)**2


        !-----------------------------------------------------------------------
        ! add an approximate capillary wave contribution to the surface tension
        !-----------------------------------------------------------------------
        st_macro(kk) = surftens(kk) / ( 1.0 + 3.0/8.0/PI *t/tc  &
             * (1.0/2.55)**2  / (0.0674*parame(1,1)+0.0045) )


        delta_st = 1.0
        IF ( kk > 1 ) delta_st = ABS( surftens(kk)-surftens(kk-1) ) / surftens(kk)

        WRITE (*,*)  '-----------------------------------------------------------'
        WRITE (*,*)  ' #      error-sum      intrinsic ST    total ST'
        WRITE (*,'(i4,3F16.8)')  kk, dev, surftens(kk), st_macro(kk)
        WRITE (*,*)  '-----------------------------------------------------------'
        WRITE (71,'(i4,4G18.8)') kk, dev, surftens(kk), st_macro(kk)
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
     WRITE(68,'(i2,6f18.8)') 0, zp(0)-zp(INT(discret/2)), rhop(0), t,  &
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
           z3t_st = PI/6.0 * parame(1,1) * ( 1.0  &   ! divided by parame(i,2)**3
                *( 1.0 - 0.12*EXP(-3.0*parame(1,3)/t) ) )**3
           d_hs = ( 1.0 - 0.12*EXP(-3.0*parame(1,3)/t) )  ! divided by parame(i,2)
           dhs_st = d_hs
           CALL rdf_matrix (msegm)
        ELSE
           EXIT Diagram_for_various_T_Loop
        END IF
        WRITE (69,'(7(f18.10))') tc, pc/1.E5, rhoc, rhoc, 0., 0., 0.
     ELSE
        EXIT Diagram_for_various_T_Loop
     END IF

  ENDDO Diagram_for_various_T_Loop

END SUBROUTINE DFT_nMFT2




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE cutoff_corrections ( i, irc, rhop, rhob, f_int_z, mu_int_z, f_int2_z, mu_int2_z )

  USE DFT_Module, ONLY: rc, NDFT
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

  cutoffz = 4.0/90.0 * rc**-9 - 1.0/3.0 * rc**-3
  cutoffz2 = 16.0/22.0/21.0 * rc**-21 - 2.0/15.0 * rc**-15 + 16.0/90.0 * rc**-9
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

END SUBROUTINE cutoff_corrections


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE aux_chain (discret,fa,dzp,d_hs,zp,rhop,rhobar,lambda)

  USE basic_variables
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: discret
  INTEGER, INTENT(IN)                    :: fa
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: d_hs
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: rhobar(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: lambda(-NDFT:NDFT)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, nn
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

  DO i = (-fa), (discret+fa)
     nn = 1
     DO j = (i-fa), (i+fa)      ! fa=NINT(1./dzp)
        IF ( zp(j+1) > (zp(i)-d_hs).AND.(zp(i)-d_hs) >= zp(j) ) THEN
           dz =  zp(j+1) - (zp(i)-d_hs)
           rx(1) = 0.0
           ry1(1) = 0.0  !  = 0.75*rhop(i-fa)*( d.**2 -d.**2 )
           ry2(1) = 0.5/d_hs*rhop(j) + ((zp(i)-d_hs)-zp(j))/dzp  &
                *( 0.5/d_hs*rhop(j+1) - 0.5/d_hs*rhop(j) )
        ELSE IF ( zp(j) > (zp(i)-d_hs) .AND. zp(j) <= (zp(i)+d_hs) ) THEN
           zz1 = ABS( zp(j)-zp(i) )          ! distance z12 between 1 and 2
           nn = nn + 1
           ry1(nn) = 0.75/d_hs**3 * rhop(j) * (d_hs**2 -zz1*zz1)
           ry2(nn) = 0.5/d_hs * rhop(j)
           rx(nn)  = rx(nn-1) + dz
           dz = dzp
           IF ( (zp(i)+d_hs) < zp(j+1) ) THEN
              dz = (zp(i)+d_hs) - zp(j)
              nn = nn + 1
              rx(nn) = rx(nn-1) + dz
              ry1(nn) = 0.0
              ry2(nn) = 0.5/d_hs*rhop(j) + ((zp(i)+d_hs)-zp(j))/dzp  &
                   *( 0.5/d_hs*rhop(j+1) - 0.5/d_hs*rhop(j) )
           END IF
        END IF
     END DO
     xl = rx(1)
     xh = rx(nn)
     CALL spline          ( rx, ry1, nn, 1.E30, 1.E30, y2 )
     CALL splint_integral ( rx, ry1, y2, nn, xl, xh, int1 )
     rhobar(i) = int1
     CALL spline          ( rx, ry2, nn, 1.E30, 1.E30, y2 )
     CALL splint_integral ( rx, ry2, y2, nn, xl, xh, int2 )
     lambda(i) = int2
  END DO

END SUBROUTINE aux_chain



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE fmt_dens ( discret, fa, dzp, d_hs, zp, rhop, ni,  &
     phi_dn0, phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5 )

  USE parameters, ONLY: PI
  USE basic_variables
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN OUT)                :: discret
  INTEGER, INTENT(IN)                    :: fa
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: d_hs
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: ni(0:5,-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn0(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn1(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn2(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn3(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn4(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: phi_dn5(-NDFT:NDFT)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, nn
  REAL                                   :: int2, int3, int5
  REAL                                   :: zz1, zs, d2, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy2, hy3, hy5
  !-----------------------------------------------------------------------------

  DO i = (-fa), (discret+fa)
     nn = 1
     d2 = d_hs / 2.0
     DO j = (i-fa/2), (i+fa/2)
        IF ( zp(j+1) > (zp(i)-d2) .AND. (zp(i)-d2) >= zp(j) ) THEN
           dz =  zp(j+1) - (zp(i)-d2)
           hx(1) = 0.0
           hy2(1) = rhop(j)  *parame(1,1) + ((zp(i)-d2)-zp(j))/dzp  &
                * ( rhop(j+1)*parame(1,1) - rhop(j)*parame(1,1) )
           hy3(1) = 0.0            ! = rhop(i-fa/2)*parame(1,1)*(0.25-(zp(i-fa/2)-zp(i))**2 )
           hy5(1) = rhop(j)*parame(1,1)*(zp(j)-zp(i)) +((zp(i)-d2)-zp(j))/dzp  &
                * ( rhop(j+1)*parame(1,1)*(zp(j+1)-zp(i))  &
                - rhop(j)  *parame(1,1)*(zp(j)-zp(i)) )
        ELSE IF ( zp(j) > (zp(i)-d2) .AND. zp(j) <= (zp(i)+d2) ) THEN
           zz1 = zp(j)-zp(i)       ! distance z12 between 1 and 2
           nn = nn + 1
           hy2(nn) = rhop(j)*parame(1,1)
           hy3(nn) = rhop(j)*parame(1,1) * ( d2**2 - zz1**2 )
           hy5(nn) = rhop(j)*parame(1,1) * zz1
           hx(nn)  = hx(nn-1) + dz
           dz = dzp
           IF ( (zp(i)+d2) < zp(j+1) ) THEN
              dz = (zp(i)+d2) - zp(j)
              nn = nn + 1
              hy2(nn) = rhop(j)*parame(1,1) + ((zp(i)+d2)-zp(j))/dzp  &
                   * ( rhop(j+1)*parame(1,1) - rhop(j)*parame(1,1) )
              hy3(nn) = 0.0
              hy5(nn) = rhop(j)*parame(1,1)*(zp(j)-zp(i)) +((zp(i)+d2)-zp(j))/dzp  &
                   * ( rhop(j+1)*parame(1,1)*(zp(j+1)-zp(i))  &
                   - rhop(j)*parame(1,1)*(zp(j)  -zp(i)) )
              hx(nn) = hx(nn-1) + dz
           END IF
        END IF
     END DO
     xl = hx(1)
     xh = hx(nn)
     CALL spline          ( hx, hy2, nn, 1.E30, 1.E30, y2 )
     CALL splint_integral ( hx, hy2, y2, nn, xl, xh, int2 )
     CALL spline          ( hx, hy3, nn, 1.E30, 1.E30, y2 )
     CALL splint_integral ( hx, hy3, y2, nn, xl, xh, int3 )
     CALL spline          ( hx, hy5, nn, 1.E30, 1.E30, y2 )
     CALL splint_integral ( hx, hy5, y2, nn, xl, xh, int5 )
     ni(2,i) = PI * d_hs *int2             ! corresponds to z2*6
     ni(3,i) = PI * int3                   ! corresponds to z3
     ni(5,i) = 2.0* PI *int5
     ni(0,i) = ni(2,i) / (PI*d_hs**2 )     ! corresponds to z0/PI*6.0
     ni(1,i) = ni(2,i) / (2.0*PI*d_hs)     ! corresponds to z1/PI*3.0
     ni(4,i) = ni(5,i) / (2.0*PI*d_hs)
  END DO

  !     derivatives of phi to ni
  DO i = (-fa), (discret+fa)
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

END SUBROUTINE fmt_dens


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE dF_FMT_dr ( i, fa, dzp, d_hs, zp, phi_dn0,  &
     phi_dn1, phi_dn2, phi_dn3, phi_dn4, phi_dn5, dF_drho_fmt )

  USE parameters, ONLY: PI
  USE basic_variables
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER, INTENT(IN)                    :: fa
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: d_hs
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn0(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn1(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn2(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn3(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn4(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: phi_dn5(-NDFT:NDFT)
  REAL, INTENT(OUT)                      :: dF_drho_fmt

  !-----------------------------------------------------------------------------
  INTEGER                                :: j, nn
  REAL                                   :: int0, int1, int2, int3, int4, int5
  REAL                                   :: zz1, d2, xl, xh, dz = 0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy0, hy1, hy2, hy3, hy4, hy5
  !-----------------------------------------------------------------------------


  nn = 1
  d2 = d_hs / 2.0
  DO j = (i-fa/2), (i+fa/2)
     IF ( zp(j+1) > (zp(i)-d2).AND.(zp(i)-d2) >= zp(j) ) THEN
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
        hy3(nn) = phi_dn3(j) * ( d2**2 -zz1**2 )
        hy4(nn) = phi_dn4(j) * zz1
        hy5(nn) = phi_dn5(j) * zz1
        hx(nn)  = hx(nn-1) + dz
        dz = dzp
        IF ( (zp(i)+d2) < zp(j+1) ) THEN
           dz = (zp(i)+d2) - zp(j)
           nn = nn + 1
           hx(nn)  = hx(nn-1) + dz
           hy0(nn) = phi_dn0(j) + ( (zp(i)+d2) - zp(j) ) / dzp  &
                * ( phi_dn0(j+1) - phi_dn0(j) )
           hy1(nn) = phi_dn1(j) + ( (zp(i)+d2) - zp(j) ) / dzp  &
                * ( phi_dn1(j+1) - phi_dn1(j) )
           hy2(nn) = phi_dn2(j) + ( (zp(i)+d2) - zp(j) ) / dzp  &
                * ( phi_dn2(j+1) - phi_dn2(j) )
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
  CALL spline         ( hx, hy0, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy0, y2, nn, xl, xh, int0 )
  CALL spline         ( hx, hy1, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy1, y2, nn, xl, xh, int1 )
  CALL spline         ( hx, hy2, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy2, y2, nn ,xl, xh, int2 )
  CALL spline         ( hx, hy3, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy3, y2, nn, xl, xh, int3 )
  CALL spline         ( hx, hy4, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy4, y2, nn, xl, xh, int4 )
  CALL spline         ( hx, hy5, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy5, y2, nn, xl, xh, int5 )
  ! write (*,*) int3,int4,int5
  ! call paus (' ')
  ! int0 = int0/d_hs
  ! int1 = 0.5*parame(1,2)*int1
  ! int2 = PI*d_hs*parame(1,2)**2 *int2
  ! int3 = PI*parame(1,2)**3 *int3
  ! int4 = parame(1,2)*int4
  ! int5 = 2.0*PI*parame(1,2)**2 *int5
  int0 = int0 / d_hs
  int1 = 0.5 * int1
  int2 = PI * d_hs * int2
  int3 = PI * int3
  int4 = int4 / d_hs
  int5 = 2.0 * PI * int5
  ! dF_drho_fmt=dF_hs/drho = dF_hs/drho_segment *m_mean
  dF_drho_fmt = (int0+int1+int2+int3+int4+int5)*parame(1,1)

END SUBROUTINE dF_FMT_dr


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE dF_chain_dr ( i, fa, dzp, d_hs, zp, rhop, lambda,  &
     rhobar, z3t_st, f_ch, dF_drho_ch )

  USE parameters, ONLY: PI
  USE basic_variables
  USE DFT_MODULE, ONLY: NDFT
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  INTEGER                                :: fa
  REAL, INTENT(IN)                       :: dzp
  REAL, INTENT(IN)                       :: d_hs
  REAL, INTENT(IN)                       :: zp(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhop(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: lambda(-NDFT:NDFT)
  REAL, INTENT(IN)                       :: rhobar(-NDFT:NDFT)
  REAL, INTENT(IN OUT)                   :: z3t_st
  REAL, INTENT(OUT)                      :: f_ch
  REAL, INTENT(OUT)                      :: dF_drho_ch

  !-----------------------------------------------------------------------------
  INTEGER                                :: j, nn
  REAL                                   :: zz1, xl, xh, dz=0.0
  REAL, DIMENSION(100)                   :: y2, hx, hy0, hy1
  REAL                                   :: ycorr, dlnydr, I_ch1, I_ch2
  !-----------------------------------------------------------------------------


  nn = 1
  DO j = (i-fa), (i+fa)         ! fa=NINT(1./dzp)
     IF ( zp(j+1) > (zp(i)-d_hs) .AND. zp(j)-1.E-11 <= (zp(i)-d_hs) ) THEN
        dz =  zp(j+1) - (zp(i)-d_hs)
        hx(1)  = 0.0
        hy0(1) = 0.0  !  = 0.75*rhop(i-fa)*DLNYDR(d_hs,rhobar(i-fa))*(d.**2 -d.**2 )
        hy1(1) = 0.5/d_hs*rhop(j)/lambda(j) + ((zp(i)-d_hs)-zp(j))/dzp  &
             * ( 0.5/d_hs*rhop(j+1)/lambda(j+1) - 0.5/d_hs*rhop(j)/lambda(j) )
     ELSE IF ( zp(j) > (zp(i)-d_hs).AND.zp(j)+1.E-11 < (zp(i)+d_hs) ) THEN
        zz1 = zp(j) - zp(i)       ! distance z12 between 1 and 2
        nn = nn + 1
        CALL cavity ( d_hs, z3t_st, rhobar(j), ycorr, dlnydr )
        hy0(nn) = 0.75/d_hs**3  *rhop(j)*dlnydr*(d_hs*d_hs-zz1*zz1)
        hy1(nn) = 0.5/d_hs*rhop(j)/lambda(j)
        hx(nn)  = hx(nn-1) + dz
        dz = dzp
        IF ( zp(j+1)+1.E-11 >= (zp(i)+d_hs) ) THEN
           dz = (zp(i)+d_hs) - zp(j)
           nn = nn + 1
           hy0(nn) = 0.0
           hy1(nn) = 0.5/d_hs*rhop(j)/lambda(j) +((zp(i)+d_hs)-zp(j))/dzp  &
                * ( 0.5/d_hs*rhop(j+1)/lambda(j+1)-0.5/d_hs*rhop(j)/lambda(j) )
           hx(nn) = hx(nn-1) + dz
        END IF
     END IF
  END DO
  xl = hx(1)
  xh = hx(nn)
  CALL spline         ( hx, hy0, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy0, y2, nn, xl, xh, I_ch1 )
  CALL spline         ( hx, hy1, nn, 1.E30, 1.E30, y2 )
  CALL splint_integral( hx, hy1, y2, nn, xl, xh, I_ch2 )
  ! write (*,*) '111',I_ch1,I_ch2

  ! DO j=1,nn
  !   write (*,*) j,hx(j),hy1(j)
  ! ENDDO

  ! I_ch1 = 0.0
  ! I_ch2 = 0.0
  ! DO j=1,nn-1
  !   I_ch1=I_ch1+(hy0(j)+hy0(j+1))/2.0*(hx(j+1)-hx(j))
  !   I_ch2=I_ch2+(hy1(j)+hy1(j+1))/2.0*(hx(j+1)-hx(j))
  ! ENDDO
  ! write (*,*) '222',I_ch1,I_ch2

  CALL cavity ( d_hs, z3t_st, rhobar(i), ycorr, dlnydr )
  f_ch = rhop(i) * ( LOG(ycorr*lambda(i)) -1.0 )
  f_ch = - ( parame(1,1) - 1.0 ) * f_ch
  f_ch = f_ch + ( parame(1,1) - 1.0 ) * rhop(i) * ( LOG(rhop(i)) - 1.0 )
  dF_drho_ch= LOG( ycorr*lambda(i) ) - 1.0 + I_ch1 + I_ch2
  dF_drho_ch= - ( parame(1,1) - 1.0 ) * dF_drho_ch
  dF_drho_ch= dF_drho_ch + ( parame(1,1) - 1.0 ) * LOG(rhop(i))

END SUBROUTINE dF_chain_dr


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE cavity ( d_hs, z3t_st, rho, ycorr, dlnydr )

  USE parameters, ONLY: PI
  USE basic_variables
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: d_hs
  REAL, INTENT(IN)                       :: z3t_st
  REAL, INTENT(IN)                       :: rho
  REAL, INTENT(OUT)                      :: ycorr
  REAL, INTENT(OUT)                      :: dlnydr

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j
  REAL                                   :: z0, z1, z2, z3, zms
  REAL, DIMENSION(nc,nc)                 :: dij_ab, gij, dgijdz
  !-----------------------------------------------------------------------------

  ncomp = 1

  !z0 = PI / 6.0 * rho * SUM( x(1:ncomp) * mseg(1:ncomp) )
  !z1 = PI / 6.0 * rho * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp) )
  !z2 = PI / 6.0 * rho * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**2 )
  !z3 = PI / 6.0 * rho * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

  z0 = PI / 6.0 * rho * parame(1,1)
  z1 = PI / 6.0 * rho * parame(1,1) * d_hs
  z2 = PI / 6.0 * rho * parame(1,1) * d_hs**2
  z3 = PI / 6.0 * rho * parame(1,1) * d_hs**3
  zms    = 1.0 - z3

  ! DO i = 1,ncomp
  !    DO j=1,ncomp
  !       dij_ab(i,j)=dhs(i)*dhs(j)/(dhs(i)+dhs(j))
  !    ENDDO
  ! END DO
  dij_ab(1,1) = d_hs / 2.0

  DO i = 1, ncomp
     DO j = 1, ncomp
        gij(i,j) = 1.0/zms + 3.0*dij_ab(i,j)*z2/zms/zms  &
             + 2.0*(dij_ab(i,j)*z2)**2 /zms**3
        dgijdz(i,j)= 1.0/zms/zms + 3.0*dij_ab(i,j)*z2*(1.0+z3)/z3/zms**3   &
             + (dij_ab(i,j)*z2/zms/zms)**2 *(4.0+2.0*z3)/z3
     END DO
  END DO

  ycorr  = gij(1,1)
  dlnydr = dgijdz(1,1) / gij(1,1) * z3t_st

END SUBROUTINE cavity



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! Module DFT_INVERSION
!
! This module contains variables for a LDA - Inversion Procedure (an algo-
! rithm for solving for the equilibrium density profile). Procedure is
! similar to R.Evans (Bristol).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module DFT_INVERSION

  implicit none
  save

  INTEGER                               :: no_step_l, no_step_h
  REAL                                  :: mu_eta_div, den_min, den_max, den_div
  REAL, DIMENSION(200)                  :: rho_array1, rho_array2
  REAL, DIMENSION(200)                  :: mu_rho_array1, mu_rho_array2

End Module DFT_INVERSION



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DENS_INV_COEFF2
!
! First step of LDA - Inversion Procedure (algorithm for solving for the
! equilibrium density profile): Determine coefficients.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE DENS_INV_COEFF2 ( rhob )

  USE basic_variables, ONLY: np, nc, dense, densta, parame
  USE DFT_MODULE
  USE DFT_INVERSION
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rhob(2)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: den_st
  REAL                                   :: lnphi(np,nc)
  !-----------------------------------------------------------------------------

  den_min = rhob(2)*z3t_st * 0.8
  den_max = rhob(1)*z3t_st * 1.1
  den_div = 0.04

  no_step_l = 80
  no_step_h = 200

  mu_eta_div = 0.0
  IF ( den_min < den_div ) THEN
     DO i = 1, no_step_l
        dense(1)  = den_min + (den_div-den_min) * REAL(i-1)/REAL(no_step_l-1)
        densta(1) = dense(1)
        rho_array1(i) = dense(1)
        CALL fugacity (lnphi)
        mu_rho_array1(i) = EXP( lnphi(1,1)+LOG(dense(1)/z3t_st / parame(1,2)**3 ) )  ! myrho = mu_res(T,rho) + ln(rho)
     END DO
     mu_eta_div = mu_rho_array1(no_step_l)
  END IF

  den_st = MAX( den_min, den_div )
  DO i = 1, no_step_h
     dense(1)  = den_st + (den_max-den_st) * REAL(i-1)/REAL(no_step_h-1)
     densta(1) = dense(1)
     rho_array2(i) = dense(1)
     CALL fugacity (lnphi)
     mu_rho_array2(i) = lnphi(1,1) + LOG(dense(1)/z3t_st /parame(1,2)**3 )  ! myrho = mu_res(T,rho) + ln(rho)
     IF ( i >= 2 ) THEN
        IF ( (mu_rho_array2(i)- mu_rho_array2(i-1)) < 0.0 ) THEN
           call paus ('DENS_INV_COEFF2: derivative negative')
        END IF
     END IF
  END DO

END SUBROUTINE DENS_INV_COEFF2


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! FUNCTION DENS_INVERS2
!
! Second step of LDA - Inversion Procedure (algorithm for solving for the
! equilibrium density profile): Determine new guess of density profile.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION DENS_INVERS2 ( rhob, mu_rho )

  USE DFT_MODULE
  USE DFT_INVERSION
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL, INTENT(IN)                       :: rhob(2)
  REAL, INTENT(IN OUT)                   :: mu_rho

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: eta_i
  LOGICAL                                :: flag_exit
  !-----------------------------------------------------------------------------

  ! den_min = rhob(2)*z3t_st * 0.8
  ! den_max = rhob(1)*z3t_st * 1.1
  ! den_div = 0.04

  ! no_step_l = 40
  ! no_step_h = 100

  IF ( EXP(mu_rho) < mu_eta_div ) THEN
     i = 2
     flag_exit = .false.
     DO WHILE (i < no_step_l .AND. .NOT. flag_exit)
        i = i + 1
        IF ( mu_rho_array1(i) > EXP(mu_rho) ) THEN
           eta_i = rho_array1(i-1)+(rho_array1(i)-rho_array1(i-1))  &
                /(mu_rho_array1(i)-mu_rho_array1(i-1))  &
                *(EXP(mu_rho)-mu_rho_array1(i-1))
           DENS_INVERS2 = eta_i/z3t_st
           flag_exit = .true.
        END IF
     END DO
     IF ( .NOT. flag_exit ) WRITE (*,*) 'error 1', EXP(mu_rho), mu_eta_div
     ! IF ( .NOT. flag_exit ) call paus (' ')
  ENDIF

  IF ( EXP(mu_rho) >= mu_eta_div ) THEN
     i = 2
     flag_exit = .false.
     DO WHILE (i < no_step_h .AND. .NOT. flag_exit )
        i = i + 1
        IF ( mu_rho_array2(i) > mu_rho ) THEN
           eta_i = rho_array2(i-1)+(rho_array2(i)-rho_array2(i-1))  &
                /(mu_rho_array2(i)-mu_rho_array2(i-1)) *(mu_rho-mu_rho_array2(i-1))
           DENS_INVERS2 = eta_i/z3t_st
           flag_exit = .true.
        END IF
     END DO
     IF ( .NOT. flag_exit ) WRITE (*,*) 'error 2', mu_rho, mu_eta_div
     ! IF ( .NOT. flag_exit ) call paus (' ')
  END IF

  IF ( DENS_INVERS2 < (rhob(2)*0.8) ) THEN
     WRITE (*,*) 'lower bound', DENS_INVERS2, mu_rho
     DENS_INVERS2 = rhob(2) * 0.8
  END IF
  IF ( DENS_INVERS2 > (rhob(1)*1.1) ) THEN
     WRITE (*,*) 'upper bound', DENS_INVERS2, mu_rho
     DENS_INVERS2 = rhob(1)
  END IF

END FUNCTION DENS_INVERS2


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE spline
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function,
! i.e., yi = f(xi), with  x1 < x2 < .. . < xN, and given values yp1 and
! ypn for the first derivative of the interpolating function at points 1
! and n, respectively, this routine returns an array y2(1:n) of length n
! which contains the second derivatives of the interpolating function at
! the tabulated points xi. If yp1 and/or ypn are equal to 1  1030 or
! larger, the routine is signaled to set  the corresponding boundary
! condition for a natural spline, with zero second derivative on that
! boundary.
! Parameter: NMAX is the largest anticipated value of n.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE spline ( x, y, n, yp1, ypn, y2 )

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n
  REAL, INTENT(IN)                       :: x(n)
  REAL, INTENT(IN)                       :: y(n)
  REAL, INTENT(IN)                       :: yp1
  REAL, INTENT(IN)                       :: ypn
  REAL, INTENT(OUT)                      :: y2(n)

  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER                     :: NMAX = 500
  INTEGER                                :: i, k
  REAL                                   :: p, qn, sig, un, u(NMAX)
  !-----------------------------------------------------------------------------

  IF ( yp1 > 0.99E30 ) THEN
     y2(1) = 0.0
     u(1)  = 0.0
  ELSE
     y2(1) = -0.5
     u(1)  = ( 3.0/(x(2)-x(1)) ) * ( (y(2)-y(1))/(x(2)-x(1))-yp1 )
  END IF
  DO  i = 2, n-1
     IF ( (x(i+1)-x(i)) == 0.0 .OR. (x(i)-x(i-1)) == 0.0 .OR. (x(i+1)-x(i-1)) == 0.0 ) THEN
        write (*,*) 'error in spline-interpolation'
        stop
     END IF
     sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
     p = sig*y2(i-1) + 2.0
     y2(i) = (sig-1.0) / p
     u(i) = ( 6.0 * ((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1))  &
          - sig * u(i-1) ) / p
  END DO
  IF ( ypn > 0.99E30 ) THEN
     qn = 0.0
     un = 0.0
  ELSE
     qn = 0.5
     un = (3.0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  END IF
  y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.0)
  DO  k = n-1, 1, -1
     y2(k) = y2(k) * y2(k+1) + u(k)
  END DO

END SUBROUTINE spline


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE splint(xa,ya,y2a,n,x,y)
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the in order), and given the array y2a(1:n), which is
! the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE splint ( xa, ya, y2a, n, x, y )

  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n
  REAL, INTENT(IN)                       :: xa(n)
  REAL, INTENT(IN)                       :: ya(n)
  REAL, INTENT(IN)                       :: y2a(n)
  REAL, INTENT(IN OUT)                   :: x
  REAL, INTENT(OUT)                      :: y

  !-----------------------------------------------------------------------------
  INTEGER                                :: k, khi, klo
  REAL                                   :: a, b, h
  !-----------------------------------------------------------------------------

  klo = 1
  khi = n
1 IF ( khi-klo > 1 ) THEN
     k = ( khi + klo ) / 2
     IF ( xa(k) > x ) THEN
        khi = k
     ELSE
        klo = k
     END IF
     GO TO 1
  END IF
  h = xa(khi) - xa(klo)
  IF ( h == 0.0 ) call paus ('bad xa input in splint')
  a = ( xa(khi) - x ) / h
  b = ( x - xa(klo) ) / h
  y = a*ya(klo) + b*ya(khi) + ( (a**3-a)*y2a(klo)+(b**3-b)*y2a(khi) )*h**2 / 6.0

END SUBROUTINE splint


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE splint_integral
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the in order), and given the array y2a(1:n), which is
! the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE splint_integral ( xa, ya, y2a, n, xlo, xhi, integral )

  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: n
  REAL, INTENT(IN)                       :: xa(n)
  REAL, INTENT(IN)                       :: ya(n)
  REAL, INTENT(IN)                       :: y2a(n)
  REAL, INTENT(IN)                       :: xlo
  REAL, INTENT(IN)                       :: xhi
  REAL, INTENT(OUT)                      :: integral

  !-----------------------------------------------------------------------------
  INTEGER                                :: k, khi_l, klo_l, khi_h, klo_h
  REAL                                   :: xl, xh = 0.0
  REAL                                   :: h, INT, x0, x1, y0, y1, y20, y21
  !-----------------------------------------------------------------------------

  integral = 0.0
  klo_l = 1
  khi_l = n
 do while ( khi_l-klo_l > 1 )
     k = ( khi_l + klo_l ) / 2
     IF ( xa(k) > xlo ) THEN
        khi_l = k
     ELSE
        klo_l = k
     END IF

  end do

  klo_h = 1
  khi_h = n
  do while ( khi_h-klo_h > 1 )
     k = ( khi_h + klo_h ) / 2
     IF ( xa(k) > xhi ) THEN
        khi_h = k
     ELSE
        klo_h = k
     END IF
  end do

  ! integration in spline pieces, the lower interval, bracketed
  ! by xa(klo_L) and xa(khi_L) is in steps shifted upward.

  ! first: determine upper integration bound
  xl = xlo
3 CONTINUE
  IF ( khi_h > khi_l ) THEN
     xh = xa(khi_l)
  ELSE IF ( khi_h == khi_l ) THEN
     xh = xhi
  ELSE
     call paus ('error in spline-integration')
  END IF

  h = xa(khi_l) - xa(klo_l)
  IF ( h == 0.0 ) call paus ('bad xa input in splint')
  x0 = xa(klo_l)
  x1 = xa(khi_l)
  y0 = ya(klo_l)
  y1 = ya(khi_l)
  y20= y2a(klo_l)
  y21= y2a(khi_l)
  ! int = -xL/h * ( (x1-.5*xL)*y0 + (0.5*xL-x0)*y1  &
  !            +y20/6.*(x1**3-1.5*xL*x1*x1+xL*xL*x1-.25*xL**3)  &
  !            -y20/6.*h*h*(x1-.5*xL)  &
  !            +y21/6.*(.25*xL**3-xL*xL*x0+1.5*xL*x0*x0-x0**3)  &
  !            -y21/6.*h*h*(.5*xL-x0) )
  ! int = int + xH/h * ( (x1-.5*xH)*y0 + (0.5*xH-x0)*y1  &
  !            +y20/6.*(x1**3-1.5*xH*x1*x1+xH*xH*x1-.25*xH**3)  &
  !            -y20/6.*h*h*(x1-.5*xH)  &
  !            +y21/6.*(.25*xH**3-xH*xH*x0+1.5*xH*x0*x0-x0**3)  &
  !            -y21/6.*h*h*(.5*xH-x0) )
  INT = -1.0/h * ( xl*((x1-.5*xl)*y0 + (0.5*xl-x0)*y1)  &
       -y20/24.*(x1-xl)**4 + y20/6.*(0.5*xl*xl-x1*xl)*h*h  &
       +y21/24.*(xl-x0)**4 - y21/6.*(0.5*xl*xl-x0*xl)*h*h )
  INT = INT + 1.0/h * ( xh*((x1-.5*xh)*y0 + (0.5*xh-x0)*y1)  &
       -y20/24.*(x1-xh)**4 + y20/6.*(0.5*xh*xh-x1*xh)*h*h  &
       +y21/24.*(xh-x0)**4 - y21/6.*(0.5*xh*xh-x0*xh)*h*h )

  integral = integral + INT
  ! write (*,*) integral,x1,xH
  klo_l = klo_l + 1
  khi_l = khi_l + 1
  xl = x1
  IF (khi_h /= (khi_l-1)) GO TO 3    ! the -1 in (khi_L-1) because khi_L was already counted up

END SUBROUTINE splint_integral
