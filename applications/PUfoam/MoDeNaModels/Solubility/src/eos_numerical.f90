MODULE EOS_NUMERICAL

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: F_NUMERICAL, PHI_NUMERICAL, P_NUMERICAL, H_NUMERICAL

CONTAINS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE F_NUMERICAL

  USE EOS_VARIABLES
  USE EOS_CONSTANTS
  use EOS_polar, only: f_polar
  USE EOS_NUMERICAL_DERIVATIVES, ONLY: ideal_gas, hard_sphere, chain_term, &
       disp_term, hb_term, LC_term, II_term, ID_term
  USE EOS_F_CONTRIBUTIONS, ONLY:F_IDEAL_GAS, F_HARD_SPHERE, F_CHAIN_TPT1,  &
                   F_CHAIN_TPT_D, F_CHAIN_HU_LIU, F_CHAIN_HU_LIU_RC, F_SPT,  &
                   F_DISP_PCSAFT, F_ASSOCIATION, F_ION_DIPOLE_TBH, F_ION_ION_PrimMSA,  &
                   F_ION_ION_nonPrimMSA, F_LC_MayerSaupe, F_pert_theory


  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j
  REAL                                   :: m_mean2
  REAL                                   :: fid, fhs, fdsp, fhc
  REAL                                   :: fhb, fdd, fqq, fdq
  REAL                                   :: fhend, fcc
  REAL                                   :: fbr, flc

  REAL                                   :: eps_kij, k_kij
  !-----------------------------------------------------------------------------

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

  if ( eta < 1.E-100 ) return


  CALL PERTURBATION_PARAMETER

  !-----------------------------------------------------------------------------
  ! overwrite the standard mixing rules by those published by Tang & Gross
  ! using an additional lij parameter
  ! WARNING : the lij parameter is set to lij = - lji in 'para_input'
  !-----------------------------------------------------------------------------
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
             ! * sig_ij(i,j)*uij(i,j)**(1.0/3.0) )**3 * lij(i,j)
             * sig_ij(i,j)*(uij(i,i)*uij(j,j))**(1.0/6.0) )**3 *lij(i,j)
     END DO
  END DO


  !-----------------------------------------------------------------------------
  ! a non-standard mixing rule scaling the hard-sphere term
  ! WARNING : the lij parameter is set to lij = - lji in 'para_input'
  ! (uses an additional lij parameter)
  !-----------------------------------------------------------------------------
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

  !---- ideal gas contribution -------------------------------------------------
  IF ( ideal_gas == 'yes' )       CALL F_IDEAL_GAS ( fid )
  !-----------------------------------------------------------------------------

  !---- hard-sphere contribution -----------------------------------------------
  IF ( hard_sphere == 'CSBM' )    CALL F_HARD_SPHERE ( m_mean2, fhs )
  !-----------------------------------------------------------------------------

  !--- chain term --------------------------------------------------------------
  IF ( chain_term == 'TPT1' )     CALL F_CHAIN_TPT1 ( fhc )
  IF ( chain_term == 'TPT2' )     CALL F_CHAIN_TPT_D ( fhc )
  IF ( chain_term == 'HuLiu' )    CALL F_CHAIN_HU_LIU ( fhc )
  IF ( chain_term == 'HuLiu_rc' ) CALL F_CHAIN_HU_LIU_RC ( fhs, fhc )
  IF ( chain_term == 'SPT' )      CALL F_SPT ( fhs, fhc )
  !-----------------------------------------------------------------------------

  !---- dispersive attraction --------------------------------------------------
  IF ( disp_term == 'PC-SAFT')    CALL F_DISP_PCSAFT ( fdsp )
  ! IF ( disp_term == 'CK')         CALL F_DISP_CKSAFT ( fdsp )
  IF ( disp_term(1:2) == 'PT')    CALL F_pert_theory ( fdsp )
  !-----------------------------------------------------------------------------

  !---- H-bonding contribution / Association -----------------------------------
  IF ( hb_term == 'TPT1_Chap')    CALL F_ASSOCIATION( eps_kij, k_kij, fhb )
  !-----------------------------------------------------------------------------

  !---- polar terms ------------------------------------------------------------
  CALL F_POLAR ( fdd, fqq, fdq )
  !-----------------------------------------------------------------------------

  !---- ion-dipole term --------------------------------------------------------
  IF ( ID_term == 'TBH')          CALL F_ION_DIPOLE_TBH ( fhend )
  !-----------------------------------------------------------------------------

  !---- ion-ion term -----------------------------------------------------------
  IF ( II_term == 'primMSA')      CALL F_ION_ION_PrimMSA ( fcc )
  IF ( II_term == 'nprMSA')       CALL F_ION_ION_nonPrimMSA ( fdd, fqq, fdq, fcc )
  !-----------------------------------------------------------------------------

  !---- liquid-crystal term ----------------------------------------------------
  IF ( LC_term == 'MSaupe')       CALL F_LC_MayerSaupe ( flc )
  !-----------------------------------------------------------------------------

  !=============================================================================
  ! SUBTRACT TERMS (local density approximation) FOR DFT
  !=============================================================================

  !IF ( subtract1 == '1PT')     CALL F_subtract_local_pert_theory ( subtract1, fdsp )
  !IF ( subtract1 == '2PT')     CALL F_subtract_local_pert_theory ( subtract1, fdsp )
  !IF ( subtract2 =='ITTpolar') CALL F_local_ITT_polar ( fdd )
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! residual Helmholtz energy F/(NkT)
  !-----------------------------------------------------------------------------
  fres = fid + fhs + fhc + fdsp + fhb + fdd + fqq + fdq + fcc + flc

  tfr = 0.0

END SUBROUTINE F_NUMERICAL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PHI_NUMERICAL

  USE EOS_VARIABLES
  USE EOS_CONSTANTS
  USE DFT_MODULE, ONLY: z_ges, fres_temp

  !-----------------------------------------------------------------------------
  INTEGER                                :: k
  REAL                                   :: zres, zges
  REAL                                   :: fres1, fres2, fres3, fres4, fres5
  REAL                                   :: delta_rho
  REAL                                   :: delta_factor, delta_rho_thrash
  REAL, DIMENSION(nc)                    :: myres
  REAL, DIMENSION(nc)                    :: rhoi, rhoi_0
  REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4, tfr_5
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! density iteration or pressure calculation
  !-----------------------------------------------------------------------------

  IF (ensemble_flag == 'tp') THEN
     CALL DENSITY_ITERATION
  ELSEIF (ensemble_flag == 'tv') THEN
     eta = eta_start
     CALL P_NUMERICAL
  ELSE
     write (*,*) 'PHI_EOS: define ensemble, ensemble_flag == (tv) or (tp)'
     stop
  END IF

  !-----------------------------------------------------------------------------
  ! compressibility factor z = p/(kT*rho)
  !-----------------------------------------------------------------------------

  zges = (p * 1.E-30) / (kbol*t*eta/z3t)
  IF ( ensemble_flag == 'tv' ) zges = (pges * 1.E-30) / (kbol*t*eta/z3t)
  zres = zges - 1.0
  z_ges = zges

  rhoi_0(1:ncomp) = x(1:ncomp) * eta/z3t
  rhoi(1:ncomp) = rhoi_0(1:ncomp)


  !-----------------------------------------------------------------------------
  ! derivative to rho_k (keeping other rho_i's constant
  !-----------------------------------------------------------------------------

  delta_factor = 1.E-2
  IF (sum(nhb_typ(1:ncomp)) /= 0) delta_rho = 1.E-3
  delta_rho_thrash = 1.E-6

  DO  k = 1, ncomp

     IF ( rhoi_0(k) > delta_rho_thrash ) THEN
        delta_rho = delta_factor * rhoi_0(k)
     ELSE
        delta_rho = delta_rho_thrash * delta_factor
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

     fres4 = 0.0
     fres5 = 0.0
     IF ( rhoi_0(k) > delta_rho ) THEN

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

     IF ( rhoi_0(k) > delta_rho ) THEN
        myres(k) = ( fres5 - 8.0*fres4 + 8.0*fres2 - fres1 ) / ( 6.0*delta_rho )
     ELSE
        myres(k) = ( -3.0*fres3 + 4.0*fres2 - fres1 ) / delta_rho
     END IF

  END DO


  !-----------------------------------------------------------------------------
  ! residual Helmholtz energy
  !-----------------------------------------------------------------------------

  fres_temp = fres

  !-----------------------------------------------------------------------------
  ! residual chemical potential
  !-----------------------------------------------------------------------------

  DO  k = 1, ncomp
     IF (ensemble_flag == 'tp') lnphi(k) = myres(k) - LOG(zges)
     IF (ensemble_flag == 'tv' .AND. eta >= 0.0) lnphi(k) = myres(k) !+LOG(rho)
     ! write (*,*) 'in',k,EXP(lnphi(k)),LOG(zges),eta
     ! IF (DFT.GE.98) write (*,*) dft
     ! write (*,*) 'lnphi',k,LNPHI(k),x(k),MYRES(k), -LOG(ZGES)
     ! pause
  END DO

END SUBROUTINE PHI_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_NUMERICAL

  USE EOS_VARIABLES
  USE utilities

  !-----------------------------------------------------------------------------
  REAL                                   :: dzetdv, eta_0, dist, fact
  REAL                                   :: fres1, fres2, fres3, fres4, fres5
  REAL                                   :: df_dr, df_drdr, pideal, dpiddz
  REAL                                   :: tfr_1, tfr_2, tfr_3, tfr_4, tfr_5
  !-----------------------------------------------------------------------------


  IF (eta > 0.1) THEN
     fact = 1.0
  ELSE IF (eta <= 0.1 .AND. eta > 0.01) THEN
     fact = 10.0
  ELSE
     fact = 10.0
  END IF
  dist = eta*3.E-3 *fact
  ! dist = eta*4.E-3 *fact

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

  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------

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

  !------------------p ideal----------------------------------------------------
  pideal = rho * (kbol*t) / 1.E-30

  !------------------p summation, p comes out in Pa ----------------------------
  pges   = pideal + pges

END SUBROUTINE P_NUMERICAL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE H_numerical

  USE PARAMETERS, ONLY: RGAS
  USE EOS_VARIABLES

  !-----------------------------------------------------------------------------
  REAL                                   :: dist, fact, rho_0
  REAL                                   :: fres1,fres2,fres3,fres4,fres5
  REAL                                   :: f_1, f_2, f_3, f_4
  REAL                                   :: cv_res, t_tmp, zges
  REAL                                   :: f_dt, f_dtdt, f_dr, f_drdr, f_drdt
  !-----------------------------------------------------------------------------


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

END MODULE EOS_NUMERICAL
