!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE GC_DATA
!
! This module contains the experimental data as well as the calculated
! properties
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
MODULE GC_DATA

  implicit none
  save

  !-----------------------------------------------------------------------------
  ! storing experimental data
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER                    :: MMX = 3500

  INTEGER                               :: nlv,nliq
  REAL, DIMENSION(MMX)                  :: plv, tlv, pliq, tliq, vliq, v_crit1, v_crit2
  CHARACTER (LEN=30), DIMENSION(MMX)    :: gc_comp1, gc_comp2

  !-----------------------------------------------------------------------------
  ! storing/transferring calculated properties
  !-----------------------------------------------------------------------------

  INTEGER                               :: n_liq(20), n_lv(20), k_nr
  REAL, DIMENSION(MMX)                  :: plvcal, v_cal
  CHARACTER (LEN=30)                    :: subst(20)

END MODULE GC_DATA





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE GC_GROUP
!
! This module contains parameters detailing a GC-group
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
MODULE GC_GROUP

  USE BASIC_VARIABLES, ONLY: nc
  implicit none
  save

  !-----------------------------------------------------------------------------
  ! N_comp: maximum number of components i considered in the GC-fit
  ! N_grs : maximum number of groups j for a component
  !         Index j runs j = 1, gc_tot_nr_groups(i) for specific comp. i
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER                    :: N_comp = nc
  INTEGER, PARAMETER                    :: N_grs = 20

  REAL, DIMENSION(N_comp,N_grs)          :: m_gc, msig3, meps, mu_gc, mm_gc

  INTEGER, DIMENSION(N_grs)             :: n_groups
  INTEGER, DIMENSION(N_comp,N_grs)      :: gc_groups
  INTEGER, DIMENSION(N_comp)            :: gc_tot_nr_groups
  CHARACTER (LEN=20), DIMENSION(N_comp,N_grs) :: gc_grouptype
  !  CHARACTER (LEN=30)                    :: comp_name(N_comp)

END MODULE GC_GROUP





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

MODULE MOD_PARA_TRANSFER
  implicit none
  save

  REAL, ALLOCATABLE                       :: para_transfer(:)

END MODULE MOD_PARA_TRANSFER






!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE GC_PARAMETER
!
! This subroutine contains the group contribution parameters of the
! individual chemical groups. It also combines the group-parameters of
! a molecule to the molecular PCP-SAFT parameters.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE GC_PARAMETER ( n_species, comp_name )

  USE BASIC_VARIABLES
  USE GC_GROUP
  USE MOD_PARA_TRANSFER
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                      :: n_species
  CHARACTER (LEN=30), INTENT(IN)           :: comp_name(N_comp)
  !-----------------------------------------------------------------------------
  INTEGER                                  :: i, j, k, jj, no
  INTEGER, DIMENSION(nc)                   :: nhb_typ
  INTEGER, DIMENSION(nc,nsite)             :: nhb_no
  REAL, DIMENSION(nc,nc,nsite,nsite)       :: eps_hb
  REAL, DIMENSION(nc,nc)                   :: kap_hb
  REAL, DIMENSION(N_comp,N_grs)            :: epsilon_hb, kappa_hb
  !-----------------------------------------------------------------------------

  n_groups  = 0
  gc_grouptype = ' '
  gc_tot_nr_groups = 0

  m_gc  = 0.0
  msig3 = 0.0
  meps  = 0.0
  mu_gc = 0.0
  mm_gc = 0.0
  epsilon_hb = 0.0
  kappa_hb = 0.0

  !-----------------------------------------------------------------------------
  ! draw structural groups of molecules from file
  !-----------------------------------------------------------------------------

  CALL GET_GC_GROUPS ( n_species, comp_name ) !, gc_tot_nr_groups, gc_groups, gc_grouptype)
  !write (*,*) n_species, comp_name

  !-----------------------------------------------------------------------------
  ! group-contribution parameters
  !-----------------------------------------------------------------------------

  DO i = 1, n_species
     DO j = 1, gc_tot_nr_groups(i)
        IF (gc_grouptype(i,j) == '-CH3') THEN  !  CH3-
           m_gc(i,j) = 0.695765
           msig3(i,j)= 34.0119
           meps(i,j) = 147.516
           mm_gc(i,j) = 15.035
        ELSE IF (gc_grouptype(i,j) == '-CH2-') THEN  ! -CH2-
           m_gc(i,j) = 0.429404
           msig3(i,j)= 25.08362
           meps(i,j) = 108.1502
           mm_gc(i,j) = 14.02675     ! ?
        ELSE IF (gc_grouptype(i,j) == '-CH<(branch)') THEN
           m_gc(i,j) = 0.143
           msig3(i,j)= 16.41211
           meps(i,j) = 49.71252
           mm_gc(i,j) = 13.015     ! ?
        ELSE IF (gc_grouptype(i,j) == '>CH<(2branch)') THEN
           m_gc(i,j) = -0.6700
           msig3(i,j)= 3.82853
           meps(i,j) = -72.1456
           mm_gc(i,j) = 12.003     ! ?
        ELSE IF (gc_grouptype(i,j) == 'CH4') THEN  ! methane
           m_gc(i,j) = 1.000
           msig3(i,j)= 3.70388767**3
           meps(i,j) = 150.033987
           mm_gc(i,j) = 16.043
        ELSE IF (gc_grouptype(i,j) == '-Cl') THEN
           m_gc(i,j) = 1.5514 / 2.0
           msig3(i,j)=  1.5514 / 2.0 * 3.3672**3
           meps(i,j) = 265.67 * 1.5514 / 2.0
           mm_gc(i,j) = 35.453
        ELSE IF (gc_grouptype(i,j) == 'N(N2)') THEN  ! nitrogen
           m_gc(i,j) = 1.205275098 / 2.0
           msig3(i,j)=  1.205275098 / 2.0 * 3.3129702**3
           meps(i,j) = 90.9606924 * 1.205275098 / 2.0
           mm_gc(i,j) = 14.005
        ELSE IF (gc_grouptype(i,j) == '-CH2-(cyclic6)') THEN  ! -CH2-  C6-cyclic
           m_gc(i,j) = 0.4217
           msig3(i,j)= m_gc(i,j) * 3.84990887**3
           meps(i,j) = m_gc(i,j) * 278.108786
           mm_gc(i,j) = 14.02675
        ELSE IF (gc_grouptype(i,j) == '-CH2-(cyclic5)') THEN  ! -CH2-  C6-cyclic
           m_gc(i,j) = 0.473
           msig3(i,j)= m_gc(i,j) * 3.71139254**3
           meps(i,j) = m_gc(i,j) * 265.828755
           mm_gc(i,j) = 14.02675
        ELSE IF (gc_grouptype(i,j) == '-CH-(aromatic)') THEN  ! -CH-  C6-aromatic
           m_gc(i,j) = 0.4108
           msig3(i,j)= 19.93845
           meps(i,j) = 118.0598
           mm_gc(i,j) = 13.015
        ELSE IF (gc_grouptype(i,j) == '-C-R(aromatic)') THEN  ! -C-R  C6-aromatic
           m_gc(i,j) = 0.026654
           msig3(i,j)= 14.1613
           meps(i,j) = 52.5285
           mm_gc(i,j) = 12.015
        ELSE IF (gc_grouptype(i,j) == '-O-ether') THEN
           m_gc(i,j) = 0.81453
           msig3(i,j)= 10.94506
           meps(i,j) = 162.6272
           mu_gc(i,j) = 1.7488
           mm_gc(i,j) = 15.998
        ELSE IF (gc_grouptype(i,j) == '-OH') THEN
           m_gc(i,j) = 0.6409
           msig3(i,j)= m_gc(i,j) * ( 2.7707 )**3
           meps(i,j) = m_gc(i,j) * 376.99
           mu_gc(i,j) = 1.7398
           mm_gc(i,j) = 17.007
           epsilon_hb(i,j) = 2374.0
           kappa_hb(i,j) = 0.00698
           m_gc(i,j) = 1.4633
           msig3(i,j)= m_gc(i,j) * ( 2.5320 )**3
           meps(i,j) = m_gc(i,j) * 198.58
           mu_gc(i,j) = 1.3843
           mm_gc(i,j) = 17.007
           epsilon_hb(i,j) =        para_transfer(1)   !2431.3
           kappa_hb(i,j) = 0.03256
        ELSE IF (gc_grouptype(i,j) == 'CO_Acetone') THEN
           m_gc(i,j) = 1.138
           msig3(i,j)= 14.0
           meps(i,j) = 200.0
           mu_gc(i,j) = 2.88
           mm_gc(i,j) = 28.010
        ELSE IF (gc_grouptype(i,j) == '-CH2-COO-CH2-') THEN  ! -CH2-COO-CH2- ester-group (from acid-side)
           m_gc(i,j) = 2.0000
           msig3(i,j)= 82.684
           meps(i,j) = 491.00
           mu_gc(i,j) = 3.674
           mm_gc(i,j) = 72.06   ! geschaetzt
           ! m_gc(i,j) = 2.5        + ( para_transfer(1) - 10.0 )
           ! msig3(i,j)= 100.0      + ( para_transfer(2) - 10.0 )*20.0
           ! meps(i,j) = 600.0      + ( para_transfer(3) - 10.0 )*100.0
           ! mu_gc(i,j) = 2.0       + ( para_transfer(4) - 10.0 )
           ! mm_gc(i,j) = 72.06   ! geschaetzt
        ELSE IF (gc_grouptype(i,j) == '-CH2-COO-CH3') THEN  ! -CH2-COO-CH3 methyl ester-group (from acid-side)
           m_gc(i,j) = 2.5386
           msig3(i,j)= 91.500
           meps(i,j) = 582.58
           mu_gc(i,j) = 3.415
           mm_gc(i,j) = 73.07   ! geschaetzt
           ! write (*,*) m_gc(i,j),msig3(i,j),meps(i,j),mu_gc(i,j)
           ! stop
        ELSE
           WRITE (*,*) 'GC_PARAMETER: GC-parameter not found !', gc_grouptype(i,j)
           STOP
        END IF
     END DO
  END DO

  !-----------------------------------------------------------------------------
  ! calculate pure component parameters (i.e. molecular parameters)
  !-----------------------------------------------------------------------------

  DO i = 1, n_species
     parame(i,1) = 0.0
     parame(i,2) = 0.0
     parame(i,3) = 0.0
     parame(i,6) = 0.0
     mm(i)       = 0.0
     DO j = 1, gc_tot_nr_groups(i)
        parame(i,1) = parame(i,1) + m_gc(i,j)  * REAL(gc_groups(i,j))
        parame(i,2) = parame(i,2) + msig3(i,j) * REAL(gc_groups(i,j))
        parame(i,3) = parame(i,3) + meps(i,j)  * REAL(gc_groups(i,j))
        parame(i,6) = parame(i,6) + mu_gc(i,j) * REAL(gc_groups(i,j))
        mm(i)       = mm(i)       + mm_gc(i,j) * REAL(gc_groups(i,j))
        if ( epsilon_hb(i,j) /= 0.0 ) then
           if ( n_species > 1 ) write (*,*) 'extend GC method'
           if ( n_species > 1 ) stop
           eps_hb = 0.0
           eps_hb(1,1,1,2) = epsilon_hb(i,j)
           eps_hb(1,1,2,1) = eps_hb(1,1,1,2)
           kap_hb(i,i)= kappa_hb(i,j)
           nhb_typ(1)  = 2
           nhb_no(1,1) = 1
           nhb_no(1,2) = 1
           IF (nhb_typ(i) /= 0) THEN
              parame(i,12) = REAL(nhb_typ(i))
              parame(i,13) = kap_hb(i,i)
              no = 0
              DO jj=1,nhb_typ(i)
                 DO k=1,nhb_typ(i)
                    parame(i,(14+no))=eps_hb(i,i,jj,k)
                    no=no+1
                 END DO
              END DO
              DO jj=1,nhb_typ(i)
                 parame(i,(14+no))=REAL(nhb_no(i,jj))
                 no=no+1
              END DO
           ELSE
              DO k=12,25
                 parame(i,k)=0.0
              END DO
           END IF
        end if
     END DO
     parame(i,2) = ( parame(i,2)/parame(i,1) )**(1.0/3.0)
     parame(i,3) =   parame(i,3)/parame(i,1)
  END DO

  DO i = 1, n_species
     IF (parame(i,1) == 0.0) THEN
        WRITE (*,*) 'error in GC_PARAMETER ',comp_name(i)
        STOP
     END IF
  END DO

  ! write (*,*) parame(1,1),parame(1,2),parame(1,3),parame(1,6),mm(1)
  ! pause

END SUBROUTINE GC_PARAMETER



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE GET_GC_GROUPS
!
! This subroutine assignes chemical groups to each species. The species
! is identified by name. The output of the subroutine is:
!
! *  gc_grouptype      : type of a group
! *  gc_groups         : numer of groups of a certain type
! *  gc_tot_nr_groups  : total number of groups of a substance i
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE GET_GC_GROUPS ( n_species, comp_name )

  USE GC_GROUP
  USE BASIC_VARIABLES, ONLY: nc
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                      :: n_species
  CHARACTER (LEN=30), INTENT(IN)           :: comp_name(N_comp)

  !-----------------------------------------------------------------------------
  INTEGER                                  :: i, j, k
  INTEGER                                  :: n_comps, tot_n_groups !, n_groups(N_grs)
  CHARACTER (LEN=30)                       :: c_name
  CHARACTER (LEN=50)                       :: gc_file
  CHARACTER (LEN=20)                       :: dum, grouptype(N_grs)
  !-----------------------------------------------------------------------------

  grouptype = ' '
  gc_groups = 0

  gc_file = './input_file/gc_method/gc_groups.inp'
  CALL file_open(gc_file,85)
  READ(85,*) dum, n_comps
  READ(85,*) dum
  READ(85,*) dum
  DO k = 1, n_comps
     READ(85,*) c_name, tot_n_groups, n_groups(1), dum, grouptype(1)
     DO j = 2, tot_n_groups
        READ(85,*) n_groups(j), dum, grouptype(j)
     END DO
     DO i = 1, n_species
        IF (c_name == comp_name(i) ) THEN
           gc_tot_nr_groups(i) = tot_n_groups
           DO j = 1, gc_tot_nr_groups(i)
              gc_groups(i,j) = n_groups(j)
              gc_grouptype(i,j) = grouptype(j) ! trim( adjustl( grouptype(j) ))
           END DO
        END IF
     END DO
  END DO
  CLOSE (85)

END SUBROUTINE GET_GC_GROUPS




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE GC_RESIDUAL
! 
! Subroutine for fitting pure component parameters. The residual between 
! calculated values and experimental data is calculated and stored in
! vector 'fvec'. The vector 'fvec' serves as the objective function for
! the parameter regression.
! GC_RESIDUAL is called by GC_FITTING via the Levenberg-Marquardt routine.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE GC_RESIDUAL (lm_m_dat, lm_n_par, para, fvec, iflag)

  USE PARAMETERS, ONLY: RGAS
  USE BASIC_VARIABLES
  USE GC_DATA
  USE GC_GROUP, ONLY: N_comp
  USE MOD_PARA_TRANSFER
  use utilities, only: dens_calc, SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! INTEGER, PARAMETER                      :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, INTENT(IN)                     :: lm_m_dat
  INTEGER, INTENT(IN)                     :: lm_n_par
  REAL, INTENT(IN)                        :: para(:)
  REAL, INTENT(IN OUT)                    :: fvec(:)
  INTEGER, INTENT(IN OUT)                 :: iflag

  !-----------------------------------------------------------------------------
  INTEGER                                 :: i, j_dat
  INTEGER                                 :: converg, crit_dat
  REAL                                    :: weigh, tc
  REAL                                    :: toterr, plv_dum, v_dum
  REAL                                    :: density(np), w(np,nc)
  REAL                                    :: rho_phas(np)
  LOGICAL                                 :: scan, vapor, nocon = .false.
  CHARACTER (LEN=30)                      :: comp_name(N_comp)
  CHARACTER (LEN=4)                       :: char_len
  !-----------------------------------------------------------------------------

  ALLOCATE( para_transfer(lm_n_par) )
  para_transfer = para

  fvec = 0.0

  parame(1,5) = 0.12

  WRITE (char_len,'(I3)') lm_n_par
  IF (iflag == 1) WRITE (*,'(a,'//char_len//'(1x,G21.14))') ' parameter ', (para(i),i=1,lm_n_par)

  !-----------------------------------------------------------------------------
  ! liquid densit data
  !-----------------------------------------------------------------------------
  DO i = 1, nliq

     comp_name(1) = gc_comp1(i)
     IF ( i == 1 ) THEN
        crit_dat = 0
        tc = 0.0
        CALL GC_PARAMETER ( 1, comp_name )
     ELSE
        IF (gc_comp1(i) /= gc_comp1(i-1)) THEN
           crit_dat = 0
           tc = 0.0
           CALL GC_PARAMETER ( 1, comp_name )
        ENDIF
     END IF


     IF ( pliq(i) > 0.0 ) THEN

        !-----------------------------------------------------------------------
        ! case 1: calculate rho to given (p,T)
        !-----------------------------------------------------------------------

        nocon   = .false.
        scan    = .false.

        j_dat   = i
        nphas   = 1
        t       = tliq(i)
        p       = pliq(i)
        xi(1,1) = 1.0
        IF (vliq(i) <= v_crit1(i)) THEN
           vapor = .false.
           densta(1) = 0.5
           ! IF (d_kond(1,i).NE.0.0) densta(1) = d_kond(1,i)
        ELSE
           vapor = .true.
           densta(1) = 0.002
           ! IF (d_kond(1,i).NE.0.0) densta(1) = d_kond(1,i)
        END IF

38      CONTINUE
        IF (scan) densta(1) = densta(1) + 0.002

        CALL dens_calc ( rho_phas )
        CALL SI_dens ( density, w )
        v_cal(i) = 1.0 / density(1)

        IF ((v_cal(i)-vliq(i))*1.d2/vliq(i) < 5.0 .AND. .NOT. scan) THEN
           d_kond(1,i) = dense(1)
        ELSE
           d_kond(1,i) = 0.0
        END IF


        IF ( vapor .AND. (1.0/density(1)) <= (v_crit1(i)/1.1) ) THEN
           IF ( densta(1) <= 0.228 ) THEN
              scan = .true.
              GO TO 38
           ELSE
              nocon = .true.
           END IF
        END IF
        IF ( dense(1) > 0.55 ) WRITE (*,*) 'density ',i,dense(1)

        weigh = 0.5
        IF (nocon) weigh = weigh / 4.0

        fvec(j_dat) = weigh * (v_cal(i)-vliq(i))/  &
             ( SQRT(ABS(vliq(i))) * SQRT(ABS(v_cal(i))) )



     ELSE

        !-----------------------------------------------------------------------
        ! case 2, when pliq(i) = 0, indicating rho_L at coexistence
        !-----------------------------------------------------------------------

        nphas = 2

        j_dat = i

        val_init(1) = 0.5
        val_init(2) = 1.E-8
        val_init(3) = tliq(i)
        val_init(4) = 10.0
        IF (d_kond(1,j_dat) /= 0.0) val_init(1)=d_kond(1,j_dat)
        IF (d_kond(2,j_dat) /= 0.0) val_init(2)=d_kond(2,j_dat)
        IF (d_kond(1,j_dat) /= 0.0 .AND. d_kond(2,j_dat) /= 0.0) val_init(4) = plv_kon(j_dat)
        val_init(5) = 0.0      ! mole fraction lnx(1,1)
        val_init(6) = 0.0      ! mole fraction lnx(2,1)

        CALL pure_equilibrium_fit ( crit_dat, tc, tliq(i), pliq(i), converg, v_cal(i), plv_dum )

        IF (converg == 1) THEN
           d_kond(1,j_dat) = dense(1)
           d_kond(2,j_dat) = dense(2)
           plv_kon(j_dat)  = val_conv(4)
        ELSE
           d_kond(1,j_dat) = 0.0
           d_kond(2,j_dat) = 0.0
           plv_kon(j_dat)  = 10.0
           plv_dum = 1.E5
           v_cal(i) = vliq(i) * 0.75
           IF (crit_dat == 1 .AND. tliq(i) > tc ) v_cal(i) = vliq(i)*1.1
        END IF

        weigh = 0.5
        IF ( nocon ) weigh = weigh / 4.0

        fvec(j_dat) = weigh * (v_cal(i)-vliq(i))/  &
             (SQRT(ABS(vliq(i)))*SQRT(ABS(v_cal(i))))

     END IF

  END DO
  !-----------------------------------------------------------------------------
  ! end: liquid density data
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  ! vapor pressure data
  !-----------------------------------------------------------------------------
  DO  i = 1, nlv

     comp_name(1) = gc_comp2(i)
     IF ( i == 1 ) THEN
        crit_dat = 0
        tc = 0.0
        CALL GC_PARAMETER ( 1, comp_name )
     ELSE
        IF(gc_comp2(i) /= gc_comp2(i-1)) THEN
           crit_dat = 0
           tc = 0.0
           CALL GC_PARAMETER ( 1, comp_name )
        ENDIF
     END IF

     nphas = 2

     j_dat = i + nliq

     val_init(1) = 0.5
     val_init(2) = 1.d-8
     val_init(3) = tlv(i)
     val_init(4) = plv(i)
     IF (d_kond(1,j_dat) /= 0.0) val_init(1) = d_kond(1,j_dat)
     IF (d_kond(2,j_dat) /= 0.0) val_init(2) = d_kond(2,j_dat)
     IF (d_kond(1,j_dat) /= 0.0 .AND. d_kond(2,j_dat) /= 0.0) val_init(4) = plv_kon(j_dat)
     val_init(5) = 0.0      ! mole fraction lnx(1,1)
     val_init(6) = 0.0      ! mole fraction lnx(2,1)

     CALL pure_equilibrium_fit ( crit_dat, tc, tlv(i), plv(i), converg, v_dum, plvcal(i) )

     IF (converg == 1) THEN
        d_kond(1,j_dat) = dense(1)
        d_kond(2,j_dat) = dense(2)
        plv_kon(j_dat)  = val_conv(4)
     ELSE
        d_kond(1,j_dat) = 0.0
        d_kond(2,j_dat) = 0.0
        plv_kon(j_dat)  = plv(i)
        IF ( crit_dat == 1 .AND. tc > tlv(i) ) THEN
           plvcal(i) = plv(i) * 0.75
           ! IF (i > 2) plvcal(i) = EXP( LOG(plvcal(i-1)) + (LOG(plvcal(i-1))-LOG(plvcal(i-2)))
           !                             /(1.0/tlv(i-1)-1.0/tlv(i-2)) * (1.0/tlv(i)-1.0/tlv(i-1)) )
           WRITE (*,'(a,i4,4G15.5)') ' loose ! ', i, tlv(i), t, tc, plv(i)
        ELSE
           plvcal(i) = plv(i)
           IF (i > 2) THEN
              IF (plvcal(i-2) > 0.0 .AND. plvcal(i-1) > 0.0 .AND.  &
                   tlv(i-1) /= tlv(i-2)) THEN
                 plvcal(i) = EXP( LOG(plvcal(i-1)) + (LOG(plvcal(i-1))-LOG(plvcal(i-2)))  &
                      /(1.0/tlv(i-1)-1.0/tlv(i-2)) *(1.0/tlv(i)-1.0/tlv(i-1)) )
              END IF
           END IF
           ! WRITE (*,*) ' supercritical at point ', i, tlv(i), tc
        END IF
     END IF

     fvec(j_dat) = 0.2 * ABS( plvcal(i)-plv(i) )/  &
          ((SQRT(ABS(plv(i)))*SQRT(ABS(plvcal(i)))))**0.85

     IF (tc < tlv(i).AND.crit_dat == 1) THEN
        fvec(j_dat) = fvec(j_dat) + 5.0 * ABS(tlv(i)-tc) + 2.0
        WRITE(*,'(a,i3,1x,a22,2(f15.5))') ' skipped ! ',i,gc_comp2(i), tlv(i),tc
     END IF

  END DO
  !-----------------------------------------------------------------------------
  ! end: vapor pressure data
  !-----------------------------------------------------------------------------


  toterr = 0.0
  DO i = 1, nlv+nliq
     toterr = toterr + fvec(i)**2
  END DO
  WRITE (*,*) '--- error ---',toterr

  IF ( j_dat /= lm_m_dat ) WRITE (*,*) 'ERROR in GC_RESIDUAL', j_dat, lm_m_dat

  DEALLOCATE( para_transfer )

END SUBROUTINE GC_RESIDUAL





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE GC_FITTING
!
! Driver for fitting group-contribution parameters of the PCP-SAFT model.
! It uses a Levenberg-Marquardt algorithm to minimize deviations of
! calculated properties to experimental data.
! Under construction: two issues are so far implemented very temprarily
! and need upgrading: (1) the components that are considered in the
! regression are specified in the following subroutine (hard wired !!!)
! and (2) the number of adjustable parameters (lm_n_par) is hard-wired in
! this subroutine !!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE GC_FITTING

  USE BASIC_VARIABLES
  USE Levenberg_Marquardt
  USE GC_DATA
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE GC_RESIDUAL(lm_m_dat, lm_n_par, para, fvec, iflag)
       IMPLICIT NONE
       INTEGER, INTENT(IN)                :: lm_m_dat, lm_n_par
       REAL, INTENT(IN)                   :: para(:)
       REAL, INTENT(IN OUT)               :: fvec(:)
       INTEGER, INTENT(IN OUT)            :: iflag
     END SUBROUTINE GC_RESIDUAL
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                :: info
  INTEGER                                :: lm_n_par         ! number of adjustable parameters
  INTEGER                                :: lm_m_dat         ! number of data points
  INTEGER, ALLOCATABLE                   :: ipvt(:)
  REAL, ALLOCATABLE                      :: para(:)
  REAL, ALLOCATABLE                      :: fvec(:)
  REAL                                   :: tol    = 1.0E-6
  REAL                                   :: epsfcn
  REAL                                   :: lm_factor

  INTEGER                                :: i, k
  INTEGER                                :: n0
  REAL                                   :: av_dev, max_dev, rms
  !-----------------------------------------------------------------------------


  ncomp     = 1

  ! numerical constants for Levenberg-Marquardt algorithm
  ! ----------------------------------------------------------------------
  scaling(1) = 1.E2
  IF ( num == 0 ) epsfcn = 1.0E-6**2  ! sqrt of relat. step size (finite differences)
  IF ( num == 1 ) epsfcn = 1.0E-5**2  ! sqrt of relat. step size (finite differences)
  lm_factor = 0.01                    ! maximum initial step size (parameter iteration)

  ! ======================================================================
  lm_n_par = 4
  ! ======================================================================


  ! initialize parameter vectors
  ! ----------------------------------------------------------------------
  ALLOCATE( fvec(lm_m_dat) )
  ALLOCATE( para(lm_n_par) )
  ALLOCATE( ipvt(lm_n_par) )

  DO i = 1, lm_n_par
     para(i) = 10.0
  END DO

  ! ----------------------------------------------------------------------
  ! read input-files
  ! ----------------------------------------------------------------------

  CALL GC_DATIN ( lm_m_dat )

  ! ----------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  ! ----------------------------------------------------------------------

  CALL lmdif1(GC_RESIDUAL,lm_m_dat,lm_n_par,para,fvec,tol,epsfcn,lm_factor,info,ipvt)



  ! ----------------------------------------------------------------------
  ! writing results to files
  ! ----------------------------------------------------------------------

  DO k = 1, k_nr

     OPEN (23,FILE='./output_file/gc_method/'//trim( adjustl( subst(k) ))//'.dat')
     WRITE(23,*) subst(k)
     WRITE(23,*) '----------------------------------------------------'
     ! --------------------------------------------------------------------
     ! writing experimental liquid density data
     ! --------------------------------------------------------------------
     WRITE(23,*) ' liquid density data'
     WRITE(23,*) ' no  T/K  p/bar rho/(kg/m**3) rho_cal/(kg/m**3) err%'
     av_dev  = 0.0
     max_dev = 0.0
     rms     = 0.0
     n0 = sum( n_liq(1:k-1) )
     DO i = n0+1 , n0+n_liq(k)
        WRITE(23,'( I3, 2X, 5(2x,G14.6) )') i-n0, tliq(i),  &
             pliq(i)/1.d5,1.0/vliq(i),1.0/v_cal(i), (vliq(i)/v_cal(i)-1.0)*100.0
        max_dev = MAX( max_dev, ABS(vliq(i)/v_cal(i)-1.0)*100.0 )
        rms     = rms + (vliq(i)/v_cal(i)-1.0)**2
        av_dev  = av_dev +  ABS(vliq(i)/v_cal(i)-1.0)*100.0
     END DO
     av_dev = av_dev / REAL(n_liq(k))
     rms    = SQRT( rms / REAL(n_liq(k)) )*100.0
     WRITE(23,*) ' '
     WRITE(23,*) 'Deviation of calculated to exp. density data'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F7.3, A/))') 'Max. Deviation:',max_dev, ' %',  &
          'Average Deviation: ',av_dev, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '

     ! --------------------------------------------------------------------
     ! writing experimental vapor pressure data
     ! --------------------------------------------------------------------
     WRITE(23,*) ' vapor pressure data'
     WRITE(23,*) ' no  T/K  psat_exp/bar  psat_cal/bar  err%'
     av_dev  = 0.0
     max_dev = 0.0
     rms     = 0.0
     n0 = sum( n_lv(1:k-1) )
     DO i = n0+1, n0+n_lv(k)
        WRITE(23,'( I3, 2X, 4(2x,G14.6) )') i-n0, tlv(i), plv(i)/1.E5,  &
             plvcal(i)/1.E5,(plvcal(i)-plv(i))/plv(i)*100.0
        max_dev = MAX( max_dev, ABS(plvcal(i)-plv(i))/plv(i)*100.0 )
        rms     = rms + ( (plvcal(i)-plv(i))/plv(i) )**2
        av_dev  = av_dev +  ABS(plvcal(i)-plv(i))/plv(i) * 100.0
     END DO
     av_dev = av_dev / REAL(n_lv(k))
     rms    = SQRT( rms / REAL(n_lv(k)) )*100.0
     WRITE(23,*) ' '
     WRITE(23,*) 'Deviation of calculated to exp. vapor pressure data'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F7.3, A/))') 'Max. Deviation:',max_dev, ' %',  &
          'Average Deviation: ',av_dev, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '


     CLOSE(23)
  END DO

  DEALLOCATE( fvec, para, ipvt )

END SUBROUTINE GC_FITTING



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE GC_DATIN ( lm_m_dat )
!
! Reading experimental data from pure component input files. Temporarily,
! the species to be considered in the regression are specified here !!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE GC_DATIN ( lm_m_dat )

  USE BASIC_VARIABLES
  USE GC_DATA
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(OUT)                     :: lm_m_dat
  !-----------------------------------------------------------------------------
  INTEGER :: datnr, i, k, n0
  REAL :: v_crit, tlv_tmp, plv_tmp
  CHARACTER (LEN=1) :: info*80
  !-----------------------------------------------------------------------------

  eos = 1

  !$$$      k_nr = 7
  !$$$      subst(1) = 'propane'
  !$$$      subst(2) = 'butane'
  !$$$      subst(3) = 'pentane'
  !$$$      subst(4) = 'hexane'
  !$$$      subst(5) = 'heptane'
  !$$$      subst(6) = 'octane'
  !$$$      subst(7) = 'decane'

  !$$$      k_nr = 4
  !$$$      subst(1) = 'toluene'
  !$$$      subst(2) = 'ethylbenzene'
  !$$$      subst(3) = 'npropylbenzene'
  !$$$      subst(4) = 'nbutylbenzene'

!!$k_nr = 8
!!$subst(1) = 'methyl-propanoate'
!!$subst(2) = 'ethyl-propanoate'
!!$subst(3) = 'npropyl-propanoate'
!!$! subst(6) = 'nbutyl-propionate'
!!$subst(4) = 'methyl-butanoate'
!!$subst(5) = 'ethyl-butanoate'
!!$subst(6) = 'npropyl-butanoate'
!!$subst(7) = 'methyl-decanoate'
!!$subst(8) = 'methyl-dodecanoate'

  k_nr = 1
  subst(1) = 'dimethyl-ether'


  nliq = 0
  nlv  = 0

  DO k = 1, k_nr

     OPEN (22,FILE='./input_file/gc_method/pure_comp/'  &
          //trim( adjustl( subst(k) ))//'.dat')

     !--------------------------------------------------------------------------
     ! reading experimental pure component data
     !--------------------------------------------------------------------------
     DO i=1,16
        READ(22,'(A80)') info
     END DO
     READ(22,*) mm(1), v_crit    ! note molecular mass mm(i) is ignored (overwritten later)

     DO i=1,15
        READ(22,'(A80)') info
     END DO
     READ(22,*) n_liq(k) ,n_lv(k)
     READ(22,'(A80)') info
     READ(22,'(A80)') info


     !--------------------------------------------------------------------------
     ! reading experimental liquid density data
     !--------------------------------------------------------------------------
     n0 = nliq
     DO i = n0+1 , n0+n_liq(k)
        READ(22,*) datnr,pliq(i),tliq(i),vliq(i)
        nliq=i
        gc_comp1(i) = subst(k)
        v_crit1(i)  = v_crit
     END DO
     READ(22,'(A80)') info
     READ(22,'(A80)') info

     !--------------------------------------------------------------------------
     ! reading experimental vapor pressure data
     !--------------------------------------------------------------------------
     n0 = nlv
     DO i = n0+1, n0+n_lv(k)
        READ(22,*) datnr,plv_tmp,tlv_tmp
        IF (plv_tmp > 20.0) THEN
           nlv = nlv + 1
           tlv(nlv) = tlv_tmp
           plv(nlv) = plv_tmp
           gc_comp2(nlv) = subst(k)
           v_crit2(nlv)  = v_crit
        END IF
     END DO
     n_lv(k) = nlv-n0

     CLOSE(22)

  END DO

  lm_m_dat = nliq + nlv

END SUBROUTINE GC_DATIN

