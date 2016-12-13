!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE pure_fit_parameters
!
! This module contains parameters and variables needed for the pure
! component parameter regression.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

Module pure_fit_parameters

  implicit none
  save

  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER          :: cont = 4           ! number of data types (e.g. p^sat., rho^liq,...) 
  INTEGER, PARAMETER          :: npmax = 40
  INTEGER, PARAMETER          :: MMX = 3500

  INTEGER                     :: nlv, nliq, n_vir, n_en
  REAL, DIMENSION(npmax)      :: rel
  REAL, DIMENSION(cont)       :: type_weight

  REAL, DIMENSION(MMX)        :: plvcal, v_cal, b_cal, u_cal
  REAL, DIMENSION(MMX)        :: plv, tlv
  REAL, DIMENSION(MMX)        :: pliq, tliq, vliq
  REAL, DIMENSION(MMX)        :: b_vir, t_vir
  REAL, DIMENSION(MMX)        :: u_en, t_en, rho_en
  REAL                        :: v_crit

  INTEGER, DIMENSION(10)      :: fit_array
  CHARACTER (LEN=25), DIMENSION(npmax)  :: aa_txt
  INTEGER, DIMENSION(cont)    :: eqa              ! flag, indicating, whether a data type is considered

  CHARACTER (LEN=50)          :: pure_fit_file

End Module pure_fit_parameters


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PURE_RESIDUAL
!
! Routine for pure component parameter fitting. The residual between 
! calculated values and experimental data is calculated and stored in
! vector 'fvec'. The vector 'fvec' serves as the objective function for
! the parameter.
! PURE_RESIDUAL is called by FITTING via the Levenberg-Marquardt routine.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PURE_RESIDUAL (lm_m_dat, lm_n_par, para, fvec, iflag)

  USE PARAMETERS, ONLY: RGAS
  USE BASIC_VARIABLES
  USE pure_fit_parameters
  USE utilities, only: SI_DENS, dens_calc, p_calc, error_message
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! INTEGER, PARAMETER                     :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, INTENT(IN)                    :: lm_m_dat
  INTEGER, INTENT(IN)                    :: lm_n_par
  REAL, INTENT(IN)                       :: para(:)
  REAL, INTENT(IN OUT)                   :: fvec(:)
  INTEGER, INTENT(IN OUT)                :: iflag

  !-----------------------------------------------------------------------------
  INTEGER                                :: converg, crit_dat
  INTEGER                                :: i, j, ps
  INTEGER                                :: j_dat = 0
  LOGICAL                                :: scan,vapor,nocon
  REAL                                   :: weigh, tc, v_dum, plv_dum
  REAL                                   :: pges,zges
  REAL                                   :: density(np), w(np,nc)
  REAL                                   :: rho_phas(np)
  CHARACTER (LEN=4)                      :: char_len
  ! REAL                                   :: tclit,pclit,pc,rhoc
  !-----------------------------------------------------------------------------

  crit_dat = 0
  tc = 0.0

  IF ( (nliq >= 1) .OR. (nlv >= 1) .OR. (n_vir >= 1) .OR. (n_en >= 1) ) THEN
     DO i=1,11
        ps = 0
        DO j = 1, lm_n_par
           IF (fit_array(j) == i) ps=j
        END DO
        IF (i == 1 .AND. ps >= 1) parame(1,1) = mm(1)*para(ps)*rel(ps)
        IF (i == 2 .AND. ps >= 1) parame(1,2) = para(ps)*rel(ps)
        IF (i == 3 .AND. ps >= 1) parame(1,3) = para(ps)*rel(ps)
        IF (i == 5 .AND. ps >= 1) parame(1,13)= para(ps)*rel(ps)
        IF (i == 6 .AND. ps >= 1) parame(1,6) = para(ps)*rel(ps)
        IF (i == 7 .AND. ps >= 1) parame(1,8) = para(ps)*rel(ps)
        IF (i == 8 .AND. ps >= 1) parame(1,7) = para(ps)*rel(ps)
        IF (i == 9 .AND. ps >= 1) parame(1,9) = para(ps)*rel(ps)
        IF (i == 11.AND. ps >= 1) parame(1,11)= para(ps)*rel(ps)
        IF (i == 4 .AND. ps >= 1) THEN
           IF ( NINT(parame(1,12)) == 2 ) THEN
              parame(1,14) = 0.0                  ! eps_hb(1,1,1,1)
              parame(1,15) = para(ps)*rel(ps)     ! eps_hb(1,1,1,2)
              parame(1,16) = para(ps)*rel(ps)     ! eps_hb(1,1,2,1)
              parame(1,17) = 0.0                  ! eps_hb(1,1,2,2)
           ELSE IF ( NINT(parame(1,12)) == 1 ) THEN
              parame(1,14) = para(ps)*rel(ps)     ! eps_hb(1,1,1,1)
           ELSE
              call error_message('extend pure_par_fit for this association case')
           END IF
        END IF
     END DO

     ! IF ( NINT(parame(1,12)) == 0) THEN        ! no association
     ! ELSEIF ( NINT(parame(1,12)) == 1) THEN    ! No of assoc sites = 1
     !    parame(1,15)=1.0
     ! ELSEIF ( NINT(parame(1,12)) == 2) THEN    ! No of assoc sites = 2
     !    parame(1,18)=1.0
     !    parame(1,19)=1.0
     ! ELSEIF ( NINT(parame(1,12)) == 3) THEN    ! No of assoc sites = 3
     !    parame(1,18)=2.0
     !    parame(1,19)=1.0
     ! ELSEIF ( NINT(parame(1,12)) == 4) THEN    ! No of assoc sites = 4
     !    parame(1,18)=2.0
     !    parame(1,19)=2.0
     ! ELSE
     !    write (*,*) ' define your association case'
     !    stop
     ! ENDIF

     parame(1,4) = 0.0
     IF (eos == 0) parame(1,4) = 10.0

  END IF


  WRITE (char_len,'(I3)') lm_n_par
  IF (iflag == 1) WRITE (*,'(a,'//char_len//'(G15.8))') ' parameter ',(para(i)*rel(i),i=1,lm_n_par)


  IF ( MINVAL( parame(1,1:25) ) < 0.0 ) THEN
     fvec(:) = fvec(:) * 2.0 ! penalty function for negative parameters
     write (*,*) 'warning negative pure component parameters encountered',parame(1,MINLOC(parame(1,1:25)))
     RETURN
  ENDIF


  fvec = 0.0


  !-----------------------------------------------------------------------------
  ! liquid densit data
  !-----------------------------------------------------------------------------
  DO i = 1, nliq

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
        IF ( vliq(i) <= v_crit ) THEN
           vapor = .false.
           densta(1) = 0.5
           IF (d_kond(1,i) /= 0.0) densta(1) = d_kond(1,i)
        ELSE
           vapor = .true.
           densta(1) = 0.002
           IF (d_kond(1,i) /= 0.0) densta(1) = d_kond(1,i)
        END IF

        vapor_scan: DO
           IF (scan) densta(1) = densta(1) + 0.002

           CALL dens_calc ( rho_phas )
           CALL SI_dens ( density, w )
           v_cal(i) = 1.0 / density(1)

           IF ((v_cal(i)-vliq(i))*1.E2/vliq(i) < 5.0 .AND. .NOT. scan) THEN
              d_kond(1,i) = dense(1)
           ELSE
              d_kond(1,i) = 0.0
           END IF


           IF ( vapor .AND. (1.0/density(1)) <= (v_crit/1.1) ) THEN
              IF ( densta(1) <= 0.228 ) THEN
                 scan = .true.
                 CYCLE vapor_scan
              ELSE
                 nocon = .true.
              END IF
           END IF
           EXIT vapor_scan
        END DO vapor_scan
        IF (dense(1) > 0.55) WRITE (*,*)'density ',i,dense(1)

        weigh = 0.5
        IF (nocon) weigh = weigh / 4.0

        fvec(j_dat) = weigh * SQRT(type_weight(1))*(v_cal(i)-vliq(i))/  &
             (SQRT(ABS(vliq(i)))*SQRT(ABS(v_cal(i))))



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
        IF (d_kond(1,j_dat) /= 0.0) val_init(1) = d_kond(1,j_dat)
        IF (d_kond(2,j_dat) /= 0.0) val_init(2) = d_kond(2,j_dat)
        IF (d_kond(1,j_dat) /= 0.0 .AND. d_kond(2,j_dat) /= 0.0) val_init(4) = plv_kon(j_dat)
        val_init(5) = 0.0      ! mole fraction lnx(1,1)
        val_init(6) = 0.0      ! mole fraction lnx(2,1)

        CALL pure_equilibrium_fit ( crit_dat, tc, tliq(i), pliq(i), converg, v_cal(i), plv_dum )

        IF ( converg == 1 ) THEN
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

        fvec(j_dat) = weigh * SQRT(type_weight(1))*(v_cal(i)-vliq(i))  &
             / (SQRT(ABS(vliq(i)))*SQRT(ABS(v_cal(i))))

     END IF

  END DO
  !-----------------------------------------------------------------------------
  ! end: liquid density data
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  ! vapor pressure data
  !-----------------------------------------------------------------------------
  DO  i = 1, nlv

     nphas = 2

     j_dat = i + nliq

     val_init(1) = 0.5
     val_init(2) = 1.E-8
     val_init(3) = tlv(i)
     val_init(4) = plv(i)
     IF (d_kond(1,j_dat) /= 0.0) val_init(1) = d_kond(1,j_dat)
     IF (d_kond(2,j_dat) /= 0.0) val_init(2) = d_kond(2,j_dat)
     IF (d_kond(1,j_dat) /= 0.0 .AND. d_kond(2,j_dat) /= 0.0) val_init(4) = plv_kon(j_dat)
     val_init(5) = 0.0      ! mole fraction lnx(1,1)
     val_init(6) = 0.0      ! mole fraction lnx(2,1)

     CALL pure_equilibrium_fit ( crit_dat, tc, tlv(i), plv(i), converg, v_dum, plvcal(i) )

     IF ( converg == 1 ) THEN
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
           WRITE (*,'(a,i3,2G14.5,G15.8)') ' loose ! ',i,tlv(i),t,plv(i)
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

     fvec(j_dat) = 0.2*SQRT(type_weight(2))*ABS(plvcal(i)-plv(i))  &
          / ((SQRT(ABS(plv(i)))*SQRT(ABS(plvcal(i)))))**0.85

     IF (tc < tlv(i) .AND. crit_dat == 1) THEN
        fvec(j_dat)=fvec(j_dat) + 5.0*(tlv(i)-tc)
        WRITE (*,'(a,i3,3(f15.5))') ' skipped ! ',i,tlv(i),tc
     END IF

  END DO
  !-----------------------------------------------------------------------------
  ! end: vapor pressure data
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! 2nd virial coefficients
  !-----------------------------------------------------------------------------

  DO  i = 1, n_vir

     j_dat = i + nliq + nlv

     nphas=1
     dense(1)= 1.E-7
     t       = t_vir(i)
     p       = 1.0
     xi(1,1)  = 1.0

     CALL p_calc (pges,zges)
     CALL si_dens (density,w)

     b_cal(i) = (zges-1.0)/density(1)*1000.0*mm(1)

     weigh = (parame(1,3)/t)**3
     fvec(j_dat) = type_weight(3)*weigh*(b_cal(i)-b_vir(i))  &
          / (SQRT(ABS(b_vir(i)))*SQRT(ABS(b_cal(i))))

  END DO



  !-----------------------------------------------------------------------------
  ! internal energy data
  !-----------------------------------------------------------------------------

  DO i = 1, n_en

     j_dat = i + nliq + nlv + n_vir

     nphas = 1
     dense(1) = rho_en(i)
     t        = t_en(i)
     xi(1,1)  = 1.0

     ! CALL P_CALC (pges,zges)
     ! CALL F_EOS (...)

     ! u_res = f_res + T*S

     weigh = 1.0
     fvec(j_dat) = type_weight(4)*weigh*(u_cal(i)-u_en(i))/ ABS(u_en(i))

  END DO

  IF (j_dat /= lm_m_dat) write (*,*) 'array length not matching!'


  !-----------------------------------------------------------------------------
  ! critical point data (work in progress)
  !-----------------------------------------------------------------------------

  !-----------critical point---------------------------------------------
  ! naphtalene
  !      tclit   = 748.210
  !      pclit   = 4100000.0
  !      rhoclit = 246.0
  ! co2
  !      tclit   = 304.210
  !      pclit   = 7382500.0
  !      rhoclit = 466.01015160
  ! propane
  !      tclit   = 370.00
  !      pclit   = 4.265d+06
  !      rhoclit = 225.00
  ! butane
  !      tclit   = 425.20
  !      pclit   = 3.796d+06
  !      rhoclit = 227.00
  !      PC-SAFT: T=432.1 P=41.95
  ! hexane
  !       tclit   = 507.40
  !       pclit   = 29.7d5
  !       rhoclit = 233.50
  ! heptane
  !      tclit   = 539.70
  !      pclit   = 2.736d+06
  !      rhoclit = 234.10
  ! octane
  !      tclit   = 569.40
  !      pclit   = 2.496d+06
  !      rhoclit = 234.70
  ! acetone
  !      tclit   = 508.15
  !      pclit   = 4.7621E+06
  !      rhoclit = 234.70

  ! tc   = tclit
  ! pc   = pclit
  ! rhoc = rhoclit

  ! CALL CRITICAL (tc,pc,rhoc)

  ! write (*,'(t5,a,f7.2,a,f7.3,a,f7.2)') 'tc =',tc,'  pc =',pc/1.E5,'  rho =',rhoc
  ! write (*,'(t5,a,f7.2,a,f7.3,a,f7.2)') 'tc =',tclit,'  pc =',pclit/1.E5,'  rho =',rhoclit

  !c   j_dat = nliq + nlv + n_vir + n_en + 1
  !c   fvec(j_dat)=0.020*SQRT(ABS(rhoc/rhoclit-1.0))
  ! j_dat = nliq + nlv + n_vir + n_en + 2
  ! fvec(j_dat) = fvec(j_dat) + 1.0 * SQRT(ABS(tc/tclit-1.0))
  ! j_dat = nliq + nlv + n_vir + n_en + 3
  ! fvec(j_dat)=1.0*SQRT(ABS(pc/pclit-1.0))

  !-----------critical point---------------------------------------------

  ! DO i = 1, j_dat
  !   write (*,*) i,fvec(i)
  ! ENDDO
  ! call pause

END SUBROUTINE PURE_RESIDUAL



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pure_equilibrium_fit
!
! input to the subroutine: val_init has to be available
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE pure_equilibrium_fit (crit_dat, tc, t_in, p_in, converg, v_cal, plvcal)

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: densav, eta_start, pges
  use utilities, only: SI_DENS, CRITICAL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN OUT)                :: crit_dat
  REAL, INTENT(IN OUT)                   :: tc
  REAL, INTENT(IN)                       :: t_in
  REAL, INTENT(IN)                       :: p_in
  INTEGER, INTENT(OUT)                   :: converg
  REAL, INTENT(OUT)                      :: v_cal
  REAL, INTENT(OUT)                      :: plvcal
  !-----------------------------------------------------------------------------
  REAL                                   :: temp, press
  REAL                                   :: pc, rhoc
  REAL                                   :: density(np), w(np,nc)
  REAL                                   :: p_high, p_low
  !-----------------------------------------------------------------------------


  n_unkw = ncomp                         ! number of quantities to be iterated
  it(1) = 'lnp'                          ! iteration of pressure
  val_conv(2) = 0.0
  densav(:) = 0.0

  temp  = t_in
  press = p_in

  plvcal = 0.0
  v_cal  = 0.0


  !-----------------------------------------------------------------------------
  ! phase equilibrium calculation
  !-----------------------------------------------------------------------------

  CALL objective_ctrl (converg)


  !-----------------------------------------------------------------------------
  ! for the case, where no convergence was found initially
  !-----------------------------------------------------------------------------

  ! --- calc critical pt. in order to assess, whether a solution is possible
  IF ( converg /= 1 .AND. crit_dat == 0 ) THEN
     tc = temp
     pc = press
     CALL critical ( tc, pc, rhoc )
     IF ( tc < 0.0 ) write (*,*) 'error: negative Tc calculated'
     IF ( tc < 0.0 ) tc = 500.0
     crit_dat = 1
  END IF

  !-----------------------------------------------------------------------------
  ! 1st step to achieve convergence: select feasible p (betw. spinodal p's)
  !-----------------------------------------------------------------------------
  IF ( converg /= 1 .AND. tc > temp ) THEN
     t = temp
     CALL PERTURBATION_PARAMETER
     eta_start = 1.E-8
     CALL PRESSURE_SPINODAL              ! vapor spinodal: max p for which a vapor exists
     p_high = pges
     eta_start = 0.45
     CALL PRESSURE_SPINODAL              ! liquid spinodal: min p for which a liquid exists
     p_low = pges
     val_init(1) = 0.5
     val_init(2) = 1.E-8
     val_init(3) = temp
     val_init(4) = 0.8 * p_high  +  0.2 * MAX( p_low, 0.0 )   ! choose p betw. the 2 spinodals
     CALL objective_ctrl (converg)
  ENDIF

  !-----------------------------------------------------------------------------
  ! 2nd step to achieve convergence: approach t from lower t in a few steps
  !-----------------------------------------------------------------------------
  IF ( converg /= 1 .AND. tc > temp ) THEN
     val_init(1) = 0.5
     val_init(2) = 1.E-6
     val_init(3) = temp - 0.1*temp
     val_init(4) = press / 10.0
     IF (val_init(4) == 0.0) val_init(4) = EXP( ( 1.0 - tc/val_init(3) )*7.0 ) * pc ! simple correspondence principle
     outp = 0           ! output to terminal (set u_out_p before turning output on)
     running = 't'      ! Temperature is running var. in PHASE_EQUILIB
     CALL phase_equilib ( temp, 4.0, converg )
     IF ( converg == 1 ) write (*,'(a,G14.5,G15.8)') ' win ! ',temp, p
  END IF

  IF ( converg == 1 ) THEN
     plvcal = val_conv(4)
     CALL SI_DENS (density,w)
     v_cal = 1.0 / density(1)
  END IF

END SUBROUTINE pure_equilibrium_fit





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE fitting
!
! Driver for the fitting of pure component parameters. This subroutine
! uses LMDIF1 (Levenberg-Marquardt MINPACK routine) in order to minimze
! the deviation of calculation results to experimental values. The
! deviations are calculated in PURE_RESIDUAL and stored in 'fvec'.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE fitting

  USE BASIC_VARIABLES
  USE Levenberg_Marquardt
  USE pure_fit_parameters, ONLY: pure_fit_file
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE PURE_RESIDUAL(lm_m_dat, lm_n_par, para, fvec, iflag)
       IMPLICIT NONE
       INTEGER, INTENT(IN)                 :: lm_m_dat, lm_n_par
       REAL, INTENT(IN)                    :: para(:)
       REAL, INTENT(IN OUT)                :: fvec(:)
       INTEGER, INTENT(IN OUT)             :: iflag
     END SUBROUTINE PURE_RESIDUAL

     SUBROUTINE paini( lm_n_par, para )
       IMPLICIT NONE
       INTEGER, INTENT(OUT)                 :: lm_n_par
       REAL, ALLOCATABLE, INTENT(OUT)       :: para(:)
     END SUBROUTINE paini

     SUBROUTINE pure_output( lm_n_par, para )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                 :: lm_n_par
       REAL, INTENT(IN)                    :: para(:)
     END SUBROUTINE pure_output
  END INTERFACE

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


  ncomp     = 1

  ! --- numerical constants for Levenberg-Marquardt algorithm ------------
  scaling(1) = 1.E2
  IF ( num == 0 ) epsfcn = 1.0E-6**2  ! sqrt of relat. step size (finite differences)
  IF ( num == 1 ) epsfcn = 1.0E-5**2  ! sqrt of relat. step size (finite differences)
  lm_factor = 0.1                     ! maximum initial step size (parameter iteration)

  IF ( num == 2 ) write (*,*) 'fitting with RGT is not yet supported.  Set num <= 1'
  IF ( num == 2 ) stop

  !-----------------------------------------------------------------------------
  ! read input-file (DATIN), read parameters (PAINI), prepare output
  !-----------------------------------------------------------------------------

  WRITE (*,*) '  SPECIFY INPUT FILE: ( ./input_file/pure_comp/ )'
  READ (*,*) pure_fit_file
  ! pure_fit_file = '3939.dat'        !JG

  CALL file_open('./input_file/pure_comp/'//pure_fit_file,22)
  OPEN (23,FILE='./output_file/pure_comp/'//pure_fit_file)
  pure_fit_file='./input_file/pure_comp/'//pure_fit_file

  CALL datin( lm_m_dat )
  CALL paini( lm_n_par, para )

  ALLOCATE( fvec(lm_m_dat) )
  ALLOCATE( ipvt(lm_n_par) )


  !-----------------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  !-----------------------------------------------------------------------------

  CALL lmdif1 (PURE_RESIDUAL,lm_m_dat,lm_n_par,para,fvec,tol,epsfcn,lm_factor,info,ipvt)


  !-----------------------------------------------------------------------------
  ! algorithm output: optimization status
  !-----------------------------------------------------------------------------

  WRITE (*,*)  ' '
  WRITE (23,*) ' '

  SELECT CASE (info)
  CASE (:-1)
     WRITE (23,*) 'SOLVER STATUS: Users FCN returned INFO = ', -info
     WRITE (*, *) 'SOLVER STATUS: Users FCN returned INFO = ', -info
  CASE (0)
     WRITE (23,*) 'SOLVER STATUS: Improper values for input parameters'
     WRITE (*, *) 'SOLVER STATUS: Improper values for input parameters'
  CASE (1:3)
     WRITE (23,*) 'SOLVER STATUS: Convergence'
     WRITE (*, *) 'SOLVER STATUS: Convergence'
  CASE (4)
     WRITE (23,*) 'SOLVER STATUS: Residuals orthogonal to the Jacobian'
     WRITE (*, *) 'SOLVER STATUS: Residuals orthogonal to the Jacobian'
     WRITE (*, *) 'There may be an error in FCN'
  CASE (5)
     WRITE (23,*) 'SOLVER STATUS: Too many calls to FCN'
     WRITE (*, *) 'SOLVER STATUS: Too many calls to FCN'
     WRITE (*, *) 'Either slow convergence, or an error in FCN'
  CASE (6:7)
     WRITE (23,*) 'SOLVER STATUS: TOL was set too small'
     WRITE (*, *) 'SOLVER STATUS: TOL was set too small'
  CASE DEFAULT
     WRITE (23,*) 'SOLVER STATUS: INFO =', info, ' ???'
     WRITE (*, *) 'SOLVER STATUS: INFO =', info, ' ???'
  END SELECT

  WRITE (*,*)  '-------------'
  WRITE (23,*) '-------------'

  CALL pure_output ( lm_n_par, para )

  CLOSE (22)
  CLOSE (23)

  DEALLOCATE( fvec, ipvt )

END SUBROUTINE fitting



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE DATIN
!
! Read experimental pure component data for the regression of pure
! component parameters. This file also writes the experimental data to
! the output file.
! So far, the experimental data is stored in vectors of fixed length
! (see module 'pure_fit_parameters'). At some point this should be
! changed - allocate the vector length here.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE datin ( lm_m_dat )

  USE BASIC_VARIABLES
  USE pure_fit_parameters
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(OUT)                   :: lm_m_dat

  !-----------------------------------------------------------------------------
  INTEGER                                :: mak, i, k
  INTEGER                                :: datatype_available(cont)
  INTEGER                                :: data_no, data_low, data_up, alldata
  REAL                                   :: mmpliq
  CHARACTER (LEN=80)                     :: info_line
  !-----------------------------------------------------------------------------

  nlv   = 0
  nliq  = 0
  n_vir = 0
  n_en  = 0

  eqa = 0

  !-----------------------------------------------------------------------------
  ! reading info-lines (input-file) and writing them (output-file)
  !-----------------------------------------------------------------------------

  DO  i = 1, 15
     READ(22,'(A80)') info_line
     WRITE(23,'(a)') info_line
  END DO

  ! --- reading molecular mass and an estimate of the crit. density ------------
  READ(22,'(A80)') info_line
  WRITE(23,'(A)') info_line
  READ(22,*) mmpliq, v_crit
  WRITE(23, '(2(2x,G13.6))') mmpliq, v_crit
  mm(1) = mmpliq

  DO i=1,14
     READ(22,'(A80)') info_line
     WRITE(23,'(A)') info_line
  END DO

  ! --- reading the number of experimental data points -------------------------
  READ(22,'(A80)') info_line
  READ(22,*) (datatype_available(i),i=1,cont)
  WRITE (*,*) ' '

  eos = 1

  !-----------------------------------------------------------------------------
  ! selection of data sets by user
  !-----------------------------------------------------------------------------

  WRITE (*,*), ' DATA SELECTION'
  WRITE (*,*) ' '
  WRITE (*,*), ' NUMBER OF DATA-SETS:'
  WRITE (*,'(T2,A40,I3)') 'PVT - data                :', datatype_available(1)
  WRITE (*,'(T2,A40,I3)') 'vapor pressure data       :', datatype_available(2)
  WRITE (*,'(T2,A40,I3)') '2nd virial coefficients   :', datatype_available(3)
  WRITE (*,'(T2,A40,I3)') 'internal res. energy data :', datatype_available(4)
  WRITE (*,*) ' '

  ! --- option to select all data sets -----------------------------------------
  WRITE (*,*), ' SELECT DATA-SET !'
  WRITE (*,*), '      consider all available data:     1'
  WRITE (*,*), '      select data-sets individually:   0'
  READ (*,*) alldata
  ! alldata = 1        !JG
  IF ( alldata == 0 ) WRITE (*,*), ' SELECT DATA-TYPES !'
  WRITE (*,*) ' '
  IF (datatype_available(1) /= 0) THEN
     eqa(1)=1
     IF (alldata == 0) WRITE (*,*), ' PVT-data ?  (0/1) '
     IF (alldata == 0) READ (*,*) eqa(1)
  END IF
  IF (datatype_available(2) /= 0) THEN
     eqa(2)=1
     IF (alldata == 0) WRITE (*,*), ' vapor pressure data ? (0/1)'
     IF (alldata == 0) READ (*,*) eqa(2)
  END IF
  IF (datatype_available(3) /= 0) THEN
     eqa(3)=1
     IF (alldata == 0) WRITE (*,*), ' 2nd virial coefficients ? (0/1)'
     IF (alldata == 0) READ (*,*) eqa(3)
  END IF
  IF (datatype_available(4) /= 0) THEN
     eqa(4)=1
     IF (alldata == 0) WRITE (*,*), ' internal res. energy data? (0/1)'
     IF (alldata == 0) READ (*,*) eqa(4)
  END IF

  WRITE (*,*) ' '

  !-----------------------------------------------------------------------------
  ! reading data sets from input file
  !-----------------------------------------------------------------------------

  DO  i = 1, cont   ! --- cont: number of data types (p^sat, rho_L, ...)

     ! --- two headers ---------------------------------------------------------
     IF (datatype_available(i) /= 0) READ(22,'(A80)') info_line
     IF (eqa(i) == 1)  WRITE(23,'(a)') info_line
     IF (datatype_available(i) /= 0) READ(22,'(A80)') info_line
     IF (eqa(i) == 1)  WRITE(23,'(a)') info_line

     ! --- reading data sets if a certain data type i was chosen ---------------
     WRITE (*,*) ' '
     IF (eqa(i) == 1) THEN
        IF (i == 1) WRITE (*,*), ' CHOOSE FROM PVT-DATA:'
        IF (i == 2) WRITE (*,*), ' CHOOSE FROM VAPOR PRESSURE DATA:'
        IF (i == 3) WRITE (*,*), ' CHOOSE FROM 2nd VIRIAL COEFFICIENT DATA:'
        IF (i == 4) WRITE (*,*), ' CHOOSE FROM INTERNAL RESID. ENERGY DATA:'

        data_low = 1
        data_up = datatype_available(i)

        ! --- user selects the considered data sets ----------------------------
        WRITE (*,'(T2,A,I3)') ' Total number of data-sets: ',datatype_available(i)
        WRITE (*,*), ' Specify the data-sets to be considered'
        WRITE (*,*), ' Lower number of data-set (e.g.  1):'
        IF (alldata == 0) READ (*,*) data_low
        WRITE (*,'(a,I3,a)') ' Upper number of data-set (e.g.',datatype_available(i),'):'
        IF (alldata == 0) READ (*,*) data_up

        ! --- check the input --------------------------------------------------
        IF ((data_low > data_up) .OR. (data_up > datatype_available(i))) THEN
           WRITE (*,*) ' Erroneous input! The lower data-set number must be smaller'
           WRITE (*,*) ' than the upper value.'
           WRITE (*,*) ' The upper data-set number must not be greater than ',datatype_available(i)
           WRITE (*,*) ' Lower number of data-set (e.g.  1):'
           READ (*,*) data_low
           WRITE (*,'(a,I3,a)') ' Upper number of data-set (e.g.',datatype_available(i),'):'
           READ (*,*) data_up
        END IF

        mak = 0

        ! --- reading data -----------------------------------------------------
        do k = 1, data_low-1
           read (22,*)
        end do
        IF (i == 1) THEN
           nliq = 0
           DO k = data_low, data_up
              nliq = nliq + 1
              READ (22,*) data_no, pliq(nliq), tliq(nliq), vliq(nliq)
              WRITE (*,'(i3,3(2x,G12.5))') data_no, pliq(nliq), tliq(nliq), vliq(nliq)
              WRITE(23,'(i3,3(2x,G12.5))') data_no, pliq(nliq), tliq(nliq), vliq(nliq)
           END DO
        ELSE IF (i == 2) THEN
           nlv = 0
           DO k = data_low, data_up
              nlv = nlv + 1
              READ (22,*) data_no, plv(nlv), tlv(nlv)
              WRITE (*,'(i3,2(2x,G12.5))') data_no, plv(nlv), tlv(nlv)
              WRITE(23,'(i3,2(2x,G12.5))') data_no, plv(nlv), tlv(nlv)
           END DO
        ELSE IF (i == 3) THEN
           n_vir = 0
           DO k = data_low, data_up
              n_vir = n_vir + 1
              READ (22,*) data_no, b_vir(n_vir), t_vir(n_vir)
              WRITE (*,'(i3,2(2x,G12.5))') data_no, b_vir(n_vir), t_vir(n_vir)
              WRITE(23,'(i3,2(2x,G12.5))') data_no, b_vir(n_vir), t_vir(n_vir)
           END DO
        ELSE IF (i == 4) THEN
           n_en = 0
           DO k = data_low, data_up
              n_en = n_en + 1
              READ (22,*) data_no, u_en(n_en), rho_en(n_en), t_en(n_en)
              WRITE (*,'(i3,3(2x,G12.5))') data_no, u_en(n_en), rho_en(n_en), t_en(n_en)
              WRITE(23,'(i3,3(2x,G12.5))') data_no, u_en(n_en), rho_en(n_en), t_en(n_en)
           END DO
        END IF
        do k = data_up+1, datatype_available(i)
           read (22,*)
        end do

     ELSE

        ! --- a certain data type was not selected for the fitting -------------
        IF (datatype_available(i) /= 0) THEN
           DO  k = 1, datatype_available(i)
              READ(22,*) data_no
           END DO
        END IF

     END IF

  END DO

  !-----------------------------------------------------------------------------
  ! total number of data sets considered in parameter fitting
  !-----------------------------------------------------------------------------

  lm_m_dat = nliq + nlv + n_vir + n_en       ! excluding crit. point
  ! lm_m_dat = lm_m_dat + 3                  ! including crit. point

  REWIND (22)

END SUBROUTINE datin



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE paini
!
! Subroutine for reading in starting values of the pure component
! parameters
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE paini ( lm_n_par, para )

  USE BASIC_VARIABLES
  USE pure_fit_parameters
  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(OUT)                   :: lm_n_par
  REAL, ALLOCATABLE, INTENT(OUT)         :: para(:)

  !-----------------------------------------------------------------------------
  INTEGER                                :: fitopt, POSITION, i_fit
  INTEGER                                :: n_sites
  INTEGER                                :: i, j, assoc, quadrup, dipole
  REAL                                   :: dummy( npmax, 3 )
  INTEGER                                :: k, no, nhb_typ(nc)
  REAL                                   :: nhb_no(nc,nsite)
  REAL                                   :: eps_hb(nc,nc,nsite,nsite)
  !-----------------------------------------------------------------------------

  type_weight = 1.0

  assoc  = 0
  dipole = 0
  quadrup= 0

  nhb_typ(:) = 0
  nhb_no(:,:) = 0.0

  !-----------------------------------------------------------------------------
  ! selecting adjustable parameters
  !-----------------------------------------------------------------------------

  WRITE (*,*) ' '
  WRITE (*,*), '*********** PARAMETER SPECIFICATION ***********'
  WRITE (*,*) ' '

  WRITE (*,*) ' You may choose from two default fitting options'
  WRITE (*,*) ' or individually specify all parameters.'
  WRITE (*,*) ' (Dipolar and quadrupolar moments are set to zero'
  WRITE (*,*) ' for both default options; for dipolar or'
  WRITE (*,*) ' quadrupolar compounds choose option 3)'
  WRITE (*,*) '  '
  WRITE (*,*) ' Choose one the following options'
  WRITE (*,*) ' Default 1: Fitting 3 parameters'
  WRITE (*,*) '            * segment number           m        [ ]'
  WRITE (*,*) '            * segment diameter         sigma    [A]'
  WRITE (*,*) '            * segment energy parameter eps/k    [K]'
  WRITE (*,*) ' Default 2: Fitting 5 parameters for associating'
  WRITE (*,*) '            compounds with 2 association sites'
  WRITE (*,*) '            * segment number           m        [ ]'
  WRITE (*,*) '            * segment diameter         sigma    [A]'
  WRITE (*,*) '            * segment energy parameter eps/k    [K]'
  WRITE (*,*) '            * assoc.-energy parameter  eps_AB/k [K]'
  WRITE (*,*) '            * assoc.-volume parameter  kappa_AB [ ]'
  WRITE (*,*) ' Option  3: Individual specification of all'
  WRITE (*,*) '            parameters to be fitted.'
  WRITE (*,*) ' Choose 1, 2, or 3.'
  READ (*,*) fitopt
  ! fitopt = 3         !JG

  IF (fitopt == 1) THEN
     lm_n_par = 3
     !parame(1,12) = 0.0                         ! number of assoc. sites (real-type !)
     nhb_typ(1)  = 0      ! 11.03.2015
     nhb_no(1,1) = 0.0    ! 11.03.2015
     aa_txt(1) = ' Segm.No.m / Molar Mass :'
     aa_txt(2) = ' Segm.Diameter sigma    :'
     aa_txt(3) = ' Segm.Energy Param eps/k:'
     fit_array(1) = 1
     fit_array(2) = 2
     fit_array(3) = 3
     DO i = 4, 10
        fit_array(i) = 0
     ENDDO
  ELSEIF (fitopt == 2) THEN
     lm_n_par = 5
     !parame(1,12) = 2.0                         ! number of assoc. sites (real-type !)
     aa_txt(1) = ' Segm.No.m / Molar Mass :'
     aa_txt(2) = ' Segm.Diameter sigma    :'
     aa_txt(3) = ' Segm.Energy Param eps/k:'
     aa_txt(4) = ' Assoc.-energy eps_AB/k :'
     aa_txt(5) = ' Assoc.-volume kappa_AB :'
     fit_array(1) = 1
     fit_array(2) = 2
     fit_array(3) = 3
     fit_array(4) = 4
     fit_array(5) = 5
     DO i = 6, 10
        fit_array(i)=0
     ENDDO
     nhb_typ(1)  = 2
     nhb_no(1,1) = 1.0
     nhb_no(1,2) = 1.0
  ELSE
     WRITE (*,*) ' The following pure component parameters can'
     WRITE (*,*) ' be fitted'
     WRITE (*,*) ' (1) segment number                  m        [ ]'
     WRITE (*,*) ' (2) segment diameter                sigma    [A]'
     WRITE (*,*) ' (3) segment energy parameter        eps/k    [K]'
     WRITE (*,*) ' (4) assoc.-energy parameter         eps_AB/k [K]'
     WRITE (*,*) ' (5) assoc.-volume parameter         kappa_AB [ ]'
     WRITE (*,*) ' (6) dipolar moment                  my       [D]'
     WRITE (*,*) ' (7) fraction of segm.with dipolar moment xp  [ ]'
     WRITE (*,*) ' (8) quadrupolar moment              Q        [D]'
     WRITE (*,*) ' (9) fraction of segm.with quadrup.moment xp  [ ]'
     WRITE (*,*) ' (11)polarizability alpha                   [A^3]'
     WRITE (*,*) ' '
     WRITE (*,*) ' Note, that the number of association sites has'
     WRITE (*,*) ' to be user-specified'
     WRITE (*,*) ' '
     WRITE (*,*) ' Choose the total number of fitting-parameters'
     READ (*,*) lm_n_par
     ! lm_n_par = 3         !JG
     WRITE (*,*) ' Choose the position-number of fitting-param.'
     WRITE (*,*) ' (see above list of parameter 1 to 9)'
     DO i = 1, lm_n_par
        WRITE (*,'(a,i1,a,i2)') '  position-number for param. ',i,' of ',lm_n_par
        READ (*,*) position
        IF (position == 1) aa_txt(i) = ' Segm.No.m / Molar Mass :'
        IF (position == 2) aa_txt(i) = ' Segm.Diameter sigma    :'
        IF (position == 3) aa_txt(i) = ' Segm.Energy Param eps/k:'
        IF (position == 4 .OR. position == 5) assoc = 1
        IF (position == 4) aa_txt(i) = ' Assoc.-energy eps_AB/k :'
        IF (position == 5) aa_txt(i) = ' Assoc.-volume kappa_AB :'
        IF (position == 6 .OR. position == 7) dipole = 1
        IF (position == 6) aa_txt(i) = ' Dipolar moment         :'
        IF (position == 7) aa_txt(i) = ' No. of segm.with dipole:'
        IF (position == 8 .OR. position == 9) quadrup = 1
        IF (position == 8) aa_txt(i) = ' Quadrupolar moment     :'
        IF (position == 9) aa_txt(i) = ' No.of segm.with quadrup:'
        IF (position == 11)aa_txt(i) = ' polarizability alpha'
        fit_array(i) = position
     ENDDO
     IF ( assoc == 1 ) THEN
        WRITE (*,*) ' Give number of association sites'
        READ (*,*) n_sites
        IF (n_sites == 1) nhb_typ(1) = 1
        IF (n_sites == 1) nhb_no(1,1)= 1.0
        IF (n_sites == 2) nhb_typ(1) = 2
        IF (n_sites == 2) nhb_no(1,1)= 1.0
        IF (n_sites == 2) nhb_no(1,2)= 1.0
        IF (n_sites == 3) nhb_typ(1) = 2
        IF (n_sites == 3) nhb_no(1,1)= 1.0
        IF (n_sites == 3) nhb_no(1,2)= 2.0
        IF (n_sites == 4) nhb_typ(1) = 2
        IF (n_sites == 4) nhb_no(1,1)= 2.0
        IF (n_sites == 4) nhb_no(1,2)= 2.0
        ! parame(1,12) = REAL( nhb_typ(1) )

     ENDIF
  ENDIF

  ALLOCATE( para(lm_n_par) )


  !-----------------------------------------------------------------------------
  ! specifying starting values of adjustable parameters
  !-----------------------------------------------------------------------------

  WRITE (*,*) ' '
  WRITE (*,*), ' Input of starting values for fitting-parameters:'
  WRITE (*,*) ' '

  WRITE (*,*), ' GIVE STARTING VALUES FOR PARAMETERS  '
  WRITE (*,*) ' '

  DO i = 1, lm_n_par
     WRITE (*,'(2a,i2,a,i2,a)') aa_txt(i),' (param.',i,' of',lm_n_par,')'
     READ (*,*) (dummy(i,j), j = 1,1)
     dummy(i,2) = dummy(i,1) / 10.0
     dummy(i,3) = dummy(i,1) / dummy(i,2)
  END DO
  WRITE (*,*) ' '

  !-----------------------------------------------------------------------------
  ! for case, where param. are selected individually, determine rest
  !-----------------------------------------------------------------------------

  n_sites = NINT( SUM( nhb_no( 1, 1:nhb_typ(1) ) ) )
  IF ( fitopt >= 3 ) THEN
     WRITE (*,*) ' '
     WRITE (*,*) ' The remaining parameters (not-fitted) now'
     WRITE (*,*) ' have to be defined:'
     WRITE (*,*) ' '
     DO i = 1,11
        i_fit = 0
        DO j = 1,lm_n_par
           IF (fit_array(j) == i) i_fit = 1
        ENDDO
        IF (i_fit == 0) THEN
           IF (i == 1) THEN
              write (*,*) ' Give segm.No.m / Molar Mass :'
              READ (*,*) parame(1,1)
              parame(1,1) = parame(1,1)*mm(1)
           ELSEIF (i == 2) THEN
              write (*,*) ' Give segm.Diameter sigma    :'
              READ (*,*) parame(1,2)
           ELSEIF (i == 3) THEN
              write (*,*) ' Give segm.Energy Param eps/k:'
              READ (*,*) parame(1,3)
           ELSEIF (i == 4) THEN
              IF ( n_sites == 0) THEN
                 write (*,*) ' Give number of association sites,'
                 write (*,*) ' choose 0 for non-associating compounds'
                 READ (*,*) n_sites
              ENDIF
              IF (n_sites /= 0) THEN
                 write (*,*) ' Give assoc.-energy eps_AB/k :'
                 IF (n_sites == 1) nhb_typ(1) = 1
                 IF (n_sites == 1) nhb_no(1,1)= 1.0
                 IF (n_sites == 1) READ (*,*) eps_hb(1,1,1,1)
                 IF (n_sites == 2) nhb_typ(1) = 2
                 IF (n_sites == 2) nhb_no(1,1)= 1.0
                 IF (n_sites == 2) nhb_no(1,2)= 1.0
                 IF (n_sites == 2) READ (*,*) eps_hb(1,1,1,2)
                 IF (n_sites == 2) eps_hb(1,1,2,1) = eps_hb(1,1,1,2)
                 IF (n_sites == 3) nhb_typ(1) = 2
                 IF (n_sites == 3) nhb_no(1,1)= 1.0
                 IF (n_sites == 3) nhb_no(1,2)= 2.0
                 IF (n_sites == 3) READ (*,*) eps_hb(1,1,1,2)
                 IF (n_sites == 3) READ (*,*) eps_hb(1,1,2,1)
                 IF (n_sites == 4) nhb_typ(1) = 2
                 IF (n_sites == 4) nhb_no(1,1)= 2.0
                 IF (n_sites == 4) nhb_no(1,2)= 2.0
                 IF (n_sites == 4) READ (*,*) eps_hb(1,1,1,2)
                 IF (n_sites == 4) READ (*,*) eps_hb(1,1,2,1)
                 ! parame(1,14) = 0.0                  ! eps_hb(1,1,1,1)
                 ! READ (*,*)     parame(1,15)         ! eps_hb(1,1,1,2)
                 ! parame(1,16) = parame(1,15)         ! eps_hb(1,1,2,1)
                 ! parame(1,17) = 0.0                  ! eps_hb(1,1,2,2)
              ENDIF
           ELSEIF (i == 5 .AND. n_sites /= 0.0) THEN
              write (*,*) ' Give assoc.-volume kappa_AB :'
              READ (*,*) parame(1,13)
              ! parame(1,13) = 0.01        !JG
           ELSEIF (i == 6) THEN
              write (*,*) ' Give dipolar moment         :'
              READ (*,*) parame(1,6)
              ! parame(1,6) = 4.886577        !JG
           ELSEIF (i == 7 .AND. (parame(1,6) /= 0.0 .OR. dipole == 1))THEN
              ! write (*,*) ' Give No. of segm.with dipole:'
              ! READ (*,*) parame(1,8)
           ELSEIF (i == 8) THEN
              write (*,*) ' Give quadrupolar moment     :'
              READ (*,*) parame(1,7)
              ! parame(1,7) = 0.0        !JG
           ELSEIF(i == 9 .AND. (parame(1,7) /= 0.0 .OR. quadrup == 1))THEN
              ! write (*,*) ' Give No. of segm.with quadrp:'
              ! READ (*,*) parame(1,9)
           ELSEIF (i == 11 .AND. (parame(1,6) /= 0.0 .OR. dipole == 1))THEN
              write (*,*) ' Specify polarizability [A**3]'
              READ (*,*) parame(1,11)
              ! parame(1,11) = 0.0        !JG
           ELSEIF(i == 11 .AND. (parame(1,7) /= 0.0 .OR. quadrup == 1))THEN
              write (*,*) ' Specify polarizability [A**3]'
              READ (*,*) parame(1,11)
           ENDIF
        ENDIF
     ENDDO
  ENDIF

  !-----------------------------------------------------------------------------
  ! write association parameters
  !-----------------------------------------------------------------------------

  IF ( nhb_typ(1)  /=  0 ) THEN
     parame(1,12) = REAL(nhb_typ(1))
     ! parame(1,13) = kap_hb(1,1)
     no = 0
     DO j = 1,nhb_typ(1)
        DO k = 1,nhb_typ(1)
           parame(1,(14+no)) = eps_hb(1,1,j,k)
           no = no + 1
        ENDDO
     ENDDO
     DO j = 1,nhb_typ(1)
        parame(1,(14+no)) = nhb_no(1,j)
        no = no + 1
     ENDDO
  ENDIF


  !-----------------------------------------------------------------------------
  ! initialize the array of adjustable parameters
  !-----------------------------------------------------------------------------

  DO  i = 1, lm_n_par
     rel(i)  = dummy(i,2)
     para(i) = dummy(i,3)
  END DO

  !-----------------------------------------------------------------------------
  ! weighting for the different data types (currently all =1.0)
  !-----------------------------------------------------------------------------

  WRITE (*,*) ' '
  WRITE (*,*), 'Input of weighting factors for the different data-sets'
  WRITE (*,*) ' '

  IF ( eqa(1) /= 0 ) THEN
     WRITE (*,*), 'weights for density data ?'
     type_weight(1) = 1.0
     write (*,*) type_weight(1)
     ! READ (*,*) type_weight(1)
  END IF
  IF (eqa(2) /= 0) THEN
     WRITE (*,*), 'weights for vapor pressure data ?'
     type_weight(2) = 1.0
     write (*,*) type_weight(2)
     ! READ (*,*) type_weight(2)
  END IF
  IF (eqa(3) /= 0) THEN
     WRITE (*,*), 'weights for 2nd virial coefficients ?'
     READ (*,*) type_weight(3)
  END IF
  IF (eqa(4) /= 0) THEN
     WRITE (*,*), 'weights for energy data ?'
     READ (*,*) type_weight(4)
  END IF

  ! --- summarizing weighting of the different data types ----------------------
  WRITE (*,*), 'the following weights were specified:'
  WRITE (*,*) ' '
  IF (eqa(1) == 1) WRITE (*,*), ' density data            : ', type_weight(1)
  IF (eqa(2) == 1) WRITE (*,*), ' vapor pressure data     : ', type_weight(2)
  IF (eqa(3) == 1) WRITE (*,*), ' 2nd virial coefficients : ', type_weight(3)
  IF (eqa(4) == 1) WRITE (*,*), ' int. res. energy data   : ', type_weight(4)

  !-----------------------------------------------------------------------------
  ! summarizing the starting values and the scaling of parameters
  !-----------------------------------------------------------------------------

  WRITE (*,*) ' '
  WRITE(23,*) ' '
  WRITE(23,*) ' '
  WRITE(23,'(a)') ' PARAMETER-TYPE           PARAM.-START   SCALING         para'
  WRITE(23,'(a)') ' --------------------------------------------------------------'
  DO i=1,lm_n_par
     WRITE(23,'(a,t26,G13.6,t41,G13.6,t56,G13.6)') aa_txt(i),dummy(i,1),rel(i),para(i) 
  END DO
  WRITE(23,*) ' '
  WRITE(23,'(a)') ' WEIGHTS'
  WRITE(23,'(a)') ' -------'
  IF (eqa(1) == 1) WRITE(23,'(A,T26,F5.2)')' density data           :',type_weight(1)
  IF (eqa(2) == 1) WRITE(23,'(A,T26,F5.2)')' vapor pressure data    :',type_weight(2)
  IF (eqa(3) == 1) WRITE(23,'(A,T26,F5.2)')' 2nd virial coeff.      :',type_weight(3)
  IF (eqa(4) == 1) WRITE(23,'(A,T26,F5.2)')' energy data            :',type_weight(4)
  WRITE(23,*) ' '

END SUBROUTINE paini


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE pure_output ( lm_n_par, para )

  USE BASIC_VARIABLES
  USE pure_fit_parameters
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: lm_n_par
  REAL, INTENT(IN)                       :: para(:)

  !-----------------------------------------------------------------------------
  INTEGER                                :: i,j,ps
  REAL                                   :: devi, devimax, devisum, deviav, rms
  REAL                                   :: p_deviav,p_devimax,r_deviav,r_devimax
  CHARACTER (LEN=4)                      :: char_len
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! writing resulting parameters to screen and to output file
  !-----------------------------------------------------------------------------

  WRITE (*,*) ' '
  WRITE (*,*) '-----------------------------------------------'
  DO i = 1, lm_n_par
     WRITE (*,'(a,2x,f12.6)') aa_txt(i), para(i)*rel(i)
  END DO
  !WRITE (*,'(a,3(2x,f15.5))') ' STD ',(stdval(i), i=1,3)
  WRITE (*,*) '-----------------------------------------------'

  WRITE (23,*) ' '
  DO i = 1, lm_n_par
     WRITE (23,'(T2,A,I2,T26,G16.9)') 'PARAMETER ', i, para(i)*rel(i)
  END DO


  !-----------------------------------------------------------------------------
  ! writing resulting parameters in the format of para_input.f90
  !-----------------------------------------------------------------------------

  WRITE(23, *) ' '
  WRITE(23, *)'-----------------------------------------------'
  WRITE(23, *)'       mm(i)       =',mm(1)
  DO i = 1, 11
     ps = 0
     DO j = 1, lm_n_par
        IF (fit_array(j) == i) ps = j
     END DO

     IF ( i == 1 .AND. ps >= 1 ) THEN
        WRITE(23, '(a,G16.8)') '!        parame(i,1) = mm(i)*',para(ps)*rel(ps)
        WRITE(23, *)'       parame(i,1) =',mm(1)*para(ps)*rel(ps)
     ELSE IF ( i == 1 ) THEN
        WRITE(23, '(a,G16.8)')'!        parame(i,1) = mm(i)*',parame(1,1)
        WRITE(23, *)'       parame(i,1) =',mm(1)*parame(1,1)
     END IF

     IF (i == 2 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,2) =',para(ps)*rel(ps)
     ELSE IF (i == 2) THEN
        WRITE(23, *)'       parame(i,2) =',parame(1,2)
     END IF

     IF (i == 3 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,3) =',para(ps)*rel(ps)
     ELSE IF (i == 3) THEN
        WRITE(23, *)'       parame(i,3) =',parame(1,3)
     END IF

     IF (parame(1,12) == 2.0) THEN
        IF (i == 4 .AND. ps >= 1) THEN
           WRITE(23, *)'       nhb_typ(i)     =   ',NINT(parame(1,12))
           WRITE(23, *)'       nhb_no(i,1)    =   ',NINT(parame(1,18))
           WRITE(23, *)'       nhb_no(i,2)    =   ',NINT(parame(1,19))
           WRITE(23, *)'       eps_hb(i,i,1,2)=   ',para(ps)*rel(ps)
           WRITE(23, *)'       eps_hb(i,i,2,1)=   ',para(ps)*rel(ps)
           WRITE(23, *)'       eps_hb(i,i,1,1)=  0.0'
           WRITE(23, *)'       eps_hb(i,i,2,2)=  0.0'
        ELSE IF (i == 4) THEN
           WRITE(23, *)'       nhb_typ(i)     =   ',NINT(parame(1,12))
           WRITE(23, *)'       nhb_no(i,1)    =   ',NINT(parame(1,18))
           WRITE(23, *)'       nhb_no(i,2)    =   ',NINT(parame(1,19))
           WRITE(23, *)'       eps_hb(i,i,1,2)=   ',parame(1,15)
           WRITE(23, *)'       eps_hb(i,i,2,1)=   ',parame(1,15)
           WRITE(23, *)'       eps_hb(i,i,1,1)=  0.0'
           WRITE(23, *)'       eps_hb(i,i,2,2)=  0.0'
        END IF
        IF (i == 5 .AND. ps >= 1) THEN
           WRITE(23, *)'       kap_hb(i,i)    =   ',para(ps)*rel(ps)
        ELSE IF (i == 5) THEN
           WRITE(23, *)'       kap_hb(i,i)    =   ',parame(1,13)
        END IF
     ELSE IF ( NINT(parame(1,12)) == 1 ) THEN
        IF (i == 4 .AND. ps >= 1) THEN
           WRITE(23, *)'       nhb_typ(i)     =   ',NINT(parame(1,12))
           WRITE(23, *)'       eps_hb(i,i,1,1)=   ',para(ps)*rel(ps)
        ELSE IF (i == 4) THEN
           WRITE(23, *)'       nhb_typ(i)     =   ',NINT(parame(1,12))
           WRITE(23, *)'       eps_hb(i,i,1,1)=   ',parame(1,14)
        END IF
        IF (i == 5 .AND. ps >= 1) THEN
           WRITE(23, *)'       kap_hb(i,i)    =    ',para(ps)*rel(ps)
        ELSE IF (i == 5) THEN
           WRITE(23, *)'       kap_hb(i,i)    =    ',parame(1,13)
        END IF
     END IF

     IF (i == 6 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,6) =      ',para(ps)*rel(ps)
     ELSE IF (i == 6 .AND. parame(1,6) /= 0.0) THEN
        WRITE(23, *)'       parame(i,6) =      ',parame(1,6)
     END IF
     IF (i == 7 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,8) =      ',para(ps)*rel(ps)
     ELSE IF (i == 7 .AND. parame(1,8) /= 0.0) THEN
        WRITE(23, *)'       parame(i,8) =      ',parame(1,8)
     END IF

     IF (i == 8 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,7) =      ',para(ps)*rel(ps)
     ELSE IF (i == 8 .AND. parame(1,7) /= 0.0) THEN
        WRITE(23, *)'       parame(i,7) =      ',parame(1,7)
     END IF
     IF (i == 9 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,9) =      ',para(ps)*rel(ps)
     ELSE IF (i == 9 .AND. parame(1,9) /= 0.0) THEN
        WRITE(23, *)'       parame(i,9) =      ',parame(1,9)
     END IF
     IF (i == 11 .AND. ps >= 1) THEN
        WRITE(23, *)'       parame(i,11) =      ',para(ps)*rel(ps)
     ELSE IF (i == 11 .AND. parame(1,11) /= 0.0) THEN
        WRITE(23, *)'       parame(i,11) =      ',parame(1,11)
     END IF
  END DO
  WRITE(23, *) '-----------------------------------------------'

  !-----------------------------------------------------------------------------
  ! analyizing result
  !-----------------------------------------------------------------------------

  !WRITE (*,*) ' '
  !WRITE (*,'(T2, A, I1, 2X, A, 1PD12.5)') 'Convergence= ', ic, 'SUMERR= ', sse
  !WRITE(23, *)
  !WRITE(23, '(A, I1, 2X, A, 1PD12.5)') 'Convergence= ', ic, 'SUMERR= ', sse
  !WRITE (*,*) ' '
  !WRITE (*,*), 'Normalized correlation of parameters:'
  !WRITE (*,*) ' '
  !WRITE(23,*)
  !WRITE(23,*) 'Normalized correlation of parameters:'
  !FORM='(T2, 2X, '//mul//'(3X, A, 5X))'
  !WRITE (*,FORM) (namp(i), i=1,lm_n_par)
  !WRITE(23, FORM) (namp(i), i=1,lm_n_par)
  !FORM='(T2, I2, '//mul//'(2X, D10.4)/)'
  !DO  i=1,lm_n_par
  !  WRITE (*,FORM) i, (corval(i,k), k=1,lm_n_par)
  !  WRITE(23, FORM) i, (corval(i,k), k=1,lm_n_par)
  !END DO
  WRITE (*,*) ' '
  WRITE(23, *)


  !-----------------------------------------------------------------------------
  ! comparing calculated values to experimental data
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! 1. density data
  !-----------------------------------------------------------------------------

  IF ( eqa(1) == 1 ) THEN
     WRITE (23,*) '------------- fluid density data --------------'
     WRITE(23,'(/A)') 'DATNR   AD (%)     T/K    v_exp/(m3/kg)   v_calc/(m3/kg) '
     devisum = 0.0
     devimax = 0.0
     devi = 0.0
     rms  = 0.0
     ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation -
     DO i = 1, nliq
        devi = ABS(v_cal(i)-vliq(i))/ (SQRT(ABS(vliq(i)))*SQRT(ABS(v_cal(i))))*100.0
        rms = rms + ((v_cal(i)-vliq(i))/vliq(i))**2
        IF (devi > devimax) devimax=devi
        devisum = devisum + devi
        WRITE(23,'(i3,2x,G10.3,3(2x,G12.5))') i, devi, tliq(i), vliq(i), v_cal(i)
     END DO
     deviav = devisum / REAL(nliq)
     rms = SQRT( rms / REAL(nliq) )*100.0
     WRITE(23,*) ' '
     WRITE(23,'(A)') 'Deviation of calculated to exp. volume data'
     WRITE(23,*) ' '
     WRITE(23,'(A, t35, F7.3, A)') 'Max. Deviation:',devimax, ' %'
     WRITE(23,'(A, t35, F7.3, A)') 'Average Deviation: ',deviav, ' %'
     WRITE(23,'(A, t35, F7.3, A)') 'RMS: ',rms, ' %'
     WRITE(23,*) ' '
     r_devimax = devimax
     r_deviav  = deviav
  END IF

  !-----------------------------------------------------------------------------
  ! 2. vapor pressure data
  !-----------------------------------------------------------------------------

  IF ( eqa(2) == 1 ) THEN
     WRITE (23,*) '------------- vapor pressure data -------------'
     WRITE(23,'(/A)') 'DATNR  AD (%)     T/K        PEXP/PA         PCAL/PA'
     devisum = 0.0
     devimax = 0.0
     devi = 0.0
     rms  = 0.0
     DO i = 1, nlv
        ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation
        devi = ABS(plvcal(i)-plv(i))/ (SQRT(ABS(plv(i)))*SQRT(ABS(plvcal(i))))*100.0
        rms = rms + ((plvcal(i)-plv(i))/plv(i))**2
        IF (devi > devimax) devimax = devi
        devisum = devisum + devi
        WRITE (23,'(i3,2x,G10.3,2x,G12.5,2(2x,G14.7))') i, devi, tlv(i), plv(i), plvcal(i)
     END DO
     deviav = devisum / REAL(nlv)
     rms = SQRT( rms / REAL(nlv)) * 100.0
     WRITE(23,*) ' '
     WRITE(23,'(A)') 'Deviation of calculated to exp. vapor pressure data'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F7.3, A/))') 'Max. Deviation:',devimax, ' %',  &
          'Average Deviation: ',deviav, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '
     p_devimax = devimax
     p_deviav  = deviav
  END IF

  !-----------------------------------------------------------------------------
  ! 3.: 2nd virial coefficient data
  !-----------------------------------------------------------------------------

  IF ( eqa(3) == 1 ) THEN
     WRITE (23,*) '---------------2nd virial coeff. --------------'
     WRITE (23,'(/A,A)') 'DATNR  AD (%)   T/K     B_exp     B_cal'
     devisum = 0.0
     devimax = 0.0
     devi = 0.0
     rms  = 0.0
     ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation -
     DO  i = 1, n_vir
        devi = ABS(b_cal(i)-b_vir(i))/ (SQRT(ABS(b_vir(i)))*SQRT(ABS(b_cal(i))))*100.0
        rms = rms + ((b_cal(i)-b_vir(i))/b_vir(i))**2
        IF (devi > devimax) devimax = devi
        devisum = devisum + devi
        WRITE(23,'(i3,2x,4(2x,G12.5))') i, devi, t_vir(i), b_vir(i), b_cal(i)
     END DO
     deviav = devisum / REAL(n_vir)
     rms = SQRT( rms / REAL(n_vir) ) * 100.0
     WRITE(23,*) ' '
     WRITE(23,'(A)') 'Deviation of calculated to exp. 2nd virial coeff. data'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F6.2, A/))') 'Max. Deviation:',devimax, ' %',  &
          'Average Deviation: ',deviav, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '
  END IF


  !-----------------------------------------------------------------------------
  ! 4.: internal energy data
  !-----------------------------------------------------------------------------

  IF ( eqa(4) == 1 ) THEN
     WRITE (23,*) '---------------- energy data ------------------'
     WRITE (23,'(/A,A)') 'DATNR  AD (%)  U_EXP      U_CAL    '
     devisum = 0.0
     devimax= 0.0
     devi = 0.0
     rms  = 0.0
     ! --- write data, calc. RMS-%, AAD-%, and maximum observed deviation -
     DO  i = 1, n_en
        devi = ABS(u_cal(i)-u_en(i))/(SQRT(ABS(u_en(i)))*SQRT(ABS(u_cal(i))))*100.0
        rms = rms + ((u_cal(i)-u_en(i))/u_en(i))**2
        IF (devi > devimax) devimax=devi
        devisum = devisum + devi
        WRITE(23,'(i3,2x,4(2x,G12.5))') i, devi, t_en(i), u_en(i), u_cal(i)
     END DO
     deviav = devisum / REAL(n_en)
     rms = SQRT( rms / REAL(n_en) ) * 100.0
     WRITE(23,*) ' '
     WRITE(23,'(A)') 'Deviation of calculated to exp. internal energy data'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F6.2, A/))') 'Max. Deviation:',devimax, ' %',  &
          'Average Deviation: ',deviav, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '
  END IF

  WRITE(23, *)' '
  WRITE (char_len,'(I3)') lm_n_par+5
  WRITE(23,'(a,'//char_len//'G15.8)') pure_fit_file, mm(1),  &
       (para(i)*rel(i),i=1,lm_n_par),r_deviav,r_devimax,p_deviav,p_devimax
  WRITE(23, *)' '


END SUBROUTINE pure_output
