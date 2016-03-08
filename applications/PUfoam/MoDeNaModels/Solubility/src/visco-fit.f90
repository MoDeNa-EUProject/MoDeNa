!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! MODULE VISC_DATA
!
! This module contains the experimental data as well as the calculated
! properties
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

MODULE VISC_DATA

  USE GC_GROUP, ONLY: N_comp
  implicit none
  save

  INTEGER, PARAMETER                    :: N_visc = 3500
  INTEGER, PARAMETER                    :: N_subst = N_comp

  INTEGER                               :: n_comps
  INTEGER, DIMENSION(N_subst)           :: ndati
  INTEGER, DIMENSION(N_visc)            :: lv_flag

  REAL, DIMENSION(N_visc)               :: t_visc, p_visc, eta_visc, s_res
  REAL, DIMENSION(N_subst)              :: eta_crit, tc_visc, mmv

  CHARACTER (LEN=30)                    :: c_name(N_comp)

  REAL, DIMENSION(N_visc,8)             :: sumarr

END MODULE VISC_DATA


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

MODULE VISC_PARA_TRANSFER
  implicit none
  save

  REAL, ALLOCATABLE                       :: para_transfer(:)

END MODULE VISC_PARA_TRANSFER





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE VISC_RESIDUAL
!
! Subroutine for fitting pure component viscosity parameters. The
! residual between calculated values and experimental data is calculated
! and stored in vector 'fvec'. The vector 'fvec' serves as the objective
! function for the parameter regression.
! VISC_RESIDUAL is called by GC_FITTING via the Levenberg-Marquardt routine.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE VISC_RESIDUAL (lm_m_dat, lm_n_par, para, fvec, iflag)

  USE BASIC_VARIABLES
  USE VISC_DATA
  USE GC_GROUP
  USE VISC_PARA_TRANSFER
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! INTEGER, PARAMETER                      :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, INTENT(IN)                     :: lm_m_dat
  INTEGER, INTENT(IN)                     :: lm_n_par
  REAL, INTENT(IN)                        :: para(:)
  REAL, INTENT(IN OUT)                    :: fvec(:)
  INTEGER, INTENT(IN OUT)                 :: iflag

  !-----------------------------------------------------------------------------
  INTEGER                                 :: i, k, k_ind

  REAL                                    :: toterr, scale, rms
  REAL                                    :: tr, visc, a_p, b_p, c_p, d_p

  REAL                                    :: ao,bo,co,doo,eo,fo,ro,so,wo,po,omega22,eta_ref,t_st
  CHARACTER (LEN=30)                      :: comp_name(N_comp)
  !-----------------------------------------------------------------------------



  ALLOCATE( para_transfer(lm_n_par) )

  para_transfer = para - 1.0

  IF (iflag == 1) WRITE (*,*) ' parameter ',(para_transfer(i), i = 1, lm_n_par)

  !-----------------------------------------------------------------------------
  ! loop over all substances considered in the fitting
  !-----------------------------------------------------------------------------
  k = 0
  DO i = 1, n_comps
     comp_name(1)   = c_name(i)
     parame      = 0.0
     parame(1,5) = 0.12

     CALL GC_PARAMETER ( 1, comp_name )
     ! compna(1)   = c_name(i)
     ! CALL PARA_INPUT            ! retriev pure comp. parameters

     CALL VISC_GC_PARAMETERS ( 1, a_p, b_p, c_p, d_p )

     !--------------------------------------------------------------------------
     ! loop over all data points ndati(i) per substance i
     !--------------------------------------------------------------------------
     DO k_ind = 1, ndati(i)
        k = k + 1

        !-----------------------------------------------------------------------
        ! calculate low density viscosity (Chapman-Enskog)
        !-----------------------------------------------------------------------
        ao = 1.16145
        bo = 0.14874
        co = 0.52487
        doo= 0.77320
        eo = 2.16178
        fo = 2.43787
        ro = -6.435E-4
        so = 18.0323
        wo = -0.76830
        po = 7.27371
        t_st = t_visc(k)/parame(1,2)

        omega22 = ao/t_st**bo + co/EXP(doo*t_st) + eo/EXP(fo*t_st)  &
             + ro*t_st**bo * SIN(so*t_st**wo -po)

        eta_ref = 26.69 / 1.e4 * (mmv(i)/parame(1,1)*t_visc(k))**0.5  &
             / parame(1,2)**2 /omega22 !/parame(1,1)


        !-----------------------------------------------------------------------
        ! viscosity
        !-----------------------------------------------------------------------
        visc = EXP( a_p + b_p*(s_res(k)/parame(1,1))  &
             + c_p*SIN( (s_res(k)/parame(1,1))*d_p ) ) * eta_ref

        !-----------------------------------------------------------------------
        ! calculate the deviation (vector of objective function)
        !-----------------------------------------------------------------------
        scale =  1.E2
        tr = t_visc(k) / tc_visc(i)
        IF (tr > 0.95.AND.tr < 1.02) scale = scale/5.0
        IF (s_res(k)/parame(1,1) < 0.5) scale = scale/10.0
        IF (s_res(k)/parame(1,1) < 0.1) scale = scale/10.0

        fvec(k) = ((visc-eta_visc(k))/eta_visc(k) ) * scale

        !-----------------------------------------------------------------------
        ! store data
        !-----------------------------------------------------------------------
        sumarr(k,1)= t_visc(k)
        sumarr(k,2)= p_visc(k)/1.E6
        sumarr(k,3)= eta_visc(k)
        sumarr(k,4)= visc
        sumarr(k,5)= (visc-eta_visc(k))/eta_visc(k)
        sumarr(k,6)=  s_res(k)/parame(1,1)
        sumarr(k,7)= LOG(eta_visc(k)/eta_ref*parame(1,1))  ! LOG(eta_visc(k)/eta_ref*parame(1,1))
        sumarr(k,8)= LOG(visc/eta_ref*parame(1,1))  ! LOG(visc/eta_ref*parame(1,1))

     END DO

  END DO


  !-----------------------------------------------------------------------------
  ! report results
  !-----------------------------------------------------------------------------
  toterr = 0.0
  rms  = 0.0
  DO i = 1, lm_m_dat
     toterr = toterr + fvec(i)
     rms = rms + sumarr(k,5)**2
  END DO
  rms = SQRT( rms/REAL(lm_m_dat) )
  WRITE (*,*) '--- error ---', toterr, rms
  WRITE(*,*)' '
  WRITE (81,'(25(e18.10))') toterr,(para_transfer(i), i = 1, lm_n_par)

  DEALLOCATE( para_transfer )

END SUBROUTINE VISC_RESIDUAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE VISC_GC_PARAMETERS ( i, a_p, b_p, c_p, d_p )

  USE BASIC_VARIABLES
  USE VISC_PARA_TRANSFER
  USE GC_GROUP
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: i
  REAL, INTENT(OUT)                      :: a_p
  REAL, INTENT(OUT)                      :: b_p
  REAL, INTENT(OUT)                      :: c_p
  REAL, INTENT(OUT)                      :: d_p

  !-----------------------------------------------------------------------------
  INTEGER :: j
  REAL :: sumgroups

  REAL :: tot_vol,coeff_1(20),coeff_2(20),coeff_3(20), coeff_4(20)
  !-----------------------------------------------------------------------------

  DO j = 1, gc_tot_nr_groups(i)

     IF (gc_grouptype(i,j) == '-CH3') THEN  !  -CH3
        !coeff_1(j)= 2.305142471923816E-003 ! 1.780022775738743E-003 (value for individual pure comp. parameters)
        !coeff_2(j)= 1.446683359706391E-002 ! 1.760271779685674E-002 (value for individual pure comp. parameters)
        coeff_1(j)= 1.852032082868327E-003
        coeff_2(j)= 1.362882416274180E-002
        coeff_3(j)= -0.12 !-4.463381242153375E-002
        coeff_4(j)= 1.432 !0.633703893096159
     ELSE IF (gc_grouptype(i,j) == '-CH2-') THEN  ! -CH2-
        !coeff_1(j)= 1.846477481626252E-003 ! 2.007261846879160E-003 (value for individual pure comp. parameters)
        !coeff_2(j)= 3.954298710657512E-002 ! 3.754281129157455E-002 (value for individual pure comp. parameters)
        coeff_1(j)= 2.057135355967743E-003
        coeff_2(j)= 4.008536597449996E-002
        coeff_3(j)= -0.30 !-0.377038448863756
        coeff_4(j)= 1.432  !2.38588281590988
     ELSE IF (gc_grouptype(i,j) == '-CH2- (cyclic)') THEN  ! -CH2-  C6-cyclic
        coeff_1(j)=para_transfer(1)
        coeff_2(j)=para_transfer(2)
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432
     ELSE IF (gc_grouptype(i,j) == '-CH- (aromatic)') THEN  ! -CH-  C6-aromatic
        coeff_1(j)= 1.495965787222886E-003 ! 4.259157289009607E-003 (value for individual pure comp. parameters)
        coeff_2(j)= 3.052614812826703E-002 ! 3.184931188943985E-002 (value for individual pure comp. parameters)
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432
     ELSE IF (gc_grouptype(i,j) == '-C-R  (aromatic)') THEN  ! -C-R  C6-aromatic
        coeff_1(j)= 1.728898150623337E-002 ! 2.080157911175284E-002 (value for individual pure comp. parameters)
        coeff_2(j)= 5.738084449918723E-002 ! 8.521275673789219E-002 (value for individual pure comp. parameters)
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432
     ELSE IF (gc_grouptype(i,j) == '-CH2-COO-CH2-') THEN  ! -CH2-COO-CH2- ester-group (from acid-side)
        coeff_1(j)= 0.0  !  1.723408759561629E-002
        coeff_2(j)= 0.0  !  0.115274566812445 !para_transfer(2) !
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432
     ELSE IF (gc_grouptype(i,j) == '-CH2-COO-CH3') THEN  ! -CH2-COO-CH3 methyl-ester-group (from acid-side)
        coeff_1(j)= 7.649246016441680E-003 ! 7.945628431905050E-003
        coeff_2(j)= 6.365329094090377E-002 ! 6.571605502893441E-002
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432 ! 2.0
     ELSE IF (gc_grouptype(i,j) == '-CH<') THEN  ! -CH<  branch
        coeff_1(j)=para_transfer(1)
        coeff_2(j)=para_transfer(2)
        coeff_3(j)= 0.0
        coeff_4(j)= 1.432
     ELSE IF (gc_grouptype(i,j) == 'individual') THEN  ! individual
        coeff_1(j)=para_transfer(1)
        coeff_2(j)=para_transfer(2)
        coeff_3(j)= -0.18
        coeff_4(j)= 1.432
     ELSE
        WRITE (*,*) 'VISC_GC_PAR: GC-parameter not found !', gc_grouptype(i,j)
        STOP
     END IF

  END DO

  sumgroups = SUM( REAL( gc_groups(i,1:gc_tot_nr_groups(i)) ) )

  tot_vol = 0.0
  c_p = 0.0
  d_p = 0.0
  DO j=1,gc_tot_nr_groups(i)
     tot_vol = tot_vol + REAL(gc_groups(i,j))*msig3(1,j)
     c_p = c_p + REAL(gc_groups(i,j))/sumgroups*coeff_3(j)
     d_p = d_p + REAL(gc_groups(i,j))/sumgroups*coeff_4(j)
  END DO

  a_p = 0.0
  b_p = 0.0
  DO j=1,8
     a_p = a_p + coeff_1(j)*msig3(1,j)*REAL(gc_groups(i,j))                ! A-param.
     b_p = b_p + coeff_2(j)*msig3(1,j)*REAL(gc_groups(i,j))/tot_vol**0.33  ! B-param.
  END DO
  a_p = -( 1.0 + a_p )
  b_p = -( 1.0 + b_p )

END SUBROUTINE VISC_GC_PARAMETERS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE viscofit_drive

  USE PARAMETERS, ONLY: RGAS
  USE BASIC_VARIABLES
  USE VISC_DATA
  USE GC_GROUP
  USE Levenberg_Marquardt
  use utilities, only: file_open, critical, SI_DENS, enthalpy_etc
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE VISC_RESIDUAL(lm_m_dat, lm_n_par, para, fvec, iflag)
       IMPLICIT NONE
       INTEGER, INTENT(IN)                :: lm_m_dat, lm_n_par
       REAL, INTENT(IN)                   :: para(:)
       REAL, INTENT(IN OUT)               :: fvec(:)
       INTEGER, INTENT(IN OUT)            :: iflag
     END SUBROUTINE VISC_RESIDUAL
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

  INTEGER                                :: i, j, k, k_ind
  INTEGER                                :: n0, converg
  REAL                                   :: av_dev, max_dev, rms

  REAL                                   :: lnphi(np,nc), density(np), w(np,nc)
  REAL                                   :: tc , pc, rhoc, end_x

  CHARACTER (LEN=20)                     :: dum
  CHARACTER (LEN=50)                     :: viscfile
  CHARACTER (LEN=30)                      :: comp_name(N_comp)
  !-----------------------------------------------------------------------------



  ncomp     = 1

  !-----------------------------------------------------------------------------
  ! numerical constants for Levenberg-Marquardt algorithm
  !-----------------------------------------------------------------------------
  scaling(1) = 1.E2
  IF ( num == 0 ) epsfcn = 1.0E-6**2  ! sqrt of relat. step size (finite differences)
  IF ( num == 1 ) epsfcn = 1.0E-5**2  ! sqrt of relat. step size (finite differences)
  lm_factor = 1.0                     ! maximum initial step size (parameter iteration)


  OPEN (81,FILE = './output_file/viscosity/visc_inter.xlo')

  viscfile = './input_file/viscosity/visc_pars.inp'
  CALL file_open ( viscfile, 82 )

  READ (82,*) lm_n_par

  ALLOCATE( para(lm_n_par) )
  ALLOCATE( ipvt(lm_n_par) )

  DO i = 1, lm_n_par
     READ (82,*) para(i)
  END DO

  viscfile = './input_file/viscosity/viscofit.inp'
  CALL file_open(viscfile,85)
  READ(85,*) dum, n_comps
  READ(85,*) dum
  READ(85,*) dum
  DO i = 1, n_comps
     READ(85,*) c_name(i)
  END DO
  CLOSE (85)
  CALL GET_GC_GROUPS ( n_comps, c_name ) !, gc_tot_nr_groups, gc_groups, gc_grouptype)


  lm_m_dat = 0
  DO i=1,n_comps
     viscfile='./input_file/viscosity/'//c_name(i)
     viscfile=trim(viscfile)//'.dat'
     CALL file_open(viscfile,85)
     READ(85,*) dum
     READ(85,*) dum
     READ(85,*) dum
     READ(85,*) dum
     READ(85,*) dum
     READ (85,*) ndati(i), eta_crit(i)
     DO k = lm_m_dat+1, lm_m_dat+ndati(i)
        READ(85,*) t_visc(k),p_visc(k),eta_visc(k),lv_flag(k)
        p_visc(k) = p_visc(k)*1.E6
     END DO
     lm_m_dat = lm_m_dat+ndati(i)
     CLOSE (85)
  END DO

  WRITE (*,*) 'Number of data points:',lm_m_dat
  ALLOCATE( fvec(lm_m_dat) )


  eos = 1
  pol = 1
  nphas = 1
  xi(1,1) = 1.0
  k = 0
  DO i = 1, n_comps
     comp_name(1)= c_name(i)
     parame   = 0.0
     parame(1,5) = 0.12
     ncomp = 1
     CALL gc_parameter ( 1, comp_name )
     !        CALL PARA_INPUT            ! retriev pure comp. parameters


     tc = 1.3*parame(1,3)*parame(1,1)
     pc = 5.E6
     mmv(i)=mm(1)
     CALL critical (tc,pc,rhoc)
     tc_visc(i) = tc
     WRITE(*,'(a,3(f16.5))') 'tc,pc,rhoc',tc, pc,rhoc
     DO k_ind = 1, ndati(i)
        k = k + 1
        p = p_visc(k)
        t = t_visc(k)
        IF (lv_flag(k) == 2) THEN
           nphas = 2
           val_init(1) = 0.5
           val_init(2) = 1.E-3
           val_init(3) = t
           val_init(4) = p
           n_unkw = ncomp          ! number of quantities to be iterated
           it(1) = 'p '            ! iteration of pressure
           val_conv(2) = 0.0
           CALL objective_ctrl (converg)
           IF (converg /= 1) THEN
              WRITE (*,*) tc,t_visc(k), c_name(i)
              val_init(3) = t_visc(k)*0.7
              val_init(4) = 100.0
              IF ((tc/t_visc(k)) > 2.0) val_init(3) = t_visc(k)*1.3
              end_x = t_visc(k)
              running = 't'         ! Temp. is running variable in PHASE_EQUILIB
              CALL phase_equilib(end_x,5.0,converg)
           END IF
           WRITE (*,*) t,val_conv(4)
           IF(converg == 1) p_visc(k) = val_conv(4)
           IF(converg == 1) p = val_conv(4)
           IF(converg /= 1)WRITE (*,*) 'pt. T=',t_visc(k),'not conv.'
           IF(converg /= 1)STOP
           nphas = 1
        END IF
        IF (eta_visc(k) >= eta_crit(i)) densta(1) =0.45
        IF (eta_visc(k) < eta_crit(i)) densta(1) =1.E-5
        CALL fugacity (lnphi)
        CALL si_dens (density,w)
        CALL enthalpy_etc
        s_res(k) = entrop(1)/RGAS -  &
             LOG( p_visc(k)/density(1)/RGAS/t_visc(k) *mmv(i)*1.E-3 )    ! s_res = S_res/(N k) for given rho,T

     END DO
  END DO



  !-----------------------------------------------------------------------------
  ! initialize parameter vectors
  !-----------------------------------------------------------------------------
  DO i=1,lm_n_par
     para(i) = para(i) + 1.0
  END DO


  ! ----------------------------------------------------------------------
  ! adjust parameters (Levenberg-Marquardt scheme)
  ! ----------------------------------------------------------------------

  CALL lmdif1(VISC_RESIDUAL,lm_m_dat,lm_n_par,para,fvec,tol,epsfcn,lm_factor,info,ipvt)


  viscfile='./output_file/viscosity/visc_pars.inp'
  CALL file_open(viscfile,83)
  DO i=1,lm_n_par
     WRITE (83,*) para(i) - 1.0
  END DO
  WRITE (83,*) ' '
  WRITE (83,*) 'standard deviations:'
  DO i=1,lm_n_par
     WRITE (83,*) i !, stdval(i)
  END DO


  k = 0
  DO i = 1, n_comps
     WRITE (81,*) ' '
     DO k_ind = 1, ndati(i)
        k = k + 1
        WRITE (81,'(8(e16.8))') (sumarr(k,j),j = 1, 5)
     END DO
  END DO

  DO k = 1, n_comps

     OPEN (23,FILE='./output_file/viscosity/'  &
          //trim( adjustl( c_name(k) ))//'.dat')
     !      WRITE(23,*) c_name(k)
     !      WRITE(23,*) '----------------------------------------------------'
     WRITE(23,'(a94)')' no T/K p/MPa sres/k LN(eta_exp/eta_ref) '  &
          //'LN(eta_cal/eta_ref) eta_exp/mPa*s eta_cal/mPa*s err%'
     av_dev  = 0.0
     max_dev = 0.0
     rms     = 0.0
     n0 = sum( ndati(1:k-1) )
     !          syntax:
     !          sumarr(k,1)= T_visc(k)
     !          sumarr(k,2)= p_visc(k)/1.E6
     !          sumarr(k,3)= eta_visc(k)
     !          sumarr(k,4)= visc
     !          sumarr(k,5)= (visc-eta_visc(k))/eta_visc(k)
     !          sumarr(k,6)=  s_res(k)/parame(1,1)
     !          sumarr(k,7)= LOG(eta_visc(k)/eta_ref)
     !          sumarr(k,8)= LOG(visc/eta_ref)
     DO i = n0+1 , n0+ndati(k)
        WRITE(23,'( I4, 2X, 8(2x,G14.6) )') i-n0, t_visc(i),  &
             p_visc(i)/1.E6, sumarr(i,6), sumarr(i,7), sumarr(i,8),  &
             eta_visc(i), sumarr(i,4), sumarr(i,5)*100.0
        max_dev = MAX( max_dev, ABS( sumarr(i,5) )*100.0 )
        rms     = rms + sumarr(i,5)**2
        av_dev  = av_dev +  ABS( sumarr(i,5) )*100.0
     END DO
     av_dev = av_dev / REAL(ndati(k))
     rms    = SQRT( rms / REAL(ndati(k)) )*100.0
     WRITE(23,*) ' '
     WRITE(23,*) 'Deviation of calculated to exp. viscosity'
     WRITE(23,*) ' '
     WRITE(23, '(3(A, t35, F7.3, A/))') 'Max. Deviation:',max_dev, ' %',  &
          'Average Deviation: ',av_dev, ' %', 'RMS: ',rms, ' %'
     WRITE(23,*) ' '

     CLOSE(23)
  END DO

END SUBROUTINE viscofit_drive
