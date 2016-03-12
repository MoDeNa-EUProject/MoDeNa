!> \file main_prog.f90
!! \brief Main program.
!!
!! 

!> \mainpage The PC-SAFT Documentation
!!
!! \section sec_intro Introduction
!! This is a program for the calculation of physical properties\n
!! and phase equilibria using the PC-SAFT equation of state.\n
!! 
!! The Code was written by \n
!! <B> Joachim Gross \n University of Stuttgart \n 
!! Institute for Thermodynamics and \n Thermal Process Engineering </B> \n
!! http://www.itt.uni-stuttgart.de
!!
!! \section sec_help Where can I get help?
!! If you stumble upon some problems while using or even extending this documentation, dont despair! \n \n
!! There are at least three simple ways to get help:
!! 1. Look at http://www.stack.nl/~dimitri/doxygen/index.html , which is a very good doxygen documentation.
!! 2. Open your favorite browser and ask <B> google </B>. Maybe like <B> doxygen bold text </B> ...
!! 3. Have a chat with Madlen. But probably she will do one of the above things anyway. ;-)

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!> \brief Starting Point of PC-SAFT program.
!!
!! Here, some basic options are set, the model constants are read \n
!! and the bubble point routine is called.
!! 
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

PROGRAM PC_SAFT

  !-----------------------------------------------------------------------------
  USE BASIC_VARIABLES
  USE EOS_CONSTANTS
  USE EOS_VARIABLES, ONLY: dd_term, qq_term, dq_term
  USE KIJ_FITTING
  USE BubblePoint
  IMPLICIT NONE

  INTEGER                                :: option
  !-----------------------------------------------------------------------------

  num = 0                     ! (num=0: using analytical derivatives of EOS)
                              ! (num=1: using numerical derivatives of EOS)
                              ! (num=2: White's renormalization group theory)

  ensemble_flag = 'tp'        ! this specifies, whether the eos-subroutines 
                              ! are executed for constant 'tp' or for constant 'tv'

  RGT_variant = 'phase_cell'  ! only for calc. including critical point renormalization:
                              ! choose variant as 'phase_cell' or 'isomorphic'

  dd_term = 'GV'              ! dipole-dipole term of Gross & Vrabec (2006)
  qq_term = 'JG'              ! quadrupole-quadrupole term of Gross (2005)
  dq_term = 'VG'              ! dipole-quadrupole term of Vrabec & Gross (2008)
  IF ( num /= 0 ) CALL SET_DEFAULT_EOS_NUMERICAL

  !-----------------------------------------------------------------------------
  ! EOS constants
  !-----------------------------------------------------------------------------
  CALL EOS_CONST ( ap, bp )


     !Perform bubble point calculation, i.e. determine 
     !the composition of a vapor phase which is in equilibrium
     !with the liquid phase that was specified in the input file 
     CALL BubblePointCalculation

END PROGRAM PC_SAFT



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE SET_DEFAULT_EOS_NUMERICAL

  USE EOS_NUMERICAL_DERIVATIVES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------

  ideal_gas   = 'no'                    ! ( yes, no )
  hard_sphere = 'CSBM'                  ! ( CSBM, no )
  chain_term  = 'TPT1'                  ! ( TPT1, HuLiu, no )
  disp_term   = 'PC-SAFT'               ! ( PC-SAFT, CK, PT1, PT2, PT_MF, PT_MIX, no )
  hb_term     = 'TPT1_Chap'             ! ( TPT1_Chap, no )
  LC_term     = 'no'                    ! ( MSaupe, OVL, no )
  branch_term = 'no'                    ! ( TPT2, no )
  II_term     = 'no'
  ID_term     = 'no'

  subtract1   = 'no'                    ! (1PT, 2PT, no)
  subtract2   = 'no'                    ! (ITTpolar, no)

END SUBROUTINE SET_DEFAULT_EOS_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE FLASH
!
! This subroutine performes a 2-phase flash calculation for multi-
! component mixtures. The feed-composition needs to be provided in
! the file INPUT.inp. The subroutine START_VAR actually perfoms the
! flash calculation.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE FLASH

  USE BASIC_VARIABLES
  USE STARTING_VALUES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                               :: converg
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2

  CALL START_VAR (converg) ! gets starting values, sets "val_init"

  IF (converg == 1) CALL OUTPUT
  IF (converg /= 1) write (*,*) ' NO PHASE SPLIT DETECTED'

END SUBROUTINE FLASH



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PURE_COMP
!
! This subroutine performes VLE calculations (option 1) for pure
! substances. It also allows the calculation of isotherms and
! isobars of pure components (option 2 and 3, respectively)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PURE_COMP

  USE BASIC_VARIABLES
  use utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, converg, calcopt, state
  REAL                                   :: steps, end_x, end_p, end_t, start_p, start_t
  REAL                                   :: lnphi(np,nc), density(np), w(np,nc), t_sav
  REAL                                   :: tc, pc, rhoc
  REAL                                   :: pges, zges, b_cal
  !-----------------------------------------------------------------------------
  steps = 20.0

  CALL READ_INPUT
  nphas = 2

  IF ( ncomp /= 1 ) THEN
     write (*,*) 'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE:'
     write (*,*) '    ./input_file/INPUT.INP'
     stop
  ENDIF

  WRITE (*,*) ' CALCULATE VLE                      (1)'
  WRITE (*,*) ' CALCULATE ISOTHERMS                (2)'
  WRITE (*,*) ' CALCULATE ISOBARS                  (3)'
  WRITE (*,*) ' CALCULATE VIRIAL COEFFICIENTS      (4)'
  WRITE (*,*) ' CRITICAL POINT                     (5)'
  READ  (*,*) calcopt

  IF (calcopt == 1) THEN

     n_unkw = ncomp           ! number of quantities to be iterated 
     it(1)   = 'p'            ! iteration of p
     running = 't'            ! Temp. is running variable in PHASE_EQUILIB

     t_sav=t

     val_init(1) = 0.5        ! starting density for a liquid phase
     val_init(2) = 1.E-6      ! starting density for a vapor phase
     val_init(3) = t_sav
     val_init(4) = 1.E-1      ! default starting value for P: 0.1 Pa
     IF (eos >= 4) val_init(4) = 1.E-3  ! default starting value for LJ-model
     val_init(5) = 0.0        ! logar. mole fraction: lnx=0, x=1, phase 1
     val_init(6) = 0.0        ! logar. mole fraction: lnx=0, x=1, phase 2

     write (*,*) '   CHOOSE END TEMPERATURE [units as in input-file]:'
     READ (*,*) end_x
     end_x = end_x + u_in_T

     outp = 1                 ! output to terminal
     CALL PHASE_EQUILIB(end_x,steps,converg)
     tc = t
     pc = p
     IF (eos < 4) CALL CRITICAL ( tc, pc, rhoc )
     ! IF (eos >= 4 .AND. eos /= 7) CALL LJ_CRITICAL (tc,pc,rhoc)
     ! IF (eos >= 4 .AND. eos /= 7) CALL eLJ_CRITICAL (tc,pc,rhoc)
     ! IF (eos >= 4 .AND. eos /= 7) write(*,*) ' eLJ routine active !!!'
     ! IF (eos >= 4 .AND. eos /= 7) write(*,*) ' eLJ routine active !!!'
     ! IF (eos == 7) CALL SW_CRITICAL (tc,pc,rhoc)
     WRITE(*,'(a,3(f16.5))') 'tc, pc, rhoc',tc-u_out_T, pc/u_out_P, rhoc
     WRITE (40,'(5(2x,f18.8))') tc-u_out_T, pc/u_out_P, rhoc, rhoc
     WRITE (40,*) ' '

  ELSEIF ( calcopt == 2 .OR. calcopt == 3 ) THEN

     nphas  = 1
     xi(1,1) = 1.0

     end_p = p
     end_t = t
     start_p = p
     start_t = t
     IF ( calcopt == 2 ) THEN
        WRITE(*,*) ' SPECIFY END-PRESSURE  [units as in input-file]'
        READ (*,*) end_p
        end_p = end_p*u_in_P
     ELSE
        WRITE(*,*) ' SPECIFY END-TEMPERATURE  [units as input-file]'
        READ (*,*) end_t
        end_t = end_t+u_in_T
     ENDIF
     WRITE (*,*) ' SPECIFY LIQUID (1) or VAPOR (2) state'
     READ (*,*) state
     densta(1) = 0.5         ! starting density for a liquid phase
     IF (state == 2) densta(1) = 1.E-5   ! density for a vapor phase

     DO i = 0, INT(steps)
        p = start_p + (end_p-start_p)*REAL(i)/steps
        t = start_t + (end_t-start_t)*REAL(i)/steps
        CALL FUGACITY (lnphi)
        IF (num /= 2) CALL ENTHALPY_ETC
        CALL SI_DENS (density,w)
        WRITE (*,'(8(2x,g16.9))') t-u_out_T,p/u_out_P,dense(1),density(1),cpres(1),enthal(1),  &
             gibbs(1), speed_of_sound(1)
        WRITE(40,'(7(2x,g16.9))') t-u_out_T,p/u_out_P,density(1),cpres(1),enthal(1),gibbs(1),  &
             speed_of_sound(1)
     ENDDO

  ELSEIF (calcopt == 4) THEN

     nphas  = 1
     xi(1,1) = 1.0
     dense(1)= 1.E-7

     start_t = t
     WRITE(*,*) ' SPECIFY END-TEMPERATURE  [units as input-file]'
     READ (*,*) end_t
     end_t = end_t + u_in_T

     WRITE (*,*) '   T    ,B / cm**3/mol'
     WRITE(40,*) '   T    ,B / cm**3/mol'

     DO i = 0, INT(steps)
        t = start_t + (end_t-start_t)*REAL(i)/steps
        CALL P_CALC (pges,zges)
        CALL SI_DENS (density,w)

        b_cal = (ZGES-1.0)/density(1)*1000.0*mm(1)
        WRITE (*,'(2(2x,f14.8))') t-u_out_T,b_cal
        WRITE(40,'(2(2x,f14.8))') t-u_out_T,b_cal
     ENDDO

  ELSE
     tc = t
     pc = p
     IF ( eos < 4) CALL CRITICAL (tc,pc,rhoc)
     WRITE (*,'(a,3(f16.5))') 'tc,pc,rhoc',tc-u_out_T, pc/u_out_P,rhoc
     WRITE (40,'(5(2x,f18.8))') tc-u_out_T, pc/u_out_P, rhoc, rhoc
     WRITE (40,*) ' '
  ENDIF

END SUBROUTINE PURE_COMP

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE TEMP

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: eta, rho, pges
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: lnphi(np,nc)
  !-----------------------------------------------------------------------------
  num = 1
  CALL SET_DEFAULT_EOS_NUMERICAL

  CALL READ_INPUT
  nphas  = 1
  ensemble_flag = 'tv'
  xi(1,1) = 1.0

  DO i = 1, 100
     !densta(1) = 0.45
     densta(1) = REAL(101-i) / 100.0 * 0.45
     dense(1) = densta(1)
     !p = 101.E5 - REAL(i) / 100.0 * 100.E5

     CALL FUGACITY (lnphi)
     write (*,'(3G18.8)') lnphi(1,1)+log(rho), pges/1.E5, eta
     !CALL ENTHALPY_ETC
  END DO

END SUBROUTINE TEMP



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE RGT_FIT_PARAMETERS
!
! ......
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE RGT_FIT_PARAMETERS

  USE BASIC_VARIABLES
  USE FITTING_RGT_PARAMETERS
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE CRIT_DEV ( fmin, optpara, n )
       INTEGER, INTENT(IN)                     :: n
       REAL, INTENT(IN)                        :: optpara(:)
       REAL, INTENT(IN OUT)                    :: fmin
     END SUBROUTINE CRIT_DEV
  END INTERFACE

  !-----------------------------------------------------------------------------
  INTEGER                                 :: nadjust, PRIN
  REAL                                    :: fmin, t0, h0, MACHEP
  REAL, ALLOCATABLE                       :: optpara(:)

  ! INTEGER                                :: IERROR !,LW,PRIN
  ! EXTERNAL CRIT_DEV,CRIT_DEV2 !,func    ! deactivated func when making a transition to f90
  ! REAL                                   :: CRIT_DEV,CRIT_DEV2 !,WVEC(14*3),GVEC(3)
  REAL                                   :: pp(3) !,func ! ,pam(4,3),yam(4),xi(3,3)
  !-----------------------------------------------------------------------------

  num = 2
  CALL READ_INPUT

  IF (ncomp /= 1) THEN
     write(*,*)'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE'
     stop
  ENDIF

  nadjust = 3
  ALLOCATE( optpara(nadjust) )

  critfit = 1

  tcf = 568.9   ! octane
  pcf = 24.9d5
  rcf = 232.0
  tcf = 373.4   ! h2s
  pcf = 89.70d5
  rcf = 388.5
  tcf = 304.21   ! co2
  pcf = 73.825d5
  rcf = 466.01
  tcf = 416.25   ! chloromethane
  pcf = 66.9d5
  rcf = 363.0
  tcf = 386.5   ! 11-difluoro-ethane (R152a)
  pcf = 45.198d5
  rcf = 367.9
  tcf = 190.555   ! methane
  pcf = 45.95d5
  rcf = 162.2
  ! tcf = 512.64   ! methanol
  ! pcf = 81.35d5
  ! rcf = 272.0
  tcf = 425.1   ! butane
  pcf = 38.E5
  rcf = 3.92*58.123
  tcf = 617.8   ! decane
  pcf = 21.1d5
  rcf = 232.
  ! tcf = 369.8   ! propane
  ! pcf = 42.5d5
  ! rcf = 225.
  tcf = 126.19   ! n2
  pcf = 33.978d5
  rcf = 312.0
  tcf = 508.15   ! acetone
  pcf = 47.62d5
  rcf = 269.0

  ! optpara(1) = (2.53  -1.0 )
  ! optpara(2) = (40.06 -20.0)  /20.0
  ! optpara(3) = (0.165 -0.12 )  *10.0
!!$ optpara(1) = (1.6  -1.0 )   *5.0
!!$ optpara(2) = (18.0 -10.0)  /20.0
!!$ optpara(3) = (0.55 -0.4 )  *10.0

  optpara(1) = (2.1  -0.0 )   *1.0
  optpara(2) = (35.0 -0.0 )   /20.0
  optpara(3) = (0.23 -0.0 )   *5.0


  ! t0 = 1.E-4
  ! optpara(1) = 2.33          +0.2
  ! optpara(2) = 30.06  /20.0 +0.2
  ! optpara(3) = 0.155  *5.0
  ! pam(1,1) = optpara(1)
  ! pam(1,2) = optpara(2)
  ! pam(1,3) = optpara(3)
  ! yam(1) = CRIT_DEV (optpara,nadjust)
  ! optpara(1) = 2.33          -0.2
  ! optpara(2) = 30.06  /20.0 -0.2
  ! optpara(3) = 0.155  *5.0
  ! pam(2,1) = optpara(1)
  ! pam(2,2) = optpara(2)
  ! pam(2,3) = optpara(3)
  ! yam(2) = CRIT_DEV (optpara,nadjust)
  ! optpara(1) = 2.33          +0.0
  ! optpara(2) = 30.06  /20.0 -0.0
  ! optpara(3) = 0.155  *5.0  +0.1
  ! pam(3,1) = optpara(1)
  ! pam(3,2) = optpara(2)
  ! pam(3,3) = optpara(3)
  ! yam(3) = CRIT_DEV (optpara,nadjust)
  ! optpara(1) = 2.33          -0.0
  ! optpara(2) = 30.06  /20.0 -0.0
  ! optpara(3) = 0.155  *5.0  -0.1
  ! pam(4,1) = optpara(1)
  ! pam(4,2) = optpara(2)
  ! pam(4,3) = optpara(3)
  ! yam(4) = CRIT_DEV (optpara,nadjust)
  ! CALL amoeba(pam,yam,4,3,nadjust,t0,CRIT_DEV,IERROR)


  ! pp(1) = optpara(1)
  ! pp(2) = optpara(2)
  ! pp(3) = optpara(3)
  ! xi(1,1)=1.0
  ! xi(2,1)=1.0
  ! xi(1,2)=-1.0
  ! xi(2,2)=1.0
  ! xi(3,3)=1.0
  ! t0 = 1.E-4
  ! CALL powell(pp,xi,nadjust,nadjust,t0,IERROR,fmin)
  ! write (*,*) 'done'
  ! fmin=func (pp)

  pp(1) = optpara(1)
  pp(2) = optpara(2)
  pp(3) = optpara(3)
  t0 = 1.E-4
  ! CALL frprmn(pp,nadjust,t0,IERROR,fmin)    ! deactivated this line when making a transition to f90
  write (*,*) 'done'
  ! fmin=func (pp)    ! deactivated this line when making a transition to f90


  t0 = 1.E-4
  h0 = 0.4
  PRIN = 0
  MACHEP = 1.E-12
  !
  call PRAXIS( t0, MACHEP, h0, nadjust, PRIN, optpara, CRIT_DEV, fmin )
  !call CRIT_DEV ( fmin, optpara, nadjust )
  !
  !! fmin = 1.E-3
  !! LW = 14*nadjust
  !! CALL TN (IERROR,nadjust,optpara,fmin,GVEC,WVEC,LW,CRIT_DEV2)


  WRITE (*,'(a,1(f16.8))') 'deviation %',fmin*100.0
  WRITE (*,*) optpara(1),optpara(2)*20.0,optpara(3)/5.0
  WRITE (*,*) ' '
  DEALLOCATE( optpara )

END SUBROUTINE RGT_FIT_PARAMETERS


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE CRIT_DEV ( fmin, optpara, nadjust )

  USE BASIC_VARIABLES
  USE FITTING_RGT_PARAMETERS
  use utilities, only: CRITICAL
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)                    :: nadjust
  REAL, INTENT(IN)                       :: optpara(nadjust)
  REAL, INTENT(IN OUT)                   :: fmin
  !-----------------------------------------------------------------------------
  REAL                                   :: tc,pc,rhoc
  !-----------------------------------------------------------------------------

  LLfit   = optpara(1)
  phifit  = optpara(2)  *20.0
  chapfit = optpara(3)  /5.0

  tc   = tcf + 40.0
  pc   = pcf
  rhoc = rcf

  WRITE (*,'(3(f18.12))') optpara(1),optpara(2)*20.0,optpara(3)/5.0
  CALL CRITICAL ( tc, pc, rhoc )

  fmin =   100.0 * ABS(1.0-tc/tcf  )**2  &
       +         ABS(1.0-pc/pcf  )**2  &
       + 1.E-1 * ABS(1.0-rhoc/rcf)**2

  ! fmin = 10.0*ABS(1.0-optpara(1)  )**2 + ABS(2.0-optpara(2)  )**1.5  &
  !          +1.E-1*ABS(3.0-optpara(3))**2.5

  WRITE (*,'(a,3(f16.5))') 'tc,pc,rhoc',tc-u_out_T, pc/u_out_P,rhoc
  WRITE (*,'(a,1(f16.8))') 'deviation %',fmin*100.0

  WRITE (40,'(5(2x,f18.6))') tc-u_out_T, pc/u_out_P, rhoc, rhoc
  WRITE (40,*) fmin,optpara(1),optpara(2)*20.0,optpara(3)/5.0

END SUBROUTINE CRIT_DEV


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE LC_COMP

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: alpha_nematic
  USE EOS_NUMERICAL_DERIVATIVES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: converg, ncompsave, phaseindex
  REAL                                   :: steps, end_x, startcon
  REAL                                   :: p_sav, t_sav
  !-----------------------------------------------------------------------------
  steps = 1.0
  num = 1

  !-----------------------------------------------------------------------------
  ! define the Helmholtz energy contributions
  !-----------------------------------------------------------------------------
  CALL SET_DEFAULT_EOS_NUMERICAL
  chain_term  = 'HuLiu'                 ! ( TPT1, HuLiu, no )
  disp_term   = 'NO'                    ! ( PC-SAFT, CK, no )
  hb_term     = 'TPT1_Chap'             ! ( TPT1_Chap, no )
  LC_term     = 'OVL'                   ! ( MSaupe, OVL, no )

  !-----------------------------------------------------------------------------
  ! read input and start by re-setting the number of sustances =1 (only the pure LC)
  !-----------------------------------------------------------------------------
  CALL READ_INPUT
  nphas = 2
  ncompsave = ncomp
  ncomp = 1

  !-----------------------------------------------------------------------------
  ! for the Mayer-Saupe theory, this used to be it(1)='t' and running='p'
  !-----------------------------------------------------------------------------
  n_unkw = ncomp      ! number of quantities to be iterated 
  it(1) = 'p'         ! iteration of pressure
  running = 't'       ! Temp. is running variable in PHASE_EQUILIB

  t_sav = t
  p_sav = p

  !-----------------------------------------------------------------------------
  ! default starting values
  !-----------------------------------------------------------------------------
  val_init(1) = 0.35
  val_init(2) = 0.35
  val_init(3) = t_sav
  val_init(4) = p_sav
  val_init(5) = 0.0  ! logar. mole fraction: lnx=0, x=1, phase 1
  val_init(6) = 0.0  ! logar. mole fraction: lnx=0, x=1, phase 2

  !-----------------------------------------------------------------------------
  ! starting value for the OVL-order parameter
  !-----------------------------------------------------------------------------
  alpha_nematic = 10.0   ! starting value for alpha of OLV-Theory

  !-----------------------------------------------------------------------------
  ! for diagrams: define end-point
  !-----------------------------------------------------------------------------
  end_x = t
  IF (ncompsave == 1) THEN
     !write(*,*)'   CHOOSE END PRESSURE [units as in input-file]:'
     !READ (*,*) end_x
     !end_x = end_x*u_in_P
  ENDIF

  !-----------------------------------------------------------------------------
  ! perform phase equilibrium calculation
  !-----------------------------------------------------------------------------
  outp = 1                         ! output to terminal
  CALL PHASE_EQUILIB(end_x,steps,converg)

  write (*,*) 'pressure=',p/1.E6,'MPa'
  write (*,*) 'eta     =',dense(1),dense(2)
  write (*,*) 'alpha    =',alpha_nematic


  !-----------------------------------------------------------------------------
  ! for mixtures: now extend the number of species
  !-----------------------------------------------------------------------------
  IF (ncompsave > 1) THEN

     ncomp = ncompsave
     n_unkw = ncomp       ! number of quantities to be iterated 
     startcon = LOG(0.000001)
     write(*,*) ' specify desired end-composition of:'
     write(*,*) ' component 2 (molefraction)'
     write(*,*) ' component 2 is:',compna(2)
     read (*,*) end_x
     ! end_x = LOG(0.01)
     ! end_x = 0.14861586
     ! end_x = 0.05489185

     val_init(5) = LOG(1.0 - EXP(startcon - 0.0))
     val_init(6) = startcon - 0.0
     val_init(7) = LOG(1.0 - EXP(startcon))
     val_init(8) = startcon

     write(*,*) ' Is this composition specified for the LC-phase (1)'
     write(*,*) ' or for the isotropic phase (2)'
     read (*,*) phaseindex

     IF (phaseindex == 1) THEN
        it(2)='x22'        ! iteration of mol fraction of comp.1 phase 2
        running = 'x12'      ! lnx(1,1) is running var. in PHASE_EQUILIB
     ELSE
        it(2)='x12'        ! iteration of mol fraction of comp.1 phase 2
        running = 'x22'      ! lnx(1,1) is running var. in PHASE_EQUILIB
     ENDIF
     it(1)='t'    ! iteration of temperature
     sum_rel(1)='x11'   ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2)='x21'   ! summation relation: x22 = 1 - sum(x2j)
     outp = 1           ! output to terminal
     steps=10.0


     CALL PHASE_EQUILIB(end_x,steps,converg)

     write(*,*)' '
     ! write(*,*) ' SPECIFY TEMPERATURE (1)'
     ! write(*,*) ' OR PRESSURE         (2)'
     ! READ (*,*) iso_x
     ! IF (iso_x == 1) THEN
     write(*,*) '  CHOOSE END TEMPERATURE [units as in input-file]:'
     READ (*,*) end_x
     end_x = end_x + u_in_T
     it(1)='p'    ! iteration of pressure
     running='t'          ! Temperature is running var. in PHASE_EQUILIB

     CALL PHASE_EQUILIB(end_x,steps,converg)

     write (77,*) val_conv
  ENDIF

END SUBROUTINE LC_COMP



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE LC_3PHASE
!
! This subroutine performes VLLE calculations for binary mixtures
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE LC_3PHASE

  USE BASIC_VARIABLES
  use utilities, only: SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, converg
  REAL                                   :: density(np), w(np,nc)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT

  write(*,*) 'we assume the phase 1 to have a target composition !!'

  nphas  = 3
  n_unkw = 4                   ! number of quantities to be iterated 
  it(1) = 'p  '                ! iteration of pressure
  it(2) = 'x11'                ! iteration of temperature
  it(3) = 'x21'                ! iteration of mol fraction of comp.1 phase 2
  it(4) = 'x31'                ! iteration of mol fraction of comp.1 phase 3
  sum_rel(1) = 'x12'           ! summation relation: x22 = 1 - sum(x2j)
  sum_rel(2) = 'x22'           ! summation relation: x22 = 1 - sum(x2j)
  sum_rel(3) = 'x32'           ! summation relation: x32 = 1 - sum(x3j)

  val_init(0) = 1.E-8          ! set starting value for vapor density (3rd phase)
  val_init(1) = 0.436817011049269
  val_init(2) = 0.435497920638304
  val_init(3) = 308.850000000000
  val_init(4) = 1858758.08530477
  val_init(5) = -0.160891853903475
  val_init(6) = -1.90639042291856
  val_init(7) = -0.173218595064226
  val_init(8) = -1.83856034172500

  ! die folgende Definition ist eine voruebergehende Kruecke:
  xi(3,1) = 0.0001             ! set starting value for 3rd phase - temporary
  xi(3,2) = 1.0 - xi(3,1)
  lnx(3,1) = LOG(xi(3,1))
  lnx(3,2) = LOG(xi(3,2))
  val_init(9) = lnx(3,1)
  val_init(10)= lnx(3,2)

  loop_condition: DO

     CALL OBJECTIVE_CTRL (converg)

     IF (converg == 1) THEN
        CALL SI_DENS (density,w)
        write(*,*) ' '
        write(*,*) '   3-phase equilibrium calc. successful'
        write(*,*) ' '
        write(*,*) '  T,  P'
        write(*,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
        write(*,*) '      x11                  x21                x31'
        write(*,*) EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
        write(*,*) EXP(val_conv(6)),EXP(val_conv(8)),EXP(val_conv(10))
        write(*,*) '      eta1 ,               eta2 ,             eta3'
        write(*,*) dense(1),dense(2),dense(3)
        write(*,*) ' '
        write(*,*) ' output written to file: ./output_file/3-phase.xlo'

        OPEN (67,file='./output_file/3-phase.xlo')
        write(67,*) '  T,  P'
        write(67,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
        write(67,*) '      x11                x21                x31'
        write(67,*)EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
        write(67,*)EXP(val_conv(6)),EXP(val_conv(8)),EXP(val_conv(10))
        write(67,*) '   rho1[kg/m3]       rho2[kg/m3]       rho3[kg/m3]'
        write(67,*) density(1),density(2),density(3)
        write(67,'(6e16.8)') val_conv(3), val_conv(4)/u_out_P, &
             1.0-EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
     ENDIF

     IF (T.lt.309.8) THEN
        DO i=1,(4+ncomp*nphas)
           val_init(i)=val_conv(i)
        ENDDO
        val_init(3) = val_init(3) + 0.2
        T = val_init(3)
     ELSE
        EXIT loop_condition
     ENDIF
  END DO loop_condition

END SUBROUTINE LC_3PHASE




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE JT_INVERSION
!
! This sub. calculates Joule-Thomson curves for pure substances
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE JT_INVERSION

  USE PARAMETERS, ONLY: RGAS
  USE BASIC_VARIABLES
  use utilities, only: dens_calc, SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: state
  REAL                                   :: zdT = 10.0
  REAL                                   :: zdT2, zges1, zges2, zges3, start_t, delta
  REAL                                   :: density(np), w(np,nc)
  REAL                                   :: rho_phas(np)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 1

  IF (ncomp /= 1) THEN
     write(*,*)'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE:'
     write(*,*)'    ./input_file/INPUT.INP'
     stop
  ENDIF

  xi(1,1) = 1.0
  delta = 1.E-1
  start_t = t
  WRITE (*,*) ' SPECIFY LIQUID (1) or VAPOR (2) state'
  READ (*,*) state
  densta(1) = 0.5
  IF (state == 2) densta(1) = 1.E-5

  DO WHILE ( ABS(zdT) > 1.E-10 )

     t = start_t  -delta
     CALL DENS_CALC(rho_phas)
     CALL SI_DENS (density,w)
     zges1 = p/(density(1)*1000.0/mm(1)*RGAS*T)
     t = start_t  +delta
     CALL DENS_CALC(rho_phas)
     CALL SI_DENS (density,w)
     zges3 = p/(density(1)*1000.0/mm(1)*RGAS*T)
     t = start_t
     CALL DENS_CALC(rho_phas)
     CALL SI_DENS (density,w)
     zges2 = p/(density(1)*1000.0/mm(1)*RGAS*T)
     zdT   = (zges3-zges1)/(2.0*delta)
     zdT2  = (zges3-2.0*zges2+zges1)/delta**2
     start_t = t - zdT  /  zdT2
     write (*,*)  start_t ,p, zdT

  END DO
  write (40,*) 'T/K              p/kPa              error=dz_dT'
  write (40,*) start_t ,p/1000.0, zdT
  write (*,*) 'result written to: output_file/output.xlo'

END SUBROUTINE JT_INVERSION



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE BINMIX
!
! This subroutine performes VLE and LLE calculations for binary
! mixtures. Allowed are regular or polymeric substances. The routine
! is designed to derive phase equilibrium diagrams of binary systems.
! Three calculation options are available:
! option 1 (ISOTHERM): gives Pxy-diagram outputs. The pressure is
!     varied in 20 steps from the starting-pressure to an end-pressure
!     defined through by the user through the terminal. The step-size
!     in pressure is adjusted if the calculation approaches the
!     critical point (or an azeotropic point)
! option 2 (ISOBAR): gives Txy-diagram outputs. Similar to option 1,
!     whereby the temperature is varied.
! option 3 (MOLE FRACTION SCAN): Either, Txy- or Pxy-diagram outputs
!     can be selected this is a suitable option for azeotropic
!     mixtures.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BINMIX

  USE BASIC_VARIABLES
  USE STARTING_VALUES, only: START_VAR
  use utilities, only: select_sum_rel
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: iso_x, converg
  REAL                                   :: steps, end_x, pinit
  LOGICAL                                :: second
  !-----------------------------------------------------------------------------

  steps = 20.0
  second = .false.
  bindiag = 1

  CALL READ_INPUT
  nphas = 2
  pinit = p

  CALL START_VAR(converg)     ! gets starting values, sets "val_init"

  n_unkw = ncomp              ! number of quantities to be iterated 
  it(2) = 'x21'               ! iteration of mol fraction of comp.1 phase 2
  sum_rel(1) = 'x12'          ! summation relation: x12 = 1 - sum(x1j)
  sum_rel(2) = 'x22'          ! summation relation: x22 = 1 - sum(x2j)


  IF (converg == 1) THEN
     write (*,*) ' '
     write(*,*)'   CHOOSE YOUR CALCULATION OPTION:'
     write(*,*)' ISOTHERM (1) or ISOBAR (2) or MOL-FRACTION-SCAN (3)'
     READ (*,*) iso_x

     other_T: DO

        IF (iso_x == 1) THEN
           write (*,*) 'Specify end pressure  [units as in input-file]'
           READ (*,*) end_x
           end_x = end_x*u_in_P
           it(1) = 'x11'      ! iteration of mol fraction of comp.1 phase 1
           running = 'p'      ! Press. is running variable in PHASE_EQUILIB
           CALL select_sum_rel (1,0,1)
           CALL select_sum_rel (2,0,2)
        ELSEIF (iso_x == 2) THEN
           write (*,*) 'Specify end temperature [units as in input-file]'
           READ (*,*) end_x
           end_x = end_x + u_in_T
           it(1) = 'x11'      ! iteration of mol fraction of comp.1 phase 1
           running = 't'      ! Temp. is running variable in PHASE_EQUILIB
           CALL select_sum_rel (1,0,1)
           CALL select_sum_rel (2,0,2)
        ELSEIF (iso_x == 3) THEN
           write (*,*) ' YOU HAVE DECIDED TO VARY THE MOLEFRACTION'
           write (*,*) ' CHOOSE AN END-CONCENTRAION, FOR EXAMPLE:'
           write (*,*) ' CALCULATE TOWARDS xi(1) = 0              (0)'
           write (*,*) ' CALCULATE TOWARDS xi(1) = 1              (1)'
           write (*,*) ' COMPOUND #1 IS:  ',compna(1)
           READ (*,*) end_x
           IF (end_x <= 0.0) THEN
              end_x = 1.E-12
           ELSE
              end_x = 1.0 - 1.E-12
           ENDIF
           running = 'x11'               ! xi(1,1) is running var. in PHASE_EQUILIB
           steps = 80.0
           write(*,*)' ISOTHERM (1)  or  ISOBAR (2)'
           READ (*,*) iso_x
           IF (iso_x == 1) it(1) = 'p'   ! iteration of pressure
           IF (iso_x == 2) it(1) = 't'   ! iteration of temperature
        ENDIF

        outp = 1                      ! output to terminal
        CALL PHASE_EQUILIB(end_x,steps,converg)

        IF ( .NOT.second .AND. iso_x == 1 .AND. running /= 'x11') THEN
           second = .true.
           p = pinit
           val_conv = 0.0
           CALL START_VAR(converg)    ! gets starting values, sets "val_init"
           end_x = 10.0               ! End-Druck f¸r zweiten Durchlauf
           IF (eos >= 4) end_x = 1.E-3
           steps = 8.0
           CALL PHASE_EQUILIB(end_x,steps,converg)
        ENDIF

        IF (iso_x == 1) THEN
           write (*,*) ' '
           write (*,*) ' CALCULATE ANOTHER TEMPERATURE ?'
           write (*,*) ' CHOOSE (0 FOR NO, 1 FOR YES)'
           READ (*,*) iso_x
           IF (iso_x == 1) THEN
              WRITE (40,*) ' '
              xiF = 0.0
              p = pinit
              second =.false.
              write (*,*) 'Specify temperature [units as in input-file]'
              READ (*,*) t
              t = t + u_in_T
              val_conv = 0.0
              CALL START_VAR(converg)    ! gets starting values, sets "val_init"
              steps = 20.0
              it(2) = 'x21'              ! iteration of mol fraction of comp.1 phase 2
              sum_rel(1) = 'x12'         ! summation relation: x12 = 1 - sum(x1j)
              sum_rel(2) = 'x22'         ! summation relation: x22 = 1 - sum(x2j)
              CYCLE other_T
           ENDIF
        ENDIF
        EXIT other_T
     END DO other_T

  ENDIF

  write (*,*) ' '
  write (*,*) 'kij(1,2)',kij(1,2)

END SUBROUTINE BINMIX



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BINCRIT

  USE BASIC_VARIABLES
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, steps_i
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2

  steps_i = 20

  OPEN (71,file='./output_file/PT-crit.xlo')

  DO i = 1, steps_i-1
     xiF(2) = REAL( i ) / REAL( steps_i )
     xiF(1) = 1.0 - xiF(2)
     dense(1) = 0.15
     CALL Heidemann_Khalil
     write (71,'(8(2x,f18.8))') dense(1), t-u_out_T,p/u_out_P, &
          xiF(1:ncomp)
  ENDDO

  write (*,*) ' '
  write (*,*) 'kij(1,2)',kij(1,2)

END SUBROUTINE BINCRIT



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE HE_MIX
!
! This subroutine calculates physical properties of binary mixtures.
! The concentration ranges from one pure compound in 40 steps to
! the other pure substance.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE HE_MIX

  USE BASIC_VARIABLES
  use utilities, only: ENTHALPY_ETC, SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, i_steps
  REAL                                   :: density(np),w(np,nc)
  REAL                                   :: lnphi(np,nc)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT

  OPEN (68,file='./output_file/caloric.xlo')
  WRITE (68,*) 'phi1 phi2 cp_res[J/molK] h_res[J/mol] g_res[J/mol] rho t[K] p[Pa] xi(1)'

  IF (ncomp /= 2) THEN
     WRITE (*,*)'CALCULATION FOR BINARY MIXTURE ONLY!'
     STOP
  ENDIF
  nphas  = 1

  densta(1) = 0.5

  i_steps = 40

  DO i= 1,(i_steps+1)
     xi(1,1) = 1.0 - REAL(i-1)/REAL(i_steps)
     xi(1,2) = 1.0 - xi(1,1)

     CALL FUGACITY (lnphi)
     CALL ENTHALPY_ETC

     CALL SI_DENS (density,w)
     WRITE (*,*) cpres(1),enthal(1),dense(1),density(1),t,p,xi(1,1)
     WRITE (*,*) gibbs(1),xi(1,1)
     WRITE (68,'(9G15.7)') EXP(lnphi(1,1)),EXP(lnphi(1,2)), &
          cpres(1),enthal(1),gibbs(1),density(1),t,p,xi(1,1)

  END DO

END SUBROUTINE HE_MIX



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE STATE_PROP
!
! This subroutine calculates the physical properties of a one-phase
! (liquid) multicomp. mixture at given concentration.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE STATE_PROP

  USE BASIC_VARIABLES
  use utilities, only: ENTHALPY_ETC, SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: density(np), w(np,nc)
  REAL                                   :: lnphi(np,nc), aver_m
  !-----------------------------------------------------------------------------

  CALL READ_INPUT

  nphas  = 1
  densta(1) = 0.5

  WRITE (*,*) ' This subroutine calculates the physical properties'
  WRITE (*,*) ' of a one-phase (liquid) multicomponent mixture at'
  WRITE (*,*) ' given concentration'
  WRITE (*,*) ' '
  WRITE (*,'(a,i3,a)') ' SPECIFY MOL FRACTION OF ALL',ncomp,' COMPONENTS'
  DO i=1,ncomp
     READ (*,*) xi(1,i)
  END DO


  OPEN (68,file='./output_file/state_point.xlo')

  CALL FUGACITY (lnphi)
  CALL ENTHALPY_ETC

  CALL SI_DENS (density,w)

  aver_m = SUM( xi(1,1:ncomp)*mm(1:ncomp) )

  WRITE (*,*) cpres(1),enthal(1),dense(1),density(1),t-u_out_T, p/u_out_P
  WRITE (*,*) ' '
  WRITE (*,*) (xi(1,i), i = 1,ncomp)
  WRITE (68,*) 'cp_res[J/molK] h_res[J/mol] g_res[J/mol] rho t p '
  WRITE (68,*) cpres(1),enthal(1),gibbs(1),density(1),t-u_out_T, p/u_out_P
  WRITE (68,*) '  '
  WRITE (68,*) 'x(i) '
  WRITE (68,*) (compna(i),i=1,ncomp)
  WRITE (68,*) (xi(1,i),i=1,ncomp)
  WRITE (68,*) '  '
  WRITE (68,*) 'w(i) '
  WRITE (68,*) ((xi(1,i)*mm(i)/aver_m),i=1,ncomp)
  WRITE (68,*) '  '
  WRITE (68,*) 'average molar mass'
  WRITE (68,*) aver_m
  write (69,*) dense(1),dense(2),p

END SUBROUTINE STATE_PROP

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE P_LINE
!
! This subroutine calculates the physical properties of a one-phase
! (liquid) multicomp. mixture at given concentration.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_LINE

  USE BASIC_VARIABLES
  use utilities, only: P_CALC
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  REAL                                   :: pges,zges
  !-----------------------------------------------------------------------------

  CALL READ_INPUT

  nphas=1
  dense(1)= 0.0
  xi(1,1)  = 1.0
  write (*,*) t,p


  DO WHILE (dense(1) < 0.25)

     dense(1) = dense(1) + 0.01
     CALL P_CALC (pges,zges)
     write (71,*) dense(1),pges

  END DO

END SUBROUTINE P_LINE


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE P_T_DIAGRAM
!
! This subroutine calculates P-T-diagrams for binary and ternary
! mixtures.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE P_T_DIAGRAM

  USE BASIC_VARIABLES
  USE STARTING_VALUES, only: START_VAR
  use utilities, only: molefrac, SWITCH
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, iso_x, converg, swtch, ncompsave
  REAL                                   :: steps, end_x, end_x1, end_x2, startcon
  REAL                                   :: w(np,nc)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2
  ncompsave = ncomp
  ncomp = 2

  CALL START_VAR(converg) ! gets starting values, sets "val_init"
  ncomp = ncompsave

  n_unkw = ncomp       ! number of quantities to be iterated 

  write (*,*) ' '
  write (*,*) ' CHOOSE A MASS FRACTION FOR COMP. #1 :',converg
  write (*,*) ' COMPOUND #1 IS ',compna(1)
  READ (*,*) w(1,1)
  IF (ncomp == 2) THEN
     w(1,2) =1.0 - w(1,1)
     CALL MOLEFRAC (w,1,0)
     end_x1 = LOG(xi(1,1))

     write (*,*) ' '
     write(*,*)' FOR APPROACHING THE END CONCENTRATION FOR'
     write(*,*)' COMPONENT # 1 CHOOSE CALCULATION PATH: '
     write(*,*)' ISOTHERM (1) or ISOBAR (2)'
     READ (*,*) iso_x
     write (*,*) ' '
     write(*,*)' END CONC. LEFT OF CRIT. POINT (yes:1 / no:0)'
     READ (*,*) swtch
     IF (swtch == 1) THEN
        CALL SWITCH
        DO i=1,8
           val_conv(i)=0.0
        ENDDO
     ENDIF

     IF (iso_x == 1) it(1)='p'    ! iteration of pressure
     IF (iso_x == 2) it(1)='t'    ! iteration of temperature
     it(2)='x21'          ! iteration of mol fraction of comp.1 phase 2
     sum_rel(1)='x12'     ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2)='x22'     ! summation relation: x22 = 1 - sum(x2j)
     outp = 1             ! output to terminal
     running = 'l11'      ! lnx(1,1) is running var. in PHASE_EQUILIB
     steps=10.0
     CALL PHASE_EQUILIB(end_x1,steps,converg)

  ELSEIF (ncomp == 3) THEN
     write (*,*) ' '
     write (*,*) ' CHOOSE A MASS FRACTION FOR COMP. #3 :'
     write (*,*) ' COMPOUND #3 IS ',compna(3)
     READ (*,*) w(1,3)
     ! xi(1,2) =1.0 - xi(1,1) - xi(1,3)
     w(1,2) =1.0 - w(1,1) - w(1,3)
     CALL MOLEFRAC (w,1,0)
     end_x1 = LOG(xi(1,1))
     end_x2 = LOG(xi(1,3))

     ! if convergence problems are encountered, play around here:
     startcon = -12.0    ! starting value for lnx(1,3)
     steps=10.0

     val_init(5) = val_conv(5)
     val_init(6) = val_conv(6)
     val_init(7) = startcon
     val_init(8) = val_conv(7)
     val_init(9) = val_conv(8)
     ! just an attempt for an initial guess for val_init(10):
     val_init(10) = startcon + (end_x2-startcon)/steps
     write (*,*) ' '
     write(*,*)' END CONC. LEFT OF CRIT. POINT (yes:1 / no:0)'
     READ (*,*) swtch
     IF (swtch == 1) CALL SWITCH

     write (*,*) ' '
     write(*,*)' FOR APPROACHING THE END CONCENTRATION FOR'
     write(*,*)' COMPONENT # 3 CHOOSE CALCULATION PATH: '
     write(*,*)' ISOTHERM (1) or ISOBAR (2)'
     READ (*,*) iso_x
     IF (iso_x == 1) it(1)='p'    ! iteration of pressure
     IF (iso_x == 2) it(1)='t'    ! iteration of temperature
     it(2)='x21'        ! iteration of mol fraction of comp.1 phase 2
     it(3)='x23'        ! iteration of mol fraction of comp.3 phase 2
     sum_rel(1)='x12'   ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2)='x22'   ! summation relation: x22 = 1 - sum(x2j)
     outp = 1           ! output to terminal
     running='l13'      ! lnx(1,3) is running var. in PHASE_EQUILIB

     CALL PHASE_EQUILIB(end_x2,steps,converg)


     write (*,*) ' '
     write(*,*)' FOR APPROACHING THE END CONCENTRATION FOR'
     write(*,*)' COMPONENT # 1 CHOOSE CALCULATION PATH: '
     write(*,*)' ISOTHERM (1) or ISOBAR (2)'
     READ (*,*) iso_x
     IF (iso_x == 1) it(1)='p'    ! iteration of pressure
     IF (iso_x == 2) it(1)='t'    ! iteration of temperature
     it(2)='x21'        ! iteration of mol fraction of comp.1 phase 2
     it(3)='x23'        ! iteration of mol fraction of comp.3 phase 2
     sum_rel(1)='x12'   ! summation relation: x12 = 1 - sum(x1j)
     sum_rel(2)='x22'   ! summation relation: x22 = 1 - sum(x2j)
     outp = 1           ! output to terminal
     running='l11'      ! lnx(1,1) is running var. in PHASE_EQUILIB
     steps=10.0
     CALL PHASE_EQUILIB(end_x1,steps,converg)

  ENDIF


  !  start calculation of P-T-diagram
  write(*,*)' '
  write(*,*) ' SPECIFY TEMPERATURE (1)'
  write(*,*) ' OR PRESSURE         (2)'
  READ (*,*) iso_x
  IF (iso_x == 1) THEN
     write(*,*) '  CHOOSE END TEMPERATURE [units as in input-file]:'
     READ (*,*) end_x
     end_x = end_x + u_in_T
     it(1)='p'    ! iteration of pressure
     running='t'          ! Temperature is running var. in PHASE_EQUILIB
  ELSEIF (iso_x == 2) THEN
     write(*,*) '  CHOOSE END PRESSURE  [units as in input-file]:'
     READ (*,*) end_x
     end_x = end_x*u_in_P
     it(1)='t'    ! iteration of temperature
     running='p'          ! Temperature is running var. in PHASE_EQUILIB
  ENDIF
  outp = 1             ! output to terminal
  steps=20.0
  CALL PHASE_EQUILIB(end_x,steps,converg)

END SUBROUTINE P_T_DIAGRAM


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!  SUBROUTINE TERN_DIAGRAM
!
! This subroutine calculates ternary diagrams for given (T and P)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE TERN_DIAGRAM

  USE BASIC_VARIABLES
  USE STARTING_VALUES, only: START_VAR
  use utilities, only: molefrac, SWITCH
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: converg, spec_ph, spec_comp
  REAL                                   :: steps, end_x
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2

  CALL START_VAR(converg) ! gets starting values, sets "val_init"

  n_unkw = ncomp       ! number of quantities to be iterated 

  write (*,*) ' '
  write (*,*) ' THE MOLE FRACTION OF OF PHASE 1 OR OF PHASE 2'
  write (*,*) ' CAN NOW BE CHANGED. SELECT:'
  write (*,*) ' PHASE 1          (1)'
  write (*,*) ' PHASE 2          (2)'
  !      write (*,*) ' FEED-PHASE       (3)'
  READ (*,*) spec_ph
  write (*,*) ' '
  write (*,*) ' SELECT A COMPONENT TO BE SPECIFIED'
  write (*,*) ' COMPONENT #1 ',compna(1),'  (1)'
  write (*,*) ' COMPONENT #2 ',compna(2),'  (2)'
  write (*,*) ' COMPONENT #3 ',compna(3),'  (3)'
  READ (*,*) spec_comp

  write (*,*) ' '
  CALL OUTPUT
  write (*,*) ' '
  write (*,*) ' CHOOSE AN END MOLE FRACTION FOR COMP. #',spec_comp
  write (*,*) ' (CAUTION: THREE-PHASE EQUIL. ARE NOT DETECTED)'
  READ (*,*) end_x
  end_x = LOG(end_x)

  steps = 10.0

  it(1) = 'x11'        ! iteration of mol fraction of comp.1 phase 1
  it(2) = 'x21'        ! iteration of mol fraction of comp.1 phase 2
  it(3) = 'x23'        ! iteration of mol fraction of comp.3 phase 2
  sum_rel(1) = 'x12'   ! summation relation: x12 = 1 - sum(x1j)
  sum_rel(2) = 'x22'   ! summation relation: x22 = 1 - sum(x2j)

  IF (spec_ph == 1) THEN
     IF (spec_comp == 1) running='l11' ! lnx(1,1) is running var. in PHASE_EQUILIB
     IF (spec_comp == 1) it(1)='x13'   ! iteration of mol fraction of comp.3 phase 1
     IF (spec_comp == 2) running='l12' ! lnx(1,2) is running var. in PHASE_EQUILIB
     IF (spec_comp == 2) sum_rel(1)='x13' ! summation relation: x12 = 1 - sum(x1j)
     IF (spec_comp == 3) running='l13' ! lnx(1,3) is running var. in PHASE_EQUILIB
  ELSEIF (spec_ph == 2) THEN
     IF (spec_comp == 1) running='l21' ! lnx(2,1) is running var. in PHASE_EQUILIB
     IF (spec_comp == 1) it(2)='x13'   ! iteration of mol fraction of comp.3 phase 1
     IF (spec_comp == 2) running='l22' ! lnx(2,2) is running var. in PHASE_EQUILIB
     IF (spec_comp == 2) it(2)='x13'   ! iteration of mol fraction of comp.3 phase 1
     IF (spec_comp == 2) sum_rel(2)='x21' ! summation relation: x21 = 1 - sum(x2j)
     IF (spec_comp == 3) running='l23' ! lnx(2,3) is running var. in PHASE_EQUILIB
     IF (spec_comp == 3) it(2)='x13'   ! iteration of mol fraction of comp.3 phase 1
     IF (spec_comp == 3) it(3)='x21'   ! iteration of mol fraction of comp.1 phase 2
  ELSEIF (spec_ph == 3) THEN
     ! it(1)='fls'         ! do a flash-calcul. for fixed feed-conc
     ! IF (spec_comp == 1) running='xF1'
     ! IF (spec_comp == 2) running='xF2'
     ! IF (spec_comp == 3) running='xF3'
  ENDIF
  outp = 1           ! output to terminal

  CALL PHASE_EQUILIB(end_x,steps,converg)

END SUBROUTINE TERN_DIAGRAM


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE BIN_3PHASE
!
! This subroutine performes VLLE calculations for binary mixtures
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE BIN_3PHASE

  USE BASIC_VARIABLES
  USE STARTING_VALUES, only: START_VAR
  use utilities, only: SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: k, ph, iso_x, converg, i
  REAL                                   :: density(np), w(np,nc)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 2

  !---  determine liquid-liquid equilibria first
  CALL START_VAR(converg) ! gets starting values, sets "val_init"

  IF (converg == 1) THEN
     write (*,*) ' '
     write (*,*) '   LLE detected for starting values'
     write (*,*) ' '
     write (*,*) '   starting 3-phase calculation'
     write (*,*) ' '
  ENDIF

  write(*,*)'   CHOOSE YOUR CALCULATION OPTION:'
  write(*,*)'   ISOTHERM (1) or ISOBAR (2)'
  READ (*,*) iso_x
  nphas = 3
  n_unkw = 4           ! number of quantities to be iterated 
  it(1)='x11'          ! iteration of mol fraction of comp.1 phase 1
  it(2)='x21'          ! iteration of mol fraction of comp.1 phase 2
  it(3)='x31'          ! iteration of mol fraction of comp.1 phase 3
  IF (iso_x == 1) it(4)='p  '          ! iteration of pressure
  IF (iso_x == 2) it(4)='t  '          ! iteration of temperature
  sum_rel(1)='x12'     ! summation relation: x12 = 1 - sum(x1j)
  sum_rel(2)='x22'     ! summation relation: x22 = 1 - sum(x2j)
  sum_rel(3)='x32'     ! summation relation: x32 = 1 - sum(x3j)

  ! die folgende Definition ist eine voruebergehende Kruecke:
  xi(3,1) = 0.01       ! set starting value for 3rd phase - temporary
  xi(3,2) = 1.0 - xi(3,1)
  lnx(3,1) = LOG(xi(3,1))
  lnx(3,2) = LOG(xi(3,2))

  val_init(0) = 1.E-8  ! set starting value for vapor density (3rd phase)

  DO ph = 1,nphas
     DO k = 1,ncomp
        val_init(4+k+(ph-1)*ncomp) = lnx(ph,k)
     ENDDO
  ENDDO

  other_condition: DO

     CALL OBJECTIVE_CTRL (converg)

     IF (converg == 1) THEN
        CALL SI_DENS (density,w)
        write(*,*) ' '
        write(*,*) '   3-phase equilibrium calc. successful'
        write(*,*) ' '
        write(*,*) '  T,  P'
        write(*,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
        write(*,*) '      x11                  x21                x31'
        write(*,*) EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
        write(*,*) EXP(val_conv(6)),EXP(val_conv(8)),EXP(val_conv(10))
        write(*,*) '      eta1 ,               eta2 ,             eta3'
        write(*,*) dense(1),dense(2),dense(3)
        write(*,*) ' '
        write(*,*) ' output written to file: ./output_file/3-phase.xlo'

        OPEN (67,file='./output_file/3-phase.xlo')
        write(67,*) '  T,  P'
        write(67,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
        write(67,*) '      x11                x21                x31'
        write(67,*)EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
        write(67,*)EXP(val_conv(6)),EXP(val_conv(8)),EXP(val_conv(10))
        write(67,*) '   rho1[kg/m3]       rho2[kg/m3]       rho3[kg/m3]'
        write(67,*) density(1),density(2),density(3)
        write(67,'(6e16.8)') val_conv(3), val_conv(4)/u_out_P, &
             1.0-EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
     ENDIF

     IF (T.lt.210.15) THEN
        DO i=1,(4+ncomp*nphas)
           val_init(i)=val_conv(i)
        ENDDO
        val_init(3) = val_init(3) + 0.025
        T = val_init(3)
     ELSE
        EXIT other_condition
     ENDIF

  END DO other_condition

END SUBROUTINE BIN_3PHASE


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE TERN_3PHASE
!
! This subroutine performes VLLE calculations for binary mixtures
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE TERN_3PHASE

  USE BASIC_VARIABLES
  use utilities, only: SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: k, ph, converg
  REAL                                   :: density(np), w(np,nc)
  !-----------------------------------------------------------------------------

  CALL READ_INPUT
  nphas = 3

  lnx(1,1) = -4.375218662
  lnx(1,2) = -0.241363748
  lnx(1,3) = -1.600186935
  lnx(2,1) = -14.98669354
  lnx(2,2) = -0.395440676
  lnx(2,3) = -1.118968702
  lnx(3,1) = -30.0
  lnx(3,2) = -0.821288895
  lnx(3,3) = -0.579576292

  t = 449.0
  p = 16.E6
  val_init(3) = t
  val_init(4) = p

  val_init(1) = 0.4  ! set starting value for vapor density (1st phase)
  val_init(2) = 0.4  ! set starting value for vapor density (2nd phase)
  val_init(0) = 1.E-5  ! set starting value for vapor density (3rd phase)

  n_unkw = 6           ! number of quantities to be iterated 
  it(1)='x21'          ! iteration of mol fraction of comp.1 phase 2
  it(2)='x22'          ! iteration of mol fraction of comp.2 phase 2
  it(3)='x31'          ! iteration of mol fraction of comp.1 phase 3
  it(4)='x32'          ! iteration of mol fraction of comp.2 phase 3
  it(5)='t  '          ! temp
  it(6)='p  '          ! pressure
  sum_rel(1)='x13'     ! summation relation: x12 = 1 - sum(x1j)
  sum_rel(2)='x23'     ! summation relation: x22 = 1 - sum(x2j)
  sum_rel(3)='x33'     ! summation relation: x32 = 1 - sum(x3j)


  DO ph = 1,nphas
     DO k = 1,ncomp
        val_init(4+k+(ph-1)*ncomp) = lnx(ph,k)
     ENDDO
  ENDDO

  ! 10 continue

  CALL OBJECTIVE_CTRL (converg)

  IF (converg == 1) THEN
     CALL SI_DENS (density,w)
     write(*,*) ' '
     write(*,*) '   3-phase equilibrium calc. successful'
     write(*,*) ' '
     write(*,*) '  T,  P'
     write(*,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
     write(*,*) '  x_phase1       x_phase2         x_phase3'
     write(*,*)EXP(val_conv(5)),EXP(val_conv(8)),EXP(val_conv(11))
     write(*,*)EXP(val_conv(6)),EXP(val_conv(9)),EXP(val_conv(12))
     write(*,*)EXP(val_conv(7)),EXP(val_conv(10)),EXP(val_conv(13))
     write(*,*) '  w_phase1       w_phase2         w_phase3'
     write(*,*)w(1,1),w(2,1),w(3,1)
     write(*,*)w(1,2),w(2,2),w(3,2)
     write(*,*)w(1,3),w(2,3),w(3,3)
     write(*,*) '      eta1 ,               eta2 ,             eta3'
     write(*,*) dense(1),dense(2),dense(3)
     write(*,*) ' '
     write(*,*) ' output written to file: ./output_file/3-phase.xlo'

     OPEN (67,file='./output_file/3-phase.xlo')
     write(67,*) '  T,  P'
     write(67,*) val_conv(3)-u_out_T,val_conv(4)/u_out_P
     write(67,*) '  x_phase1       x_phase2         x_phase3'
     write(67,*)EXP(val_conv(5)),EXP(val_conv(8)),EXP(val_conv(11))
     write(67,*)EXP(val_conv(6)),EXP(val_conv(9)),EXP(val_conv(12))
     write(67,*)EXP(val_conv(7)),EXP(val_conv(10)),EXP(val_conv(13))
     write(67,*) '   rho1[kg/m3]       rho2[kg/m3]       rho3[kg/m3]'
     write(67,*) density(1),density(2),density(3)
     write(67,'(6e16.8)') val_conv(3), val_conv(4)/u_out_P, &
          1.0-EXP(val_conv(5)),EXP(val_conv(7)),EXP(val_conv(9))
  ENDIF

  write (*,*) ' '
  write (*,*) parame(1,13),parame(1,16)
  write (*,*) parame(1,3),kij(1,2),kij(1,3)
  write (*,*) 'modify parameter eps(1)/k'
  read(*,*)   parame(1,3)
  write (*,*) 'modify parameter eps_assoc'
  read(*,*)   parame(1,13)
  ! parame(1,16) = parame(1,15)
  write (*,*) 'modify parameter kij(1,2)'
  read(*,*)   kij(1,2)
  kij(2,1) = kij(1,2)
  write (*,*) 'modify parameter kij(1,3)'
  read(*,*)   kij(1,3)
  kij(3,1) = kij(1,3)

  ! IF (t.lt.473.15) THEN
  ! DO i=1,(4+ncomp*nphas)
  !   val_init(i)=val_conv(i)
  ! ENDDO
  ! val_init(3) = val_init(3) + 5.0
  ! t = val_init(3)
  ! GOTO 10
  ! ENDIF

END SUBROUTINE TERN_3PHASE


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PARITER_CRIT
!
! This subroutine calculates the three pure component parameters of
! PC-SAFT for non-polar and non-associating components from the
! critical temp., the critical pressur and the acentric factor
! OMEGA.
! This procedure sacrifices the density in the liquid phase and
! sacrifices therefore also caloric properties, like cp and enthalpy
! (of vaporization).
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE PARITER_CRIT

  USE BASIC_VARIABLES
  USE utilities, only: CRITICAL, paus
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i, j, converg
  REAL                                   :: sig_new, eps_new, fvek1, fvek2
  REAL                                   :: sig_save, eps_save
  REAL                                   :: df1dx, df2dx, df1dy, df2dy
  REAL                                   :: betrag, attenu
  REAL                                   :: tc, pc, rhoc, omega, seg_new, par1_sv, par2_sv
  REAL                                   :: par3_sv, omeg_sv, omeg_nw, err_sv, err_nw, derrdr
  REAL                                   :: steps, end_x
  !-----------------------------------------------------------------------------


  OPEN (10,file = './output_file/crit-par.out')

  write (*,*) ' '
  write (*,*) ' ENTER MOLECULAR MASS         M  / g/mol :'
  READ (*,*) mm(1)
  ! mm(1)=16.043
  write (*,*) ' ENTER CRITICAL TEMPERATURE   Tc / K    :'
  READ (*,*) tc
  ! tc = - 82.59
  write (*,*) ' ENTER CRITICAL PRESSURE      Pc / bar   :'
  READ (*,*) pc
  ! pc = 45.95
  write (*,*) ' ENTER ACENTRIC FACTOR        omega      :'
  READ (*,*) omega
  ! omega = 0.01
  ! tc = tc + 273.15
  pc = pc * 1.E5

  eos   = 1
  ncomp = 1
  nphas = 2

  !----------------------------------------------------------
  sig_new = 4.2
  eps_new = 240.0
  seg_new = mm(1)*0.05
  !----------------------------------------------------------
  parame(1,4) = 0.0
  parame(1,5) = 0.12
  parame(1,6) = 0.0
  parame(1,7) = 0.0
  DO j=8,25
     parame(1,j) = 0.0
  ENDDO
  scaling(1)   = 1.E6

  attenu = 1.0

  crit_loop: DO

     DO i = 1,2

        j = 1
        DO WHILE ( (ABS(fvek1)+ABS(fvek2)) > 2.6E-4 )
           j = j + 1

           parame(1,2) = sig_new
           parame(1,3) = eps_new

           !------------------------------
           xi(1,1)      = 1.0
           parame(1,1) = seg_new + 0.01*REAL(i-1)
           t  = tc
           p  = pc

           CALL critical (tc,pc,rhoc)

           fvek1= (t-tc)/tc
           fvek2= (p-pc)/pc

           ! Bildung der numerischen Ableitung nach sigma:
           sig_save = parame(1,2)
           parame(1,2) = sig_save*1.0002

           CALL critical (tc,pc,rhoc)

           df1dx = (fvek1-(t-tc)/tc )/( sig_save-parame(1,2) )
           df2dx = (fvek2-(p-pc)/pc )/( sig_save-parame(1,2) )
           parame(1,2) = sig_save



           ! Bildung der numerischen Ableitung nach epsilon:
           eps_save = parame(1,3)
           parame(1,3) = eps_save*1.0002

           CALL critical (tc,pc,rhoc)

           df1dy = (fvek1-(t-tc)/tc )/( eps_save-parame(1,3) )
           df2dy = (fvek2-(p-pc)/pc )/( eps_save-parame(1,3) )
           parame(1,3) = eps_save

           betrag = df1dx*df2dy - df1dy*df2dx
           sig_new = parame(1,2) - attenu / betrag * (df2dy*fvek1+(-df1dy)*fvek2)
           eps_new = parame(1,3) - attenu / betrag * ((-df2dx)*fvek1+df1dx*fvek2)

           write (10,*) parame(1,2),parame(1,3),fvek1,fvek2,i
           write (*,*)  parame(1,2),parame(1,3),fvek1,fvek2,i
           ! if (j == 15) stop

           ! IF (((ABS(fvek1) > 1.E-4).OR.(ABS(fvek2) > 1.E-4)).
        END DO


        ! --- calculate vapor pressure at 0.7*tc (for acentric factor)
        n_unkw = ncomp         ! number of quantities to be iterated 
        it(1) = 'p'            ! iteration of mol fraction of comp.1 phase 1
        running = 't'          ! Temp. is running variable in PHASE_EQUILIB

        val_init(1) = 0.5
        val_init(2) = 1.E-6
        val_init(3) = 0.7*tc
        val_init(4) = 1.E+2    ! default starting value for P: 100 Pa
        val_init(5) = 0.0      ! logar. mole fraction: lnx=0, x=1, phase 1
        val_init(6) = 0.0      ! logar. mole fraction: lnx=0, x=1, phase 2

        end_x = 0.7*tc

        steps=1.0
        outp = 0               ! output to terminal
        CALL PHASE_EQUILIB(end_x,steps,converg)
        IF (converg /= 1) THEN
           call paus ('PARITER_CRIT: equilibrium calculation failed')
        ENDIF

        IF (i == 1) THEN
           par1_sv = parame(1,2)
           par2_sv = parame(1,3)
           par3_sv = parame(1,1)
           omeg_sv = -1.0 - LOG10(val_conv(4)/pc)
           err_sv  = omeg_sv-omega
           write (10,*) 'error',i,err_sv
           write (*,*) 'error',i,err_sv
        ELSE
           omeg_nw = -1.0 - LOG10(val_conv(4)/pc)
           err_nw  = omeg_nw-omega
           write (10,*) 'error',i,err_nw
           write (*,*) 'error',i,err_nw
        ENDIF
        write (10,*) 'omega',i,-1.0 - LOG10(val_conv(4)/pc)
        write (*,*) 'omega',i,-1.0 - LOG10(val_conv(4)/pc)

     END DO

     derrdr =(err_sv - err_nw)/(par3_sv-parame(1,1))
     seg_new = par3_sv - 0.9 * err_sv/derrdr
     write (10,*) 'derivative',derrdr
     write (10,*) 'new value for m',seg_new
     write (*,*) 'derivative',derrdr
     write (*,*) 'new value for m',seg_new

     IF (ABS(err_sv) <= 3.E-5) EXIT crit_loop

  END DO crit_loop

  write (*,*) '  '
  write (*,*) '  Pure component parameters :'
  write (*,*) '  '
  write (*,*) '     sigma = ',par1_sv
  write (*,*) '     eps/k = ',par2_sv
  write (*,*) '     m     = ',par3_sv
  write (10,*) '  '
  write (10,*) '  Pure component parameters :'
  write (10,*) '  '
  write (10,*) '        mm(i)       = ',mm(1)
  write (10,*) '        parame(i,1) = ',par3_sv
  write (10,*) '        parame(i,2) = ',par1_sv
  write (10,*) '        parame(i,3) = ',par2_sv

  CLOSE (10)

END SUBROUTINE PARITER_CRIT


