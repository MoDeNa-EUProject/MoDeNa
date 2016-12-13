!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE DFT

  USE PARAMETERS, ONLY: RGAS, PI
  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  USE DFT_MODULE
  use utilities
  USE EOS_NUMERICAL_DERIVATIVES, ONLY: subtract1
  IMPLICIT NONE

  LOGICAL diagram
  INTEGER                                :: i, j, kk, converg, discret, fa, outpt, irc
  INTEGER                                :: ITHICK, IDEL, IAB1, IAB2         ! ,stepno
  REAL, DIMENSION(-NDFT:NDFT)            :: zp
  REAL, DIMENSION(-NDFT:NDFT)            :: rhop, rhop_o
  REAL, DIMENSION(-NDFT:NDFT)            :: mf_att
  REAL ,DIMENSION(2)                     :: rhob
  REAL                                   :: end_x, steps, my0
  REAL                                   :: zz1
  REAL                                   :: integrand, dev
  REAL                                   :: omega, free, phs1, om1, pbulk, surftens
  REAL                                   :: zges(np), inte2, f01, rm
  REAL                                   :: tc, pc, rhoc, tsav     ! ,attan ,myrho(-5000:5000)
  REAL                                   :: density(np), w(np,nc), lnphi(np,nc)
  REAL                                   :: cploc, DENS_INV, DUMMY, CPCP
  REAL                                   :: intgrid(0:5000)  ! ,utri(5000),intf(0:5000)
  !-----------------------------------------------------------------------------

  num = 1
  CALL SET_DEFAULT_EOS_NUMERICAL

  WCA = .false.
  diagram =.true.

  CALL READ_INPUT

  ! T_st=t/parame(1,3)
  ! IF (WCA)      CALL Barker(T_st,d_hs)
  ! IF (.NOT.WCA) CALL Barker_org(T_st,d_hs)
  z3t = PI/6.0* parame(1,1) * ( 1.0  &  ! divided by parame(i,2)**3
       *(1.0-0.12*EXP(-3.0*parame(1,3)/t)) )**3

  IF (ncomp /= 1) THEN
     write(*,*)'SPECIFY ONLY ONE COMPONENT IN THE INPUT-FILE:'
     write(*,*)'    ./input_file/INPUT.INP'
     stop
  ENDIF

  !-----------------------------------------------------------------------------
  ! The truncated and shifted LJ-potential is considered
  ! The attractive part of this potential is defined by:
  !    u(r) = - epsilon                  - TAU  for  0 < r < rm
  !    u(r) = 4*epsilon*(r**-12 - r**-6) - TAU  for  rm< r < rc
  !    u(r) = 0.                                for      r > rc
  ! rm=2**(1/6) corresponds to the minimum of full LJ-potential
  ! rc=7.5 is the cut-off distance (both values dimensionless)
  ! (-TAU) is the value by which the potetial is shifted upwards,
  ! with TAU=4*epsilon*(rc**-12 - rc**-6)
  !-----------------------------------------------------------------------------
  rc = 4.0
  rm = 2.0**(1.0/6.0)
  tau_cut = 4.0*(rc**(-12)-rc**(-6))

  !-----------------------------------------------------------------------------
  ! basic definitions for calculating density profile,
  !-----------------------------------------------------------------------------
  discret= 1200                ! number of spacial discretizations for profile
  fa     = 40                  ! size of the box is discret/fa * sigma
  ! (with sigma: LJ diameter)
  ! fa is here dimensionless, however
  dzp     = 1.0/REAL(fa)       ! dimensionless grid-distance dzp*=dzp/sigma
  ! grid-size  dzp = zp(1)-zp(0)    
  outpt  = 10                  ! number output-points = discret/outpt

  !-----------------------------------------------------------------------------
  ! prepare for phase equilibrium calculation for given T
  !-----------------------------------------------------------------------------
  diagram_condition: DO
     nphas = 2
     n_unkw = ncomp            ! number of quantities to be iterated 
     it(1)='p'                 ! iteration of pressure

     val_init(1) = 0.45        ! starting value for liquid density
     val_init(2) = 1.E-6       ! starting value for vapor density
     val_init(3) = t           ! value of temperature
     val_init(4) = 1.0         ! default starting value for P
     IF (eos >= 4) val_init(4) = 1.E-3  ! default starting value for LJ-model
     val_init(5) = 0.0         ! logar. mole fraction: lnx=0, x=1, phase 1
     val_init(6) = 0.0         ! logar. mole fraction: lnx=0, x=1, phase 2

     running = 't'             ! T is running variable in PHASE_EQUILIB - here T=const.
     end_x  = t                ! end temperature - here T=const.
     steps = 1.0               ! number of steps of running var. - here 1, since T=const.

     outp = 0                  ! output to terminal
     !--------------------------------------------------------------------------
     ! calculate phase equilibrium for given T
     !--------------------------------------------------------------------------
     find_equil: DO
        ensemble_flag = 'tp'             ! this flag is for: 'regular calculation'
        CALL PHASE_EQUILIB(end_x,steps,converg)
        IF (converg == 1) THEN
           EXIT find_equil
        ELSE
           val_init(3) = t - 10.0        ! value of temperature
           steps = 2.0                   ! number of steps of running var. - here 1, since T=const.
           write (*,*) 'no VLE found'
        ENDIF
     END DO find_equil
     CALL SI_DENS (density,w)
     rhob(1) = dense(1)/z3t              ! coexisting bulk density L
     rhob(2) = dense(2)/z3t              ! coexisting bulk density V
     write (*,*) 'densities       ',rhob(1),rhob(2)

     ! (re-)calculate residual chemical potential of both phases
     ensemble_flag = 'tv'                ! this flag is for: my_res=my_res(T,rho)
     CALL FUGACITY (lnphi)
     my0 = lnphi(1,1)
     zges(1) = p / (RGAS *T *density(1)) * mm(1)/1000.0
     zges(2) = p / (RGAS *T *density(2)) * mm(1)/1000.0

     pbulk = zges(1)*rhob(1)             ! pressure  p/kT (= Z*rho)
     write (*,*) 'chem. potentials',lnphi(1,1),lnphi(2,1)
     ! write (*,*) 'pressures  p/kT ',zges(1)*rhob(1),zges(2)*rhob(2)
     ! write (*,*) 'pressures  p   ',p
     write (*,*) ' '

     !--------------------------------------------------------------------------
     ! define initial density profile rhop(i)
     ! and dimensionless space coordinates zp(i)
     !
     ! discret  : number of grid-points within "the box"
     ! irc      : number of grid-points extending the the box to left and right
     !            the box is extended in order to allow for numerical integration
     !            'irc' is determined by the cut-off distance 'rc'
     ! IAB1,IAB2: help-integers allowing a linear initial density profile in
     !            vicinity of the box-center
     !--------------------------------------------------------------------------
     irc = NINT(rc/dzp) + 1
     ITHICK = discret / 2
     idel   = discret / 10
     ITHICK = ITHICK + 1
     IAB2 = ITHICK + idel/2
     IAB1 = ITHICK - idel/2
     do i = -irc, iab1
        zp(i) = REAL(i) * dzp
        rhop(i) = rhob(1)
     END do
     do i = iab1+1, iab2
        zp(i) = REAL(i) * dzp
        rhop(i) = rhob(1) - (rhob(1) - rhob(2)) *REAL(I-IAB1) / REAL(IAB2-IAB1)
     end do
     do i = iab2, discret+irc
        zp(i) = REAL(i) * dzp
        rhop(i) =rhob(2)
     end do

     !--------------------------------------------------------------------------
     ! Initialize the DENS_INV subroutine
     !--------------------------------------------------------------------------
     nphas = 1

     CPCP = 1.E7
     ensemble_flag = 'tv'                ! this flag is for: my_res=my_res(T,rho)
     subtract1 = '1PT'                   ! this flag is for: 'minus 1st order pert. theory' for chem. pot.

     DUMMY = DENS_INV(CPCP)


     !--------------------------------------------------------------------------
     ! Start iterating the density profile
     !--------------------------------------------------------------------------
     kk=1
     converge_profile: DO

        DO i=0,discret
           ! dense(1)  = rhop(i)*z3t
           ! densta(1) = rhop(i)*z3t

           ! -------------------------------------------------------------------
           ! the residual chemical pot. that is treated with a LDA is here calculated.
           ! this part together with 'my_att' gives the functional derivative of Fres to rho(z)
           ! -------------------------------------------------------------------
           ! CALL FUGACITY (lnphi)
           ! myrho(i) = lnphi(1,1)

           ! inte = 0.0
           inte2 = 0.0
           f01 = 0.0
           DO j=-irc,discret+irc
              zz1 = ABS(zp(j)-zp(i))
              IF (WCA) THEN
                 IF (zz1 > rm .AND. zz1 <= rc)THEN
                    integrand= 0.4*(zz1**(-10)-rc**(-10)) &
                         -(zz1**(-4 )-rc**(-4 )) &
                         +tau_cut/2.0*(zz1**2 -rc**2 )
                 ELSEIF (zz1 <= rm) THEN
                    integrand= 0.40*(rm**(-10)- rc**(-10)) &
                         - (rm**(-4) - rc**(-4)) &
                         +0.5*(zz1**2- rm**2 ) &
                         +tau_cut/2.0*(zz1**2 -rc**2 )
                 ELSE
                    integrand= 0.0
                 ENDIF
              ELSE
                 IF (zz1 > 1.0 .AND. zz1 <= rc) THEN
                    integrand= 0.4*(zz1**(-10)-rc**(-10)) &
                         -(zz1**(-4 )-rc**(-4 )) &
                         +tau_cut/2.0*(zz1**2-rc**2)
                 ELSEIF (zz1 <= 1.0) THEN
                    integrand= 0.4*(1.0- rc**(-10)) &
                         - (1.0 - rc**(-4)) &
                         +tau_cut/2.0*(1.0-rc**2)
                 ELSE
                    integrand= 0.0
                 ENDIF
              ENDIF
              intgrid(j+irc) = rhop(j)* integrand
              inte2=inte2+dzp*(intgrid(j+irc)+f01)/2.0
              f01 = intgrid(j+irc)
           ENDDO

           !--------------------------------------------------------------------
           ! The integration of the DFT-integral is done with cubic splines
           !--------------------------------------------------------------------
           ! stepno = discret+irc*2
           ! CALL SPLINE_PARA (dzp,intgrid,utri,stepno)
           ! CALL SPLINE_INT (inte,dzp,intgrid,utri,stepno)


           !--------------------------------------------------------------------
           ! the attractive part of the chemical potential
           !      mf_att(i) = d(F_att)/d(rho*)
           ! where F_att is attractive part of the intrinsic Helmholtz energy
           ! and where rho*=rhop(i) denotes the dimensionless density.
           ! parame(1,3) is epsilon/k (LJ-energy parameter)
           ! parame(1,2) is sigma (LJ-size parameter)
           !--------------------------------------------------------------------
           mf_att(i) = -2.0*PI*parame(1,3)/t*parame(1,1)**2 *inte2
           ! if (REAL(i)/REAL(outpt) == REAL(INT(i/outpt)))  &
           ! write(*,*) i,rhop(i), -( myrho(i)-mf_att(i)-my0 ),inte/rhop(i)
        ENDDO


        dev = 0.0
        DO i=0,discret
           rhop_o(i) = rhop(i)
           !--------------------------------------------------------------------
           ! LDA - Inversion Procedure according to R.Evans, Bristol
           !--------------------------------------------------------------------
           cploc   = my0 + mf_att(i)
           rhop(i)   = DENS_INV(cploc)
           ! write (*,*) i,rhop_o(i),rhop(i)
           ! read (*,*)
           ! attan = 0.5
           ! rhop(i) = rhop_o(i)*(1.0-attan)+rhop(i)*attan
           dev = dev + (rhop(i)-rhop_o(i))**2
        ENDDO


        kk=kk+1
        write (*,*) kk,dev

        IF (dev >= 4.E-8 .AND. kk <= 50) CYCLE converge_profile
        IF (kk > 50) write (*,*) ' no convergence in 50 steps'
        EXIT converge_profile
     END DO converge_profile

     !--------------------------------------------------------------------------
     ! write resulting density profile
     !--------------------------------------------------------------------------
     write (68,'(i3,66(f18.8))') 0,zp(0)-zp(INT(discret/2)),rhop(0),T, &
          val_conv(4)/1.E5,density(1),density(2)
     DO i=0,discret
        IF (REAL(i)/REAL(outpt) == REAL(INT(i/outpt))) &
             write (68,'(i6,3(f18.10))') i,zp(i)-zp(INT(discret/2)),rhop(i)
        ! IF (REAL(i)/REAL(outpt) == REAL(INT(i/outpt)))  &
        ! write (*,'(i6,3(f18.12))') i,zp(i)-zp(INT(discret/2)),rhop(i)
     ENDDO
     write (68,*) ' '

     !--------------------------------------------------------------------------
     ! calculate surface tension
     !--------------------------------------------------------------------------
     free =0.0
     DO i=1,discret
        dense(1)  = rhop(i)*z3t
        densta(1)  = rhop(i)*z3t
        ! call CPHS(ro1,fmh,zmh,hsmy)
        CALL FUGACITY (lnphi)
        phs1 = dense(1)/z3t*z_ges               ! phs must be p/kT
        om1  = lnphi(1,1) - my0 - mf_att(i)/2.0
        omega = Pbulk - phs1 + dense(1)/z3t*om1 ! Pbulk must be p/kT
        free  = free  + omega*dzp
     END DO
     surftens = 1.380662* t *free /parame(1,2)**2
     write (*,*) ' '
     write (*,*) 'Temp. [K], Pressure [bar]    ',T,val_conv(4)/1.E5
     write (*,*) 'Density [kg/m**3]            ',density(1),density(2)
     write (*,*) 'Dimensionless Density (rho*) ',rhob(1),rhob(2)
     write (*,*) 'Excess Grand Potential       ',free
     write (*,*) 'Interfacial Tension [mN/m] = ',surftens
     write (*,*) ' '
     write (69,'(7(f18.10))') T,val_conv(4)/1.E5, &
          density(1),density(2),surftens,free,dev

     !--------------------------------------------------------------------------
     ! calculate phase diagram and diagram of surface tension
     !--------------------------------------------------------------------------
     IF (diagram) THEN
        tc=t
        ensemble_flag = 'tp'
        subtract1     = 'no'
        tsav=t
        IF (eos < 4 .AND. val_conv(4) > 1.E6) CALL CRITICAL (tc,pc,rhoc)
        ! IF (eos >= 4 .AND. eos /= 7) CALL LJ_CRITICAL (tc,pc,rhoc)
        t=tsav
        IF (val_conv(4) > 1.E6) THEN
           write (*,'(a,3(f16.4))') 'critical point',tc,pc/1.E5,rhoc
           write (*,*) ' '
        ENDIF
        IF (val_conv(4) < 1.E6) tc=1000.0
        IF ((t+2.5) <= tc) THEN
           t=t+2.0
           IF ((t+8.0) <= tc) t=t+3.0
           IF ((t+25.0) <= tc) t=t+15.0
           val_init(3) = t       ! value of temperature
           end_x  = t    ! end temperature - here T=const.
           steps = 1.0  ! number of steps of running var. - here 1, since T=const.
           nphas = 2
           z3t = PI/6.0* parame(1,1) * ( 1.0  &  ! divided by parame(i,2)**3
                *(1.0-0.12*EXP(-3.0*parame(1,3)/t)) )**3
           CYCLE diagram_condition
        ENDIF
        IF(tc /= 1.E3)write(69,'(6(f18.10))')tc,pc/1.E5,rhoc,rhoc,0.,0.
     ENDIF
     EXIT diagram_condition
  END DO diagram_condition

END SUBROUTINE DFT


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE TRIDIAG (a,b,c,r,u,n)

  USE utilities
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: n
  REAL                                   :: a(n),b(n),c(n),r(n),u(n)

  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER                     :: NMAX = 5000
  INTEGER                                :: j
  REAL                                   :: bet, gam(NMAX)
  !-----------------------------------------------------------------------------

  bet = b(1)
  u(1) = r(1)/bet
  DO j = 2,n
     gam(j) = c(j-1)/bet
     bet = b(j)-a(j)*gam(j)
     IF (bet == 0.0) call paus ('TRIDIAG failed')
     u(j) = (r(j)-a(j)*u(j-1))/bet
  END DO
  DO j = n-1,1,-1
     u(j) = u(j) - gam(j+1)*u(j+1)
  END DO

END SUBROUTINE TRIDIAG


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE SPLINE_PARA
!
! Considered are cubic splines of the form
! P(i)  =   alpha(i) 
!         + beta(i) *(x-x0(i))
!         + gamma(i)*(x-x0(i))**2
!         + delta(i)*(x-x0(i))**3    for i= (0 , stepno-1 )
!
! parameters of the tridiagonal matrix ( M * utri = rtri )
! with matrix-elements: m(j+1,j) = ctri(j)  fuer alle j
!                       m(j,j)   = btri(j)  fuer alle j
!                       m(j-1,j) = atri(j)  fuer j >= 2
! where u(i) = gamma(i) (defined for i=(1,stepno-2))  and gamma(0)=0
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE SPLINE_PARA (dzp,intgrid,utri,stepno)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: stepno,j
  REAL                                   :: dzp
  REAL, DIMENSION(5000)                  :: atri, btri, ctri, rtri, utri
  REAL, DIMENSION(0:5000)                :: intgrid
  !-----------------------------------------------------------------------------


  ! intgrid(index) is the function-value of the spline-grid at the
  ! grid-point 'index'
  ! dzp is the distance of the grid-points
  DO j = 1, (stepno-1)
     atri(j) = dzp
     btri(j) = 4.0*dzp
     ctri(j) = dzp
     rtri(j) = 3.0/dzp*(intgrid(j+1)-2.0*intgrid(j)+intgrid(j-1))
  ENDDO

  CALL TRIDIAG (atri,btri,ctri,rtri,utri,(stepno-1))

END SUBROUTINE SPLINE_PARA


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE SPLINE_COEFF(beta,gamma,delta,dzp,intgrid,utri,stepno)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: stepno, j
  REAL                                   :: dzp, utri(5000)
  REAL, DIMENSION(0:5000)                :: intgrid, beta, gamma, delta
  !-----------------------------------------------------------------------------

  beta(0)  = (intgrid(1)-intgrid(0))/dzp - dzp/3.0*utri(1)
  gamma(0) = 0.0
  delta(0) = (utri(1)-0.0)/(3.0*dzp)
  DO j = 1, (stepno-2)
     gamma(j)= utri(j)
     beta(j) = (intgrid(j+1)-intgrid(j))/dzp - dzp/3.0*(utri(j+1)+2.0*utri(j))
     delta(j)= (utri(j+1)-utri(j))/(3.0*dzp)
  END DO
  gamma(stepno-1)= utri(stepno-1)
  beta(stepno-1) = (intgrid(stepno)-intgrid(stepno-1))/dzp - dzp/3.0*(0.0+2.0*utri(stepno-1))
  delta(stepno-1)= (0.0-utri(stepno-1))/(3.0*dzp)

  ! Integration over splines
  ! I1=0.0
  ! DO j= 0,(stepno-1)
  !   I1 = I1 + intgrid(j)*dzp +0.5*beta(j)*dzp**2 +1.0/3.0*gamma(j)*dzp**3  +0.25*delta(j)*dzp**4
  ! ENDDO

END SUBROUTINE SPLINE_COEFF


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE SPLINE_INT(INTEGRAL,dzp,intgrid,utri,stepno)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER stepno,j
  REAL                                   :: dzp, INTEGRAL, utri(5000)
  REAL, DIMENSION(0:5000)                :: intgrid, beta, gamma, delta
  !-----------------------------------------------------------------------------

  beta(0)  = (intgrid(1)-intgrid(0))/dzp - dzp/3.0*utri(1)
  gamma(0) = 0.0
  delta(0) = (utri(1)-0.0)/(3.0*dzp)
  DO j = 1, (stepno-2)
     gamma(j)= utri(j)
     beta(j) = (intgrid(j+1)-intgrid(j))/dzp - dzp/3.0*(utri(j+1)+2.0*utri(j))
     delta(j)= (utri(j+1)-utri(j))/(3.0*dzp)
  END DO
  gamma(stepno-1)= utri(stepno-1)
  beta(stepno-1) = (intgrid(stepno)-intgrid(stepno-1))/dzp  &
       - dzp/3.0*(0.0+2.0*utri(stepno-1))
  delta(stepno-1)= (0.0-utri(stepno-1))/(3.0*dzp)

  ! Integration over splines
  INTEGRAL = 0.0
  DO j= 0,(stepno-1)
     INTEGRAL = INTEGRAL + intgrid(j)*dzp +0.5*beta(j)*dzp**2  &
          +1.0/3.0*gamma(j)*dzp**3  +0.25*delta(j)*dzp**4
  ENDDO

END SUBROUTINE SPLINE_INT


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! FUNCTION DENS_INV
!
! INVERSION OF HARD SPHERE CHEMICAL POTENTIAL
! routine: R. Evans, Bristol
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

REAL FUNCTION DENS_INV(CPCP)

  USE BASIC_VARIABLES
  USE EOS_VARIABLES, ONLY: z3t
  IMPLICIT None

  !-----------------------------------------------------------------------------
  INTEGER                                :: L, N, I, II, J, K, KK
  INTEGER                                :: NN(3)
  REAL                                   :: DL(3), DM(3), DMM(3)
  REAL                                   :: B(10), C(10), D(10)
  REAL                                   :: A(10,3)
  REAL                                   :: CPCP, CP, DF, DI, S, hsmy, lnphi(np,nc)

  COMMON/HELP/ A,B,C,D,DMM
  DATA DM / 0.0, 0.0, 10.0 /
  DATA DL / 0.1, 0.6, 0.95 /
  DATA NN / 8, 10, 8 /
  !-----------------------------------------------------------------------------

  DL(1) = 0.1  /parame(1,1)
  DL(2) = 0.6  /parame(1,1)
  DL(3) = 0.95 /parame(1,1)

  ! EXTERNAL CPHS
  CP = CPCP
  IF (CPCP >= 1.E6) THEN
     t=1.E-6 /parame(1,1)
     DO L=1,3
        N=NN(L)
        DI=DF
        DF=DL(L)
        !-----------
        dense(1)  = DI*z3t
        densta(1)  = DI*z3t
        CALL FUGACITY (lnphi)
        hsmy=lnphi(1,1)
        !-----------
        ! call cphs(DI,fhs,zhs,hsmy)
        DMM(L)=hsmy
        DO I=1,N
           D(I)=DI+(DF-DI)/(REAL(N)-1.0)*(REAL(I)-1.0)
           !-----------
           dense(1)  = D(I)*z3t
           densta(1)  = D(I)*z3t
           CALL FUGACITY (lnphi)
           hsmy=lnphi(1,1)
           !-----------
           ! call cphs(d(i),fhs,zhs,hsmy)
           C(I)=hsmy-DM(L)
           IF (L == 1) C(I)=EXP(C(I))
           A(I,L)=0.0
        END DO
        DO I=1,N
           S=1.0
           B(2:N)=0.0
           B(1)=1.0
           DO J=1,N
              IF (J == I) EXIT
              S=S/(C(I)-C(J))
              DO KK=2,N
                 K=N-KK+2
                 B(K)=B(K-1)-C(J)*B(K)
              END DO
              B(1)=-C(J)*B(1)
           END DO
           DO K=1,N
              A(K,L)=A(K,L)+D(I)*S*B(K)
           END DO
        END DO
     END DO
     DENS_INV=0.0
  ELSE
     DENS_INV=0.0
     IF (CP >= DMM(1)) THEN
        L=1
        IF (CP > DMM(2)) L=2
        IF (CP > DMM(3)) L=3
        IF (L == 1) CP=EXP(CP)
        CP=CP-DM(L)
        N=NN(L)
        DENS_INV=A(N,L)
        DO I=2,N
           II=N-I+1
           DENS_INV=DENS_INV*CP+A(II,L)
        END DO
     END IF
  END IF

END FUNCTION DENS_INV

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE Barker_org(T_st,d_hs)
!
! calculation of BARKER-HENDERSON effecive diameter dbh
! using the direct numerical integration
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE Barker_org(T_st,d_hs)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: kint, j
  REAL                                   :: gint, summ, xsig, eux, add, xinv
  REAL                                   :: x6, u0, u0r, d_hs, epst, T_st, rm
  !-----------------------------------------------------------------------------

  epst = 1.0 / T_st
  rm = 1.0
  kint = 200
  gint = REAL(kint)
  summ = 1.0 / (2.0 * gint)
  DO j = 1, kint
     xsig = REAL(j) /gint * rm
     IF ( xsig < 0.5 ) then
        eux = 0.0
     ELSE
        xinv = 1.0 / xsig
        x6  = xinv**6
        u0  = 4.0 * x6 * (x6-1.0)
        u0r = u0 * epst
        IF ( u0r > 20.0 ) THEN
           eux =0.0
        ELSE
           eux =EXP( -u0r )
        END IF
     END IF
     add = (1.0-eux) / gint
     summ = summ + add
  END DO
  summ = summ - add/2.0
  d_hs  = rm*summ
  ! dbh  = gam*sig

END SUBROUTINE Barker_org


