!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE INPUT
!
! This subroutine draws user input data from the file
!     INPUT.INP      which is located:    ./INPUT_FILE/INPUT.INP
! This subroutine also provides all pure component parameters as
! well as kij's. (At the end of this routine, the SUBROUTINE
! PARA_INPUT is called, which delivers these pure component
! parameters.)
!
! SUMMERY OF OUTPUT VARIABLES:
! eos           type of equation of state ( 1 = PC-SAFT )
! T, P, ncomp   temperature, pressure, number of components
! xiF           array of feed molefractions. For binary mixtures
!               most calculation options don't require a definition
!               of feed molefractions. Rather, values are generated
!               automatically if all feed molefractions are set to
!               zero in the INPUT-file.
! parame        array of pure component parameters. See
!               SUBROUTINE PARA_INPUT (in file: para_input.f)
! compna, mm    name (string) and molec.mass of pure component
! kij           binary interaction parameter
! uinT,uinP     units of T and p for input and all outputs
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE READ_INPUT

  USE BASIC_VARIABLES
  use utilities, only: file_open
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  REAL                                   :: reading2,reading3,sumfeed
  CHARACTER (LEN=4)                      :: uoutp, uinp
  CHARACTER (LEN=1)                      :: uoutt, uint
  CHARACTER (LEN=50)                     :: filename
  CHARACTER (LEN=30)                     :: reading1
  !-----------------------------------------------------------------------------

  filename = './in.txt'
  CALL file_open(filename,30)
!  READ (30,*) eos, pol 
  eos = 1
  pol = 1

  !READ (30,*) t, uint, p, uinp
  READ(30,*) t
  uint = 'K'
  p = 1.
  uinp = 'bar'



!   ncomp = 0
!   i = 0
!   sumfeed = 0.0
!   read_loop: DO
!      READ (30,*) reading1,reading2,reading3
!      IF (reading1 == 'end') EXIT read_loop
!      ncomp = ncomp + 1
!      i = i + 1
!      compna(i)= reading1    ! comp.name
!      mm(i)    = reading2    ! molec.mass (mandatory only for polymers)
!      xif(i)   = reading3
!      sumfeed  = sumfeed + xif(i)
!   ENDDO read_loop

  READ(30,*) ncomp      !read number of components present in system
  READ(30,*) compna(1)  !read name of component which the solubility is calculated for
  sumfeed = 0.
  Do i = 1,ncomp        !read composition of liquid phase
     mm(i) = 0.     
     READ(30,*) xif(i)
     sumfeed = sumfeed + xif(i)
  End Do


  !Determine wheter two component system (-> solubility in polyurethane polymer)
  !Or three component system (-> solubility in reaction mixture consisting of polyol and MDI)
  If(ncomp == 2) Then
      compna(2) = 'pu'
  Else If (ncomp == 3) Then
      compna(2) = 'polyol'!'hexane'
      compna(3) = 'mdi'!'butane'
  Else
      write(*,*) 'Solubility Code: This code is only for systems with 2 or 3 components!'
      Stop 5 
  End If


  CLOSE (30)

  IF ( sumfeed /= 0.0 .AND. sumfeed /= 1.0 ) THEN
     xif(1:ncomp) = xif(1:ncomp)/sumfeed
  END IF

  uoutt = uint
  uoutp = uinp
  IF (uint == 'C' ) THEN
     u_in_t = 273.15
  ELSE
     u_in_t = 0.0
  END IF
  IF ( uinp == 'bar' ) THEN
     u_in_p = 1.E5
  ELSE IF ( uinp == 'mbar' ) THEN
     u_in_p = 1.E2
  ELSE IF ( uinp == 'MPa' ) THEN
     u_in_p = 1.E6
  ELSE IF ( uinp == 'kPa' ) THEN
     u_in_p = 1.E3
  ELSE
     u_in_p = 1.E0
  END IF

  IF ( uoutt == 'C' ) THEN
     u_out_t = 273.15
  ELSE
     u_out_t = 0.0
  END IF
  IF ( uoutp == 'bar' ) THEN
     u_out_p = 1.E5
  ELSE IF ( uoutp == 'mbar' ) THEN
     u_out_p = 1.E2
  ELSE IF ( uoutp == 'MPa' ) THEN
     u_out_p = 1.E6
  ELSE IF ( uoutp == 'kPa' ) THEN
     u_out_p = 1.E3
  ELSE
     u_out_p = 1.0
  END IF

  t = t + u_in_t
  p = p * u_in_p

  CALL para_input            ! retriev pure comp. parameters

  IF (ncomp == 1) THEN
     WRITE (40,*) '   T       P      rho_1      rho_2      h_LV'
  ELSE IF (ncomp == 2) THEN
     ! WRITE (40,*) ' x2_phase1 x2_phase2  w1_phase1 w2_phase2  T  P  rho1  rho2'
     WRITE (40,*) ' '
  ELSE IF (ncomp == 3) THEN
     WRITE (40,*) ' x1_ph1 x2_ph1 x3_ph1  x1_ph2  x2_ph2  x3_ph2  T  P  rho1  rho2'
  END IF

END SUBROUTINE READ_INPUT



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE OUTPUT
!
! The purpose of this subroutine is obvious.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

SUBROUTINE OUTPUT

  USE BASIC_VARIABLES
  use utilities, only: SI_DENS
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  INTEGER                                :: i
  CHARACTER (LEN=4)                      :: t_ind
  CHARACTER (LEN=4)                      :: p_ind
  CHARACTER (LEN=4)                      :: char_ncomp
  REAL                                   :: density(np),w(np,nc)
  !-----------------------------------------------------------------------------

  CALL SI_DENS ( density, w )

  IF ( u_in_p == 1.E5 ) THEN
     p_ind = ' bar'
  ELSE IF ( u_in_p == 1.E2 ) THEN
     p_ind = 'mbar'
  ELSE IF ( u_in_p == 1.E6 ) THEN
     p_ind = ' MPa'
  ELSE IF ( u_in_p == 1.E3 ) THEN
     p_ind = ' kPa'
  ELSE
     p_ind = ' Pa'
  END IF
  IF ( u_in_t == 273.15 ) THEN
     t_ind = ' C'
  ELSE
     t_ind = ' K'
  END IF

  WRITE(*,*) '--------------------------------------------------'
  WRITE (char_ncomp,'(I3)') ncomp
  WRITE(*,'(t2,a,f7.2,2a,f9.4,a)') ' T =',t-u_out_t,t_ind  &
       ,'   P =',p/u_out_p,p_ind
  WRITE(*,'(t15,4(a12,1x),10x,a)') (compna(i),i=1,ncomp)
  WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  w', (w(1,i),i=1,ncomp)
  WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II w', (w(2,i),i=1,ncomp)
  WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  x', (EXP(lnx(1,i)),i=1,ncomp)
  WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II x', (EXP(lnx(2,i)),i=1,ncomp)
  WRITE(*,'(2x,a,2(g13.6,1x))')                'DENSITY   ', density(1),density(2)

  !-----------------------------------------------------------------------------
  ! output to files
  !-----------------------------------------------------------------------------
  IF ( ncomp == 1 ) THEN
     WRITE (40,'(9(2x,f18.8))') t-u_out_t, p/u_out_p,  &
          density(1), density(2),h_lv,cpres(1),cpres(2),  &
          speed_of_sound(1),speed_of_sound(2)
     !     &      ,(enthal(2)-enthal(1))/mm(1)
     !        WRITE (40,'(4(2x,f15.8))') t, p, 0.3107*dense(1)
     !     &  /0.1617*(0.689+0.311*(T/1.328)**0.3674),0.3107
     !     &  *dense(2)/0.1617*(0.689+0.311*(T/1.328)**0.3674)
  ELSE IF ( ncomp == 2 ) THEN
     WRITE (40,'(12(2x,G15.8))') 1.0-xi(1,1),1.0-xi(2,1),  &
          w(1,1),w(2,1),t-u_out_t, p/u_out_p, density(1),density(2)  &
          ,enthal(1),enthal(2),cpres(1),cpres(2)
  ELSE IF ( ncomp == 3 ) THEN
     WRITE (40,'(10(2x,f15.8))') xi(1,1),xi(1,2),xi(1,3),xi(2,1),xi(2,2),  &
          xi(2,3),t-u_out_t, p/u_out_p, density(1),density(2)
  END IF

END SUBROUTINE OUTPUT

