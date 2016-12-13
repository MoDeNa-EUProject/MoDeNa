
!> This file contains auxiliary subroutines.



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE eos_const
!
! This subroutine provides the constants of the PC-SAFT EOS.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE eos_const (ap,bp,dnm)
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: ap(0:6,3)
 REAL, INTENT(OUT)                      :: bp(0:6,3)
 REAL, INTENT(OUT)                      :: dnm(4,9)
! ----------------------------------------------------------------------


! --- dispersion term constants ----------------------------------------
ap(0,1) =  0.91056314451539
ap(0,2) = -0.30840169182720
ap(0,3) = -0.09061483509767
ap(1,1) =  0.63612814494991
ap(1,2) =  0.18605311591713
ap(1,3) =  0.45278428063920
ap(2,1) =  2.68613478913903
ap(2,2) = -2.50300472586548
ap(2,3) =  0.59627007280101
ap(3,1) = -26.5473624914884
ap(3,2) =  21.4197936296668
ap(3,3) = -1.72418291311787
ap(4,1) =  97.7592087835073
ap(4,2) = -65.2558853303492
ap(4,3) = -4.13021125311661
ap(5,1) = -159.591540865600
ap(5,2) =  83.3186804808856
ap(5,3) =  13.7766318697211
ap(6,1) =  91.2977740839123
ap(6,2) = -33.7469229297323
ap(6,3) = -8.67284703679646

bp(0,1) =  0.72409469413165
bp(0,2) = -0.57554980753450
bp(0,3) =  0.09768831158356
bp(1,1) =  1.11913959304690  *2.0
bp(1,2) =  0.34975477607218  *2.0
bp(1,3) = -0.12787874908050  *2.0
bp(2,1) = -1.33419498282114  *3.0
bp(2,2) =  1.29752244631769  *3.0
bp(2,3) = -3.05195205099107  *3.0
bp(3,1) = -5.25089420371162  *4.0
bp(3,2) = -4.30386791194303  *4.0
bp(3,3) =  5.16051899359931  *4.0
bp(4,1) =  5.37112827253230  *5.0
bp(4,2) =  38.5344528930499  *5.0
bp(4,3) = -7.76088601041257  *5.0
bp(5,1) =  34.4252230677698  *6.0
bp(5,2) = -26.9710769414608  *6.0
bp(5,3) =  15.6044623461691  *6.0
bp(6,1) = -50.8003365888685  *7.0
bp(6,2) = -23.6010990650801  *7.0
bp(6,3) = -4.23812936930675  *7.0


! square-well fluid
!      ap(1,1)=  0.79152347258784
!      ap(1,2)= -0.62269805320654
!      ap(1,3)= -0.06798823934067
!      ap(2,1)=  1.07120982251709
!      ap(2,2)=  0.48628215731716
!      ap(2,3)=  0.02837828512515
!      ap(3,1)=  0.92084839459226
!      ap(3,2)=  1.11652038059747
!      ap(3,3)=  0.09713202077943
!      ap(4,1)= -7.84708350369249
!      ap(4,2)= -2.04200599876547
!      ap(4,3)=  0.06475764015088
!      ap(5,1)= 25.90284137818050
!      ap(5,2)=  9.27791640100603
!      ap(5,3)=  0.07729792971827
!      ap(6,1)= -57.1528726997640
!      ap(6,2)= -16.8377999920957
!      ap(6,3)=  0.24883598436184
!      ap(7,1)= 42.02314637860930
!      ap(7,2)=  7.62432635016420
!      ap(7,3)= -0.72472024688888

!      bp(1,1)=  0.79152347258784
!      bp(1,2)= -0.62269805320654
!      bp(1,3)= -0.06798823934067
!      bp(2,1)=  1.07120982251709  *2.0
!      bp(2,2)=  0.48628215731716  *2.0
!      bp(2,3)=  0.02837828512515  *2.0
!      bp(3,1)=  0.92084839459226  *3.0
!      bp(3,2)=  1.11652038059747  *3.0
!      bp(3,3)=  0.09713202077943  *3.0
!      bp(4,1)= -7.84708350369249  *4.0
!      bp(4,2)= -2.04200599876547  *4.0
!      bp(4,3)=  0.06475764015088  *4.0
!      bp(5,1)= 25.90284137818050  *5.0
!      bp(5,2)=  9.27791640100603  *5.0
!      bp(5,3)=  0.07729792971827  *5.0
!      bp(6,1)= -57.1528726997640  *6.0
!      bp(6,2)= -16.8377999920957  *6.0
!      bp(6,3)=  0.24883598436184  *6.0
!      bp(7,1)= 42.02314637860930  *7.0
!      bp(7,2)=  7.62432635016420  *7.0
!      bp(7,3)= -0.72472024688888  *7.0


dnm(1,1) = -8.8043
dnm(1,2) = +4.1646270
dnm(1,3) = -48.203555
dnm(1,4) = +140.43620
dnm(1,5) = -195.23339
dnm(1,6) = +113.51500
dnm(2,1) = +2.9396
dnm(2,2) = -6.0865383
dnm(2,3) = +40.137956
dnm(2,4) = -76.230797
dnm(2,5) = -133.70055
dnm(2,6) = +860.25349
dnm(2,7) = -1535.3224
dnm(2,8) = +1221.4261
dnm(2,9) = -409.10539
dnm(3,1) = -2.8225
dnm(3,2) = +4.7600148
dnm(3,3) = +11.257177
dnm(3,4) = -66.382743
dnm(3,5) = +69.248785
dnm(4,1) = +0.3400
dnm(4,2) = -3.1875014
dnm(4,3) = +12.231796
dnm(4,4) = -12.110681

END SUBROUTINE eos_const



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE dq_const
!
! This subr. provides the constants of the dipole-quadrupole term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE dq_const ( dqp2,dqp3,dqp4 )
!
 USE PARAMETERS, ONLY: nc
 USE EOS_VARIABLES, ONLY: ncomp, parame
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: dqp2(nc,nc,0:8)
 REAL, INTENT(OUT)                      :: dqp3(nc,nc,nc,0:8)
 REAL, INTENT(OUT)                      :: dqp4(nc,nc,0:8)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k
 REAL                                   :: mdq(nc)
 REAL                                   :: mf1, mf2, msegij
! ----------------------------------------------------------------------

DO i=1,ncomp
  mdq(i) = parame(i,1)
  IF (mdq(i) > 2.0) mdq(i) = 2.0
END DO


DO i=1,ncomp
  DO j=1,ncomp
    
    msegij=(mdq(i)*mdq(j))**0.5
    mf1 = (msegij-1.0)/msegij
    mf2 = mf1*(msegij-2.0)/msegij
    
    dqp2(i,j,0) =  0.697094963 + mf1*(-0.673459279) + mf2*0.670340770
    dqp2(i,j,1) = -0.633554144 + mf1*(-1.425899106) + mf2*(-4.338471826)
    dqp2(i,j,2) =  2.945509028 + mf1 * 4.19441392   + mf2*7.234168360
    dqp2(i,j,3) = -1.467027314 + mf1 * 1.0266216
    dqp2(i,j,4) = 0.0
    
    dqp4(i,j,0) = -0.484038322 + mf1 * 0.67651011   + mf2*(-1.167560146)
    dqp4(i,j,1) =  1.970405465 + mf1*(-3.013867512) + mf2*2.13488432
    dqp4(i,j,2) = -2.118572671 + mf1 * 0.46742656
    dqp4(i,j,3) = 0.0
    dqp4(i,j,4) = 0.0
    
    
    DO k=1,ncomp
      msegij=(mdq(i)*mdq(j)*mdq(k))**(1.0/3.0)
      mf1 = (msegij-1.0)/msegij
      mf2 = (msegij-2.0)/msegij
      dqp3(i,j,k,0) = 0.795009692 + mf1*(-2.099579397)
      dqp3(i,j,k,1) = 3.386863396 + mf1*(-5.941376392)
      dqp3(i,j,k,2) = 0.475106328 + mf1*(-0.178820384)
      dqp3(i,j,k,3) = 0.0
      dqp3(i,j,k,4) = 0.0
    END DO
    
  END DO
END DO

END SUBROUTINE dq_const


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE dd_const
!
! This subroutine provides the constants of the dipole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE dd_const ( ddp2,ddp3,ddp4 )
!
 USE PARAMETERS, ONLY: nc, PI
 USE EOS_VARIABLES, ONLY: ncomp, parame
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: ddp2(nc,nc,0:8)
 REAL, INTENT(OUT)                      :: ddp3(nc,nc,nc,0:8)
 REAL, INTENT(OUT)                      :: ddp4(nc,nc,0:8)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k
 REAL                                   :: pardd(nc)
 REAL                                   :: mf1,mf2,msegij,sin2t
! ----------------------------------------------------------------------

sin2t = SIN( 0.0 * PI / 180.0 )
sin2t = sin2t*sin2t

DO i = 1, ncomp
  pardd(i) = parame(i,1)
  IF (pardd(i) > 2.0) pardd(i) = 2.0
END DO

DO i=1,ncomp
  DO j=1,ncomp
!      IF (parame(i,6).NE.0.0.AND.parame(j,6).NE.0.0) THEN
    
    msegij=(pardd(i)*pardd(j))**0.5
    mf1 = (msegij-1.0)/msegij
    mf2 = mf1*(msegij-2.0)/msegij
    
    ddp2(i,j,0) =  0.30435038064 + mf1*(0.95346405973+0.201436*sin2t)  &
                  + mf2*(-1.16100802773-1.74114*sin2t)
    ddp2(i,j,1) = -0.13585877707 + mf1*(-1.83963831920+1.31649*sin2t)  &
                  + mf2*4.52586067320
    ddp2(i,j,2) =  1.44933285154 + mf1 * 2.01311801180  + mf2*0.97512223853
    ddp2(i,j,3) =  0.35569769252 + mf1*(-7.37249576667) + mf2*(-12.2810377713)
    ddp2(i,j,4) = -2.06533084541 + mf1 * 8.23741345333  + mf2*5.93975747420
    
    ddp4(i,j,0) =  0.21879385627 + mf1*(-0.58731641193) + mf2*3.48695755800
    ddp4(i,j,1) = -1.18964307357 + mf1 * 1.24891317047  + mf2*(-14.9159739347)
    ddp4(i,j,2) =  1.16268885692 + mf1*(-0.50852797392) + mf2*15.3720218600
    ddp4(i,j,3) =  0.0
    ddp4(i,j,4) =  0.0
    
    DO k=1,ncomp
!      IF (parame(k,6).NE.0.0) THEN
      msegij=(pardd(i)*pardd(j)*pardd(k))**(1.0/3.0)
      mf1 = (msegij-1.0)/msegij
      mf2 = mf1*(msegij-2.0)/msegij
      ddp3(i,j,k,0) = -0.06467735252 + mf1*(-0.95208758351+0.28503*sin2t)  &
                      + mf2*(-0.62609792333+2.2195*sin2t)
      ddp3(i,j,k,1) =  0.19758818347 + mf1 * 2.99242575222  + mf2*1.29246858189
      ddp3(i,j,k,2) = -0.80875619458 + mf1*(-2.38026356489) + mf2*1.65427830900
      ddp3(i,j,k,3) =  0.69028490492 + mf1*(-0.27012609786) + mf2*(-3.43967436378)
      ddp3(i,j,k,4) =  0.0
      
!      ENDIF
    END DO
    
!      ENDIF
  END DO
END DO

END SUBROUTINE dd_const


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE qq_const
!
! This subroutine provides the constants of the quadrupole-term.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE qq_const ( qqp2,qqp3,qqp4 )
!
 USE PARAMETERS, ONLY: nc
 USE EOS_VARIABLES, ONLY: ncomp, parame
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: qqp2(nc,nc,0:8)
 REAL, INTENT(OUT)                      :: qqp3(nc,nc,nc,0:8)
 REAL, INTENT(OUT)                      :: qqp4(nc,nc,0:8)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, k
 REAL                                   :: mqq(nc)
 REAL                                   :: mf1, mf2, msegij
! ----------------------------------------------------------------------

DO i = 1,ncomp
  mqq(i) = parame(i,1)
  IF (mqq(i) > 2.0) mqq(i) = 2.0
END DO

DO i = 1,ncomp
  DO j = 1,ncomp
    IF (parame(i,7) /= 0.0 .AND. parame(j,7) /= 0.0) THEN
      
      msegij=(mqq(i)*mqq(j))**0.5
!      msegij=(parame(i,1)*parame(j,1))**0.50
      mf1 = (msegij-1.0)/msegij
      mf2 = mf1*(msegij-2.0)/msegij
      
      qqp2(i,j,0) =  1.237830788 + mf1 * 1.285410878  + mf2*1.794295401
      qqp2(i,j,1) =  2.435503144 + mf1*(-11.46561451) + mf2*0.769510293
      qqp2(i,j,2) =  1.633090469 + mf1 *22.08689285   + mf2*7.264792255
      qqp2(i,j,3) = -1.611815241 + mf1 * 7.46913832   + mf2*94.48669892
      qqp2(i,j,4) =  6.977118504 + mf1*(-17.19777208) + mf2*(-77.1484579)
      
      qqp4(i,j,0) =  0.454271755 + mf1*(-0.813734006) + mf2*6.868267516
      qqp4(i,j,1) = -4.501626435 + mf1 * 10.06402986  + mf2*(-5.173223765)
      qqp4(i,j,2) =  3.585886783 + mf1*(-10.87663092) + mf2*(-17.2402066)
      qqp4(i,j,3) =  0.0
      qqp4(i,j,4) =  0.0
      
      DO k = 1,ncomp
        IF (parame(k,7) /= 0.0) THEN
          msegij=(mqq(i)*mqq(j)*mqq(k))**(1.0/3.0)
!      msegij=(parame(i,1)*parame(j,1)*parame(k,1))**(1.0/3.0)
          mf1 = (msegij-1.0)/msegij
          mf2 = mf1*(msegij-2.0)/msegij
          qqp3(i,j,k,0) = -0.500043713 + mf1 * 2.000209381 + mf2*3.135827145
          qqp3(i,j,k,1) =  6.531869153 + mf1*(-6.78386584) + mf2*7.247588801
          qqp3(i,j,k,2) = -16.01477983 + mf1 * 20.38324603 + mf2*3.075947834
          qqp3(i,j,k,3) =  14.42597018 + mf1*(-10.89598394)
          qqp3(i,j,k,4) =  0.0
        END IF
      END DO
      
    END IF
  END DO
END DO

END SUBROUTINE qq_const

 





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
 SUBROUTINE SET_DEFAULT_EOS_NUMERICAL
!
 USE EOS_NUMERICAL_DERIVATIVES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------

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









! SUBROUTINE READ_INPUT
! !
!  USE BASIC_VARIABLES
!  IMPLICIT NONE
! !
! ! ----------------------------------------------------------------------
!  INTEGER                                :: i
!  REAL                                   :: reading2,reading3,sumfeed
!  CHARACTER (LEN=4)                      :: uoutp, uinp
!  CHARACTER (LEN=1)                      :: uoutt, uint
!  CHARACTER (LEN=50)                     :: filename
!  CHARACTER (LEN=30)                     :: reading1
! ! ----------------------------------------------------------------------
! 
!  filename='./input_file/INPUT.INP'
!  CALL file_open(filename,30)
!  READ (30,*) eos, pol		!J: specify by numbers! eos(1=pcsaft, 2=SRK,...) pol (=polar) yes(1) no(0)
!  READ (30,*) t, uint, p, uinp	!J: t: value of temp, uint: unit of temp, p: value of pressure, uinp: unit of pressure
! 
!  ncomp = 0
!  i = 0
!  sumfeed = 0.0
!  read_loop: DO
!     READ (30,*) reading1,reading2,reading3
!     IF (reading1 == 'end') EXIT read_loop
!     ncomp = ncomp + 1
!     i = i + 1
!     compna(i)= reading1    ! comp.name
!     mm(i)    = reading2    ! molec.mass (mandatory only for polymers)
!     xif(i)   = reading3	   !J: molefractions
!     sumfeed = sumfeed + xif(i)
!  ENDDO read_loop
! 
!  CLOSE (30)
! 
!  IF (sumfeed /= 0.0 .AND. sumfeed /= 1.0) THEN	!J: in case mole fractions dont sum up to 1??
!     xif(1:ncomp) = xif(1:ncomp)/sumfeed
!  END IF
! 
!  uoutt = uint
!  uoutp = uinp
!  IF (uint == 'C') THEN			!J: unit stuff
!     u_in_t = 273.15
!  ELSE
!     u_in_t = 0.0
!  END IF
!  IF (uinp == 'bar') THEN
!     u_in_p = 1.E5
!  ELSE IF (uinp == 'mbar') THEN
!     u_in_p = 1.E2
!  ELSE IF (uinp == 'MPa') THEN
!     u_in_p = 1.E6
!  ELSE IF (uinp == 'kPa') THEN
!     u_in_p = 1.E3
!  ELSE
!     u_in_p = 1.E0
!  END IF
! 
!  IF (uoutt == 'C') THEN
!     u_out_t = 273.15
!  ELSE
!     u_out_t = 0.0
!  END IF
!  IF (uoutp == 'bar') THEN
!     u_out_p = 1.E5
!  ELSE IF (uoutp == 'mbar') THEN
!     u_out_p = 1.E2
!  ELSE IF (uoutp == 'MPa') THEN
!     u_out_p = 1.E6
!  ELSE IF (uoutp == 'kPa') THEN
!     u_out_p = 1.E3
!  ELSE
!     u_out_p = 1.0
!  END IF
! 
!  t = t + u_in_t			!J: calculate temp in Kelvin	
!  p = p * u_in_p			!J: calculate pressure in Pascal
! 
!  CALL para_input            ! retriev pure comp. parameters
! 
!  
!  END SUBROUTINE READ_INPUT





!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE file_open
!
! This subroutine opens files for reading. Beforehand, it checks
! whether this file is available.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE file_open(filename,file_number)
!
! ----------------------------------------------------------------------
 CHARACTER (LEN=50)                       :: filename
 INTEGER                                  :: file_number
 LOGICAL                                  :: filefound
! ----------------------------------------------------------------------

INQUIRE (FILE=filename, EXIST = filefound)
IF (filefound) THEN
  OPEN (file_number, FILE = filename)
ELSE
  write (*,*) ' FOLLOWING FILE CAN NOT BE OPENED', filename
  stop
END IF

END SUBROUTINE file_open



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE para_input
!
! This subroutine provides pure component parameters and kij parameters.
! The following syntax applies:
!
! compna(i)                  component name
! parame(i,k)                pure comp. parameter:
!              parame(i,1):  segment number  [/]
!              parame(i,2):  segment diameter "sigma" [Angstrom]
!              parame(i,3):  segment energy param. epsilon/k [K]
!              parame(i,4):  model parameter; not used for PC-SAFT (=0)
!                            it is 10K most of the time for SAFT [K]
!              parame(i,5):  Param. for T-dependent segment diameter [/]
!              parame(i,6):  dipolar moment [debye]
!              parame(i,7):  quadrupolar moment [debye]
!              parame(i,8):  number of segments that are part of a branching 4-mer [/]
!              parame(i,9):
!              parame(i,10): ionic charge number (positiv or negativ) [/]
!              parame(i,11): polarizability [A**3]
!              parame(i,12): number of association sites [/]
!              parame(i,13): (=kap_hb, see below) [/]
!              parame(i,14 to 25): (=eps_hb, see below) [K]
! nhb_typ(i)                 number of different types of association sites (comp. i)
! nhb_no(i,k)                number of association sites of type k
! eps_hb                     depth of association potential [K]
! kap_hb                     effective width of assoc. potential (angle-averg.)
! mm                         molec. mass
! scaling                    param. for roughly scaling the set of objective functions
!
! As opposed to low-molec mass compounds, the molecular mass of a
! polymer is not obtained from this routine. Rather, it is a
! user-specification given in the file INPUT.INP
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE para_input
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
!----------------------------------------------------------------------
 INTEGER                                :: i
!----------------------------------------------------------------------

IF (eos == 1) THEN
  CALL pcsaft_par
ELSE IF (eos == 4 .OR. eos == 5 .OR. eos == 6 .OR. eos == 8) THEN
  ! CALL lj_par
  write (*,*) 'deactivated this line when making a transition to f90'
  stop
ELSE IF (eos == 7) THEN
  ! CALL sw_par
  write (*,*) 'deactivated this line when making a transition to f90'
  stop
ELSE IF (eos == 10) THEN
  i = 1
  IF (compna(i) == 'LC_generic' .AND. ncomp == 1 ) THEN
    mm(i)       = 1.0
    parame(i,1) = 7.0
    parame(i,2) = 1.0
    parame(i,3) = 0.0
  ELSE
    write (*,*) 'PARA_INPUT: define the component !'
    stop
  ENDIF
ELSE
  !CALL saft_par
END IF

DO i = 1, ncomp
  IF ( mm(i) >= 1.0       .AND. mm(i) < 45.0  ) THEN
    scaling(i) = 10000.0
  ELSE IF( mm(i) >= 45.0  .AND. mm(i) < 90.0  ) THEN
    scaling(i) = 1000.0
  ELSE IF( mm(i) >= 90.0  .AND. mm(i) < 150.0 ) THEN
    scaling(i) = 100.0
  ELSE IF( mm(i) >= 150.0 .AND. mm(i) < 250.0 ) THEN
    scaling(i) = 10.0
  ELSE
    scaling(i) = 1.0
  END IF
  IF (parame(i,10) /= 0.0) scaling(i) = scaling(i) / 1.E4   ! Electrolytes
END DO

END SUBROUTINE para_input




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE pcsaft_par
!
! pure component parameters and kij parameters
! (as described in SUBROUTINE para_input)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE pcsaft_par
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
!----------------------------------------------------------------------
 INTEGER                                :: i, j, k, no
 INTEGER, DIMENSION(nc)                 :: nhb_typ
 INTEGER, DIMENSION(nc,nsite)           :: nhb_no
 REAL, DIMENSION(nc,nc,nsite,nsite)     :: eps_hb
 REAL, DIMENSION(nc,nc)                 :: kap_hb
!----------------------------------------------------------------------


DO  i = 1, ncomp
  parame(i,4) = 0.0     ! T correct. required for SAFT, not PC-SAFT
  parame(i,5) = 0.12    ! Param. for T-dependent segment diameter
  parame(i,6) = 0.0     ! dipolar moment
  parame(i,7) = 0.0     ! quadrupolar moment
  parame(i,8) = 0.0     ! number of segments that are part of a branching 4-mer
  parame(i,9) = 0.0
  parame(i,10)= 0.0     ! ionic charge number
  parame(i,11)= 0.0     ! polarizability
  lli(i)      = 0.0
  phi_criti(i)= 0.0
  chap(i)     = 0.0
  
  nhb_typ(i)  = 0
  kap_hb(i,i) = 0.0
  ! irgendwann sollten nhb_typ und kap_hb durch parame(i,12) und (i,13)
  ! ersetzt werden.
  
   IF (compna(i) == '14-butandiol') THEN
    mm(i)       = 90.12
    parame(i,1) = 4.35923557
    parame(i,2) = 3.02947364
    parame(i,3) = 197.11998863
  
    ELSE IF (compna(i) == 'air') THEN
    mm(i)       = 28.899	!n2 and o2 according to mole fractions
    parame(i,1) = 1.18938	!n2 and o2 according to mole fractions (weighted artihm. avg)
    parame(i,2) = 3.28694	!n2 and o2 according to mole fractions (weighted artihm. avg)
    parame(i,3) = 95.672	!n2 and o2 according to mole fractions (weighted artihm. avg)

   Else IF(compna(i) == 'surfactant') THEN
     mm(i) =  2655.24078
     parame(i,1) = 78.5859962
     parame(i,2) = 4.17006833
     parame(i,3) = 230.284526
     parame(i,6) = 17.9645
    
    
    
  
   Else IF(compna(i) == 'mdi') THEN
     mm(i) =  2.50252E+02
     parame(i,1) = mm(i)*0.030769
     parame(i,2) = 2.886003
     parame(i,3) = 283.052778
  
     Else IF(compna(i) == 'po') THEN
    mm(i)       = 90.12
    parame(i,1) = 4.35923557
    parame(i,2) = 3.02947364
    parame(i,3) = 197.11998863
    Else IF(compna(i) == 'pu') THEN
!      mm(i) =  2042.22 !pu n = 5
!      parame(i,1) = mm(i)*0.008845
!      parame(i,2) = 5.680270
!      parame(i,3) = 497.997594
     mm(i) =  340.37 !pu n = 0
     parame(i,1) = mm(i)*0.043312
     parame(i,2) = 3.008359
     parame(i,3) = 273.445205
!     mm(i) =  680.74 !pu n = 1
!     parame(i,1) = mm(i)*0.024106
!     parame(i,2) = 3.744327
!     parame(i,3) = 321.486386
!       mm(i) =  1021.11 !pu n = 2
!       parame(i,1) = mm(i)*0.015076
!       parame(i,2) = 4.537837
!       parame(i,3) = 400.036950




   Else IF(compna(i) == 'tpg') THEN
     mm(i) =  192.25
     parame(i,1) = mm(i)*0.01239
     parame(i,2) =  4.549
     parame(i,3) = 148.678
     parame(i,6) = 0.41

    nhb_typ(i)  = 2             ! no. of different association sites
    nhb_no(i,1) = 2            ! no. of sites of type 1
    nhb_no(i,2) = 2            ! no. of sites of type 2

    eps_hb(i,i,1,2)= 5597.844
    eps_hb(i,i,2,1)= 5597.844
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 0.03

    
    
  
  
  ELSE IF (compna(i) == 'ps') THEN
    parame(i,1) = mm(i)*1.9E-2
    parame(i,2) = 4.10705961
    parame(i,3) = 267.0
  ELSE IF (compna(i) == 'pg2') THEN !Polyglycerol 2
    mm(i) = 2000.0
    parame(i,1) = mm(i)*2.37E-2 ! from figure 5 PCSAFT paper
    parame(i,2) = 3.8           ! from figure 5 PCSAFT paper
    parame(i,3) = 270.0         ! starting value for iteration
    ! this is the extra parameter
    parame(i,8) = mm(i)*2.37E-2
    
    nhb_typ(i)  = 2             ! no. of different association sites
    nhb_no(i,1) = 27            ! no. of sites of type 1
    nhb_no(i,2) = 27            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2544.6     ! taken from butanol (same M/OH)
    eps_hb(i,i,2,1)= 2544.6
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)= .00489087833   ! taken from butanol (same M/OH)
  ELSE IF (compna(i) == 'peva') THEN
    parame(i,1) = mm(i)*2.63E-2
    ! --   0 Gew.% VA-------------
    ! parame(i,2) = 4.021767
    ! parame(i,3) = 249.5
    ! -- 7.5 Gew.% VA-------------
    ! parame(i,2) = 4.011
    ! parame(i,3) = 248.1864
    ! parame(i,3) = 247.6286
    ! ---12.7 Gew.% VA------------
    ! parame(i,2) = 4.0028
    ! parame(i,3) = 247.2075
    ! parame(i,3) = 246.24454
    ! ---27.3 Gew.% VA------------
    ! parame(i,2) = 3.9762
    ! parame(i,3) = 244.114
    ! parame(i,3) = 241.9345
    ! ---31.8 Gew.% VA------------
    parame(i,2) = 3.9666
    parame(i,3) = 243.0436
    ! parame(i,3) = 240.46
    ! ---42.7 Gew.% VA------------
    ! parame(i,2) = 3.9400
    ! parame(i,3) = 240.184
    ! parame(i,3) = 236.62
    ! ---------------
  ELSE IF (compna(i) == 'pp') THEN
    parame(i,1) = mm(i)*2.2E-2
    parame(i,2) = 4.2
    parame(i,3) = 220.0
    
    parame(i,1) = mm(i)*0.0230487701
    parame(i,2) = 4.1
    parame(i,3) = 217.0
  ELSE IF (compna(i) == 'pe') THEN
    parame(i,1) = mm(i)*2.622E-2
    parame(i,2) = 4.021767
    parame(i,3) = 252.0
    ! HDPE: extrapolated from pure comp. param. of n-alkane series!
    ! parame(i,1) = mm(i)*2.4346E-2
    ! parame(i,2) = 4.07182
    ! parame(i,3) = 269.67
    !!  parame(i,3) = 252.5
  ELSE IF (compna(i) == 'ldpe') THEN
    parame(i,1) = mm(i)*2.63E-2
    parame(i,2) = 4.021767
    parame(i,3) = 249.5
  ELSE IF (compna(i) == 'pba') THEN
    parame(i,1) = mm(i)*2.5872E-2
    parame(i,2) = 3.95
    parame(i,3) = 229.0
  ELSE IF (compna(i) == 'dextran') THEN
    parame(i,1) = mm(i)*2.E-2
    parame(i,2) = 4.0
    parame(i,3) = 300.0
  ELSE IF (compna(i) == 'glycol-ethers') THEN
    ! mm(i) = 218.0
    ! parame(i,1) = 7.4044
    ! parame(i,2) = 3.61576
    ! parame(i,3) = 244.0034598
    mm(i) = 222.0
    parame(i,1) = 7.994
    parame(i,2) = 3.445377778
    parame(i,3) = 234.916506
  ELSE IF (compna(i) == 'LJ') THEN
    mm(i)       = 1.0
    parame(i,1)  = 1.0
    parame(i,2) = 1.0
    parame(i,3) = 1.0
  ELSE IF (compna(i) == 'LJ1205') THEN
    mm(i)       = 1.0
    parame(i,1)  = 1.0
    parame(i,2) = 1.0
    parame(i,3) = 140.0
  ELSE IF (compna(i) == 'adamantane') THEN
    mm(i)       =   136.235000000000     
    parame(i,1) =   4.81897145432221     
    parame(i,2) =   3.47128575274660     
    parame(i,3) =   266.936967922521     
  ELSE IF (compna(i) == 'methane') THEN
    mm(i)       = 16.043
    parame(i,1)  = 1.0
    parame(i,2) = 3.70388767
    parame(i,3) = 150.033987
    ! LLi(i)      = 1.185*parame(i,2)
    ! phi_criti(i)= 11.141
    ! chap(i)     = 0.787
    lli(i)      = 1.398*parame(i,2)
    phi_criti(i)= 16.01197
    chap(i)     = 0.6
    IF (pol == 2) parame(i,11)= 2.593
    ! --- adjusted to Tc, Pc und omega ---
    ! mm(i)       = 16.0430000000000
    ! parame(i,1) = 1.03353666429362
    ! parame(i,2) = 3.64824920605089
    ! parame(i,3) = 147.903953522994
    lli(i)      = 2.254442763775*parame(i,2)
    phi_criti(i)= 42.060975627454
    chap(i)     = 0.704895924
    lli(i)      = 1.935801125833*parame(i,2)
    phi_criti(i)= 26.363325937261
    chap(i)     = 0.700112854298
    lli(i)      = 2.610103087662*parame(i,2)
    phi_criti(i)= 38.192854403173
    chap(i)     = 0.812100472735
    ! 2.122960316503   34.937141524804    0.734513223627
    ! 2.082897379591   33.036391564859    0.877578492999
  ELSE IF (compna(i) == 'ethane') THEN
    mm(i)       = 30.070
    parame(i,1) =mm(i)*  .0534364758
    parame(i,2) = 3.5205923
    parame(i,3) = 191.423815
    lli(i)      = 1.40*parame(i,2)
    phi_criti(i)= 15.38
    chap(i)     = 0.520
    IF (pol == 2) parame(i,11)= 4.3
    ! --- adjusted to Tc, Pc und omega ---
    ! mm(i)       = 30.069
    ! parame(i,1) = 1.74034548122
    ! parame(i,2) = 3.4697441893134
    ! parame(i,3) = 181.90770083591
    IF (pol >= 1) mm(i)       = 30.0700000000000
    IF (pol >= 1) parame(i,1) = mm(i)*  5.341907666260094E-002
    IF (pol >= 1) parame(i,2) = 3.52104466654628
    IF (pol >= 1) parame(i,3) = 191.449300423694
    IF (pol >= 1) parame(i,7) = 0.650000000000000
    IF (pol >= 1) lli(i)      = 0.0
    IF (pol >= 1) phi_criti(i)= 0.0
    IF (pol >= 1) chap(i)     = 0.0
  ELSE IF (compna(i) == 'propane') THEN
    mm(i)       = 44.096
    parame(i,1) = mm(i)*  .0453970622
    parame(i,2) = 3.61835302
    parame(i,3) = 208.110116
    lli(i)      = 1.8*parame(i,2)
    phi_criti(i)= 21.0
    chap(i)     = 1.0
    lli(i)      = 1.63*parame(i,2)
    phi_criti(i)= 20.37
    chap(i)     = 0.397
    IF (pol == 2) parame(i,11)= 6.29
  ELSE IF (compna(i) == 'butane_debug') THEN
    mm(i)       = 58.123
    parame(i,1) = 2.3374
    parame(i,2) = 3.6655
    parame(i,3) = 214.805
  ELSE IF (compna(i) == 'butane') THEN
    mm(i)       = 58.123
    parame(i,1) = mm(i)*  .0401146927
    parame(i,2) = 3.70860139
    parame(i,3) = 222.877405
    lli(i)      = 1.75*parame(i,2)
    phi_criti(i)= 23.43
    chap(i)     = 0.304
    ! LLi(i)      = 1.942079633622*parame(i,2)
    ! phi_criti(i)= 24.527323443155
    ! chap(i)     = 0.734064026277
    ! LLi(i)      = 1.515115760477*parame(i,2)
    ! phi_criti(i)= 17.682929717796
    ! chap(i)     = 0.335848717079
    IF (pol == 2) parame(i,11)= 8.2
    ! --- adjusted to Tc, Pc und omega ---
    ! mm(i)       = 58.1230000000
    ! parame(i,1) = 2.45352304112
    ! parame(i,2) = 3.74239117802
    ! parame(i,3) = 214.185157925
  ELSE IF (compna(i) == 'pentane') THEN
    mm(i)       = 72.146
    parame(i,1) = mm(i)*  .03727896
    parame(i,2) = 3.77293174
    parame(i,3) = 231.197015
    IF (pol == 2) parame(i,11)= 9.99
  ELSE IF (compna(i) == 'hexane') THEN
    mm(i)       = 86.177
    parame(i,1) = mm(i)*  .0354812325
    parame(i,2) = 3.79829291
    parame(i,3) = 236.769054
    lli(i)      = 2.24*parame(i,2)
    phi_criti(i)= 33.25
    chap(i)     = 0.205
    IF (pol == 2) parame(i,11)= 11.9
  ELSE IF (compna(i) == 'heptane') THEN
    mm(i)       = 100.203
    parame(i,1) = mm(i)*  .034762384
    parame(i,2) = 3.80487025
    parame(i,3) = 238.400913
    lli(i)      = 2.35*parame(i,2)
    phi_criti(i)= 38.10
    chap(i)     = 0.173
    IF (pol == 2) parame(i,11)= 13.61
  ELSE IF (compna(i) == 'octane') THEN
    mm(i)       = 114.231
    parame(i,1) = mm(i)*  .0334228038
    parame(i,2) = 3.83732677
    parame(i,3) = 242.775853
    ! LLi(i)      = 2.0*parame(i,2)
    ! phi_criti(i)= 18.75
    ! chap(i)     = 1.0
    lli(i)      = 2.63*parame(i,2)
    phi_criti(i)= 42.06
    chap(i)     = 0.155
    IF (pol == 2) parame(i,11)= 15.9
  ELSE IF (compna(i) == 'nonane') THEN
    mm(i)       = 128.25
    parame(i,1) = mm(i)*  .0328062594
    parame(i,2) = 3.84483643
    parame(i,3) = 244.508457
  ELSE IF (compna(i) == 'decane') THEN
    mm(i)       = 142.285
    parame(i,1) = mm(i)*  .03277373
    parame(i,2) = 3.8384498
    parame(i,3) = 243.866074
    lli(i)      = 1.845*parame(i,2)
    phi_criti(i)= 21.27
    chap(i)     = 1.0
    lli(i)      = 2.68*parame(i,2)
    phi_criti(i)= 45.0
    chap(i)     = 0.15
    IF (pol == 2) parame(i,11)= 19.1
    ! --- adjusted to Tc, Pc und omega ---
    !  parame(i,1) = 4.794137228322
    !  parame(i,2) = 4.030446690586
    !  parame(i,3) = 236.5884493386
  ELSE IF (compna(i) == 'dodecane') THEN
    mm(i)       = 170.338
    parame(i,1) = mm(i)*  .0311484156
    parame(i,2) = 3.89589236
    parame(i,3) = 249.214532
  ELSE IF (compna(i) == 'hexadecane') THEN
    mm(i)       = 226.446
    parame(i,1) = mm(i)*  .0293593045
    parame(i,2) = 3.95516743
    parame(i,3) = 254.700131
  ELSE IF (compna(i) == 'octadecane') THEN
    mm(i)       = 254.5
    parame(i,1) = 7.3271
    parame(i,2) = 3.9668
    parame(i,3) = 256.20
    IF (pol == 2) parame(i,11)= 30.2
    ! --- adjusted to Tc, Pc und omega ---
    !  mm(i)       = 226.446000000000
    !  parame(i,1) = 6.66976520488694
    !  parame(i,2) = 4.25025597912511
    !  parame(i,3) = 249.582941976119
  ELSE IF (compna(i) == 'eicosane') THEN
    mm(i)       = 282.553
    parame(i,1) = mm(i)*  .0282572812
    parame(i,2) = 3.98692612
    parame(i,3) = 257.747939
  ELSE IF (compna(i) == 'triacontane') THEN
    ! mm(i)       = 422.822             ! polyethylene parameters
    ! parame(i,1) = mm(i)*2.622E-2
    ! parame(i,2) = 4.021767
    ! parame(i,3) = 252.0
    mm(i)       = 422.822             ! param. by extrapolation of n-alkanes
    parame(i,1) = mm(i)*  0.026922527
    parame(i,2) = 4.007608009
    parame(i,3) = 262.28622
  ELSE IF (compna(i) == 'octaeicosane') THEN
    mm(i)       = 395.0             ! param. by extrapolation of n-alkanes (sloppy!!)
    parame(i,1) = mm(i)*  0.026922527
    parame(i,2) = 4.007608009
    parame(i,3) = 262.28622
  ELSE IF (compna(i) == 'tetracontane') THEN
    ! mm(i)       = 563.1             ! polyethylene parameters
    ! parame(i,1) = mm(i)*2.622E-2
    ! parame(i,2) = 4.021767
    ! parame(i,3) = 252.0
    mm(i)       = 563.1             ! param. by extrapolation of n-alkanes
    parame(i,1) = mm(i)*0.026287593
    parame(i,2) = 4.023277
    parame(i,3) = 264.10466
  ELSE IF (compna(i) == 'isobutane') THEN
    mm(i)       = 58.123
    parame(i,1) = mm(i)*  .0389105395
    parame(i,2) = 3.75735249
    parame(i,3) = 216.528584
  ELSE IF (compna(i) == 'isopentane') THEN
    mm(i)       = 72.15
    parame(i,1) = 2.5620
    parame(i,2) = 3.8296
    parame(i,3) = 230.75
  ELSE IF (compna(i) == '2-methylpentane') THEN
    mm(i)       = 86.177
    parame(i,1) = mm(i)*  .0340166994
    parame(i,2) = 3.85354665
    parame(i,3) = 235.5801
  ELSE IF (compna(i) == '23-dimethylbutane') THEN
    mm(i)       = 86.177
    parame(i,1) = mm(i)*  .0311599207
    parame(i,2) = 3.9544545
    parame(i,3) = 246.068188
  ELSE IF (compna(i) == 'ethylene') THEN
    mm(i)       = 28.05
    parame(i,1) = mm(i)*  .0567939013
    parame(i,2) = 3.44499904
    parame(i,3) = 176.468725
    IF (pol == 2) parame(i,11)= 4.252
! eigener 3-ter Anlauf.
    IF (pol >= 1) parame(i,1) = mm(i)*  5.574644443117726E-002
    IF (pol >= 1) parame(i,2) = 3.43281482228714
    IF (pol >= 1) parame(i,3) = 178.627308564610
    IF (pol >= 1) parame(i,7) = 1.56885870200446
    IF (pol == 2) parame(i,11)= 4.252
  ELSE IF (compna(i) == 'propylene') THEN
    mm(i)       = 42.081
    parame(i,1) = mm(i)*  .0465710324
    parame(i,2) = 3.53559831
    parame(i,3) = 207.189309
    ! --- adjusted to Tc, Pc und omega ---
    ! mm(i)       = 42.081
    ! parame(i,1) = 2.086735327675
    ! parame(i,2) = 3.536779407969
    ! parame(i,3) = 198.3529810625
  ELSE IF (compna(i) == '1-butene') THEN
    mm(i)       = 56.107
    parame(i,1) = mm(i)*  .0407524782
    parame(i,2) = 3.64305136
    parame(i,3) = 222.002756
    IF (pol == 2) parame(i,11)= 7.97
  ELSE IF (compna(i) == '1-pentene') THEN
    mm(i)       = 70.134
    parame(i,1) = 2.6006
    parame(i,2) = 3.7399
    parame(i,3) = 231.99
  ELSE IF (compna(i) == '1-hexene') THEN
    mm(i)       = 84.616
    parame(i,1) = mm(i)*  .0352836857
    parame(i,2) = 3.77529612
    parame(i,3) = 236.810973
  ELSE IF (compna(i) == '1-octene') THEN
    mm(i)       = 112.215
    parame(i,1) = mm(i)*  .033345175
    parame(i,2) = 3.81329011
    parame(i,3) = 243.017587
  ELSE IF (compna(i) == 'cyclopentane') THEN
    mm(i)       = 70.13
    parame(i,1) = mm(i)*  .0337262571
    parame(i,2) = 3.71139254
    parame(i,3) = 265.828755
  ELSE IF (compna(i) == 'cyclohexane') THEN
    mm(i)       = 84.147
    parame(i,1) = mm(i)*  .0300695505
    parame(i,2) = 3.84990887
    parame(i,3) = 278.108786
    IF (pol == 2) parame(i,11)= 10.87
  ELSE IF (compna(i) == 'toluene') THEN
    mm(i)       = 92.141
    parame(i,1) = mm(i)*  .0305499338
    parame(i,2) = 3.71689689
    parame(i,3) = 285.68996
    IF (pol == 2) parame(i,11)= 11.8
    ! --- adjusted to Tc, Pc und omega ---
    ! mm(i)       = 92.141
    ! parame(i,1) = 3.002119827762
    ! parame(i,2) = 3.803702734224
    ! parame(i,3) = 271.9428642880
  ELSE IF (compna(i) == 'm-xylene') THEN
    mm(i)       = 106.167
    parame(i,1) = mm(i)*  .030011086
    parame(i,2) = 3.75625585
    parame(i,3) = 283.977525
  ELSE IF (compna(i) == 'o-xylene') THEN
    mm(i)     =           106.167
    parame(i,1)   = mm(i)*  .0295409161
    parame(i,2) =           3.76000631
    parame(i,3) =           291.049123
  ELSE IF (compna(i) == 'thf') THEN
    mm(i)       =   72.1057000000000     
    ! parame(i,1) = mm(i)*  0.34311391E-01
    parame(i,1) =   2.47404685540709     
    parame(i,2) =   3.51369375633677     
    parame(i,3) =   274.181927093696     
    parame(i,6) =   1.63100000000000     
  ELSE IF (compna(i) == 'co2') THEN
    mm(i)       = 44.01
    parame(i,1) = mm(i)*  .0470968503
    parame(i,2) = 2.7851954
    parame(i,3) = 169.207418
    IF (pol >= 1) parame(i,1) = mm(i)*  3.438191426159075E-002
    IF (pol >= 1) parame(i,2) = 3.18693935424469
    IF (pol >= 1) parame(i,3) = 163.333232725156
    IF (pol >= 1) parame(i,7) = 4.400000000000
    IF (pol >= 1) lli(i)      = 1.472215*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 17.706567
    IF (pol >= 1) chap(i)     = 0.5
    IF (pol == 2) parame(i,11)= 2.911
  ELSE IF (compna(i) == 'co') THEN
    IF (pol /= 1) write (*,*) 'parameters for co missing'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       = 28.01
    IF (pol >= 1) parame(i,1) = mm(i)*  5.126059746332587E-002  ! 1.43580933494776     
    IF (pol >= 1) parame(i,2) = 3.13556624711756     
    IF (pol >= 1) parame(i,3) = 87.7191028693595     
    IF (pol >= 1) parame(i,6) = 0.1098
  ELSE IF (compna(i) == 'n2') THEN
    mm(i)       = 28.01
    parame(i,1) = mm(i)*  .0430301713
    parame(i,2) = 3.3129702
    parame(i,3) = 90.9606924
    IF (pol >= 1) parame(i,1) = mm(i)*  3.971157114787596E-002
    IF (pol >= 1) parame(i,2) = 3.42116853868336
    IF (pol >= 1) parame(i,3) = 92.3972606842862
    IF (pol >= 1) parame(i,7) = 1.52000000000000
    IF (pol >= 1) lli(i)      = 1.5188*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 19.9247
    IF (pol >= 1) chap(i)     = 0.375
    ! better RGT-results came later, with: 1.5822   21.201   0.3972
  ELSE IF (compna(i) == 'o2') THEN
    mm(i)       = 32.05
    parame(i,1) = mm(i)*  .0353671563
    parame(i,2) = 3.19465166
    parame(i,3) = 114.430197
  ELSE IF (compna(i) == 'hydrogen') THEN
    mm(i)       = 2.016
    parame(i,1) = mm(i)*  .258951975
    parame(i,2) = 4.43304935
    parame(i,3) = 29.6509579
    
    mm(i)       = 2.016
    parame(i,1) = 1.0
    parame(i,2) = 2.915
    parame(i,3) = 38.0
    
    ! mm(i)       = 2.016  !  Ghosh et al. 2003
    ! parame(i,1) = 1.0
    ! parame(i,2) = 2.986
    ! parame(i,3) = 19.2775
  ELSE IF (compna(i) == 'argon') THEN
    ! mm(i)       = 39.948  ! adjusted m !!
    ! parame(i,1) = 0.9285
    ! parame(i,2) = 3.4784
    ! parame(i,3) = 122.23
    mm(i)       = 39.948 ! enforced m=1 !!
    parame(i,1) = 1.0
    parame(i,2) = 3.3658
    parame(i,3) = 118.34
    IF (pol == 2) parame(i,11)= 1.6411
  ELSE IF (compna(i) == 'xenon') THEN
    mm(i)       =   131.29
    parame(i,1) =   1.0
    parame(i,2) =   3.93143
    parame(i,3) =   227.749
  ELSE IF (compna(i) == 'chlorine') THEN  ! Cl2
    mm(i)       =   70.906
    parame(i,1) =   1.5514
    parame(i,2) =   3.3672
    parame(i,3) =   265.67
  ELSE IF (compna(i) == 'SF6') THEN
    mm(i)       = 146.056  ! adjusted m !!
    parame(i,1) = 2.48191
    parame(i,2) = 3.32727
    parame(i,3) = 161.639
    ! mm(i)       = 146.056 ! enforced m=1 !!
    ! parame(i,1) = 1.0
    ! parame(i,2) = 4.55222
    ! parame(i,3) = 263.1356
  ELSE IF (compna(i) == 'benzene') THEN
    mm(i)       = 78.114
    parame(i,1) = mm(i)*  .0315590546
    parame(i,2) = 3.64778975
    parame(i,3) = 287.354574
    IF (pol >= 1) mm(i)       = 78.114   ! PCP-SAFT with m=2 in QQ term
    IF (pol >= 1) parame(i,1) = mm(i)*  2.932783311E-2 !  = 2.29091435590515
    IF (pol >= 1) parame(i,2) =         3.7563854
    IF (pol >= 1) parame(i,3) =         294.06253
    IF (pol >= 1) parame(i,7) = 5.5907
  ELSE IF (compna(i) == 'ethylbenzene') THEN
    mm(i)       = 106.167
    parame(i,1) = mm(i)*  .0290120497
    parame(i,2) = 3.79741116
    parame(i,3) = 287.348098
    IF (pol == 2) parame(i,11)= 13.3
  ELSE IF (compna(i) == 'propylbenzene') THEN
    mm(i)     =           120.194
    parame(i,1)   = mm(i)*  .0278171627
    parame(i,2) =           3.8437772
    parame(i,3) =           288.128269
  ELSE IF (compna(i) == 'n-butylbenzene') THEN
    mm(i)       = 134.221
    parame(i,1) = mm(i)*  .0280642225
    parame(i,2) = 3.87267961
    parame(i,3) = 283.072331
  ELSE IF (compna(i) == 'tetralin') THEN
    mm(i)       = 132.205
    parame(i,1) = mm(i)*  .0250640795
    parame(i,2) = 3.87498866
    parame(i,3) = 325.065688
  ELSE IF (compna(i) == 'methylcyclohexane') THEN
    mm(i)       = 98.182
    parame(i,1) = mm(i)*  .0271259953
    parame(i,2) = 3.99931892
    parame(i,3) = 282.334148
    IF (pol == 2) parame(i,11)= 13.1
  ELSE IF (compna(i) == 'methylcyclopentane') THEN
    mm(i)       = 84.156
    parame(i,1) = mm(i)*  .0310459009
    parame(i,2) = 3.82534693
    parame(i,3) = 265.122799
  ELSE IF (compna(i) == 'acetone') THEN
    mm(i)       = 58.0800000000000     ! PC-SAFT
    parame(i,1) = mm(i)*  4.870380408159182E-002  ! =2.82871694105885
    parame(i,2) = 3.24969003020675
    parame(i,3) = 250.262241927379
    lli(i)      = 2.0021*parame(i,2)
    phi_criti(i)= 21.336
    chap(i)     = 0.24931
    IF (pol >= 1) mm(i)       = 58.0800000000000     ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.725811736856114E-002  ! =2.74475145676603
    IF (pol >= 1) parame(i,2) = 3.27423145271184
    IF (pol >= 1) parame(i,3) = 232.990879135326
    IF (pol >= 1) parame(i,6) = 2.88000000000000
    IF (pol >= 1) lli(i)      = 2.0641*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 28.1783
    IF (pol >= 1) chap(i)     = 0.22695
    IF (pol >= 2) mm(i)       = 58.0800000000000     ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.902301475689938E-002  !  =2.84725669708072
    IF (pol >= 2) parame(i,2) = 3.23880349104868
    IF (pol >= 2) parame(i,3) = 220.884202656054
    IF (pol >= 2) parame(i,6) = 2.88000000000000
    IF (pol == 2) parame(i,11)= 6.40000000000000
  ELSE IF (compna(i) == 'butanone') THEN
    mm(i)       = 72.1066               !  PC-SAFT
    parame(i,1) = mm(i)*  4.264192830122321E-002  !  =3.07476446724498
    parame(i,2) = 3.39324011060028
    parame(i,3) = 252.267273608975
    IF (pol >= 1) mm(i)       = 72.1066               !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.137668924230600E-002  !  =2.98353238051926
    IF (pol >= 1) parame(i,2) = 3.42393701353423
    IF (pol >= 1) parame(i,3) = 244.994381354681
    IF (pol >= 1) parame(i,6) = 2.78000000000000
    IF (pol >= 2) mm(i)       = 72.1066               !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.254697075199448E-002  !  =3.06791740122577
    IF (pol >= 2) parame(i,2) = 3.39138375903252
    IF (pol >= 2) parame(i,3) = 236.527763837528
    IF (pol >= 2) parame(i,6) = 2.78000000000000
    IF (pol == 2) parame(i,11)= 8.13000000000000
  ELSE IF (compna(i) == '2-pentanone') THEN
    ! mm(i)       = 86.134                !  PC-SAFT
    ! parame(i,1) = mm(i)*  3.982654501296355E-002  !  =3.43041962814660
    ! parame(i,2) = 3.46877976946838
    ! parame(i,3) = 249.834724442656
    ! mm(i)       = 86.134                !  PCP-SAFT
    ! parame(i,1) = mm(i)*  3.893594769994072E-002  !  =3.35370891918669
    ! parame(i,2) = 3.49417356096593
    ! parame(i,3) = 246.656329096835
    ! parame(i,6) = 2.70000000000000
    mm(i)       = 86.134                !  PCIP-SAFT
    parame(i,1) = mm(i)*  3.973160761515879E-002  !  =3.42224229032409
    parame(i,2) = 3.46827593107280
    parame(i,3) = 240.904278156822
    parame(i,6) = 2.70000000000000
    IF (pol == 2) parame(i,11)= 9.93000000000000
  ELSE IF (compna(i) == '3-pentanone') THEN
    mm(i)       = 86.134                !  PC-SAFT
    parame(i,1) = 3.36439508013322
    parame(i,2) = 3.48770251979329
    parame(i,3) = 252.695415552376
    IF (pol >= 1) mm(i)       = 86.134                !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = 3.27863398611842
    IF (pol >= 1) parame(i,2) = 3.51592571835030
    IF (pol >= 1) parame(i,3) = 248.690775540981
    IF (pol >= 1) parame(i,6) = 2.82000000000000
    IF (pol == 2) mm(i)       = 86.134                !  PCIP-SAFT
    IF (pol == 2) parame(i,1) = 3.34821857026283
    IF (pol == 2) parame(i,2) = 3.48903345340516
    IF (pol == 2) parame(i,3) = 242.314578558329
    IF (pol == 2) parame(i,6) = 2.82000000000000
    IF (pol == 2) parame(i,11)= 9.93000000000000
  ELSE IF (compna(i) == 'cyclohexanone') THEN          ! from DIPPR
    ! IF (pol.GE.1) mm(i)       = 98.1430                !  PCP-SAFT
    ! IF (pol.GE.1) parame(i,1) = 3.084202
    ! IF (pol.GE.1) parame(i,2) = 3.613681
    ! IF (pol.GE.1) parame(i,3) = 286.15865
    ! IF (pol.GE.1) parame(i,6) = 3.087862
    IF (pol >= 1) mm(i)       = 98.1500000000000
    IF (pol >= 1) parame(i,1) = 2.72291913132818
    IF (pol >= 1) parame(i,2) = 3.79018433908522
    IF (pol >= 1) parame(i,3) = 314.772193827344
    IF (pol >= 1) parame(i,6) = 3.24600000000000
    IF (pol /= 1) WRITE (*,*) 'no non-polar param. for cyclohexanone'
    IF (pol /= 1) STOP
  ELSE IF (compna(i) == 'propanal') THEN
    mm(i)       = 58.08                             !  PC-SAFT
    parame(i,1) = 2.67564746980910
    parame(i,2) = 3.26295953984941
    parame(i,3) = 251.888982765626
    IF (pol >= 1) mm(i)       = 58.08               !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = 2.60007872084995
    IF (pol >= 1) parame(i,2) = 3.28720732189761
    IF (pol >= 1) parame(i,3) = 235.205188090107
    IF (pol >= 1) parame(i,6) = 2.72000000000000
    IF (pol >= 2) mm(i)       = 58.08               !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 2.72471167411028
    IF (pol >= 2) parame(i,2) = 3.24781643022922
    IF (pol >= 2) parame(i,3) = 221.642071811094
    IF (pol >= 2) parame(i,6) = 2.72000000000000
    IF (pol >= 2) parame(i,11)= 6.50000000000000
  ELSE IF (compna(i) == 'butanal') THEN
    mm(i)       = 72.1066000000000                  ! PC-SAFT
    parame(i,1) = 2.96824823599784
    parame(i,2) = 3.44068916025889
    parame(i,3) = 253.929404992884
    IF (pol >= 1) mm(i)       = 72.1066000000000    ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = 2.86783706423953
    IF (pol >= 1) parame(i,2) = 3.47737904036296
    IF (pol >= 1) parame(i,3) = 247.543312127310
    IF (pol >= 1) parame(i,6) = 2.72000000000000
  ELSE IF (compna(i) == 'dmso') THEN
    mm(i)       = 78.1300000000000                  ! PC-SAFT
    parame(i,1) = 2.92225114054231
    parame(i,2) = 3.27780791606297
    parame(i,3) = 355.688793038512
    IF (pol >= 1) mm(i)       = 78.1300000000000    ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = 3.02433694138348
    IF (pol >= 1) parame(i,2) = 3.24270742566613
    IF (pol >= 1) parame(i,3) = 309.357476696679
    IF (pol >= 1) parame(i,6) = 3.96000000000000
    IF (pol >= 2) mm(i)       = 78.1300000000000    ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 3.19078234633277
    IF (pol >= 2) parame(i,2) = 3.19778269816832
    IF (pol >= 2) parame(i,3) = 286.337981216861
    IF (pol >= 2) parame(i,6) = 3.96000000000000
    IF (pol >= 2) parame(i,11)= 7.97000000000000
  ELSE IF (compna(i) == 'acetone_JC') THEN  ! Jog-Chapman
    ! mm(i)       = 58.0800000000000   ! Dominik et al.2005
    ! parame(i,1) = 2.221
    ! parame(i,2) = 3.607908
    ! parame(i,3) = 259.99
    ! parame(i,6) = 2.7
    ! parame(i,8) = 0.2258
    ! mm(i)       = 58.0800000000000
    ! parame(i,1) = mm(i)*  3.556617369195472E-002
    ! parame(i,2) = 3.58780367502515
    ! parame(i,3) = 273.025100470307
    ! parame(i,6) = 2.70000000000000
    ! parame(i,8) = 0.229800000000000
    
    mm(i)       = 58.08     !  Tumakaka Sadowski 2004
    parame(i,1) = mm(i)*  3.766E-2
    parame(i,2) = 3.6028
    parame(i,3) = 245.49
    parame(i,6) = 2.72
    parame(i,8) = 0.2969
    ! mm(i)       = 58.0800000000000     ! no adjust. DD-param.
    ! parame(i,1) = 1.87041620247774
    ! parame(i,2) = 3.79783535570774
    ! parame(i,3) = 208.885730881588
    ! parame(i,6) = 2.88000000000000
    ! parame(i,8) = 1.0/parame(i,1)
    WRITE (*,*) 'caution: parame(i,8) is now used for branching'
    STOP
    kij(1,2) = -0.005
    kij(2,1)=kij(1,2)
  ELSE IF (compna(i) == 'acetone_SF') THEN  ! Saager-Fischer
    mm(i)       = 58.08
    parame(i,1) = mm(i)*  4.603296414764944E-002
    parame(i,2) = 3.29454924451643
    parame(i,3) = 221.052649057645
    parame(i,6) = 2.70000000000000
    parame(i,8) = 0.625410000000000
    mm(i)       = 58.08 ! form as expected from me - no DD-param adjusted.dat
    parame(i,1) = mm(i)*  4.364264724158790E-002  !  =2.53476495179143
    parame(i,2) = 3.37098670735567
    parame(i,3) = 254.366379701851
    parame(i,6) = 2.88000000000000
    ! mm(i)       = 58.08  ! form as expected but w/ sumseg/1.5 - no DD-param adjusted.dat
    ! parame(i,1) = mm(i)*  4.694644361257521E-002  !  =2.72664944501837
    ! parame(i,2) = 3.27842292595463
    ! parame(i,3) = 238.398883501772
    ! parame(i,6) = 2.88000000000000
    ! mm(i)       = 58.08  ! form as expected but w/ sumseg/1.5 and fdd*sumseg- no DD-param adjusted.dat
    ! parame(i,1) = mm(i)*  4.458214655521766E-002  !  =2.58933107192704
    ! parame(i,2) = 3.32050824493493
    ! parame(i,3) = 218.285994651271
    ! parame(i,6) = 2.88000000000000
    WRITE (*,*) 'caution: parame(i,8) is now used for branching'
    STOP
    kij(1,2) = 0.035
    kij(2,1)=kij(1,2)
  ELSE IF (compna(i) == 'ethylacetate_JC') THEN  ! Jog-Chapman
    ! mm(i)       = 88.11
    ! parame(i,1) = 2.7481
    ! parame(i,2) = 3.6511
    ! parame(i,3) = 236.99
    ! parame(i,6) = 1.84
    ! parame(i,8) = 0.5458
    mm(i)       = 88.1060000000000
    parame(i,1) = mm(i)*  0.03117 ! 2.74626402
    parame(i,2) = 3.6493
    parame(i,3) = 236.75
    parame(i,6) = 1.8
    parame(i,8) = 0.5462
  ELSE IF (compna(i) == 'ethylacetate_SF') THEN  ! Saager-Fischer
    mm(i)       = 88.106
    parame(i,1) = mm(i)*  3.564165384763394E-002
    parame(i,2) = 3.447379322
    parame(i,3) = 226.0930487
    parame(i,6) = 1.8
    parame(i,8) = 0.849967000000000
    WRITE (*,*) 'caution: parame(i,8) is now used for branching'
    STOP
  ELSE IF (compna(i) == '12po_JC') THEN  ! Jog-Chapman
    mm(i)       = 58.08
    parame(i,1) = 2.0105
    parame(i,2) = 3.6095
    parame(i,3) = 258.82
    parame(i,6) = 2.0
    parame(i,8) = 0.3979
    WRITE (*,*) 'caution: parame(i,8) is now used for branching'
    STOP
  ELSE IF (compna(i) == '12po_SF') THEN  ! Saager-Fischer
    mm(i)       = 58.08
    parame(i,1) = 2.1341
    parame(i,2) = 3.4739
    parame(i,3) = 252.95
    parame(i,6) = 2.0
    parame(i,8) = 0.916
    WRITE (*,*) 'caution: parame(i,8) is now used for branching'
    STOP
  ELSE IF (compna(i) == 'acrylonitrile') THEN
    IF (pol >= 2) mm(i)       = 53.06     ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 2.168
    IF (pol >= 2) parame(i,2) = 3.575
    IF (pol >= 2) parame(i,3) = 214.83
    IF (pol >= 2) parame(i,6) = 3.91
    IF (pol == 2) parame(i,11)= 8.04
    IF (pol >= 2) mm(i)       = 53.0000000000000  ! second parameter set ??
    IF (pol >= 2) parame(i,1) = 2.45403467006041
    IF (pol >= 2) parame(i,2) = 3.41276825781723
    IF (pol >= 2) parame(i,3) = 195.194353082408
    IF (pol >= 2) parame(i,6) = 3.91000000000000
    IF (pol == 2) parame(i,11)= 8.04000000000000
  ELSE IF (compna(i) == 'butyronitrile') THEN
    ! mm(i)       = 69.11
    ! parame(i,1) = 2.860
    ! parame(i,2) = 3.478
    ! parame(i,3) = 253.21
    ! parame(i,6) = 4.07
    mm(i)       = 69.11
    parame(i,1) = 2.989
    parame(i,2) = 3.441
    parame(i,3) = 234.04
    parame(i,6) = 4.07
    IF (pol == 2) parame(i,11)= 8.4
  ELSE IF (compna(i) == 'propionitrile') THEN
    mm(i)       = 55.079                !  PC-SAFT
    parame(i,1) = 2.66211021227108
    parame(i,2) = 3.34032231132738
    parame(i,3) = 294.078737359580
    IF (pol >= 1) mm(i)       = 55.079                !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = 2.50958981615666
    IF (pol >= 1) parame(i,2) = 3.39806320429568
    IF (pol >= 1) parame(i,3) = 239.152759066148
    IF (pol >= 1) parame(i,6) = 4.05000000000000
    IF (pol >= 2) mm(i)       = 55.079                !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 2.54684827683436
    IF (pol >= 2) parame(i,2) = 3.41240089912190
    IF (pol >= 2) parame(i,3) = 218.299491580335
    IF (pol >= 2) parame(i,6) = 4.05000000000000
    IF (pol == 2) parame(i,11)= 6.24000000000000
    ! IF (pol.GE.2) mm(i)       = 55.079                !  PCIP-SAFT my_DD adjusted
    ! IF (pol.GE.2) parame(i,1) = 2.61175151088002
    ! IF (pol.GE.2) parame(i,2) = 3.37194293181453
    ! IF (pol.GE.2) parame(i,3) = 233.346110749402
    ! IF (pol.GE.2) parame(i,6) = 3.74682245835235
    ! IF (pol.EQ.2) parame(i,11)= 6.24000000000000
  ELSE IF (compna(i) == 'nitromethane') THEN
    mm(i)       = 61.04                 !  PC-SAFT
    parame(i,1) = mm(i)*  4.233767489308791E-002  !  =2.58429167547409
    parame(i,2) = 3.10839592337018
    parame(i,3) = 310.694151426943
    IF (pol >= 1) mm(i)       = 61.04                 !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.191475020685036E-002  !  =2.55847635262615
    IF (pol >= 1) parame(i,2) = 3.10129282495975
    IF (pol >= 1) parame(i,3) = 256.456941430554
    IF (pol >= 1) parame(i,6) = 3.46000000000000
    IF (pol >= 2) mm(i)       = 61.04                 !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.394323357988009E-002  !  =2.68229497771588
    IF (pol >= 2) parame(i,2) = 3.10654492320028
    IF (pol >= 2) parame(i,3) = 225.973607468282
    IF (pol >= 2) parame(i,6) = 3.46000000000000
    IF (pol >= 2) parame(i,11)= 7.37000000000000
  ELSE IF (compna(i) == 'nitroethane') THEN
    mm(i)       = 75.067                !  PC-SAFT
    parame(i,1) = mm(i)*  4.019977215251163E-002  !  =3.01767629617259
    parame(i,2) = 3.21364231060938
    parame(i,3) = 286.571650044235
    IF (pol >= 1) mm(i)       = 75.067                !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  3.928506808347654E-002  !  =2.94901220582233
    IF (pol >= 1) parame(i,2) = 3.23117331990738
    IF (pol >= 1) parame(i,3) = 265.961000131109
    IF (pol >= 1) parame(i,6) = 3.23000000000000
    IF (pol >= 2) mm(i)       = 75.067                !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.117677400894779E-002  !  =3.09101689452968
    IF (pol >= 2) parame(i,2) = 3.19364569858756
    IF (pol >= 2) parame(i,3) = 246.676040248662
    IF (pol >= 2) parame(i,6) = 3.23000000000000
    IF (pol >= 2) parame(i,11)= 9.63000000000000
  ELSE IF (compna(i) == 'acetonitrile') THEN
    mm(i)       = 41.052                !  PC-SAFT
    parame(i,1) = mm(i)*  5.673187410405271E-002  !  =2.32895689571957
    parame(i,2) = 3.18980108373791
    parame(i,3) = 311.307486044181
    IF (pol >= 1) mm(i)       = 41.052                !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  5.254832931037250E-002  !  =2.15721401484941
    IF (pol >= 1) parame(i,2) = 3.27301469369132
    IF (pol >= 1) parame(i,3) = 216.888948676921
    IF (pol >= 1) parame(i,6) = 3.92520000000000
    IF (pol >= 2) mm(i)       = 41.052                !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  5.125846581157176E-002  !  =2.10426253849664
    IF (pol >= 2) parame(i,2) = 3.39403305120647
    IF (pol >= 2) parame(i,3) = 199.070191065791
    IF (pol >= 2) parame(i,6) = 3.92520000000000
    IF (pol >= 2) parame(i,11)= 4.40000000000000
    ! IF (pol >= 2) mm(i)       = 41.052                !  PCIP-SAFT my_DD adjusted
    ! IF (pol >= 2) parame(i,1) = mm(i)*  5.755347845863738E-002  !  =2.36268539768398
    ! IF (pol >= 2) parame(i,2) = 3.18554306395900
    ! IF (pol >= 2) parame(i,3) = 225.143934506015
    ! IF (pol >= 2) parame(i,6) = 3.43151866932598
    ! IF (pol >= 2) parame(i,11)= 4.40000000000000
    ! mm(i)       = 41.053                   ! PCP-SAFT dipole and quadrupole
    ! parame(i,1) = 1.79993
    ! parame(i,2) = 3.47366
    ! parame(i,3) = 211.001
    ! parame(i,6) = 3.93800
    ! parame(i,7) = 2.44000
    ! parame(i,11)= 0.0
  ELSE IF (compna(i) == 'dmf') THEN
    mm(i)       = 73.09              !  PC-SAFT
    parame(i,1) = 2.388
    parame(i,2) = 3.658
    parame(i,3) = 363.77
    IF (pol >= 1) mm(i)       = 73.09              !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = 2.269
    IF (pol >= 1) parame(i,2) = 3.714
    IF (pol >= 1) parame(i,3) = 331.56
    IF (pol >= 1) parame(i,6) = 3.82
    IF (pol >= 2) mm(i)       = 73.09              !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 2.375
    IF (pol >= 2) parame(i,2) = 3.667
    IF (pol >= 2) parame(i,3) = 308.42
    IF (pol >= 2) parame(i,6) = 3.82
    IF (pol >= 2) parame(i,11)= 7.81
  ELSE IF (compna(i) == 'chloroform') THEN
    mm(i)       = 119.377              !  PCIP-SAFT
    parame(i,1) = 2.5957
    parame(i,2) = 3.4299
    parame(i,3) = 264.664
    parame(i,6) = 1.04
    IF (pol == 2) parame(i,11)= 8.23
  ELSE IF (compna(i) == 'dimethyl-ether') THEN
    mm(i)       = 46.069              ! PC-SAFT
    parame(i,1) = mm(i)*  0.049107715  !  =2.26234331
    parame(i,2) = 3.276640534
    parame(i,3) = 212.9343244
    IF (pol >= 1) mm(i)       = 46.0690000000000     ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  0.048170452  !  =2.219164566
    IF (pol >= 1) parame(i,2) = 3.296939638
    IF (pol >= 1) parame(i,3) = 212.1048888
    IF (pol >= 1) parame(i,6) = 1.30000000000000
    IF (pol >= 2) mm(i)       = 46.0690000000000     ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.939183716945787E-002  !  =2.27543254655976
    IF (pol >= 2) parame(i,2) = 3.26584718800835
    IF (pol >= 2) parame(i,3) = 206.904551967059
    IF (pol >= 2) parame(i,6) = 1.30000000000000
    IF (pol == 2) parame(i,11)= 5.29000000000000
  ELSE IF (compna(i) == 'methyl-ethyl-ether') THEN
    mm(i)     =           60.096
    parame(i,1) = mm(i)*  .0442404671
    parame(i,2) =           3.37282595
    parame(i,3) =           216.010217
    IF (pol >= 1) mm(i)       =         60.096            ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.3971676124088D-002  !  =2.64252184835325
    IF (pol >= 1) parame(i,2) =         3.37938465390
    IF (pol >= 1) parame(i,3) =         215.787173860
    IF (pol >= 1) parame(i,6) =         1.17000000000
    IF (pol >= 2) mm(i)       =         60.096            ! PICP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.4580196137984D-002  !  =2.67909146710834
    IF (pol >= 2) parame(i,2) =         3.36105342286
    IF (pol >= 2) parame(i,3) =         212.871911999
    IF (pol >= 2) parame(i,6) =         1.17000000000
    IF (pol >= 2) parame(i,11) =         7.93000000000
  ELSE IF (compna(i) == 'diethyl-ether') THEN
    mm(i)       = 74.123            ! PC-SAFT
    parame(i,1) = mm(i)*  .0409704089
    parame(i,2) = 3.48569553
    parame(i,3) = 217.64113
    IF (pol >= 1) mm(i)       = 74.123            ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.0103121403686E-2  !  =2.97256367
    IF (pol >= 1) parame(i,2) = 3.51268687697978
    IF (pol >= 1) parame(i,3) = 219.527376572135
    IF (pol >= 1) parame(i,6) = 1.15000000000000
    IF (pol >= 2) mm(i)       = 74.123            ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.04144179873E-2    !  =2.9956379
    IF (pol >= 2) parame(i,2) = 3.501724569
    IF (pol >= 2) parame(i,3) = 217.8941822
    IF (pol >= 2) parame(i,6) = 1.15
    IF (pol == 2) parame(i,11)= 8.73
  ELSE IF (compna(i) == 'vinylacetate') THEN
    mm(i)       = 86.092
    parame(i,1) = mm(i)*  .0374329292
    parame(i,2) = 3.35278602
    parame(i,3) = 240.492049
  ELSE IF (compna(i) == 'chloromethane') THEN    ! R40
    mm(i)       = 50.488        ! PC-SAFT
    parame(i,1) = mm(i)*  0.039418879  !  1.9902
    parame(i,2) = 3.1974
    parame(i,3) = 237.27
    IF (pol >= 1) mm(i)       = 50.488        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  0.035790801  !  1.8070
    IF (pol >= 1) parame(i,2) = 3.3034
    IF (pol >= 1) parame(i,3) = 229.97
    IF (pol >= 1) parame(i,6) = 1.8963
    IF (pol >= 1) lli(i)      = 1.67703*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 20.75417
    IF (pol >= 1) chap(i)     = 0.5
    IF (pol >= 2) mm(i)       = 50.488        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  3.68559992E-2  !  1.86078
    IF (pol >= 2) parame(i,2) = 3.275186
    IF (pol >= 2) parame(i,3) = 216.4621
    IF (pol >= 2) parame(i,6) = 1.8963
    IF (pol == 2) parame(i,11)= 4.72
  ELSE IF (compna(i) == 'fluoromethane') THEN    ! R41
    IF (pol /= 1) write (*,*) 'non-polar parameters missing for fluoromethane'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       =   34.0329000000000     
    IF (pol >= 1) parame(i,1) =   1.94494757526896     
    IF (pol >= 1) parame(i,2) =   2.96858005012635     
    IF (pol >= 1) parame(i,3) =   168.938697391009     
    IF (pol >= 1) parame(i,6) =         1.57823038894029     
  ELSE IF (compna(i) == 'dichloromethane') THEN   ! R30
    mm(i)       = 84.932        ! PC-SAFT
    parame(i,1) = 2.3117
    parame(i,2) = 3.3161
    parame(i,3) = 270.98
    IF (pol >= 1) mm(i)       = 84.932        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = 2.2687
    IF (pol >= 1) parame(i,2) = 3.3373
    IF (pol >= 1) parame(i,3) = 269.08
    IF (pol >= 1) parame(i,6) = 1.6
    IF (pol >= 2) mm(i)       = 84.932        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = 2.3435
    IF (pol >= 2) parame(i,2) = 3.2987
    IF (pol >= 2) parame(i,3) = 260.66
    IF (pol >= 2) parame(i,6) = 1.6
    IF (pol == 2) parame(i,11)= 6.48
  ELSE IF (compna(i) == 'difluoromethane') THEN   ! R32
    IF (pol /= 1) write (*,*) 'non-polar parameters missing for difluoromethane'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       = 52.0236             ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.814700934384165E-002  !  2.50478075530028
    IF (pol >= 1) parame(i,2) = 2.79365980535456
    IF (pol >= 1) parame(i,3) = 160.893555378523
    IF (pol >= 1) parame(i,6) = 1.97850000000000
  ELSE IF (compna(i) == 'trifluoromethane') THEN   ! R23
    IF (pol /= 1) write (*,*) 'non-polar parameters missing for trifluoromethane'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       =   70.0138000000000     
    IF (pol >= 1) parame(i,1) =   2.66039274225485     
    IF (pol >= 1) parame(i,2) =   2.82905884530501     
    IF (pol >= 1) parame(i,3) =   149.527709542333     
    IF (pol >= 1) parame(i,6) =   1.339963415253999E-002
  ELSE IF (compna(i) == 'tetrachloromethane') THEN   ! R10
    mm(i)       = 153.822
    parame(i,1) = mm(i)*  .0150432213
    parame(i,2) = 3.81801454
    parame(i,3) = 292.838632
  ELSE IF (compna(i) == 'trichlorofluoromethane') THEN   ! R11
    IF (pol /= 1) write (*,*) 'non-polar parameters missing for trichlorofluoromethane'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       =   137.368000000000     
    IF (pol >= 1) parame(i,1) =   2.28793359008803     
    IF (pol >= 1) parame(i,2) =   3.69013104930876     
    IF (pol >= 1) parame(i,3) =   248.603173885090     
    IF (pol >= 1) parame(i,6) =   0.23225538492979
  ELSE IF (compna(i) == 'chlorodifluoromethane') THEN   ! R22   ( CHClF2 or CHF2Cl)
    IF (pol /= 1) write (*,*) 'non-polar parameters missing for chlorodifluoromethane'
    IF (pol /= 1) stop
    IF (pol >= 1) mm(i)       =   86.4684000000000     
    IF (pol >= 1) parame(i,1) =   2.47218586047893     
    IF (pol >= 1) parame(i,2) =   3.13845692489930     
    IF (pol >= 1) parame(i,3) =   187.666355083434     
    IF (pol >= 1) parame(i,6) =         1.04954264812860     
  ELSE IF (compna(i) == 'chloroethane') THEN
    mm(i)       = 64.514
    parame(i,1) = mm(i)*  .0350926868
    parame(i,2) = 3.41602397
    parame(i,3) = 245.42626
  ELSE IF (compna(i) == '11difluoroethane') THEN
    ! mm(i)       = 66.0500000000000         !  PC-SAFT
    ! parame(i,1) = mm(i)*  4.109944338817734E-002
    ! parame(i,2) = 3.10257444633546
    ! parame(i,3) = 192.177159144029
    ! mm(i)      = 66.05                             !  PC-SAFT assoc
    ! parame(i,1)= 2.984947188
    ! parame(i,2)= 2.978630027
    ! parame(i,3)= 137.8192282
    ! nhb_typ(i) = 2
    ! nhb_no(i,1)= 1
    ! nhb_no(i,2)= 1
    ! eps_hb(i,i,1,2)= 823.3478288
    ! eps_hb(i,i,2,1)= 823.3478288
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 0.96345994
    IF (pol >= 1) mm(i)       = 66.0500000000000        !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  3.949665745363346E-002  !  =2.60875422481249
    IF (pol >= 1) parame(i,2) = 3.13758353925036
    IF (pol >= 1) parame(i,3) = 179.517952627836
    IF (pol >= 1) parame(i,6) = 2.27000000000000
    IF (pol >= 1) lli(i)      = 2.03907*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 26.5
    IF (pol >= 1) chap(i)     = 0.4
    IF (pol >= 2) mm(i)       = 66.0500000000000        !  PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.093647666154238E-002  !  =2.70385428349487
    IF (pol >= 2) parame(i,2) = 3.10437129415885
    IF (pol >= 2) parame(i,3) = 170.464400902455
    IF (pol >= 2) parame(i,6) = 2.27000000000000
    IF (pol == 2) parame(i,11)= 5.01000000000000
  ELSE IF (compna(i) == '1-chlorobutane') THEN
    mm(i)       = 92.568
    parame(i,1) = mm(i)*  .0308793201
    parame(i,2) = 3.64240187
    parame(i,3) = 258.655298
  ELSE IF (compna(i) == 'chlorobenzene') THEN
    ! mm(i)       = 112.558
    ! parame(i,1) = mm(i)*  .0235308686
    ! parame(i,2) = 3.75328494
    ! parame(i,3) = 315.039018
    mm(i)       = 112.558        ! PCIP-SAFT
    parame(i,1) = mm(i)*  0.023824167  !  =2.6816
    parame(i,2) = 3.7352
    parame(i,3) = 308.82
    parame(i,6) = 1.69
    IF (pol == 2) parame(i,11)= 14.1
  ELSE IF (compna(i) == 'styrene') THEN
    mm(i)       = 104.150
    parame(i,1) = mm(i)*  2.9124104853E-2
    parame(i,2) = 3.760233548
    parame(i,3) = 298.51287564
  ELSE IF (compna(i) == 'methylmethanoate') THEN
    mm(i)     = 60.053
    parame(i,1) = mm(i)*  .0446000264
    parame(i,2) = 3.08753499
    parame(i,3) = 242.626755
    IF (pol >= 1) mm(i)       = 60.053        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.366991153963102E-002  !  =2.62250919768946
    IF (pol >= 1) parame(i,2) = 3.10946396964
    IF (pol >= 1) parame(i,3) = 239.051951942
    IF (pol >= 1) parame(i,6) = 1.77
    IF (pol >= 2) mm(i)       = 60.053        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.492572388931002E-2  !  2.69792449
    IF (pol >= 2) parame(i,2) = 3.078467837
    IF (pol >= 2) parame(i,3) = 232.1842551
    IF (pol >= 2) parame(i,6) = 1.77
    IF (pol == 2) parame(i,11)= 5.05
  ELSE IF (compna(i) == 'ethylmethanoate') THEN
    mm(i)     =           74.079        ! PC-SAFT
    parame(i,1) = mm(i)*  .03898009
    parame(i,2) =           3.31087192
    parame(i,3) =           246.465646
    IF (pol >= 1) mm(i)       =         74.079        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  3.825407152074255E-002  !  =2.83382336418509
    IF (pol >= 1) parame(i,2) =         3.33160046679
    IF (pol >= 1) parame(i,3) =         244.495680932
    IF (pol >= 1) parame(i,6) =         1.93000000000
  ELSE IF (compna(i) == 'propylmethanoate') THEN
    mm(i)     =           88.106
    parame(i,1) = mm(i)*  .0364206062
    parame(i,2) =           3.41679642
    parame(i,3) =           246.457732
    IF (pol >= 1) mm(i)       =         88.106
    IF (pol >= 1) parame(i,1) = mm(i)*  3.60050739149E-2  !  =3.17226304235232
    IF (pol >= 1) parame(i,2) =         3.42957609309
    IF (pol >= 1) parame(i,3) =         245.637644107
    IF (pol >= 1) parame(i,6) =         1.89
  ELSE IF (compna(i) == 'methylacetate') THEN
    mm(i)       = 74.079        ! PC-SAFT
    parame(i,1) = mm(i)*  4.286817177E-2  !  =3.175631296
    parame(i,2) = 3.18722021277843
    parame(i,3) = 234.106931032456
    IF (pol >= 1) mm(i)       = 74.079        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.228922065E-2  !  =3.132743176
    IF (pol >= 1) parame(i,2) = 3.2011401688
    IF (pol >= 1) parame(i,3) = 233.17562886
    IF (pol >= 1) parame(i,6) = 1.72
    IF (pol >= 2) mm(i)       = 74.079        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  4.298900538E-2  !  =3.18458252
    IF (pol >= 2) parame(i,2) = 3.180642322
    IF (pol >= 2) parame(i,3) = 229.3132680
    IF (pol >= 2) parame(i,6) = 1.72
    IF (pol == 2) parame(i,11)= 6.94
  ELSE IF (compna(i) == 'ethylacetate') THEN
    mm(i)       = 88.106        ! PC-SAFT
    parame(i,1) = mm(i)*  .0401464427   !  =3.537142481
    parame(i,2) = 3.30789258
    parame(i,3) = 230.800689
    IF (pol >= 1) mm(i)       = 88.106        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  0.039792575   !  =3.505964572
    IF (pol >= 1) parame(i,2) = 3.317655188
    IF (pol >= 1) parame(i,3) = 230.2434769
    IF (pol >= 1) parame(i,6) = 1.78
    IF (pol >= 2) mm(i)       = 88.106        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  0.040270267   !  =3.548052143
    IF (pol >= 2) parame(i,2) = 3.302097562
    IF (pol >= 2) parame(i,3) = 227.5026191
    IF (pol >= 2) parame(i,6) = 1.78
    IF (pol == 2) parame(i,11)= 8.62
  ELSE IF (compna(i) == 'ethyl-propanoate') THEN
    mm(i)       = 102.133
    parame(i,1) = mm(i)*  .0375692464
    parame(i,2) = 3.40306071
    parame(i,3) = 232.778374
  ELSE IF (compna(i) == 'propyl-ethanoate') THEN
    mm(i)       = 102.133
    parame(i,1) = mm(i)*  .0370721275
    parame(i,2) = 3.42272266
    parame(i,3) = 235.758378
    IF (pol >= 1) mm(i)       = 102.133        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  3.687149995200072E-2  !  =3.76579690459769
    IF (pol >= 1) parame(i,2) = 3.4289353421006
    IF (pol >= 1) parame(i,3) = 235.41679442910
    IF (pol >= 1) parame(i,6) = 1.78
    ! IF (pol.EQ.2) parame(i,11)= 10.41
  ELSE IF (compna(i) == 'nbutyl-ethanoate') THEN
    mm(i)       = 116.16        ! PC-SAFT
    parame(i,1) = mm(i)*  .03427004
    parame(i,2) = 3.54269638
    parame(i,3) = 242.515768
    IF (pol >= 1) mm(i)       = 116.16        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  3.411585209773470E-002  !  =3.96289737967286
    IF (pol >= 1) parame(i,2) = 3.54821589228130
    IF (pol >= 1) parame(i,3) = 242.274388267447
    IF (pol >= 1) parame(i,6) = 1.87000000000000
    IF (pol >= 2) mm(i)       = 116.16        ! PCIP-SAFT
    IF (pol >= 2) parame(i,1) = mm(i)*  3.442139015733717E-002  !  =3.99838868067629
    IF (pol >= 2) parame(i,2) = 3.53576054452119
    IF (pol >= 2) parame(i,3) = 240.154409609249
    IF (pol >= 2) parame(i,6) = 1.87000000000000
    IF (pol == 2) parame(i,11)= 14.2000000000000
  ELSE IF (compna(i) == 'methyl-octanoate') THEN
    mm(i)       = 158.24        ! PC-SAFT
    parame(i,1) = 5.2074
    parame(i,2) = 3.6069
    parame(i,3) = 244.12
  ELSE IF (compna(i) == 'methyl-decanoate') THEN
    mm(i)       = 186.2912        ! PC-SAFT
    parame(i,1) = 5.8402
    parame(i,2) = 3.6871
    parame(i,3) = 248.27
    
    mm(i)       = 186.2912        ! PC-SAFT from GC-method Tim
    parame(i,1) = 7.716
    parame(i,2) = 3.337303029
    parame(i,3) = 204.250907
    
    mm(i)       = 186.2912        ! PC-SAFT from GC-method (tightly fit) Tim
    parame(i,1) = 7.728
    parame(i,2) = 3.334023322
    parame(i,3) = 206.9099379
    
    ! mm(i)       = 186.2912        ! PC-SAFT from fit to DIPPR
    ! parame(i,1) = 6.285005
    ! parame(i,2) = 3.594888
    ! parame(i,3) = 236.781461
    ! ! parame(i,6) = 2.08056

    ! mm(i)       =         186.291000000000
    ! parame(i,1) =    6.28500485898895
    ! parame(i,2) =         3.59488828061149
    ! parame(i,3) =         236.781461491921
    ! parame(i,6) =         2.08055996894836
    ! parame(i,8) =         1.00000000000000
    mm(i)       =         186.291000000000
    parame(i,1) =    6.14436331493372
    parame(i,2) =         3.61977264863944
    parame(i,3) =         242.071887817656
    
  ELSE IF (compna(i) == 'methyl-dodecanoate') THEN
    mm(i)       = 214.344        ! PC-SAFT
    parame(i,1) = 6.5153
    parame(i,2) = 3.7406
    parame(i,3) = 250.7
  ELSE IF (compna(i) == 'methyl-tetradecanoate') THEN
    mm(i)       = 242.398        ! PC-SAFT
    parame(i,1) = 7.1197
    parame(i,2) = 3.7968
    parame(i,3) = 253.77
  ELSE IF (compna(i) == 'methyl-hexadecanoate') THEN
    mm(i)       = 270.451        ! PC-SAFT
    parame(i,1) = 7.891
    parame(i,2) = 3.814
    parame(i,3) = 253.71
  ELSE IF (compna(i) == 'methyl-octadecanoate') THEN
    mm(i)       = 298.504        ! PC-SAFT
    parame(i,1) = 8.8759
    parame(i,2) = 3.7932
    parame(i,3) = 250.81
  ELSE IF (compna(i) == 'CH2F2') THEN
    mm(i)       = 52.02
    parame(i,1) = 3.110084171
    parame(i,2) = 2.8145230485
    parame(i,3) = 158.98060151
  ELSE IF (compna(i) == 'naphthalene') THEN
    ! mm(i)       = 128.174000000
    ! parame(i,1) = mm(i)*  2.4877834216412E-2
    ! parame(i,2) = 3.82355815011
    ! parame(i,3) = 341.560675334
    
    mm(i)       = 128.17400000000
    parame(i,1) = mm(i)*  2.6400924157729E-2
    parame(i,2) = 3.8102186020014
    parame(i,3) = 328.96792935903
  ELSE IF (compna(i) == 'h2s') THEN
    mm(i)       =         34.0820000000000     ! PC-SAFT
    parame(i,1) = mm(i)*  4.838886696385162E-002   ! =    1.64918936386199
    parame(i,2) =         3.05478289838459
    parame(i,3) =         229.838873939562
    nhb_typ(i)     =              2
    nhb_no(i,1)    =              1
    nhb_no(i,2)    =              1
    eps_hb(i,i,1,2)=      536.634834731413
    eps_hb(i,i,2,1)=      536.634834731413
    eps_hb(i,i,1,1)=  0.0
    eps_hb(i,i,2,2)=  0.0
    kap_hb(i,i)    =     1.000000000000000E-003
! PC-SAFT from Xiaohua
    mm(i)       =         34.082     ! PC-SAFT
    parame(i,1) = 1.63677
    parame(i,2) =         3.06565
    parame(i,3) =         230.2121
    nhb_typ(i)     =              2
    nhb_no(i,1)    =              1
    nhb_no(i,2)    =              1
    eps_hb(i,i,1,2)=      275.1088
    eps_hb(i,i,2,1)=      275.1088
    eps_hb(i,i,1,1)=  0.0
    eps_hb(i,i,2,2)=  0.0
    kap_hb(i,i)    =     1.E-2
    ! IF (pol.GE.1) mm(i)       =         34.082       ! PCP-SAFT with quadrupole
    ! IF (pol.GE.1) parame(i,1) = mm(i)*  3.03171032558E-2  !  =1.03326751316478
    ! IF (pol.GE.1) parame(i,2) =         3.6868189914018
    ! IF (pol.GE.1) parame(i,3) =         246.862831266172
    ! IF (pol.GE.1) nhb_typ(i)     =              2
    ! IF (pol.GE.1) nhb_no(i,1)    =              1
    ! IF (pol.GE.1) nhb_no(i,2)    =              1
    ! IF (pol.GE.1) eps_hb(i,i,1,2)=      987.4927232
    ! IF (pol.GE.1) eps_hb(i,i,2,1)=      987.4927232
    ! IF (pol.GE.1) eps_hb(i,i,1,1)=  0.0
    ! IF (pol.GE.1) eps_hb(i,i,2,2)=  0.0
    ! IF (pol.GE.1) kap_hb(i,i)    =     5.5480659623d-4
    ! IF (pol.GE.1) parame(i,6) =        0.97833
    ! IF (pol.GE.1) parame(i,7) =         3.8623
    ! IF (pol.GE.1) LLi(i)      = 1.2737*parame(i,2)
    ! IF (pol.GE.1) phi_criti(i)= 14.316
    ! IF (pol.GE.1) chap(i)     = 0.4473
    IF (pol >= 1) mm(i)       =         34.0820000000000       ! PCP-SAFT no quadrupoLE
    IF (pol >= 1) parame(i,1) = mm(i)*  4.646468487062725E-002 ! 1.58360938976072
    IF (pol >= 1) parame(i,2) =         3.10111012646306
    IF (pol >= 1) parame(i,3) =         230.243457544889
    IF (pol >= 1) nhb_typ(i)     =              2
    IF (pol >= 1) nhb_no(i,1)    =              1
    IF (pol >= 1) nhb_no(i,2)    =              1
    IF (pol >= 1) eps_hb(i,i,1,2)=      584.367708701411
    IF (pol >= 1) eps_hb(i,i,2,1)=      584.367708701411
    IF (pol >= 1) eps_hb(i,i,1,1)=  0.0
    IF (pol >= 1) eps_hb(i,i,2,2)=  0.0
    IF (pol >= 1) kap_hb(i,i)    =     1.000000000000000E-003
    IF (pol >= 1) parame(i,6) =        0.978330000000000
    
    IF (pol >= 1) lli(i)      = 1.2737*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 14.316
    IF (pol >= 1) chap(i)     = 0.4473
    
    
    IF (pol == 2) parame(i,7) =        0.0
    IF (pol == 2) mm(i)       =         34.0820000000000     ! PCIP-SAFT
    IF (pol == 2) parame(i,1) = mm(i)*  4.806418212963168E-002 ! 1.63812345534211
    IF (pol == 2) parame(i,2) =         3.06556006883749
    IF (pol == 2) parame(i,3) =         221.746622243054
    IF (pol == 2) nhb_typ(i)     =              2
    IF (pol == 2) nhb_no(i,1)    =              1
    IF (pol == 2) nhb_no(i,2)    =              1
    IF (pol == 2) eps_hb(i,i,1,2)=      672.164783984789
    IF (pol == 2) eps_hb(i,i,2,1)=      672.164783984789
    IF (pol == 2) eps_hb(i,i,1,1)=  0.0
    IF (pol == 2) eps_hb(i,i,2,2)=  0.0
    IF (pol == 2) kap_hb(i,i)    =     1.000000000000000E-003
    IF (pol == 2) parame(i,6) =        0.978330000000000
    IF (pol == 2) parame(i,11) =         3.60200000000000
    IF (pol == 2) parame(i,7) =        0.0
    
    IF (pol >= 1)mm(i)       =    34.0820000000000     !PCP-SAFT D&Q
    IF (pol >= 1)parame(i,1) = mm(i)*  3.974667896078737E-002  !  = 1.35464631234155
    IF (pol >= 1)parame(i,2) =         3.30857082333438
    IF (pol >= 1)parame(i,3) =         234.248947273191
    IF (pol >= 1)nhb_typ(i)     =              2
    IF (pol >= 1)nhb_no(i,1)    =              1
    IF (pol >= 1)nhb_no(i,2)    =              1
    IF (pol >= 1)eps_hb(i,i,1,2)=      780.770936834770
    IF (pol >= 1)eps_hb(i,i,2,1)=      780.770936834770
    IF (pol >= 1)eps_hb(i,i,1,1)=  0.0
    IF (pol >= 1)eps_hb(i,i,2,2)=  0.0
    IF (pol >= 1)kap_hb(i,i)    =     1.000000000000000E-003
    IF (pol >= 1)parame(i,6) =        0.978330000000000
    IF (pol >= 1)parame(i,7) =         2.93750500000000
    
  ELSE IF (compna(i) == 'methanol') THEN
    mm(i)       = 32.042        ! PC-SAFT
    parame(i,1) = mm(i)*  .0476100379
    parame(i,2) = 3.23000005
    parame(i,3) = 188.904644
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2899.49055
    eps_hb(i,i,2,1)= 2899.49055
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0351760892
    IF (pol >= 1) mm(i)       = 32.042       ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  7.213091821E-2  !  =2.31121888139672
    IF (pol >= 1) parame(i,2) = 2.8270129502
    IF (pol >= 1) parame(i,3) = 176.3760515
    IF (pol >= 1) nhb_typ(i)     = 2
    IF (pol >= 1) nhb_no(i,1)    = 1
    IF (pol >= 1) nhb_no(i,2)    = 1
    IF (pol >= 1) eps_hb(i,i,1,2)= 2332.5845803
    IF (pol >= 1) eps_hb(i,i,2,1)= 2332.5845803
    IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
    IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
    IF (pol >= 1) kap_hb(i,i)    = 8.9248658086E-2
    IF (pol >= 1) parame(i,6) = 1.7
    IF (pol >= 1) lli(i)      = 1.75*parame(i,2)
    IF (pol >= 1) phi_criti(i)= 23.43
    IF (pol >= 1) chap(i)     = 0.304
    IF (pol == 2) mm(i)       = 32.042       ! PCIP-SAFT
    IF (pol == 2) parame(i,1) = 2.0693
    IF (pol == 2) parame(i,2) = 2.9547
    IF (pol == 2) parame(i,3) = 174.51
    IF (pol == 2) nhb_typ(i)     = 2
    IF (pol == 2) nhb_no(i,1)    = 1
    IF (pol == 2) nhb_no(i,2)    = 1
    IF (pol == 2) eps_hb(i,i,1,2)= 2418.5
    IF (pol == 2) eps_hb(i,i,2,1)= 2418.5
    IF (pol == 2) eps_hb(i,i,1,1)= 0.0
    IF (pol == 2) eps_hb(i,i,2,2)= 0.0
    IF (pol == 2) kap_hb(i,i)    = 0.06319
    IF (pol == 2) parame(i,6) = 1.7
    IF (pol == 2) parame(i,11)= 3.29
    ! mm(i)       =         32.0420000000000     ! PCP-SAFT with adjusted QQ
    ! parame(i,1) = mm(i)*  6.241807629559099E-002
    ! ! parame(i,1) =    2.00000000066333
    ! parame(i,2) =         2.97610169698593
    ! parame(i,3) =         163.268505098639
    ! nhb_typ(i)     =              2
    ! nhb_no(i,1)    =              1
    ! nhb_no(i,2)    =              1
    ! eps_hb(i,i,1,2)=      2449.55621933612
    ! eps_hb(i,i,2,1)=      2449.55621933612
    ! eps_hb(i,i,1,1)=  0.0
    ! eps_hb(i,i,2,2)=  0.0
    ! kap_hb(i,i)    =     8.431015160393653E-002
    ! parame(i,6) =         1.72000000000000
    ! parame(i,7) =         1.59810028824523
  ELSE IF (compna(i) == 'ethanol') THEN
    mm(i)       = 46.069
    parame(i,1) = mm(i)*  .0517195521
    parame(i,2) = 3.17705595
    parame(i,3) = 198.236542
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2653.38367
    eps_hb(i,i,2,1)= 2653.38367
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0323840159
    IF (pol >= 1) mm(i)       = 46.0690000000000
    IF (pol >= 1) parame(i,1) = mm(i)*  4.753626908781145E-002  !  =2.18994838060639
    IF (pol >= 1) parame(i,2) = 3.30120000000000
    IF (pol >= 1) parame(i,3) = 209.824555801706
    IF (pol >= 1) nhb_typ(i)     = 2
    IF (pol >= 1) nhb_no(i,1)    = 1
    IF (pol >= 1) nhb_no(i,2)    = 1
    IF (pol >= 1) eps_hb(i,i,1,2)= 2584.53116785767
    IF (pol >= 1) eps_hb(i,i,2,1)= 2584.53116785767
    IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
    IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
    IF (pol >= 1) kap_hb(i,i)    = 2.349382956935725E-002
    IF (pol >= 1) parame(i,6) = 1.69000000000000
    ! mm(i)       = 46.0690000000000
    ! parame(i,1) = mm(i)*  5.117957752785066E-002   !  =2.357791957
    ! parame(i,2) = 3.24027031244304
    ! parame(i,3) = 175.657110615456
    ! nhb_typ(i)     = 2
    ! nhb_no(i,1)    = 1
    ! nhb_no(i,2)    = 1
    ! eps_hb(i,i,1,2)= 2273.62670516146
    ! eps_hb(i,i,2,1)= 2273.62670516146
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 7.030279197039477E-002
    ! parame(i,6) = 1.69000000000000
    ! parame(i,7) = 3.63701294195013
    IF (pol == 2) mm(i)       = 46.0690000000000
    IF (pol == 2) parame(i,1) = mm(i)*  4.733436280008321E-002  !  =2.18064676
    IF (pol == 2) parame(i,2) = 3.31260000000000
    IF (pol == 2) parame(i,3) = 207.594119926613
    IF (pol == 2) nhb_typ(i)  = 2
    IF (pol == 2) nhb_no(i,1) = 1
    IF (pol == 2) nhb_no(i,2) = 1
    IF (pol == 2) eps_hb(i,i,1,2)= 2589.68311382661
    IF (pol == 2) eps_hb(i,i,2,1)= 2589.68311382661
    IF (pol == 2) eps_hb(i,i,1,1)= 0.0
    IF (pol == 2) eps_hb(i,i,2,2)= 0.0
    IF (pol == 2) kap_hb(i,i) = 2.132561218631547E-002
    IF (pol == 2) parame(i,6) = 1.69000000000000
    IF (pol == 2) parame(i,7) = 0.0
    IF (pol == 2) parame(i,11)= 5.11000000000000
  ELSE IF (compna(i) == '1-propanol') THEN
    mm(i)       = 60.096
    parame(i,1) = mm(i)*  .0499154461
    parame(i,2) = 3.25221234
    parame(i,3) = 233.396705
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2276.77867
    eps_hb(i,i,2,1)= 2276.77867
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0152683094
  ELSE IF (compna(i) == '1-butanol') THEN
    mm(i)       = 74.123
    parame(i,1) = mm(i)*  .0341065046
    parame(i,2) = 3.72361538
    parame(i,3) = 269.798048
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2661.37119
    eps_hb(i,i,2,1)= 2661.37119
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .00489087833
    mm(i)       = 74.1230000000000
    parame(i,1) = mm(i)*  3.329202420547412E-002  !  =2.46770471018236
    parame(i,2) = 3.76179376417092
    parame(i,3) = 270.237284242002
    nhb_typ(i)     = 2
    nhb_no(i,1)    = 1
    nhb_no(i,2)    = 1
    eps_hb(i,i,1,2)= 2669.28754983370
    eps_hb(i,i,2,1)= 2669.28754983370
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 4.855584122733399E-003
    parame(i,6) = 1.66000000000000
  ELSE IF (compna(i) == '1-pentanol') THEN
    mm(i)       = 88.15               ! PC-SAFT
    parame(i,1) = mm(i)*  .041134139
    parame(i,2) = 3.45079143
    parame(i,3) = 247.278748
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2252.09237
    eps_hb(i,i,2,1)= 2252.09237
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0103189939
    IF (pol >= 1) mm(i)       = 88.1500000000000        ! PCP-SAFT
    IF (pol >= 1) parame(i,1) = mm(i)*  4.138903382168521E-002  !  =3.64844333138155
    IF (pol >= 1) parame(i,2) = 3.44250118689142
    IF (pol >= 1) parame(i,3) = 246.078034174947
    IF (pol >= 1) nhb_typ(i)     = 2
    IF (pol >= 1) nhb_no(i,1)    = 1
    IF (pol >= 1) nhb_no(i,2)    = 1
    IF (pol >= 1) eps_hb(i,i,1,2)= 2236.72830142446
    IF (pol >= 1) eps_hb(i,i,2,1)= 2236.72830142446
    IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
    IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
    IF (pol >= 1) kap_hb(i,i)    = 1.040067895187016E-002
    IF (pol >= 1) parame(i,6) = 1.70000000000000
    IF (pol == 2) mm(i)       = 88.1500000000000        ! PCIP-SAFT
    IF (pol == 2) parame(i,1) = mm(i)*  4.161521814399406E-002  !  =3.66838147939308
    IF (pol == 2) parame(i,2) = 3.43496654431777
    IF (pol == 2) parame(i,3) = 244.177313808431
    IF (pol == 2) nhb_typ(i)     = 2
    IF (pol == 2) nhb_no(i,1)    = 1
    IF (pol == 2) nhb_no(i,2)    = 1
    IF (pol == 2) eps_hb(i,i,1,2)= 2241.27880639096
    IF (pol == 2) eps_hb(i,i,2,1)= 2241.27880639096
    IF (pol == 2) eps_hb(i,i,1,1)= 0.0
    IF (pol == 2) eps_hb(i,i,2,2)= 0.0
    IF (pol == 2) kap_hb(i,i)    = 1.049516309928397E-002
    IF (pol == 2) parame(i,6) = 1.70000000000000
    IF (pol == 2) parame(i,11)= 10.8000000000000
  ELSE IF (compna(i) == '1-octanol') THEN
    mm(i)       = 130.23
    parame(i,1) = mm(i)*  .0334446084
    parame(i,2) = 3.714535
    parame(i,3) = 262.740637
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2754.77272
    eps_hb(i,i,2,1)= 2754.77272
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .00219656803
  ELSE IF (compna(i) == '1-nonanol') THEN
    mm(i)       = 144.257
    parame(i,1) = mm(i)*  .0324692669
    parame(i,2) = 3.72924286
    parame(i,3) = 263.636673
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2941.9231
    eps_hb(i,i,2,1)= 2941.9231
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .00142696883
  ELSE IF (compna(i) == '2-propanol') THEN
    mm(i)       = 60.096
    parame(i,1) = mm(i)*  .0514663133
    parame(i,2) = 3.20845858
    parame(i,3) = 208.420809
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2253.91064
    eps_hb(i,i,2,1)= 2253.91064
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0246746934
  ELSE IF (compna(i) == '2-methyl-2-butanol') THEN
    mm(i)       = 88.15
    parame(i,1) = mm(i)*  .0289135026
    parame(i,2) = 3.90526707
    parame(i,3) = 266.011828
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2618.80124
    eps_hb(i,i,2,1)= 2618.80124
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .00186263367
  ELSE IF (compna(i) == 'acetic-acid') THEN
    mm(i)       = 60.053
    parame(i,1) = mm(i)*  .0227076949
    parame(i,2) = 3.79651163
    parame(i,3) = 199.225066
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 3092.40109
    eps_hb(i,i,2,1)= 3092.40109
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0870093874
    
    
    mm(i)       = 60.053
    parame(i,1) = mm(i)*  .0181797646
    parame(i,2) = 4.13711044
    parame(i,3) = 207.552969
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 3198.84362
    eps_hb(i,i,2,1)= 3198.84362
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0586552968
    
! mit gesetztem Dipol-Moment
    mm(i)       = 60.0530000000000
    parame(i,1) = mm(i)*  1.736420143637533E-002
    parame(i,2) = 4.25220708070687
    parame(i,3) = 190.957247854820
    parame(i,6) = 3.50000000000000
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 3096.36190957945
    eps_hb(i,i,2,1)= 3096.36190957945
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 6.154307094782551E-002
    
  ELSE IF (compna(i) == 'propionic-acid') THEN
    mm(i)       = 74.0800000000000
    parame(i,1) = mm(i)*  2.359519915877884E-002
    parame(i,2) = 3.99339224153844
    parame(i,3) = 295.947729838532
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2668.97826430874
    eps_hb(i,i,2,1)= 2668.97826430874
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 3.660242292423115E-002
  ELSE IF (compna(i) == 'acrylic-acid') THEN
    mm(i)       = 72.0636
    parame(i,1) = mm(i)*  .0430585424
    parame(i,2) = 3.0545415
    parame(i,3) = 164.115604
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 3065.40667
    eps_hb(i,i,2,1)= 3065.40667
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .336261669
  ELSE IF (compna(i) == 'caproic-acid') THEN
    mm(i)       =         116.16
    parame(i,1) = 5.87151
    parame(i,2) =         3.0694697
    parame(i,3) =         241.4569
    nhb_typ(i)     =      1
    eps_hb(i,i,1,1)=      2871.37
    kap_hb(i,i)    =      3.411613D-3
  ELSE IF (compna(i) == 'aniline') THEN
    mm(i)       = 93.13
    parame(i,1) = mm(i)*  .0285695992
    parame(i,2) = 3.70214085
    parame(i,3) = 335.471062
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 1351.64234
    eps_hb(i,i,2,1)= 1351.64234
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0748830615
    
    mm(i)       = 93.1300000000000
    parame(i,1) = mm(i)*  2.834372610192228E-002  !  =2.63965121187202
    parame(i,2) = 3.71326867619433
    parame(i,3) = 332.253796842009
    nhb_typ(i)     = 2
    nhb_no(i,1)    = 1
    nhb_no(i,2)    = 1
    eps_hb(i,i,1,2)= 1392.14266886674
    eps_hb(i,i,2,1)= 1392.14266886674
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 7.424612087328866E-002
    parame(i,6) = 1.55000000000000
    IF (pol == 2) parame(i,11)= 12.1000000000000
  ELSE IF (compna(i) == 'HF') THEN
    ! mm(i)       = 20.006       !  PC-SAFT
    ! parame(i,1) = 0.87731
    ! parame(i,2) = 3.0006
    ! parame(i,3) = 112.24
    ! nhb_typ(i)     = 2
    ! nhb_no(i,1)    = 1
    ! nhb_no(i,2)    = 1
    ! eps_hb(i,i,1,2)= 2208.22
    ! eps_hb(i,i,2,1)= 2208.22
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 0.71265
    mm(i)       = 20.0060000000000    !  PCP-SAFT
    parame(i,1) = 1.00030000000000
    parame(i,2) = 3.17603622195029
    parame(i,3) = 331.133373208282
    nhb_typ(i)     = 2
    nhb_no(i,1)    = 1
    nhb_no(i,2)    = 1
    eps_hb(i,i,1,2)= 348.251433080979
    eps_hb(i,i,2,1)= 348.251433080979
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 2.868887975449893E-002
    parame(i,6) = 1.82600000000000
  ELSE IF (compna(i) == 'HCl') THEN
    ! mm(i)       = 36.4610000000000
    ! parame(i,1) = mm(i)*  3.922046741026943E-002
    ! parame(i,2) = 3.08731180727493
    ! parame(i,3) = 203.898845304388
    ! nhb_typ(i)     = 2
    ! nhb_no(i,1)    = 1
    ! nhb_no(i,2)    = 1
    ! eps_hb(i,i,1,2)= 245.462773177367
    ! eps_hb(i,i,2,1)= 245.462773177367
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 0.256454187330899
    mm(i)       = 36.461        ! PCIP-SAFT
    parame(i,1) = 1.6335
    parame(i,2) = 2.9066
    parame(i,3) = 190.17
    parame(i,6) = 1.1086
    IF (pol == 2) parame(i,11)= 2.63
  ELSE IF (compna(i) == 'gen') THEN
    mm(i)       = 302.0
    parame(i,1) = 8.7563
    parame(i,2) = 3.604243
    parame(i,3) = 255.6434
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 0.0
    eps_hb(i,i,2,1)= 0.0
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = 0.02
  ELSE IF (compna(i) == 'h2o') THEN
    mm(i)       = 18.015
    parame(i,1) = mm(i)*  .05915
    parame(i,2) = 3.00068335
    parame(i,3) = 366.512135
    nhb_typ(i)  = 2
    nhb_no(i,1) = 1            ! no. of sites of type 1
    nhb_no(i,2) = 1            ! no. of sites of type 2
    eps_hb(i,i,1,2)= 2500.6706
    eps_hb(i,i,2,1)= 2500.6706
    eps_hb(i,i,1,1)= 0.0
    eps_hb(i,i,2,2)= 0.0
    kap_hb(i,i)    = .0348679836
    
    ! mm(i)       = 18.015
    ! parame(i,1) = 1.706
    ! parame(i,2) = 2.429
    ! parame(i,3) = 242.19
    ! nhb_typ(i)  = 2
    ! nhb_no(i,1) = 1            ! no. of sites of type 1
    ! nhb_no(i,2) = 1            ! no. of sites of type 2
    ! eps_hb(i,i,1,2)= 2644.2
    ! eps_hb(i,i,2,1)= 2644.2
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 0.153
    
    ! mm(i)       = 18.015
    ! parame(i,1) = mm(i)*  .0588185709
    ! parame(i,2) = 3.02483704
    ! parame(i,3) = 382.086672
    !c! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    ! nhb_typ(i)  = 2
    ! nhb_no(i,1) = 1            ! no. of sites of type 1
    ! nhb_no(i,2) = 2            ! no. of sites of type 2
    !c! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    ! eps_hb(i,i,1,2)= 2442.49782
    ! eps_hb(i,i,2,1)= 2442.49782
    ! eps_hb(i,i,1,1)=0.0
    ! eps_hb(i,i,2,2)=0.0
    ! kap_hb(i,i)    = .0303754635
    
    ! mit gefittetem Dipol-Moment - Haarlem-night
    ! mm(i)       = 18.015
    ! parame(i,1) = mm(i)*  7.0037160952278E-2
    ! parame(i,2) = 2.79276650240763
    ! parame(i,3) = 270.970053834558
    ! nhb_typ(i)  = 2
    ! nhb_no(i,1) = 1            ! no. of sites of type 1
    ! nhb_no(i,2) = 1            ! no. of sites of type 2
    ! eps_hb(i,i,1,2)= 1427.8287
    ! eps_hb(i,i,2,1)= 1427.8287
    ! eps_hb(i,i,1,1)=0.0
    ! eps_hb(i,i,2,2)=0.0
    ! kap_hb(i,i)    = 4.335167238E-2
    ! parame(i,6) = 3.968686856378
    
    IF (pol >= 1) mm(i)       = 18.015       !  PCP-SAFT
    IF (pol >= 1) parame(i,1) = 0.922688825223317
    IF (pol >= 1) parame(i,2) = 3.17562052023944
    IF (pol >= 1) parame(i,3) = 388.462197714696
    IF (pol >= 1) nhb_typ(i)     = 2
    IF (pol >= 1) nhb_no(i,1)    = 1
    IF (pol >= 1) nhb_no(i,2)    = 1
    IF (pol >= 1) eps_hb(i,i,1,2)= 2000.67247409031
    IF (pol >= 1) eps_hb(i,i,2,1)= 2000.67247409031
    IF (pol >= 1) eps_hb(i,i,1,1)= 0.0
    IF (pol >= 1) eps_hb(i,i,2,2)= 0.0
    IF (pol >= 1) kap_hb(i,i)    = 2.040614952751225E-003
    IF (pol >= 1) parame(i,6) = 1.85500000000000
    IF (pol >= 1) parame(i,7) = 2.00000000000000
    ! IF (pol.EQ.2) mm(i)       = 18.015       !  PCIP-SAFT - DQ with my=my_RPT
    ! IF (pol.EQ.2) parame(i,1) = 1.0
    ! IF (pol.EQ.2) parame(i,2) = 3.14540664928026
    ! IF (pol.EQ.2) parame(i,3) = 320.283823615925
    ! IF (pol.EQ.2) nhb_typ(i)     = 2
    ! IF (pol.EQ.2) nhb_no(i,1)    = 2
    ! IF (pol.EQ.2) nhb_no(i,2)    = 2
    ! IF (pol.EQ.2) eps_hb(i,i,1,2)= 1335.72887678032
    ! IF (pol.EQ.2) eps_hb(i,i,2,1)= 1335.72887678032
    ! IF (pol.EQ.2) eps_hb(i,i,1,1)= 0.0
    ! IF (pol.EQ.2) eps_hb(i,i,2,2)= 0.0
    ! IF (pol.EQ.2) kap_hb(i,i)    = 4.162619960844732E-003
    ! IF (pol.EQ.2) parame(i,6) = 1.85500000000000
    ! IF (pol.EQ.2) parame(i,7) = 2.00000000000000
    ! IF (pol.EQ.2) parame(i,11)= 1.45000000000000
    ! mm(i)       = 18.0150000000000
    ! parame(i,1) = 1.0
    ! parame(i,2) = 3.11505069470915
    ! parame(i,3) = 320.991387913502
    ! nhb_typ(i)     = 2
    ! nhb_no(i,1)    = 1
    ! nhb_no(i,2)    = 1
    ! eps_hb(i,i,1,2)= 2037.76329812542
    ! eps_hb(i,i,2,1)= 2037.76329812542
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 3.763148982832804E-003
    ! parame(i,6) = 1.85500000000000
    ! parame(i,7) = 2.00000000000000
    ! IF (pol.EQ.2) parame(i,11)= 1.45000000000000
    IF (pol == 2) mm(i)       = 18.015       !  PCIP-SAFT - DQ with my=my_0
    IF (pol == 2) parame(i,1) = 1.0
    IF (pol == 2) parame(i,2) = 3.11574491885322
    IF (pol == 2) parame(i,3) = 322.699984283499
    IF (pol == 2) nhb_typ(i)     = 2
    IF (pol == 2) nhb_no(i,1)    = 1
    IF (pol == 2) nhb_no(i,2)    = 1
    IF (pol == 2) eps_hb(i,i,1,2)= 2033.87777692450
    IF (pol == 2) eps_hb(i,i,2,1)= 2033.87777692450
    IF (pol == 2) eps_hb(i,i,1,1)= 0.0
    IF (pol == 2) eps_hb(i,i,2,2)= 0.0
    IF (pol == 2) kap_hb(i,i)    = 3.815764667176484E-003
    IF (pol == 2) parame(i,6) = 1.85500000000000
    IF (pol == 2) parame(i,7) = 2.00000000000000
    IF (pol == 2) parame(i,11)= 1.45000000000000
    ! mm(i)       = 18.015          ! Dortmund
    ! parame(i,1) = 0.11065254*mm(i)
    ! parame(i,2) = 2.36393615
    ! parame(i,3) = 300.288589
    ! nhb_typ(i)     = 2
    ! nhb_no(i,1)    = 1
    ! nhb_no(i,2)    = 1
    ! eps_hb(i,i,1,2)= 1193.45585
    ! eps_hb(i,i,2,1)= 1193.45585
    ! eps_hb(i,i,1,1)= 0.0
    ! eps_hb(i,i,2,2)= 0.0
    ! kap_hb(i,i)    = 0.091203519
    ! parame(i,6) = 1.8546
    ! parame(i,7) = 0.0
    ! parame(i,11)= 0.0
  ELSE IF (compna(i) == 'MBBA') THEN
    mm(i)     =         267.37
    parame(i,1) =        12.194
    parame(i,2) =         3.064
    parame(i,3) =       270.7
    e_lc(i,i)   =        13.7      !Hino & Prausnitz
    s_lc(i,i)   =         0.176    !Hino & Prausnitz
  ELSE IF (compna(i) == 'PCH5') THEN
    mm(i)     =         255.41
    parame(i,1) =        11.6
    parame(i,2) =         3.2
    parame(i,3) =       270.7
    ! E_LC(i,i)   =        16.7      !Hino & Prausnitz
    ! S_LC(i,i)   =         0.176    !Hino & Prausnitz
    e_lc(i,i)   =        8.95
    s_lc(i,i)   =         0.2
    
    ! mm(i)     =         255.41
    ! parame(i,1) =        11.6
    ! parame(i,2) =         3.2
    ! parame(i,3) =       290.7
    ! E_LC(i,i)   =        7.18
    ! S_LC(i,i)   =         0.2
    
  ELSE IF (compna(i) == 'Li') THEN
    mm(i)       = 6.9
    parame(i,1) = 1.0
    parame(i,2) = 1.4
    parame(i,3) = 96.83
    parame(i,10)= 1.0
!   the self-association is set to zero in routine F_EXPL for ions
    ! nhb_typ(i)  = 1
    ! nhb_no(i,1) = 3                      ! no. of sites of type 1
    ! eps_hb(i,i,1,1)= 2000.0
    ! kap_hb(i,i)    = 0.008
  ELSE IF (compna(i) == 'Na') THEN
    mm(i)       = 23.0
    parame(i,1) = 1.0
    parame(i,2) = 1.9
    parame(i,3) = 147.38
    parame(i,10)= 1.0
!   the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 3                       ! no. of sites of type 1
    eps_hb(i,i,1,1)= 8946.28257     ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.001648933
  ELSE IF (compna(i) == 'Ka') THEN
    mm(i)       = 39.1
    parame(i,1) = 1.0
    parame(i,2) = 2.66
    parame(i,3) = 221.44
    parame(i,10)= 1.0
    ! the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 3                      ! no. of sites of type 1
    eps_hb(i,i,1,1)= 3118.336216     ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'Cs') THEN
    mm(i)       = 132.9
    parame(i,1) = 1.0
    parame(i,2) = 3.38
    parame(i,3) = 523.28
    parame(i,10)= 1.0
    ! the self-association is set to zero in routine F_EXPL for ions
    ! nhb_typ(i)  = 1
    ! nhb_no(i,1) = 3                      ! no. of sites of type 1
    ! eps_hb(i,i,1,1)= 2000.0
    ! kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'Cl') THEN
    mm(i)       = 35.5
    parame(i,1) = 1.0
    parame(i,2) = 3.62
    parame(i,3) = 225.44
    parame(i,10)= -1.0
!   the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 4                      ! no. of sites of type 1
    eps_hb(i,i,1,1)= 6744.12509     ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.00155252
  ELSE IF (compna(i) == 'Br') THEN
    mm(i)       = 79.9
    parame(i,1) = 1.0
    parame(i,2) = 3.9
    parame(i,3) = 330.82
    parame(i,10)= -1.0
!   the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 4                      ! no. of sites of type 1
    eps_hb(i,i,1,1)= 4516.033227     ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'Io') THEN
    mm(i)       = 126.9
    parame(i,1) = 1.0
    parame(i,2) = 4.4
    parame(i,3) = 380.60
    parame(i,10)= -1.0
!   the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 4                      ! no. of sites of type 1
    eps_hb(i,i,1,1)= 1631.203342     ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'OH') THEN
    mm(i)       = 17.0
    parame(i,1) = 1.0
    parame(i,2) = 3.06
    parame(i,3) = 217.26
    parame(i,10)= -1.0
    ! the self-association is set to zero in routine F_EXPL for ions
    nhb_typ(i)  = 1
    nhb_no(i,1) = 4                      ! no. of sites of type 1
    eps_hb(i,i,1,1)= 14118.68089    ! 25C, 3 sites, dG-ref-1
    kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'NO3') THEN
    mm(i)       = 62.0
    parame(i,1) = 1.0
    parame(i,2) = 4.12
    parame(i,3) = 239.48
    parame(i,10)= -1.0
    ! the self-association is set to zero in routine F_EXPL for ions
    ! nhb_typ(i)  = 1
    ! nhb_no(i,1) = 4                      ! no. of sites of type 1
    ! eps_hb(i,i,1,1)= 2000.0
    ! kap_hb(i,i)    = 0.00200559
  ELSE IF (compna(i) == 'bf4') THEN
    mm(i)       = 86.8
    parame(i,1) = 1.0
    parame(i,2) = 4.51 ! *0.85
    parame(i,3) = 164.7
    parame(i,10)= -1.0
  ELSE IF (compna(i) == 'pf6') THEN
    mm(i)       = 145.0
    parame(i,1) = 1.0
    parame(i,2) = 5.06
    parame(i,3) = 224.9
    parame(i,10)= -1.0
  ELSE IF (compna(i) == 'emim') THEN
    mm(i)       = 111.16
    parame(i,1) = 3.11
    parame(i,2) = 4.0
    parame(i,3) = 250.0
    parame(i,10)= 1.0
  ELSE IF (compna(i) == 'bmim') THEN
    mm(i)       = 139.21
    ! parame(i,1) = 2.81
    ! parame(i,2) = 3.5
    parame(i,1) = 3.81
    parame(i,2) = 4.0
    parame(i,3) = 250.0
    parame(i,6) = 0.0
    parame(i,10)= 1.0
  ELSE IF (compna(i) == 'hmim') THEN
    mm(i)       = 167.27
    parame(i,1) = 4.53
    parame(i,2) = 4.0
    parame(i,3) = 250.0
    parame(i,10)= 1.0
  ELSE IF (compna(i) == 'omim') THEN
    mm(i)       = 195.32
    parame(i,1) = 5.30
    parame(i,2) = 4.0
    parame(i,3) = 250.0
    parame(i,10)= 1.0
  ELSE IF (compna(i) == 'sw') THEN
    parame(i,1) = 1.0
    parame(i,2) = 1.0
    parame(i,3) = 100.0
    parame(i,4) = 0.0
    parame(i,5) = 0.0
    mm(i)       = 1.0
    parame(i,6) = 0.1175015839*2.0
    ! use Temp. in Kelvin in the input-file. For dimensionless quantities
    ! (P*=P*sig**3/epsilon, T*=T*kBol/epsilon, rho*=rho*sig**3) calculate
    ! P* = P *1E5 * (1.e-10)^3 / (100*8.31441/6.022045E+23)
    ! T* = (T+273.15)/100
    ! for rho* go to utilities.f (subroutine SI_DENS) and write
    ! density(ph) = dense(ph)*6.0/PI
  ELSE IF (compna(i) == 'c8-sim') THEN
    mm(i)       =         114.231
    parame(i,1) = mm(i)*  4.095944E-2  !  =4.67883774717337
    parame(i,2) =         3.501769
    parame(i,3) =         163.8606
    ! mm(i)       =         114.231000000000
    ! parame(i,1) = mm(i)*  3.547001476437745E-002  !  =4.05177525654960
    ! parame(i,2) =         3.70988567055814
    ! parame(i,3) =         192.787548176343
  ELSE IF (compna(i) == 'argon_ge') THEN
    mm(i)       = 39.948
    parame(i,1) = mm(i)*0.030327
    parame(i,2) = 3.149910
    parame(i,3) = 100.188975
  ELSE IF (compna(i) == 'argon_ge2') THEN
    mm(i)       = 39.948
    parame(i,1) = mm(i)*0.030327
    parame(i,2) = 3.149910
    parame(i,3) = 0.8*100.188975
  ELSE
    WRITE (*,*) ' pure component parameters missing for ',compna(i)
    STOP
  END IF
  
  IF (pol == 2.AND.parame(i,11) == 0.0) THEN
    WRITE (*,*) ' polarizability missing for comp. ',compna(i)
    STOP
  END IF
  
  IF (nhb_typ(i) /= 0) THEN
    parame(i,12) = DBLE(nhb_typ(i))
    parame(i,13) = kap_hb(i,i)
    no = 0
    DO j=1,nhb_typ(i)
      DO k=1,nhb_typ(i)
        parame(i,(14+no))=eps_hb(i,i,j,k)
        no=no+1
      END DO
    END DO
    DO j=1,nhb_typ(i)
      parame(i,(14+no))=DBLE(nhb_no(i,j))
      no=no+1
    END DO
  ELSE
    DO k=12,25
      parame(i,k)=0.0
    END DO
  END IF
  
END DO


DO  i = 1,ncomp
  DO  j = 1,ncomp
    IF (compna(i) == 'ps'.AND.compna(j) == 'cyclohexane')THEN
      kij(i,j) = 0.0075
    ELSE IF(compna(i) == 'peva'.AND.compna(j) == 'ethylene')THEN
! --   0 Gew.% VA-------------
      ! kij(i,j) = 0.039
! -- 7.5 Gew.% VA-------------
      ! kij(i,j) = 0.0377325
      ! kij(i,j) = 0.0374867
! ---12.7 Gew.% VA------------
      ! kij(i,j) = 0.036854
      ! kij(i,j) = 0.0366508
! ---27.3 Gew.% VA------------
      ! kij(i,j) = 0.034386
      ! kij(i,j) = 0.0352375
! ---31.8 Gew.% VA------------
      kij(i,j) = 0.033626
      ! kij(i,j) = 0.0350795
! ---42.7 Gew.% VA------------
      ! kij(i,j) = 0.031784
      ! kij(i,j) = 0.035239
    ELSE IF(compna(i) == 'gen'.AND.compna(j) == 'h2o')THEN
      kij(i,j) = - 0.2
    ELSE IF(compna(i) == 'peva'.AND.compna(j) == 'vinylacetate')THEN
      kij(i,j) = 0.019
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'co2') THEN
      IF ( pol == 0 ) kij(i,j) = 0.195
      IF ( pol == 1 ) kij(i,j) = 0.06
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'acetone') THEN
      kij(i,j) = 0.021
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'hexane') THEN
      kij(i,j) = 0.012
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'pentane')THEN
      kij(i,j) = 0.005
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'methylcyclohexane') THEN
      kij(i,j) = 0.0073
    ELSE IF(compna(i) == 'ps'.AND.compna(j) == 'ethylbenzene')THEN
      kij(i,j) = 0.008
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.181
      kij(i,j) = 0.088
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'propane') THEN
      kij(i,j) = 0.0206
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'butane') THEN
      kij(i,j) = 0.01
    ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'argon') THEN
      kij(i,j) = 0.01
    ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'butane') THEN
      kij(i,j) = 0.026
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'pentane') THEN
      ! kij(i,j) = -0.0195
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'hexane') THEN
      ! kij(i,j) = 0.008
      kij(i,j) = 0.004
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'ethylene') THEN
      kij(i,j) = 0.0404
      ! kij(i,j) = 0.0423
      ! kij(i,j) = 0.044
    ELSE IF(compna(i) == 'pe'.AND.compna(j) == 'cyclohexane') THEN
      kij(i,j) = -0.1
    ELSE IF(compna(i) == 'ldpe'.AND.compna(j) == 'cyclopentane')THEN
      kij(i,j) = -0.016
    ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'propane') THEN
      kij(i,j) = 0.0242
    ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'pentane') THEN
      kij(i,j) = 0.0137583176
    ELSE IF(compna(i) == 'pp'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.1767         ! without quadrupol-term
      kij(i,j) = 0.063         ! with quadrupol-term
    ELSE IF(compna(i) == 'pba'.AND.compna(j) == 'ethylene') THEN
      kij(i,j) = 0.026
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'ethane') THEN
      kij(i,j) = 0.017
    ELSE IF(compna(i) == 'n2'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = -0.04
    ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.051875  ! PC-SAFT
      IF (pol == 1) kij(i,j) = -0.0353125  ! PCP-SAFT
    ELSE IF(compna(i) == 'methane'.AND.compna(j) == 'co') THEN
      ! IF (pol == 1) kij(i,j) = -0.003  ! PCP-SAFT 
      IF (pol == 1) kij(i,j) = 0.018  ! PCP-SAFT 
    ELSE IF(compna(i) == 'ethane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.095
      kij(i,j) = 0.021
      ! kij(i,j) = 0.024
    ELSE IF(compna(i) == 'propane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.042
    ELSE IF(compna(i) == 'argon_ge'.AND.compna(j) == 'argon_ge2') THEN
      read (*,*) kij(i,j)
    ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.115
      ! kij(i,j) = 0.048
      kij(i,j) = 0.036
    ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.143         ! without quadrupol-term
      kij(i,j) = 0.0            ! with quadrupol-term
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.125         ! without quadrupol-term
      kij(i,j) = 0.0495
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.11         ! without quadrupol-term
      ! kij(i,j) = 0.05
      ! kij(i,j) = 0.039         ! with quadrupol-term
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.128         ! without quadrupol-term
      kij(i,j) = 0.053         ! with quadrupol-term
    ELSE IF(compna(i) == 'dodecane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.12         ! without quadrupol-term
      kij(i,j) = 0.0508         ! with quadrupol-term
    ELSE IF(compna(i) == 'benzene'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.087968750000    ! without quadrupol-term
      ! kij(i,j) =  0.008203125   ! only co2 quadrupol
      kij(i,j) = 0.042     ! both quadrupol
      ! kij(i,j) = 0.003     ! both quadrupol
    ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'co2') THEN
      ! kij(i,j) = 0.110784912  ! without quadrupol-term
      kij(i,j) = 0.0305        ! with quadrupol-term
    ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'co2') THEN
      kij(i,j) = 0.13
      lij(i,j) = - 0.00
      ! kij(i,j) = 0.045
    ELSE IF(compna(i) == 'chloromethane'.AND.compna(j) == 'co2')THEN
      ! kij(i,j) = 0.04  ! PC-SAFT
      kij(i,j) = 0.025  ! PCP-SAFT
      ! kij(i,j) = 0.083  ! PCIP-SAFT
    ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'n2')THEN
      kij(i,j) = 0.035211  ! PCP-SAFT
      lij(i,j) = + 0.013225  ! PCP-SAFT
      kij(i,j) = 0.023  ! PCP-SAFT
      lij(i,j) = + 0.013225  ! PCP-SAFT
      !kij(i,j) = 1.722238535635467E-002  ! PCP-SAFT
      !lij(i,j) = 2.815974678394451E-003  ! PCP-SAFT
      !kij(i,j) = 1.931522058164026E-002  ! PCP-SAFT
      !lij(i,j) = 0.0  ! PCP-SAFT
      !kij(i,j) =  1.641053794134795E-002  ! PCP-SAFT
      !lij(i,j) = -5.850421759950764E-003  ! PCP-SAFT
      if ( num == 0 ) write (*,*) 'calculation with lij only possible with num=1'
      if ( num == 0 ) stop
    ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'co2')THEN
      kij(i,j) = 0.015  ! PC-SAFT
      IF (pol == 1) kij(i,j) = -0.02  ! PCP-SAFT
      IF (pol == 2) kij(i,j) = -0.005  ! PCIP-SAFT  where DQ with my=my_vacuum
      ! IF (pol.EQ.2) kij(i,j) = 0.0  ! PCIP-SAFT  where DQ with my=my_RPT
    ELSE IF(compna(i) == 'methanol'.AND.compna(j) == 'co2')THEN
      ! kij(i,j) =  0.0288  ! PC-SAFT
      ! kij(i,j) = - 0.035  ! PCP-SAFT for co2 and PC-SAFT methanol
      ! kij(i,j) = - 0.035  ! PCP-SAFT
      ! lij(i,j) =  0.3  ! PCP-SAFT
    ELSE IF(compna(i) == 'dimethyl-ether'.AND.compna(j) == 'co2')THEN
      kij(i,j) = 0.00896894  ! PC-SAFT
      ! kij(i,j) = - 0.015  ! PCP-SAFT
    ELSE IF(compna(i) == 'dimethyl-ether'.AND.compna(j) == 'h2o')THEN
      ! kij(i,j) = -0.134  ! PC-SAFT
    ELSE IF(compna(i) == 'dichloromethane'.AND.compna(j) == 'co2')THEN
      ! kij(i,j) = 0.06881725  ! PC-SAFT
      ! kij(i,j) = 0.05839145  ! PCP-SAFT
      kij(i,j) = -0.00944346  ! PCP-SAFT co2 dichloromethane PC-SAFT
      ! kij(i,j) = 0.06  ! PCIP-SAFT
    ELSE IF(compna(i) == 'h2s'.AND.compna(j) == 'methane')THEN
      ! kij(i,j) =  0.0414 ! PC-SAFT
      kij(i,j) =  0.0152 ! PCP-SAFT Dipole momnet (d with Q)
    ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'h2s')THEN
      kij(i,j) = -0.002  ! PCP-SAFT
    ELSE IF(compna(i) == 'methanol'.AND.compna(j) == 'h2s')THEN
      kij(i,j) = 0.0  ! PCP-SAFT
      lij(i,j) = 0.0  ! PCP-SAFT
    ELSE IF(compna(i) == 'co2'.AND.compna(j) == 'hydrogen') THEN
      ! lij(i,j) = -0.08 !!!!! Lij not kij !!!!
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'n2') THEN
      kij(i,j) = 0.11
    ELSE IF(compna(i) == 'propane'.AND.compna(j) == 'n2') THEN
      kij(i,j) = 0.0251171875
    ELSE IF(compna(i) == 'co2'.AND.compna(j) == 'hexadecane') THEN
      ! kij(i,j) = 0.1194   ! PC-SAFT ohne QQ
      kij(i,j) = 0.0588
    ELSE IF(compna(i) == 'ethane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.065  ! no DD
      kij(i,j) = 0.038   ! DD non-polarizable
      ! kij(i,j) = 0.025   ! DD polarizable
    ELSE IF(compna(i) == 'butane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.065  ! no DD
      kij(i,j) = 0.037  ! DD non-polarizable
      ! kij(i,j) = 0.025  ! DD polarizable
    ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.072  ! no DD
      ! kij(i,j) = 0.041  ! DD non-polarizable
      kij(i,j) = 0.039  ! DD polarizable
      ! kij(i,j) = 0.035  ! DD polarizable
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.063
      kij(i,j) = 0.038  ! PCP-SAFT
      ! kij(i,j) = 0.036  ! PCIP-SAFT
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'acetone') THEN
      kij(i,j) = 0.035  ! PCP-SAFT
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.059  ! no DD
      ! kij(i,j) =  0.03281250  ! DD non-polarizable
      kij(i,j) = 0.028  ! DD polarizable
    ELSE IF(compna(i) == 'hexadecane'.AND.compna(j) == 'acetone') THEN
      ! kij(i,j) = 0.07  ! no DD
      ! kij(i,j) =  0.0415  ! DD non-polarizable
      kij(i,j) = 0.035  ! DD polarizable
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'butanone')THEN
      kij(i,j) = 0.027   ! PCP-SAFT
      ! kij(i,j) = 0.033   ! PCP-SAFT with lij
      ! lij(i,j) = 0.13    ! PCP-SAFT
      ! kij(i,j) = 0.042    ! PC-SAFT
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'butanone')THEN
      kij(i,j) = 0.042  ! no DD
      ! kij(i,j) = 0.027  ! DD non-polarizable
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == '2-pentanone')THEN
      kij(i,j) = 0.041  ! no DD
      ! kij(i,j) = 0.031  ! DD non-polarizable
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == '3-pentanone')THEN
      kij(i,j) = 0.0
    ELSE IF(compna(i) == 'pentane'.AND.compna(j) == 'propanal')THEN
      ! kij(i,j) = 0.055  ! no DD
      kij(i,j) = 0.027  ! DD non-polarizable
      ! kij(i,j) = 0.026  ! DD polarizable 22
    ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'propanal')THEN
      ! kij(i,j) = 0.06  ! no DD
      kij(i,j) = 0.036  ! DD non-polarizable
      kij(i,j) = 0.035  ! DD polarizable 22
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'butanal')THEN
      kij(i,j) = 0.041  ! no DD
      ! kij(i,j) = 0.025  ! DD non-polarizable
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'thf')THEN
      kij(i,j) = 0.012  ! PCP-SAFT
    ELSE IF(compna(i) == 'octane'.AND.compna(j) == 'thf')THEN
      kij(i,j) = 0.012  ! PCP-SAFT
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'thf')THEN
      kij(i,j) = 0.012  ! PCP-SAFT
    ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'dmso')THEN
      ! kij(i,j) = 0.025  ! no DD
      kij(i,j) = - 0.0105  ! DD non-polarizable
      ! kij(i,j) = - 0.019  ! DD polarizable
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == 'acrylonitrile')THEN
      kij(i,j) = - 0.05  ! DD polarizable
    ELSE IF(compna(i) == 'heptane' .AND. compna(j) == 'butyronitrile')THEN
      kij(i,j) = - 0.002  ! DD polarizable 11
      kij(i,j) = 0.002  ! DD polarizable 22
    ELSE IF(compna(i) == '1-butene'.AND.compna(j) == 'dmf')THEN
      ! kij(i,j) = 0.04  ! no DD
      ! kij(i,j) = 0.004  ! DD non-polarizable
      kij(i,j) = 0.005  ! DD polarizable 22
    ELSE IF(compna(i) == 'cyclohexane'.AND.compna(j) == 'dmf')THEN
      kij(i,j) = 0.0135  ! DD polarizable 11
      kij(i,j) = 0.022  ! DD polarizable 22
    ELSE IF(compna(i) == 'ethylene'.AND.compna(j) == 'dmf')THEN
      ! kij(i,j) = - 0.0215  ! DD polarizable 11
      kij(i,j) = - 0.01  ! DD polarizable 22
    ELSE IF(compna(i) == 'nbutyl-ethanoate'.AND.compna(j) == 'dmf')THEN
      ! kij(i,j) = 0.016  ! no DD
      ! kij(i,j) = -0.01  ! DD non-polarizable
      kij(i,j) = - 0.015  ! DD polarizable 22
    ELSE IF(compna(i) == 'methylacetate' .AND. compna(j) == 'cyclohexane')THEN
      kij(i,j) = 0.066  ! PC-SAFT
      ! kij(i,j) = 0.061  ! PCP-SAFT
      ! kij(i,j) = 0.0625  ! PCIP-SAFT
    ELSE IF(compna(i) == 'methylacetate'.AND.compna(j) == 'decane')THEN
      kij(i,j) = 0.0625  ! PCIP-SAFT
    ELSE IF(compna(i) == 'methylacetate' .AND. compna(j) == 'methanol')THEN
      ! kij(i,j) = -0.07  ! PCIP-SAFT
    ELSE IF(compna(i) == 'pentane' .AND. compna(j) == 'propionitrile')THEN
      kij(i,j) =  0.0498
      IF (pol >= 1) kij(i,j) = -0.01
      IF (pol >= 2) kij(i,j) = -0.027
    ELSE IF(compna(i) == 'hexane' .AND. compna(j) == 'propionitrile')THEN
      kij(i,j) = 0.05
      IF (pol >= 1) kij(i,j) = 0.0
      IF (pol >= 2) kij(i,j) = -0.03
    ELSE IF(compna(i) == 'octane' .AND. compna(j) == 'propionitrile')THEN
      kij(i,j) = 0.0  ! DD polarizable 22
    ELSE IF(compna(i) == 'cyclohexane' .AND. compna(j) == 'nitromethane')THEN
      kij(i,j) =  0.14  ! no DD
      ! kij(i,j) =   0.07  ! DD non-polarizable
      ! kij(i,j) =   0.055  ! DD polarizable 22
    ELSE IF(compna(i) == 'cyclohexane' .AND. compna(j) == 'nitroethane')THEN
      ! kij(i,j) =  0.06  ! no DD
      kij(i,j) =   0.03  ! DD polarizable 22
    ELSE IF(compna(i) == 'acetone' .AND. compna(j) == 'nitromethane')THEN
      ! kij(i,j) = - 0.017  ! no DD
      kij(i,j) = - 0.021  ! DD non-polarizable
      ! kij(i,j) = - 0.023  ! DD polarizable 22
    ELSE IF(compna(i) == 'acetone' .AND. compna(j) == 'h2o')THEN
      kij(i,j) = - 0.2    ! PCP-SAFT (no cross-association)
    ELSE IF(compna(i) == 'methylcyclohexane' .AND. compna(j) == 'acetonitrile')THEN
      ! kij(i,j) =   0.09  ! no DD
      ! kij(i,j) =   0.033  ! DD non-polarizable
      ! kij(i,j) =   0.025  ! DD polarizable 22
      kij(i,j) =   0.04  ! DD polarizable 22 und my angepasst
    ELSE IF(compna(i) == 'ethylacetate' .AND. compna(j) == 'acetonitrile')THEN
      kij(i,j) =   0.007  ! no DD
      ! kij(i,j) =   -0.045  ! DD polarizable 22
    ELSE IF(compna(i) == 'dimethyl-ether' .AND. compna(j) == 'propane')THEN
      ! kij(i,j) = 0.03  ! no DD
      kij(i,j) =  0.0225  ! DD non-polarizable
    ELSE IF(compna(i) == 'benzene' .AND. compna(j) == 'pentane')THEN
      ! kij(i,j) =  0.012968750   ! ohne QQ
      ! kij(i,j) =  0.004921875   ! mit QQ=5.6D (gefittet)
      ! kij(i,j) = -0.006406250   ! mit QQ=7.45D (Literatur)
    ELSE IF(compna(i) == 'benzene' .AND. compna(j) == 'heptane')THEN
      kij(i,j) = 0.01328125
      ! kij(i,j) = 0.0038
    ELSE IF(compna(i) == 'benzene' .AND. compna(j) == '1-hexene')THEN
      kij(i,j) = 0.0067
    ELSE IF(compna(i) == 'ethylene' .AND. compna(j) == 'co2') THEN
      kij(i,j) = 0.04
      kij(i,j) = -0.029
    ELSE IF(compna(i) == 'ethylene' .AND. compna(j) == 'vinylacetate') THEN
      kij(i,j) = - 0.013847
    ELSE IF(compna(i) == 'triacontane' .AND. compna(j) == 'ethylene') THEN
      kij(i,j) =  0.028
    ELSE IF(compna(i) == 'triacontane' .AND. compna(j) == 'methane')THEN
      ! kij(i,j) =  0.061953125   ! polyethylene parameters
      kij(i,j) =  0.039609375   ! param. by extrapolation of n-alkanes
    ELSE IF(compna(i) == 'tetracontane' .AND. compna(j) == 'methane')THEN
      ! kij(i,j) =  0.06515625    ! polyethylene parameters
      kij(i,j) =  0.04453125    ! param. by extrapolation of n-alkanes
      ! --- kij and lij adjusted -------
      ! kij(i,j) =  0.045786119   ! param. by extrapolation of n-alkanes
      ! lij(i,j) =  +8.53466437d-4   ! param. by extrapolation of n-alkanes
    ELSE IF(compna(i) == 'eicosane' .AND. compna(j) == 'methane')THEN
      kij(i,j) =  0.0360134457445
    ELSE IF(compna(i) == 'tetracontane' .AND. compna(j) == 'methane')THEN
      kij(i,j) =  0.0360134457445  ! assumed equal to eicosane-C1
    ELSE IF(compna(i) == 'chlorobenzene' .AND. compna(j) == 'cyclohexane') THEN
      kij(i,j) = 0.013
    ELSE IF(compna(i) == 'chloroethane' .AND. compna(j) == 'butane')THEN
      kij(i,j) = 0.025
    ELSE IF(compna(i) == 'tetrachloromethane' .AND. compna(j) == '2-methylpentane') THEN
      kij(i,j) = 0.0070105
    ELSE IF(compna(i) == 'tetrachloromethane' .AND. compna(j) == 'hexane') THEN
      kij(i,j) = 0.004
    ELSE IF(compna(i) == 'hydrogen' .AND. compna(j) == 'hexane') THEN
      kij(i,j) = 0.1501
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'co2') THEN
      ! kij(i,j) = -0.08
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'propane') THEN
      kij(i,j) = - 0.07
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'ethane') THEN
      kij(i,j) = 0.0
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'butane') THEN
      kij(i,j) = 0.028
      kij(i,j) = 0.016
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == 'cyclohexane')THEN
      kij(i,j) = 0.037
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == '2-methylpentane') THEN
      kij(i,j) = 0.028
    ELSE IF(compna(i) == 'ethanol' .AND. compna(j) == '1-octanol')THEN
      kij(i,j) = 0.06
    ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'cyclohexane') THEN
      kij(i,j) = 0.0508
      ! kij(i,j) = 0.03
    ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'heptane')THEN
      kij(i,j) = 0.034
    ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'decane')THEN
      ! kij(i,j) = 0.042   !  PC-SAFT
      ! kij(i,j) = 0.011   !  PCP-SAFT
      kij(i,j) = 0.000    !  PCIP-SAFT
    ELSE IF(compna(i) == 'methanol' .AND. compna(j) == 'isobutane') THEN
      kij(i,j) = 0.051
    ELSE IF(compna(i) == 'methanol' .AND. compna(j) == '1-octanol') THEN
      kij(i,j) = 0.02
    ELSE IF(compna(i) == '1-butanol' .AND. compna(j) == 'butane') THEN
      kij(i,j) = 0.015
    ELSE IF(compna(i) == '1-nonanol' .AND. compna(j) == 'co2') THEN
      kij(i,j) = 0.076
      kij(i,j) = 0.01443
    ELSE IF(compna(i) == '1-propanol' .AND. compna(j) == 'ethylbenzene') THEN
      kij(i,j) = 0.0251757813
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'ethanol') THEN
      kij(i,j) = 0.085
    ELSE IF(compna(i) == 'hexane'.AND.compna(j) == '1-chlorobutane') THEN
      kij(i,j) = 0.017
    ELSE IF(compna(i) == 'aniline'.AND.compna(j) == 'methylcyclopentane') THEN
      kij(i,j) = 0.0153
    ELSE IF(compna(i) == 'heptane'.AND.compna(j) == 'nbutyl-ethanoate')THEN
      kij(i,j) = 0.027
    ELSE IF(compna(i) == '1-hexene'.AND.compna(j) == 'ethyl-ethanoate')THEN
      kij(i,j) = 0.026
    ELSE IF(compna(i) == 'co2'.AND.compna(j) == '1-butanol')THEN
      ! kij(i,j) = 0.075
      kij(i,j) = 0.0
    ELSE IF(compna(i) == 'acrylic-acid'.AND.compna(j) == 'co2')THEN
      kij(i,j) = -0.1
    ELSE IF(compna(i) == 'bmim'.AND.compna(j) == 'hydrogen')THEN
      lij(i,j) =  0.55   !!!!! Lij not kij !!!!
    ELSE IF(compna(i) == 'bf4'.AND.compna(j) == 'hydrogen')THEN
      lij(i,j) =  0.55   !!!!! Lij not kij !!!!
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'butane')THEN !K.R.Hall FPE 2007 254 112-125 kij=0.1850
      kij(i,j) = -0.07
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == '1-butanol')THEN
      kij(i,j) = -0.12
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'aniline')THEN
      kij(i,j) = 0.0
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'methanol') THEN
      kij(i,j) = -0.02
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'ethanol') THEN
      kij(i,j) = -0.027
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'styrene') THEN
      kij(i,j) = 0.1
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'propyl-ethanoate') THEN
      kij(i,j) = -0.205
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'ethyl-propanoate') THEN
      kij(i,j) = 0.0
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == '1-pentanol') THEN
      kij(i,j) = 0.0165
      ! kij(i,j) = 0.0294
      ! kij(i,j) = -0.082
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'methane') THEN
      ! kij(i,j) = +0.06
      kij(i,j) = -0.08
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'propane') THEN
      kij(i,j) = 0.02
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'hexane') THEN
      kij(i,j) = -0.02
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'acetic-acid') THEN
      kij(i,j) = -0.0
    ELSE IF(compna(i) == 'h2o'.AND.compna(j) == 'co2') THEN
      if (pol == 0) kij(i,j) = 0.0030625  ! for T=50C, according to X.Tang
      stop                                ! very T-dependent
    ELSE IF(compna(i) == 'toluene'.AND.compna(j) == 'acetic-acid') THEN
      kij(i,j) = -0.1
    ELSE IF(compna(i) == 'caproic-acid'.AND.compna(j) == 'cyclohexane') THEN
      kij(i,j) = 0.041531
    ELSE IF(compna(i) == '1-octanol'.AND.compna(j) == 'h2o')THEN
      kij(i,j) = 0.07
    ELSE IF(compna(i) == 'acetone'.AND.compna(j) == 'benzene')THEN
      kij(i,j) = 0.02132466    ! PC-SAFT
      ! kij(i,j) = 0.01495148    ! PCP-SAFT
      ! kij(i,j) = -0.00735575    ! PCP-SAFT but non-polar benzene
    ELSE IF(compna(i) == '1-propanol'.AND.compna(j) == 'benzene')THEN
      kij(i,j) = 0.02017
    ELSE IF(compna(i) == '2-propanol'.AND.compna(j) == 'benzene')THEN
      kij(i,j) = 0.021386
    ELSE IF(compna(i) == '1-pentanol'.AND.compna(j) == 'benzene')THEN
      kij(i,j) = 0.0129638671875
    ELSE IF(compna(i) == 'CH2F2' .AND. compna(j) == 'co2') THEN
      kij(i,j) =  2.2548828125E-2
    ELSE IF(compna(i) == 'dmso' .AND. compna(j) == 'co2') THEN
      kij(i,j) =  -0.00
    ELSE IF(compna(i) == 'dmf'.AND.compna(j) == 'h2o')THEN
      kij(i,j) =  -0.0
    ELSE IF(compna(i) == 'decane'.AND.compna(j) == 'h2o')THEN
      kij(i,j) =  0.11
    ELSE IF(compna(i) == '11difluoroethane'.AND.compna(j) == 'HF')THEN
      kij(i,j) =  0.032 !  association:   eps_kij = 0.16
    ELSE IF(compna(i) == '11difluoroethane'.AND.compna(j) == 'co2')THEN
      kij(i,j) =  -0.004 !  PCP-SAFT (taken from simulation)
    ELSE IF(compna(i) == 'difluoromethane'.AND.compna(j) == 'HF')THEN
      kij(i,j) = -0.02
    ELSE IF(compna(i) == 'naphthalene'.AND.compna(j) == 'co2')THEN
      kij(i,j) =  0.137  ! angepasst an SVLE-Linie (tradition.CO2-Parameter)
      kij(i,j) =  0.09  ! angepasst an SVLE-Linie (tradition.CO2-Parameter)
    ELSE IF(compna(i) == 'pg2'.AND.compna(j) == 'methanol')THEN
      kij(i,j) =  0.03
    ELSE IF(compna(i) == 'pg2'.AND.compna(j) == 'co2')THEN
      ! kij(i,j) =  0.05
    ELSE IF(compna(i) == 'PCH5'.AND.compna(j) == 'co2')THEN
      ! kij(i,j) =  -0.047
      kij(i,j) =  +0.055
      ! kij(i,j) =  +0.036
    ELSE
    END IF
    kij(j,i) = kij(i,j)
    lij(j,i) = -lij(i,j)
    
  END DO
END DO

END SUBROUTINE pcsaft_par


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE init_vars
!
! This subroutine writes variables from an array "val_init" to the
! system-variables. Those variables
! include some specifications but also some starting values for a
! phase equilibrium calculation. (val_init stands for 'initial value')

! The array comprises the following elements
! element 0              density of third phase
! element 1              density of first phase
! element 2              density of second phase
! element 3              temperature [K]
! element 4              pressure [Pa]
! element 5,..(5+NCOMP)  mole-fraction of comp. i in log. from (phase 1)
! element ...            mole-fraction of comp. i in log. from (phase 2)
! element ...            mole-fraction of comp. i in log. from (phase 3)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE init_vars
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
 INTEGER                                :: i, ph
! ----------------------------------------------------------------------

densta(3)=val_init(0)
densta(1)=val_init(1)
densta(2)=val_init(2)
t = val_init(3)
p = val_init(4)
DO ph = 1,nphas
  DO i = 1,ncomp
    lnx(ph,i) = val_init(4+i+(ph-1)*ncomp)
  END DO
END DO

END SUBROUTINE init_vars


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE converged
!
! Once a converged solution for an equilibrium calculation is found,
! this subroutine writes variables to an array "val_conv".
! (= short for converged values)
! The array comprises the following elements
! element 0              density of third phase
! element 1              density of first phase
! element 2              density of second phase
! element 3              temperature [K]
! element 4              pressure [Pa]
! element 5,..(4+NCOMP)  mole-fraction of comp. i in log. from (phase 1)
! element ...            mole-fraction of comp. i in log. from (phase 2)
! element ...            mole-fraction of comp. i in log. from (phase 3)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE converged
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
 INTEGER                                :: i, ph
! ----------------------------------------------------------------------

val_conv(0)  = dense(3)
val_conv(1)  = dense(1)
val_conv(2)  = dense(2)
val_conv(3)  = t
val_conv(4)  = p
DO ph = 1,nphas
  DO i = 1,ncomp
    val_conv(4+i+(ph-1)*ncomp) = lnx(ph,i)
  END DO
END DO

END SUBROUTINE converged

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE PERTURBATION_PARAMETER
!
! calculates density-independent parameters of the equation of state.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE PERTURBATION_PARAMETER
!
 USE PARAMETERS, ONLY: PI, KBOL, RGAS, NAV, TAU
 USE EOS_VARIABLES
 USE EOS_CONSTANTS
 IMPLICIT NONE
!
! --- local variables --------------------------------------------------
 INTEGER                                :: i, j, k, l, m, no
 LOGICAL                                :: assoc, qudpole, dipole
 REAL                                   :: m_mean
 REAL, DIMENSION(nc)                    :: v00, v0, d00, u
 REAL, DIMENSION(nc,nc,nsite,nsite)     :: eps_hb
 REAL, DIMENSION(nc,nc)                 :: kap_hb
 REAL                                   :: zmr, nmr
 REAL                                   :: eps_kij, k_kij
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! pure component parameters
! ----------------------------------------------------------------------
DO  i = 1, ncomp
  u(i)     = parame(i,3) * (1.0 + parame(i,4)/t)
  mseg(i)  = parame(i,1)
  IF (eos == 0) THEN
    v00(i) = parame(i,2)
    v0(i)  = v00(i)*(1.0-0.12*EXP(-3.0*parame(i,3)/t))**3
    d00(i) = (1.d30/1.d6 *tau *v00(i)*6.0/PI /NAV)**(1.0/3.0)
    dhs(i) = d00(i)*(1.0-0.12*EXP(-3.0*parame(i,3)/t))
  ELSE
    dhs(i) = parame(i,2)*(1.0-0.12*EXP(-3.0*parame(i,3)/t))
    d00(i) = parame(i,2)
  END IF
END DO


! ----------------------------------------------------------------------
! combination rules
! ----------------------------------------------------------------------
DO  i = 1, ncomp
  DO  j = 1, ncomp
    sig_ij(i,j) = 0.5 * ( d00(i) + d00(j) )
    uij(i,j) = ( 1.0 - kij(i,j) ) * ( u(i)*u(j) )**0.5
    vij(i,j) = ( 0.5*( v0(i)**(1.0/3.0) + v0(j)**(1.0/3.0) ) )**3 
  END DO
END DO


! ----------------------------------------------------------------------
! abbreviations
! ----------------------------------------------------------------------
z0t = PI / 6.0 * SUM( x(1:ncomp) * mseg(1:ncomp) )
z1t = PI / 6.0 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp) )
z2t = PI / 6.0 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**2 )
z3t = PI / 6.0 * SUM( x(1:ncomp) * mseg(1:ncomp) * dhs(1:ncomp)**3 )

m_mean = z0t/(PI/6.0)

DO i = 1, ncomp
  DO j = 1, ncomp
    dij_ab(i,j) = dhs(i)*dhs(j) / ( dhs(i) + dhs(j) )
  END DO
END DO

! ----------------------------------------------------------------------
! dispersion term parameters for chain molecules
! ----------------------------------------------------------------------
DO m = 0, 6
  apar(m) = ap(m,1) + (1.0-1.0/m_mean)*ap(m,2)  &
            + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*ap(m,3)
  bpar(m) = bp(m,1) + (1.0-1.0/m_mean)*bp(m,2)  &
            + (1.0-1.0/m_mean)*(1.0-2.0/m_mean)*bp(m,3)
END DO


! ----------------------------------------------------------------------
! van der Waals mixing rules for perturbation terms
! ----------------------------------------------------------------------
order1 = 0.0
order2 = 0.0
DO i = 1, ncomp
  DO j = 1, ncomp
    order1 = order1 + x(i)*x(j)* mseg(i)*mseg(j)*sig_ij(i,j)**3 * uij(i,j)/t
    order2 = order2 + x(i)*x(j)* mseg(i)*mseg(j)*sig_ij(i,j)**3 * (uij(i,j)/t)**2
  END DO
END DO


! ----------------------------------------------------------------------
! SAFT parameters
! ----------------------------------------------------------------------
IF (eos == 0) THEN
  zmr = 0.0
  nmr = 0.0
  DO  i = 1, ncomp
    DO  j = 1, ncomp
      zmr = zmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)*uij(i,j)
      nmr = nmr + x(i)*x(j)*mseg(i)*mseg(j)*vij(i,j)
    END DO
  END DO
  um = zmr / nmr
END IF



! ----------------------------------------------------------------------
! association and polar parameters
! ----------------------------------------------------------------------
assoc   = .false.
qudpole = .false.
dipole  = .false.
DO i = 1, ncomp
  IF (NINT(parame(i,12)) /= 0) assoc  = .true.
  IF (parame(i,7) /= 0.0)     qudpole = .true.
  IF (parame(i,6) /= 0.0)     dipole  = .true.
END DO

! --- dipole and quadrupole constants ----------------------------------
IF (qudpole) CALL qq_const ( qqp2, qqp3, qqp4 )
IF (dipole)  CALL dd_const ( ddp2, ddp3, ddp4 )
IF (dipole .AND. qudpole) CALL dq_const ( dqp2, dqp3, dqp4 )


! --- TPT-1-association parameters -------------------------------------
IF (assoc) THEN
  
  eps_kij = 0.0
  k_kij   = 0.0
  
  DO i = 1, ncomp
    IF (NINT(parame(i,12)) /= 0) THEN
      nhb_typ(i)  = NINT(parame(i,12))
      kap_hb(i,i) = parame(i,13)
      no = 0
      DO j = 1, nhb_typ(i)
        DO k = 1, nhb_typ(i)
          eps_hb(i,i,j,k) = parame(i,(14+no))
          no=no+1
        END DO
      END DO
      DO j = 1, nhb_typ(i)
        nhb_no(i,j) = parame(i,(14+no))
        no=no+1
      END DO
    ELSE
      nhb_typ(i) = 0
      kap_hb(i,i)= 0.0
      DO k = 1, nsite
        DO l = 1, nsite
          eps_hb(i,i,k,l) = 0.0
        END DO
      END DO
    END IF
  END DO
  
  DO i = 1, ncomp
    DO j = 1, ncomp
      IF (i /= j .AND. (nhb_typ(i) /= 0 .AND. nhb_typ(j) /= 0)) THEN
        kap_hb(i,j)= (kap_hb(i,i)*kap_hb(j,j))**0.5  &
                     *((parame(i,2)*parame(j,2))**3 )**0.5  &
                     /(0.5*(parame(i,2)+parame(j,2)))**3 
        kap_hb(i,j)= kap_hb(i,j)*(1.0-k_kij)
        DO k = 1, nhb_typ(i)
          DO l = 1, nhb_typ(j)
            IF (k /= l) THEN
              eps_hb(i,j,k,l) = (eps_hb(i,i,k,l)+eps_hb(j,j,l,k))/2.0
              eps_hb(i,j,k,l) = eps_hb(i,j,k,l)*(1.0-eps_kij)
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
  IF (nhb_typ(1) == 3) THEN
!        write(*,*)'eps_hb manuell eingegeben'
    eps_hb(1,2,3,1) = 0.5*(eps_hb(1,1,3,1)+eps_hb(2,2,1,2))
    eps_hb(2,1,1,3) = eps_hb(1,2,3,1)
  END IF
  IF (nhb_typ(2) == 3) THEN
    eps_hb(2,1,3,1) = 0.5*(eps_hb(2,2,3,1)+eps_hb(1,1,1,2))
    eps_hb(1,2,1,3) = eps_hb(2,1,3,1)
  END IF
  
  DO i = 1, ncomp
    DO k = 1, nhb_typ(i)
      DO j = 1, ncomp
        DO l = 1, nhb_typ(j)
          ass_d(i,j,k,l) = kap_hb(i,j) *sig_ij(i,j)**3 *(EXP(eps_hb(i,j,k,l)/t)-1.0)
        END DO
      END DO
    END DO
  END DO
  
END IF

END SUBROUTINE PERTURBATION_PARAMETER




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE OUTPUT
!
! The purpose of this subroutine is obvious.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE OUTPUT
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i
 CHARACTER (LEN=4)                      :: t_ind
 CHARACTER (LEN=4)                      :: p_ind
 CHARACTER (LEN=4)                      :: char_ncomp
 REAL                                   :: density(np),w(np,nc)
! ----------------------------------------------------------------------

 CALL si_dens (density,w)

 IF (u_in_p == 1.E5) THEN
    p_ind = ' bar'
 ELSE IF (u_in_p == 1.E2) THEN
    p_ind = 'mbar'
 ELSE IF (u_in_p == 1.E6) THEN
    p_ind = ' MPa'
 ELSE IF (u_in_p == 1.E3) THEN
    p_ind = ' kPa'
 ELSE
    p_ind = ' Pa'
 END IF
 IF (u_in_t == 273.15) THEN
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

 !-----output to files--------------------------------------------------
 IF (ncomp == 1) THEN
    WRITE (40,'(7(2x,f18.8))') t-u_out_t, p/u_out_p,  &
         density(1), density(2),h_lv,cpres(1),cpres(2)
    !     &      ,(enthal(2)-enthal(1))/mm(1)
    !        WRITE (40,'(4(2x,f15.8))') t, p, 0.3107*dense(1)
    !     &  /0.1617*(0.689+0.311*(T/1.328)**0.3674),0.3107
    !     &  *dense(2)/0.1617*(0.689+0.311*(T/1.328)**0.3674)
 ELSE IF (ncomp == 2) THEN
    WRITE (40,'(12(2x,G15.8))') 1.0-xi(1,1),1.0-xi(2,1),  &
         w(1,1),w(2,1),t-u_out_t, p/u_out_p, density(1),density(2)  &
         ,enthal(1),enthal(2),cpres(1),cpres(2)
 ELSE IF (ncomp == 3) THEN
    WRITE (40,'(10(2x,f15.8))') xi(1,1),xi(1,2),xi(1,3),xi(2,1),xi(2,2),  &
         xi(2,3),t-u_out_t, p/u_out_p, density(1),density(2)
 END IF

 END SUBROUTINE OUTPUT



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE neutr_charge
!
! This subroutine is called for electrolye solutions, where a
! neutral overall-charge has to be enforced in all phases. The basic
! philosophy is similar to the above described routine X_SUMMATION.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE neutr_charge(i)
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: i
!
! ----------------------------------------------------------------------
 INTEGER                                :: comp_e, ph_e
 REAL                                   :: sum_c
 CHARACTER (LEN=2)                      :: phasno
 CHARACTER (LEN=2)                      :: compno
! ----------------------------------------------------------------------

phasno = sum_rel(i)(2:2)
READ(phasno,*) ph_e
compno = sum_rel(i)(3:3)
READ(compno,*) comp_e

sum_c = 0.0
write (*,*) 'there must be an error in neutr_charge'
stop
! there is an error in the following passage. The index i is an
! argument to this subroutine - I guess it is INTENT(IN), so the
! index in the following loop can not be i.
!
! I have commented the loop until I check the code.
!DO i=1,ncomp
!  IF ( comp_e /= i .AND. parame(i,10) /= 0.0)  &
!      sum_c = sum_c + xi(ph_e,i)*parame(i,10)
!END DO

xi(ph_e,comp_e) = - sum_c
IF (xi(ph_e,comp_e) < 0.0) xi(ph_e,comp_e)=0.0
IF (xi(ph_e,comp_e) /= 0.0) THEN
  lnx(ph_e,comp_e) = LOG(xi(ph_e,comp_e))
ELSE
  lnx(ph_e,comp_e) = -100000.0
END IF

! xi(2,1) = xi(2,2)
! IF (xi(2,1).NE.0.0) lnx(2,1) = LOG(xi(2,1))

END SUBROUTINE neutr_charge


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE flash_sum
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
 INTEGER                                :: i, j, ph_i, phase1, phase2
! ----------------------------------------------------------------------

phase1=0
phase2=0
DO j=1,ncomp
  IF (it(j)(2:2) == '1') phase1=phase1+1
  IF (it(j)(2:2) == '2') phase2=phase2+1
END DO

IF (phase1 == ncomp-1) THEN
  ph_i = 1
ELSE IF (phase2 == ncomp-1) THEN
  ph_i = 2
ELSE
  WRITE (*,*) ' FLASH_SUM: undefined flash-case'
  STOP
END IF



IF (ph_i == 1) THEN
  DO i=1,ncomp
    IF (alpha > DMIN1(1.0,xif(i)/xi(1,i),  &
          (xif(i)-1.0)/(xi(1,i)-1.0),alpha)) THEN
      WRITE (*,*) ' FLASH_SUM: exeeded 1st alpha-bound'
      alpha=DMIN1(1.0,xif(i)/xi(1,i),(xif(i)-1.0)/(xi(1,i)-1.0))
    END IF
  END DO
  DO i=1,ncomp
    xi(2,i) = ( xif(i) - alpha*xi(1,i) ) / (1.0-alpha)
    IF (xi(2,i) > 0.0) THEN
      lnx(2,i) = LOG(xi(2,i))
    ELSE
      xi(2,i) = 0.0
      lnx(2,i) = -100000.0
    END IF
  END DO
ELSE IF (ph_i == 2) THEN
  DO i=1,ncomp
    IF (alpha > DMAX1(0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)),  &
          1.0-xif(i)/xi(2,i),alpha)) THEN
      WRITE (*,*) ' FLASH_SUM: exeeded 2nd alpha-bound'
      WRITE (*,*) 0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)), 1.0-xif(i)/xi(2,i)
      alpha=DMAX1(0.0,(xif(i)-xi(2,i))/(1.0-xi(2,i)), 1.0-xif(i)/xi(2,i))
    END IF
  END DO
  DO i=1,ncomp
    xi(1,i) = ( xif(i) - (1.0-alpha)*xi(2,i) ) / alpha
!    write (*,*) 'x1,i',xi(1,i),xi(2,i),alpha
    IF (xi(1,i) > 0.0) THEN
      lnx(1,i) = LOG(xi(1,i))
    ELSE
      xi(1,i) = 0.0
      lnx(1,i) = -100000.0
    END IF
  END DO
END IF

! pause

RETURN
END SUBROUTINE flash_sum

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE flash_alpha
!
! This subroutine calculates all molefractions of one phase
! from a component balance. What is needed for this calculation
! are all molefractions of the other phase (nphas=2, so far)
! and the phase fraction alpha.
! Alpha is calculated from the mole fraction
! of component {sum_rel(j)(3:3)}. If for example sum_rel(2)='fl3',
! then the alpha is determined from the molefraction of comp. 3 and
! all molefractions of one phase are calculated using that alpha-value.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE flash_alpha
!
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, j, comp_i, phase1, phase2
 CHARACTER (LEN=2)                      :: compno
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! first calculate the phase fraction alpha from a known composition
! of component sum_rel(j)(3:3).
! ----------------------------------------------------------------------

DO j = 1, nphas
  IF ( sum_rel(j)(1:2) == 'fl' ) THEN
    compno = sum_rel(j)(3:3)
    READ(compno,*) comp_i
    IF ( (xi(1,comp_i)-xi(2,comp_i)) /= 0.0 ) THEN
      alpha = (xif(comp_i)-xi(2,comp_i)) / (xi(1,comp_i)-xi(2,comp_i))
      write (*,*) 'flsh',(xif(comp_i)-xi(2,comp_i)),(xi(1,comp_i)-xi(2,comp_i))
    ELSE
      alpha = 0.5
      WRITE (*,*) 'FLASH_ALPHA:error in calc. of phase fraction',comp_i
    END IF
    ! IF (alpha <= 0.0 .OR. alpha >= 1.0) WRITE(*,*) 'FLASH_ALPHA: error',alpha
    IF (alpha > 1.0) alpha = 1.0
    IF (alpha < 0.0) alpha = 0.0
  END IF
END DO

! ----------------------------------------------------------------------
! determine which phase is fully determined by iterated molefractions (+ summation relation)
! ----------------------------------------------------------------------
phase1 = 0
phase2 = 0
DO i = 1, ncomp
  IF ( it(i)(2:2) == '1' ) phase1 = phase1 + 1
  IF ( it(i)(2:2) == '2' ) phase2 = phase2 + 1
END DO
DO i = 1, ncomp
  IF ( sum_rel(i)(2:2) == '1' ) phase1 = phase1 + 1
  IF ( sum_rel(i)(2:2) == '2' ) phase2 = phase2 + 1
END DO


IF ( phase1 == ncomp ) THEN
  ! --------------------------------------------------------------------
  ! phase 1 is defined by iterated molefractions + summation relation
  ! phase 2 is determined from componennt balance (using alpha)
  ! --------------------------------------------------------------------
  IF ( alpha == 1.0 ) alpha = 1.0 - 1.0E-10
  DO i=1,ncomp
    xi(2,i) = ( xif(i) - alpha*xi(1,i) ) / (1.0-alpha)
    IF ( xi(2,i) < 0.0 ) xi(2,i) = 0.0
    IF ( xi(2,i) > 1.0 ) xi(2,i) = 1.0
    IF ( xi(2,i) /= 0.0 ) THEN
      lnx(2,i) = LOG( xi(2,i) )
    ELSE
      lnx(2,i) = -100000.0
    END IF
    write (*,*) 'fl_cal ph=2',i,lnx(2,i),xi(2,i)
  END DO
ELSE IF ( phase2 == ncomp ) THEN
  ! --------------------------------------------------------------------
  ! phase 2 is defined by iterated molefractions + summation relation
  ! phase 1 is determined from componennt balance (using alpha)
  ! --------------------------------------------------------------------
  DO i = 1, ncomp
    xi(1,i) = ( xif(i) - (1.0-alpha)*xi(2,i) ) /alpha
    IF ( xi(1,i) < 0.0 ) xi(1,i) = 0.0
    IF ( xi(1,i) > 1.0 ) xi(1,i) = 1.0
    IF ( xi(1,i) /= 0.0 ) THEN
      lnx(1,i) = LOG( xi(1,i) )
    ELSE
      lnx(1,i) = -100000.0
    END IF
    write (*,*) 'fl_cal ph=1',i,lnx(1,i),xi(1,i),alpha
  END DO
ELSE
  WRITE (*,*) ' FLASH_ALPHA: undefined flash-case'
  STOP
END IF

END SUBROUTINE flash_alpha


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE SI_DENS (density,w)
!
! This subroutine calculates the (macroskopic) fluid-density in
! units [kg/m3] from the dimensionless density (eta=zeta3).
! Further, mass fractions are calculated from mole fractions.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE SI_DENS (density,w)
!
 USE PARAMETERS, ONLY: pi, nav, tau
 USE BASIC_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(OUT)                      :: density(np)
 REAL, INTENT(OUT)                      :: w(np,nc)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, ph
 REAL                                   :: mm_mean, rho, z3t
 REAL                                   :: dhs(nc), d00(nc), t_p, pcon, l_st
! ----------------------------------------------------------------------


DO i = 1,ncomp
  IF (eos == 1) THEN
    dhs(i) = parame(i,2) * ( 1.0 - 0.12 *EXP( -3.0*parame(i,3)/t ) )
  ELSE IF (eos == 0) THEN
    d00(i) = ( 1.E30/1.E6*tau*parame(i,2)*6.0/pi/nav )**(1.0/3.0)
    dhs(i) = d00(i) * ( 1.0 - 0.12 *EXP( -3.0*parame(i,3)/t ) )
  ELSE IF (eos == 4) THEN
    dhs(i) = (  0.1617/0.3107 / ( 0.689+0.311*(t/parame(i,3)/1.328)**0.3674 )  &
        / ( pi/6.0 )  )**(1.0/3.0) * parame(i,2)
  ELSE IF (eos == 5.OR.eos == 6) THEN
    l_st = parame(1,25)
    IF (ncomp /= 1) write (*,*) ' ERROR for EOS = 5'
    t_p  =((34.037352+17.733741*l_st) /(1.0+0.53237307*l_st+12.860239*l_st**2 ))**0.5
    IF (l_st == 0.0) t_p = t_p/4.0
    IF (eos == 5 .AND. l_st /= 0.0) t_p = t_p/4.0*parame(1,1)**2 
    t_p = t/parame(i,3)/t_p
    pcon =0.5256+3.2088804*l_st**2 -3.1499114*l_st**2.5 +0.43049357*l_st**4 
    dhs(i) = (  pcon/(0.67793+0.32207*(t_p)**0.3674) /(pi/6.0)  )**(1.0/3.0) *parame(i,2)
  ELSE IF (eos == 8) THEN
    dhs(i) = parame(i,2)*(1.0+0.2977*t/parame(i,3))  &
        /(1.0+0.33163*t/parame(i,3) +1.0477E-3*(t/parame(i,3))**2 )
  ELSE
    write (*,*) 'define EOS in subroutine: SI_DENS'
    stop
  END IF
END DO

DO ph = 1,nphas
  mm_mean = 0.0
  z3t = 0.0
  DO i = 1, ncomp
    mm_mean = mm_mean + xi(ph,i)*mm(i)
    z3t = z3t + xi(ph,i) * parame(i,1) * dhs(i)**3
  END DO
  z3t = pi/6.0 * z3t
  rho = dense(ph)/z3t
  density(ph) = rho * mm_mean * 1.E27 /nav
  DO i = 1, ncomp
    w(ph,i) = xi(ph,i) * mm(i)/mm_mean
  END DO
!  write (*,*) density(ph),rho,mm_mean,1.d27 /NAV
END DO

END SUBROUTINE SI_DENS



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 REAL FUNCTION F_STABILITY ( optpara, n )
!
 USE BASIC_VARIABLES
 USE STARTING_VALUES
 USE EOS_VARIABLES, ONLY: dhs, PI
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN OUT)                   :: optpara(n)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, stabil
 REAL                                   :: rhoi(nc),gradterm
 REAL                                   :: fden,punish
 REAL                                   :: dens
! ----------------------------------------------------------------------

COMMON /stabil / stabil

punish = 0.0
stabil = 1

DO i = 1, n
  IF ( optpara(i) < 0.5 ) rhoi(i) = EXP(optpara(i) )
  IF ( optpara(i) >= 0.5) rhoi(i) = EXP(0.5)
END DO

dens = PI/6.0 * SUM( rhoi(1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )

IF (dens > 0.65) THEN
  punish = punish + (dens-0.65)*10000.0
  rhoi(1:n) = rhoi(1:n)*0.65/dens
END IF

CALL fden_calc (fden, rhoi)

gradterm = sum( grad_fd(1:n) * ( rhoi(1:n) - rhoif(1:n) ) )

f_stability = fden - fdenf - gradterm + punish

! write (*,'(5G16.8)') F_STABILITY,(rhoi(i),i=1,n)
! pause

stabil = 0

END FUNCTION F_STABILITY




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE p_calc (pges_transfer, zges)
!
! This subroutine serves as an iterface to the EOS-routines. The
! system pressure corresponding to given (desity,T,xi) is calculated.
! (Note: the more common interface is SUBROUTINE FUGACITY. This
! routine is only used for one-phase systems, e.g. calculation of
! virial coefficients)
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE p_calc (pges_transfer, zges)
!
 USE BASIC_VARIABLES
 USE EOS_VARIABLES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: pges_transfer
 REAL, INTENT(OUT)                      :: zges
! ----------------------------------------------------------------------

IF (nphas /= 1 ) THEN
  write (*,*) 'P_CALC: can only be called for single phases'
  stop
ENDIF

IF (eos < 2) THEN

  phas = 1
  eta = dense(1)
  x(1:ncomp) = xi(1,1:ncomp)
  
  CALL PERTURBATION_PARAMETER
  IF (num == 0) CALL P_EOS
  IF(num == 1) CALL P_NUMERICAL
 !! IF(num == 2) CALL F_EOS_RN

  pges_transfer = pges
  rho = eta/z3t
  zges = (pges * 1.E-30) / (kbol*t*rho)
  
ELSE
  write (*,*) ' SUBROUTINE P_CALC not available for cubic EOS'
  stop
END IF

END SUBROUTINE p_calc



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
SUBROUTINE ONLY_ONE_TERM_EOS_NUMERICAL ( only_term, type_of_term )
!
 USE EOS_NUMERICAL_DERIVATIVES
 IMPLICIT NONE
!
  character (LEN=9)                        :: only_term, type_of_term
! ----------------------------------------------------------------------

  save_eos_terms(1) = ideal_gas
  save_eos_terms(2) = hard_sphere
  save_eos_terms(3) = chain_term
  save_eos_terms(4) = disp_term
  save_eos_terms(5) = hb_term
  save_eos_terms(6) = LC_term
  save_eos_terms(7) = branch_term
  save_eos_terms(8) = II_term
  save_eos_terms(9) = ID_term

  ideal_gas   = 'no'
  hard_sphere = 'no'
  chain_term  = 'no'
  disp_term   = 'no'
  hb_term     = 'no'
  LC_term     = 'no'
  branch_term = 'no'
  II_term     = 'no'
  ID_term     = 'no'

  IF ( only_term == 'ideal_gas' )   ideal_gas   = trim( adjustl( type_of_term ) )
  IF ( only_term == 'hard_sphere' ) hard_sphere = trim( adjustl( type_of_term ) )
  IF ( only_term == 'chain_term' )  chain_term  = trim( adjustl( type_of_term ) )
  IF ( only_term == 'disp_term' )   disp_term   = trim( adjustl( type_of_term ) )
  IF ( only_term == 'hb_term' )     hb_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'LC_term' )     LC_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'branch_term' ) branch_term = trim( adjustl( type_of_term ) )
  IF ( only_term == 'II_term' )     II_term     = trim( adjustl( type_of_term ) )
  IF ( only_term == 'ID_term' )     ID_term     = trim( adjustl( type_of_term ) )

END SUBROUTINE ONLY_ONE_TERM_EOS_NUMERICAL


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
SUBROUTINE RESTORE_PREVIOUS_EOS_NUMERICAL
!
 USE EOS_NUMERICAL_DERIVATIVES
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
  ideal_gas   = trim( adjustl( save_eos_terms(1) ) )
  hard_sphere = trim( adjustl( save_eos_terms(2) ) )
  chain_term  = trim( adjustl( save_eos_terms(3) ) )
  disp_term   = trim( adjustl( save_eos_terms(4) ) )
  hb_term     = trim( adjustl( save_eos_terms(5) ) )
  LC_term     = trim( adjustl( save_eos_terms(6) ) )
  branch_term = trim( adjustl( save_eos_terms(7) ) )
  II_term     = trim( adjustl( save_eos_terms(8) ) )
  ID_term     = trim( adjustl( save_eos_terms(9) ) )

END SUBROUTINE RESTORE_PREVIOUS_EOS_NUMERICAL




