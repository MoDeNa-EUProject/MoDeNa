module besselj
    implicit none
    private
    public ZBESJ
contains
!*****************************************************************
!* EVALUATE A J-BESSEL FUNCTION OF COMPLEX ARGUMENT (FIRST KIND) *
!* ------------------------------------------------------------- *
!* SAMPLE RUN:                                                   *
!* (Evaluate J0 to J4 for argument Z=(1.0,2.0) ).                *
!*                                                               *
!* zr(0) =   1.586259                                            *
!* zi(0) =  -1.391602                                            *
!* zr(1) =   1.291848                                            *
!* zi(1) =   1.010488                                            *
!* zr(2) =  -0.261130                                            *
!* zi(2) =   0.762320                                            *
!* zr(3) =  -0.281040                                            *
!* zi(3) =   0.017175                                            *
!* zr(4) =  -0.034898                                            *
!* zi(4) =  -0.067215                                            *
!* NZ =          0                                               *
!* Error code:           0                                       *
!*                                                               *
!* ------------------------------------------------------------- *
!* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77    *
!*       [BIBLI 18].                                             *
!*                                                               *
!*                        F90 Release 1.0 By J-P Moreau, Paris   *
!*                                 (www.jpmoreau.fr)             *
!*****************************************************************
!Note: to link with Complex.f90 and Utilit.f90.
!---------------------------------------------
!PROGRAM TEST_ZBESJ
!real*8 zr,zi, cyr(10), cyi(10)

!  n=5
!  zr=1.d0; zi=2.d0
!  call ZBESJ(zr,zi,0,1,n,cyr,cyi,nz,ierr)

!  print *,' '
!  do i=1, n
!    write(*,10) i-1, cyr(i)
!    write(*,11) i-1, cyi(i)
!  end do
!  print *,' NZ=', NZ
!  print *,' Error code:', ierr
!  print *,' '

!10 format('  zr(',I1,') = ',F10.6)
!11 format('  zi(',I1,') = ',F10.6)
!
!END

SUBROUTINE ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
USE UTILIT  !For I1MACH,D1MACH.
USE COMPLEX !For ZABS
!***BEGIN PROLOGUE  ZBESJ
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
!             BESSEL FUNCTION OF FIRST KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
!***DESCRIPTION
!
!                      ***A DOUBLE PRECISION ROUTINE***
!         ON KODE=1, CBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
!         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
!         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
!         -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESJ RETURNS THE SCALED
!         FUNCTIONS
!
!         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
!
!         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
!         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
!         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
!         (REF. 1).
!
!         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
!           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
!           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0D0
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE= 1  RETURNS
!                             CY(I)=J(FNU+I-1,Z), I=1,...,N
!                        = 2  RETURNS
!                             CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
!           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
!
!         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
!           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
!                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
!                    CY(I)=J(FNU+I-1,Z)  OR
!                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N
!                    DEPENDING ON KODE, Y=AIMAG(Z).
!           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!                    NZ= 0   , NORMAL RETURN
!                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE
!                             TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
!                              I = N-NZ+1,...,N
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
!                            TOO LARGE ON KODE=1
!                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
!                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
!                            ACCURACY
!                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
!                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
!                            CANCE BY ARGUMENT REDUCTION
!                    IERR=5, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!***LONG DESCRIPTION
!
!         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
!
!         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
!
!         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
!
!         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
!
!         FOR NEGATIVE ORDERS,THE FORMULA
!
!              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
!
!         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
!         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
!         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
!         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
!         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
!         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
!         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
!         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
!         LARGE MEANS FNU.GT.CABS(Z).
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
!         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
!         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
!         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
!         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
!         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
!         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
!         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
!         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
!         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
!         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
!         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
!         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
!         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
!         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
!         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
!         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!         OR -PI/2+P.
!
!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
!                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
!                 COMMERCE, 1955.
!
!               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
!
!               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
!
!               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!***ROUTINES CALLED  ZABS,ZBINU,I1MACH,D1MACH
!***END PROLOGUE  ZBESJ
!
!     COMPLEX CI,CSGN,CY,Z,ZN
      DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG, &
      ELIM, FNU, FNUL, HPI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,       &
      BB, FN, AZ
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ
      DIMENSION CYR(1), CYI(1)
      DATA HPI /1.57079632679489662D0/

!***FIRST EXECUTABLE STATEMENT  ZBESJ
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
!-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
!-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+DBLE(FLOAT(N-1))
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
!-----------------------------------------------------------------------
!     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
      CII = 1.0D0
      INU = INT(SNGL(FNU))
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-DBLE(FLOAT(INU-IR)))*HPI
      CSGNR = DCOS(ARG)
      CSGNI = DSIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE
!-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      IF (ZI.GE.0.0D0) GO TO 50
      ZNR = -ZNR
      ZNI = -ZNI
      CSGNI = -CSGNI
      CII = -CII
   50 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, &
      ELIM, ALIM)
      IF (NZ.LT.0) GO TO 130
      NL = N - NZ
      IF (NL.EQ.0) RETURN
      DO 60 I=1,NL
        STR = CYR(I)*CSGNR - CYI(I)*CSGNI
        CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
        CYR(I) = STR
        STR = -CSGNI*CII
        CSGNI = CSGNR*CII
        CSGNR = STR
   60 CONTINUE
      RETURN
  130 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 140
      NZ = 0
      IERR = 2
      RETURN
  140 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END



SUBROUTINE ZUCHK(YR, YI, NZ, ASCLE, TOL)
!***BEGIN PROLOGUE  ZUCHK
!***REFER TO ZSERI,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
!
!      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
!      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
!      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
!      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
!      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
!      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
!      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZUCHK
!
!     COMPLEX Y
      DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
      INTEGER NZ
      NZ = 0
      WR = DABS(YR)
      WI = DABS(YI)
      ST = DMIN1(WR,WI)
      IF (ST.GT.ASCLE) RETURN
      SS = DMAX1(WR,WI)
      ST = ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
END

SUBROUTINE ZBINU(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM)
USE COMPLEX
!***BEGIN PROLOGUE  ZBINU
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY

!   ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE

!***ROUTINES CALLED  ZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK
!***END PROLOGUE  ZBINU
      DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU, &
      FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CYR(1), CYI(1), CWR(2), CWI(2)
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /

      NZ = 0
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + DBLE(FLOAT(N-1))
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES
!-----------------------------------------------------------------------
      CALL ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      INW = IABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + DBLE(FLOAT(NN-1))
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0D0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!-----------------------------------------------------------------------
   30 CONTINUE
      CALL ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0D0) GO TO 70
   50 CONTINUE
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+DBLE(FLOAT(NN-1))
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!-----------------------------------------------------------------------
      CALL ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CYR(I) = ZEROR
        CYI(I) = ZEROI
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
      NUI = INT(SNGL(FNUL-DFNU)) + 1
      NUI = MAX0(NUI,0)
      CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL, &
      TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END

SUBROUTINE ZSERI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZSERI
!***REFER TO  ZBESI,ZBESK
!
!     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
!     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
!     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
!
!***ROUTINES CALLED  DGAMLN,D1MACH,ZUCHK,ZABS,ZDIV,ZLOG,ZMLT
!***END PROLOGUE  ZSERI
!     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
      DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL, &
      AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU, &
      ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,  &
      STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,  &
      ZR!, DGAMLN
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
      DIMENSION YR(1), YI(1), WR(2), WI(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /

      NZ = 0
      AZ = ZABS(ZR,ZI)
      IF (AZ.EQ.0.0D0) GO TO 160
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = DSQRT(ARM)
      CRSCR = 1.0D0
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 150
      HZR = 0.5D0*ZR
      HZI = 0.5D0*ZI
      CZR = ZEROR
      CZI = ZEROI
      IF (AZ.LE.RTR1) GO TO 10
      CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
      ACZ = ZABS(CZR,CZI)
      NN = N
      CALL ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
      DFNU = FNU + DBLE(FLOAT(NN-1))
      FNUP = DFNU + 1.0D0
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
      AK1R = CKR*DFNU
      AK1I = CKI*DFNU
      AK = DGAMLN(FNUP,IDUM)
      AK1R = AK1R - AK
      IF (KODE.EQ.2) AK1R = AK1R - ZR
      IF (AK1R.GT.(-ELIM)) GO TO 40
   30 CONTINUE
      NZ = NZ + 1
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      IF (ACZ.GT.DFNU) GO TO 190
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 20
   40 CONTINUE
      IF (AK1R.GT.(-ALIM)) GO TO 50
      IFLAG = 1
      SS = 1.0D0/TOL
      CRSCR = TOL
      ASCLE = ARM*SS
   50 CONTINUE
      AA = DEXP(AK1R)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEFR = AA*DCOS(AK1I)
      COEFI = AA*DSIN(AK1I)
      ATOL = TOL*ACZ/FNUP
      IL = MIN0(2,NN)
      DO 90 I=1,IL
        DFNU = FNU + DBLE(FLOAT(NN-I))
        FNUP = DFNU + 1.0D0
        S1R = CONER
        S1I = CONEI
        IF (ACZ.LT.TOL*FNUP) GO TO 70
        AK1R = CONER
        AK1I = CONEI
        AK = FNUP + 2.0D0
        S = FNUP
        AA = 2.0D0
   60   CONTINUE
        RS = 1.0D0/S
        STR = AK1R*CZR - AK1I*CZI
        STI = AK1R*CZI + AK1I*CZR
        AK1R = STR*RS
        AK1I = STI*RS
        S1R = S1R + AK1R
        S1I = S1I + AK1I
        S = S + AK
        AK = AK + 2.0D0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 60
   70   CONTINUE
        S2R = S1R*COEFR - S1I*COEFI
        S2I = S1R*COEFI + S1I*COEFR
        WR(I) = S2R
        WI(I) = S2I
        IF (IFLAG.EQ.0) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 30
   80   CONTINUE
        M = NN - I + 1
        YR(M) = S2R*CRSCR
        YI(M) = S2I*CRSCR
        IF (I.EQ.IL) GO TO 90
        CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
        COEFR = STR*DFNU
        COEFI = STI*DFNU
   90 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = DBLE(FLOAT(K))
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IF (IFLAG.EQ.1) GO TO 120
      IB = 3
  100 CONTINUE
      DO 110 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
  110 CONTINUE
      RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD WITH SCALED VALUES
!-----------------------------------------------------------------------
  120 CONTINUE
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
!     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
!-----------------------------------------------------------------------
      S1R = WR(1)
      S1I = WI(1)
      S2R = WR(2)
      S2I = WI(2)
      DO 130 L=3,NN
        CKR = S2R
        CKI = S2I
        S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
        S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
        S1R = CKR
        S1I = CKI
        CKR = S2R*CRSCR
        CKI = S2I*CRSCR
        YR(K) = CKR
        YI(K) = CKI
        AK = AK - 1.0D0
        K = K - 1
        IF (ZABS(CKR,CKI).GT.ASCLE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 100
  150 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160 CONTINUE
      YR(1) = ZEROR
      YI(1) = ZEROI
      IF (FNU.NE.0.0D0) GO TO 170
      YR(1) = CONER
      YI(1) = CONEI
  170 CONTINUE
      IF (N.EQ.1) RETURN
      DO 180 I=2,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  180 CONTINUE
      RETURN
!-----------------------------------------------------------------------
!     RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
!     THE CALCULATION IN CBINU WITH N=N-IABS(NZ)
!-----------------------------------------------------------------------
  190 CONTINUE
      NZ = -NZ
      RETURN
      END

SUBROUTINE ZASYI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZASYI
!***REFER TO  ZBESI,ZBESK

!     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
!     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
!     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
!
!***ROUTINES CALLED  D1MACH,ZABS,ZDIV,ZEXP,ZMLT,ZSQRT
!***END PROLOGUE  ZASYI
!     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
      DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,   &
      AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI, &
      CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,    &
      P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,  &
      S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION YR(1), YI(1)
      DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /

      NZ = 0
      AZ = ZABS(ZR,ZI)
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = DSQRT(ARM)
      IL = MIN0(2,N)
      DFNU = FNU + DBLE(FLOAT(N-IL))
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      AK1R = RTPI*STR*RAZ
      AK1I = RTPI*STI*RAZ
      CALL ZSQRT(AK1R, AK1I, AK1R, AK1I)
      CZR = ZR
      CZI = ZI
      IF (KODE.NE.2) GO TO 10
      CZR = ZEROR
      CZI = ZI
   10 CONTINUE
      IF (DABS(CZR).GT.ELIM) GO TO 100
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((DABS(CZR).GT.ALIM) .AND. (N.GT.2)) GO TO 20
      KODED = 0
      CALL ZEXP(CZR, CZI, STR, STI)
      CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20 CONTINUE
      FDN = 0.0D0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZR = ZR*8.0D0
      EZI = ZI*8.0D0
!-----------------------------------------------------------------------
!     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
!     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
!     EXPANSION FOR THE IMAGINARY PART.
!-----------------------------------------------------------------------
      AEZ = 8.0D0*AZ
      S = TOL/AEZ
      JL = INT(SNGL(RL+RL)) + 2
      P1R = ZEROR
      P1I = ZEROI
      IF (ZI.EQ.0.0D0) GO TO 30
!-----------------------------------------------------------------------
!     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!     SIGNIFICANCE WHEN FNU OR N IS LARGE
!-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*PI
      INU = INU + N - IL
      AK = -DSIN(ARG)
      BK = DCOS(ARG)
      IF (ZI.LT.0.0D0) BK = -BK
      P1R = AK
      P1I = BK
      IF (MOD(INU,2).EQ.0) GO TO 30
      P1R = -P1R
      P1I = -P1I
   30 CONTINUE
      DO 70 K=1,IL
        SQK = FDN - 1.0D0
        ATOL = S*DABS(SQK)
        SGN = 1.0D0
        CS1R = CONER
        CS1I = CONEI
        CS2R = CONER
        CS2I = CONEI
        CKR = CONER
        CKI = CONEI
        AK = 0.0D0
        AA = 1.0D0
        BB = AEZ
        DKR = EZR
        DKI = EZI
        DO 40 J=1,JL
          CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI)
          CKR = STR*SQK
          CKI = STI*SQK
          CS2R = CS2R + CKR
          CS2I = CS2I + CKI
          SGN = -SGN
          CS1R = CS1R + CKR*SGN
          CS1I = CS1I + CKI*SGN
          DKR = DKR + EZR
          DKI = DKI + EZI
          AA = AA*DABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0D0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 50
   40   CONTINUE
        GO TO 110
   50   CONTINUE
        S2R = CS1R
        S2I = CS1I
        IF (ZR+ZR.GE.ELIM) GO TO 60
        TZR = ZR + ZR
        TZI = ZI + ZI
        CALL ZEXP(-TZR, -TZI, STR, STI)
        CALL ZMLT(STR, STI, P1R, P1I, STR, STI)
        CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI)
        S2R = S2R + STR
        S2I = S2I + STI
   60   CONTINUE
        FDN = FDN + 8.0D0*DFNU + 4.0D0
        P1R = -P1R
        P1I = -P1I
        M = N - IL + K
        YR(M) = S2R*AK1R - S2I*AK1I
        YI(M) = S2R*AK1I + S2I*AK1R
   70 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = DBLE(FLOAT(K))
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IB = 3
      DO 80 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
   80 CONTINUE
      IF (KODED.EQ.0) RETURN
      CALL ZEXP(CZR, CZI, CKR, CKI)
      DO 90 I=1,NN
        STR = YR(I)*CKR - YI(I)*CKI
        YI(I) = YR(I)*CKI + YI(I)*CKR
        YR(I) = STR
   90 CONTINUE
      RETURN
  100 CONTINUE
      NZ = -1
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END

SUBROUTINE ZUOIK(ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZUOIK
!***REFER TO  ZBESI,ZBESK,ZBESH
!
!     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
!     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
!     (IN LOGARITHMI! FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
!     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
!     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
!     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
!     MULTIPLIERS (IN LOGARITHMI! FORM) IS MADE BASED ON ELIM. HERE
!     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
!     EXP(-ELIM)/TOL
!
!     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
!          =2 MEANS THE K SEQUENCE IS TESTED
!     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
!         =-1 MEANS AN OVERFLOW WOULD OCCUR
!     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
!             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
!     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
!     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
!             ANOTHER ROUTINE
!
!***ROUTINES CALLED  ZUCHK,ZUNHJ,ZUNIK,D1MACH,ZABS,ZLOG
!***END PROLOGUE  ZUOIK
!     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
!    *ZR
      DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR, &
      ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,  &
      FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,   &
      YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,  &
      ZNI, ZNR, ZR, ZRI, ZRR
      INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION YR(1), YI(1), CWRKR(16), CWRKI(16)
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
      DATA AIC / 1.265512123484645396D+00 /
      NUF = 0
      NN = N
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      ZBR = ZRR
      ZBI = ZRI
      AX = DABS(ZR)*1.7321D0
      AY = DABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = DMAX1(FNU,1.0D0)
      IF (IKFLG.EQ.1) GO TO 20
      FNN = DBLE(FLOAT(NN))
      GNN = FNU + FNN - 1.0D0
      GNU = DMAX1(GNN,FNN)
   20 CONTINUE
!-----------------------------------------------------------------------
!     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
!     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
!     THE SIGN OF THE IMAGINARY PART CORRECT.
!-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 30
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII, &
      ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 50
   30 CONTINUE
      ZNR = ZRI
      ZNI = -ZRR
      IF (ZI.GT.0.0D0) GO TO 40
      ZNR = -ZNR
   40 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, &
      ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
   50 CONTINUE
      IF (KODE.EQ.1) GO TO 60
      CZR = CZR - ZBR
      CZI = CZI - ZBI
   60 CONTINUE
      IF (IKFLG.EQ.1) GO TO 70
      CZR = -CZR
      CZI = -CZI
   70 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 210
      IF (RCZ.LT.ALIM) GO TO 80
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 210
      GO TO 130
   80 CONTINUE
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 90
      IF (RCZ.GT.(-ALIM)) GO TO 130
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 110
   90 CONTINUE
      DO 100 I=1,NN
        YR(I) = ZEROR
        YI(I) = ZEROI
  100 CONTINUE
      NUF = NN
      RETURN
  110 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 120
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  120 CONTINUE
      AX = DEXP(RCZ)/TOL
      AY = CZI
      CZR = AX*DCOS(AY)
      CZI = AX*DSIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 90
  130 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOWS ON I SEQUENCE
!-----------------------------------------------------------------------
  140 CONTINUE
      GNU = FNU + DBLE(FLOAT(NN-1))
      IF (IFORM.EQ.2) GO TO 150
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII, &
      ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 160
  150 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, &
      ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
  160 CONTINUE
      IF (KODE.EQ.1) GO TO 170
      CZR = CZR - ZBR
      CZI = CZI - ZBI
  170 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
      IF (RCZ.LT.(-ELIM)) GO TO 180
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + DLOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*DLOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 190
  180 CONTINUE
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 140
  190 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 200
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  200 CONTINUE
      AX = DEXP(RCZ)/TOL
      AY = CZI
      CZR = AX*DCOS(AY)
      CZI = AX*DSIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 180
      RETURN
  210 CONTINUE
      NUF = -1
      RETURN
      END

SUBROUTINE ZBUNI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST, FNUL, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZBUNI
!***REFER TO  ZBESI,ZBESK
!
!     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
!     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
!
!***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,D1MACH
!***END PROLOGUE  ZBUNI
!     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
      DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,    &
      ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,  &
      S2I, S2R, TOL, YI, YR, ZI, ZR, ASCLE, BRY, C1R, C1I, C1M
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION YR(1), YI(1), CYR(2), CYI(2), BRY(3)
      NZ = 0
      AX = DABS(ZR)*1.7321D0
      AY = DABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = DBLE(FLOAT(NUI))
      DFNU = FNU + DBLE(FLOAT(N-1))
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3.LE.ARG(Z).LE.PI/3
!-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      STR = ZABS(CYR(1),CYI(1))
!----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!----------------------------------------------------------------------
      BRY(1)=1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      CSCLR = 1.0D0
      IF (STR.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      CSCLR = 1.0D0/TOL
      GO TO 25
   21 CONTINUE
      IF (STR.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE=BRY(3)
      CSCLR = TOL
   25 CONTINUE
      CSCRR = 1.0D0/CSCLR
      S1R = CYR(2)*CSCLR
      S1I = CYI(2)*CSCLR
      S2R = CYR(1)*CSCLR
      S2I = CYI(1)*CSCLR
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      DO 30 I=1,NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0D0
        IF (IFLAG.GE.3) GO TO 30
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   30 CONTINUE
      YR(N) = S2R*CSCRR
      YI(N) = S2I*CSCRR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = DBLE(FLOAT(NL))
      K = NL
      DO 40 I=1,NL
        STR = S2R
        STI = S2I
        S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        YR(K) = STR
        YI(K) = STI
        FNUI = FNUI - 1.0D0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        C1R = DABS(STR)
        C1I = DABS(STI)
        C1M = DMAX1(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3.LE.ARG(Z).LE.PI/3
!-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
      GO TO 80
   70 CONTINUE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
!     AND HPI=PI/2
!-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END

SUBROUTINE ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZUNI1
!***REFER TO  ZBESI,ZBESK
!
!     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***ROUTINES CALLED  ZUCHK,ZUNIK,ZUOIK,D1MACH,ZABS
!***END PROLOGUE  ZUNI1
!     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
!    *S2,Y,Z,ZETA1,ZETA2
      DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONEI, CONER, CRSC,  &
      CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN, &
      FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,   &
      SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,  &
      ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), YR(1), YI(1), CWRKR(16), CWRKI(16), CSSR(3), &
      CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /

      NZ = 0
      ND = N
      NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0D0)
      INIT = 0
      CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R, &
      ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      IF (KODE.EQ.1) GO TO 10
      STR = ZR + ZETA2R
      STI = ZI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 20
   10 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   20 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN0(2,ND)
      DO 80 I=1,NN
        FN = FNU + DBLE(FLOAT(ND-I))
        INIT = 0
        CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R, &
        ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
        IF (KODE.EQ.1) GO TO 40
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        GO TO 50
   40   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   50   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 60
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        RS1 = RS1 + DLOG(APHI)
        IF (DABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 IF CABS(S1).LT.ASCLE
!-----------------------------------------------------------------------
        S2R = PHIR*SUMR - PHII*SUMI
        S2I = PHIR*SUMI + PHII*SUMR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 70
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        CYR(I) = S2R
        CYI(I) = S2I
        M = ND - I + 1
        YR(M) = S2R*CSRR(IFLAG)
        YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RAST = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAST
      STI = -ZI*RAST
      RZR = (STR+STR)*RAST
      RZI = (STI+STI)*RAST
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = DBLE(FLOAT(K))
      DO 90 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 90
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + DBLE(FLOAT(ND-1))
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      NZ = N
      DO 140 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  140 CONTINUE
      RETURN
      END

SUBROUTINE ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZUNI2
!***REFER TO  ZBESI,ZBESK
!
!     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!***ROUTINES CALLED  ZAIRY,ZUCHK,ZUNHJ,ZUOIK,D1MACH,ZABS
!***END PROLOGUE  ZUNI2
!     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
!    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,     &
      ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,  &
      CONEI, CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,  &
      DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI, &
      RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI, &
      ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR, CYI
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,  &
      NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), YR(1), YI(1), CIPR(4), CIPI(4), CSSR(3),       &
      CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),    &
      CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
      DATA HPI, AIC  /1.57079632679489662D+00,  1.265512123484645396D+00/

      NZ = 0
      ND = N
      NLAST = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      ZBR = ZR
      ZBI = ZI
      CIDI = -CONER
      INU = INT(SNGL(FNU))
      ANG = HPI*(FNU-DBLE(FLOAT(INU)))
      C2R = DCOS(ANG)
      C2I = DSIN(ANG)
      IN = INU + N - 1
      IN = MOD(IN,4) + 1
      STR = C2R*CIPR(IN) - C2I*CIPI(IN)
      C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
      C2R = STR
      IF (ZI.GT.0.0D0) GO TO 10
      ZNR = -ZNR
      ZBI = -ZBI
      CIDI = -CIDI
      C2I = -C2I
   10 CONTINUE
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0D0)
      CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R, &
      ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      IF (KODE.EQ.1) GO TO 20
      STR = ZBR + ZETA2R
      STI = ZBI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 30
   20 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   30 CONTINUE
      RS1 = S1R
      IF (DABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN0(2,ND)
      DO 90 I=1,NN
        FN = FNU + DBLE(FLOAT(ND-I))
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI, &
        ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
        IF (KODE.EQ.1) GO TO 50
        STR = ZBR + ZETA2R
        STI = ZBI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + DABS(ZI)
        GO TO 60
   50   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   60   CONTINUE
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
        RS1 = S1R
        IF (DABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (DABS(RS1).LT.ALIM) GO TO 70
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        AARG = ZABS(ARGR,ARGI)
        RS1 = RS1 + DLOG(APHI) - 0.25D0*DLOG(AARG) - AIC
        IF (DABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
        CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR - DAII*BSUMI
        STI = DAIR*BSUMI + DAII*BSUMR
        STR = STR + (AIR*ASUMR-AII*ASUMI)
        STI = STI + (AIR*ASUMI+AII*ASUMR)
        S2R = PHIR*STR - PHII*STI
        S2I = PHIR*STI + PHII*STR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (ZI.LE.0.0D0) S2I = -S2I
        STR = S2R*C2R - S2I*C2I
        S2I = S2R*C2I + S2I*C2R
        S2R = STR
        CYR(I) = S2R
        CYI(I) = S2I
        J = ND - I + 1
        YR(J) = S2R*CSRR(IFLAG)
        YI(J) = S2I*CSRR(IFLAG)
        STR = -C2I*CIDI
        C2I = C2R*CIDI
        C2R = STR
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = DBLE(FLOAT(K))
      DO 100 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 100
        STR = DABS(C2R)
        STI = DABS(C2I)
        C2M = DMAX1(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + DBLE(FLOAT(ND-1))
      IF (FN.LT.FNUL) GO TO 130
      FN = CIDI
      J = NUF + 1
      K = MOD(J,4) + 1
      S1R = CIPR(K)
      S1I = CIPI(K)
      IF (FN.LT.0.0D0) S1I = -S1I
      STR = C2R*S1R - C2I*S1I
      C2I = C2R*S1I + C2I*S1R
      C2R = STR
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
      NZ = N
      DO 160 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  160 CONTINUE
      RETURN
      END

SUBROUTINE ZWRSK(ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZWRSK
!***REFER TO  ZBESI,ZBESK
!
!     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
!
!***ROUTINES CALLED  D1MACH,ZBKNU,ZRATI,ZABS
!***END PROLOGUE  ZWRSK
!     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
      DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI, &
      CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,     &
      STI, STR, TOL, YI, YR, ZRI, ZRR
      INTEGER I, KODE, N, NW, NZ
      DIMENSION YR(1), YI(1), CWR(2), CWI(2)
!-----------------------------------------------------------------------
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!-----------------------------------------------------------------------
      NZ = 0
      CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
!-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z)=Y(J),  J=1,...,N
!-----------------------------------------------------------------------
      CINUR = 1.0D0
      CINUI = 0.0D0
      IF (KODE.EQ.1) GO TO 10
      CINUR = DCOS(ZRI)
      CINUI = DSIN(ZRI)
   10 CONTINUE
!-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
!     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
!     THE RESULT IS ON SCALE.
!-----------------------------------------------------------------------
      ACW = ZABS(CWR(2),CWI(2))
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CSCLR = 1.0D0
      IF (ACW.GT.ASCLE) GO TO 20
      CSCLR = 1.0D0/TOL
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0D0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCLR = TOL
   30 CONTINUE
      C1R = CWR(1)*CSCLR
      C1I = CWI(1)*CSCLR
      C2R = CWR(2)*CSCLR
      C2I = CWI(2)*CSCLR
      STR = YR(1)
      STI = YI(1)
!-----------------------------------------------------------------------
!     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
!-----------------------------------------------------------------------
      PTR = STR*C1R - STI*C1I
      PTI = STR*C1I + STI*C1R
      PTR = PTR + C2R
      PTI = PTI + C2I
      CTR = ZRR*PTR - ZRI*PTI
      CTI = ZRR*PTI + ZRI*PTR
      ACT = ZABS(CTR,CTI)
      RACT = 1.0D0/ACT
      CTR = CTR*RACT
      CTI = -CTI*RACT
      PTR = CINUR*RACT
      PTI = CINUI*RACT
      CINUR = PTR*CTR - PTI*CTI
      CINUI = PTR*CTI + PTI*CTR
      YR(1) = CINUR*CSCLR
      YI(1) = CINUI*CSCLR
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        PTR = STR*CINUR - STI*CINUI
        CINUI = STR*CINUI + STI*CINUR
        CINUR = PTR
        STR = YR(I)
        STI = YI(I)
        YR(I) = CINUR*CSCLR
        YI(I) = CINUI*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
END

SUBROUTINE ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZBKNU
!***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
!
!     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
!
!***ROUTINES CALLED  DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,ZABS,ZDIV,
!                    ZEXP,ZLOG,ZMLT,ZSQRT
!***END PROLOGUE  ZBKNU
!
      DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ, &
      CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER, &
      CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOI, CTWOR,   &
      CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS, &
      FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI,&
      PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,    &
      RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,  &
      TOL, TTH, T1, T2, YI, YR, ZI, ZR, ELM, CELMR, ZDR, ZDI, &
	  AS, ALAS, HELIM, CYR, CYI!, DGAMLN
      INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,  &
      IDUM, J, IC, INUB, NW
      DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),&
      CYI(2)
!     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
!     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK

      DATA KMAX / 30 /
      DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,CTWOI,R1/ &
        0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 0.0D0 , 2.0D0 /
      DATA DPI, RTHPI, SPI ,HPI, FPI, TTH /  &
           3.14159265358979324D0,       1.25331413731550025D0, &
           1.90985931710274403D0,       1.57079632679489662D0, &
           1.89769999331517738D0,       6.66666666666666666D-01/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/ &
           5.77215664901532861D-01,    -4.20026350340952355D-02,   &
          -4.21977345555443367D-02,     7.21894324666309954D-03,   &
          -2.15241674114950973D-04,    -2.01348547807882387D-05,   &
           1.13302723198169588D-06,     6.11609510448141582D-09/

      CAZ = ZABS(ZR,ZI)
      CSCLR = 1.0D0/TOL
      CRSCR = TOL
      CSSR(1) = CSCLR
      CSSR(2) = 1.0D0
      CSSR(3) = CRSCR
      CSRR(1) = CRSCR
      CSRR(2) = 1.0D0
      CSRR(3) = CSCLR
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RCAZ = 1.0D0/CAZ
      STR = ZR*RCAZ
      STI = -ZI*RCAZ
      RZR = (STR+STR)*RCAZ
      RZI = (STI+STI)*RCAZ
      INU = INT(SNGL(FNU+0.5D0))
      DNU = FNU - DBLE(FLOAT(INU))
      IF (DABS(DNU).EQ.0.5D0) GO TO 110
      DNU2 = 0.0D0
      IF (DABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
!-----------------------------------------------------------------------
!     SERIES FOR CABS(Z).LE.R1
!-----------------------------------------------------------------------
      FC = 1.0D0
      CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
      FMUR = SMUR*DNU
      FMUI = SMUI*DNU
      CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
      IF (DNU.EQ.0.0D0) GO TO 10
      FC = DNU*DPI
      FC = FC/DSIN(FC)
      SMUR = CSHR/DNU
      SMUI = CSHI/DNU
   10 CONTINUE
      A2 = 1.0D0 + DNU
!-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
!-----------------------------------------------------------------------
      T2 = DEXP(-DGAMLN(A2,IDUM))
      T1 = 1.0D0/(T2*FC)
      IF (DABS(DNU).GT.0.1D0) GO TO 40
!-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!-----------------------------------------------------------------------
      AK = 1.0D0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (DABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      FR = FC*(CCHR*G1+SMUR*G2)
      FI = FC*(CCHI*G1+SMUI*G2)
      CALL ZEXP(FMUR, FMUI, STR, STI)
      PR = 0.5D0*STR/T2
      PI = 0.5D0*STI/T2
      CALL ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
      QR = PTR/T1
      QI = PTI/T1
      S1R = FR
      S1I = FI
      S2R = PR
      S2I = PI
      AK = 1.0D0
      A1 = 1.0D0
      CKR = CONER
      CKI = CONEI
      BK = 1.0D0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
!-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
!-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   60 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF (KODED.EQ.1) RETURN
      CALL ZEXP(ZR, ZI, STR, STI)
      CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
      RETURN
!-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   90 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      STR = PR - FR*AK
      STI = PI - FI*AK
      S2R = CKR*STR - CKI*STI + S2R
      S2I = CKR*STI + CKI*STR + S2I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      A1 = FNU + 1.0D0
      AK = A1*DABS(SMUR)
      IF (AK.GT.ALIM) KFLAG = 3
      STR = CSSR(KFLAG)
      P2R = S2R*STR
      P2I = S2I*STR
      CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
      S1R = S1R*STR
      S1I = S1I*STR
      IF (KODED.EQ.1) GO TO 210
      CALL ZEXP(ZR, ZI, FR, FI)
      CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
      CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
      GO TO 210
!-----------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!-----------------------------------------------------------------------
  110 CONTINUE
      CALL ZSQRT(ZR, ZI, STR, STI)
      CALL ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (ZR.GT.ALIM) GO TO 290
!     BLANK LINE
      STR = DEXP(-ZR)*CSSR(KFLAG)
      STI = -STR*DSIN(ZI)
      STR = STR*DCOS(ZI)
      CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
      IF (DABS(DNU).EQ.0.5D0) GO TO 300
!-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR CABS(Z).GT.R1
!-----------------------------------------------------------------------
      AK = DCOS(DPI*DNU)
      AK = DABS(AK)
      IF (AK.EQ.CZEROR) GO TO 300
      FHS = DABS(0.25D0-DNU2)
      IF (FHS.EQ.CZEROR) GO TO 300
!-----------------------------------------------------------------------
!     COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
!     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!-----------------------------------------------------------------------
      T1 = DBLE(FLOAT(I1MACH(14)-1))
      T1 = T1*D1MACH(5)*3.321928094D0
      T1 = DMAX1(T1,12.0D0)
      T1 = DMIN1(T1,60.0D0)
      T2 = TTH*T1 - 6.0D0
      IF (ZR.NE.0.0D0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = DATAN(ZI/ZR)
      T1 = DABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
!-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
!-----------------------------------------------------------------------
      ETEST = AK/(DPI*CAZ*TOL)
      FK = CONER
      IF (ETEST.LT.CONER) GO TO 180
      FKS = CTWOR
      CKR = CAZ + CAZ + CTWOR
      P1R = CZEROR
      P2R = CONER
      DO 150 I=1,KMAX
        AK = FHS/FKS
        CBR = CKR/(FK+CONER)
        PTR = P2R
        P2R = CBR*P2R - P1R*AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = DABS(P2R)*FK
        IF (ETEST.LT.STR) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*DSQRT(T2/CAZ)
      FHS = DABS(0.25D0-DNU2)
      GO TO 180
  170 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
!-----------------------------------------------------------------------
      A2 = DSQRT(CAZ)
      AK = FPI*AK/(TOL*DSQRT(A2))
      AA = 3.0D0*T1/(1.0D0+CAZ)
      BB = 14.7D0*T1/(28.0D0+CAZ)
      AK = (DLOG(AK)+CAZ*DCOS(AA)/(1.0D0+0.008D0*CAZ))/DCOS(BB)
      FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
      K = INT(SNGL(FK))
      FK = DBLE(FLOAT(K))
      FKS = FK*FK
      P1R = CZEROR
      P1I = CZEROI
      P2R = TOL
      P2I = CZEROI
      CSR = P2R
      CSI = P2I
      DO 190 I=1,K
        A1 = FKS - FK
        AK = (FKS+FK)/(A1+FHS)
        RAK = 2.0D0/(FK+CONER)
        CBR = (FK+ZR)*RAK
        CBI = ZI*RAK
        PTR = P2R
        PTI = P2I
        P2R = (PTR*CBR-PTI*CBI-P1R)*AK
        P2I = (PTI*CBR+PTR*CBI-P1I)*AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = A1 - FK + CONER
        FK = FK - CONER
  190 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
!     SCALING
!-----------------------------------------------------------------------
      TM = ZABS(CSR,CSI)
      PTR = 1.0D0/TM
      S1R = P2R*PTR
      S1I = P2I*PTR
      CSR = CSR*PTR
      CSI = -CSI*PTR
      CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
      CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
!-----------------------------------------------------------------------
      TM = ZABS(P2R,P2I)
      PTR = 1.0D0/TM
      P1R = P1R*PTR
      P1I = P1I*PTR
      P2R = P2R*PTR
      P2I = -P2I*PTR
      CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
      STR = DNU + 0.5D0 - PTR
      STI = -PTI
      CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
      STR = STR + 1.0D0
      CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
!-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!-----------------------------------------------------------------------
  210 CONTINUE
      STR = DNU + 1.0D0
      CKR = STR*RZR
      CKI = STR*RZI
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.GT.1) GO TO 215
      S1R = S2R
      S1I = S2I
  215 CONTINUE
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF(IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        IF (KFLAG.GE.3) GO TO 230
        P2R = S2R*P1R
        P2I = S2I*P1R
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  230 CONTINUE
      IF (N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
  240 CONTINUE
      STR = CSRR(KFLAG)
      YR(1) = S1R*STR
      YI(1) = S1I*STR
      IF (N.EQ.1) RETURN
      YR(2) = S2R*STR
      YI(2) = S2I*STR
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2R = S2R
        P2I = S2I
        S2R = CKR*P2R - CKI*P2I + S1R
        S2I = CKI*P2R + CKR*P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R*P1R
        P2I = S2I*P1R
        YR(I) = P2R
        YI(I) = P2I
        IF (KFLAG.GE.3) GO TO 260
        STR = DABS(P2R)
        STI = DABS(P2I)
        P2M = DMAX1(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  260 CONTINUE
      RETURN
!-----------------------------------------------------------------------
!     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5D0*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ASCLE = BRY(1)
      ZDR = ZR
      ZDI = ZI
      IC = -1
      J = 2
      DO 262 I=1,INU
        STR = S2R
        STI = S2I
        S2R = STR*CKR-STI*CKI+S1R
        S2I = STI*CKR+STR*CKI+S1I
        S1R = STR
        S1I = STI
        CKR = CKR+RZR
        CKI = CKI+RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        P2R = -ZDR+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        CALL ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR+STR
        P2I = -ZDI+STI
        P2M = DEXP(P2R)/TOL
        P1R = P2M*DCOS(P2I)
        P1I = P2M*DSIN(P2I)
        CALL ZUCHK(P1R,P1I,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J = 3 - J
        CYR(J) = P1R
        CYI(J) = P1I
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        ZDR = ZDR-ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
  262 CONTINUE
      IF(N.NE.1) GO TO 270
      S1R = S2R
      S1I = S2I
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2R = CYR(J)
      S2I = CYI(J)
      J = 3 - J
      S1R = CYR(J)
      S1I = CYI(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
      GO TO 240
  270 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF(N.EQ.1) GO TO 280
      YR(2) = S2R
      YI(2) = S2I
  280 CONTINUE
      ASCLE = BRY(1)
      CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1R = YR(KK)
      S1I = YI(KK)
      YR(KK) = S1R*CSRR(1)
      YI(KK) = S1I*CSRR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2R = YR(KK)
      S2I = YI(KK)
      YR(KK) = S2R*CSRR(1)
      YI(KK) = S2I*CSRR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + DBLE(FLOAT(KK-1))
      CKR = T2*RZR
      CKI = T2*RZI
      KFLAG = 1
      GO TO 250
  290 CONTINUE
!-----------------------------------------------------------------------
!     SCALE BY DEXP(Z), IFLAG = 1 CASES
!-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
!-----------------------------------------------------------------------
!     FNU=HALF ODD INTEGER CASE, DNU=-0.5
!-----------------------------------------------------------------------
  300 CONTINUE
      S1R = COEFR
      S1I = COEFI
      S2R = COEFR
      S2I = COEFI
      GO TO 210

  310 CONTINUE
      NZ=-2
      RETURN
END

SUBROUTINE ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
USE COMPLEX
!***BEGIN PROLOGUE  ZKSCL
!***REFER TO  ZBESK
!
!     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
!     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
!
!***ROUTINES CALLED  ZUCHK,ZABS,ZLOG
!***END PROLOGUE  ZKSCL
!     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI, &
      CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I, S2R,    &
      TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZDR, ZDI, CELMR,     &
	  ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(1), YI(1), CYR(2), CYI(2)
      DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /

      NZ = 0
      IC = 0
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZABS(S1R,S1I)
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 10
        CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5D0*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI

!     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
!     S2 GETS LARGER THAN EXP(ELIM/2)

      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 25
        CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        IF (IC.EQ.KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
END

SUBROUTINE ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)
!***BEGIN PROLOGUE  ZSHCH
!***REFER TO  ZBESK,ZBESH
!
!     ZSHCH COMPUTES THE COMPLEX HYPERBOLI! FUNCTIONS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZSHCH
!
  DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR, DCOSH, DSINH
      SH = DSINH(ZR)
      CH = DCOSH(ZR)
      SN = DSIN(ZI)
      CN = DCOS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
  RETURN
END

SUBROUTINE ZMLRI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZMLRI
!***REFER TO  ZBESI,ZBESK
!
!     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
!
!***ROUTINES CALLED  DGAMLN,D1MACH,ZABS,ZEXP,ZLOG,ZMLT
!***END PROLOGUE  ZMLRI
!     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
      DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,     &
      CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I, &
      P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,  &
      SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR!, DGAMLN
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
      DIMENSION YR(1), YI(1)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      SCLE = D1MACH(1)/TOL
      NZ=0
      AZ = ZABS(ZR,ZI)
      IAZ = INT(SNGL(AZ))
      IFNU = INT(SNGL(FNU))
      INU = IFNU + N - 1
      AT = DBLE(FLOAT(IAZ)) + 1.0D0
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      ACK = (AT+1.0D0)*RAZ
      RHO = ACK + DSQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKI*PTR+CKR*PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0D0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!-----------------------------------------------------------------------
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      AT = DBLE(FLOAT(INU)) + 1.0D0
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      ACK = AT*RAZ
      TST = DSQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKR*PTI+CKI*PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = ZABS(CKR,CKI)
        FLAM = ACK + DSQRT(ACK*ACK-1.0D0)
        FKAP = AP/ZABS(P1R,P1I)
        RHO = DMIN1(FLAM,FKAP)
        TST = TST*DSQRT(RHO/(RHO*RHO-1.0D0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
!-----------------------------------------------------------------------
      K = K + 1
      KK = MAX0(I+IAZ,K+INU)
      FKK = DBLE(FLOAT(KK))
      P1R = ZEROR
      P1I = ZEROI
!-----------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!-----------------------------------------------------------------------
      P2R = SCLE
      P2I = ZEROI
      FNF = FNU - DBLE(FLOAT(IFNU))
      TFNF = FNF + FNF
      BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) - &
      DGAMLN(TFNF+1.0D0,IDUM)
      BK = DEXP(BK)
      SUMR = ZEROR
      SUMI = ZEROI
      KM = KK - INU
      DO 50 I=1,KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   50 CONTINUE
      YR(N) = P2R
      YI(N) = P2I
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
        M = N - I + 1
        YR(M) = P2R
        YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
      PTR = ZR
      PTI = ZI
      IF (KODE.EQ.2) PTR = ZEROR
      CALL ZLOG(RZR, RZI, STR, STI, IDUM)
      P1R = -FNF*STR + PTR
      P1I = -FNF*STI + PTI
      AP = DGAMLN(1.0D0+FNF,IDUM)
      PTR = P1R - AP
      PTI = P1I
!-----------------------------------------------------------------------
!     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!-----------------------------------------------------------------------
      P2R = P2R + SUMR
      P2I = P2I + SUMI
      AP = ZABS(P2R,P2I)
      P1R = 1.0D0/AP
      CALL ZEXP(PTR, PTI, STR, STI)
      CKR = STR*P1R
      CKI = STI*P1R
      PTR = P2R*P1R
      PTI = -P2I*P1R
      CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
      DO 100 I=1,N
        STR = YR(I)*CNORMR - YI(I)*CNORMI
        YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
        YR(I) = STR
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
END

DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
USE UTILIT
!***BEGIN PROLOGUE  DGAMLN
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5F
!***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
!***DESCRIPTION
!
!               **** A DOUBLE PRECISION ROUTINE ****
!         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
!         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
!         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
!         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
!         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
!         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
!         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
!
!         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
!         VALUES IS USED FOR SPEED OF EXECUTION.
!
!     DESCRIPTION OF ARGUMENTS
!
!         INPUT      Z IS D0UBLE PRECISION
!           Z      - ARGUMENT, Z.GT.0.0D0
!
!         OUTPUT      DGAMLN IS DOUBLE PRECISION
!           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
!           IERR    - ERROR FLAG
!                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
!                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
!
!
!***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
!***ROUTINES CALLED  I1MACH,D1MACH
!***END PROLOGUE  DGAMLN
      DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST, &
	  T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ
      INTEGER I, IERR, I1M, K, MZ, NZ
      DIMENSION CF(22), GLN(100)
!       LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),       &
           GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),  &
           GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),         &
           GLN(21), GLN(22)/                                             &
           0.00000000000000000D+00,     0.00000000000000000D+00,         &
           6.93147180559945309D-01,     1.79175946922805500D+00,         &
           3.17805383034794562D+00,     4.78749174278204599D+00,         &
           6.57925121201010100D+00,     8.52516136106541430D+00,         &
           1.06046029027452502D+01,     1.28018274800814696D+01,         &
           1.51044125730755153D+01,     1.75023078458738858D+01,         &
           1.99872144956618861D+01,     2.25521638531234229D+01,         &
           2.51912211827386815D+01,     2.78992713838408916D+01,         &
           3.06718601060806728D+01,     3.35050734501368889D+01,         &
           3.63954452080330536D+01,     3.93398841871994940D+01,         &
           4.23356164607534850D+01,     4.53801388984769080D+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),         &
           GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),         &
           GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),         &
           GLN(41), GLN(42), GLN(43), GLN(44)/                           &
           4.84711813518352239D+01,     5.16066755677643736D+01,         &
           5.47847293981123192D+01,     5.80036052229805199D+01,         &
           6.12617017610020020D+01,     6.45575386270063311D+01,         &
           6.78897431371815350D+01,     7.12570389671680090D+01,         &
           7.46582363488301644D+01,     7.80922235533153106D+01,         &
           8.15579594561150372D+01,     8.50544670175815174D+01,         &
           8.85808275421976788D+01,     9.21361756036870925D+01,         &
           9.57196945421432025D+01,     9.93306124547874269D+01,         &
           1.02968198614513813D+02,     1.06631760260643459D+02,         &
           1.10320639714757395D+02,     1.14034211781461703D+02,         &
           1.17771881399745072D+02,     1.21533081515438634D+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),         &
           GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),         &
           GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),         &
           GLN(63), GLN(64), GLN(65), GLN(66)/                           &
           1.25317271149356895D+02,     1.29123933639127215D+02,         &
           1.32952575035616310D+02,     1.36802722637326368D+02,         &
           1.40673923648234259D+02,     1.44565743946344886D+02,         &
           1.48477766951773032D+02,     1.52409592584497358D+02,         &
           1.56360836303078785D+02,     1.60331128216630907D+02,         &
           1.64320112263195181D+02,     1.68327445448427652D+02,         &
           1.72352797139162802D+02,     1.76395848406997352D+02,         &
           1.80456291417543771D+02,     1.84533828861449491D+02,         &
           1.88628173423671591D+02,     1.92739047287844902D+02,         &
           1.96866181672889994D+02,     2.01009316399281527D+02,         &
           2.05168199482641199D+02,     2.09342586752536836D+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),         &
           GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),         &
           GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),         &
           GLN(85), GLN(86), GLN(87), GLN(88)/                           &
           2.13532241494563261D+02,     2.17736934113954227D+02,         &
           2.21956441819130334D+02,     2.26190548323727593D+02,         &
           2.30439043565776952D+02,     2.34701723442818268D+02,         &
           2.38978389561834323D+02,     2.43268849002982714D+02,         &
           2.47572914096186884D+02,     2.51890402209723194D+02,         &
           2.56221135550009525D+02,     2.60564940971863209D+02,         &
           2.64921649798552801D+02,     2.69291097651019823D+02,         &
           2.73673124285693704D+02,     2.78067573440366143D+02,         &
           2.82474292687630396D+02,     2.86893133295426994D+02,         &
           2.91323950094270308D+02,     2.95766601350760624D+02,         &
           3.00220948647014132D+02,     3.04686856765668715D+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),         &
           GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/        &
           3.09164193580146922D+02,     3.13652829949879062D+02,         &
           3.18152639620209327D+02,     3.22663499126726177D+02,         &
           3.27185287703775217D+02,     3.31717887196928473D+02,         &
           3.36261181979198477D+02,     3.40815058870799018D+02,         &
           3.45379407062266854D+02,     3.49954118040770237D+02,         &
           3.54539085519440809D+02,     3.59134205369575399D+02/
!             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),       &
           CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),        &
           CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/       &
           8.33333333333333333D-02,    -2.77777777777777778D-03,         &
           7.93650793650793651D-04,    -5.95238095238095238D-04,         &
           8.41750841750841751D-04,    -1.91752691752691753D-03,         &
           6.41025641025641026D-03,    -2.95506535947712418D-02,         &
           1.79644372368830573D-01,    -1.39243221690590112D+00,         &
           1.34028640441683920D+01,    -1.56848284626002017D+02,         &
           2.19310333333333333D+03,    -3.61087712537249894D+04,         &
           6.91472268851313067D+05,    -1.52382215394074162D+07,         &
           3.82900751391414141D+08,    -1.08822660357843911D+10,         &
           3.47320283765002252D+11,    -1.23696021422692745D+13,         &
           4.88788064793079335D+14,    -2.13203339609193739D+16/

!             LN(2*PI)
      DATA CON                    /     1.83787706640934548D+00/

!***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0D0) GO TO 70
      IF (Z.GT.101.0D0) GO TO 10
      NZ = INT(SNGL(Z))
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0D0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = DMAX1(WDTOL,0.5D-18)
      I1M = I1MACH(14)
      RLN = D1MACH(5)*FLOAT(I1M)
      FLN = DMIN1(RLN,20.0D0)
      FLN = DMAX1(FLN,3.0D0)
      FLN = FLN - 3.0D0
      ZM = 1.8000D0 + 0.3875D0*FLN
      MZ = INT(SNGL(ZM)) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0D0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0D0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (DABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0D0) GO TO 50
      TLG = DLOG(Z)
      DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0D0
      NZ = INT(SNGL(ZINC))
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = DLOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0D0) - DLOG(ZP) + 0.5D0*(CON-TLG) + S
      RETURN

   70 CONTINUE
      IERR=1
      RETURN
END

SUBROUTINE ZUNIK(ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,  &
       PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
USE COMPLEX
!***BEGIN PROLOGUE  ZUNIK
!***REFER TO  ZBESI,ZBESK
!
!        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
!        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
!        RESPECTIVELY BY
!
!        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
!
!        WHERE       ZETA=-ZETA1 + ZETA2       OR
!                          ZETA1 - ZETA2
!
!        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
!        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
!        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
!        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
!        ZETA1,ZETA2.
!
!***ROUTINES CALLED  ZDIV,ZLOG,ZSQRT
!***END PROLOGUE  ZUNIK
!     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
!    *ZETA2,ZN,ZR
      DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,  &
      CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,   &
      SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R, &
      ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR
      INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      DATA CON(1), CON(2)  /  &
       3.98942280401432678D-01,  1.25331413731550025D+00 /
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
           C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),      &
           C(19), C(20), C(21), C(22), C(23), C(24)/                    &
           1.00000000000000000D+00,    -2.08333333333333333D-01,        &
           1.25000000000000000D-01,     3.34201388888888889D-01,        &
          -4.01041666666666667D-01,     7.03125000000000000D-02,        &
          -1.02581259645061728D+00,     1.84646267361111111D+00,        &
          -8.91210937500000000D-01,     7.32421875000000000D-02,        &
           4.66958442342624743D+00,    -1.12070026162229938D+01,        &
           8.78912353515625000D+00,    -2.36408691406250000D+00,        &
           1.12152099609375000D-01,    -2.82120725582002449D+01,        &
           8.46362176746007346D+01,    -9.18182415432400174D+01,        &
           4.25349987453884549D+01,    -7.36879435947963170D+00,        &
           2.27108001708984375D-01,     2.12570130039217123D+02,        &
          -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),      &
           C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),      &
           C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/      &
          -6.99579627376132541D+02,     2.18190511744211590D+02,        &
          -2.64914304869515555D+01,     5.72501420974731445D-01,        &
          -1.91945766231840700D+03,     8.06172218173730938D+03,        &
          -1.35865500064341374D+04,     1.16553933368645332D+04,        &
          -5.30564697861340311D+03,     1.20090291321635246D+03,        &
          -1.08090919788394656D+02,     1.72772750258445740D+00,        &
           2.02042913309661486D+04,    -9.69805983886375135D+04,        &
           1.92547001232531532D+05,    -2.03400177280415534D+05,        &
           1.22200464983017460D+05,    -4.11926549688975513D+04,        &
           7.10951430248936372D+03,    -4.93915304773088012D+02,        &
           6.07404200127348304D+00,    -2.42919187900551333D+05,        &
           1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),      &
           C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),      &
           C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/      &
           3.76327129765640400D+06,    -2.81356322658653411D+06,        &
           1.26836527332162478D+06,    -3.31645172484563578D+05,        &
           4.52187689813627263D+04,    -2.49983048181120962D+03,        &
           2.43805296995560639D+01,     3.28446985307203782D+06,        &
          -1.97068191184322269D+07,     5.09526024926646422D+07,        &
          -7.41051482115326577D+07,     6.63445122747290267D+07,        &
          -3.75671766607633513D+07,     1.32887671664218183D+07,        &
          -2.78561812808645469D+06,     3.08186404612662398D+05,        &
          -1.38860897537170405D+04,     1.10017140269246738D+02,        &
          -4.93292536645099620D+07,     3.25573074185765749D+08,        &
          -9.39462359681578403D+08,     1.55359689957058006D+09,        &
          -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),      &
           C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),      &
           C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/      &
          -4.95889784275030309D+08,     1.42062907797533095D+08,        &
          -2.44740627257387285D+07,     2.24376817792244943D+06,        &
          -8.40054336030240853D+04,     5.51335896122020586D+02,        &
           8.14789096118312115D+08,    -5.86648149205184723D+09,        &
           1.86882075092958249D+10,    -3.46320433881587779D+10,        &
           4.12801855797539740D+10,    -3.30265997498007231D+10,        &
           1.79542137311556001D+10,    -6.56329379261928433D+09,        &
           1.55927986487925751D+09,    -2.25105661889415278D+08,        &
           1.73951075539781645D+07,    -5.49842327572288687D+05,        &
           3.03809051092238427D+03,    -1.46792612476956167D+10,        &
           1.14498237732025810D+11,    -3.99096175224466498D+11,        &
           8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104), &
           C(105), C(106), C(107), C(108), C(109), C(110), C(111),      &
           C(112), C(113), C(114), C(115), C(116), C(117), C(118)/      &
           1.00815810686538209D+12,    -6.45364869245376503D+11,        &
           2.87900649906150589D+11,    -8.78670721780232657D+10,        &
           1.76347306068349694D+10,    -2.16716498322379509D+09,        &
           1.43157876718888981D+08,    -3.87183344257261262D+06,        &
           1.82577554742931747D+04,     2.86464035717679043D+11,        &
          -2.40629790002850396D+12,     9.10934118523989896D+12,        &
          -2.05168994109344374D+13,     3.05651255199353206D+13,        &
          -3.16670885847851584D+13,     2.33483640445818409D+13,        &
          -1.23204913055982872D+13,     4.61272578084913197D+12,        &
          -1.19655288019618160D+12,     2.05914503232410016D+11,        &
          -2.18229277575292237D+10,     1.24700929351271032D+09/
      DATA C(119), C(120)/  &
          -2.91883881222208134D+07,     1.18838426256783253D+05/

      IF (INIT.NE.0) GO TO 40
!-----------------------------------------------------------------------
!     INITIALIZE ALL VARIABLES
!-----------------------------------------------------------------------
      RFN = 1.0D0/FNU
      TR = ZRR*RFN
      TI = ZRI*RFN
      SR = CONER + (TR*TR-TI*TI)
      SI = CONEI + (TR*TI+TI*TR)
      CALL ZSQRT(SR, SI, SRR, SRI)
      STR = CONER + SRR
      STI = CONEI + SRI
      CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI)
      CALL ZLOG(ZNR, ZNI, STR, STI, IDUM)
      ZETA1R = FNU*STR
      ZETA1I = FNU*STI
      ZETA2R = FNU*SRR
      ZETA2I = FNU*SRI
      CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
      SRR = TR*RFN
      SRI = TI*RFN
      CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
      PHIR = CWRKR(16)*CON(IKFLG)
      PHII = CWRKI(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
      CWRKR(1) = CONER
      CWRKI(1) = CONEI
      CRFNR = CONER
      CRFNI = CONEI
      AC = 1.0D0
      L = 1
      DO 20 K=2,15
        SR = ZEROR
        SI = ZEROI
        DO 10 J=1,K
          L = L + 1
          STR = SR*T2R - SI*T2I + C(L)
          SI = SR*T2I + SI*T2R
          SR = STR
   10   CONTINUE
        STR = CRFNR*SRR - CRFNI*SRI
        CRFNI = CRFNR*SRI + CRFNI*SRR
        CRFNR = STR
        CWRKR(K) = CRFNR*SR - CRFNI*SI
        CWRKI(K) = CRFNR*SI + CRFNI*SR
        AC = AC*RFN
        TEST = DABS(CWRKR(K)) + DABS(CWRKI(K))
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE I FUNCTION
!-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      DO 50 I=1,INIT
        SR = SR + CWRKR(I)
        SI = SI + CWRKI(I)
   50 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(1)
      PHII = CWRKI(16)*CON(1)
      RETURN
   60 CONTINUE
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE K FUNCTION
!-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      TR = CONER
      DO 70 I=1,INIT
        SR = SR + TR*CWRKR(I)
        SI = SI + TR*CWRKI(I)
        TR = -TR
   70 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(2)
      PHII = CWRKI(16)*CON(2)
      RETURN
END

SUBROUTINE ZUNHJ(ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,  &
       ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
USE COMPLEX
!***BEGIN PROLOGUE  ZUNHJ
!***REFER TO  ZBESI,ZBESK
!
!     REFERENCES
!         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
!         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
!
!         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
!         PRESS, N.Y., 1974, PAGE 420
!
!     ABSTRACT
!         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
!         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
!         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
!
!         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
!
!         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
!         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
!
!               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
!
!         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
!         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
!
!         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
!         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
!         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
!
!***ROUTINES CALLED  ZABS,ZDIV,ZLOG,ZSQRT
!***END PROLOGUE  ZUNHJ
!     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
!    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
!    *ZETA2,ZTH
      DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,      &
       ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,   &
       CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,    &
       PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,  &
       RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,  &
       SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,    &
       TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,   &
       ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,     &
       ZETA2R, ZI, ZR, ZTHI, ZTHR
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,   &
       LRP1, L1, L2, M, IDUM
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),  &
       AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),       &
       DRR(14), DRI(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),  &
           AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/           &
           1.00000000000000000D+00,     1.04166666666666667D-01,    &
           8.35503472222222222D-02,     1.28226574556327160D-01,    &
           2.91849026464140464D-01,     8.81627267443757652D-01,    &
           3.32140828186276754D+00,     1.49957629868625547D+01,    &
           7.89230130115865181D+01,     4.74451538868264323D+02,    &
           3.20749009089066193D+03,     2.40865496408740049D+04,    &
           1.98923119169509794D+05,     1.79190200777534383D+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),  &
           BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/           &
           1.00000000000000000D+00,    -1.45833333333333333D-01,    &
          -9.87413194444444444D-02,    -1.43312053915895062D-01,    &
          -3.17227202678413548D-01,    -9.42429147957120249D-01,    &
          -3.51120304082635426D+00,    -1.57272636203680451D+01,    &
          -8.22814390971859444D+01,    -4.92355370523670524D+02,    &
          -3.31621856854797251D+03,    -2.48276742452085896D+04,    &
          -2.04526587315129788D+05,    -1.83844491706820990D+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
           C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),      &
           C(19), C(20), C(21), C(22), C(23), C(24)/                    &
           1.00000000000000000D+00,    -2.08333333333333333D-01,        &
           1.25000000000000000D-01,     3.34201388888888889D-01,        &
          -4.01041666666666667D-01,     7.03125000000000000D-02,        &
          -1.02581259645061728D+00,     1.84646267361111111D+00,        &
          -8.91210937500000000D-01,     7.32421875000000000D-02,        &
           4.66958442342624743D+00,    -1.12070026162229938D+01,        &
           8.78912353515625000D+00,    -2.36408691406250000D+00,        &
           1.12152099609375000D-01,    -2.82120725582002449D+01,        &
           8.46362176746007346D+01,    -9.18182415432400174D+01,        &
           4.25349987453884549D+01,    -7.36879435947963170D+00,        &
           2.27108001708984375D-01,     2.12570130039217123D+02,        &
          -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),      &
           C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),      &
           C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/      &
          -6.99579627376132541D+02,     2.18190511744211590D+02,        &
          -2.64914304869515555D+01,     5.72501420974731445D-01,        &
          -1.91945766231840700D+03,     8.06172218173730938D+03,        &
          -1.35865500064341374D+04,     1.16553933368645332D+04,        &
          -5.30564697861340311D+03,     1.20090291321635246D+03,        &
          -1.08090919788394656D+02,     1.72772750258445740D+00,        &
           2.02042913309661486D+04,    -9.69805983886375135D+04,        &
           1.92547001232531532D+05,    -2.03400177280415534D+05,        &
           1.22200464983017460D+05,    -4.11926549688975513D+04,        &
           7.10951430248936372D+03,    -4.93915304773088012D+02,        &
           6.07404200127348304D+00,    -2.42919187900551333D+05,        &
           1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),      &
           C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),      &
           C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/      &
           3.76327129765640400D+06,    -2.81356322658653411D+06,        &
           1.26836527332162478D+06,    -3.31645172484563578D+05,        &
           4.52187689813627263D+04,    -2.49983048181120962D+03,        &
           2.43805296995560639D+01,     3.28446985307203782D+06,        &
          -1.97068191184322269D+07,     5.09526024926646422D+07,        &
          -7.41051482115326577D+07,     6.63445122747290267D+07,        &
          -3.75671766607633513D+07,     1.32887671664218183D+07,        &
          -2.78561812808645469D+06,     3.08186404612662398D+05,        &
          -1.38860897537170405D+04,     1.10017140269246738D+02,        &
          -4.93292536645099620D+07,     3.25573074185765749D+08,        &
          -9.39462359681578403D+08,     1.55359689957058006D+09,        &
          -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),      &
           C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),      &
           C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/      &
          -4.95889784275030309D+08,     1.42062907797533095D+08,        &
          -2.44740627257387285D+07,     2.24376817792244943D+06,        &
          -8.40054336030240853D+04,     5.51335896122020586D+02,        &
           8.14789096118312115D+08,    -5.86648149205184723D+09,        &
           1.86882075092958249D+10,    -3.46320433881587779D+10,        &
           4.12801855797539740D+10,    -3.30265997498007231D+10,        &
           1.79542137311556001D+10,    -6.56329379261928433D+09,        &
           1.55927986487925751D+09,    -2.25105661889415278D+08,        &
           1.73951075539781645D+07,    -5.49842327572288687D+05,        &
           3.03809051092238427D+03,    -1.46792612476956167D+10,        &
           1.14498237732025810D+11,    -3.99096175224466498D+11,        &
           8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104), &
           C(105)/                                                      &
           1.00815810686538209D+12,    -6.45364869245376503D+11,        &
           2.87900649906150589D+11,    -8.78670721780232657D+10,        &
           1.76347306068349694D+10,    -2.16716498322379509D+09,        &
           1.43157876718888981D+08,    -3.87183344257261262D+06,        &
           1.82577554742931747D+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),        &
           ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),     &
           ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),  &
           ALFA(19), ALFA(20), ALFA(21), ALFA(22)/                      &
          -4.44444444444444444D-03,    -9.22077922077922078D-04,        &
          -8.84892884892884893D-05,     1.65927687832449737D-04,        &
           2.46691372741792910D-04,     2.65995589346254780D-04,        &
           2.61824297061500945D-04,     2.48730437344655609D-04,        &
           2.32721040083232098D-04,     2.16362485712365082D-04,        &
           2.00738858762752355D-04,     1.86267636637545172D-04,        &
           1.73060775917876493D-04,     1.61091705929015752D-04,        &
           1.50274774160908134D-04,     1.40503497391269794D-04,        &
           1.31668816545922806D-04,     1.23667445598253261D-04,        &
           1.16405271474737902D-04,     1.09798298372713369D-04,        &
           1.03772410422992823D-04,     9.82626078369363448D-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),  &
           ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),  &
           ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),  &
           ALFA(41), ALFA(42), ALFA(43), ALFA(44)/                      &
           9.32120517249503256D-05,     8.85710852478711718D-05,        &
           8.42963105715700223D-05,     8.03497548407791151D-05,        &
           7.66981345359207388D-05,     7.33122157481777809D-05,        &
           7.01662625163141333D-05,     6.72375633790160292D-05,        &
           6.93735541354588974D-04,     2.32241745182921654D-04,        &
          -1.41986273556691197D-05,    -1.16444931672048640D-04,        &
          -1.50803558053048762D-04,    -1.55121924918096223D-04,        &
          -1.46809756646465549D-04,    -1.33815503867491367D-04,        &
          -1.19744975684254051D-04,    -1.06184319207974020D-04,        &
          -9.37699549891194492D-05,    -8.26923045588193274D-05,        &
          -7.29374348155221211D-05,    -6.44042357721016283D-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),  &
           ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),  &
           ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),  &
           ALFA(63), ALFA(64), ALFA(65), ALFA(66)/                      &
          -5.69611566009369048D-05,    -5.04731044303561628D-05,        &
          -4.48134868008882786D-05,    -3.98688727717598864D-05,        &
          -3.55400532972042498D-05,    -3.17414256609022480D-05,        &
          -2.83996793904174811D-05,    -2.54522720634870566D-05,        &
          -2.28459297164724555D-05,    -2.05352753106480604D-05,        &
          -1.84816217627666085D-05,    -1.66519330021393806D-05,        &
          -1.50179412980119482D-05,    -1.35554031379040526D-05,        &
          -1.22434746473858131D-05,    -1.10641884811308169D-05,        &
          -3.54211971457743841D-04,    -1.56161263945159416D-04,        &
           3.04465503594936410D-05,     1.30198655773242693D-04,        &
           1.67471106699712269D-04,     1.70222587683592569D-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),  &
           ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),  &
           ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),  &
           ALFA(85), ALFA(86), ALFA(87), ALFA(88)/                      &
           1.56501427608594704D-04,     1.36339170977445120D-04,        &
           1.14886692029825128D-04,     9.45869093034688111D-05,        &
           7.64498419250898258D-05,     6.07570334965197354D-05,        &
           4.74394299290508799D-05,     3.62757512005344297D-05,        &
           2.69939714979224901D-05,     1.93210938247939253D-05,        &
           1.30056674793963203D-05,     7.82620866744496661D-06,        &
           3.59257485819351583D-06,     1.44040049814251817D-07,        &
          -2.65396769697939116D-06,    -4.91346867098485910D-06,        &
          -6.72739296091248287D-06,    -8.17269379678657923D-06,        &
          -9.31304715093561232D-06,    -1.02011418798016441D-05,        &
          -1.08805962510592880D-05,    -1.13875481509603555D-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),  &
           ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100), &
           ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),       &
           ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/       &
          -1.17519675674556414D-05,    -1.19987364870944141D-05,        &
           3.78194199201772914D-04,     2.02471952761816167D-04,        &
          -6.37938506318862408D-05,    -2.38598230603005903D-04,        &
          -3.10916256027361568D-04,    -3.13680115247576316D-04,        &
          -2.78950273791323387D-04,    -2.28564082619141374D-04,        &
          -1.75245280340846749D-04,    -1.25544063060690348D-04,        &
          -8.22982872820208365D-05,    -4.62860730588116458D-05,        &
          -1.72334302366962267D-05,     5.60690482304602267D-06,        &
           2.31395443148286800D-05,     3.62642745856793957D-05,        &
           4.58006124490188752D-05,     5.24595294959114050D-05,        &
           5.68396208545815266D-05,     5.94349820393104052D-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),       &
           ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),       &
           ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),       &
           ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/       &
           6.06478527578421742D-05,     6.08023907788436497D-05,        &
           6.01577894539460388D-05,     5.89199657344698500D-05,        &
           5.72515823777593053D-05,     5.52804375585852577D-05,        &
           5.31063773802880170D-05,     5.08069302012325706D-05,        &
           4.84418647620094842D-05,     4.60568581607475370D-05,        &
          -6.91141397288294174D-04,    -4.29976633058871912D-04,        &
           1.83067735980039018D-04,     6.60088147542014144D-04,        &
           8.75964969951185931D-04,     8.77335235958235514D-04,        &
           7.49369585378990637D-04,     5.63832329756980918D-04,        &
           3.68059319971443156D-04,     1.88464535514455599D-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),       &
           ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),       &
           ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),       &
           ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/       &
           3.70663057664904149D-05,    -8.28520220232137023D-05,        &
          -1.72751952869172998D-04,    -2.36314873605872983D-04,        &
          -2.77966150694906658D-04,    -3.02079514155456919D-04,        &
          -3.12594712643820127D-04,    -3.12872558758067163D-04,        &
          -3.05678038466324377D-04,    -2.93226470614557331D-04,        &
          -2.77255655582934777D-04,    -2.59103928467031709D-04,        &
          -2.39784014396480342D-04,    -2.20048260045422848D-04,        &
          -2.00443911094971498D-04,    -1.81358692210970687D-04,        &
          -1.63057674478657464D-04,    -1.45712672175205844D-04,        &
          -1.29425421983924587D-04,    -1.14245691942445952D-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),       &
           ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),       &
           ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),       &
           ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/       &
           1.92821964248775885D-03,     1.35592576302022234D-03,        &
          -7.17858090421302995D-04,    -2.58084802575270346D-03,        &
          -3.49271130826168475D-03,    -3.46986299340960628D-03,        &
          -2.82285233351310182D-03,    -1.88103076404891354D-03,        &
          -8.89531718383947600D-04,     3.87912102631035228D-06,        &
           7.28688540119691412D-04,     1.26566373053457758D-03,        &
           1.62518158372674427D-03,     1.83203153216373172D-03,        &
           1.91588388990527909D-03,     1.90588846755546138D-03,        &
           1.82798982421825727D-03,     1.70389506421121530D-03,        &
           1.55097127171097686D-03,     1.38261421852276159D-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),       &
           ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/       &
           1.20881424230064774D-03,     1.03676532638344962D-03,        &
           8.71437918068619115D-04,     7.16080155297701002D-04,        &
           5.72637002558129372D-04,     4.42089819465802277D-04,        &
           3.24724948503090564D-04,     2.20342042730246599D-04,        &
           1.28412898401353882D-04,     4.82005924552095464D-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),        &
           BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),     &
           BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),  &
           BETA(19), BETA(20), BETA(21), BETA(22)/                      &
           1.79988721413553309D-02,     5.59964911064388073D-03,        &
           2.88501402231132779D-03,     1.80096606761053941D-03,        &
           1.24753110589199202D-03,     9.22878876572938311D-04,        &
           7.14430421727287357D-04,     5.71787281789704872D-04,        &
           4.69431007606481533D-04,     3.93232835462916638D-04,        &
           3.34818889318297664D-04,     2.88952148495751517D-04,        &
           2.52211615549573284D-04,     2.22280580798883327D-04,        &
           1.97541838033062524D-04,     1.76836855019718004D-04,        &
           1.59316899661821081D-04,     1.44347930197333986D-04,        &
           1.31448068119965379D-04,     1.20245444949302884D-04,        &
           1.10449144504599392D-04,     1.01828770740567258D-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),  &
           BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),  &
           BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),  &
           BETA(41), BETA(42), BETA(43), BETA(44)/                      &
           9.41998224204237509D-05,     8.74130545753834437D-05,        &
           8.13466262162801467D-05,     7.59002269646219339D-05,        &
           7.09906300634153481D-05,     6.65482874842468183D-05,        &
           6.25146958969275078D-05,     5.88403394426251749D-05,        &
          -1.49282953213429172D-03,    -8.78204709546389328D-04,        &
          -5.02916549572034614D-04,    -2.94822138512746025D-04,        &
          -1.75463996970782828D-04,    -1.04008550460816434D-04,        &
          -5.96141953046457895D-05,    -3.12038929076098340D-05,        &
          -1.26089735980230047D-05,    -2.42892608575730389D-07,        &
           8.05996165414273571D-06,     1.36507009262147391D-05,        &
           1.73964125472926261D-05,     1.98672978842133780D-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),  &
           BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),  &
           BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),  &
           BETA(63), BETA(64), BETA(65), BETA(66)/                      &
           2.14463263790822639D-05,     2.23954659232456514D-05,        &
           2.28967783814712629D-05,     2.30785389811177817D-05,        &
           2.30321976080909144D-05,     2.28236073720348722D-05,        &
           2.25005881105292418D-05,     2.20981015361991429D-05,        &
           2.16418427448103905D-05,     2.11507649256220843D-05,        &
           2.06388749782170737D-05,     2.01165241997081666D-05,        &
           1.95913450141179244D-05,     1.90689367910436740D-05,        &
           1.85533719641636667D-05,     1.80475722259674218D-05,        &
           5.52213076721292790D-04,     4.47932581552384646D-04,        &
           2.79520653992020589D-04,     1.52468156198446602D-04,        &
           6.93271105657043598D-05,     1.76258683069991397D-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),  &
           BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),  &
           BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),  &
           BETA(85), BETA(86), BETA(87), BETA(88)/                      &
          -1.35744996343269136D-05,    -3.17972413350427135D-05,        &
          -4.18861861696693365D-05,    -4.69004889379141029D-05,        &
          -4.87665447413787352D-05,    -4.87010031186735069D-05,        &
          -4.74755620890086638D-05,    -4.55813058138628452D-05,        &
          -4.33309644511266036D-05,    -4.09230193157750364D-05,        &
          -3.84822638603221274D-05,    -3.60857167535410501D-05,        &
          -3.37793306123367417D-05,    -3.15888560772109621D-05,        &
          -2.95269561750807315D-05,    -2.75978914828335759D-05,        &
          -2.58006174666883713D-05,    -2.41308356761280200D-05,        &
          -2.25823509518346033D-05,    -2.11479656768912971D-05,        &
          -1.98200638885294927D-05,    -1.85909870801065077D-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),  &
           BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100), &
           BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),       &
           BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/       &
          -1.74532699844210224D-05,    -1.63997823854497997D-05,        &
          -4.74617796559959808D-04,    -4.77864567147321487D-04,        &
          -3.20390228067037603D-04,    -1.61105016119962282D-04,        &
          -4.25778101285435204D-05,     3.44571294294967503D-05,        &
           7.97092684075674924D-05,     1.03138236708272200D-04,        &
           1.12466775262204158D-04,     1.13103642108481389D-04,        &
           1.08651634848774268D-04,     1.01437951597661973D-04,        &
           9.29298396593363896D-05,     8.40293133016089978D-05,        &
           7.52727991349134062D-05,     6.69632521975730872D-05,        &
           5.92564547323194704D-05,     5.22169308826975567D-05,        &
           4.58539485165360646D-05,     4.01445513891486808D-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),       &
           BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),       &
           BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),       &
           BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/       &
           3.50481730031328081D-05,     3.05157995034346659D-05,        &
           2.64956119950516039D-05,     2.29363633690998152D-05,        &
           1.97893056664021636D-05,     1.70091984636412623D-05,        &
           1.45547428261524004D-05,     1.23886640995878413D-05,        &
           1.04775876076583236D-05,     8.79179954978479373D-06,        &
           7.36465810572578444D-04,     8.72790805146193976D-04,        &
           6.22614862573135066D-04,     2.85998154194304147D-04,        &
           3.84737672879366102D-06,    -1.87906003636971558D-04,        &
          -2.97603646594554535D-04,    -3.45998126832656348D-04,        &
          -3.53382470916037712D-04,    -3.35715635775048757D-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),       &
           BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),       &
           BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),       &
           BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/       &
          -3.04321124789039809D-04,    -2.66722723047612821D-04,        &
          -2.27654214122819527D-04,    -1.89922611854562356D-04,        &
          -1.55058918599093870D-04,    -1.23778240761873630D-04,        &
          -9.62926147717644187D-05,    -7.25178327714425337D-05,        &
          -5.22070028895633801D-05,    -3.50347750511900522D-05,        &
          -2.06489761035551757D-05,    -8.70106096849767054D-06,        &
           1.13698686675100290D-06,     9.16426474122778849D-06,        &
           1.56477785428872620D-05,     2.08223629482466847D-05,        &
           2.48923381004595156D-05,     2.80340509574146325D-05,        &
           3.03987774629861915D-05,     3.21156731406700616D-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),       &
           BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),       &
           BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),       &
           BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/       &
          -1.80182191963885708D-03,    -2.43402962938042533D-03,        &
          -1.83422663549856802D-03,    -7.62204596354009765D-04,        &
           2.39079475256927218D-04,     9.49266117176881141D-04,        &
           1.34467449701540359D-03,     1.48457495259449178D-03,        &
           1.44732339830617591D-03,     1.30268261285657186D-03,        &
           1.10351597375642682D-03,     8.86047440419791759D-04,        &
           6.73073208165665473D-04,     4.77603872856582378D-04,        &
           3.05991926358789362D-04,     1.60315694594721630D-04,        &
           4.00749555270613286D-05,    -5.66607461635251611D-05,        &
          -1.32506186772982638D-04,    -1.90296187989614057D-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),       &
           BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),       &
           BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),       &
           BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/       &
          -2.32811450376937408D-04,    -2.62628811464668841D-04,        &
          -2.82050469867598672D-04,    -2.93081563192861167D-04,        &
          -2.97435962176316616D-04,    -2.96557334239348078D-04,        &
          -2.91647363312090861D-04,    -2.83696203837734166D-04,        &
          -2.73512317095673346D-04,    -2.61750155806768580D-04,        &
           6.38585891212050914D-03,     9.62374215806377941D-03,        &
           7.61878061207001043D-03,     2.83219055545628054D-03,        &
          -2.09841352012720090D-03,    -5.73826764216626498D-03,        &
          -7.70804244495414620D-03,    -8.21011692264844401D-03,        &
          -7.65824520346905413D-03,    -6.47209729391045177D-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),       &
           BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),       &
           BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),       &
           BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/       &
          -4.99132412004966473D-03,    -3.45612289713133280D-03,        &
          -2.01785580014170775D-03,    -7.59430686781961401D-04,        &
           2.84173631523859138D-04,     1.10891667586337403D-03,        &
           1.72901493872728771D-03,     2.16812590802684701D-03,        &
           2.45357710494539735D-03,     2.61281821058334862D-03,        &
           2.67141039656276912D-03,     2.65203073395980430D-03,        &
           2.57411652877287315D-03,     2.45389126236094427D-03,        &
           2.30460058071795494D-03,     2.13684837686712662D-03,        &
           1.95896528478870911D-03,     1.77737008679454412D-03,        &
           1.59690280765839059D-03,     1.42111975664438546D-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),        &
           GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),     &
           GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),  &
           GAMA(19), GAMA(20), GAMA(21), GAMA(22)/                      &
           6.29960524947436582D-01,     2.51984209978974633D-01,        &
           1.54790300415655846D-01,     1.10713062416159013D-01,        &
           8.57309395527394825D-02,     6.97161316958684292D-02,        &
           5.86085671893713576D-02,     5.04698873536310685D-02,        &
           4.42600580689154809D-02,     3.93720661543509966D-02,        &
           3.54283195924455368D-02,     3.21818857502098231D-02,        &
           2.94646240791157679D-02,     2.71581677112934479D-02,        &
           2.51768272973861779D-02,     2.34570755306078891D-02,        &
           2.19508390134907203D-02,     2.06210828235646240D-02,        &
           1.94388240897880846D-02,     1.83810633800683158D-02,        &
           1.74293213231963172D-02,     1.65685837786612353D-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),  &
           GAMA(29), GAMA(30)/                                          &
           1.57865285987918445D-02,     1.50729501494095594D-02,        &
           1.44193250839954639D-02,     1.38184805735341786D-02,        &
           1.32643378994276568D-02,     1.27517121970498651D-02,        &
           1.22761545318762767D-02,     1.18338262398482403D-02/
      DATA EX1, EX2, HPI, GPI, THPI /                                   &
           3.33333333333333333D-01,     6.66666666666666667D-01,        &
           1.57079632679489662D+00,     3.14159265358979324D+00,        &
           4.71238898038468986D+00/
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /

      RFNU = 1.0D0/FNU
      ZBR = ZR*RFNU
      ZBI = ZI*RFNU
      RFNU2 = RFNU*RFNU
!-----------------------------------------------------------------------
!     COMPUTE IN THE FOURTH QUADRANT
!-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = 1.0D0/FN13
      W2R = CONER - ZBR*ZBR + ZBI*ZBI
      W2I = CONEI - ZBR*ZBI - ZBR*ZBI
      AW2 = ZABS(W2R,W2I)
      IF (AW2.GT.0.25D0) GO TO 130
!-----------------------------------------------------------------------
!     POWER SERIES FOR CABS(W2).LE.0.25D0
!-----------------------------------------------------------------------
      K = 1
      PR(1) = CONER
      PI(1) = CONEI
      SUMAR = GAMA(1)
      SUMAI = ZEROI
      AP(1) = 1.0D0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
        PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
        SUMAR = SUMAR + PR(K)*GAMA(K)
        SUMAI = SUMAI + PI(K)*GAMA(K)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETAR = W2R*SUMAR - W2I*SUMAI
      ZETAI = W2R*SUMAI + W2I*SUMAR
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI)
      CALL ZSQRT(W2R, W2I, STR, STI)
      ZETA2R = STR*FNU
      ZETA2I = STI*FNU
      STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
      STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
      ZETA1R = STR*ZETA2R - STI*ZETA2I
      ZETA1I = STR*ZETA2I + STI*ZETA2R
      ZAR = ZAR + ZAR
      ZAI = ZAI + ZAI
      CALL ZSQRT(ZAR, ZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
!-----------------------------------------------------------------------
!     SUM SERIES FOR ASUM AND BSUM
!-----------------------------------------------------------------------
      SUMBR = ZEROR
      SUMBI = ZEROI
      DO 30 K=1,KMAX
        SUMBR = SUMBR + PR(K)*BETA(K)
        SUMBI = SUMBI + PI(K)*BETA(K)
   30 CONTINUE
      ASUMR = ZEROR
      ASUMI = ZEROI
      BSUMR = SUMBR
      BSUMI = SUMBI
      L1 = 0
      L2 = 30
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      ATOL = TOL
      PP = 1.0D0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMAR = ZEROR
        SUMAI = ZEROI
        DO 40 K=1,KMAX
          M = L1 + K
          SUMAR = SUMAR + PR(K)*ALFA(M)
          SUMAI = SUMAI + PI(K)*ALFA(M)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUMR = ASUMR + SUMAR*PP
        ASUMI = ASUMI + SUMAI*PP
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMBR = ZEROR
        SUMBI = ZEROI
        DO 70 K=1,KMAX
          M = L2 + K
          SUMBR = SUMBR + PR(K)*BETA(M)
          SUMBI = SUMBI + PI(K)*BETA(M)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUMR = BSUMR + SUMBR*PP
        BSUMI = BSUMI + SUMBI*PP
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUMR = ASUMR + CONER
      PP = RFNU*RFN13
      BSUMR = BSUMR*PP
      BSUMI = BSUMI*PP
  120 CONTINUE
      RETURN
!-----------------------------------------------------------------------
!     CABS(W2).GT.0.25D0
!-----------------------------------------------------------------------
  130 CONTINUE
      CALL ZSQRT(W2R, W2I, WR, WI)
      IF (WR.LT.0.0D0) WR = 0.0D0
      IF (WI.LT.0.0D0) WI = 0.0D0
      STR = CONER + WR
      STI = WI
      CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
      CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
      IF (ZCI.LT.0.0D0) ZCI = 0.0D0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0D0) ZCR = 0.0D0
      ZTHR = (ZCR-WR)*1.5D0
      ZTHI = (ZCI-WI)*1.5D0
      ZETA1R = ZCR*FNU
      ZETA1I = ZCI*FNU
      ZETA2R = WR*FNU
      ZETA2I = WI*FNU
      AZTH = ZABS(ZTHR,ZTHI)
      ANG = THPI
      IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0D0) GO TO 140
      ANG = DATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*DCOS(ANG)
      ZETAI = PP*DSIN(ANG)
      IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
      CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
      TZAR = ZAR + ZAR
      TZAI = ZAI + ZAI
      CALL ZSQRT(TZAR, TZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      RAW = 1.0D0/DSQRT(AW2)
      STR = WR*RAW
      STI = -WI*RAW
      TFNR = STR*RFNU*RAW
      TFNI = STI*RFNU*RAW
      RAZTH = 1.0D0/AZTH
      STR = ZTHR*RAZTH
      STI = -ZTHI*RAZTH
      RZTHR = STR*RAZTH*RFNU
      RZTHI = STI*RAZTH*RFNU
      ZCR = RZTHR*AR(2)
      ZCI = RZTHI*AR(2)
      RAW2 = 1.0D0/AW2
      STR = W2R*RAW2
      STI = -W2I*RAW2
      T2R = STR*RAW2
      T2I = STI*RAW2
      STR = T2R*C(2) + C(3)
      STI = T2I*C(2)
      UPR(2) = STR*TFNR - STI*TFNI
      UPI(2) = STR*TFNI + STI*TFNR
      BSUMR = UPR(2) + ZCR
      BSUMI = UPI(2) + ZCI
      ASUMR = ZEROR
      ASUMI = ZEROI
      IF (RFNU.LT.TOL) GO TO 220
      PRZTHR = RZTHR
      PRZTHI = RZTHI
      PTFNR = TFNR
      PTFNI = TFNI
      UPR(1) = CONER
      UPI(1) = CONEI
      PP = 1.0D0
      BTOL = TOL*(DABS(BSUMR)+DABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
!-----------------------------------------------------------------------
!     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
!     NEXT SUMA AND SUMB
!-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZAR = C(L)
          ZAI = ZEROI
          DO 150 J=2,KP1
            L = L + 1
            STR = ZAR*T2R - T2I*ZAI + C(L)
            ZAI = ZAR*T2I + ZAI*T2R
            ZAR = STR
  150     CONTINUE
          STR = PTFNR*TFNR - PTFNI*TFNI
          PTFNI = PTFNR*TFNI + PTFNI*TFNR
          PTFNR = STR
          UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
          UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
          CRR(KS) = PRZTHR*BR(KS+1)
          CRI(KS) = PRZTHI*BR(KS+1)
          STR = PRZTHR*RZTHR - PRZTHI*RZTHI
          PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
          PRZTHR = STR
          DRR(KS) = PRZTHR*AR(KS+2)
          DRI(KS) = PRZTHI*AR(KS+2)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMAR = UPR(LRP1)
        SUMAI = UPI(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
          SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170   CONTINUE
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = DABS(SUMAR) + DABS(SUMAI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
        SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
          SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190   CONTINUE
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = DABS(SUMBR) + DABS(SUMBI)
        IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUMR = ASUMR + CONER
      STR = -BSUMR*RFN13
      STI = -BSUMI*RFN13
      CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
      GO TO 120
END

SUBROUTINE ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)
USE COMPLEX
!***BEGIN PROLOGUE  ZRATI
!***REFER TO  ZBESI,ZBESK,ZBESH
!
!     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
!     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
!     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
!     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
!     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
!     BY D. J. SOOKNE.
!
!***ROUTINES CALLED  ZABS,ZDIV
!***END PROLOGUE  ZRATI
!     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR, &
      CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU, &
      FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,  &
      RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(1), CYI(1)
      DATA CZEROR,CZEROI,CONER,CONEI,RT2  &
      /0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0/
      AZ = ZABS(ZR,ZI)
      INU = INT(SNGL(FNU))
      IDNU = INU + N - 1
      MAGZ = INT(SNGL(AZ))
      AMAGZ = DBLE(FLOAT(MAGZ+1))
      FDNU = DBLE(FLOAT(IDNU))
      FNUP = DMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0D0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
!-----------------------------------------------------------------------
!     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
!     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
!     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
!     PREMATURELY.
!-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = DSQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5D0
      FLAM = AK + DSQRT(AK*AK-1.0D0)
      RHO = DMIN1(AP2/AP1,FLAM)
      TEST = TEST1*DSQRT(RHO/(RHO*RHO-1.0D0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = DBLE(FLOAT(KK))
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + DBLE(FLOAT(N-1))
      P1R = 1.0D0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = DBLE(FLOAT(K))
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
END

SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZAIRY
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
!***DESCRIPTION
!
!                      ***A DOUBLE PRECISION ROUTINE***
!         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
!         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
!         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
!         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
!         -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
!         PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
!
!         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTI! IN
!         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
!         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
!         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
!         MATHEMATICAL FUNCTIONS (REF. 1).
!
!         INPUT      ZR,ZI ARE DOUBLE PRECISION
!           ZR,ZI  - Z=CMPLX(ZR,ZI)
!           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE= 1  RETURNS
!                             AI=AI(Z)                ON ID=0 OR
!                             AI=DAI(Z)/DZ            ON ID=1
!                        = 2  RETURNS
!                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
!                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
!                             ZTA=(2/3)*Z*CSQRT(Z)
!
!         OUTPUT     AIR,AII ARE DOUBLE PRECISION
!           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
!                    KODE
!           NZ     - UNDERFLOW INDICATOR
!                    NZ= 0   , NORMAL RETURN
!                    NZ= 1   , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN
!                              -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
!                            TOO LARGE ON KODE=1
!                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
!                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
!                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
!                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
!                            REDUCTION
!                    IERR=5, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!***LONG DESCRIPTION
!
!         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
!         FUNCTIONS BY
!
!            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
!                           C=1.0/(PI*SQRT(3.0))
!                            ZTA=(2/3)*Z**(3/2)
!
!         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
!         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
!         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
!         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
!         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
!         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
!         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
!         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
!         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!         LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
!         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
!         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
!         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
!         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
!         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
!         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
!         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
!         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
!         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
!         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
!         MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
!         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!         OR -PI/2+P.
!
!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
!                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
!                 COMMERCE, 1955.
!
!               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
!
!               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,I1MACH,D1MACH
!***END PROLOGUE  ZAIRY
!     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,  &
      CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,   &
      DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,       &
      S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI, &
      ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, ALAZ, BB
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ
      DIMENSION CYR(1), CYI(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667D-01, &
      3.55028053887817240D-01,2.58819403792806799D-01, &
      1.83776298473930683D-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
!***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = DMAX1(D1MACH(4),1.0D-18)
      FID = DBLE(FLOAT(ID))
      IF (AZ.GT.1.0D0) GO TO 70
!-----------------------------------------------------------------------
!     POWER SERIES FOR CABS(Z).LE.1.
!-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 170
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = DMIN1(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0D0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
!-----------------------------------------------------------------------
!     CASE FOR CABS(Z).GT.1.0
!-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      ALAZ = DLOG(AZ)
!-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA=DMIN1(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
!-----------------------------------------------------------------------
!     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -DABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0) GO TO 90
      IF (ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      IF (KODE.EQ.2) GO TO 100
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 100
      AA = -AA + 0.25D0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 270
  100 CONTINUE
!-----------------------------------------------------------------------
!     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL, &
      ELIM, ALIM)
      IF (NN.LT.0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      IF (KODE.EQ.2) GO TO 120
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 120
      AA = -AA - 0.25D0*ALAZ
      IFLAG = 2
      SFAC = 1.0D0/TOL
      IF (AA.LT.(-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM, ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      IF (IFLAG.NE.0) GO TO 150
      IF (ID.EQ.1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      IF (ID.EQ.1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0D+3*D1MACH(1)
      S1R = ZEROR
      S1I = ZEROI
      IF (ID.EQ.1) GO TO 190
      IF (AZ.LE.AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0D0
      AA = DSQRT(AA)
      IF (AZ.LE.AA) GO TO 200
      S1R = 0.5D0*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      IF(NN.EQ.(-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
END

SUBROUTINE ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL, ELIM, ALIM)
USE UTILIT
USE COMPLEX
!***BEGIN PROLOGUE  ZACAI
!***REFER TO  ZAIRY
!
!     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
!     IS CALLED FROM ZAIRY.
!
!***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,ZABS
!***END PROLOGUE  ZACAI
!     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,     &
      CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,  &
      RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(1), YI(1), CYR(2), CYI(2)
      DATA PI / 3.14159265358979324D0 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + DBLE(FLOAT(N-1))
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
!-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
      CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 80
      GO TO 40
   30 CONTINUE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      IF(NW.LT.0) GO TO 80
   40 CONTINUE
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      FMR = DBLE(FLOAT(MR))
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0D0
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*DSIN(YY)
      CSGNI = CSGNI*DCOS(YY)
   50 CONTINUE
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-DBLE(FLOAT(INU)))*SGN
      CSPNR = DCOS(ARG)
      CSPNI = DSIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      IF (KODE.EQ.1) GO TO 70
      IUF = 0
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
END

SUBROUTINE ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM, IUF)
USE COMPLEX
!***BEGIN PROLOGUE  ZS1S2
!***REFER TO  ZBESK,ZAIRY
!
!     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
!     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
!     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
!     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
!     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
!     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
!     PRECISION ABOVE THE UNDERFLOW LIMIT.
!
!***ROUTINES CALLED  ZABS,ZEXP,ZLOG
!***END PROLOGUE  ZS1S2
!     COMPLEX CZERO,C1,S1,S1D,S2,ZR
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI, &
      S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR
      INTEGER IUF, IDUM, NZ
      DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
      NZ = 0
      AS1 = ZABS(S1R,S1I)
      AS2 = ZABS(S2R,S2I)
      IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
      IF (AS1.EQ.0.0D0) GO TO 10
      ALN = -ZRR - ZRR + DLOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      IF (ALN.LT.(-ALIM)) GO TO 10
      CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZEXP(C1R, C1I, S1R, S1I)
      AS1 = ZABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = DMAX1(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
END
end module besselj
!end of file tzbesj.f90
