! ---------------------------------------------------------------------
!      Utility subroutines used by any program from Numath library
!      with (not intrinsic) complex numbers z = (zr,zi).
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               F90 Release 1.0 By J-P Moreau, Paris
!                                         (www.jpmoreau.fr)
! ---------------------------------------------------------------------
MODULE COMPLEX

CONTAINS

SUBROUTINE ZSQRT(AR, AI, BR, BI)
!***BEGIN PROLOGUE  ZSQRT
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
!
!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZSQRT
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DATA DRT , DPI / 7.071067811865475244008443621D-1, &
                       3.141592653589793238462643383D+0/
      ZM = ZABS(AR,AI)
      ZM = DSQRT(ZM)
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.GT.0.0D+0) GO TO 60
      IF (AI.LT.0.0D+0) GO TO 70
      BR = 0.0D+0
      BI = 0.0D+0
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = 0.0D+0
      BI = DSQRT(DABS(AR))
      RETURN
   30 BR = DSQRT(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
      BR = ZM*DCOS(DTHETA)
      BI = ZM*DSIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
END SUBROUTINE ZSQRT

SUBROUTINE ZEXP(AR, AI, BR, BI)
!***BEGIN PROLOGUE  ZEXP
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
      ZM = DEXP(AR)
      CA = ZM*DCOS(AI)
      CB = ZM*DSIN(AI)
      BR = CA
      BI = CB
RETURN
END SUBROUTINE ZEXP

SUBROUTINE ZMLT(AR, AI, BR, BI, CR, CI)
!***BEGIN PROLOGUE  ZMLT
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZMLT
DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
RETURN
END SUBROUTINE ZMLT

SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)
!***BEGIN PROLOGUE  ZDIV
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.

!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZDIV
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      BM = 1.0D0/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
RETURN
END SUBROUTINE ZDIV

SUBROUTINE ZLOG(AR, AI, BR, BI, IERR)
!***BEGIN PROLOGUE  ZLOG
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY

!     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
!     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZLOG
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      DATA DPI , DHPI  / 3.141592653589793238462643383D+0, &
                         1.570796326794896619231321696D+0/

      IERR=0
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0D+0) GO TO 60
      BI = DHPI
      BR = DLOG(DABS(AI))
      IF (AI.LT.0.0D+0) BI = -BI
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = DLOG(DABS(AR))
      BI = DPI
      RETURN
   30 BR = DLOG(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = DLOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
RETURN
END SUBROUTINE ZLOG

DOUBLE PRECISION FUNCTION ZABS(ZR, ZI)
!***BEGIN PROLOGUE  ZABS
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY

!     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
!     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
      U = DABS(ZR)
      V = DABS(ZI)
      S = U + V
!-----------------------------------------------------------------------
!     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CD! MACHINES INTO A
!     TRUE FLOATING ZERO
!-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZABS = V*DSQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS = U*DSQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS = 0.0D+0
      RETURN
END FUNCTION ZABS

END MODULE COMPLEX
! end of file Complex.f90
