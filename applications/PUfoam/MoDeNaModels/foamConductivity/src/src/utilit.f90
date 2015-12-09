! ---------------------------------------------------------------------
!      Utility subroutines used by any program from Numath library
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               F90 Release 1.0 By J-P Moreau, Paris
!                                         (www.jpmoreau.fr)
! ---------------------------------------------------------------------
MODULE UTILIT

CONTAINS

!REAL*8 FUNCTION D1MACH(I)
!!***BEGIN PROLOGUE  D1MACH
!!***DATE WRITTEN   750101   (YYMMDD)
!!***REVISION DATE  860501   (YYMMDD)
!!***CATEGORY NO.  R1
!!***KEYWORDS  MACHINE CONSTANTS
!!***AUTHOR  FOX, P. A., (BELL LABS)
!!           HALL, A. D., (BELL LABS)
!!           SCHRYER, N. L., (BELL LABS)
!!***PURPOSE  RETURN DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS.
!!***DESCRIPTION
!
!!     D1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!!     SUBPROGRAM WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!!     AS FOLLOWS, FOR EXAMPLE
!
!!          D = D1MACH(I)
!
!!     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF D ABOVE IS
!!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
!
!!  DOUBLE-PRECISION MACHINE CONSTANTS
!!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!!  D1MACH( 5) = LOG10(B)
!!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!!***ROUTINES CALLED  XERROR
!!***END PROLOGUE  D1MACH
!
!      INTEGER SMALL(4)
!      INTEGER LARGE(4)
!      INTEGER RIGHT(4)
!      INTEGER DIVER(4)
!      INTEGER LOG10(4)
!
!      DOUBLE PRECISION DMACH(5)
!
!      EQUIVALENCE (DMACH(1),SMALL(1))
!      EQUIVALENCE (DMACH(2),LARGE(1))
!      EQUIVALENCE (DMACH(3),RIGHT(1))
!      EQUIVALENCE (DMACH(4),DIVER(1))
!      EQUIVALENCE (DMACH(5),LOG10(1))
!
!!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!      DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
!      DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!      DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
!      DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
!      DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /

!     MACHINE CONSTANTS FOR THE APOLLO DNXXXX SERIES,

!     DATA DMACH(1) / 2.22559D-308/
!     DATA DMACH(2) / 1.79728D308/
!     DATA DMACH(3) / 1.11048D-16  /
!     DATA DMACH(4) / 2.22096D-16  /
!     DATA DMACH(5) / .301029995663981198D0  /

!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.

!     DATA SMALL(1) / ZC00800000 /
!     DATA SMALL(2) / Z000000000 /

!     DATA LARGE(1) / ZDFFFFFFFF /
!     DATA LARGE(2) / ZFFFFFFFFF /

!     DATA RIGHT(1) / ZCC5800000 /
!     DATA RIGHT(2) / Z000000000 /

!     DATA DIVER(1) / ZCC6800000 /
!     DATA DIVER(2) / Z000000000 /

!     DATA LOG10(1) / ZD00E730E7 /
!     DATA LOG10(2) / ZC77800DC0 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.

!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O0000000000000000 /

!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O0007777777777777 /

!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /

!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /

!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.

!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O7770000000000000 /

!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O7777777777777777 /

!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /

!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /

!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /

!     MACHINE CONSTANTS FOR THE CD! 6000/7000 SERIES.
!     FOR FTN4

!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /

!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /

!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /

!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /

!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /

!     MACHINE CONSTANTS FOR THE CD! 6000/7000 SERIES.
!     FOR FTN5

!     DATA SMALL(1) / O"00564000000000000000" /
!     DATA SMALL(2) / O"00000000000000000000" /

!     DATA LARGE(1) / O"37757777777777777777" /
!     DATA LARGE(2) / O"37157777777777777777" /

!     DATA RIGHT(1) / O"15624000000000000000" /
!     DATA RIGHT(2) / O"00000000000000000000" /

!     DATA DIVER(1) / O"15634000000000000000" /
!     DATA DIVER(2) / O"00000000000000000000" /

!     DATA LOG10(1) / O"17164642023241175717" /
!     DATA LOG10(2) / O"16367571421742254654" /

!     MACHINE CONSTANTS FOR THE CRAY 1

!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /

!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777774B /

!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /

!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /

!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /

!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200

!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATI! DMACH(5)

!     DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
!     DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
!     DATA LOG10/40423K,42023K,50237K,74776K/

!     MACHINE CONSTANTS FOR THE HARRIS 220

!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
!     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /

!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.

!     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /

!     MACHINE CONSTANTS FOR THE HP 2100
!     THREE WORD DOUBLE PRECISION OPTION WITH FTN4

!     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /

!     MACHINE CONSTANTS FOR THE HP 2100
!     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4

!     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!     DATA SMALL(3), SMALL(4) /       0,       1 /
!     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!     DATA DIVER(3), DIVER(4) /       0,    227B /
!     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!     DATA LOG10(3), LOG10(4) /  76747B, 176377B /

!     MACHINE CONSTANTS FOR THE HP 9000

!     D1MACH(1) = 2.8480954D-306
!     D1MACH(2) = 1.40444776D+306
!     D1MACH(3) = 2.22044605D-16
!     D1MACH(4) = 4.44089210D-16
!     D1MACH(5) = 3.01029996D-1

!     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
!     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
!     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
!     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
!     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /

!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).

!     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /

!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).

!     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).

!     DATA SMALL(1), SMALL(2) /    8388608,           0 /
!     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
!     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
!     DATA DIVER(1), DIVER(2) /  620756992,           0 /
!     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /

!     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
!     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
!     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
!     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
!     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).

!     DATA SMALL(1), SMALL(2) /    128,      0 /
!     DATA SMALL(3), SMALL(4) /      0,      0 /

!     DATA LARGE(1), LARGE(2) /  32767,     -1 /
!     DATA LARGE(3), LARGE(4) /     -1,     -1 /

!     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
!     DATA RIGHT(3), RIGHT(4) /      0,      0 /

!     DATA DIVER(1), DIVER(2) /   9472,      0 /
!     DATA DIVER(3), DIVER(4) /      0,      0 /

!     DATA LOG10(1), LOG10(2) /  16282,   8346 /
!     DATA LOG10(3), LOG10(4) / -31493, -12296 /

!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA SMALL(3), SMALL(4) / O000000, O000000 /

!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA LARGE(3), LARGE(4) / O177777, O177777 /

!     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
!     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /

!     DATA DIVER(1), DIVER(2) / O022400, O000000 /
!     DATA DIVER(3), DIVER(4) / O000000, O000000 /

!     DATA LOG10(1), LOG10(2) / O037632, O020232 /
!     DATA LOG10(3), LOG10(4) / O102373, O147770 /

!     MACHINE CONSTANTS FOR THE UNIVA! 1100 SERIES. FTN COMPILER

!     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /

!     MACHINE CONSTANTS FOR VAX 11/780
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***

!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /

!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /

!     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
!     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***

!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /

!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /

!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)

!     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
!     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
!     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
!     DATA LOG10(1), DIVER(2) / '3FD34413'X,'509F79FF'X /

!     MACHINE CONSTANTS FOR THE IBM PC - MICROSOFT FORTRAN

!     DATA SMALL(1), SMALL(2) / £00000000, £00100000 /
!     DATA LARGE(1), LARGE(2) / £FFFFFFFF, £7FEFFFFF /
!     DATA RIGHT(1), RIGHT(2) / £00000000, £3CA00000 /
!     DATA DIVER(1), DIVER(2) / £00000000, £3CB00000 /
!     DATA LOG10(1), LOG10(2) / £509F79FF, £3FD34413 /

!     MACHINE CONSTANTS FOR THE IBM PC - PROFESSIONAL FORTRAN
!     AND LAHEY FORTRAN

!      DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000' /
!      DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF' /
!      DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000' /
!      DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000' /
!      DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413' /
!
!!***FIRST EXECUTABLE STATEMENT  D1MACH
!      IF (I .LT. 1  .OR.  I .GT. 5)  &
!         CALL XERROR( 'D1MACH -- I OUT OF BOUNDS',25,1,2)
!
!      D1MACH = DMACH(I)
!      RETURN
!
!END FUNCTION D1MACH
!DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)
!***END PROLOGUE  D1MACH
!
      X = 1.0D0
      B = RADIX(X)
      SELECT CASE (I)
        CASE (1)
          D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
      RETURN
      END
! ---------------------------------------------------------------------
INTEGER FUNCTION I1MACH(I)
!***BEGIN PROLOGUE  I1MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  890313   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS)
!***PURPOSE  Return integer machine dependent constants.
!***DESCRIPTION

!     I1MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function
!     subroutine with one (input) argument, and can be called
!     as follows, for example

!          K = I1MACH(I)

!     where I=1,...,16.  The (output) value of K above is
!     determined by the (input) value of I.  The results for
!     various values of I are discussed below.

!  I/O unit numbers.
!    I1MACH( 1) = the standard input unit.
!    I1MACH( 2) = the standard output unit.
!    I1MACH( 3) = the standard punch unit.
!    I1MACH( 4) = the standard error message unit.

!  Words.
!    I1MACH( 5) = the number of bits per integer storage unit.
!    I1MACH( 6) = the number of characters per integer storage unit.

!  Integers.
!    assume integers are represented in the S-digit, base-A form

!               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )

!               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!    I1MACH( 7) = A, the base.
!    I1MACH( 8) = S, the number of base-A digits.
!    I1MACH( 9) = A**S - 1, the largest magnitude.

!  Floating-Point Numbers.
!    Assume floating-point numbers are represented in the T-digit,
!    base-B form
!               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )

!               where 0 .LE. X(I) .LT. B for I=1,...,T,
!               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, the base.

!  Single-Precision
!    I1MACH(11) = T, the number of base-B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.

!  Double-Precision
!    I1MACH(14) = T, the number of base-B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.

!  To alter this function for a particular environment,
!  the desired set of DATA statements should be activated by
!  removing the ! from column 1.  Also, the values of
!  I1MACH(1) - I1MACH(4) should be checked for consistency
!  with the local operating system.

!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  I1MACH

      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)

!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT COMPILER

!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    5 /
!     DATA IMACH(4) /    6 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) / 2147483647 /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -126 /
!     DATA IMACH(13)/  127 /
!     DATA IMACH(14)/   53 /
!     DATA IMACH(15)/ -1022 /
!     DATA IMACH(16)/  1023 /

!     MACHINE CONSTANTS FOR THE APOLLO

!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    6 /
!     DATA IMACH(4) /    6 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) / 2147483647 /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -125 /
!     DATA IMACH(13)/  129 /
!     DATA IMACH(14)/   53 /
!     DATA IMACH(15)/ -1021 /
!     DATA IMACH(16)/  1025 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM

!     DATA IMACH( 1) /    7 /
!     DATA IMACH( 2) /    2 /
!     DATA IMACH( 3) /    2 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -256 /
!     DATA IMACH(13) /  255 /
!     DATA IMACH(14) /   60 /
!     DATA IMACH(15) / -256 /
!     DATA IMACH(16) /  255 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM

!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  48 /
!     DATA IMACH( 6) /   6 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /   8 /
!     DATA IMACH(11) /  13 /
!     DATA IMACH(12) / -50 /
!     DATA IMACH(13) /  76 /
!     DATA IMACH(14) /  26 /
!     DATA IMACH(15) / -50 /
!     DATA IMACH(16) /  76 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS

!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  48 /
!     DATA IMACH( 6) /   6 /
!     DATA IMACH( 7) /   2 /
!     DATA IMACH( 8) /  39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /   8 /
!     DATA IMACH(11) /  13 /
!     DATA IMACH(12) / -50 /
!     DATA IMACH(13) /  76 /
!     DATA IMACH(14) /  26 /
!     DATA IMACH(15) / -32754 /
!     DATA IMACH(16) /  32780 /

!     MACHINE CONSTANTS FOR THE CD! 170/180 SERIES USING NOS/VE

!     DATA IMACH( 1) /     5 /
!     DATA IMACH( 2) /     6 /
!     DATA IMACH( 3) /     7 /
!     DATA IMACH( 4) /     6 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -4095 /
!     DATA IMACH(13) /  4094 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -4095 /
!     DATA IMACH(16) /  4094 /

!     MACHINE CONSTANTS FOR THE CD! 6000/7000 SERIES

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    7 /
!     DATA IMACH( 4) /6LOUTPUT/
!     DATA IMACH( 5) /   60 /
!     DATA IMACH( 6) /   10 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   47 /
!     DATA IMACH(12) / -929 /
!     DATA IMACH(13) / 1070 /
!     DATA IMACH(14) /   94 /
!     DATA IMACH(15) / -929 /
!     DATA IMACH(16) / 1069 /

!     MACHINE CONSTANTS FOR THE CELERITY C1260

!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    6 /
!     DATA IMACH(4) /    0 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) / Z'7FFFFFFF' /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -126 /
!     DATA IMACH(13)/  127 /
!     DATA IMACH(14)/   53 /
!     DATA IMACH(15)/ -1022 /
!     DATA IMACH(16)/  1023 /

!     MACHINE CONSTANTS FOR THE CONVEX C-1

!     DATA IMACH( 1) /     5/
!     DATA IMACH( 2) /     6/
!     DATA IMACH( 3) /     7/
!     DATA IMACH( 4) /     6/
!     DATA IMACH( 5) /    32/
!     DATA IMACH( 6) /     4/
!     DATA IMACH( 7) /     2/
!     DATA IMACH( 8) /    31/
!     DATA IMACH( 9) /2147483647/
!     DATA IMACH(10) /     2/
!     DATA IMACH(11) /    24/
!     DATA IMACH(12) /  -128/
!     DATA IMACH(13) /   127/
!     DATA IMACH(14) /    53/
!     DATA IMACH(15) / -1024/
!     DATA IMACH(16) /  1023/

!     MACHINE CONSTANTS FOR THE CRAY-1
!     USING THE 46 BIT INTEGER COMPILER OPTION

!     DATA IMACH( 1) /   100 /
!     DATA IMACH( 2) /   101 /
!     DATA IMACH( 3) /   102 /
!     DATA IMACH( 4) /   101 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    46 /
!     DATA IMACH( 9) /  1777777777777777B /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -8189 /
!     DATA IMACH(13) /  8190 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -8099 /
!     DATA IMACH(16) /  8190 /

!     MACHINE CONSTANTS FOR THE CRAY-1
!     USING THE 64 BIT INTEGER COMPILER OPTION

!     DATA IMACH( 1) /   100 /
!     DATA IMACH( 2) /   101 /
!     DATA IMACH( 3) /   102 /
!     DATA IMACH( 4) /   101 /
!     DATA IMACH( 5) /    64 /
!     DATA IMACH( 6) /     8 /
!     DATA IMACH( 7) /     2 /
!     DATA IMACH( 8) /    63 /
!     DATA IMACH( 9) /  777777777777777777777B /
!     DATA IMACH(10) /     2 /
!     DATA IMACH(11) /    47 /
!     DATA IMACH(12) / -8189 /
!     DATA IMACH(13) /  8190 /
!     DATA IMACH(14) /    94 /
!     DATA IMACH(15) / -8099 /
!     DATA IMACH(16) /  8190 /

!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200

!     DATA IMACH( 1) /   11 /
!     DATA IMACH( 2) /   12 /
!     DATA IMACH( 3) /    8 /
!     DATA IMACH( 4) /   10 /
!     DATA IMACH( 5) /   16 /
!     DATA IMACH( 6) /    2 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   15 /
!     DATA IMACH( 9) /32767 /
!     DATA IMACH(10) /   16 /
!     DATA IMACH(11) /    6 /
!     DATA IMACH(12) /  -64 /
!     DATA IMACH(13) /   63 /
!     DATA IMACH(14) /   14 /
!     DATA IMACH(15) /  -64 /
!     DATA IMACH(16) /   63 /

!     MACHINE CONSTANTS FOR THE ELXSI 6400

!     DATA IMACH( 1) /     5/
!     DATA IMACH( 2) /     6/
!     DATA IMACH( 3) /     6/
!     DATA IMACH( 4) /     6/
!     DATA IMACH( 5) /    32/
!     DATA IMACH( 6) /     4/
!     DATA IMACH( 7) /     2/
!     DATA IMACH( 8) /    32/
!     DATA IMACH( 9) /2147483647/
!     DATA IMACH(10) /     2/
!     DATA IMACH(11) /    24/
!     DATA IMACH(12) /  -126/
!     DATA IMACH(13) /   127/
!     DATA IMACH(14) /    53/
!     DATA IMACH(15) / -1022/
!     DATA IMACH(16) /  1023/

!     MACHINE CONSTANTS FOR THE HARRIS 220

!     DATA IMACH( 1) /       5 /
!     DATA IMACH( 2) /       6 /
!     DATA IMACH( 3) /       0 /
!     DATA IMACH( 4) /       6 /
!     DATA IMACH( 5) /      24 /
!     DATA IMACH( 6) /       3 /
!     DATA IMACH( 7) /       2 /
!     DATA IMACH( 8) /      23 /
!     DATA IMACH( 9) / 8388607 /
!     DATA IMACH(10) /       2 /
!     DATA IMACH(11) /      23 /
!     DATA IMACH(12) /    -127 /
!     DATA IMACH(13) /     127 /
!     DATA IMACH(14) /      38 /
!     DATA IMACH(15) /    -127 /
!     DATA IMACH(16) /     127 /

!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /   43 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    6 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   63 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /

!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4

!     DATA IMACH(1) /      5/
!     DATA IMACH(2) /      6 /
!     DATA IMACH(3) /      4 /
!     DATA IMACH(4) /      1 /
!     DATA IMACH(5) /     16 /
!     DATA IMACH(6) /      2 /
!     DATA IMACH(7) /      2 /
!     DATA IMACH(8) /     15 /
!     DATA IMACH(9) /  32767 /
!     DATA IMACH(10)/      2 /
!     DATA IMACH(11)/     23 /
!     DATA IMACH(12)/   -128 /
!     DATA IMACH(13)/    127 /
!     DATA IMACH(14)/     39 /
!     DATA IMACH(15)/   -128 /
!     DATA IMACH(16)/    127 /

!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4

!     DATA IMACH(1) /      5 /
!     DATA IMACH(2) /      6 /
!     DATA IMACH(3) /      4 /
!     DATA IMACH(4) /      1 /
!     DATA IMACH(5) /     16 /
!     DATA IMACH(6) /      2 /
!     DATA IMACH(7) /      2 /
!     DATA IMACH(8) /     15 /
!     DATA IMACH(9) /  32767 /
!     DATA IMACH(10)/      2 /
!     DATA IMACH(11)/     23 /
!     DATA IMACH(12)/   -128 /
!     DATA IMACH(13)/    127 /
!     DATA IMACH(14)/     55 /
!     DATA IMACH(15)/   -128 /
!     DATA IMACH(16)/    127 /

!     MACHINE CONSTANTS FOR THE HP 9000

!     DATA IMACH(1)  /    5 /
!     DATA IMACH(2)  /    6 /
!     DATA IMACH(3)  /    6 /
!     DATA IMACH(3)  /    7 /
!     DATA IMACH(5)  /   32 /
!     DATA IMACH(6)  /    4 /
!     DATA IMACH(7)  /    2 /
!     DATA IMACH(8)  /   32 /
!     DATA IMACH(9)  /2147483647 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -126 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   53 /
!     DATA IMACH(15) /-1015 /
!     DATA IMACH(16) / 1017 /

!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.

!     DATA IMACH( 1) /   5 /
!     DATA IMACH( 2) /   6 /
!     DATA IMACH( 3) /   7 /
!     DATA IMACH( 4) /   6 /
!     DATA IMACH( 5) /  32 /
!     DATA IMACH( 6) /   4 /
!     DATA IMACH( 7) /  16 /
!     DATA IMACH( 8) /  31 /
!     DATA IMACH( 9) / Z7FFFFFFF /
!     DATA IMACH(10) /  16 /
!     DATA IMACH(11) /   6 /
!     DATA IMACH(12) / -64 /
!     DATA IMACH(13) /  63 /
!     DATA IMACH(14) /  14 /
!     DATA IMACH(15) / -64 /
!     DATA IMACH(16) /  63 /

!     MACHINE CONSTANTS FOR THE IBM PC

      DATA IMACH( 1) /     5 /
      DATA IMACH( 2) /     6 /
      DATA IMACH( 3) /     0 /
      DATA IMACH( 4) /     0 /
      DATA IMACH( 5) /    32 /
      DATA IMACH( 6) /     4 /
      DATA IMACH( 7) /     2 /
      DATA IMACH( 8) /    31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /     2 /
      DATA IMACH(11) /    24 /
      DATA IMACH(12) /  -125 /
      DATA IMACH(13) /   127 /
      DATA IMACH(14) /    53 /
      DATA IMACH(15) / -1021 /
      DATA IMACH(16) /  1023 /

!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    5 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   54 /
!     DATA IMACH(15) / -101 /
!     DATA IMACH(16) /  127 /

!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    5 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   62 /
!     DATA IMACH(15) / -128 /
!     DATA IMACH(16) /  127 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETIC.

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   32 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   56 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGER ARITHMETIC.

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    5 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   16 /
!     DATA IMACH( 6) /    2 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   15 /
!     DATA IMACH( 9) / 32767 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   24 /
!     DATA IMACH(12) / -127 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   56 /
!     DATA IMACH(15) / -127 /
!     DATA IMACH(16) /  127 /

!     MACHINE CONSTANTS FOR THE SUN

!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    6 /
!     DATA IMACH(4) /    6 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) /2147483647 /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -125 /
!     DATA IMACH(13)/  128 /
!     DATA IMACH(14)/   53 /
!     DATA IMACH(15)/ -1021 /
!     DATA IMACH(16)/  1024 /

!     MACHINE CONSTANTS FOR THE UNIVA! 1100 SERIES FTN COMPILER

!     DATA IMACH( 1) /    5 /
!     DATA IMACH( 2) /    6 /
!     DATA IMACH( 3) /    1 /
!     DATA IMACH( 4) /    6 /
!     DATA IMACH( 5) /   36 /
!     DATA IMACH( 6) /    4 /
!     DATA IMACH( 7) /    2 /
!     DATA IMACH( 8) /   35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /    2 /
!     DATA IMACH(11) /   27 /
!     DATA IMACH(12) / -128 /
!     DATA IMACH(13) /  127 /
!     DATA IMACH(14) /   60 /
!     DATA IMACH(15) /-1024 /
!     DATA IMACH(16) / 1023 /

!     MACHINE CONSTANTS FOR THE VAX 11/780

!     DATA IMACH(1) /    5 /
!     DATA IMACH(2) /    6 /
!     DATA IMACH(3) /    5 /
!     DATA IMACH(4) /    6 /
!     DATA IMACH(5) /   32 /
!     DATA IMACH(6) /    4 /
!     DATA IMACH(7) /    2 /
!     DATA IMACH(8) /   31 /
!     DATA IMACH(9) /2147483647 /
!     DATA IMACH(10)/    2 /
!     DATA IMACH(11)/   24 /
!     DATA IMACH(12)/ -127 /
!     DATA IMACH(13)/  127 /
!     DATA IMACH(14)/   56 /
!     DATA IMACH(15)/ -127 /
!     DATA IMACH(16)/  127 /

!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR

!     DATA IMACH( 1) /     1/
!     DATA IMACH( 2) /     1/
!     DATA IMACH( 3) /     0/
!     DATA IMACH( 4) /     1/
!     DATA IMACH( 5) /    16/
!     DATA IMACH( 6) /     2/
!     DATA IMACH( 7) /     2/
!     DATA IMACH( 8) /    15/
!     DATA IMACH( 9) / 32767/
!     DATA IMACH(10) /     2/
!     DATA IMACH(11) /    24/
!     DATA IMACH(12) /  -127/
!     DATA IMACH(13) /   127/
!     DATA IMACH(14) /    56/
!     DATA IMACH(15) /  -127/
!     DATA IMACH(16) /   127/

!***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10

      I1MACH = IMACH(I)
      RETURN

   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')

!     CALL FDUMP

      STOP
END FUNCTION I1MACH
! ---------------------------------------------------------------------
SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
!***BEGIN PROLOGUE  XERROR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERROR-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Process an error (diagnostic) message.
!***DESCRIPTION

!     Abstract
!        XERROR processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)

!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed, containing
!                no more than 72 characters.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.

!     Examples
!        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
!        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
!    1                43,2,1)
!        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
!    1ULLY COLLAPSED.',65,3,0)
!        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)

!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  XERRWV
!***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
END SUBROUTINE XERROR
! ---------------------------------------------------------------------
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
!***BEGIN PROLOGUE  XERRWV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  890531   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERRWV-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Process an error message allowing 2 integer and 2 real
!            values to be included in the message.
!***DESCRIPTION

!     Abstract
!        XERRWV processes a diagnostic message, in a manner
!        determined by the value of LEVEL and the current value
!        of the library error control flag, KONTRL.
!        (See subroutine XSETF for details.)
!        In addition, up to two integer values and two real
!        values may be printed along with the message.

!     Description of Parameters
!      --Input--
!        MESSG - the Hollerith message to be processed.
!        NMESSG- the actual number of characters in MESSG.
!        NERR  - the error number associated with this message.
!                NERR must not be zero.
!        LEVEL - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (I.e., it is
!                   non-fatal if XSETF has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!        NI    - number of integer values to be printed. (0 to 2)
!        I1    - first integer value.
!        I2    - second integer value.
!        NR    - number of real values to be printed. (0 to 2)
!        R1    - first real value.
!        R2    - second real value.

!     Examples
!        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
!    1   1,NUM,0,0,0.,0.)
!        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
!    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)

!     Latest revision ---  1 August 1985
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
!                    XGETUA
!***END PROLOGUE  XERRWV

! ----------------------------------------------------------------------

!  Change record:
!     89-05-31   Changed all specific intrinsics to generic.  (WRB)

! ----------------------------------------------------------------------

      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
!     GET FLAGS
!***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
!     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.  &
          (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.', &
         29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
!     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
!     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
!     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX(-2,MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN(1,MAXMES)))  &
      .OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))         &
      .OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))  &
      .OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
!           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT  &
      ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT  &
            ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
!        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(REAL(I1MACH(9))) + 1.0
         ISIZEF = LOG10(REAL(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',  &
               I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
!              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
!        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))  &
      IFATAL = 1
!     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX(1,MAXMES))) GO TO 120
!        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT  &
         ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT  &
         ('JOB ABORT DUE TO FATAL ERROR.',29)
!        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
!     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
END SUBROUTINE XERRWV
! ---------------------------------------------------------------------
FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
!***BEGIN PROLOGUE  J4SAVE
!***REFER TO  XERROR
!***ROUTINES CALLED  (NONE)
!***DESCRIPTION

!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.

!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                 = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                 = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                 = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                 = 6 Refers to the 2nd unit for error messages
!                 = 7 Refers to the 3rd unit for error messages
!                 = 8 Refers to the 4th unit for error messages
!                 = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.

!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!    Adapted from Bell Laboratories PORT Library Error Handler
!     Latest revision ---  1 August 1985
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END FUNCTION J4SAVE
! ---------------------------------------------------------------------
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
!***BEGIN PROLOGUE  XERSAV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERSAV-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Record that an error has occurred.
!***DESCRIPTION

!     Abstract
!        Record that this error occurred.

!     Description of Parameters
!     --Input--
!       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
!       except that when NMESSG=0 the tables will be
!       dumped and cleared, and when NMESSG is less than zero the
!       tables will be dumped and not cleared.
!     --Output--
!       ICOUNT will be the number of times this message has
!       been seen, or zero if the table has overflowed and
!       does not contain this message specifically.
!       When NMESSG=0, ICOUNT will not be altered.

!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  1 August 1985
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I1MACH,XGETUA
!***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
!     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
!     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),  &
           KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)  &
           /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
!***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
!     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
!        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/  &
            51H MESSAGE START             NERR     LEVEL     COUNT)
!           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
!           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
!        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
!     PROCESS A MESSAGE...
!     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
!     THREE POSSIBLE CASES...
!     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
!     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
!     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
END SUBROUTINE XERSAV
! ---------------------------------------------------------------------
SUBROUTINE XGETUA(IUNITA,N)
!***BEGIN PROLOGUE  XGETUA
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XGETUA-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Return unit number(s) to which error messages are being
!            sent.
!***DESCRIPTION

!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.

!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.

!     Latest revision ---  19 MAR 1980
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
END SUBROUTINE XGETUA
! ---------------------------------------------------------------------
SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
!***BEGIN PROLOGUE  XERCTL
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERCTL-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Allow user control over handling of errors.
!***DESCRIPTION

!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCTL.
!        If the user has provided his own version of XERCTL, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.

!     Description of Parameters

!      --Input--
!        MESSG1 - the first word (only) of the error message.
!        NMESSG - same as in the call to XERROR or XERRWV.
!        NERR   - same as in the call to XERROR or XERRWV.
!        LEVEL  - same as in the call to XERROR or XERRWV.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.

!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
!***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END SUBROUTINE XERCTL
! ---------------------------------------------------------------------
      SUBROUTINE XERPRT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERPRT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  890531   (YYMMDD)
!***CATEGORY NO.  R3
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERPRT-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Print error messages.
!***DESCRIPTION

!     Abstract
!        Print the Hollerith message in MESSG, of length NMESSG,
!        on each file indicated by XGETUA.
!     Latest revision ---  1 August 1985
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  I1MACH,XGETUA
!***END PROLOGUE  XERPRT

! ----------------------------------------------------------------------

!  Change record:
!     89-05-31   Changed all specific intrinsics to generic.  (WRB)

! ----------------------------------------------------------------------

      INTEGER LUN(5)
      CHARACTER*(*) MESSG
!     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
!***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
END SUBROUTINE XERPRT
! ---------------------------------------------------------------------
SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(FDUMP-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Symbolic dump (should be locally written).
!***DESCRIPTION

!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.

!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END SUBROUTINE FDUMP
! ---------------------------------------------------------------------
      SUBROUTINE XERABT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERABT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERABT-A),ERROR
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Abort program execution and print error message.
!***DESCRIPTION

!     Abstract
!        ***Note*** machine dependent routine
!        XERABT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.

!     Description of Parameters
!        MESSG and NMESSG are as in XERROR, except that NMESSG may
!        be zero, in which case no message is being supplied.

!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!     Latest revision ---  1 August 1982
!***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
!                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
!                 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
END SUBROUTINE XERABT
! ---------------------------------------------------------------------
END MODULE UTILIT
! end of file Utilit.f90
