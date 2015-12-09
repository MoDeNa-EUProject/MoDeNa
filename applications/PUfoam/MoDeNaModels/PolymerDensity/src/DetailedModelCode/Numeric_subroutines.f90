!>WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!! This file contains subroutines for subtasks like spline interpolations,
!! spline integration and function minimization
!!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW





 SUBROUTINE bicub_derivative ( ya, x1a, x2a, y1a, y2a, y12a, i_max, k_max )
!
 USE DFT_MODULE, ONLY: NDFT, r_grid
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
 REAL, INTENT(IN)                       :: x1a(r_grid)
 REAL, INTENT(IN)                       :: x2a(NDFT)
 REAL, INTENT(OUT)                      :: y1a(r_grid,NDFT)
 REAL, INTENT(OUT)                      :: y2a(r_grid,NDFT)
 REAL, INTENT(OUT)                      :: y12a(r_grid,NDFT)
 INTEGER, INTENT(IN)                    :: i_max
 INTEGER, INTENT(IN)                    :: k_max
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, k
! ----------------------------------------------------------------------


DO i = 2, i_max-1
  DO k = 2, k_max-1
    y1a(i,k) = (ya(i+1,k)-ya(i-1,k)) / (x1a(i+1)-x1a(i-1))
    y2a(i,k) = (ya(i,k+1)-ya(i,k-1)) / (x2a(k+1)-x2a(k-1))
    y12a(i,k)= (ya(i+1,k+1)-ya(i+1,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
               /((x1a(i+1)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
  END DO
END DO

i = 1
DO k = 1, k_max
  y1a(i,k) = 0.0
  y2a(i,k) = 0.0
  y12a(i,k)= 0.0
END DO

k = 1
DO i = 1, i_max
  y1a(i,k) = 0.0
  y2a(i,k) = 0.0
  y12a(i,k)= 0.0
END DO


i = i_max
DO k = 2, k_max-1
  y1a(i,k) = (ya(i,k)-ya(i-1,k)) / (x1a(i)-x1a(i-1))
  y2a(i,k) = (ya(i,k+1)-ya(i,k-1)) / (x2a(k+1)-x2a(k-1))
  y12a(i,k)= (ya(i,k+1)-ya(i,k-1)-ya(i-1,k+1)+ya(i-1,k-1))  &
             /((x1a(i)-x1a(i-1))*(x2a(k+1)-x2a(k-1)))
END DO


k = k_max
DO i = 2, i_max-1
  y1a(i,k) = (ya(i+1,k)-ya(i-1,k)) / (x1a(i+1)-x1a(i-1))
  y2a(i,k) = (ya(i,k)-ya(i,k-1)) / (x2a(k)-x2a(k-1))
  y12a(i,k)= (ya(i+1,k)-ya(i+1,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
             /((x1a(i+1)-x1a(i-1))*(x2a(k)-x2a(k-1)))
END DO

k = k_max
i = i_max
y1a(i,k) = (ya(i,k)-ya(i-1,k)) / (x1a(i)-x1a(i-1))
y2a(i,k) = (ya(i,k)-ya(i,k-1)) / (x2a(k)-x2a(k-1))
y12a(i,k)= (ya(i,k)-ya(i,k-1)-ya(i-1,k)+ya(i-1,k-1))  &
           /((x1a(i)-x1a(i-1))*(x2a(k)-x2a(k-1)))

END SUBROUTINE bicub_derivative


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE bicub_c ( ya, x1a, x2a, y1a, y2a, y12a, c_bicub, i_max, k_max )
!
 USE DFT_MODULE, ONLY: NDFT, r_grid
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
 REAL, INTENT(IN)                       :: x1a(r_grid)
 REAL, INTENT(IN)                       :: x2a(NDFT)
 REAL, INTENT(IN)                       :: y1a(r_grid,NDFT)
 REAL, INTENT(IN)                       :: y2a(r_grid,NDFT)
 REAL, INTENT(IN)                       :: y12a(r_grid,NDFT)
 REAL, INTENT(OUT)                      :: c_bicub(r_grid,NDFT,4,4)
 INTEGER, INTENT(IN)                    :: i_max
 INTEGER, INTENT(IN)                    :: k_max
!
! ----------------------------------------------------------------------
 INTEGER                                :: i, k, m, n
 REAL                                   :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
 REAL                                   :: c(4,4)
! ----------------------------------------------------------------------

DO i = 1, i_max-1
  DO k = 1, k_max-1
    y(1)=ya(i,k)
    y(2)=ya(i+1,k)
    y(3)=ya(i+1,k+1)
    y(4)=ya(i,k+1)
    
    y1(1)=y1a(i,k)
    y1(2)=y1a(i+1,k)
    y1(3)=y1a(i+1,k+1)
    y1(4)=y1a(i,k+1)
    
    y2(1)=y2a(i,k)
    y2(2)=y2a(i+1,k)
    y2(3)=y2a(i+1,k+1)
    y2(4)=y2a(i,k+1)
    
    y12(1)=y12a(i,k)
    y12(2)=y12a(i+1,k)
    y12(3)=y12a(i+1,k+1)
    y12(4)=y12a(i,k+1)
    
    x1l=x1a(i)
    x1u=x1a(i+1)
    x2l=x2a(k)
    x2u=x2a(k+1)
    
    CALL bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
    DO m=1,4
      DO n=1,4
        c_bicub(i,k,m,n)=c(m,n)
      END DO
    END DO
    
  END DO
END DO

END SUBROUTINE bicub_c

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN)                       :: y(4)
 REAL, INTENT(IN)                       :: y1(4)
 REAL, INTENT(IN)                       :: y2(4)
 REAL, INTENT(IN)                       :: y12(4)
 REAL, INTENT(IN)                       :: d1
 REAL, INTENT(IN)                       :: d2
 REAL, INTENT(OUT)                      :: c(4,4)
!
! ----------------------------------------------------------------------
 INTEGER                                :: i,j,k,l
 REAL                                   :: d1d2,xx,cl(16),wt(16,16),x(16)
 SAVE wt
 DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*  &
    0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,  &
    1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,  &
    -6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,  &
    10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,  &
    -2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,  &
    2,-2,2*0,-1,1/
! ----------------------------------------------------------------------

d1d2 = d1 * d2
DO  i = 1, 4
  x(i)    = y(i)
  x(i+4)  = y1(i)*d1
  x(i+8)  = y2(i)*d2
  x(i+12) = y12(i)*d1d2
END DO
DO i = 1, 16
  xx = 0.0
  DO  k = 1, 16
    xx = xx + wt(i,k) * x(k)
  END DO
  cl(i) = xx
END DO
l = 0
DO i = 1, 4
  DO j = 1, 4
    l = l + 1
    c(i,j) = cl(l)
  END DO
END DO

END SUBROUTINE bcucof


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE BI_CUB_SPLINE ( rho_rdf, xg, ya, x1a, x2a, y1a, y2a, y12a,  &
                            c_bicub, rdf, dg_drho, dg_dr, i_max, ih, k )
!
 USE DFT_MODULE, ONLY: NDFT, r_grid
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: rho_rdf
 REAL, INTENT(IN OUT)                   :: xg
 REAL, INTENT(IN)                       :: ya(r_grid,NDFT)
 REAL, INTENT(IN)                       :: x1a(r_grid)
 REAL, INTENT(IN)                       :: x2a(NDFT)
 REAL, INTENT(IN)                       :: y1a(r_grid,NDFT)
 REAL, INTENT(IN)                       :: y2a(r_grid,NDFT)
 REAL, INTENT(IN)                       :: y12a(r_grid,NDFT)
 REAL, INTENT(IN)                       :: c_bicub(r_grid,NDFT,4,4)
 REAL, INTENT(OUT)                      :: rdf
 REAL, INTENT(OUT)                      :: dg_drho
 REAL, INTENT(OUT)                      :: dg_dr
 INTEGER, INTENT(IN OUT)                :: i_max
 !INTEGER, INTENT(IN OUT)                :: k_max
 INTEGER, INTENT(OUT)                   :: ih
 INTEGER, INTENT(IN)                    :: k
!
! ----------------------------------------------------------------------
 INTEGER                                :: m, n

 REAL                                   :: y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u
 REAL                                   :: c(4,4)
! ----------------------------------------------------------------------

IF ( rho_rdf < x1a(1) ) THEN
  dg_drho = 0.0
  dg_dr = 0.0
  rdf = 1.0
  RETURN
END IF
IF ( x1a(ih) <= rho_rdf .AND. rho_rdf < x1a(ih+1) ) GO TO 10
IF ( ih > 2 ) THEN
  IF ( x1a(ih-1) <= rho_rdf .AND. rho_rdf < x1a(ih) ) THEN
    ih = ih - 1
    GO TO 10
  END IF
END IF
! write (*,*) 'in ',ih
CALL hunt ( x1a, i_max, rho_rdf, ih )
! write (*,*) 'out',ih
10   CONTINUE
IF ( x2a(k) /= xg ) THEN
! write (*,*) 'error bi-cubic-spline',k,x2a(k),xg
! DO k=1,NDFT
!   write (*,*) k,x2a(k)
! ENDDO
! stop
END IF



y(1) = ya(ih,k)
y(2) = ya(ih+1,k)
y(3) = ya(ih+1,k+1)
y(4) = ya(ih,k+1)

y1(1) = y1a(ih,k)
y1(2) = y1a(ih+1,k)
y1(3) = y1a(ih+1,k+1)
y1(4) = y1a(ih,k+1)

y2(1) = y2a(ih,k)
y2(2) = y2a(ih+1,k)
y2(3) = y2a(ih+1,k+1)
y2(4) = y2a(ih,k+1)

y12(1) = y12a(ih,k)
y12(2) = y12a(ih+1,k)
y12(3) = y12a(ih+1,k+1)
y12(4) = y12a(ih,k+1)

x1l = x1a(ih)
x1u = x1a(ih+1)
x2l = x2a(k)
x2u = x2a(k+1)

DO m = 1, 4
  DO n = 1, 4
    c(m,n) = c_bicub( ih, k, m, n )
  END DO
END DO
CALL bcuint ( x1l, x1u, x2l, x2u, rho_rdf, xg, c, rdf, dg_drho, dg_dr )

END SUBROUTINE BI_CUB_SPLINE

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE hunt
!
! Given an array xx(1:n), and given a value x, returns a value jlo
! such that x is between xx(jlo) and xx(jlo+1). xx(1:n) must be
! monotonic, either increasing or decreasing. jlo=0 or jlo=n is
! returned to indicate that x is out of range. jlo on input is taken
! as the initial guess for jlo on output.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE hunt ( xx, n, x, jlo )
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: n
 INTEGER, INTENT(OUT)                   :: jlo
 REAL, INTENT(IN)                       :: xx(n)
 REAL                                   :: x
!
! ----------------------------------------------------------------------
 INTEGER                                :: inc,jhi,jm
 LOGICAL                                :: ascnd
! ----------------------------------------------------------------------

ascnd = xx(n) >= xx(1)
IF( jlo <= 0 .OR. jlo > n ) THEN
  jlo = 0
  jhi = n + 1
  GO TO 3
END IF
inc = 1
IF( x >= xx(jlo) .EQV. ascnd ) THEN
1 jhi = jlo + inc
  IF ( jhi > n ) THEN
    jhi = n + 1
  ELSE IF ( x >= xx(jhi) .EQV. ascnd ) THEN
    jlo = jhi
    inc = inc + inc
    GO TO 1
  END IF
ELSE
  jhi = jlo
2 jlo = jhi - inc
  IF ( jlo < 1 ) THEN
    jlo = 0
  ELSE IF ( x < xx(jlo) .EQV. ascnd ) THEN
    jhi = jlo
    inc = inc + inc
    GO TO 2
  END IF
END IF
3  IF (jhi-jlo == 1 ) THEN
  IF ( x == xx(n)) jlo = n - 1
  IF ( x == xx(1) ) jlo = 1
  RETURN
END IF
jm = ( jhi + jlo ) / 2
IF ( x >= xx(jm) .EQV. ascnd ) THEN
  jlo = jm
ELSE
  jhi = jm
END IF
GO TO 3
END SUBROUTINE hunt



!**********************************************************************
!
!**********************************************************************
!
 !SUBROUTINE bcuint ( y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, c, ansy, ansy1, ansy2 )
 SUBROUTINE bcuint ( x1l, x1u, x2l, x2u, x1, x2, c, ansy, ansy1, ansy2 )
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 !REAL, INTENT(IN OUT)                   :: y(4)
 !REAL, INTENT(IN OUT)                   :: y1(4)
 !REAL, INTENT(IN OUT)                   :: y2(4)
 !REAL, INTENT(IN OUT)                   :: y12(4)
 REAL, INTENT(IN OUT)                   :: x1l
 REAL, INTENT(IN OUT)                   :: x1u
 REAL, INTENT(IN OUT)                   :: x2l
 REAL, INTENT(IN OUT)                   :: x2u
 REAL, INTENT(IN OUT)                   :: x1
 REAL, INTENT(IN OUT)                   :: x2
 REAL, INTENT(IN)                       :: c(4,4)
 REAL, INTENT(OUT)                      :: ansy
 REAL, INTENT(OUT)                      :: ansy1
 REAL, INTENT(OUT)                      :: ansy2
!
! ----------------------------------------------------------------------
 !U    USES bcucof
 INTEGER                                :: i
 REAL                                   :: t, u
! ----------------------------------------------------------------------

! call bcucof ( y, y1, y2, y12, x1u-x1l, x2u-x2l, c )

IF ( x1u == x1l .OR. x2u == x2l ) PAUSE 'bad input in bcuint'
t = (x1-x1l) / (x1u-x1l)
u = (x2-x2l) / (x2u-x2l)
ansy  = 0.0
ansy2 = 0.0
ansy1 = 0.0
DO  i = 4, 1, -1
  ansy  = t *ansy  + ( (c(i,4)*u + c(i,3))*u+c(i,2) )*u + c(i,1)
  ansy2 = t *ansy2 + ( 3.*c(i,4)*u+2.*c(i,3) )*u + c(i,2)
  ansy1 = u *ansy1 + ( 3.*c(4,i)*t+2.*c(3,i) )*t + c(2,i)
END DO
ansy1 = ansy1 / (x1u-x1l)
ansy2 = ansy2 / (x2u-x2l)

END SUBROUTINE bcuint


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE spline
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function,
! i.e., yi = f(xi), with  x1 < x2 < .. . < xN, and given values yp1 and
! ypn for the first derivative of the interpolating function at points 1
! and n, respectively, this routine returns an array y2(1:n) of length n
! which contains the second derivatives of the interpolating function at
! the tabulated points xi. If yp1 and/or ypn are equal to 1  1030 or
! larger, the routine is signaled to set  the corresponding boundary
! condition for a natural spline, with zero second derivative on that
! boundary.
! Parameter: NMAX is the largest anticipated value of n.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE spline ( x, y, n, yp1, ypn, y2 )
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN)                       :: x(n)
 REAL, INTENT(IN)                       :: y(n)
 REAL, INTENT(IN)                       :: yp1
 REAL, INTENT(IN)                       :: ypn
 REAL, INTENT(OUT)                      :: y2(n)
!
! ----------------------------------------------------------------------
 INTEGER, PARAMETER                     :: NMAX = 500
 INTEGER                                :: i, k
 REAL                                   :: p, qn, sig, un, u(NMAX)
! ----------------------------------------------------------------------

 IF ( yp1 > 0.99E30 ) THEN
    y2(1) = 0.0
    u(1)  = 0.0
 ELSE
    y2(1) = -0.5
    u(1)  = ( 3.0/(x(2)-x(1)) ) * ( (y(2)-y(1))/(x(2)-x(1))-yp1 )
 END IF
 DO  i = 2, n-1
    IF ( (x(i+1)-x(i)) == 0.0 .OR. (x(i)-x(i-1)) == 0.0 .OR. (x(i+1)-x(i-1)) == 0.0 ) THEN
       write (*,*) 'error in spline-interpolation'
       stop
    END IF
    sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    p = sig*y2(i-1) + 2.0
    y2(i) = (sig-1.0) / p
    u(i) = ( 6.0 * ((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1))  &
         - sig * u(i-1) ) / p
 END DO
 IF ( ypn > 0.99E30 ) THEN
    qn = 0.0
    un = 0.0
 ELSE
    qn = 0.5
    un = (3.0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
 END IF
 y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.0)
 DO  k = n-1, 1, -1
    y2(k) = y2(k) * y2(k+1) + u(k)
 END DO

END SUBROUTINE spline

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SUBROUTINE splint_integral
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the in order), and given the array y2a(1:n), which is
! the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE splint_integral ( xa, ya, y2a, n, xlo, xhi, integral )
!
 IMPLICIT NONE
!
! ----------------------------------------------------------------------
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN)                       :: xa(n)
 REAL, INTENT(IN)                       :: ya(n)
 REAL, INTENT(IN)                       :: y2a(n)
 REAL, INTENT(IN)                       :: xlo
 REAL, INTENT(IN)                       :: xhi
 REAL, INTENT(OUT)                      :: integral
!
! ----------------------------------------------------------------------
 INTEGER                                :: k, khi_l, klo_l, khi_h, klo_h
 REAL                                   :: xl, xh, h, INT, x0, x1, y0, y1, y20, y21
! ----------------------------------------------------------------------

 integral = 0.0
 klo_l = 1
 khi_l = n
1 IF ( khi_l-klo_l > 1 ) THEN
    k = ( khi_l + klo_l ) / 2
    IF ( xa(k) > xlo ) THEN
       khi_l = k
    ELSE
       klo_l = k
    END IF
    GO TO 1
 END IF

 klo_h = 1
 khi_h = n
2 IF ( khi_h-klo_h > 1 ) THEN
    k = ( khi_h + klo_h ) / 2
    IF ( xa(k) > xhi ) THEN
       khi_h = k
    ELSE
       klo_h = k
    END IF
    GO TO 2
 END IF

 ! integration in spline pieces, the lower interval, bracketed
 ! by xa(klo_L) and xa(khi_L) is in steps shifted upward.

 ! first: determine upper integration bound
 xl = xlo
3 CONTINUE
 IF ( khi_h > khi_l ) THEN
    xh = xa(khi_l)
 ELSE IF ( khi_h == khi_l ) THEN
    xh = xhi
 ELSE
    WRITE (*,*) 'error in spline-integration'
    PAUSE
 END IF

 h = xa(khi_l) - xa(klo_l)
 IF ( h == 0.0 ) PAUSE 'bad xa input in splint'
 x0 = xa(klo_l)
 x1 = xa(khi_l)
 y0 = ya(klo_l)
 y1 = ya(khi_l)
 y20= y2a(klo_l)
 y21= y2a(khi_l)
 ! int = -xL/h * ( (x1-.5*xL)*y0 + (0.5*xL-x0)*y1  &
 !            +y20/6.*(x1**3-1.5*xL*x1*x1+xL*xL*x1-.25*xL**3)  &
 !            -y20/6.*h*h*(x1-.5*xL)  &
 !            +y21/6.*(.25*xL**3-xL*xL*x0+1.5*xL*x0*x0-x0**3)  &
 !            -y21/6.*h*h*(.5*xL-x0) )
 ! int = int + xH/h * ( (x1-.5*xH)*y0 + (0.5*xH-x0)*y1  &
 !            +y20/6.*(x1**3-1.5*xH*x1*x1+xH*xH*x1-.25*xH**3)  &
 !            -y20/6.*h*h*(x1-.5*xH)  &
 !            +y21/6.*(.25*xH**3-xH*xH*x0+1.5*xH*x0*x0-x0**3)  &
 !            -y21/6.*h*h*(.5*xH-x0) )
 INT = -1.0/h * ( xl*((x1-.5*xl)*y0 + (0.5*xl-x0)*y1)  &
      -y20/24.*(x1-xl)**4 + y20/6.*(0.5*xl*xl-x1*xl)*h*h  &
      +y21/24.*(xl-x0)**4 - y21/6.*(0.5*xl*xl-x0*xl)*h*h )
 INT = INT + 1.0/h * ( xh*((x1-.5*xh)*y0 + (0.5*xh-x0)*y1)  &
      -y20/24.*(x1-xh)**4 + y20/6.*(0.5*xh*xh-x1*xh)*h*h  &
      +y21/24.*(xh-x0)**4 - y21/6.*(0.5*xh*xh-x0*xh)*h*h )

 integral = integral + INT
 ! write (*,*) integral,x1,xH
 klo_l = klo_l + 1
 khi_l = khi_l + 1
 xl = x1
 IF (khi_h /= (khi_l-1)) GO TO 3    ! the -1 in (khi_L-1) because khi_L was already counted up

END SUBROUTINE splint_integral




!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 REAL FUNCTION praxis( t0, machep, h0, n, prin, x, f, fmin )

 IMPLICIT NONE

! ----------------------------------------------------------------------
 REAL, INTENT(IN OUT)                   :: t0
 REAL, INTENT(IN)                       :: machep
 REAL, INTENT(IN)                       :: h0
 INTEGER                                :: n
 INTEGER, INTENT(IN OUT)                :: prin
 REAL, INTENT(IN OUT)                   :: x(n)
 REAL                                   :: f
 REAL, INTENT(IN OUT)                   :: fmin
! ----------------------------------------------------------------------

EXTERNAL f

!     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
!     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
!     NOT REQUIRED.

!     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
!     "ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
!     CALCULATING DERIVATIVES" BY RICHARD P BRENT.

!     THE PARAMETERS ARE:
!     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
!              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
!              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
!     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
!              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
!              2.22D-16) FOR REAL*8 ARITHMETIC ON THE IBM 360.
!     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
!              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
!              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
!              CONVERGENCE MAY BE SLOW.)
!     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
!              THE FUNCTION DEPENDS.
!     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
!              IF PRIN=0, NOTHING IS PRINTED.
!              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
!              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
!              PRINTED ONLY IF N IS AT MOST 4.
!              IF PRIN=2, THE SCALE FACTORS AND THE PRINCIPAL VALUES OF
!              THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
!              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
!              MINIMIZATIONS.
!              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
!              QUADRATIC FORM ARE ALSO PRINTED.
!     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
!              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
!     F(X,N)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE A REAL*8
!              FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
!     FMIN     IS AN ESTIMATE OF THE MINIMUM, USED ONLY IN PRINTING
!              INTERMEDIATE RESULTS.
!     THE APPROXIMATING QUADRATIC FORM IS
!              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
!     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
!              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
!     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
!     OF SECOND DIFFERENCES).  IF F HAS CONTINUOUS SECOND DERIVATIVES
!     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.

!     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
!     TO ZERO.
!     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
!     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.

 LOGICAL :: illc
 INTEGER :: nl,nf,kl,kt,ktm,idim,i,j,k,k2,km1,klmk,ii,im1
 REAL :: s,sl,dn,dmin,fx,f1,lds,ldt,t,h,sf,df,qf1,qd0, qd1,qa,qb,qc
 REAL :: m2,m4,small,vsmall,large,vlarge,scbd,ldfac,t2, dni,value
 REAL :: random

!.....IF N>20 OR IF N<20 AND YOU NEED MORE SPACE, CHANGE '20' TO THE
!     LARGEST VALUE OF N IN THE NEXT CARD, IN THE CARD 'IDIM=20', AND
!     IN THE DIMENSION STATEMENTS IN SUBROUTINES MINFIT,MIN,FLIN,QUAD.

 REAL :: d(20),y(20),z(20),q0(20),q1(20),v(20,20)
 COMMON /global/ fx,ldt,dmin,nf,nl /q/ v,q0,q1,qa,qb,qc,qd0,qd1,qf1

! ---------------------------------
!     introduced by Joachim........
 idim = n
! ---------------------------------



!.....INITIALIZATION.....
!     MACHINE DEPENDENT NUMBERS:

small = machep*machep
vsmall = small*small
large = 1.d0/small
vlarge = 1.d0/vsmall
m2 = SQRT(machep)
m4 = SQRT(m2)

!     HEURISTIC NUMBERS:
!     IF THE AXES MAY BE BADLY SCALED (WHICH IS TO BE AVOIDED IF
!     POSSIBLE), THEN SET SCBD=10.  OTHERWISE SET SCBD=1.
!     IF THE PROBLEM IS KNOWN TO BE ILL-CONDITIONED, SET ILLC=TRUE.
!     OTHERWISE SET ILLC=FALSE.
!     KTM IS THE NUMBER OF ITERATIONS WITHOUT IMPROVEMENT BEFORE THE
!     ALGORITHM TERMINATES.  KTM=4 IS VERY CAUTIOUS; USUALLY KTM=1
!     IS SATISFACTORY.

scbd = 1.0
illc = .false.
ktm = 1

ldfac = 0.01
IF (illc) ldfac = 0.1
kt = 0
nl = 0
nf = 1
fx = f(x,n)
qf1 = fx
t = small+ABS(t0)
t2 = t
dmin = small
h = h0
IF (h < 100*t) h = 100*t
ldt = h
!.....THE FIRST SET OF SEARCH DIRECTIONS V IS THE IDENTITY MATRIX.....
DO  i = 1,n
  DO  j = 1,n
    v(i,j) = 0.0
  END DO
  v(i,i) = 1.0
END DO
d(1) = 0.0
qd0 = 0.0
DO  i = 1,n
  q0(i) = x(i)
  q1(i) = x(i)
END DO
IF (prin > 0) CALL PRINT(n,x,prin,fmin)

!.....THE MAIN LOOP STARTS HERE.....
40    sf=d(1)
d(1)=0.d0
s=0.d0

!.....MINIMIZE ALONG THE FIRST DIRECTION V(*,1).
!     FX MUST BE PASSED TO MIN BY VALUE.
value=fx
CALL MIN(n,1,2,d(1),s,value,.false.,f,x,t,machep,h)
IF (s > 0.d0) GO TO 50
DO  i=1,n
  v(i,1)=-v(i,1)
END DO
50    IF (sf > 0.9D0*d(1).AND.0.9D0*sf < d(1)) GO TO 70
DO  i=2,n
  d(i)=0.d0
END DO

!.....THE INNER LOOP STARTS HERE.....
70    DO  k=2,n
  DO  i=1,n
    y(i)=x(i)
  END DO
  sf=fx
  IF (kt > 0) illc=.true.
  80         kl=k
  df=0.d0
  
!.....A RANDOM STEP FOLLOWS (TO AVOID RESOLUTION VALLEYS).
!     PRAXIS ASSUMES THAT RANDOM RETURNS A RANDOM NUMBER UNIFORMLY
!     DISTRIBUTED IN (0,1).
  
  IF(.NOT.illc) GO TO 95
  DO  i=1,n
    s=(0.1D0*ldt+t2*(10**kt))*(random(n)-0.5D0)
    z(i)=s
    DO  j=1,n
      x(j)=x(j)+s*v(j,i)
    END DO
  END DO
  fx=f(x,n)
  nf=nf+1
  
!.....MINIMIZE ALONG THE "NON-CONJUGATE" DIRECTIONS V(*,K),...,V(*,N)
  
  95         DO  k2=k,n
    sl=fx
    s=0.d0
    value=fx
    CALL MIN(n,k2,2,d(k2),s,value,.false.,f,x,t,machep,h)
    IF (illc) GO TO 97
    s=sl-fx
    GO TO 99
    97              s=d(k2)*((s+z(k2))**2)
    99              IF (df > s) CYCLE
    df=s
    kl=k2
  END DO
  IF (illc.OR.(df >= ABS((100*machep)*fx))) GO TO 110
  
!.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
!     ILLC=TRUE AND START THE INNER LOOP AGAIN.....
  
  illc=.true.
  GO TO 80
  110        IF (k == 2.AND.prin > 1) CALL vcprnt(1,d,n)
  
!.....MINIMIZE ALONG THE "CONJUGATE" DIRECTIONS V(*,1),...,V(*,K-1)
  
  km1=k-1
  DO  k2=1,km1
    s=0
    value=fx
    CALL MIN(n,k2,2,d(k2),s,value,.false.,f,x,t,machep,h)
  END DO
  f1=fx
  fx=sf
  lds=0
  DO  i=1,n
    sl=x(i)
    x(i)=y(i)
    sl=sl-y(i)
    y(i)=sl
    lds=lds+sl*sl
  END DO
  lds=SQRT(lds)
  IF (lds <= small) GO TO 160
  
!.....DISCARD DIRECTION V(*,KL).
!     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE "NON-CONJUGATE"
!     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
  
  klmk=kl-k
  IF (klmk < 1) GO TO 141
  DO  ii=1,klmk
    i=kl-ii
    DO  j=1,n
      v(j,i+1)=v(j,i)
    END DO
    d(i+1)=d(i)
  END DO
  141        d(k)=0
  DO  i=1,n
    v(i,k)=y(i)/lds
  END DO
  
!.....MINIMIZE ALONG THE NEW "CONJUGATE" DIRECTION V(*,K), WHICH IS
!     THE NORMALIZED VECTOR:  (NEW X) - (0LD X).....
  
  value=f1
  CALL MIN(n,k,4,d(k),lds,value,.true.,f,x,t,machep,h)
  IF (lds > 0.d0) GO TO 160
  lds=-lds
  DO  i=1,n
    v(i,k)=-v(i,k)
  END DO
  160        ldt=ldfac*ldt
  IF (ldt < lds) ldt=lds
  IF (prin > 0) CALL PRINT(n,x,prin,fmin)
  t2=0.d0
  DO  i=1,n
    t2=t2+x(i)**2
  END DO
  t2=m2*SQRT(t2)+t
  
!.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
!     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
  
  IF (ldt > (0.5*t2)) kt=-1
  kt=kt+1
  IF (kt > ktm) GO TO 400
END DO
!.....THE INNER LOOP ENDS HERE.

!     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.

CALL quad(n,f,x,t,machep,h)
dn=0.d0
DO  i=1,n
  d(i)=1.d0/SQRT(d(i))
  IF (dn < d(i)) dn=d(i)
END DO
IF (prin > 3) CALL maprnt(1,v,idim,n)
DO  j=1,n
  s=d(j)/dn
  DO  i=1,n
    v(i,j)=s*v(i,j)
  END DO
END DO

!.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....

IF (scbd <= 1.d0) GO TO 200
s=vlarge
DO  i=1,n
  sl=0.d0
  DO  j=1,n
    sl=sl+v(i,j)*v(i,j)
  END DO
  z(i)=SQRT(sl)
  IF (z(i) < m4) z(i)=m4
  IF (s > z(i)) s=z(i)
END DO
DO  i=1,n
  sl=s/z(i)
  z(i)=1.d0/sl
  IF (z(i) <= scbd) GO TO 189
  sl=1.d0/scbd
  z(i)=scbd
  189             DO  j=1,n
    v(i,j)=sl*v(i,j)
  END DO
END DO

!.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
!     THE MAIN LOOP.
!     FIRST TRANSPOSE V FOR MINFIT:

200   DO  i=2,n
  im1=i-1
  DO  j=1,im1
    s=v(i,j)
    v(i,j)=v(j,i)
    v(j,i)=s
  END DO
END DO

!.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
!     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
!     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
!     NUMBER.....

CALL minfit(idim,n,machep,vsmall,v,d)

!.....UNSCALE THE AXES.....

IF (scbd <= 1.d0) GO TO 250
DO  i=1,n
  s=z(i)
  DO  j=1,n
    v(i,j)=s*v(i,j)
  END DO
END DO
DO  i=1,n
  s=0.d0
  DO  j=1,n
    s=s+v(j,i)**2
  END DO
  s=SQRT(s)
  d(i)=s*d(i)
  s=1/s
  DO  j=1,n
    v(j,i)=s*v(j,i)
  END DO
END DO

250   DO  i=1,n
  dni=dn*d(i)
  IF (dni > large) GO TO 265
  IF (dni < small) GO TO 260
  d(i)=1/(dni*dni)
  CYCLE
  260             d(i)=vlarge
  CYCLE
  265        d(i)=vsmall
END DO

!.....SORT THE EIGENVALUES AND EIGENVECTORS.....

CALL sort(idim,n,d,v)
dmin=d(n)
IF (dmin < small) dmin=small
illc=.false.
IF (m2*d(1) > dmin) illc=.true.
IF (prin > 1.AND.scbd > 1.d0) CALL vcprnt(2,z,n)
IF (prin > 1) CALL vcprnt(3,d,n)
IF (prin > 3) CALL maprnt(2,v,idim,n)
!.....THE MAIN LOOP ENDS HERE.....

GO TO 40

!.....RETURN.....

400   IF (prin > 0) CALL vcprnt(4,x,n)
praxis=fx

END FUNCTION praxis


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE minfit(m,n,machep,tol,ab,q)

 IMPLICIT NONE

 INTEGER, INTENT(IN OUT)                :: m
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN)                       :: machep
 REAL, INTENT(IN OUT)                   :: tol
 REAL, INTENT(IN OUT)                   :: ab(m,n)
 REAL, INTENT(OUT)                      :: q(n)
 INTEGER                                :: i,j,k,l, kk,kt,l2,ll2,ii,lp1
!      IMPLICIT REAL (A-H,O-Z)


REAL                                    :: x,eps,e(20),g,s, f,h,y,c,z,temp
!...AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
!   RESTRICTED TO M=N,P=0.
!   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
!   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
!   WHERE U IS ANOTHER ORTHOGONAL MATRIX.

!...HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM...
IF (n == 1) GO TO 200
eps = machep
g = 0.d0
x = 0.d0
DO  i=1,n
  e(i) = g
  s = 0.d0
  l = i + 1
  DO  j=i,n
    s = s + ab(j,i)**2
  END DO
  g = 0.d0
  IF (s < tol) GO TO 4
  f = ab(i,i)
  g = SQRT(s)
  IF (f >= 0.d0) g = -g
  h = f*g - s
  ab(i,i)=f-g
  IF (l > n) GO TO 4
  DO  j=l,n
    f = 0.d0
    DO  k=i,n
      f = f + ab(k,i)*ab(k,j)
    END DO
    f = f/h
    DO  k=i,n
      ab(k,j) = ab(k,j) + f*ab(k,i)
    END DO
  END DO
  4        q(i) = g
  s = 0.d0
  IF (i == n) GO TO 6
  DO  j=l,n
    s = s + ab(i,j)*ab(i,j)
  END DO
  6        g = 0.d0
  IF (s < tol) GO TO 10
  IF (i == n) GO TO 16
  f = ab(i,i+1)
  16          g = SQRT(s)
  IF (f >= 0.d0) g = -g
  h = f*g - s
  IF (i == n) GO TO 10
  ab(i,i+1) = f - g
  DO  j=l,n
    e(j) = ab(i,j)/h
  END DO
  DO  j=l,n
    s = 0.d0
    DO  k=l,n
      s = s + ab(j,k)*ab(i,k)
    END DO
    DO  k=l,n
      ab(j,k) = ab(j,k) + s*e(k)
    END DO
  END DO
  10       y = ABS(q(i)) + ABS(e(i))
  IF (y > x) x = y
END DO
!...ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS...
ab(n,n) = 1.d0
g = e(n)
l = n
DO  ii=2,n
  i = n - ii + 1
  IF (g == 0.d0) GO TO 23
  h = ab(i,i+1)*g
  DO  j=l,n
    ab(j,i) = ab(i,j)/h
  END DO
  DO  j=l,n
    s = 0.d0
    DO  k=l,n
      s = s + ab(i,k)*ab(k,j)
    END DO
    DO  k=l,n
      ab(k,j) = ab(k,j) + s*ab(k,i)
    END DO
  END DO
  23       DO  j=l,n
    ab(i,j) = 0.d0
    ab(j,i) = 0.d0
  END DO
  ab(i,i) = 1.d0
  g = e(i)
  l = i
END DO
!...DIAGONALIZATION OF THE BIDIAGONAL FORM...
eps = eps*x
DO  kk=1,n
  k = n - kk + 1
  kt = 0
  101      kt = kt + 1
  IF (kt <= 30) GO TO 102
  e(k) = 0.d0
  WRITE (6,1000)
  1000        FORMAT (' QR FAILED')
  102      DO  ll2=1,k
    l2 = k - ll2 + 1
    l = l2
    IF (ABS(e(l)) <= eps) GO TO 120
    IF (l == 1) CYCLE
    IF (ABS(q(l-1)) <= eps) EXIT
  END DO
!...CANCELLATION OF E(L) IF L>1...
  c = 0.d0
  s = 1.d0
  DO  i=l,k
    f = s*e(i)
    e(i) = c*e(i)
    IF (ABS(f) <= eps) GO TO 120
    g = q(i)
!...Q(I) = H = SQRT(G*G + F*F)...
    IF (ABS(f) < ABS(g)) GO TO 113
    IF (f == 0.0) THEN
      GO TO   111
    ELSE
      GO TO   112
    END IF
    111         h = 0.d0
    GO TO 114
    112         h = ABS(f)*SQRT(1 + (g/f)**2)
    GO TO 114
    113         h = ABS(g)*SQRT(1 + (f/g)**2)
    114         q(i) = h
    IF (h /= 0.d0) GO TO 115
    g = 1.d0
    h = 1.d0
    115         c = g/h
    s = -f/h
  END DO
!...TEST FOR CONVERGENCE...
  120      z = q(k)
  IF (l == k) GO TO 140
!...SHIFT FROM BOTTOM 2*2 MINOR...
  x = q(l)
  y = q(k-1)
  g = e(k-1)
  h = e(k)
  f = ((y - z)*(y + z) + (g - h)*(g + h))/(2*h*y)
  g = SQRT(f*f + 1.0D0)
  temp = f - g
  IF (f >= 0.d0) temp = f + g
  f = ((x - z)*(x + z) + h*(y/temp - h))/x
!...NEXT QR TRANSFORMATION...
  c = 1.d0
  s = 1.d0
  lp1 = l + 1
  IF (lp1 > k) GO TO 133
  DO  i=lp1,k
    g = e(i)
    y = q(i)
    h = s*g
    g = g*c
    IF (ABS(f) < ABS(h)) GO TO 123
    IF (f == 0.0) THEN
      GO TO   121
    ELSE
      GO TO   122
    END IF
    121         z = 0.d0
    GO TO 124
    122         z = ABS(f)*SQRT(1 + (h/f)**2)
    GO TO 124
    123         z = ABS(h)*SQRT(1 + (f/h)**2)
    124         e(i-1) = z
    IF (z /= 0.d0) GO TO 125
    f = 1.d0
    z = 1.d0
    125         c = f/z
    s = h/z
    f = x*c + g*s
    g = -x*s + g*c
    h = y*s
    y = y*c
    DO  j=1,n
      x = ab(j,i-1)
      z = ab(j,i)
      ab(j,i-1) = x*c + z*s
      ab(j,i) = -x*s + z*c
    END DO
    IF (ABS(f) < ABS(h)) GO TO 129
    IF (f == 0.0) THEN
      GO TO   127
    ELSE
      GO TO   128
    END IF
    127         z = 0.d0
    GO TO 130
    128         z = ABS(f)*SQRT(1 + (h/f)**2)
    GO TO 130
    129         z = ABS(h)*SQRT(1 + (f/h)**2)
    130         q(i-1) = z
    IF (z /= 0.d0) GO TO 131
    f = 1.d0
    z = 1.d0
    131         c = f/z
    s = h/z
    f = c*g + s*y
    x = -s*g + c*y
  END DO
  133      e(l) = 0.d0
  e(k) = f
  q(k) = x
  GO TO 101
!...CONVERGENCE:  Q(K) IS MADE NON-NEGATIVE...
  140      IF (z >= 0.d0) CYCLE
  q(k) = -z
  DO  j=1,n
    ab(j,k) = -ab(j,k)
  END DO
END DO
RETURN
200   q(1) = ab(1,1)
ab(1,1) = 1.d0

END SUBROUTINE minfit

 


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE MIN(n,j,nits,d2,x1,f1,fk,f,x,t,machep,h)

 IMPLICIT NONE

 INTEGER, INTENT(IN)                    :: n
 INTEGER                                :: j
 INTEGER                                :: nits
 REAL, INTENT(IN OUT)                   :: d2
 REAL, INTENT(IN OUT)                   :: x1
 REAL, INTENT(IN OUT)                   :: f1
 LOGICAL                                :: fk
 REAL                                   :: f
 REAL, INTENT(IN OUT)                   :: x(n)
 REAL, INTENT(IN)                       :: t
 REAL, INTENT(IN)                       :: machep
 REAL, INTENT(IN)                       :: h

 INTEGER                                :: i,k
 EXTERNAL f


 REAL                                   :: flin   ! function
 REAL                                   :: small,sf1,sx1,s,temp, xm,x2,f2,d1
 REAL                                   :: fm,f0,t2
!----------------------------------------------
 INTEGER                                :: nf,nl
 REAL                                   :: fx,ldt,dmin
 REAL :: v(20,20),q0(20),q1(20),qa,qb,qc,qd0,qd1,qf1
 COMMON /global/ fx,ldt,dmin,nf,nl /q/ v,q0,q1,qa,qb,qc,qd0,qd1,qf1
!----------------------------------------------
!...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
!   J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
!   DEFINED BY Q0,Q1,X.
!   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
!   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
!   ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
!   FOUND.
!   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
!   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
!   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
!   THE INTERVAL. 
 LOGICAL                                :: dz
 REAL                                   :: m2,m4

small = machep**2
m2 = SQRT(machep)
m4 = SQRT(m2)
sf1 = f1
sx1 = x1
k = 0
xm = 0.d0
fm = fx
f0 = fx
dz = d2 < machep
!...FIND THE STEP SIZE...
s = 0.d0
DO  i=1,n
  s = s + x(i)**2
END DO
s = SQRT(s)
temp = d2
IF (dz) temp = dmin
t2 = m4*SQRT(ABS(fx)/temp + s*ldt) + m2*ldt
s = m4*s + t
IF (dz.AND.t2 > s) t2 = s
t2 = DMAX1(t2,small)
t2 = DMIN1(t2,.01D0*h)
IF (.NOT.fk.OR.f1 > fm) GO TO 2
xm = x1
fm = f1
2     IF (fk.AND.ABS(x1) >= t2) GO TO 3
temp=1.d0
IF (x1 < 0.d0) temp=-1.d0
x1=temp*t2
f1 = flin(n,j,x1,f,x,nf)
3     IF (f1 > fm) GO TO 4
xm = x1
fm = f1
4     IF (.NOT.dz) GO TO 6
!...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
x2 = -x1
IF (f0 >= f1) x2 = 2.d0*x1
f2 = flin(n,j,x2,f,x,nf)
IF (f2 > fm) GO TO 5
xm = x2
fm = f2
5     d2 = (x2*(f1 - f0)-x1*(f2 - f0))/((x1*x2)*(x1 - x2))
!...ESTIMATE THE FIRST DERIVATIVE AT 0...
6     d1 = (f1 - f0)/x1 - x1*d2
dz = .true.
!...PREDICT THE MINIMUM...
IF (d2 > small) GO TO 7
x2 = h
IF (d1 >= 0.d0) x2 = -x2
GO TO 8
7        x2 = (-.5D0*d1)/d2
8     IF (ABS(x2) <= h) GO TO 11
IF (x2 > 0.0) THEN
  GO TO    10
END IF
x2 = -h
GO TO 11
10       x2 = h
!...EVALUATE F AT THE PREDICTED MINIMUM...
11    f2 = flin(n,j,x2,f,x,nf)
IF (k >= nits.OR.f2 <= f0) GO TO 12
!...NO SUCCESS, SO TRY AGAIN...
k = k + 1
IF (f0 < f1.AND.(x1*x2) > 0.d0) GO TO 4
x2 = 0.5D0*x2
GO TO 11
!...INCREMENT THE ONE-DIMENSIONAL SEARCH COUNTER...
12    nl = nl + 1
IF (f2 <= fm) GO TO 13
x2 = xm
GO TO 14
13    fm = f2
!...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
14    IF (ABS(x2*(x2 - x1)) <= small) GO TO 15
d2 = (x2*(f1-f0) - x1*(fm-f0))/((x1*x2)*(x1 - x2))
GO TO 16
15       IF (k > 0) d2 = 0.d0
16    IF (d2 <= small) d2 = small
x1 = x2
fx = fm
IF (sf1 >= fx) GO TO 17
fx = sf1
x1 = sx1
!...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
17    IF (j == 0) RETURN
DO  i=1,n
  x(i) = x(i) + x1*v(i,j)
END DO

END SUBROUTINE MIN


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE vcprnt(option,v,n)

 IMPLICIT NONE
 INTEGER                                :: option
 REAL, INTENT(IN OUT)                   :: v(n)
 INTEGER                                :: n

 INTEGER                                :: i

SELECT CASE ( option )
  CASE (    1)
    GO TO 1
  CASE (    2)
    GO TO 2
  CASE (    3)
    GO TO 3
  CASE (    4)
    GO TO 4
END SELECT

1     WRITE (6,101) (v(i),i=1,n)
RETURN
2     WRITE (6,102) (v(i),i=1,n)
RETURN
3     WRITE (6,103) (v(i),i=1,n)
RETURN
4     WRITE (6,104) (v(i),i=1,n)
RETURN
101   FORMAT (/' THE SECOND DIFFERENCE ARRAY D(*) IS:'/ (e32.14,4E25.14))
102   FORMAT (/' THE SCALE FACTORS ARE:'/(e32.14,4E25.14))
103   FORMAT (/' THE APPROXIMATING QUADR. FORM HAS PRINCIPAL VALUES:'/  &
    (e32.14,4E25.14))
104   FORMAT (/' X IS:',e26.14/(e32.14))
END SUBROUTINE vcprnt


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE PRINT(n,x,prin,fmin)

 IMPLICIT NONE
 INTEGER                                :: n
 REAL, INTENT(IN OUT)                   :: x(n)
 INTEGER, INTENT(IN OUT)                :: prin
 REAL, INTENT(IN OUT)                   :: fmin

 INTEGER                                :: i
 REAL                                   :: ln
!----------------------------------------------
INTEGER :: nf,nl
REAL :: fx,ldt,dmin
COMMON /global/ fx,ldt,dmin,nf,nl
!----------------------------------------------
WRITE (6,101) nl,nf,fx

IF (fx <= fmin) GO TO 1
ln = LOG10(fx-fmin)
WRITE (6,102) fmin,ln
GO TO 2
1     WRITE (6,103) fmin
2     IF (n > 4.AND.prin <= 2) RETURN
WRITE (6,104) (x(i),i=1,n)
RETURN
101   FORMAT (/' AFTER',i6,  &
    ' LINEAR SEARCHES, THE FUNCTION HAS BEEN EVALUATED',i6,  &
    ' TIMES.  THE SMALLEST VALUE FOUND IS F(X) = ',e21.14)
102   FORMAT (' LOG (F(X)-',e21.14,') = ',e21.14)
103   FORMAT (' LOG (F(X)-',e21.14,') IS UNDEFINED.')
104   FORMAT (' X IS:',e26.14/(e32.14))
END SUBROUTINE PRINT


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE maprnt(option,v,m,n)

 IMPLICIT NONE
 INTEGER                                :: option
 REAL, INTENT(IN OUT)                   :: v(m,n)
 INTEGER, INTENT(IN OUT)                :: m
 INTEGER, INTENT(IN)                    :: n
 INTEGER                                :: i,j

 INTEGER                                :: low,upp
!...THE SUBROUTINE MAPRNT PRINTS THE COLUMNS OF THE NXN MATRIX V
!   WITH A HEADING AS SPECIFIED BY OPTION.
!   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM...
low = 1
upp = 5
SELECT CASE ( option )
  CASE (    1)
    GO TO 1
  CASE (    2)
    GO TO 2
END SELECT
1    WRITE (6,101)
101  FORMAT (/' THE NEW DIRECTIONS ARE:')
GO TO 3
2    WRITE (6,102)
102  FORMAT (' AND THE PRINCIPAL AXES:')
3    IF (n < upp) upp = n
DO  i=1,n
  WRITE (6,104) (v(i,j),j=low,upp)
END DO
low = low + 5
IF (n < low) RETURN
upp = upp + 5
WRITE (6,103)
GO TO 3
103  FORMAT (' ')
104  FORMAT (e32.14,4E25.14)
END SUBROUTINE maprnt


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 REAL FUNCTION random(naught)

 IMPLICIT NONE
 INTEGER, INTENT(IN OUT)                :: naught

 REAL                                   :: ran1,ran3(127),half
 INTEGER                                :: i,j,ran2,q,r
 LOGICAL :: init
 DATA init/.false./
 SAVE init,ran2,ran1,ran3

IF (init) GO TO 3
r = MOD(naught,8190) + 1
ran2 = 128
DO  i=1,127
  ran2 = ran2 - 1
  ran1 = -2.d0**55
  DO  j=1,7
    r = MOD(1756*r,8191)
    q = r/32
    ran1 = (ran1 + q)*(1.0D0/256)
  END DO
  ran3(ran2) = ran1
END DO
init = .true.
3     IF (ran2 == 1) ran2 = 128
ran2 = ran2 - 1
ran1 = ran1 + ran3(ran2)
half = .5D0
IF (ran1 >= 0.d0) half = -half
ran1 = ran1 + half
ran3(ran2) = ran1
random = ran1 + .5D0

END FUNCTION random



!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 REAL FUNCTION flin (n,j,l,f,x,nf)
 IMPLICIT NONE

 INTEGER, INTENT(IN)                    :: n
 INTEGER, INTENT(IN OUT)                :: j
 REAL, INTENT(IN)                       :: l
 REAL                                   :: f
 REAL, INTENT(IN)                       :: x(n)
 INTEGER, INTENT(OUT)                   :: nf

 INTEGER                                :: i
 REAL                                   :: t(20)

 EXTERNAL f

!...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
!   BY THE SUBROUTINE MIN...
!----------------------------------------------
 REAL :: v(20,20),q0(20),q1(20),qa,qb,qc,qd0,qd1,qf1
 COMMON /q/ v,q0,q1,qa,qb,qc,qd0,qd1,qf1
!----------------------------------------------
IF (j == 0) GO TO 2
!...THE SEARCH IS LINEAR...
DO  i=1,n
  t(i) = x(i) + l*v(i,j)
END DO
GO TO 4
!...THE SEARCH IS ALONG A PARABOLIC SPACE CURVE...
2     qa = (l*(l - qd1))/(qd0*(qd0 + qd1))
qb = ((l + qd0)*(qd1 - l))/(qd0*qd1)
qc = (l*(l + qd0))/(qd1*(qd0 + qd1))
DO  i=1,n
  t(i) = (qa*q0(i) + qb*x(i)) + qc*q1(i)
END DO
!...THE FUNCTION EVALUATION COUNTER NF IS INCREMENTED...
4     nf = nf + 1
flin = f(t,n)

END FUNCTION flin


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE sort(m,n,d,v)
 IMPLICIT NONE
!
 INTEGER, INTENT(IN OUT)                :: m
 INTEGER, INTENT(IN)                    :: n
 REAL, INTENT(IN OUT)                   :: d(n)
 REAL, INTENT(IN OUT)                   :: v(m,n)

 INTEGER                                :: i,j,k,nm1,ip1
 REAL                                   :: s
!...SORTS THE ELEMENTS OF D(N) INTO DESCENDING ORDER AND MOVES THE
!   CORRESPONDING COLUMNS OF V(N,N).
!   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM.
IF (n == 1) RETURN
nm1 = n - 1
DO  i = 1,nm1
  k=i
  s = d(i)
  ip1 = i + 1
  DO  j = ip1,n
    IF (d(j) <= s) CYCLE
    k = j
    s = d(j)
  END DO
  IF (k <= i) CYCLE
  d(k) = d(i)
  d(i) = s
  DO  j = 1,n
    s = v(j,i)
    v(j,i) = v(j,k)
    v(j,k) = s
  END DO
END DO
END SUBROUTINE sort


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!
 SUBROUTINE quad(n,f,x,t,machep,h)
 IMPLICIT NONE

 INTEGER, INTENT(IN)                    :: n
 REAL                                   :: f
 REAL, INTENT(IN OUT)                   :: x(n)
 REAL, INTENT(IN OUT)                   :: t
 REAL                                   :: machep
 REAL, INTENT(IN OUT)                   :: h
!      IMPLICIT REAL (A-H,O-Z)
 EXTERNAL f

!...QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X...
 INTEGER                                :: i
 REAL                                   :: l
 REAL                                   :: s,value
!----------------------------------------------
 INTEGER                                :: nf,nl
 REAL                                   :: fx,ldt,dmin

REAL :: v(20,20),q0(20),q1(20),qa,qb,qc,qd0,qd1,qf1
COMMON /global/ fx,ldt,dmin,nf,nl /q/ v,q0,q1,qa,qb,qc,qd0,qd1,qf1
!----------------------------------------------
s = fx
fx = qf1
qf1 = s
qd1 = 0.d0
DO  i=1,n
  s = x(i)
  l = q1(i)
  x(i) = l
  q1(i) = s
  qd1 = qd1 + (s-l)**2
END DO
qd1 = SQRT(qd1)
l = qd1
s = 0.d0
IF (qd0 <= 0.d0 .OR. qd1 <= 0.d0 .OR. nl < 3*n*n) GO TO 2
value=qf1
CALL MIN(n,0,2,s,l,value,.true.,f,x,t,machep,h)
qa = (l*(l-qd1))/(qd0*(qd0+qd1))
qb = ((l+qd0)*(qd1-l))/(qd0*qd1)
qc = (l*(l+qd0))/(qd1*(qd0+qd1))
GO TO 3
2    fx = qf1
qa = 0.d0
qb = qa
qc = 1.d0
3    qd0 = qd1
DO  i=1,n
  s = q0(i)
  q0(i) = x(i)
  x(i) = (qa*s + qb*x(i)) + qc*q1(i)
END DO
END SUBROUTINE quad

