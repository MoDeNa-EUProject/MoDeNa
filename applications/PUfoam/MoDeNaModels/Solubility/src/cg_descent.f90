!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!      ________________________________________________________________
!     |      A conjugate gradient method with guaranteed descent       |
!     |                                                                |
!     |             Version 1.1  (December 10, 2004)                   |
!     |             Version 1.2  (June 4, 2005)                        |
!     |             Version 1.3  (October 6, 2005)                     |
!     |             Version 1.4  (November 14, 2005)                   |
!     |                                                                |
!     |           William W. Hager    and   Hongchao Zhang             |
!     |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
!     |                   Department of Mathematics                    |
!     |                     University of Florida                      |
!     |                 Gainesville, Florida 32611 USA                 |
!     |                      352-392-0281 x 244                        |
!     |                                                                |
!     |              Copyright 2004 by William W. Hager                |
!     |                                                                |
!     |This program is free software; you can redistribute it and/or   |
!     |modify it under the terms of the GNU General Public License as  |
!     |published by the Free Software Foundation; either version 2 of  |
!     |the License, or (at your option) any later version.             |
!     |This program is distributed in the hope that it will be useful, |
!     |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
!     |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
!     |GNU General Public License for more details.                    |
!     |                                                                |
!     |You should have received a copy of the GNU General Public       |
!     |License along with this program; if not, write to the Free      |
!     |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
!     |MA  02110-1301  USA                                             |
!     |                                                                |
!     |          http://www.math.ufl.edu/~hager/papers/CG              |
!     |                                                                |
!     |    INPUT:                                                      |
!     |                                                                |
!     |(double) grad_tol-- StopRule = T: |g|_infty <= max (grad_tol,   |
!     |                          StopFac*initial |g|_infty) [default]  |
!     |                    StopRule = F: |g|_infty <= grad_tol(1+|f|)  |
!     |                                                                |
!     |(double) x       --starting guess (length n)                    |
!     |                                                                |
!     |(int)    dim     --problem dimension (also denoted n)           |
!     |                                                                |
!     |         cg_value--name of cost evaluation subroutine           |
!     |                  (external in main program, cg_value(f, x, n)  |
!     |                   puts value of cost function at x in f        |
!     |                   f is double precision scalar and x is        |
!     |                   double precision array of length n)          |
!     |                                                                |
!     |         cg_grad --name gradient evaluation subroutine          |
!     |                  (external in main program, cg_grad (g, x, n)  |
!     |                   puts gradient at x in g, g and x are         |
!     |                   double precision arrays of length n)         |
!     |                                                                |
!     |(double) gnorm   --if the parameter Step in cg_descent.parm is  |
!     |                   .true., then gnorm contains the initial step |
!     |                   used at iteration 0 in the line search       |
!     |                                                                |
!     |(double) d       --direction (work array, length n)             |
!     |                                                                |
!     |(double) g       --gradient (work array, length n)              |
!     |                                                                |
!     |(double) xtemp   --work array (work array, length n)            |
!     |                                                                |
!     |(double) gtemp   --work array (work array, length n)            |
!     |                                                                |
!     |    OUTPUT:                                                     |
!     |                                                                |
!     |(int)    status  -- 0 (convergence tolerance satisfied)         |
!     |                    1 (change in func <= feps*|f|)              |
!     |                    2 (total iterations exceeded maxit)         |
!     |                    3 (slope always negative in line search)    |
!     |                    4 (number secant iterations exceed nsecant) |
!     |                    5 (search direction not a descent direction)|
!     |                    6 (line search fails in initial interval)   |
!     |                    7 (line search fails during bisection)      |
!     |                    8 (line search fails during interval update)|
!     |                                                                |
!     |(double) gnorm   --max abs component of gradient                |
!     |                                                                |
!     |(double) f       --function value at solution                   |
!     |                                                                |
!     |(double) x       --solution (length n)                          |
!     |                                                                |
!     |(int)    iter    --number of iterations                         |
!     |                                                                |
!     |(int)    nfunc   --number of function evaluations               |
!     |                                                                |
!     |(int)    ngrad   --number of gradient evaluations               |
!     |                                                                |
!     |Note: The file cg_descent.parm must be placed in the directory  |
!     |      where the code is run                                     |
!     |________________________________________________________________|
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

MODULE cg_minimization


  implicit none
  save

  integer                               :: n, n5, n6, nf, ng
  integer                               :: info
  integer                               :: nrestart, nexpand, nsecant
  integer                               :: maxit
  real                                  :: delta, sigma, eps
  real                                  :: gamma, rho, tol
  real                                  :: eta, fpert, f0
  real                                  :: ck, qdecay
  real                                  :: wolfe_hi, wolfe_lo, awolfe_hi
  real                                  :: quadcutoff, stopfac, awolfefac
  real                                  :: zero, feps
  real                                  :: psi0, psi1, psi2
  logical                               :: pertrule, quadok
  logical                               :: quadstep, printlevel
  logical                               :: printfinal, stoprule
  logical                               :: awolfe
  logical                               :: step
  logical                               :: debug

  real                                  :: f_best
  real, dimension(:), allocatable       :: x_best

  PRIVATE
  PUBLIC  :: cg_descent

CONTAINS


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE cg_descent (grad_tol, x, dim, cg_value, cg_grad,  &
       STATUS, gnorm, f, iter, nfunc, ngrad, d, g, xtemp, gtemp)

    implicit none


    INTEGER, INTENT(IN)                      :: dim        ! (= n) dimensionality of problem
    REAL, INTENT(IN)                         :: grad_tol
    REAL, INTENT(IN OUT)                     :: x(dim)       ! initial guess / for output: solution (length n)
    INTEGER, INTENT(OUT)                     :: STATUS
    REAL, INTENT(IN OUT)                     :: gnorm      ! max abs component of gradient
    REAL, INTENT(IN OUT)                     :: f          ! function value at solution
    INTEGER, INTENT(OUT)                     :: iter       ! number of iterations
    INTEGER, INTENT(OUT)                     :: nfunc      ! number of function evaluations
    INTEGER, INTENT(OUT)                     :: ngrad      ! number of gradient evaluations
    REAL, INTENT(IN OUT)                     :: d(dim)       ! direction (work array, length n)
    REAL, INTENT(IN OUT)                     :: g(dim)       ! gradient (work array, length n)
    REAL, INTENT(IN OUT)                     :: xtemp(dim)   ! work array (length n)
    REAL, INTENT(IN OUT)                     :: gtemp(dim)   ! work array (length n)

    INTERFACE

       SUBROUTINE cg_value (f, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: f
       END SUBROUTINE cg_value

       SUBROUTINE cg_grad (g, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: g(:)
       END SUBROUTINE cg_grad

    END INTERFACE

    integer                               :: i, j, i1, i2, i3, i4
    real                                  :: delta2, eta_sq, qk
    real                                  :: ftemp, xnorm, gnorm2
    real                                  :: dnorm2, denom
    real                                  :: t, t1, t2, t3, t4
    real                                  :: dphi, dphi0
    real                                  :: alpha, talpha
    real                                  :: yk, yk2, ykgk, dkyk
    real                                  :: beta
    !---------------------------------------------------------------------------

    allocate( x_best( dim ) )
    f_best = 1.E30

    ! initialize the parameters

    CALL cg_init (grad_tol, dim)
    IF ( step ) THEN
       alpha = gnorm
    END IF
    delta2 = 2.0 * delta - 1.0
    eta_sq = eta * eta
    iter = 0
    ck = 0
    qk = 0

    ! initial function and gradient evaluations, initial direction

    CALL cg_value (f, x, n)
    if ( f < f_best ) x_best(1:n) = x(1:n)
    if ( f < f_best ) f_best = f
    nf = nf + 1
    CALL cg_grad (g, x, n)
    ng = ng + 1
    f0 = f + f
    gnorm = zero
    xnorm = zero
    gnorm2 = zero
    DO i = 1, n5
       xnorm = MAX(xnorm, ABS(x(i)))
       t = g(i)
       d(i) = -t
       gnorm = MAX(gnorm, ABS(t))
       gnorm2 = gnorm2 + t*t
    END DO
    DO i = n6, n, 5
       xnorm = MAX(xnorm, ABS(x(i)))
       t = g(i)
       gnorm = MAX(gnorm, ABS(t))
       d(i)   = -t
       j = i + 1
       t1 = g(j)
       d(j) = -t1
       gnorm = MAX(gnorm, ABS(t1))
       xnorm = MAX(xnorm, ABS(x(j)))
       j = i + 2
       t2 = g(j)
       d(j) = -t2
       gnorm = MAX(gnorm, ABS(t2))
       xnorm = MAX(xnorm, ABS(x(j)))
       j = i + 3
       t3 = g(j)
       d(j) = -t3
       gnorm = MAX(gnorm, ABS(t3))
       xnorm = MAX(xnorm, ABS(x(j)))
       j = i + 4
       t4 = g(j)
       d(j) = -t4
       gnorm = MAX(gnorm, ABS(t4))
       xnorm = MAX(xnorm, ABS(x(j)))
       gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
    END DO

    IF ( stoprule ) THEN
       tol = MAX(gnorm*stopfac, tol)
    END IF

    IF ( printlevel ) THEN
       WRITE (*, 10) iter, f, gnorm, awolfe
10     FORMAT ('iter: ', i5, ' f= ', e14.6,  &
            ' gnorm= ', e14.6, ' AWolfe= ', l2)
    END IF

    IF ( cg_tol(f, gnorm) ) GO TO 100

    dphi0 = -gnorm2
    IF ( .NOT.step ) THEN
       alpha = psi0*xnorm/gnorm
       IF ( xnorm == zero ) THEN
          IF ( f /= zero ) THEN
             alpha = psi0*ABS(f)/gnorm2
          ELSE
             alpha = 1.0
          END IF
       END IF
    END IF

    ! start the conjugate gradient iteration


    !   alpha starts as old step, ends as initial step for next iteration
    !   f is function value for alpha = 0
    !   QuadOK = .true. means that a quadratic step was taken


    DO iter = 1, maxit
       quadok = .false.
       alpha = psi2*alpha
       IF ( quadstep ) THEN
          IF ( f /= zero ) THEN
             t = ABS((f-f0)/f)
          ELSE
             t = 1.0
          END IF
          IF ( t > quadcutoff ) THEN
             talpha = psi1*alpha
             CALL cg_step(xtemp, x, d, talpha)
             CALL cg_value(ftemp, xtemp, n)
             if ( ftemp < f_best ) x_best(1:n) = xtemp(1:n)
             if ( ftemp < f_best ) f_best = ftemp
             nf = nf + 1
             IF ( ftemp < f ) THEN
                denom = 2.0D0*(((ftemp-f)/talpha)-dphi0)
                IF ( denom > zero ) THEN
                   quadok = .true.
                   alpha = -dphi0*talpha/denom
                END IF
             END IF
          END IF
       END IF
       f0 = f

       IF ( printlevel ) THEN
          WRITE (*, 20) quadok, alpha, f0, dphi0
20        FORMAT ('QuadOK:', l2, ' initial a:',  &
               e14.6,' f0:', e14.6, ' dphi', e14.6)
       END IF

       ! parameters in Wolfe and approximiate Wolfe conditions, and in update

       qk = qdecay*qk + 1.
       ck = ck + (ABS(f) - ck)/qk

       IF ( pertrule ) THEN
          fpert = f + eps*ck
       ELSE
          fpert = f + eps
       END IF

       wolfe_hi = delta*dphi0
       wolfe_lo = sigma*dphi0
       awolfe_hi = delta2*dphi0
       IF ( awolfe ) THEN
          CALL cg_line (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,  &
               cg_value, cg_grad)
       ELSE
          CALL cg_linew (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,  &
               cg_value, cg_grad)
       END IF

       IF ( info > 0 ) GO TO 100

       ! Test for convergence to within machine epsilon
       ! (set feps to zero to remove this test)

       IF ( -alpha*dphi0 <= feps*ABS(f) ) THEN
          info = 1
          GO TO 100
       END IF

       ! compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g,

       IF ( MOD(iter, nrestart) /= 0 ) THEN
          gnorm = zero
          dnorm2 = zero
          yk2 = zero
          ykgk = zero
          DO i = 1, n5
             x(i) = xtemp(i)
             t = gtemp(i)
             yk = t - g(i)
             yk2 = yk2 + yk**2
             ykgk = ykgk + yk*t
             g(i) = t
             gnorm = MAX(gnorm, ABS(t))
             dnorm2 = dnorm2 + d(i)**2
          END DO
          DO i = n6, n, 5
             x(i) = xtemp(i)
             t = gtemp(i)
             yk = t - g(i)
             yk2 = yk2 + yk**2
             ykgk = ykgk + yk*t
             i1 = i + 1
             x(i1) = xtemp(i1)
             t1 = gtemp(i1)
             i2 = i + 2
             x(i2) = xtemp(i2)
             t2 = gtemp(i2)
             i3 = i + 3
             x(i3) = xtemp(i3)
             t3 = gtemp(i3)
             i4 = i + 4
             x(i4) = xtemp(i4)
             t4 = gtemp(i4)
             yk2 = yk2 + (t1-g(i1))**2 + (t2-g(i2))**2  &
                  + (t3-g(i3))**2 + (t4-g(i4))**2
             ykgk = ykgk + (t1-g(i1))*t1 + (t2-g(i2))*t2  &
                  + (t3-g(i3))*t3 + (t4-g(i4))*t4
             g(i) = t
             gnorm = MAX(gnorm, ABS(t))
             g(i1) = t1
             gnorm = MAX(gnorm, ABS(t1))
             g(i2) = t2
             gnorm = MAX(gnorm, ABS(t2))
             g(i3) = t3
             gnorm = MAX(gnorm, ABS(t3))
             g(i4) = t4
             gnorm = MAX(gnorm, ABS(t4))
             dnorm2 = dnorm2 + d(i)**2 + d(i1)**2 + d(i2)**2  &
                  + d(i3)**2 + d(i4)**2
          END DO
          IF ( cg_tol(f, gnorm) ) GO TO 100
          dkyk = dphi - dphi0
          beta = (ykgk - 2.0*dphi*yk2/dkyk)/dkyk

          !   faster: initialize dnorm2 = gnorm2 at start, then
          !             dnorm2 = gnorm2 + beta**2*dnorm2 - 2.0*beta*dphi
          !             gnorm2 = ||g_{k+1}||^2
          !             dnorm2 = ||d_{k+1}||^2
          !             dpi = g_{k+1}' d_k

          beta = MAX(beta, -1.0/SQRT(MIN(eta_sq, gnorm2)*dnorm2))

          !     update search direction d = -g + beta*dold

          gnorm2 = zero
          DO i = 1, n5
             t = g(i)
             d(i) = -t + beta*d(i)
             gnorm2 = gnorm2 + t*t
          END DO
          DO i = n6, n, 5
             d(i) = -g(i) + beta*d(i)
             i1 = i + 1
             d(i1) = -g(i1) + beta*d(i1)
             i2 = i + 2
             d(i2) = -g(i2) + beta*d(i2)
             i3 = i + 3
             d(i3) = -g(i3) + beta*d(i3)
             i4 = i + 4
             d(i4) = -g(i4) + beta*d(i4)
             gnorm2 = gnorm2 + g(i)**2 + g(i1)**2 + g(i2)**2  &
                  + g(i3)**2 + g(i4)**2
          END DO
          dphi0 = -gnorm2 + beta*dphi

       ELSE

          !     search direction d = -g

          IF ( printlevel ) THEN
             WRITE (*, *) "RESTART CG"
          END IF
          gnorm = zero
          gnorm2 = zero
          DO i = 1, n5
             x(i) = xtemp(i)
             t = gtemp(i)
             g(i) = t
             d(i) = -t
             gnorm = MAX(gnorm, ABS(t))
             gnorm2 = gnorm2 + t*t
          END DO
          DO i = n6, n, 5
             x(i) = xtemp(i)
             t = gtemp(i)
             g(i) = t
             d(i) = -t
             gnorm = MAX(gnorm, ABS(t))
             j = i + 1
             x(j) = xtemp(j)
             t1 = gtemp(j)
             g(j) = t1
             d(j) = -t1
             gnorm = MAX(gnorm, ABS(t1))
             j = i + 2
             x(j) = xtemp(j)
             t2 = gtemp(j)
             g(j) = t2
             d(j) = -t2
             gnorm = MAX(gnorm, ABS(t2))
             j = i + 3
             x(j) = xtemp(j)
             t3 = gtemp(j)
             g(j) = t3
             d(j) = -t3
             gnorm = MAX(gnorm, ABS(t3))
             j = i + 4
             x(j) = xtemp(j)
             t4 = gtemp(j)
             g(j) = t4
             d(j) = -t4
             gnorm = MAX(gnorm, ABS(t4))
             gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
          END DO
          IF ( cg_tol(f, gnorm) ) GO TO 100
          dphi0 = -gnorm2
       END IF
       IF ( .NOT.awolfe ) THEN
          IF ( ABS(f-f0) < awolfefac*ck ) THEN
             awolfe = .true.
          END IF
       END IF

       IF ( printlevel ) THEN
          WRITE (*, 10) iter, f, gnorm, awolfe
       END IF

       IF ( debug ) THEN
          IF ( f > f0 + 1.e-10*ck ) THEN
             WRITE (*, 270)
             WRITE (*, 50) f, f0
50           FORMAT (' new value:', e30.16, 'old value:', e30.16)
             STOP
          END IF
       END IF

       IF ( dphi0 > zero ) THEN
          info = 5
          GO TO 100
       END IF
    END DO
    info = 2
100 nfunc = nf
    ngrad = ng
    STATUS = info
    IF ( info > 2 ) THEN
       gnorm = zero
       DO i = 1, n
          x(i) = xtemp(i)
          g(i) = gtemp(i)
          gnorm = MAX(gnorm, ABS(g(i)))
       END DO
    END IF
    IF ( printfinal ) THEN
       WRITE (6, *) 'Termination status:', STATUS
       IF ( STATUS == 0 ) THEN
          WRITE (6, 200)
       ELSE IF ( STATUS == 1 ) THEN
          WRITE (6, 210)
       ELSE IF ( STATUS == 2 ) THEN
          WRITE (6, 220) maxit
          WRITE (6, 300)
          WRITE (6, 400) grad_tol
       ELSE IF ( STATUS == 3 ) THEN
          WRITE (6, 230)
          WRITE (6, 300)
          WRITE (6, 430)
          WRITE (6, 410)
       ELSE IF ( STATUS == 4 ) THEN
          WRITE (6, 240)
          WRITE (6, 300)
          WRITE (6, 400) grad_tol
       ELSE IF ( STATUS == 5 ) THEN
          WRITE (6, 250)
       ELSE IF ( STATUS == 6 ) THEN
          WRITE (6, 260)
          WRITE (6, 300)
          WRITE (6, 400) grad_tol
          WRITE (6, 410)
          WRITE (6, 420)
       ELSE IF ( STATUS == 7 ) THEN
          WRITE (6, 260)
          WRITE (6, 300)
          WRITE (6, 400) grad_tol
       ELSE IF ( STATUS == 8 ) THEN
          WRITE (6, 260)
          WRITE (6, 300)
          WRITE (6, 400) grad_tol
          WRITE (6, 410)
          WRITE (6, 420)
       END IF
       WRITE (6, 500) gnorm
       WRITE (6, *) 'function value:', f
       WRITE (6, *) 'cg iterations:', iter
       WRITE (6, *) 'function evaluations:', nfunc
       WRITE (6, *) 'gradient evaluations:', ngrad
    END IF

    if ( f_best < f ) then
       f = f_best
       x(1:dim) = x_best(1:dim)
    end if

    deallocate( x_best )

    RETURN

200 FORMAT (' Convergence tolerance for gradient satisfied')
210 FORMAT (' Terminating since change in function value <= feps*|f|')
220 FORMAT (' Total number of iterations exceed max allow:', i10)
230 FORMAT (' Slope always negative in line search')
240 FORMAT (' Line search fails, too many secant steps')
250 FORMAT (' Search direction not a descent direction')
260 FORMAT (' Line search fails')
270 FORMAT (' Debugger is on, function value does not improve')
300 FORMAT (' Possible causes of this error message:')
400 FORMAT ('   - your tolerance (grad_tol = ', d11.4,  &
         ') may be too strict')
410 FORMAT ('   - your gradient routine has an error')
420 FORMAT ('   - parameter epsilon in cg_descent.parm is too small')
430 FORMAT ('   - your cost function has an error')
500 FORMAT (' absolute largest component of gradient: ', d11.4)

  END SUBROUTINE cg_descent




  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !     PARAMETERS:

  !     delta - range (0, .5), used in the Wolfe conditions
  !     sigma - range [delta, 1), used in the Wolfe conditions
  !     eps - range [0, infty), used to compute line search perturbation
  !     gamma - range (0,1), determines when to perform bisection step
  !     rho   - range (1, infty), growth factor when finding initial interval
  !     eta   - range (0, infty), used in lower bound for beta
  !     psi0  - range (0, 1), factor used in very initial starting guess
  !     psi1  - range (0, 1), factor previous step multiplied by in QuadStep
  !     psi2  - range (1, infty), factor previous step is multipled by for startup
  !     QuadCutOff - perform QuadStep if relative change in f > QuadCutOff
  !     StopFac - used in StopRule
  !     AWolfeFac - used to decide when to switch from Wolfe to AWolfe
  !     restart_fac - range (0, infty) restart cg when iter = n*restart
  !     maxit_fac - range (0, infty) terminate in maxit = maxit_fac*n iterations
  !     feps - stop when -alpha*dphi0 (est. change in value) <= feps*|f|
  !            (feps = 0 removes this test, example: feps = eps*1.e-5
  !             where eps is machine epsilon)
  !     tol   - range (0, infty), convergence tolerance
  !     nexpand - range [0, infty), number of grow/shrink allowed in bracket
  !     nsecant - range [0, infty), maximum number of secant steps
  !     PertRule - gives the rule used for the perturbation in f
  !                 F => fpert = eps
  !                 T => fpert = eps*Ck, Ck is an average of prior |f|
  !                             Ck is an average of prior |f|
  !     QuadStep- .true. (use quadratic step) .false. (no quadratic step)
  !     PrintLevel- .false. (no printout) .true. (print intermediate results)
  !     PrintFinal- .false. (no printout) .true. (print messages, final error)
  !     StopRule - .true. (max abs grad <= max (tol, StopFac*initial abs grad))
  !                .false. (... <= tol*(1+|f|))
  !     AWolfe - .false. (use standard Wolfe initially)
  !            - .true. (use approximate + standard Wolfe)
  !     Step - .false. (program computing starting step at iteration 0)
  !          - .true. (user provides starting step in gnorm argument of cg_descent
  !     debug - .false. (no debugging)
  !           - .true. (check that function values do not increase)
  !     info  - same as status

  !     DEFAULT PARAMETER VALUES:

  !         delta : 0.1
  !         sigma : 0.9
  !         eps : 1.e-6
  !         gamma : 0.66
  !         rho   : 5.0
  !         restart: 1.0
  !         eta   : 0.01
  !         psi0  : 0.01
  !         psi1  : 0.1
  !         psi2  : 2.0
  !         QuadCutOff: 1.E-12
  !         StopFac: 0.0
  !         AWolfeFac: 1.E-3
  !         tol   : grad_tol
  !         nrestart: n (restart_fac = 1)
  !         maxit : 500*n (maxit_fac = 500)
  !         feps : 0.0
  !         Qdecay : 0.7
  !         nexpand: 50
  !         nsecant: 50
  !         PertRule: .true.
  !         QuadStep: .true.
  !         PrintLevel: .false.
  !         PrintFinal: .true.
  !         StopRule: .true.
  !         AWolfe: .false.
  !         Step: .false.
  !         debug: .false.
  !         info  : 0
  !         feps  : 0.0


  !      (double) grad_tol-- used in stopping rule
  !      (int)    dim     --problem dimension (also denoted n)
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE cg_init (grad_tol, dim)

    !  use define_cg_variables
    implicit none

    real, INTENT(IN)                      :: grad_tol
    integer, INTENT(IN)                   :: dim

    real                                  :: restart_fac, maxit_fac
    !---------------------------------------------------------------------------

    n = dim
    tol = grad_tol
    OPEN (10, FILE='cg_descent_f.parm')
    READ (10, *) delta
    READ (10, *) sigma
    READ (10, *) eps
    READ (10, *) gamma
    READ (10, *) rho
    READ (10, *) eta
    READ (10, *) psi0
    READ (10, *) psi1
    READ (10, *) psi2
    READ (10, *) quadcutoff
    READ (10, *) stopfac
    READ (10, *) awolfefac
    READ (10, *) restart_fac
    READ (10, *) maxit_fac
    READ (10, *) feps
    READ (10, *) qdecay
    READ (10, *) nexpand
    READ (10, *) nsecant
    READ (10, *) pertrule
    READ (10, *) quadstep
    READ (10, *) printlevel
    READ (10, *) printfinal
    READ (10, *) stoprule
    READ (10, *) awolfe
    READ (10, *) step
    READ (10, *) debug
    nrestart = n*restart_fac
    maxit = n*maxit_fac
    zero = 0.0
    info = 0
    n5 = MOD(n, 5)
    n6 = n5 + 1
    nf = 0
    ng = 0
    CLOSE (10)

  END SUBROUTINE cg_init



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! check whether the Wolfe or the approximate Wolfe conditions
  !     are satisfied

  !      (double) alpha   -- stepsize
  !      (double) f       -- function value associated with stepsize alpha
  !      (double) dphi    -- derivative value associated with stepsize alpha
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  LOGICAL FUNCTION cg_wolfe (alpha, f, dphi)

    !  use define_cg_variables
    implicit none

    real, INTENT(IN)             :: alpha
    real, INTENT(IN)             :: f
    real, INTENT(IN)             :: dphi
    !---------------------------------------------------------------------------

    IF ( dphi >= wolfe_lo ) THEN

       ! test original Wolfe conditions

       IF ( f-f0 <= alpha*wolfe_hi ) THEN
          cg_wolfe = .true.
          IF ( printlevel ) THEN
             WRITE (*, 10) f, f0, alpha*wolfe_hi, dphi
10           FORMAT (' wolfe f:', e14.6, ' f0: ',  &
                  e14.6, e14.6, ' dphi:', e14.6)
          END IF
          RETURN

          ! test approximate Wolfe conditions

       ELSE IF ( awolfe ) THEN
          IF ( (f <= fpert).AND.(dphi <= awolfe_hi) ) THEN
             cg_wolfe = .true.
             IF ( printlevel ) THEN
                WRITE (*, 20) f, fpert, dphi, awolfe_hi
20              FORMAT ('f:', e14.6, ' fpert:', e14.6,  &
                     ' dphi: ', e14.6, ' fappx:', e14.6)
             END IF
             RETURN
          END IF
       END IF
    END IF
    cg_wolfe = .false.

  END FUNCTION cg_wolfe



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! check for convergence of the cg iterations
  !      (double) f       -- function value associated with stepsize
  !      (double) gnorm   -- gradient (infinity) norm
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  LOGICAL FUNCTION cg_tol (f, gnorm)

    !  use define_cg_variables
    implicit none

    real, INTENT(IN OUT)         :: f
    real, INTENT(IN)             :: gnorm
    !---------------------------------------------------------------------------

    IF ( stoprule ) THEN
       IF ( gnorm <= tol ) THEN
          cg_tol = .true.
          RETURN
       END IF
    ELSE
       IF ( gnorm <= tol*(1.0 + ABS(f)) ) THEN
          cg_tol = .true.
          RETURN
       END IF
    END IF
    cg_tol = .false.

  END FUNCTION cg_tol



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! compute dot product of x and y, vectors of length n
  !      (double) x       -- first vector
  !      (double) y       -- second vector
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  real FUNCTION cg_dot (x, y)

    !  use define_cg_variables
    implicit none

    real, INTENT(IN)                      :: x(:)
    real, INTENT(IN)                      :: y(:)
    real                                  :: t

    integer                               :: i
    !---------------------------------------------------------------------------

    t = zero
    DO i = 1, n5
       t = t + x(i)*y(i)
    END DO
    DO i = n6, n, 5
       t = t + x(i)*y(i) + x(i+1)*y(i+1) + x(i+2)*y(i+2)  &
            + x(i+3)*y(i+3) + x(i+4)*y(i+4)
    END DO
    cg_dot = t

  END FUNCTION cg_dot


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  compute xtemp = x + alpha d

  !     (double) xtemp   -- output vector
  !     (double) x       -- initial vector
  !     (double) d       -- search direction vector
  !     (double) alpha   -- stepsize along search direction vector
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE cg_step (xtemp, x, d, alpha)

    !  use define_cg_variables
    implicit none

    real, INTENT(OUT)            :: xtemp(:)
    real, INTENT(IN)             :: x(:)
    real, INTENT(IN)             :: d(:)
    real, INTENT(IN)             :: alpha

    integer :: i, j
    !---------------------------------------------------------------------------

    DO i = 1, n5
       xtemp(i) = x(i) + alpha*d(i)
    END DO
    DO i = n6, n, 5
       xtemp(i) = x(i) + alpha*d(i)
       j = i + 1
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 2
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 3
       xtemp(j) = x(j) + alpha*d(j)
       j = i + 4
       xtemp(j) = x(j) + alpha*d(j)
    END DO
  END SUBROUTINE cg_step



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !      (double) alpha   -- stepsize along search direction vector
  !      (double) phi     -- function value for step alpha
  !      (double) dphi    -- function derivative for step alpha
  !      (double) dphi0   -- function derivative at starting point (alpha = 0)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE cg_line (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp,  &
       cg_value, cg_grad)

    !  use define_cg_variables
    implicit none


    real, INTENT(IN OUT)         :: alpha
    real, INTENT(IN OUT)         :: phi
    real, INTENT(OUT)            :: dphi
    real, INTENT(IN)             :: dphi0
    real, INTENT(IN OUT)         :: x(:)
    real, INTENT(IN OUT)         :: xtemp(:)
    real, INTENT(IN OUT)         :: d(:)
    real, INTENT(IN OUT)         :: gtemp(:)

    INTERFACE

       SUBROUTINE cg_value (f, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: f
       END SUBROUTINE cg_value

       SUBROUTINE cg_grad (g, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: g(:)
       END SUBROUTINE cg_grad

    END INTERFACE

    real                                  :: a, dphia, b, dphib, c
    real                                  :: a0, da0, b0, db0
    real                                  :: width, fquad

    integer                               :: ngrow, nshrink, iter, flag
    !---------------------------------------------------------------------------

    CALL cg_step (xtemp, x, d, alpha)
    CALL cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)

    ! Find initial interval [a,b] such that dphia < 0, dphib >= 0,
    !        and phia <= phi0 + feps*dabs (phi0)

    a = zero
    dphia = dphi0
    ngrow = 0
    nshrink = 0
    DO WHILE ( dphi < zero )
       CALL cg_value (phi, xtemp, n)
       nf = nf + 1

       ! if quadstep in effect and quadratic conditions hold, check wolfe condition

       IF ( quadok ) THEN
          IF ( ngrow == 0 ) fquad = MIN(phi, f0)
          IF ( phi <= fquad ) THEN
             IF ( printlevel ) THEN
                WRITE (*, 10) alpha, phi, fquad
10              FORMAT ('alpha:', e14.6, ' phi:', e14.6,  &
                     ' fquad:', e14.6)
             END IF
             IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
          END IF
       END IF
       IF ( phi <= fpert ) THEN
          a = alpha
          dphia = dphi
       ELSE

          ! contraction phase

          b = alpha
          DO WHILE ( .true. )
             alpha = .5D0*(a+b)
             nshrink = nshrink + 1
             IF ( nshrink > nexpand ) THEN
                info = 6
                RETURN
             END IF
             CALL cg_step(xtemp, x, d, alpha)
             CALL cg_grad(gtemp, xtemp, n)
             ng = ng + 1
             dphi = cg_dot(gtemp, d)
             IF ( dphi >= zero ) GO TO 100
             CALL cg_value(phi, xtemp, n)
             nf = nf + 1
             IF ( printlevel ) THEN
                WRITE (6, 20) a, b, alpha, phi, dphi
20              FORMAT ('contract, a:', e14.6,  &
                     ' b:', e14.6, ' alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
             END IF
             IF ( quadok .AND. (phi <= fquad) ) THEN
                IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
             END IF
             IF ( phi <= fpert ) THEN
                a = alpha
                dphia = dphi
             ELSE
                b = alpha
             END IF
          END DO
       END IF

       ! expansion phase

       ngrow = ngrow + 1
       IF ( ngrow > nexpand ) THEN
          info = 3
          RETURN
       END IF
       alpha = rho*alpha
       CALL cg_step(xtemp, x, d, alpha)
       CALL cg_grad(gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       IF ( printlevel ) THEN
          WRITE (*, 30) a, alpha, phi, dphi
30        FORMAT ('expand,   a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       END IF
    END DO
100 CONTINUE
    b = alpha
    dphib = dphi
    IF ( quadok ) THEN
       CALL cg_value(phi, xtemp, n)
       nf = nf + 1
       IF ( ngrow + nshrink == 0 ) fquad = MIN(phi, f0)
       IF ( phi <= fquad ) THEN
          IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
       END IF
    END IF
    DO iter = 1, nsecant
       IF ( printlevel ) THEN
          WRITE (*, 40) a, b, dphia, dphib
40        FORMAT ('secant, a:', e14.6, ' b:', e14.6,  &
               ' da:', e14.6, ' db:', e14.6)
       END IF
       width = gamma*(b - a)
       IF ( -dphia <= dphib ) THEN
          alpha = a - (a-b)*(dphia/(dphia-dphib))
       ELSE
          alpha = b - (a-b)*(dphib/(dphia-dphib))
       END IF
       c = alpha
       a0 = a
       b0 = b
       da0 = dphia
       db0 = dphib
       flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
            dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
       IF ( flag > 0 ) THEN
          RETURN
       ELSE IF ( flag == 0 ) THEN
          IF ( c == a ) THEN
             IF ( dphi > da0 ) THEN
                alpha = c - (c-a0)*(dphi/(dphi-da0))
             ELSE
                alpha = a
             END IF
          ELSE
             IF ( dphi < db0 ) THEN
                alpha = c - (c-b0)*(dphi/(dphi-db0))
             ELSE
                alpha = b
             END IF
          END IF
          IF ( (alpha > a) .AND. (alpha < b) ) THEN
             IF ( printlevel ) WRITE (*, *) "2nd secant"
             flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
                  dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
             IF ( flag > 0 ) RETURN
          END IF
       END IF

       !    bisection iteration

       IF ( (b-a) >= width ) THEN
          alpha = .5D0*(b+a)
          IF ( printlevel ) WRITE (*, *) "bisection"
          flag = cg_update(a, dphia, b, dphib, alpha, phi,  &
               dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
          IF ( flag > 0 ) RETURN
       ELSE
          IF ( b <= a ) THEN
             info = 7
             RETURN
          END IF
       END IF
    END DO
    info = 4

  END SUBROUTINE cg_line


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  This routine is identical to cg_line except that the function
  !  psi (a) = phi (a) - phi (0) - a*delta*dphi(0) is miniminized instead of
  !  the function phi

  !      (double) alpha   -- stepsize along search direction vector
  !      (double) phi     -- function value for step alpha
  !      (double) dphi    -- function derivative for step alpha
  !      (double) dphi0   -- function derivative at starting point (alpha = 0)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE cg_linew (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    real, INTENT(IN OUT)         :: alpha
    real, INTENT(IN OUT)             :: phi
    real, INTENT(OUT)            :: dphi
    real, INTENT(IN)             :: dphi0
    real, INTENT(IN OUT)         :: x(:)
    real, INTENT(IN OUT)         :: xtemp(:)
    real, INTENT(IN OUT)         :: d(:)
    real, INTENT(IN OUT)         :: gtemp(:)

    INTERFACE

       SUBROUTINE cg_value (f, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: f
       END SUBROUTINE cg_value

       SUBROUTINE cg_grad (g, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: g(:)
       END SUBROUTINE cg_grad

    END INTERFACE

    real :: a, dpsia, b, dpsib, c,  &
         a0, da0, b0, db0, width, fquad, psi, dpsi

    INTEGER :: ngrow, nshrink, iter, flag
    !---------------------------------------------------------------------------


    CALL cg_step (xtemp, x, d, alpha)
    CALL cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)
    dpsi = dphi - wolfe_hi

    ! Find initial interval [a,b] such that dpsia < 0, dpsib >= 0,
    !        and psia <= phi0 + feps*dabs(phi0)

    a = zero
    dpsia = dphi0 - wolfe_hi
    ngrow = 0
    nshrink = 0
    DO WHILE ( dpsi < zero )
       CALL cg_value(phi, xtemp, n)
       psi = phi - alpha*wolfe_hi

       nf = nf + 1

       ! if quadstep in effect and quadratic conditions hold, check wolfe condition

       IF ( quadok ) THEN
          IF ( ngrow == 0 ) fquad = MIN(phi, f0)
          IF ( phi <= fquad ) THEN
             IF ( printlevel ) THEN
                WRITE (*, 10) alpha, phi, fquad
10              FORMAT ('alpha:', e14.6, ' phi:', e14.6,  &
                     ' fquad:', e14.6)
             END IF
             IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
          END IF
       END IF
       IF ( psi <= fpert ) THEN
          a = alpha
          dpsia = dpsi
       ELSE

          ! contraction phase

          b = alpha
          DO WHILE ( .true. )
             alpha = .5D0*(a+b)
             nshrink = nshrink + 1
             IF ( nshrink > nexpand ) THEN
                info = 6
                RETURN
             END IF
             CALL cg_step (xtemp, x, d, alpha)
             CALL cg_grad (gtemp, xtemp, n)
             ng = ng + 1
             dphi = cg_dot(gtemp, d)
             dpsi = dphi - wolfe_hi
             IF ( dpsi >= zero ) GO TO 100
             CALL cg_value (phi, xtemp, n)
             psi = phi - alpha*wolfe_hi
             nf = nf + 1
             IF ( printlevel ) THEN
                WRITE (6, 20) a, b, alpha, phi, dphi
20              FORMAT ('contract, a:', e14.6,  &
                     ' b:', e14.6, ' alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
             END IF
             IF ( quadok .AND. (phi <= fquad) ) THEN
                IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
             END IF
             IF ( psi <= fpert ) THEN
                a = alpha
                dpsia = dpsi
             ELSE
                b = alpha
             END IF
          END DO
       END IF

       ! expansion phase

       ngrow = ngrow + 1
       IF ( ngrow > nexpand ) THEN
          info = 3
          RETURN
       END IF
       alpha = rho*alpha
       CALL cg_step (xtemp, x, d, alpha)
       CALL cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       dpsi = dphi - wolfe_hi
       IF ( printlevel ) THEN
          WRITE (*, 30) a, alpha, phi, dphi
30        FORMAT ('expand,   a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
          WRITE (6, *) "expand, alpha:", alpha, "dphi:", dphi
       END IF
    END DO
100 CONTINUE
    b = alpha
    dpsib = dpsi
    IF ( quadok ) THEN
       CALL cg_value(phi, xtemp, n)
       nf = nf + 1
       IF ( ngrow + nshrink == 0 ) fquad = MIN(phi, f0)
       IF ( phi <= fquad ) THEN
          IF ( cg_wolfe(alpha, phi, dphi) ) RETURN
       END IF
    END IF
    DO iter = 1, nsecant
       IF ( printlevel ) THEN
          WRITE (*, 40) a, b, dpsia, dpsib
40        FORMAT ('secant, a:', e14.6, ' b:', e14.6,  &
               ' da:', e14.6, ' db:', e14.6)
       END IF
       width = gamma*(b - a)
       IF ( -dpsia <= dpsib ) THEN
          alpha = a - (a-b)*(dpsia/(dpsia-dpsib))
       ELSE
          alpha = b - (a-b)*(dpsib/(dpsia-dpsib))
       END IF
       c = alpha
       a0 = a
       b0 = b
       da0 = dpsia
       db0 = dpsib
       flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
            phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
       IF ( flag > 0 ) THEN
          RETURN
       ELSE IF ( flag == 0 ) THEN
          IF ( c == a ) THEN
             IF ( dpsi > da0 ) THEN
                alpha = c - (c-a0)*(dpsi/(dpsi-da0))
             ELSE
                alpha = a
             END IF
          ELSE
             IF ( dpsi < db0 ) THEN
                alpha = c -(c-b0)*(dpsi/(dpsi-db0))
             ELSE
                alpha = b
             END IF
          END IF
          IF ( (alpha > a) .AND. (alpha < b) ) THEN
             IF ( printlevel ) WRITE (*, *) "2nd secant"
             flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
                  phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
             IF ( flag > 0 ) RETURN
          END IF
       END IF

       !    bisection iteration

       IF ( (b-a) >= width ) THEN
          alpha = .5D0*(b+a)
          IF ( printlevel ) WRITE (*, *) "bisection"
          flag = cg_updatew(a, dpsia, b, dpsib, alpha,  &
               phi, dphi, dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)
          IF ( flag > 0 ) RETURN
       ELSE
          IF ( b <= a ) THEN
             info = 7
             RETURN
          END IF
       END IF
    END DO
    info = 4

  END SUBROUTINE cg_linew


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  ! update returns 1 if Wolfe condition is satisfied or too many iterations
  !        returns  0 if the interval updated successfully
  !        returns -1 if search done

  !      (double) a       -- left side of bracketting interval
  !      (double) dphia   -- derivative at a
  !      (double) b       -- right side of bracketting interval
  !      (double) dphib   -- derivative at b
  !      (double) alpha   -- trial step (between a and b)
  !      (double) phi     -- function value at alpha (returned)
  !      (double) dphi    -- function derivative at alpha (returned)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  INTEGER FUNCTION cg_update (a, dphia, b, dphib, alpha, phi,  &
       dphi, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    real, INTENT(OUT)            :: a
    real, INTENT(OUT)            :: dphia
    real, INTENT(OUT)            :: b
    real, INTENT(OUT)            :: dphib
    real, INTENT(IN OUT)         :: alpha
    real, INTENT(IN OUT)         :: phi
    real, INTENT(OUT)            :: dphi
    real, INTENT(IN OUT)         :: x(:)
    real, INTENT(IN OUT)         :: xtemp(:)
    real, INTENT(IN OUT)         :: d(:)
    real, INTENT(IN OUT)         :: gtemp(:)

    INTERFACE

       SUBROUTINE cg_value (f, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: f
       END SUBROUTINE cg_value

       SUBROUTINE cg_grad (g, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: g(:)
       END SUBROUTINE cg_grad

    END INTERFACE


    INTEGER :: nshrink
    !---------------------------------------------------------------------------

    CALL cg_step (xtemp, x, d, alpha)
    CALL cg_value (phi, xtemp, n)
    nf = nf + 1
    CALL cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot(gtemp, d)
    IF ( printlevel ) THEN
       WRITE (*, 10) alpha, phi, dphi
10     FORMAT ('update alpha:', e14.6, ' phi:', e14.6, ' dphi:', e14.6)
    END IF
    cg_update = 0
    IF ( cg_wolfe(alpha, phi, dphi) ) THEN
       cg_update = 1
       GO TO 110
    END IF
    IF ( dphi >= zero ) THEN
       b = alpha
       dphib = dphi
       GO TO 110
    ELSE
       IF ( phi <= fpert ) THEN
          a = alpha
          dphia = dphi
          GO TO 110
       END IF
    END IF
    nshrink = 0
    b = alpha
    DO WHILE ( .true. )
       alpha = .5D0*(a+b)
       nshrink = nshrink + 1
       IF ( nshrink > nexpand ) THEN
          info = 8
          cg_update = 1
          GO TO 110
       END IF
       CALL cg_step (xtemp, x, d, alpha)
       CALL cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       CALL cg_value (phi, xtemp, n)
       nf = nf + 1
       IF ( printlevel ) THEN
          WRITE (6, 20) a, alpha, phi, dphi
20        FORMAT ('contract, a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       END IF
       IF ( cg_wolfe(alpha, phi, dphi) ) THEN
          cg_update = 1
          GO TO 110
       END IF
       IF ( dphi >= zero ) THEN
          b = alpha
          dphib = dphi
          GO TO 100
       END IF
       IF ( phi <= fpert ) THEN
          IF ( printlevel ) THEN
             WRITE (6, *) "update a:", alpha, "dphia:", dphi
          END IF
          a = alpha
          dphia = dphi
       ELSE
          b = alpha
       END IF
    END DO
100 CONTINUE
    cg_update = -1
110 CONTINUE
    IF ( printlevel ) THEN
       WRITE (*, 200) a, b, dphia, dphib, cg_update
200    FORMAT ('UP a:', e14.6, ' b:', e14.6,  &
            ' da:', e14.6, ' db:', e14.6, ' up:', i2)
    END IF

  END FUNCTION cg_update


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !  This routine is identical to cg_update except that the function
  !  psi (a) = phi(a) - phi(0) - a*delta*dphi(0) is miniminized instead of
  !  the function phi

  ! update returns 1 if Wolfe condition is satisfied or too many iterations
  !        returns  0 if the interval updated successfully
  !        returns -1 if search done

  !      (double) a       -- left side of bracketting interval
  !      (double) dpsia   -- derivative at a
  !      (double) b       -- right side of bracketting interval
  !      (double) dpsib   -- derivative at b
  !      (double) alpha   -- trial step (between a and b)
  !      (double) phi     -- function value at alpha(returned)
  !      (double) dphi    -- derivative of phi at alpha(returned)
  !      (double) dpsi    -- derivative of psi at alpha(returned)
  !      (double) x       -- current iterate
  !      (double) xtemp   -- x + alpha*d
  !      (double) d       -- current search direction
  !      (double) gtemp   -- gradient at x + alpha*d
  !      (external) cg_value -- routine to evaluate function value
  !      (external) cg_grad  -- routine to evaluate function gradient
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  INTEGER FUNCTION cg_updatew (a, dpsia, b, dpsib, alpha, phi, dphi,  &
       dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)

    !  use define_cg_variables
    implicit none

    real, INTENT(OUT)            :: a
    real, INTENT(OUT)            :: dpsia
    real, INTENT(OUT)            :: b
    real, INTENT(OUT)            :: dpsib
    real, INTENT(IN OUT)         :: alpha
    real, INTENT(IN OUT)         :: phi
    real, INTENT(OUT)            :: dphi
    real, INTENT(OUT)            :: dpsi
    real, INTENT(IN OUT)         :: x(:)
    real, INTENT(IN OUT)         :: xtemp(:)
    real, INTENT(IN OUT)         :: d(:)
    real, INTENT(IN OUT)         :: gtemp(:)

    INTERFACE

       SUBROUTINE cg_value (f, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: f
       END SUBROUTINE cg_value

       SUBROUTINE cg_grad (g, x, n)
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: n
         REAL, INTENT(IN)           :: x(:)
         REAL, INTENT(IN OUT)       :: g(:)
       END SUBROUTINE cg_grad

    END INTERFACE

    real :: psi

    INTEGER :: nshrink
    !---------------------------------------------------------------------------

    CALL cg_step (xtemp, x, d, alpha)
    CALL cg_value (phi, xtemp, n)
    psi = phi - alpha*wolfe_hi
    nf = nf + 1
    CALL cg_grad (gtemp, xtemp, n)
    ng = ng + 1
    dphi = cg_dot (gtemp, d)
    dpsi = dphi - wolfe_hi
    IF ( printlevel ) THEN
       WRITE (*, 10) alpha, psi, dpsi
10     FORMAT ('update alpha:', e14.6, ' psi:', e14.6, ' dpsi:', e14.6)
    END IF
    cg_updatew = 0
    IF ( cg_wolfe(alpha, phi, dphi) ) THEN
       cg_updatew = 1
       GO TO 110
    END IF
    IF ( dpsi >= zero ) THEN
       b = alpha
       dpsib = dpsi
       GO TO 110
    ELSE
       IF ( psi <= fpert ) THEN
          a = alpha
          dpsia = dpsi
          GO TO 110
       END IF
    END IF
    nshrink = 0
    b = alpha
    DO WHILE ( .true. )
       alpha = .5D0*(a+b)
       nshrink = nshrink + 1
       IF ( nshrink > nexpand ) THEN
          info = 8
          cg_updatew = 1
          GO TO 110
       END IF
       CALL cg_step (xtemp, x, d, alpha)
       CALL cg_grad (gtemp, xtemp, n)
       ng = ng + 1
       dphi = cg_dot(gtemp, d)
       dpsi = dphi - wolfe_hi
       CALL cg_value (phi, xtemp, n)
       psi = phi - alpha*wolfe_hi
       nf = nf + 1
       IF ( printlevel ) THEN
          WRITE (6, 20) a, alpha, phi, dphi
20        FORMAT ('contract, a:', e14.6, ' alpha:', e14.6,  &
               ' phi:', e14.6, ' dphi:', e14.6)
       END IF
       IF ( cg_wolfe(alpha, phi, dphi) ) THEN
          cg_updatew = 1
          GO TO 110
       END IF
       IF ( dpsi >= zero ) THEN
          b = alpha
          dpsib = dpsi
          GO TO 100
       END IF
       IF ( psi <= fpert ) THEN
          IF ( printlevel ) THEN
             WRITE (6, *) "update a:", alpha, "dpsia:", dpsi
          END IF
          a = alpha
          dpsia = dpsi
       ELSE
          b = alpha
       END IF
    END DO
100 CONTINUE
    cg_updatew = -1
110 CONTINUE
    IF ( printlevel ) THEN
       WRITE (*, 200) a, b, dpsia, dpsib, cg_updatew
200    FORMAT ('UP a:', e14.6, ' b:', e14.6,  &
            ' da:', e14.6, ' db:', e14.6, ' up:', i2)
    END IF

  END FUNCTION cg_updatew
  ! Version 1.2 Changes:

  !   1. Fix problem with user specified initial step (overwriting step)
  !   2. Change dphi to dpsi at lines 1228 and 1234 in cg_lineW
  !   3. Add comment about how to compute dnorm2 by an update of previous dnorm2
  !   4. In comment statements for cg_lineW and cg_updateW, insert "delta"
  !      in definition of psi (a)
  !   5. In dimension statements, change "(1)" to "(*)"

  ! Version 1.3 Changes:
  !   1. Remove extraneous write in line 985 (same thing written out twice)
  !   2. Remove the parameter theta from cg_descent.parm and from the code
  !      (we use theta = .5 in the cg_update)

  ! Version 1.4 Change:
  !   1. The variable dpsi needs to be included in the argument list for
  !      subroutine updateW (update of a Wolfe line search)

END MODULE cg_minimization
